##' @include abstracts.R parameters.R markov.annualHomogeneous.R
##' @export

nstates = 2
parameter.length=c(nstates,nstates,1,1,1,1,1,1,nstates-1)
markov.annualNonHomogeneous <- setClass(
  # Set the name for the class
  "markov.annualNonHomogeneous",

  package='hydroState',

  contains=c('markov.annualHomogeneous'),

  # Define the slots
  slots = c(
    inputForcing= "matrix",
    nstates = "numeric",  # Add a slot for nstates
    parameter.length = "numeric",  # Add a slot for parameter.length to be set later
    do.Logistic.Displacement = "logical"
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    transition.graph = matrix(TRUE,2,2),
    transition.graph.parameter.index = matrix(c(1,-1,2,-1),nrow = 2,ncol=2),
    parameter.length=c(nstates,nstates,1,1,1,1,1,1,nstates-1),
    parameters = new('parameters',c('transition.prob.amp','transition.prob.lower',
                                    'transition.prob.slope', 'transition.prob.center','transition.prob.disp','transition.prob.pearson.A',
                                    'transition.prob.pearson.B','transition.prob.pearson.N',
                                    'initial.state.prob'),parameter.length),
    inputForcing = matrix(NA,10,3),
    do.Logistic.Displacement = F

  )
)


# Valid object?
validObject <- function(object) {
  return(TRUE)
}
setValidity("markov.annualNonHomogeneous", validObject)

# Initialise object
#setGeneric(name="initialize",def=function(.Object,transition.graph) {standardGeneric("initialize")})
setMethod(f="initialize",
          signature="markov.annualNonHomogeneous",
          definition = function(.Object, transition.graph = matrix(TRUE, 2, 2), do.Logistic.Displacement=F,inputForcing = matrix(NA, 10, 3))
          {
            .Object@transition.graph <- transition.graph

            # Set nstates
            nstates = ncol(.Object@transition.graph)

            # Chec nstates ==2
            if (nstates!=2)
              stop('Only two states can be examined using the non-homogenous model.')

            # Set te number of transition parameters and a matrxi of indexes to them.
            #
            #.Object <- initialize.TransitionParameters(.Object)
            parameter.length=c(nstates,nstates,1,1,1,1,1,1,nstates-1)
            .Object@parameters = new('parameters',c('transition.prob.amp','transition.prob.lower',
                                                    'transition.prob.slope', 'transition.prob.center', 'transition.prob.disp','transition.prob.pearson.A',
                                                    'transition.prob.pearson.B','transition.prob.pearson.N',
                                                    'initial.state.prob'),parameter.length)

            # Set the initial stat probs.
            .Object <- initialize.StateProbParameters(.Object)

            .Object@inputForcing = inputForcing
            .Object@do.Logistic.Displacement = do.Logistic.Displacement
            validObject(.Object)
            return(.Object)
          }
)





#' # Set the number of transition parameters and how they relate to the transition graph.
setMethod(
  f="getTransitionForcing",
  signature ="markov.annualNonHomogeneous",
  definition = function(.Object) {

    P <- .Object@inputForcing[,2]
    PET <- .Object@inputForcing[,3]

    # P <- Qhat@input.data[["precipitation"]]
    # PET <- Qhat@input.data[["PET"]]


    d <- P - PET
    meand <- mean(d)
    d <- c(rep(meand, 50), d)

    parameters = getParameters(.Object@parameters)

    b <- parameters$transition.prob.pearson.B
    k <- parameters$transition.prob.pearson.N
    a <- parameters$transition.prob.pearson.A




    N <- length(d)
    f <- numeric(N)
    xpeak <- ((k - 1) * 1 / b)
    fpeak <- (a * (b^k) * (xpeak^(k - 1)) * exp(-b * xpeak))

    for (i in 1:N) {
      NN <- seq(i, 1, -1)
      w <- (a * (b^k) * (NN^(k - 1)) * exp(-b * NN)) / gamma(k)

      if (i == 1) {
        f[i] <- w[i] * d[i]
      } else if (i == 2) {
        f[i] <- w[1] * d[1] + w[i] * d[i]
      } else if (i == 3) {
        f[i] <- w[1] * d[1] + w[i] * d[i] + 3 * w[2] * d[2]
      } else {
        i_end <- ifelse(i %% 3 == 0, i - 3, i - i %% 3)
        ind3 <- seq(3, i_end, by = 3)
        ind12 <- 1:i
        ind12 <- ind12[!(ind12 %in% ind3)]
        ind12 <- ind12[!(ind12 %in% c(1, i))]

        f[i] <- w[1] * d[1] + w[i] * d[i] + sum(2 * w[ind3] * d[ind3]) + sum(3 * w[ind12] * d[ind12])
      }

      f[i] <- 3/8 * f[i]
    }

    f <- f[51:N]
    #message(paste('f:',f[51]))
    return(f)
  }
)

setGeneric("constraintx", def = function(.Object) {standardGeneric("constraintx")})

setMethod("constraintx", signature = "markov.annualNonHomogeneous",definition = function(.Object) {


  parameters = getParameters(.Object@parameters)

  forcingValues = getTransitionForcing(.Object)

  b <- params@transitionProbPearsonB
  k <- params@transitionProbPearsonN
  a <- params@transitionProbPearsonA
  xpeak <- (k - 1) / b

  if (xpeak > 10) {
    forcingValues[]<-Inf
  }
}
)

setGeneric("constraintdm", def = function(.Object) {standardGeneric("constraintdm")})

setMethod("constraintdm", signature = "markov.annualNonHomogeneous",definition = function(.Object) {

  parameters = getParameters(.Object@parameters)

  forcingValues = getTransitionForcing(.Object)

  b <- params@transitionProbPearsonB
  k <- params@transitionProbPearsonN
  a <- params@transitionProbPearsonA
  xpeak <- (k-1) / b

  if (b < 1/30*k) {
    forcingValues[]<-Inf
  }
}
)





setGeneric(
  "applyLogistic",
  def = function(.Object) {
    standardGeneric("applyLogistic")
  }
)


setMethod(
  "applyLogistic",
  signature = "markov.annualNonHomogeneous", # or "HydroModelParameters" if you're reusing that
  definition = function(.Object) {
    do.Logistic.Displacement = .Object@do.Logistic.Displacement
    forcingValues = getTransitionForcing(.Object)
    forcingRange<-(max(forcingValues, na.rm =TRUE) - min(forcingValues, na.rm = TRUE))
    parameters = getParameters(.Object@parameters)
    # Extract specific parameters
    cmdx <- parameters$transition.prob.center*forcingRange+min(forcingValues, na.rm = TRUE)
    slp <- parameters$transition.prob.slope

    N <- length(forcingValues)
    L <- numeric(N)

    # Apply the logistic equation
    L <- 1 / (1 + exp(-slp * (forcingValues - cmdx)))




    if (do.Logistic.Displacement==T){
      disp<-c(1,L[1:length(L)-1])
      cmdt <- parameters$transition.prob.disp*forcingRange
      cmd <- cmdx + cmdt*(1-disp)
      LD <- 1 / (1+exp(-slp*((forcingValues - cmd))))


    } else{

      LD <- L


    }

    Logi<- data.frame(L,LD)

    # Return logistic curve values
    return(Logi)
  }
)



#setGeneric( name = "getTransitionProbabilities", def = function(.Object, data) {standardGeneric("getTransitionProbabilities")})
# Get transition matrix with no input data.
setMethod(f="getTransitionProbabilities",
          signature="markov.annualNonHomogeneous",
          definition=function(.Object)
          {
            # Get number of states
            nStates = getNumStates(.Object)
            # Handle 1 state model.
            if (nStates==1) {
              return(c(1))
            }

            # Get object parameter vector
            parameters = getParameters(.Object@parameters)
            # Get weighted forcing
            forcingValues = getTransitionForcing(.Object)
            # Get logistic values and its length
            logisticValuesL = applyLogistic(.Object)$L
            logisticValuesLD = applyLogistic(.Object)$LD
            nT = length(forcingValues)
            do.Logistic.Displacement = .Object@do.Logistic.Displacement
            # # Check amp and lower sum to <=1
            if (any(parameters$transition.prob.amp + parameters$transition.prob.lower > 1)) {
              return(array(Inf,c(nStates,nStates, nT)))
            }

            # if (any(parameters$transition.prob.amp + parameters$transition.prob.lower > 1)) {
            #   parameters$transition.prob.amp[] <- Inf
            #   parameters$transition.prob.lower[] <- Inf
            #
            # }

            # Initials output
            Tprob = array(0.0,c(nStates,nStates, nT))

            # Cal. vector of transition prob for first time step from amp and lower for state 1.
            Tprob[1,1,] = logisticValuesL * parameters$transition.prob.amp[1] + parameters$transition.prob.lower[1]

            # Cal. vector of transition prob for first time step from amp and lower for state 2.
            # Tprob[1,2,] = (1-logisticValues) * parameters$transition.prob.amp[1] + parameters$transition.prob.lower[1]
            Tprob[1,2,] = 1 - Tprob[1,1,]

            # Cal. vector of transition prob for first time step from amp and lower for state 1.
            # TP: If displacement set by user, apply here to T22 and T21. This might require re-calling logistic fn

            if (do.Logistic.Displacement == T) {
              Tprob[2, 2, ] = (1 - logisticValuesLD) * parameters$transition.prob.amp[2] + parameters$transition.prob.lower[2]
            } else {
              Tprob[2, 2, ] = (1 - logisticValuesL) * parameters$transition.prob.amp[2] + parameters$transition.prob.lower[2]
            }
            Tprob[2, 1, ] = 1 - Tprob[2, 2, ]





            return(Tprob)
          }
)


# Get the log likelihood for the input data.
setMethod(f="getLogLikelihood", signature=c("markov.annualNonHomogeneous","data.frame","matrix"),
          definition=function(.Object, data, emission.probs)
          {
            # Check all the emmision probs. are finite.
            if (any(is.infinite(emission.probs)))
              return(Inf)

            # Get number of states
            nStates = getNumStates(.Object)

            # Check required fields exist
            if (!any(names(data)=='Qhat.flow'))
              stop('"data" must be a data.frame with the field "Qhat.flow".')

            # Built filter for non NAs.
            filt <- is.finite(data$Qhat.flow)

            # Handle 1 state model.
            if (nStates==1) {
              return(sum(log(emission.probs[filt])))
            }

            # Get the initial state probabilites
            alpha = getInitialStateProbabilities(.Object)
            # Check the initial state probs are sum to 1 and are all b/w 0 and 1.
            if (abs(sum(alpha)-1)>sqrt(.Machine$double.eps) || any(alpha<0) || any(alpha>1))
              return(Inf)

            # Get the transition matrix.
            Tprob = getTransitionProbabilities(.Object)
            if (any (Tprob[,,]  == Inf))
              return(Inf)

            if (length(dim(Tprob)) < 3) {
              stop("Tprob does not have the expected three dimensions")
            }

            # Only accept the parameters if the forward probabilities from the first to next time step show
            # that the state has not switched state after the first time step.
            P.forward <- getLogForwardProbabilities(.Object, data[1:2,], as.matrix(emission.probs[1:2,], ncol=nStates))
            if (all(!is.na(P.forward))) {
              if (which(max(P.forward[,1])==P.forward[,1]) != which(max(P.forward[,2])==P.forward[,2]))
                return(Inf)
            }


            # Set Emmision probs 1.0 if no obs data.
            emission.probs[!filt,] = 1

            # Adapted from Zucchini, McDonald and Langrock, 2016, Hidden Markov Models for Time Series, 2nd Edision, Appendix A.
            #---------------------------------------------------------------------

            # Get alpha at first time step and log liklihood
            alpha      <- alpha * emission.probs[1,]
            sumalpha   <- sum(alpha)
            lscale   <- log(sumalpha)
            alpha      <- alpha/sumalpha

            # Loop through 2+ time steps
            for (i in 2:nrow(data)) {
              alpha    <- alpha %*% Tprob[,,i] * as.vector(emission.probs[i,])
              sumalpha <- sum(alpha)
              lscale <- lscale+log(sumalpha)
              alpha    <- alpha/sumalpha
            }

            return(lscale)

          }
)

# Get the log forward probabilities for the input data.
#' @exportMethod getLogForwardProbabilities
#setGeneric(name="getLogForwardProbabilities",def=function(.Object, data, emission.probs) {standardGeneric("getLogForwardProbabilities")})
setMethod(f="getLogForwardProbabilities", signature=c("markov.annualNonHomogeneous","data.frame","matrix"),
          definition=function(.Object, data, emission.probs)
          {

            # Get number of states
            nStates = getNumStates(.Object)

            # Built filter for non NAs.
            filt <- !is.na(data$Qhat.flow)

            # Get the initial state probabilites
            alpha = getInitialStateProbabilities(.Object)

            # Get the transition matrix.
            Tprob = getTransitionProbabilities(.Object)
            f=getTransitionForcing(.Object)
            if (length(dim(Tprob)) < 3) {
              stop("Tprob does not have the expected three dimensions")
            }
            # Set Emmision probs 1.0 if no obs data.
            emission.probs[!filt,] = 1

            # Adapted from Zucchini, McDonald and Langrock, 2016, Hidden Markov Models for Time Series, 2nd Edision, Appendix A.
            #---------------------------------------------------------------------
            n             <- nrow(data)
            lalpha        <- matrix(NA,nStates,n)
            foo           <- alpha * as.vector(emission.probs[1,])
            sumfoo        <- sum(foo)
            lscale        <- log(sumfoo)
            foo           <- foo/sumfoo
            lalpha[,1]    <- lscale+log(foo)
            for (i in 2:n)
            {
              foo          <- foo %*% Tprob[,,i] * as.vector(emission.probs[i,])
              sumfoo       <- sum(foo)
              lscale       <- lscale+log(sumfoo)
              foo          <- foo/sumfoo
              lalpha[,i]   <- log(foo)+lscale
            }
            return(lalpha)
            #---------------------------------------------------------------------
          }
)

# Get the log backward probabilities for the input data.
#setGeneric(name="getLogBackwardProbabilities",def=function(.Object, data, emission.probs) {standardGeneric("getLogBackwardProbabilities")})
setMethod(f="getLogBackwardProbabilities", signature=c("markov.annualNonHomogeneous","data.frame","matrix"),
          definition=function(.Object, data, emission.probs)
          {

            # Get number of states
            nStates = getNumStates(.Object)

            # Built filter for non NAs.
            filt <- !is.na(data$Qhat.flow)

            # Get the initial state probabilites
            alpha = getInitialStateProbabilities(.Object)

            # Get the transition matrix.
            Tprob = getTransitionProbabilities(.Object)
            if (length(dim(Tprob)) < 3) {
              stop("Tprob does not have the expected three dimensions")
            }
            # Set Emmision probs 1.0 if no obs data.
            emission.probs[!filt,] = 1

            # Adapted from Zucchini, McDonald and Langrock, 2016, Hidden Markov Models for Time Series, 2nd Edision, Appendix A.
            #---------------------------------------------------------------------
            n          <- nrow(data)
            m          <- nStates
            lbeta      <- matrix(NA,m,n)
            lbeta[,n]  <- rep(0,m)
            foo        <- rep(1/m,m)
            lscale     <- log(m)
            for (i in (n-1):1)
            {
              foo        <- Tprob[,,i+1] %*% (as.vector(emission.probs[i+1,])*foo)
              lbeta[,i]  <- log(foo)+lscale
              sumfoo     <- sum(foo)
              foo        <- foo/sumfoo
              lscale     <- lscale+log(sumfoo)
            }
            return(lbeta)
            #---------------------------------------------------------------------
          }
)


# Get the conditional probability of a givn Qhat observation at time t given all other observations.
setMethod(f="getConditionalProbabilities", signature="markov.annualNonHomogeneous",
          definition=function(.Object, data, emission.probs, cumprob.atQhatIncrements)
          {
            # Get number of states
            nStates = getNumStates(.Object)

            # Get the initial state probabilites
            alpha = getInitialStateProbabilities(.Object)

            # Get the transition matrix.
            Tprob = getTransitionProbabilities(.Object)
            if (length(dim(Tprob)) < 3) {
              stop("Tprob does not have the expected three dimensions")
            }
            # Adapted from Zucchini, McDonald and Langrock, 2016, Hidden Markov Models for Time Series, 2nd Edision, Appendix A.
            #---------------------------------------------------------------------

            n          <- nrow(data)
            m          <- nStates
            nxc       <- dim(cumprob.atQhatIncrements)[3]
            message(paste('cumprob.atQhatIncrements :',cumprob.atQhatIncrements))
            message(paste('nxc :',nxc))
            dxc       <- matrix(NA,nrow=nxc,ncol=n)
            message(paste('dxc :',dxc))

            la        <- getLogForwardProbabilities(.Object, data, emission.probs)
            lb        <- getLogBackwardProbabilities(.Object, data, emission.probs)
            message(paste('la :',la))
            message(paste('lb :',lb))
            la        <- cbind(log(alpha),la)
            message(paste('la2 :',la))
            lafact    <- apply(la,2,max)
            lbfact    <- apply(lb,2,max)
            for (i in 1:n)
            {
              foo      <- (exp(la[,i]-lafact[i]) %*% Tprob[,,i])*exp(lb[,i]-lbfact[i])
              foo      <- foo/sum(foo)
              message(paste('foo :',foo))
              # Note, Zucchini, McDonald and Langrock, 2016 code calculates the probability of a given
              # Xt outside this for-loop because the emmision probs in their model are independent of
              # time. However, here the time-varying means and auto-regressive terms make the emmision probs
              # time-dependent and so here a 3D array is passed to this function with the depth dimension
              # containing the cumulative probs. at pre-defined incrementts of Qhat for a given state and time point.
              Px <- cumprob.atQhatIncrements[i,,]

              dxc[,i]  <- as.vector(foo%*%Px)

            }

            return(dxc)
            #---------------------------------------------------------------------
          }
)


setMethod(f="generate.sample.states",signature="markov.annualNonHomogeneous",definition=function(.Object, data)
{

  # Generate a synthtic series from the HMM states.
  # Note, the sampling uses the transition prob. matrix to estimate the long
  # term prob. of being in each state.
  nStates <- getNumStates(.Object)
  nSamples = nrow(data)
  states = rep(0,nSamples)

  # Get the initial state probabilites
  alpha = getInitialStateProbabilities(.Object)

  # Get the transition matrix.
  Tprob = getTransitionProbabilities(.Object)

  states[1] <- sample(1:nStates,1,prob=alpha)
  for (i in 2:nSamples) {
    alpha <- alpha %*% Tprob[,,i]
    states[i] = sample(1:nStates,1,prob=alpha)
  }

  return(states)
}
)


