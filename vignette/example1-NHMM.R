# This HydroState examples builds and calibrates one two-state annual model and plots the results.
#----------------
rm(list=ls())
library(hydroState)
library(DEoptim)
library(truncnorm)

# Load flow data
#data(streamflow_annual)
# Load flow data
#data(streamflow_annual)
streamflow_annual<- read.csv(file = "G:/My Drive/Annual19502022/Annualrunoff2021-22.csv",TRUE,",")

# Extract one catchment
gaugeID = 238208
streamflow_annual = streamflow_annual[streamflow_annual$gauge==gaugeID,]
PET_ex<-streamflow_annual$pet

# Convert to format for hydroState
streamflow_annual = data.frame(year = streamflow_annual$hy_year, flow=streamflow_annual$q, precipitation=streamflow_annual$p)
filt = streamflow_annual$year>=1900
streamflow_annual = streamflow_annual[filt,]

# Build input objects.
# Note a 2-state linear model with first-order serial correlation is built. A truncated normal distribution is used for
# each Markov state.
transition.graph=matrix(TRUE,2,2)
Qhat = new('Qhat.log', input.data=streamflow_annual)
QhatModel = new('QhatModel.homo.gamma.linear.AR1', input.data=streamflow_annual, transition.graph=transition.graph)


PET_ex_adjusted <- PET_ex[(length(PET_ex) - nrow(streamflow_annual) + 1):length(PET_ex)]
markov_input <- cbind(year = streamflow_annual$year,
                      precipitation = streamflow_annual$precipitation,
                      PET = PET_ex_adjusted)
markov = new('markov.annualNonHomogeneous', transition.graph=transition.graph, inputForcing=markov_input,do.Logistic.Displacement=T)
#markov = new('markov.annualHomogeneous', transition.graph=transition.graph)
#markov@inputForcing <- as.matrix(markov_input)

# Build HydroState object
model = new('hydroState',input.data=streamflow_annual, Qhat.object=Qhat, QhatModel.object=QhatModel, markov.model.object=markov)

# Fit the model using DEoptim
model <- hydroState::fit(model,pop.size.perParameter = 10, max.generations=150)

# Name the states names with 1990 being defined as a 'norma' runoff year.
model <- setStateNames(model, 1990)

# Plot Viterbi states
viterbi(model)

# Plot pseduo residuals
#check.PseudoResiduals(model)

# # Build all seasonal models for this gauge, the claibration each model, select the best (byt AIC) and plot
# #-------------------------------------------
# all.Models <- new('hydroState.allModels',as.character(gaugeID), streamflow_annual, allow.flickering=F)
# all.Models <- fit(all.Models, pop.size.perParameter=5, max.generations=25, doParallel=F)
# best.model = getAIC.bestModel(all.Models)
