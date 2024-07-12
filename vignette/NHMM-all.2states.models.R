
rm(list=ls())
library(hydroState)
library(DEoptim)
library(truncnorm)

# Load flow data

#data(streamflow_annual)
streamflow_annual<- read.csv(file = "G:/My Drive/Annual19502022/Annualrunoff2021-22.csv",TRUE,",")

# Extract one catchment
gaugeID = 415206
streamflow_annual = streamflow_annual[streamflow_annual$gauge==gaugeID,]
PET_ex<-streamflow_annual$pet

# Convert to format for hydroState
streamflow_annual = data.frame(year = streamflow_annual$hy_year, flow=streamflow_annual$q, precipitation=streamflow_annual$p)
filt = streamflow_annual$year>=1900
streamflow_annual = streamflow_annual[filt,]


PET_ex_adjusted <- PET_ex[(length(PET_ex) - nrow(streamflow_annual) + 1):length(PET_ex)]
markov_input <- cbind(year = streamflow_annual$year,
                      precipitation = streamflow_annual$precipitation,
                      PET = PET_ex_adjusted)

# Build input objects.
# Note a 2-state linear model with first-order serial correlation is built. A truncated normal distribution is used for
# each Markov state.

# Build all combinations of annual models for this gauge.

all.Models <- new('hydroState.All2St.NHMM',as.character(gaugeID), streamflow_annual, allow.flickering=F)

# Calibrate each of the models.
# NOTE: comment out line 30 and uncomment line 31 to apply more robust calibration settings.
all.Models <- fit(all.Models, pop.size.perParameter=10, max.generations=1500, doParallel=F)
#all.Models <- fit(all.Models,pop.size.perParameter = 75,max.generations = 10000,reltol=1e-8,steptol=50, doParallel=F)


# Select the best model (byt AIC)
best.model = getAIC.bestModel(all.Models)

# Name the states names with 1990 being defined as a 'normal' runoff year.
best.model1 <- setStateNames(best.model[["model"]], c(1990, 1989, 1991, 1988, 1992))
best.model1<- setStateNames(all.Models@models[["model.2State.BC.AR1.S"]], c(1990, 1989, 1991, 1988, 1992))
# Plot Viterbi states
viterbi(best.model1)


# # Select the best model (byt AIC)
# all.Models = getAIC.bestModel(all.Models)
#
# # Name the states names with 1990 being defined as a 'normal' runoff year.
# # The years after 1990 are the options used if there is no flow data for 1990.
# all.Models <- setStateNames(all.Models, c(1990, 1989, 1991, 1988, 1992))
#
# # Plot Viterbi states
# viterbi(all.Models)
#
# # Plot pseduo residuals
# check.PseudoResiduals(all.Models)
