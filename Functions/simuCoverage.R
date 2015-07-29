## ----------------------------------------------------------------------------
## Name : simuCoverage.R
## ----------------------------------------------------------------------------
## Authors: Cathy Y. Y. Lee 
## ----------------------------------------------------------------------------
## Last updated: 29 JULY 2015
## ----------------------------------------------------------------------------
## Description: Obtain coverage percentage summaries of simulation results
## ----------------------------------------------------------------------------

## Read in simulation results:

simRes <- read.table(filenameCov,header=TRUE)

## Obtain the coverage percentages:

simRes <- as.data.frame(simRes)

covPercs <- round(100*colSums(simRes[,-1])/nrow(simRes))

print(covPercs)

############ End of simuCoverage.R ############

