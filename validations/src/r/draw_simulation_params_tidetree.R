## ---------------------------
##
## Script name: draw_simulation_params_tidetree
##
## Purpose of script: Draw simulation parameters from prior
## distribution for simulation study
##
## Author: Antoine Zwaans (based on S.Seidel's original simulation study)
##
## Date Created: 2025-10-16
##
## Copyright (c) Antoine Zwaans
## Email: zwaansantoine@gmail.com
##

## specify output location
outputDir = "./parameters/"

if (! dir.exists(outputDir)){
  dir.create(outputDir)
}


for (seed in 1:100){

# Substitution Model
#raw scarring rates from Exponential distribution
nScarringOutcomes = 13

set.seed(seed)
scarringRates = round(rexp(n = nScarringOutcomes, rate = 2), digits = 2)
scarringRates[nScarringOutcomes] = 1 #set last scarringRate to 1. This will be fixed during the inference as well.
hdi_limits <- hdi(qexp, credMass = 0.95, rate=2)
hdi_limits


clockRate = round(rlnorm(n = 1, meanlog = -4, sdlog = 0.5), digits = 5)
hdi_limits <- hdi(qlnorm, credMass = 0.95, meanlog = -4, sdlog = 0.5)
hdi_limits

lossRate = round(rlnorm(n = 1, meanlog = -1, sdlog = 1.0), digits = 5)
hdi_limits <- hdi(qlnorm, credMass = 0.95, meanlog = -1, sdlog = 1.0)
hdi_limits

dropoutProb = round(rbeta(1,3,7), digits = 5)

data.frame(matrix(nrow = 1, ncol = nScarringOutcomes +3, data = 0))
dataForFile = data.frame(matrix(nrow = 1, ncol = nScarringOutcomes +3, data = 0))
colnames(dataForFile) = c(paste("scarringRate_", 1:nScarringOutcomes), "clockRate", "dropoutProb", "lossRate")

dataForFile[1, ] = c(scarringRates, clockRate, dropoutProb, lossRate)
write.csv(x = dataForFile, file = paste0(outputDir, "simParams_", seed, ".csv"), quote = F, row.names = F)
}

