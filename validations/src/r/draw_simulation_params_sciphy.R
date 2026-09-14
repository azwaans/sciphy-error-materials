## ---------------------------
##
## Script name: draw_simulation_params_sciphy
##
## Purpose of script: Draw simulation parameters from prior
## distribution for simulation study
##
## Author: Antoine Zwaans
##
## Date Created: 2025-10-16
##
## Copyright (c) Antoine Zwaans
## Email: zwaansantoine@gmail.com
##
library(extraDistr)
library(HDInterval)

## specify output location
outputDir = paste0("parameters/")

if (! dir.exists(outputDir)){
  dir.create(outputDir)
}


# Calculate origin time to end up with identical tree size distribution 
# but different total sampling proportion, given that the growth rate is 0.6
#exp(25*0.6)*0.0005 = exp(t*0.6)*0.05
#exp(t*0.6 - 25*0.6) = 0.0005/0.05
#0.6(t - 25) = log(0.0005/0.05)
#t - 25 = log(0.0005/0.05)/0.6
t0005 = 25
t05 = t0005 + log(0.0005/0.05)/0.6
#17.32472 
t0005 = 25
t1 = t0005 + log(0.0005)/0.6
#12.33183 

for (seed in 1:100){

  # 1. Substitution Model
  set.seed(seed)

  edit_probabilities = extraDistr::rdirichlet(n = 1, alpha = rep(1.5, 13))
  samples <- extraDistr::rdirichlet(n = 100000, alpha = rep(1.5, 13))
  samples <- extraDistr::rdirichlet(n = 100000, alpha = rep(1.5, 42))
  hdi_limits <- hdi(samples, credMass = 0.95)
  hdi_limits

  # Calculate the prior mean for the clock rates to match expected average mutations per tape
  # time * clock_rate = average_mutation
  clock_rate0005 = exp(-2)
  average_mutation = t0005*clock_rate0005
  clock_rate05 = average_mutation / t05
 
  # Calculate the prior mean for the miss rates to match expected loss per tape
  # expected loss = 1 - exp(-miss_rate* time)
  desired_loss_percent_heritable = 0.5
  miss_rate0005 = -log(1- desired_loss_percent_heritable)/t0005
  miss_rate05 = -log(1- desired_loss_percent_heritable)/t05
  
  clockRate = round(rlnorm(n = 1, meanlog = log(clock_rate0005), sdlog = 0.5), digits = 3)
  hdi_limits <- hdi(qlnorm, credMass = 0.95, meanlog = log(clock_rate0005), sdlog = 0.5)
  hdi_limits
  missingRate = round(rlnorm(n = 1, meanlog =  log(miss_rate0005), sdlog = 0.5), digits = 3)
  hdi_limits <- hdi(qlnorm, credMass = 0.95, meanlog =  log(miss_rate0005), sdlog = 0.5)
  hdi_limits
  missingProb = 0.0
 
  
  
  # 2. Write to file
  data_for_file = c(edit_probabilities, clockRate,missingRate,missingProb)
  names(data_for_file) = c( paste0("edit_prob_", 1:13), "clock_rate","missingRate","missingProb")

  write.csv(x = data_for_file, file = paste0(outputDir, "simParams_50_heritable_0.0005_", seed, ".csv"),

            quote = F, row.names = T,)
}


for (seed in 1:100){
  
  # 1. Substitution Model
  set.seed(seed)
  
  edit_probabilities = extraDistr::rdirichlet(n = 1, alpha = rep(1.5, 13))
  
  # Calculate the prior mean for the clock rates to match expected average mutations per tape
  # time * clock_rate = average_mutation
  clock_rate0005 = exp(-2)
  average_mutation = t0005*clock_rate0005
  clock_rate05 = average_mutation / t05
  
  # Calculate the prior mean for the miss rates to match expected loss per tape
  # expected loss = 1 - exp(-miss_rate* time)
  desired_loss_percent_heritable = 0.5
  miss_rate0005 = -log(1- desired_loss_percent_heritable)/t0005
  miss_rate05 = -log(1- desired_loss_percent_heritable)/t05
  
  clockRate = round(rlnorm(n = 1, meanlog = log(clock_rate05), sdlog = 0.5), digits = 3)
  hdi_limits <- hdi(qlnorm, credMass = 0.95, meanlog = log(clock_rate05), sdlog = 0.5)
  hdi_limits
  
  missingRate = round(rlnorm(n = 1, meanlog =  log(miss_rate05), sdlog = 0.5), digits = 3)
  hdi_limits <- hdi(qlnorm, credMass = 0.95, meanlog =  log(miss_rate05), sdlog = 0.5)
  hdi_limits
  
  missingProb = 0.0
  
  # 2. Write to file
  data_for_file = c(edit_probabilities, clockRate,missingRate,missingProb)
  names(data_for_file) = c( paste0("edit_prob_", 1:13), "clock_rate","missingRate","missingProb")
  
  write.csv(x = data_for_file, file = paste0(outputDir, "simParams_50_heritable_0.05_", seed, ".csv"),
            
            quote = F, row.names = T,)
}


for (seed in 1:100){
  
  # 1. Substitution Model
  set.seed(seed)
  
  edit_probabilities = extraDistr::rdirichlet(n = 1, alpha = rep(1.5, 13))
  
  # Calculate the prior mean for the clock rates to match expected average mutations per tape
  # time * clock_rate = average_mutation
  clock_rate0005 = exp(-2)
  average_mutation = t0005*clock_rate0005
  clock_rate1 = average_mutation / t1
  
  # Calculate the prior mean for the miss rates to match expected loss per tape
  # expected loss = 1 - exp(-miss_rate* time)
  desired_loss_percent_heritable = 0.5
  miss_rate0005 = -log(1- desired_loss_percent_heritable)/t0005
  miss_rate1 = -log(1- desired_loss_percent_heritable)/t1
  
  clockRate = round(rlnorm(n = 1, meanlog = log(clock_rate1), sdlog = 0.5), digits = 3)
  hdi_limits <- hdi(qlnorm, credMass = 0.95, meanlog = log(clock_rate1), sdlog = 0.5)
  hdi_limits
  
  
  missingRate = round(rlnorm(n = 1, meanlog =  log(miss_rate1), sdlog = 0.5), digits = 3)
  hdi_limits <- hdi(qlnorm, credMass = 0.95, meanlog = log(miss_rate1), sdlog = 0.5)
  hdi_limits
  
  
  missingProb = 0.0
  
  # 2. Write to file
  data_for_file = c(edit_probabilities, clockRate,missingRate,missingProb)
  names(data_for_file) = c( paste0("edit_prob_", 1:13), "clock_rate","missingRate","missingProb")
  
  write.csv(x = data_for_file, file = paste0(outputDir, "simParams_50_heritable_1_", seed, ".csv"),
            
            quote = F, row.names = T,)
}

for (seed in 1:100){
  
  # 1. Substitution Model
  set.seed(seed)
  
  edit_probabilities = extraDistr::rdirichlet(n = 1, alpha = rep(1.5, 13))
  
  # Calculate the prior mean for the clock rates to match expected average mutations per tape
  # time * clock_rate = average_mutation
  clock_rate0005 = exp(-2)
  average_mutation = t0005*clock_rate0005
  clock_rate05 = average_mutation / t05
  
  # Calculate the prior mean for the miss rates to match expected loss per tape
  # expected loss = 1 - exp(-miss_rate* time)
  desired_loss_percent_heritable = 0.25
  miss_rate0005 = -log(1- desired_loss_percent_heritable)/t0005
  miss_rate05 = -log(1- desired_loss_percent_heritable)/t05
  
  clockRate = round(rlnorm(n = 1, meanlog = log(clock_rate0005), sdlog = 0.5), digits = 3)
  
  missingRate = round(rlnorm(n = 1, meanlog =  log(miss_rate0005), sdlog = 0.5), digits = 3)
  hdi_limits <- hdi(qlnorm, credMass = 0.95, meanlog = log(miss_rate0005), sdlog = 0.5)
  hdi_limits
  
  
  desired_total_loss = 0.5
  #total_loss = desired_loss_percent_heritable + (1 - desired_loss_percent_heritable)*dropout_loss
  #(1 - desired_loss_percent_heritable)*dropout_loss = total_loss - desired_loss_percent_heritable
  dropout_loss = (desired_total_loss - desired_loss_percent_heritable) / (1 - desired_loss_percent_heritable)
  
  missingProb = round(rbeta(1,10*dropout_loss,10*(1-dropout_loss)), digits = 3)
  hdi_limits <- hdi(qbeta, credMass = 0.95,shape1= 10*dropout_loss,shape2=10*(1-dropout_loss))
  hdi_limits
  
  # 2. Write to file
  data_for_file = c(edit_probabilities, clockRate,missingRate,missingProb)
  names(data_for_file) = c( paste0("edit_prob_", 1:13), "clock_rate","missingRate","missingProb")
  
  write.csv(x = data_for_file, file = paste0(outputDir, "simParams_both_0.0005_", seed, ".csv"),
            
            quote = F, row.names = T,)
}

for (seed in 1:100){
  
  # 1. Substitution Model
  set.seed(seed)
  
  edit_probabilities = extraDistr::rdirichlet(n = 1, alpha = rep(1.5, 13))
  
  # Calculate the prior mean for the clock rates to match expected average mutations per tape
  # time * clock_rate = average_mutation
  clock_rate0005 = exp(-2)
  average_mutation = t0005*clock_rate0005
  clock_rate05 = average_mutation / t05
  
  # Calculate the prior mean for the miss rates to match expected loss per tape
  # expected loss = 1 - exp(-miss_rate* time)
  desired_loss_percent_heritable = 0.25
  miss_rate0005 = -log(1- desired_loss_percent_heritable)/t0005
  miss_rate05 = -log(1- desired_loss_percent_heritable)/t05
  
  clockRate = round(rlnorm(n = 1, meanlog = log(clock_rate05), sdlog = 0.5), digits = 3)
  
  missingRate = round(rlnorm(n = 1, meanlog =  log(miss_rate05), sdlog = 0.5), digits = 3)
  hdi_limits <- hdi(qlnorm, credMass = 0.95, meanlog = log(miss_rate0005), sdlog = 0.5)
  hdi_limits
  
  desired_total_loss = 0.5
  #total_loss = desired_loss_percent_heritable + (1 - desired_loss_percent_heritable)*dropout_loss
  #(1 - desired_loss_percent_heritable)*dropout_loss = total_loss - desired_loss_percent_heritable
  dropout_loss = (desired_total_loss - desired_loss_percent_heritable) / (1 - desired_loss_percent_heritable)
  
  missingProb = round(rbeta(1,10*dropout_loss,10*(1-dropout_loss)), digits = 3)
  hdi_limits <- hdi(qbeta, credMass = 0.95,shape1= 10*dropout_loss,shape2=10*(1-dropout_loss))
  hdi_limits
  
  # 2. Write to file
  data_for_file = c(edit_probabilities, clockRate,missingRate,missingProb)
  names(data_for_file) = c( paste0("edit_prob_", 1:13), "clock_rate","missingRate","missingProb")
  
  write.csv(x = data_for_file, file = paste0(outputDir, "simParams_both_0.05_", seed, ".csv"),
            
            quote = F, row.names = T,)
}


for (seed in 1:100){
  
  # 1. Substitution Model
  set.seed(seed)
  
  edit_probabilities = extraDistr::rdirichlet(n = 1, alpha = rep(1.5, 13))
  
  # Calculate the prior mean for the clock rates to match expected average mutations per tape
  # time * clock_rate = average_mutation
  clock_rate0005 = exp(-2)
  average_mutation = t0005*clock_rate0005
  clock_rate1 = average_mutation / t1
  
  # Calculate the prior mean for the miss rates to match expected loss per tape
  # expected loss = 1 - exp(-miss_rate* time)
  desired_loss_percent_heritable = 0.25
  miss_rate0005 = -log(1- desired_loss_percent_heritable)/t0005
  miss_rate1 = -log(1- desired_loss_percent_heritable)/t1
  
  clockRate = round(rlnorm(n = 1, meanlog = log(clock_rate1), sdlog = 0.5), digits = 3)
  
  missingRate = round(rlnorm(n = 1, meanlog =  log(miss_rate1), sdlog = 0.5), digits = 3)
  hdi_limits <- hdi(qlnorm, credMass = 0.95, meanlog = log(miss_rate1), sdlog = 0.5)
  hdi_limits
  
  
  desired_total_loss = 0.5
  #total_loss = desired_loss_percent_heritable + (1 - desired_loss_percent_heritable)*dropout_loss
  #(1 - desired_loss_percent_heritable)*dropout_loss = total_loss - desired_loss_percent_heritable
  dropout_loss = (desired_total_loss - desired_loss_percent_heritable) / (1 - desired_loss_percent_heritable)
  
  missingProb = round(rbeta(1,10*dropout_loss,10*(1-dropout_loss)), digits = 3)
  
  # 2. Write to file
  data_for_file = c(edit_probabilities, clockRate,missingRate,missingProb)
  names(data_for_file) = c( paste0("edit_prob_", 1:13), "clock_rate","missingRate","missingProb")
  
  write.csv(x = data_for_file, file = paste0(outputDir, "simParams_both_1_", seed, ".csv"),
            
            quote = F, row.names = T,)
}



for (seed in 1:100){
  
  # 1. Substitution Model
  set.seed(seed)
  
  edit_probabilities = extraDistr::rdirichlet(n = 1, alpha = rep(1.5, 13))
  
  # Calculate the prior mean for the clock rates to match expected average mutations per tape
  # time * clock_rate = average_mutation
  clock_rate0005 = exp(-2)
  average_mutation = t0005*clock_rate0005
  clock_rate05 = average_mutation / t05
  
  # Calculate the prior mean for the miss rates to match expected loss per tape
  # expected loss = 1 - exp(-miss_rate* time)
  desired_loss_percent_heritable = 0.0
  miss_rate0005 = -log(1- desired_loss_percent_heritable)/t0005
  miss_rate05 = -log(1- desired_loss_percent_heritable)/t05
  
  clockRate = round(rlnorm(n = 1, meanlog = log(clock_rate0005), sdlog = 0.5), digits = 3)
  
  missingRate = 0.0
  
  desired_total_loss = 0.5
  #total_loss = desired_loss_percent_heritable + (1 - desired_loss_percent_heritable)*dropout_loss
  #(1 - desired_loss_percent_heritable)*dropout_loss = total_loss - desired_loss_percent_heritable
  dropout_loss = (desired_total_loss - desired_loss_percent_heritable) / (1 - desired_loss_percent_heritable)
  
  missingProb = round(rbeta(1,10*dropout_loss,10*(1-dropout_loss)), digits = 3)
  hdi_limits <- hdi(qbeta, credMass = 0.95,shape1= 10*dropout_loss,shape2=10*(1-dropout_loss))
  hdi_limits
  
  
  # 2. Write to file
  data_for_file = c(edit_probabilities, clockRate,missingRate,missingProb)
  names(data_for_file) = c( paste0("edit_prob_", 1:13), "clock_rate","missingRate","missingProb")
  
  write.csv(x = data_for_file, file = paste0(outputDir, "simParams_50_dropout_0.0005_", seed, ".csv"),
            
            quote = F, row.names = T,)
}

for (seed in 1:100){
  
  # 1. Substitution Model
  set.seed(seed)
  
  edit_probabilities = extraDistr::rdirichlet(n = 1, alpha = rep(1.5, 13))
  
  # Calculate the prior mean for the clock rates to match expected average mutations per tape
  # time * clock_rate = average_mutation
  clock_rate0005 = exp(-2)
  average_mutation = t0005*clock_rate0005
  clock_rate05 = average_mutation / t05
  
  # Calculate the prior mean for the miss rates to match expected loss per tape
  # expected loss = 1 - exp(-miss_rate* time)
  desired_loss_percent_heritable = 0.0
  miss_rate0005 = -log(1- desired_loss_percent_heritable)/t0005
  miss_rate05 = -log(1- desired_loss_percent_heritable)/t05
  
  clockRate = round(rlnorm(n = 1, meanlog = log(clock_rate05), sdlog = 0.5), digits = 3)
  
  missingRate = 0.0
  
  desired_total_loss = 0.5
  #total_loss = desired_loss_percent_heritable + (1 - desired_loss_percent_heritable)*dropout_loss
  #(1 - desired_loss_percent_heritable)*dropout_loss = total_loss - desired_loss_percent_heritable
  dropout_loss = (desired_total_loss - desired_loss_percent_heritable) / (1 - desired_loss_percent_heritable)
  
  missingProb = round(rbeta(1,10*dropout_loss,10*(1-dropout_loss)), digits = 3)
  
  
  # 2. Write to file
  data_for_file = c(edit_probabilities, clockRate,missingRate,missingProb)
  names(data_for_file) = c( paste0("edit_prob_", 1:13), "clock_rate","missingRate","missingProb")
  
  write.csv(x = data_for_file, file = paste0(outputDir, "simParams_50_dropout_0.05_", seed, ".csv"),
            
            quote = F, row.names = T,)
}


for (seed in 1:100){
  
  # 1. Substitution Model
  set.seed(seed)
  
  edit_probabilities = extraDistr::rdirichlet(n = 1, alpha = rep(1.5, 13))
  
  # Calculate the prior mean for the clock rates to match expected average mutations per tape
  # time * clock_rate = average_mutation
  clock_rate0005 = exp(-2)
  average_mutation = t0005*clock_rate0005
  clock_rate1 = average_mutation / t1
  
  # Calculate the prior mean for the miss rates to match expected loss per tape
  # expected loss = 1 - exp(-miss_rate* time)
  desired_loss_percent_heritable = 0.0
  miss_rate0005 = -log(1- desired_loss_percent_heritable)/t0005
  miss_rate1 = -log(1- desired_loss_percent_heritable)/t1
  
  clockRate = round(rlnorm(n = 1, meanlog = log(clock_rate1), sdlog = 0.5), digits = 3)
  
  missingRate = 0.0
  
  desired_total_loss = 0.5
  #total_loss = desired_loss_percent_heritable + (1 - desired_loss_percent_heritable)*dropout_loss
  #(1 - desired_loss_percent_heritable)*dropout_loss = total_loss - desired_loss_percent_heritable
  dropout_loss = (desired_total_loss - desired_loss_percent_heritable) / (1 - desired_loss_percent_heritable)
  
  missingProb = round(rbeta(1,10*dropout_loss,10*(1-dropout_loss)), digits = 3)
  
  
  # 2. Write to file
  data_for_file = c(edit_probabilities, clockRate,missingRate,missingProb)
  names(data_for_file) = c( paste0("edit_prob_", 1:13), "clock_rate","missingRate","missingProb")
  
  write.csv(x = data_for_file, file = paste0(outputDir, "simParams_50_dropout_1_", seed, ".csv"),
            
            quote = F, row.names = T,)
}

