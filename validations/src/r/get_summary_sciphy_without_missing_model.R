#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
}


analysis_name <- args[1]
sampling_proportion  <- args[2]
targets <- args[3]

library(tracerer)
library(HDInterval)
library(ggplot2)

## ---------------------------
log_dir = "log/"

simulation_dir = "parameters/"

parameters_of_interest = c("clockRate", paste0("insertRate.", 1:13))


inc = function(x){
  eval.parent(substitute(x <- x+1))
}

is_parameter_in_hpd = function(hpd_lower, hpd_upper, true_parameter){
  if(true_parameter >= hpd_lower && true_parameter <= hpd_upper)
    return(TRUE)
  else
    return(FALSE)
}

nr_converged_chains = 0

clock_rate_array_5_inference = data.frame(seed=1:50, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F)

insert_rate_inference = data.frame(seed=rep(1:50, each = 13), hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F,
                                   insertRate= rep(1:13,100))

birth_rate_array_5_inference = data.frame(seed=1:50, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F)

death_rate_array_5_inference = data.frame(seed=1:50, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F)


growth_rate_array_5_inference = data.frame(seed=1:50, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F)


for (seed in  1:50){

  print(seed)
  # get inference log
  
  log_file = paste0("infer_filtered_simulated_",analysis_name,"_",sampling_proportion,"_",targets,"targets.seed=", seed, ".log")
  if(file.exists(paste0(log_dir, log_file))) {
    
    log_data = parse_beast_tracelog_file(paste0(log_dir, log_file))
    log_data_wo_burnin = remove_burn_ins(log_data, burn_in_fraction = 0.3)
    
    # get true parameter from simulation
    simulation_parameter_file = paste0(simulation_dir, "simParams_",analysis_name,"_",sampling_proportion,"_",seed, ".csv")
    simulation_parameters = read.csv(simulation_parameter_file)
    colnames(simulation_parameters) = c("parameter", "value")
    
    #check that chains has converged
    esses = calc_esses(log_data_wo_burnin, sample_interval = 1000)
    esses = esses[!colnames(esses) %in% c("prior")]
    esses = esses[! names(esses) %in% c("samplingProportion", "treeHeight.t.alignment", "treeLength.t.alignment",paste0("treeLikelihood_",1:20))]
    esses
    
    if (any(esses < 200)){
      problematic_ess = esses[which(esses < 200)]
      print(paste("Ess for log file ", log_file, " has ESS < 200: for parameters", paste(names(problematic_ess), collapse = ",")))
      next

    }else{
    inc(nr_converged_chains)
    }
    
    clock_rate_array_5_hpd = hdi(object = log_data_wo_burnin$clockRate, credMass = 0.95)
    
    true_clock_rate = simulation_parameters[14, 2]
    
    ## check for clock rate for alignments with array length 5
    recovered = is_parameter_in_hpd(hpd_lower = clock_rate_array_5_hpd["lower"], hpd_upper = clock_rate_array_5_hpd["upper"],
                                    true_parameter = true_clock_rate)
    clock_rate_array_5_inference[seed, ] = c(seed, clock_rate_array_5_hpd, median(log_data_wo_burnin$clockRate), true_clock_rate, recovered)
    
    
    clock_rate_array_5_hpd = hdi(object = log_data_wo_burnin$clockRate, credMass = 0.95)
    
    for (insert_rate_nr in 1:13){
      
      insert_rate = paste0("insertRates.", insert_rate_nr)
      insert_rate_hpd = hdi(object = log_data_wo_burnin[, insert_rate], credMass = 0.95)
      median_insert_rate = median(log_data_wo_burnin[, insert_rate])
      true_insert_rate = simulation_parameters[insert_rate_nr, "value"]
      recovered = is_parameter_in_hpd(hpd_lower = insert_rate_hpd["lower"], hpd_upper = insert_rate_hpd["upper"],
                                      true_parameter = true_insert_rate)
      
      
      row_index = which(insert_rate_inference$seed == seed & insert_rate_inference$insertRate == insert_rate_nr)
      insert_rate_inference[row_index, ] = c(seed, as.numeric(insert_rate_hpd), median_insert_rate, true_insert_rate,
                                             recovered, insert_rate_nr)
    }
    
    
    ### for the birthrate
    birth_rate_array_5_hpd = hdi(object = log_data_wo_burnin$birthRate, credMass = 0.95)
    
    ###
    true_birth_rate = 0.8
    print(true_birth_rate)
    
    ## check for clock rate for alignments with array length 5
    recovered = is_parameter_in_hpd(hpd_lower = birth_rate_array_5_hpd["lower"], hpd_upper = birth_rate_array_5_hpd["upper"],
                                    true_parameter = true_birth_rate)
    print(recovered)
    birth_rate_array_5_inference[seed, ] = c(seed, birth_rate_array_5_hpd, median(log_data_wo_burnin$birthRate), true_birth_rate, recovered)
    
    
    ### for the birthrate
    death_rate_array_5_hpd = hdi(object = log_data_wo_burnin$deathRate, credMass = 0.95)
    
    ###
    true_death_rate = 0.2
    print(true_death_rate)
    
    ## check for clock rate for alignments with array length 5
    recovered = is_parameter_in_hpd(hpd_lower = death_rate_array_5_hpd["lower"], hpd_upper = death_rate_array_5_hpd["upper"],
                                    true_parameter = true_death_rate)
    print(recovered)
    death_rate_array_5_inference[seed, ] = c(seed, death_rate_array_5_hpd, median(log_data_wo_burnin$deathRate), true_death_rate, recovered)
    
    log_data_wo_burnin$growthRate = log_data_wo_burnin$birthRate - log_data_wo_burnin$deathRate
    
    growth_rate_array_5_hpd = hdi(object = log_data_wo_burnin$growthRate, credMass = 0.95)
    
    
    true_growth_rate = 0.6
    print(true_growth_rate)
    
    ## check for clock rate for alignments with array length 5
    recovered = is_parameter_in_hpd(hpd_lower = growth_rate_array_5_hpd["lower"], hpd_upper = growth_rate_array_5_hpd["upper"],
                                    true_parameter = true_growth_rate)
    print(recovered)
    growth_rate_array_5_inference[seed, ] = c(seed, growth_rate_array_5_hpd, median(log_data_wo_burnin$growthRate), true_growth_rate, recovered) 
    
    
    
  }
  
  
  
  
}

### clock rate of alignment with length 5 targets
clock_rate_array_5_inference = clock_rate_array_5_inference[order(clock_rate_array_5_inference$true_value), ]
clock_rate_array_5_inference$orderedSeed = 1:50


birth_rate_array_5_inference = birth_rate_array_5_inference[order(birth_rate_array_5_inference$true_value), ]
birth_rate_array_5_inference$orderedSeed = 1:50


death_rate_array_5_inference = death_rate_array_5_inference[order(death_rate_array_5_inference$true_value), ]
death_rate_array_5_inference$orderedSeed = 1:50

growth_rate_array_5_inference = growth_rate_array_5_inference[order(growth_rate_array_5_inference$true_value), ]
growth_rate_array_5_inference$orderedSeed = 1:50



#coverages per insert rate
coverages_per_insert <- c()
for(i in 1:13) {
  coverage <- sum(insert_rate_inference[which(insert_rate_inference$insertRate == i),"recovered"])/nr_converged_chains
  coverages_per_insert <- c(coverages_per_insert,coverage)
}
coverages_per_insert

#coverage for the clock rate
sum(clock_rate_array_5_inference$recovered)/nr_converged_chains
sum(birth_rate_array_5_inference$recovered)/nr_converged_chains
sum(death_rate_array_5_inference$recovered)/nr_converged_chains


write.csv(x = clock_rate_array_5_inference, paste0(analysis_name,"_",sampling_proportion,"_clock_rate_inference_no_error.csv"),quote = F, row.names = F)
write.csv(x = birth_rate_array_5_inference, paste0(analysis_name,"_",sampling_proportion,"_","birth_rate_inference_no_error.csv"),quote = F, row.names = F)
write.csv(x = death_rate_array_5_inference, paste0(analysis_name,"_",sampling_proportion,"_","death_rate_inference_no_error.csv"),quote = F, row.names = F)
write.csv(x = growth_rate_array_5_inference, paste0(analysis_name,"_",sampling_proportion,"_","growth_rate_inference_no_error.csv"),quote = F, row.names = F)


