## ---------------------------
##
## Script name: get summary data from all log files (based on S. Seidel's scripts)
##
## -----------------------------

library(tracerer)
library(HDInterval)
library(ggplot2)
library(beastio)
library(stringr)
library(tidyverse)
library(phytools)
library(treebalance)

## ---------------------------
log_dir = "log/"

simulation_dir = "parameters/"

parameters_of_interest = c("clockRate", paste0("editRate.", 1:13),"lossRate","dropoutProb")


inc = function(x){
  eval.parent(substitute(x <- x+1))
}

is_parameter_in_hpd = function(hpd_lower, hpd_upper, true_parameter){
  if(true_parameter >= hpd_lower && true_parameter <= hpd_upper)
    return(TRUE)
  else
    return(FALSE)
}

num_tips <- c()
tree_heights <- c()
tree_lengths <- c()
tree_B1 <- c()
for(i in 1:100) {
  file_name <- paste0("data/simulate_with_both_sampling=0.00003.",i,".newick")
  tree <- readLines(file_name)
  tree <- str_remove_all(tree,"\\[.........\\]")
  tree <- str_remove_all(tree,"\\[..........\\]")
  tree <- str_remove_all(tree,"\\[...........\\]")
  tree <- paste0(tree,";")
  p_tree <- ape::read.tree(text = tree)
  tree_heights <- c(tree_heights,max(nodeHeights(p_tree)))
  tree_lengths <- c(tree_lengths,sum(p_tree$edge.length))
  tree_B1 <- c(tree_B1,treebalance::B1I(as.phylo(p_tree)))
  num_tips <- c(num_tips,length(p_tree$tip.label))
  
}

true_tree_stats <- data.frame(seed=1:100,ntips=num_tips,true_heights=tree_heights,true_lengths=tree_lengths,true_b1=tree_B1)


nr_converged_chains = 0

clock_rate_array_5_inference = data.frame(seed=1:100, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F)

insert_rate_inference = data.frame(seed=rep(1:100, each = 13), hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F,
                                   insertRate= rep(1:13,100))

miss_rate_array_5_inference = data.frame(seed=1:100, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F)

miss_prob_array_5_inference = data.frame(seed=1:100, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F)

birth_rate_array_5_inference = data.frame(seed=1:100, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F)

death_rate_array_5_inference = data.frame(seed=1:100, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F)


growth_rate_array_5_inference = data.frame(seed=1:100, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F)

tree_height_inference = data.frame(seed=1:100, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F)

tree_length_inference = data.frame(seed=1:100, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F)


for (seed in  1:100){

  
  print(seed)
  # get inference log
  log_file = paste0("infer_with_error_model_both_sampling=0.00003.", seed, ".log")
  if(file.exists(paste0(log_dir, log_file))) {
    

  
  
  log_data = parse_beast_tracelog_file(paste0(log_dir, log_file))
  log_data_wo_burnin = remove_burn_ins(log_data, burn_in_fraction = 0.3)

  # get true parameter from simulation
  simulation_parameter_file = paste0(simulation_dir, "simParams_", seed, ".csv")
  simulation_parameters = read.csv(simulation_parameter_file)
  colnames(simulation_parameters) = c("parameter", "value")

  #check that chains has converged
  esses = calc_esses(log_data_wo_burnin, sample_interval = 1000)
  esses = esses[!colnames(esses) %in% c("prior")]
  esses = esses[! names(esses) %in% c("editRate.13","rho_BDSKY_Contemp.t.input_alignment")]
  esses

  if (any(esses < 200)){
    problematic_ess = esses[which(esses < 200)]
    print(paste("Ess for log file ", log_file, " has ESS < 200: for parameters", paste(names(problematic_ess), collapse = ",")))
    next

   }else{
    inc(nr_converged_chains)
   }

  clock_rate_array_5_hpd = hdi(object = log_data_wo_burnin$clockRate, credMass = 0.95)

  true_clock_rate = simulation_parameters[[14]]

  ## check for clock rate for alignments with array length 5
    recovered = is_parameter_in_hpd(hpd_lower = clock_rate_array_5_hpd["lower"], hpd_upper = clock_rate_array_5_hpd["upper"],
                      true_parameter = true_clock_rate)
  clock_rate_array_5_inference[seed, ] = c(seed, clock_rate_array_5_hpd, median(log_data_wo_burnin$clockRate), true_clock_rate, recovered)


  ### for the missRate
  miss_rate_array_5_hpd = hdi(object = log_data_wo_burnin$lossRate, credMass = 0.95)
  
  ### extract simulation value
  true_miss_rate = simulation_parameters[[16]]
  
  ## check for clock rate for alignments with array length 5
  recovered = is_parameter_in_hpd(hpd_lower = miss_rate_array_5_hpd["lower"], hpd_upper = miss_rate_array_5_hpd["upper"],
                                  true_parameter = true_miss_rate)
  miss_rate_array_5_inference[seed, ] = c(seed, miss_rate_array_5_hpd, median(log_data_wo_burnin$lossRate), true_miss_rate, recovered)
  
  
  ### for the missProb
  miss_prob_array_5_hpd = hdi(object = log_data_wo_burnin$dropoutProb, credMass = 0.95)
  
  ###
  true_miss_prob = simulation_parameters[[15]]
  print(true_miss_prob)
  
  ## check for clock rate for alignments with array length 5
  recovered = is_parameter_in_hpd(hpd_lower = miss_prob_array_5_hpd["lower"], hpd_upper = miss_prob_array_5_hpd["upper"],
                                  true_parameter = true_miss_prob)
  print(recovered)
  miss_prob_array_5_inference[seed, ] = c(seed, miss_prob_array_5_hpd, median(log_data_wo_burnin$dropoutProb), true_miss_prob, recovered)
  
  

  
  ####
  # TREE HEIGHT
  ###
  
  ### extrac tree height hpd
  tree_height_hpd = hdi(object = log_data_wo_burnin$'treeHeight.t.h3n2_2deme', credMass = 0.95)
  
  ### extract simulation value
  true_tree_height = true_tree_stats[true_tree_stats$seed == seed, "true_heights"]
  
  ## check recovery
  recovered = is_parameter_in_hpd(hpd_lower = tree_height_hpd["lower"], hpd_upper = tree_height_hpd["upper"],
                                  true_parameter = true_tree_height)
  
  ### add all stats to table
  tree_height_inference[seed, ] = c(seed, tree_height_hpd, median(log_data_wo_burnin$'treeHeight.t.h3n2_2deme'), true_tree_height, recovered)
  
  ####
  # TREE LENGTH
  ###
  
  ## extract tree length hpd
  tree_length_hpd = hdi(object = log_data_wo_burnin$'treeLength.t.h3n2_2deme', credMass = 0.95)
  
  ### extract simulation value
  true_tree_length = true_tree_stats[true_tree_stats$seed == seed, "true_lengths"]
  
  ## check recovery
  recovered = is_parameter_in_hpd(hpd_lower = tree_length_hpd["lower"], hpd_upper = tree_length_hpd["upper"],
                                  true_parameter = true_tree_length)
  
  ### add all stats to table
  tree_length_inference[seed, ] = c(seed, tree_length_hpd, median(log_data_wo_burnin$'treeLength.t.h3n2_2deme'), true_tree_length, recovered)
  
  
  
  for (insert_rate_nr in 1:13){

    insert_rate = paste0("editRate.", insert_rate_nr)
    insert_rate_hpd = hdi(object = log_data_wo_burnin[, insert_rate] * log_data_wo_burnin$clockRate, credMass = 0.95)
    median_insert_rate = median(log_data_wo_burnin[, insert_rate] * log_data_wo_burnin$clockRate)
    true_insert_rate = simulation_parameters[[insert_rate_nr]] * simulation_parameters[[14]]
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
clock_rate_array_5_inference$orderedSeed = 1:100

### miss rate of alignment with length 5 targets
miss_rate_array_5_inference = miss_rate_array_5_inference[order(miss_rate_array_5_inference$true_value), ]
miss_rate_array_5_inference$orderedSeed = 1:100

### miss prob of alignment with length 5 targets
miss_prob_array_5_inference = miss_prob_array_5_inference[order(miss_prob_array_5_inference$true_value), ]
miss_prob_array_5_inference$orderedSeed = 1:100


birth_rate_array_5_inference = birth_rate_array_5_inference[order(birth_rate_array_5_inference$true_value), ]
birth_rate_array_5_inference$orderedSeed = 1:100


death_rate_array_5_inference = death_rate_array_5_inference[order(death_rate_array_5_inference$true_value), ]
death_rate_array_5_inference$orderedSeed = 1:100

growth_rate_array_5_inference = growth_rate_array_5_inference[order(growth_rate_array_5_inference$true_value), ]
growth_rate_array_5_inference$orderedSeed = 1:100


### tree height
tree_height_inference = tree_height_inference[order(tree_height_inference$true_value), ]
tree_height_inference$orderedSeed = 1:100

### miss prob of alignment with length 5 targets
tree_length_inference = tree_length_inference[order(tree_length_inference$true_value), ]
tree_length_inference$orderedSeed = 1:100


#coverages per insert rate
coverages_per_insert <- c()
for(i in 1:13) {
  coverage <- sum(insert_rate_inference[which(insert_rate_inference$insertRate == i),"recovered"])/nr_converged_chains
  coverages_per_insert <- c(coverages_per_insert,coverage)
}
coverages_per_insert

#coverage for the clock rate
sum(clock_rate_array_5_inference$recovered)/nr_converged_chains
sum(miss_rate_array_5_inference$recovered)/nr_converged_chains
sum(miss_prob_array_5_inference$recovered)/nr_converged_chains

sum(birth_rate_array_5_inference$recovered)/nr_converged_chains
sum(death_rate_array_5_inference$recovered)/nr_converged_chains

#coverage for the tree height
sum(tree_height_inference$recovered)/nr_converged_chains
#coverage for the tree length
sum(tree_length_inference$recovered)/nr_converged_chains

write.csv(x = insert_rate_inference, "scarring_rate_product_inference.csv",quote = F, row.names = F)
write.csv(x = clock_rate_array_5_inference, "clock_rate_inference.csv",quote = F, row.names = F)
write.csv(x = miss_rate_array_5_inference, "miss_rate_inference.csv",quote = F, row.names = F)
write.csv(x = miss_prob_array_5_inference, "miss_prob_inference.csv",quote = F, row.names = F)

write.csv(x = birth_rate_array_5_inference, "birth_rate_inference.csv",quote = F, row.names = F)
write.csv(x = death_rate_array_5_inference, "death_rate_inference.csv",quote = F, row.names = F)
write.csv(x = growth_rate_array_5_inference, "growth_rate_inference.csv",quote = F, row.names = F)

write.csv(x = tree_height_inference, "tree_height_inference.csv",quote = F, row.names = F)
write.csv(x = tree_length_inference, "tree_length_inference.csv",quote = F, row.names = F)



