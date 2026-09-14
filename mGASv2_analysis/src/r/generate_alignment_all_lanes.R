library(stringr)
library(dplyr)
library(ggplot2)

pick_code <- function(edit,code_map) {
  code <- code_map[which(code_map$insert == edit), "integer"]
  #print(code)
  return(code)
  
}

cell_has_selected_barcodes = function(cell, dat, selected_barcodes){
  
  cell_barcodes = dat[dat$Cell == cell, "TargetBC"]
  
  if (all(selected_barcodes %in% cell_barcodes)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

cell_has_at_least_one_of_the_selected_barcodes = function(cell, dat, selected_barcodes, at_least){
  
  cell_barcodes = dat[dat$Cell == cell, "TargetBC"]
  
  if (sum(selected_barcodes %in% cell_barcodes) >= at_least){
    return(TRUE)
  }else{
    return(FALSE)
  }
}


create_nexus_alignments <- function(edit_table,analysis_name,tbcs_length_1,tbcs_length_2,tbcs_length_3,tbcs_length_4){
  
  
  list_of_cells <- unique(edit_table$Cell)
  taxa <- length(list_of_cells)
  targets <- unique(edit_table$TargetBC)
  loss <- c()
  for(i in 1:length(targets)){
    
    target = targets[i]
    #dir.create(file.path(paste0("data/",analysis_name,"alignment.targetBC=",target,".sciphy")), recursive = TRUE)
    sink(paste0("runs_on_all_lanes_final/data/",analysis_name,".alignment.targetBC=",i,".sciphy"))
    cat("#NEXUS")
    cat("\n")
    cat("\n")
    cat("begin taxa;")
    cat("\n")
    cat("\t")
    cat(paste0("dimensions ntax=",taxa,";"))
    cat("\n")
    cat("\t")
    cat(paste0("taxlabels ",paste(list_of_cells,collapse=" "),";"))
    cat("\n")
    cat("end;")
    cat("\n")
    cat("\n")
    cat("begin characters;")
    cat("\n")
    cat("\t")
    if(target %in% tbcs_length_1) {
      cat("dimensions nchar=1;")
    }
    if(target %in% tbcs_length_2) {
      cat("dimensions nchar=2;")
    }
    if(target %in% tbcs_length_3) {
      cat("dimensions nchar=3;")
    }
    else if(target  %in% tbcs_length_4) {
      cat("dimensions nchar=4;")
    }
    else {
      cat("dimensions nchar=5;")
    }
    
    cat("\n")
    cat("\t")
    cat("format datatype=integer;")
    cat("\n")
    cat("\t")
    cat("matrix")
    counter <- 0
    for(cellID in list_of_cells) {
      cat("\n")
      cat("\t")
      cat("\t")
      cat(cellID)
      cat(" ")
      #concatenate all sites into one sequence for particular cellID: 
      
      if(paste((edit_table[(edit_table$Cell == cellID) & (edit_table$TargetBC == target),3:7]),collapse = ",") == "character(0),character(0),character(0),character(0),character(0)") {
        if(target %in% tbcs_length_1) {
          cat("?")
        }
        else if(target %in% tbcs_length_2) {
          cat("?,?")
        }
        else if(target %in% tbcs_length_3) {
          cat("?,?,?")
        }
        else if(target  %in% tbcs_length_4) {
          cat("?,?,?,?")
        }
        else {
          cat("?,?,?,?,?")
        }
        
        counter <- counter +1 
      }
      else {
        #target length 2
        if(target %in% tbcs_length_1) {
          cat(paste((edit_table[(edit_table$Cell == cellID) & (edit_table$TargetBC == target),3:3]),collapse = ","))
        }
        else if(target %in% tbcs_length_2) {
          cat(paste((edit_table[(edit_table$Cell == cellID) & (edit_table$TargetBC == target),3:4]),collapse = ","))
        }
        else if(target %in% tbcs_length_3) {
          cat(paste((edit_table[(edit_table$Cell == cellID) & (edit_table$TargetBC == target),3:5]),collapse = ","))
        }
        else if(target  %in% tbcs_length_4) {
          cat(paste((edit_table[(edit_table$Cell == cellID) & (edit_table$TargetBC == target),3:6]),collapse = ","))
        }
        else {
          cat(paste((edit_table[(edit_table$Cell == cellID) & (edit_table$TargetBC == target),3:7]),collapse = ","))
        }
      }
      
      
    }
    
    loss <- c(loss,counter)
    cat(";")
    cat("\n")
    cat("end;")
    
    
    sink()
  }
  
  return(data.frame(losses=loss,targetBCs=targets))
}

# get inserts across all cells, arrange in descrending freq,
# map each insert to integer;
# assumes that the inserts are already tri-nucleotides
get_insert_to_integer_map = function(edit_table, columns){
  
  insert_to_integer_map <- data.frame(table(unlist(edit_table[, columns]))) # %>% arrange(desc(Freq))

  insert_to_integer_map = insert_to_integer_map %>% arrange(desc(Freq))
  
  colnames(insert_to_integer_map)[1] = "insert"
  
  index_unedited <- which(insert_to_integer_map$insert == "NA")
  
  without_unedited <- insert_to_integer_map[-index_unedited,]
  
  corrected_order_insert_to_integer_map <- rbind(insert_to_integer_map[index_unedited,],insert_to_integer_map[-index_unedited,])
  
  print("where are we")
  corrected_order_insert_to_integer_map$integer <- as.character(seq(0,dim(corrected_order_insert_to_integer_map)[1] -1))
  
  
  return(corrected_order_insert_to_integer_map)
  
}

#lane_number <- 2
#clustering_k <- 8

#setwd("/Users/azwaans/Documents/Projects/TYPEWRITER/missing_data/final_mGASv2_missing")


import_and_clean_table <- function(lane_number) { 
  
  
  input_file = paste0("unfiltered_data/mGASv2_Lane",lane_number,"_CellByTape_10X_bamExtractV2_t3_collapse.csv")
  cat("\n",'Reading-in file :', input_file)
  edit_table = read.csv(paste0("unfiltered_data/","mGASv2_Lane",lane_number,"_CellByTape_10X_bamExtractV2_t3_collapse.csv"), stringsAsFactors = F,
                        header = T, na.strings=c("","NA"))
  
  cell_filters = read.csv(paste0("unfiltered_data/mGASv2_Lane",lane_number,"_filtered_feature_bc_matrix_barcodes.tsv"), header = F)
  
  edit_table <- filter(edit_table, Cell %in% cell_filters$V1)
  
  cat("\n",'Total number of cells in unfiltered data')
  cat("\n",length(unique(edit_table$Cell))) 
  
  cat("\n",'Filtering step 1 - remove NA TBcs and cBcs')
  ### filter out entries that have no TargetBC or UMI
  edit_table = edit_table[which(!(is.na(edit_table$TargetBC))), ]
  edit_table = edit_table[which(!(is.na(edit_table$Cell))), ]
  cat("\n",'Total number of cells after filtering')
  cat("\n",length(unique(edit_table$Cell))) 
  
  
  cat("\n",'Filtering step 2 - remove barcodes that violate ordering')
  ### filter out entries that have not correctly ordered edits
  #### pairwise comparisons, remove all entries without edit at site 1 and edit at site 2-5
  edit_table = edit_table[! (is.na(edit_table[, 3]) &
                               !(is.na(edit_table[, 4]))), ]
  length(unique(edit_table$Cell)) 
  edit_table = edit_table[! (is.na(edit_table[, 3]) &
                               !(is.na(edit_table[, 5]))), ]
  length(unique(edit_table$Cell)) 
  edit_table = edit_table[! (is.na(edit_table[, 3]) &
                               !(is.na(edit_table[, 6]))), ]
  length(unique(edit_table$Cell)) 
  edit_table = edit_table[! (is.na(edit_table[, 3]) &
                               !(is.na(edit_table[, 7]))), ]
  
  
  #### same for site 2
  edit_table = edit_table[! (is.na(edit_table[, 4]) &
                               !(is.na(edit_table[, 5]))), ]
  edit_table = edit_table[! (is.na(edit_table[, 4]) &
                               !(is.na(edit_table[, 6]))), ]
  edit_table = edit_table[! (is.na(edit_table[, 4]) &
                               !(is.na(edit_table[, 7]))), ]
  
  #### same for site 3
  edit_table = edit_table[! (is.na(edit_table[, 5]) &
                               !(is.na(edit_table[, 6]))), ]
  edit_table = edit_table[! (is.na(edit_table[, 5]) &
                               !(is.na(edit_table[, 7]))), ]
  
  #### same for site 4
  edit_table = edit_table[! (is.na(edit_table[, 6]) &
                               !(is.na(edit_table[, 7]))), ]
  
  cat("\n",'Total number of cells after filtering')
  cat("\n",length(unique(edit_table$Cell))) 
  
  #The filtering was likely done previously, as n_cells is still= 58554
  frequency_of_target_bcs = as.data.frame(table(edit_table$TargetBC, useNA = "ifany"))
  frequency_of_target_bcs = frequency_of_target_bcs[order(frequency_of_target_bcs$Freq, decreasing = T), ]

  
  return(edit_table)
}


convert_table_to_integers <- function(original_csv,lane_number,clustering_level,cluster_number) {
  
   site_columns = 3:7

  # # 1st Convert 6-mers to 3-mers
  # # by removing the last 3 nucleotides that do not contain information
  # # (they are always GGA); replace "None by "NA"
  for(site_column in site_columns) {

    original_csv[original_csv[, site_column] == "None", site_column] <- "NA"
    original_csv[,site_column] <- substring(original_csv[, site_column], 1,3)
    

  }

  cat("\n",'Create trinucleotide to integer mapping by frequency (most frequent = 0, second most frequent = 1 ... )')

  # # 2nd, create a map from the trinucleotides to an integer
  insert_to_integer_map = get_insert_to_integer_map( original_csv, columns=site_columns)
  
  if(insert_to_integer_map$insert[1] != "NA") {
    stop("the most frequent edit isn't non-edited")
  }
  #insert_to_integer_map_full = get_insert_to_integer_map(original_csv, columns=site_columns, "cell_culture")

  write.csv(insert_to_integer_map, file = paste0("runs_on_all_lanes_final/conversion/insert_to_integer_map_","mGASv2_Lane",lane_number,"_Clustering",clustering_level,"_Cluster",cluster_number,".csv"))
  cat("\n",'Mapping saved as :', paste0("insert_to_integer_map","mGASv2_Lane",lane_number,"_Clustering",clustering_level,"_Cluster",cluster_number,".csv"))


  cat("\n",'Convert the whole table with mapping')
  # convert these trinucleotides to integer in the edit table
  for(site_column in site_columns) {
    original_csv[, site_column] <- unlist(lapply(original_csv[, site_column], function(x) {pick_code(x, insert_to_integer_map)}))

  }

  return(list(table=original_csv,dim=nrow(insert_to_integer_map)))
  
  
}







create_stats_per_lane_per_cluster <- function(lane_number, clustering_level) {
  
  
  edit_table <- import_and_clean_table(lane_number)
  
  # cat("\n",'Converting hexanucleotides to trinucleotides')
  
   original_csv <- edit_table
  #
  
  cat("\n",'Import tbcs filters obtained with mGASv2_TargetBC_filtering_v2.R, with entropy')
  #select sets of cells with a set of TBCs identified to be from the same group of celss
  
  for(cluster_number in 1:clustering_level) { 
  selected_barcode_file = paste0("runs_on_all_lanes_final/clustering/filtered_tbcs_Lane",lane_number,"_clustering",clustering_level,"_group",cluster_number,".csv")
  selected_barcodes_dat = read.csv(selected_barcode_file)
  selected_barcodes = selected_barcodes_dat$TargetBC
  cat("\n",'For Lane ', lane_number, ", group ",cluster_number, " " )
  cat("\n","We use: ", selected_barcodes)
  
  cells_in_dat = unique(original_csv$Cell)
  
  cat("\n",'Identify a set of cells possessing all selected barcodes')
  cells_filtered =  sapply(cells_in_dat, FUN = function(x){
    cell_has_selected_barcodes(cell = x, dat = original_csv,
                               selected_barcodes = selected_barcodes)
    
  })
  
  cells_IDs_filtered<- names(which(cells_filtered == TRUE))
  

  cat("\n",'Save cell IDS')
  write.csv(cells_IDs_filtered,file=paste0("runs_on_all_lanes_final/filtering/mGASv2_","Lane",lane_number,"_Clustering",clustering_level,"_Cluster",cluster_number,"_cell_IDs.csv"))
  
  
  
  cat("\n",'Subset of cells has: ', length(cells_IDs_filtered), " cells")
  

  dat_filtered_cells = original_csv[which(original_csv$Cell %in% cells_IDs_filtered), ]
  
  
  cat("\n",'Get number of cells per target')
  num <- c()
  for(target in unique(dat_filtered_cells$TargetBC)) { 
    #print(target)
    num <- c(num,length(dat_filtered_cells$Cell[which(dat_filtered_cells$TargetBC==target)]))
    
    
  }
  filtered_data <- data.frame(nums=num, targetBCS=unique(dat_filtered_cells$TargetBC))
  
  
  
  
  cat("\n",'Filter out all TBCs recovered for less than 2 cells')
  at_least_2_seq_TBCs <- filtered_data$targetBCS[which(filtered_data$nums > 1)]
  
  
  cat("\n",'Calculate entropy per TBCs')
  TBC_all_entropy <- data.frame()
  TBC_samples_all <- data.frame()
  for (TBC in at_least_2_seq_TBCs){
    TBC_temp <- filter(dat_filtered_cells, TargetBC == TBC) %>%
      mutate(Sites = paste(Site1,Site2,Site3,Site4,Site5,Site6))
    TBC_freq <- table(TBC_temp$Sites)/dim(TBC_temp)[1]
    TBC_entropy <- (-1)*sum(TBC_freq*log(TBC_freq))
    #TBC_all_entropy[TBC] <- TBC_entropy
    TBC_samples <- arrange(data.frame(table(TBC_temp$Sites)), desc(Freq))[1:10,]
    TBC_samples$TargetBC <- TBC
    TBC_samples$Entropy <- TBC_entropy
    TBC_num <- length(dat_filtered_cells$Cell[which(dat_filtered_cells$TargetBC==TBC)])
    TBC_all_entropy <- bind_rows(TBC_all_entropy, data.frame(TBC, TBC_entropy,TBC_num))
    TBC_samples_all <- bind_rows(TBC_samples_all, TBC_samples)
  }
  colnames(TBC_all_entropy) <- c('TargetBC','Entropy','Num. cells')
  
  cat("\n",'Filter out TBCs with entropy 0 (all sequences identical!)')
  filtered_tbcs <- TBC_all_entropy[TBC_all_entropy$Entropy != 0.0,]
  
  cat("\n",'Save stats for all usable TBCs sequences')
  write.csv(filtered_tbcs,file=paste0("runs_on_all_lanes_final/filtering/mGASv2_","Lane",lane_number,"_Clustering",clustering_level,"_Cluster",cluster_number,"_filtering_stats.csv"))
  
  cat("\n",'Save corresponding plot  usable TBCs sequences')
  edit_stats_plot <- ggplot(filtered_tbcs, aes(x=Entropy,y=`Num. cells`)) + geom_point() + ggtitle(paste0("Editing and recovery for Lane ", lane_number, ", Cluster ", cluster_number )) + 
    ylab("Num. cells with TargetBC") + xlab("TargetBC alignment entropy (-1*sum(sequence_frequencies*log(sequence_frequencies))")
  
  
  ggsave( paste0("runs_on_all_lanes_final/plots/edit_stats/mGASv2_edit_stats_plot_Lane", lane_number,"_Clustering",clustering_level, "_Cluster", cluster_number,".pdf"),edit_stats_plot)
  
  #extract the desired TBCS
  
  
  
  }
  
}


####################################################
# select which target bcs tp make into alignments. #
####################################################


convert_filter_and_write <- function(lane_number,clustering_level) {
#create a txt file that can be used as input for BEAST2
stats <- paste0("clustering",";","cluster_number",";","first_cellID",";","sampling_proportion",";","dimension_insertBCs",";","array_lenghts",";","initial_insertProb",";","frequencies",";","range_tbcs")
write(stats,file=paste0("runs_on_all_lanes_final/mGASv2_","Lane",lane_number,"_Clustering",clustering_level,"_filtering_stats.txt"))
#function input

#for loop 1 to number of clusters
for(cluster_number in 1:clustering_level) { 

filtered_tbcs_all <- read.csv(,file=paste0("runs_on_all_lanes_final/filtering/mGASv2_","Lane",lane_number,"_Clustering",clustering_level,"_Cluster",cluster_number,"_filtering_stats.csv"))



filtered_tbcs_all$percent_recovery <- filtered_tbcs_all$Num..cells / max(filtered_tbcs_all$Num..cells)

filtered_tbcs_all$product_entropy_cells <- filtered_tbcs_all$Entropy * filtered_tbcs_all$Num..cells

filtered_tbcs_all <- arrange(filtered_tbcs_all, desc(product_entropy_cells))

#create a threshold of entropty and recovery for selection
filtered_tbcs_final <- filter(filtered_tbcs_all,(Entropy>2) & (Num..cells > max(Num..cells)/5)) 

#top 15 tbcs with high entropy and coverage
if(length(filtered_tbcs_final$TargetBC) > 10) {
  filtered_tbcs_final <- filtered_tbcs_final[1:10,]
}


#convert to dataframe
filtered_tbcs_all <- as.data.frame(filtered_tbcs_all)

#keep order
filtered_tbcs_all$X <- factor(filtered_tbcs_all$X, levels=filtered_tbcs_all$X)

#color by discarded or kept 
filtered_tbcs_all$col <- "Discarded"
filtered_tbcs_all$col[filtered_tbcs_all$TargetBC %in% filtered_tbcs_final$TargetBC] <- "Selected"

filtered_tbcs_all$label <- filtered_tbcs_all$TargetBC
filtered_tbcs_all$label[!(filtered_tbcs_all$TargetBC %in% filtered_tbcs_final$TargetBC)] <- ""


selected_tbcs_plot_product <- ggplot(data=filtered_tbcs_all,aes(x=as.numeric(X), y=product_entropy_cells, label=label, col=col)) + geom_point() + geom_text(hjust=-0.1) +  theme_bw() + xlab("TargetBC") + ylab("Entropy x Num. cells")
selected_tbcs_plot_entropy <- ggplot(data=filtered_tbcs_all,aes(x=as.numeric(X), y=Entropy, label=label, col=col)) + geom_point() + geom_text(hjust=-0.1) +  theme_bw() + xlab("TargetBC") + ylab("Shannon entropy") + geom_hline(yintercept=2)
selected_tbcs_plot_num_cell <- ggplot(data=filtered_tbcs_all,aes(x=as.numeric(X), y=Num..cells, label=label, col=col)) + geom_point() + geom_text(hjust=-0.1) +  theme_bw() + xlab("TargetBC") + ylab("Num. cells")  + geom_hline(yintercept=max(filtered_tbcs_all$Num..cells)/5)


selected_tbcs_plot_numxentropy <- ggplot(data=filtered_tbcs_all,aes(x=Entropy, y=Num..cells, label=label, col=col)) + geom_point() + geom_text(hjust=-0.1) +  theme_bw() + xlab("Entropy") + ylab("Num. cells") + geom_vline(xintercept=2) + geom_hline(yintercept=max(filtered_tbcs_all$Num..cells)/5)


ggsave( paste0("runs_on_all_lanes_final/plots/edit_stats/mGASv2_selected_tbcs_plot_Lane", lane_number,"_Clustering",clustering_level, "_Cluster", cluster_number,".pdf"),cowplot::plot_grid(selected_tbcs_plot_product,selected_tbcs_plot_entropy,selected_tbcs_plot_num_cell,selected_tbcs_plot_numxentropy),width = 25,height=15)

edit_table <- import_and_clean_table(lane_number)


selected_tbcs <- filtered_tbcs_final$TargetBC

cells_IDs_filtered <- read.csv(paste0("runs_on_all_lanes_final/filtering/mGASv2_","Lane",lane_number,"_Clustering",clustering_level,"_Cluster",cluster_number,"_cell_IDs.csv"),header = T)

cells_IDs_filtered <- cells_IDs_filtered$x


dat_filtered_cells_tbcs = edit_table[which((edit_table$Cell %in% cells_IDs_filtered) & (edit_table$TargetBC %in% selected_tbcs) ), ]

if(length(dat_filtered_cells_tbcs$Cell) == 0) {
  print("=================NO CELLS IN FILTER===================")
  next
  
}

#remove the 6th column 
if(all(dat_filtered_cells_tbcs$Site6 == "None")) { 
  
  dat_filtered_cells_tbcs <- dat_filtered_cells_tbcs[,colnames(dat_filtered_cells_tbcs) != "Site6"   ]
  
}


converted_table <- convert_table_to_integers(dat_filtered_cells_tbcs,lane_number,clustering_level,cluster_number)$table
dimension_insertBCs <- convert_table_to_integers(dat_filtered_cells_tbcs,lane_number,clustering_level,cluster_number)$dim

#additionally, select ones with some editing
unedited_datasets <- c()
edited_datasets <- c()

truncated_1 <-c()
truncated_2 <- c()
truncated_3 <- c()
truncated_4 <- c()
full <- c()
array_lengths <- c()

write.csv(cbind(unique(converted_table$TargetBC)),paste0("runs_on_all_lanes_final/conversion/mGASv2_","Lane",lane_number,"_Clustering",clustering_level,"_Cluster",cluster_number,"_TargetBC_names_to_number.csv"))


for(targetBC in unique(converted_table$TargetBC)) {
  
  all_data_per_target <- converted_table[ (converted_table$TargetBC %in% targetBC), ]
  if((length(which(all_data_per_target$Site1 == 0))/length(all_data_per_target$Site1 == 0))==1) {
    unedited_datasets <-c(unedited_datasets,targetBC)
  }
  else { 
    
    print("truncation at 1")
    print((length(which(all_data_per_target$Site2 == 0))/length(all_data_per_target$Site2 == 0)))
    print("truncation at 2")
    print((length(which(all_data_per_target$Site3 == 0))/length(all_data_per_target$Site3 == 0)))
    print("truncation at 3")
    print((length(which(all_data_per_target$Site4 == 0))/length(all_data_per_target$Site4 == 0)))
    print("truncation at 4")
    print((length(which(all_data_per_target$Site5 == 0))/length(all_data_per_target$Site5 == 0)))
    
    edited_datasets <- c(edited_datasets,targetBC)
    if((length(which(all_data_per_target$Site2 == 0))/length(all_data_per_target$Site2 == 0))==1) {
     
      truncated_1 <-c(truncated_1,targetBC)
      array_lengths <- c(array_lengths,1)
    }
    else if((length(which(all_data_per_target$Site3 == 0))/length(all_data_per_target$Site3 == 0))==1) {
      
      truncated_2 <-c(truncated_2,targetBC)
      array_lengths <- c(array_lengths,2)
    }
    else if((length(which(all_data_per_target$Site4 == 0))/length(all_data_per_target$Site4 == 0))==1) {

      truncated_3 <-c(truncated_3,targetBC)
      array_lengths <- c(array_lengths,3)
    }
    else if((length(which(all_data_per_target$Site5 == 0))/length(all_data_per_target$Site5 == 0))==1) {
      
      truncated_4 <-c(truncated_4,targetBC)
      array_lengths <- c(array_lengths,4)
    }
    else { 
      full <- c(full,targetBC)
      array_lengths <- c(array_lengths,5)
    }
    
    } 
  
}


if(length(unique(converted_table$Cell)) > 700) { 
cell_sample <- sample(converted_table$Cell,700)
converted_table <- filter(converted_table,Cell %in% cell_sample)
    
}


#calculate the effective sampling proportion based on these filters (40000 or 35000:
sampling_proportion = length(unique(converted_table$Cell))/35000

#add filler tape lengths
missing_lengths <- 10 - length(array_lengths)
array_lengths_padded <- c(array_lengths,rep(0,missing_lengths))

#save filtering stats that can be read in for BEAST job submission
stats <- paste0(clustering_level,";",cluster_number,";",unique(converted_table$Cell)[1],";",sampling_proportion,";",dimension_insertBCs,";",paste(array_lengths_padded,collapse=";"),";",1/dimension_insertBCs,";",paste(c(1.0,rep(0,dimension_insertBCs -1 )),collapse = " "),";",paste(1:length(array_lengths),collapse=","))
write(stats,file=paste0("runs_on_all_lanes_final/mGASv2_","Lane",lane_number,"_Clustering",clustering_level,"_filtering_stats.txt"),append=TRUE)

create_nexus_alignments(edit_table = converted_table,paste0("mGASv2_Lane",lane_number,"_Clustering",clustering_level,"_Cluster",cluster_number),truncated_1,truncated_2,truncated_3,truncated_4)

}
}

#check editing in the sixth site
table_lane1 <- import_and_clean_table(1)
all(table_lane1$Site6 == "None")
table_lane2 <- import_and_clean_table(2)
all(table_lane2$Site6 == "None")
table_lane3 <- import_and_clean_table(3)
all(table_lane3$Site6 == "None")

#the "ideal" clustering level can be found with 
# in folder:clustering scripts Lane2_TargetBC_clustering_k3.pdf

create_stats_per_lane_per_cluster(lane_number=1,clustering_level=3)
convert_filter_and_write(1,3)

create_stats_per_lane_per_cluster(lane_number=2,clustering_level=3)
convert_filter_and_write(2,3)

create_stats_per_lane_per_cluster(lane_number=3,clustering_level=5)
convert_filter_and_write(3,5)



create_stats_per_lane_per_cluster(lane_number=1,clustering_level=8)
convert_filter_and_write(1,8)
create_stats_per_lane_per_cluster(lane_number=2,clustering_level=8)
convert_filter_and_write(2,8)
create_stats_per_lane_per_cluster(lane_number=3,clustering_level=8)
convert_filter_and_write(3,8)


