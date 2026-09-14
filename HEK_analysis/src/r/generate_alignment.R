library(stringr)
library(dplyr)

pick_code <- function(edit,code_map) {
  code <- code_map[which(code_map$insert == edit), "integer"]
  print(code)
  return(code)
  
}


# get inserts across all cells, arrange in descrending freq,
# map each insert to integer;
# assumes that the inserts are already tri-nucleotides
get_insert_to_integer_map = function(edit_table, columns, dataset){
  
  insert_to_integer_map <- data.frame(table(unlist(edit_table[, columns]))) # %>% arrange(desc(Freq))

  insert_to_integer_map = insert_to_integer_map %>% arrange(desc(Freq))
   
  colnames(insert_to_integer_map)[1] = "insert"
  
  insert_to_integer_map$integer <- as.character(seq(0,19))
  
  
  return(insert_to_integer_map)
  
}

create_nexus_alignments <- function(edit_table,analysis_name){
  
  
  list_of_cells <- unique(edit_table$Cell)
  taxa <- length(list_of_cells)
  targets <- unique(edit_table$TargetBC)
  loss <- c()
  for(target in targets){
    #dir.create(file.path(paste0("data/",analysis_name,"alignment.targetBC=",target,".sciphy")), recursive = TRUE)
    sink(paste0("data/",analysis_name,".alignment.targetBC=",target,".sciphy"))
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
    
    if(target %in% tbcs_length_2) {
      cat("dimensions nchar=2;")
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
        
        if(target %in% tbcs_length_2) {
          cat("?,?")
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
        if(target %in% tbcs_length_2) {
          cat(paste((edit_table[(edit_table$Cell == cellID) & (edit_table$TargetBC == target),3:4]),collapse = ","))
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


edit_table = read.csv("Supplementary_File_2_DataTableMOI19.csv", stringsAsFactors = F,
                      header = T, na.strings=c("","NA"))



## filtering
### filter out entries that have no TargetBC or UMI
edit_table = edit_table[which(!(is.na(edit_table$TargetBC))), ]
edit_table = edit_table[which(!(is.na(edit_table$nUMI))), ]
edit_table = edit_table[which(!(is.na(edit_table$Cell))), ]
#length(unique(edit_table$Cell)) == 16810


### filter out entries that have not correctly ordered edits
#### pairwise comparisons, remove all entries without edit at site 1 and edit at site 2-5
edit_table = edit_table[! (is.na(edit_table[, 3]) &
                             !(is.na(edit_table[, 4]))), ]
edit_table = edit_table[! (is.na(edit_table[, 3]) &
                             !(is.na(edit_table[, 5]))), ]
edit_table = edit_table[! (is.na(edit_table[, 3]) &
                             !(is.na(edit_table[, 6]))), ]
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



site_columns = 3:7


frequency_of_target_bcs = as.data.frame(table(edit_table$TargetBC, useNA = "ifany"))
frequency_of_target_bcs = frequency_of_target_bcs[order(frequency_of_target_bcs$Freq, decreasing = T), ]
#frequent_target_bcs = frequency_of_target_bcs[1:17, "Var1"]

#write.csv(x = frequent_target_bcs, file = "./src/exploratory/13_frequent_TargetBCs.csv")


targetBCs_from_paper = c("ATGGTAAG", "ATTTATAT",
                         "ATTTGGTT", "GCAGGGTG",
                         "GTAAAGAT", "TAGATTTT",
                         "TGCGATTT", "TGGACGAC",
                         "TGGTTTTG", "TTAGATTG",
                         "TTCACGTA", "TTGAGGTG",
                         "TTTCGTGA")


site_columns = 3:7

original_csv <- edit_table
# 1st Convert 6-mers to 3-mers
# by removing the last 3 nucleotides that do not contain information
# (they are always GGA); replace <Na> by "NA"

for(site_column in site_columns) {
  
  original_csv[,site_column] <- substring(original_csv[, site_column], 1,3)
  original_csv[is.na(original_csv[, site_column]), site_column] <- "NA"
  
  
  
}

truncated_1 <-c()
truncated_2 <- c()
truncated_3 <- c()
truncated_4 <- c()
full <- c()

for(targetBC in unique(original_csv$TargetBC)) {
  
  all_data_per_target <- original_csv[ (original_csv$TargetBC %in% targetBC), ]
  
  if((length(which(all_data_per_target$Site2 == 0))/length(all_data_per_target$Site2 == 0))>0.99) {
    truncated_1 <-c(truncated_1,targetBC)
  }
  else if((length(which(all_data_per_target$Site3 == 0))/length(all_data_per_target$Site3 == 0))>0.99) {
    truncated_2 <-c(truncated_2,targetBC)
  }
  else if((length(which(all_data_per_target$Site4 == 0))/length(all_data_per_target$Site4 == 0))>0.99) {
    truncated_3 <-c(truncated_3,targetBC)
  }
  else if((length(which(all_data_per_target$Site5 == 0))/length(all_data_per_target$Site5 == 0))>0.99) {
    truncated_4 <-c(truncated_4,targetBC)
  }
  else { 
    full <- c(full,targetBC)
    }
  
  
  
  
}




# 2nd, create a map from the trinucleotides to an integer
insert_to_integer_map = get_insert_to_integer_map( original_csv, columns=site_columns, "cell_culture")
#insert_to_integer_map_full = get_insert_to_integer_map(original_csv, columns=site_columns, "cell_culture")

write.csv(insert_to_integer_map, file = "insert_to_integer_map.csv")

# convert these trinucleotides to integer in the edit table
for(site_column in site_columns) {
  original_csv[, site_column] <- unlist(lapply(original_csv[, site_column], function(x) {pick_code(x, insert_to_integer_map)}))
  
}


length(unique(original_csv$Cell))

set.seed(1)
cell_sample <- sample(unique(original_csv$Cell),1000)
sample <- original_csv[ (original_csv$Cell %in% cell_sample), ]


#create the BEAST2 input lines: 
cat(paste(unique(original_csv$TargetBC),collapse=","))
#ATTTGGTT,TTCACGTA,TAGATTTT,TTGTTTAC,ATTTATAT,TTAGATTG,TGGACGAC,TGGTTTTG,TGCGATTT,ACCTCGTG,TTTCGTGA,GTAAAGAT,TTGAGGTG,ATGGTAAG,GCAGGGTG,ATTAGTCA,ATCACGCT
cat(sample$Cell[1])
#AAGGTAAGTCATAACC-1
# targetBC used to initialise the tree, on the the following
unique(original_csv$TargetBC)
#"TTCACGTA"
cat(paste(truncated_2,collapse=","))
#TTCACGTA
cat(paste(truncated_4,collapse=","))
#TTGTTTAC,TGGACGAC,TGGTTTTG,ACCTCGTG,TTTCGTGA
cat(paste(full,collapse=","))
#ATTTGGTT,TAGATTTT,ATTTATAT,TTAGATTG,TGCGATTT,GTAAAGAT,TTGAGGTG,ATGGTAAG,GCAGGGTG,ATTAGTCA,ATCACGCT

# -D target_bcs="ATTTGGTT,TTCACGTA,TAGATTTT,TTGTTTAC,ATTTATAT,TTAGATTG,TGGACGAC,TGGTTTTG,TGCGATTT,ACCTCGTG,TTTCGTGA,GTAAAGAT,TTGAGGTG,ATGGTAAG,GCAGGGTG,ATTAGTCA,ATCACGCT" -D first_cell="AAGGTAAGTCATAACC-1" -D first_targetBC="TTCACGTA" -D range_length_2="TTCACGTA" -D range_length_4="TTGTTTAC,TGGACGAC,TGGTTTTG,ACCTCGTG,TTTCGTGA" -D range_length_5="TTTTGGTT,TAGATTTT,ATTTATAT,TTAGATTG,TGCGATTT,GTAAAGAT,TTGAGGTG,ATGGTAAG,GCAGGGTG,ATTAGTCA,ATCACGCT"


tbcs_length_4 <- truncated_4
tbcs_length_2 <- truncated_2

losses_stats <- create_nexus_alignments(sample,"HEK_missing")
stats <- data.frame(losses_stats)
write.csv(stats,file=paste0("HEK_missing","_filtering_stats.csv"))