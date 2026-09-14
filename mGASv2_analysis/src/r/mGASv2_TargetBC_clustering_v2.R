options(stringsAsFactors = F)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)
library(RColorBrewer)
library(dplyr)
library(ape)
library(phangorn)
library(tidyverse)

#read-in the extracted targetBC sequences
TargetBC_raw = read.csv('unfiltered_data/mGASv2_Lane3_CellByTape_10X_bamExtractV2_t3.csv', 
                         stringsAsFactors = F, header = T, na.strings=c("","NA"))

#filter out artifactual TBcs ("GGGGGGGGGGGG") and  ("AAAAAAAAAAAA"), and retain only TBcs with cumulative nreads > 100000
TargetBC_list <- filter(TargetBC_raw, mBC != "GGGGGGGGGGGG") %>%
  filter(mBC != "AAAAAAAAAAAA") %>%
  group_by(mBC) %>%
  summarise(read_total = sum(n_reads))

hist(TargetBC_list$read_total, xlim = c(0,10000000),breaks = 200)

#take the ~50 TBCs with the highest total reads?
TargetBC_list  <- filter(TargetBC_list,read_total > 100000)

hist(TargetBC_list$read_total, xlim = c(100000,10000000),breaks = 200)

TargetBC_list_UMI <- filter(TargetBC_raw, mBC != "GGGGGGGGGGGG") %>%
  filter(mBC != "AAAAAAAAAAAA") %>%
  group_by(mBC) %>%
  summarise(UMI_total = sum(n_UMI)) %>%
  filter(UMI_total > 10000)


hist(TargetBC_list_UMI$UMI_total,xlim = c(0,1000000),breaks = 400)

#Find a UMI threshold for the top 31 TBCs (matching the number of retained TBCs from the nreads based approach)
TargetBC_list_UMI  <- filter(TargetBC_list_UMI,UMI_total > 35000)

hist(TargetBC_list_UMI$UMI_total,xlim = c(35000,1000000),breaks = 400)

#check whether the top 31 UMI tbcs are also the top 31 nreads selected tbcs
length(intersect(TargetBC_list_UMI$mBC,TargetBC_list$mBC))
#28 out of 31 are identical. 

nTargetBC =  nrow(TargetBC_list)
# plot hist of reads x per sequences in the unfiltered set
hist(TargetBC_raw$n_reads, xlim = c(0,2000),breaks = 200)

#total sequences before nread
unfiltered_total <- nrow(TargetBC_raw)

# filter from the list of TBCs with total nreads> 100000, individual sequences with n_reads>1000
CellxTargetBC <- filter(TargetBC_raw, mBC %in% TargetBC_list$mBC) %>%
  filter(n_reads > 1000)

#find out the total sampling proportion from the reads based filtering
summed_and_individual_filtered <- nrow(CellxTargetBC)
summed_and_individual_filtered/unfiltered_total

# filter from the list of TBCs with total n_UMI > 100000, individual sequences with n_reads>1000
CellxTargetBC_UMI <- filter(TargetBC_raw, mBC %in% TargetBC_list_UMI$mBC) %>% filter(n_UMI > 15)

#find out the total sampling proportion from this filtering approach
summed_and_individual_filtered_UMI <- nrow(CellxTargetBC_UMI)
summed_and_individual_filtered_UMI/unfiltered_total

hist(CellxTargetBC$n_reads, xlim = c(1000,5000),breaks = 200)
plot(TargetBC_raw$n_reads,TargetBC_raw$n_UMI)

ggplot(medians,aes(x=median)) + geom_histogram(binwidth = 1)

#create an cell barcode overlap per targetBC 
TargetBC_overlaps <- matrix(0, nTargetBC, nTargetBC)
colnames(TargetBC_overlaps) <- TargetBC_list$mBC
rownames(TargetBC_overlaps) <- TargetBC_list$mBC

# extract all cBCs from the unfiltered set
CellBC_list_100Kreads = unique(unlist(TargetBC_raw$cBC))

#for each individual cell, go through the list of TBCs (if there is more than 1), increment the overlap by 1 for all of the TBCs. The matrix is a symmetric
for (CellBC in CellBC_list_100Kreads){
  ByCell_df <- filter(CellxTargetBC, cBC == CellBC)
  TBC_list <- unique(ByCell_df$mBC)
  if (length(TBC_list) > 1){
    for (ii in 1:(length(TBC_list)-1)){
      for (jj in (ii):length(TBC_list)){
        TargetBC_overlaps[TBC_list[ii],TBC_list[jj]] <- TargetBC_overlaps[TBC_list[ii],TBC_list[jj]] + 1
        TargetBC_overlaps[TBC_list[jj],TBC_list[ii]] <- TargetBC_overlaps[TBC_list[ii],TBC_list[jj]]
      }
    }
  }
}


#process the overlap matrix into a distance matrix 
TargetBC_DM <-1 / (TargetBC_overlaps+1)
diag(TargetBC_DM) <- 0
#build tree from that distance matrix
TargetBC_tree <- hclust(as.dist(TargetBC_DM), "complete")
#color palette
color_set3 <- brewer.pal(n = , name = "Set1")
#extract the clusters from the tree
member <- as.data.frame(cutree(TargetBC_tree, k = 8))
member$TargetBC <- rownames(member)
colnames(member) <- c('group','TargetBC')

TargetBC_tree <- as.phylo(hclust(as.dist(TargetBC_DM),"complete"))
plot(TargetBC_tree,tip.col=color_set3[member$group])

umap_DM <- uwot::umap(as.dist(TargetBC_DM))
plot(umap_DM, pch = 1)

pdf(file = "preprocessing_scripts/Lane3_TargetBC_clustering_k8.pdf")
plot(TargetBC_tree,tip.col=color_set3[member$group])
dev.off()

write.csv(member,file = "preprocessing_scripts/Lane3_TargetBC_cell_assignment_k8.csv")
