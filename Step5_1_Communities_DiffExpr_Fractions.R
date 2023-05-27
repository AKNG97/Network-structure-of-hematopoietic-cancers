#Step 5.1: Get the fractions of differentially expressed genes in each community of a network.

#This script uses ans inputs a community file and a table containing the differentially expressed genes of a phenotype.

library(dplyr)

communities <- as.data.frame(read.csv("communities_TALL_MI_8th.csv", row.names = 1))
communities[communities == "" | communities == " "] <- NA
DE <- read.table("resLFC_TALL_vs_NormalBM.tsv")
DE <- DE %>% select(baseMean, log2FoldChange)

for(i in 1:nrow(communities)) {

  query <- DE[rownames(DE) %in% communities[i,], ]
  
  CommNumGenes <- length(which(!is.na(communities[i,])))
  Num_Diff_up <- sum(query$log2FoldChange > 0)
  Num_Diff_down <- sum(query$log2FoldChange < 0)
  
  DE_Comm <- data.frame(Community_name = rownames(communities[i,]), 
                      Diff_up = Num_Diff_up/CommNumGenes, 
                      Diff_down = Num_Diff_down/CommNumGenes,
                      Num_Diff_up = Num_Diff_up,
                      Num_Diff_down = Num_Diff_down,
                      NA_DEG = (1 - ((Num_Diff_up + Num_Diff_down)/CommNumGenes)),
                      CommNumGenes = CommNumGenes,
                      DENumGenes = nrow(query))
  
     if(i==1){
          x <- DE_Comm
          } else {
                    x <- rbind(x, DE_Comm)
          } 
}

write.table(x, file = "TALL_Comm_DiffExpr.csv", row.names = FALSE, col.names = TRUE, quote = FALSE)
