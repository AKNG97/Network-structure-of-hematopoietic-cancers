# Step 4: Functional enrichment

#This scprit performs the functional enrichment of each community detected by HiDeF in Cytoscape. HiDeF was run using default parameters 
#(Weight column: None, Max resolution parameter: 50, Consensus threshold: 75, Persistent threshold: 5, Algorith: Louvain). The files
#of the communities in each network is found in the Supplementary material folder.

library(clusterProfiler)
library(GOSemSim)
library(DOSE)
library(dplyr)

communities <- read.csv("../Communities_Diff_Expr/MM_GenesEnsembl_Communities.csv", row.names = 1)
communities_t <- t(communities)

for(i in 1:ncol(communities_t)) {
    
FuncEnrich <- enrichGO(gene = communities_t[,i], OrgDb = "org.Hs.eg.db", ont = "BP",
                       keyType = "ENSEMBL", pvalueCutoff=0.0000000001)

  #We applied the redundancy reduction function simplify to remove similar biological processes
  
  if(nrow(FuncEnrich) != 0) {
    simpl_FuncEnrich <- simplify(FuncEnrich, cutoff=0.7, by="p.adjust", select_fun=min) %>% as.data.frame()
    
    network <- data.frame(Community = rep(colnames(communities_t)[i], each=nrow(simpl_FuncEnrich)),
                          Interaction = rep(1, each=nrow(simpl_FuncEnrich)),
                          GO_BP = simpl_FuncEnrich$Description,
                          ID = simpl_FuncEnrich$ID,
                          GeneRatio = simpl_FuncEnrich$GeneRatio,
                          BgRatio = simpl_FuncEnrich$BgRatio,
                          pvalue = simpl_FuncEnrich$pvalue,
                          p.adjust = simpl_FuncEnrich$p.adjust,
                          qvalue = simpl_FuncEnrich$qvalue,
                          geneID = simpl_FuncEnrich$geneID,
                          Count = simpl_FuncEnrich$Count,
                          row.names=NULL)
    
    if(!exists("x")){
        x <- network
    } else {
      x <- rbind(x, network)
      
    }
  }

}

write.table(x, file = "NormalBMvsMM_GenesEnsembl_Communities_Enriched.sif", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)


