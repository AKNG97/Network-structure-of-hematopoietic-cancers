#Step 5.5: Analyze the common biological processes among the four hematopoietic cancer networks.

#This script uses the results of the Step 4, it identifies the common BPs and the genes that are responsible
#for the each enrichment in the communities of the networks. This information is summarized in a heatmap. The names of 
#the different branches in the dendrograms were added in a image editor.


library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(stringr)

AML <- read.delim("/Users/kenzuke/Documents/Paper_1/AML_GenesEnsembl_Communities_Enriched.sif", header = FALSE)
BALL <- read.delim("/Users/kenzuke/Documents/Paper_1/BALL_GenesEnsembl_Communities_Enriched.sif", header = FALSE)
TALL <- read.delim("/Users/kenzuke/Documents/Paper_1/TALL_GenesEnsembl_Communities_Enriched.sif", header = FALSE)
MM <- read.delim("/Users/kenzuke/Documents/Paper_1/MM_GenesEnsembl_Communities_Enriched.sif", header = FALSE)

colnames(AML)<-c("Community", "val", "Description", "ID", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "count")
colnames(BALL)<-c("Community", "val", "Description", "ID", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "count")
colnames(TALL)<-c("Community", "val", "Description", "ID", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "count")
colnames(MM)<-c("Community", "val", "Description", "ID", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "count")

AML_cutoff <- AML[AML$p.adjust <= 1e-10, ] %>% select("Community", "p.adjust", "Description", "geneID")
BALL_cutoff <- BALL[BALL$p.adjust <= 1e-10, ] %>% select("Community", "p.adjust", "Description","geneID")
TALL_cutoff <- TALL[TALL$p.adjust <= 1e-10, ] %>% select("Community", "p.adjust", "Description","geneID")
MM_cutoff <- MM[MM$p.adjust <= 1e-10, ] %>% select("Community", "p.adjust", "Description","geneID")

intersection <- Reduce(intersect, list(AML_cutoff$Description, BALL_cutoff$Description, TALL_cutoff$Description, MM_cutoff$Description))

AML_InterProccs <- AML_cutoff[AML_cutoff$Description %in% intersection, ]
BALL_InterProccs <- BALL_cutoff[BALL_cutoff$Description %in% intersection, ]
TALL_InterProccs <- TALL_cutoff[TALL_cutoff$Description %in% intersection, ]
MM_InterProccs <- MM_cutoff[MM_cutoff$Description %in% intersection, ]

communts_procc <- function(x){
  for(i in 1:length(intersection)) {
    
    
    communities <- x$Community[x$Description==intersection[i]]
    table <- data.frame(process = intersection[i],
                        community = paste(unlist(communities), collapse =" "))
    
    if(i==1){
      y <- table
    } else {
      y <- rbind(y, table)
      
    }
    
  }
  return(y)
}

AML_communts_procc <- communts_procc(AML_InterProccs)
BALL_communts_procc <- communts_procc(BALL_InterProccs)
TALL_communts_procc <- communts_procc(TALL_InterProccs)
MM_communts_procc <- communts_procc(MM_InterProccs)

AML_communities <- read.csv("AML_communities_Ensembl.csv", row.names = 1)
BALL_communities <- read.csv("BALL_communities_Ensembl.csv")
TALL_communities <- read.csv("TALL_communities_Ensembl.csv")
MM_communities <- read.csv("MM_communities_Ensembl.csv")

gene_set <- function(x){
  for(i in 1:length(intersection)) { 
    a <- x[x$Description == intersection[i],]
    for(j in 1:nrow(a)) {
      b <- str_split(a$geneID[j], pattern = "/")[[1]]
      if(j == 1) {
        w <- b 
      } else {
        w <- union(w, b)
      }
      df <- data.frame(BP = intersection[i], Genes = paste(w, collapse = " "))}
    if(i == 1) {
      y <- df 
    } else {
      y <- rbind(y, df)
    }
  }
  return(y)
}

AML_genes_processes <- gene_set(AML_InterProccs)
BALL_genes_processes <- gene_set(BALL_InterProccs)
TALL_genes_processes <- gene_set(TALL_InterProccs)
MM_genes_processes <- gene_set(MM_InterProccs)

write.csv(AML_genes_processes, file = "AML_genes_processes.csv")
write.csv(BALL_genes_processes, file = "BALL_genes_processes.csv")
write.csv(TALL_genes_processes, file = "TALL_genes_processes.csv")
write.csv(MM_genes_processes, file = "MM_genes_processes.csv")

Split_AML_genes_processes <- read.csv("AML_genes_processes.csv", row.names = 1)
Split_AML_genes_processes[Split_AML_genes_processes == ""]  <- NA

Split_BALL_genes_processes <- read.csv("BALL_genes_processes.csv", row.names = 1)
Split_BALL_genes_processes[Split_BALL_genes_processes == ""]  <- NA

Split_TALL_genes_processes <- read.csv("TALL_genes_processes.csv", row.names = 1)
Split_TALL_genes_processes[Split_TALL_genes_processes == ""]  <- NA

Split_MM_genes_processes <- read.csv("MM_genes_processes.csv", row.names = 1)
Split_MM_genes_processes[Split_MM_genes_processes == ""]  <- NA

#Counting up/down genes

ensemble_to_symbol <- function(x) {
  i = 1
  while(i <= nrow(x)) {
    
    for(j in 1:ncol(x)) {
      
      if(is.na(x[i,j])){
        
        x[i,j] <- NA
        
      } else {
        x[i,j] <- annot[annot$ensembl_gene_id == x[i,j] ,]$external_gene_name
      }
    }
    
    i <- i +1
    
  }
  return(x)
  
}

rownames(Split_AML_genes_processes) <- Split_AML_genes_processes$BP
rownames(Split_TALL_genes_processes) <- Split_TALL_genes_processes$BP
rownames(Split_BALL_genes_processes) <- Split_BALL_genes_processes$BP
rownames(Split_MM_genes_processes) <- Split_MM_genes_processes$BP

Split_AML_genes_processes <- Split_AML_genes_processes %>% select(!BP)
Split_TALL_genes_processes <- Split_TALL_genes_processes %>% select(!BP)
Split_BALL_genes_processes <- Split_BALL_genes_processes %>% select(!BP)
Split_MM_genes_processes <- Split_MM_genes_processes %>% select(!BP)

AML_symbol <- ensemble_to_symbol(Split_AML_genes_processes)
BALL_symbol <- ensemble_to_symbol(Split_BALL_genes_processes)
TALL_symbol <- ensemble_to_symbol(Split_TALL_genes_processes)
MM_symbol <- ensemble_to_symbol(Split_MM_genes_processes)

count_up_down <- function(x,y) {
  for(i in 1:nrow(x)) {
    
    count_up = 0
    count_down = 0
    
    # Count gene set for the BP
    
    for(j in 1:ncol(x)) {
      DEi <- y[x[i,j],]$log2FoldChange
      
      if(!is.na(DEi)) {
        if(DEi > 0) {
          count_up = count_up + 1
        } else if (DEi < 0) {
          count_down = count_down + 1
        }
      }
      
      BP_expr_up <- count_up/(count_up + count_down)
      BP_expr_down <- count_down/(count_up + count_down)
      
    }
    
    df_2 <- data.frame(BP = intersection[i], 
                       DE_up = BP_expr_up, 
                       DE_down = BP_expr_down,
                       DE_sign = BP_expr_up - BP_expr_down)
    
    if(i == 1) {
      df_DE_2 <- df_2 
    } else {
      df_DE_2 <- rbind(df_DE_2, df_2)
    }
    
  }
  return(df_DE_2)
}

AML_DE_count <- count_up_down(AML_symbol, DE_AML)
BALL_DE_count <- count_up_down(BALL_symbol, DE_BALL)
TALL_DE_count <- count_up_down(TALL_symbol, DE_TALL)
MM_DE_count <- count_up_down(MM_symbol, DE_MM)

df_Global_DE_count <- data.frame(AML=AML_DE_count$DE_sign, BALL=BALL_DE_count$DE_sign, 
                                 TALL=TALL_DE_count$DE_sign, MM=MM_DE_count$DE_sign, row.names = AML_DE_count$BP)

abc_10 <- Heatmap(as.matrix(df_Global_DE_count), , width = unit(5, "cm"), row_names_gp = gpar(fontsize = 7.5), border = TRUE, col = c("#0072b2", "#FFFFFF", "#E62E00"), 
                  row_names_side = "left", row_dend_side = "right", name = "Differential \nExpression Trend", row_km = 5, row_km_repeats = 100, row_title_side = "right",
                  row_title_rot = 0, column_names_side = "top", row_title = c("Adaptive immune system", "Innate immune system", "Cell cycle",  "RNA processing", 
                                                                              "Adaptive immune system"))

draw(abc_10,merge_legend = TRUE,annotation_legend_side = "bottom")

tiff("Shared_BP.tiff", units="in", width=20, height=10, res=300)
draw(abc_10,merge_legend = TRUE,annotation_legend_side = "bottom")
dev.off()


