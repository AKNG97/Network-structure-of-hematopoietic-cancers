library(tidyr)
library(dplyr)
library(ggplot2)

#### Get an annotation file ####

TALL <- read.table("/Users/kenzuke/Documents/R_projects/BulkEdges_Files/TALL_MI_10M.sif", header = TRUE)
NormalBMvsTALL <- read.table("/Users/kenzuke/Documents/R_projects/BulkEdges_Files/NormalBMvsTALL_MI_10M.sif", header = TRUE)
AML <- read.table("/Users/kenzuke/Documents/R_projects/BulkEdges_Files/AML_MI_10M.sif", header = TRUE)
NormalBMvsAML <- read.table("/Users/kenzuke/Documents/R_projects/BulkEdges_Files/NormalBMvsAML_MI_10M.sif", header = TRUE)
MM <- read.table("/Users/kenzuke/Documents/R_projects/BulkEdges_Files/MM_MI_10M.sif", header = TRUE)
NormalBMvsMM <- read.table("/Users/kenzuke/Documents/R_projects/BulkEdges_Files/NormalBMvsMM_MI_10M.sif", header = TRUE)
BALL <- read.table("/Users/kenzuke/Documents/R_projects/BulkEdges_Files/BALL_MI_10M.sif", header = TRUE)
NormalBMvsBALL <- read.table("/Users/kenzuke/Documents/R_projects/BulkEdges_Files/NormalBMvsBALL_MI_10M.sif", header = TRUE)

annot <- read.csv("/Users/kenzuke/Documents/R_projects/BulkEdges_Files/annot_genes_MI.csv")
annot$name <- gsub("-", ".", annot$name)

Edges_Chromosomes <- function(x,y) {
  Sources <- merge(x, y, by.x = "source", by.y = "name")
  Sources <- Sources %>% dplyr::select("source", "mi", "target", "Chromosome.scaffold.name")
  colnames(Sources) <- c("source", "mi", "target", "Chromosome_Source")
  Sources <- merge(Sources, y, by.x = "target", by.y = "name")
  Sources <- Sources %>% dplyr::select("source", "mi", "target", "Chromosome_Source", "Chromosome.scaffold.name")
  colnames(Sources) <- c("source", "mi", "target", "Chromosome_Source", "Chromosome_Target")
  
  return(Sources[order(abs(Sources$mi), decreasing = TRUE),])
}

TALL_chromosomes <- Edges_Chromosomes(TALL, annot)
AML_chromosomes <- Edges_Chromosomes(AML, annot)
MM_chromosomes <- Edges_Chromosomes(MM, annot)
BALL_chromosomes <- Edges_Chromosomes(BALL, annot)
DLBCL_chromosomes <- Edges_Chromosomes(DLBCL, annot)
BL_chromosomes <- Edges_Chromosomes(BL, annot)

NormalBMvsTALL_chromosomes <- Edges_Chromosomes(NormalBMvsTALL, annot)
NormalBMvsAML_chromosomes <- Edges_Chromosomes(NormalBMvsAML, annot)
NormalBMvsMM_chromosomes <- Edges_Chromosomes(NormalBMvsMM, annot)
NormalBMvsBALL_chromosomes <- Edges_Chromosomes(NormalBMvsBALL, annot)

cis_counting <- function(x) {
  
  h <- 1
  j <- 1
  y <- 1
  cis_count <- 0

      while(h*10^(j) <= 1e+07) {
            for(i in y:nrow(x[1:(h*10^(j)),])) {
              if(x$Chromosome_Source[i] == x$Chromosome_Target[i]) { 
                if(h*10^(j) == 10) {
                  cis_count = 1
                } else {
                  cis_count = 1 + cis_count
                }
              } else {
                  cis_count = cis_count
                }
              }
        if(h*10^(j) == 10) {
          cis_matrix <- tibble(total_edges = h*10^(j), cis = (cis_count), cis_proportion = (cis_count/(h*10^(j))))
          } else {
          a <- tibble(total_edges = h*10^(j), cis = (cis_count), cis_proportion = (cis_count/(h*10^(j))))
          cis_matrix <- rbind(cis_matrix, a)
          
        }
        if(h != 9) {
          y = h*10^(j) + 1
          h = h + 1
        } else {
          y = h*10^(j) + 1
          h = 1
          j = j + 1
        }
      }
  return(cis_matrix)
}

TALL_cis_count <- cis_counting(TALL_chromosomes)
AML_cis_count <- cis_counting(AML_chromosomes)
MM_cis_count <- cis_counting(MM_chromosomes)
BALL_cis_count <- cis_counting(BALL_chromosomes)
DLBCL_cis_count <- cis_counting(DLBCL_chromosomes)
BL_cis_count <- cis_counting(BL_chromosomes)

NormalBMvsTALL_cis_count <- cis_counting(NormalBMvsTALL_chromosomes)
NormalBMvsAML_cis_count <- cis_counting(NormalBMvsAML_chromosomes)
NormalBMvsMM_cis_count <- cis_counting(NormalBMvsMM_chromosomes)
NormalBMvsBALL_cis_count <- cis_counting(NormalBMvsBALL_chromosomes)

conditions <- c(rep("TALL", 55), rep("AML", 55), rep("MM", 55), rep("BALL", 55),
  rep("NormalBMvsTALL", 55), rep("NormalBMvsAML", 55), rep("NormalBMvsMM", 55), rep("NormalBMvsBALL", 55))

global_cis_count <- data.frame(condition=conditions, total_interactions=rep(TALL_cis_count$total_edges, 440), 
                               cis_fraction = c(TALL_cis_count$cis_proportion, AML_cis_count$cis_proportion,
                                                   MM_cis_count$cis_proportion, BALL_cis_count$cis_proportion,
                                                    NormalBMvsTALL_cis_count$cis_proportion, NormalBMvsAML_cis_count$cis_proportion,
                                                   NormalBMvsMM_cis_count$cis_proportion, NormalBMvsBALL_cis_count$cis_proportion), 
                               phenotype=c(rep("Cancer", 220), rep("Normal", 220)))


#### Graphs ####
cols <- c("BALL" = "#33a02c", "NormalBMvsBALL" = "#b2df8a",
          "NormalBMvsTALL" = "#fb9a99", "TALL" = "#e31a1c", "NormalBMvsMM" = "#fdbf6f", "MM" = "#ff7f00", "NormalBMvsAML" = "#cab2d6", 
          "AML" = "#6a3d9a")
cis_graph <- ggplot(global_cis_count[], aes(x=total_interactions, y=cis_fraction, color=condition)) + geom_point() + 
  ylim(0, 1) + scale_colour_manual(values = cols) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ggtitle("Fraction of Intra-chromosomal interactions calculated \nin Mutual information networks") +
  labs(y = "Fraction of intra-chromocomal interactions", x = "Total edges") + geom_line() +  scale_x_log10() + theme(plot.title = element_text(face = "bold")) +
  guides(col=guide_legend("Phenotypes"))



