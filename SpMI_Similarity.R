#Similarity between the MI and Sp networks.

# This is script has as inputs the networks on the cut-off point of 10,000,000.


library(tidyr)
library(dplyr)
library(ggplot2)

S_T_unite <- function(x) {
#S_T_unite gets an ID for each coexpression interaction using the names of the genes. 
  
  
  x$S_T <- paste(pmin(x$source,x$target), 
                 pmax(x$source,x$target),sep="-")
  return(x)
}

####################### Get MI Networks #################################################
TALL_MI_10M <- read.table("TALL_MI_10M.sif", header = TRUE)
AML_MI_10M <- read.table("AML_MI_10M.sif", header = TRUE)
BALL_MI_10M <- read.table("BALL_MI_10M.sif", header = TRUE)
MM_MI_10M <- read.table("MM_MI_10M.sif", header = TRUE)

Normal_TALL_MI_10M <- read.table("NormalBMvsTALL_MI_10M.sif", header = TRUE)
Normal_AML_MI_10M <- read.table("NormalBMvsAML_MI_10M.sif", header = TRUE)
Normal_BALL_MI_10M <- read.table("NormalBMvsBALL_MI_10M.sif", header = TRUE)
Normal_MM_MI_10M <- read.table("NormalBMvsMM_MI_10M.sif", header = TRUE)

TALL_MI_interactions_10M <- S_T_unite(TALL_MI_10M)
AML_MI_interactions_10M <- S_T_unite(AML_MI_10M)
BALL_MI_interactions_10M <- S_T_unite(BALL_MI_10M)
MM_MI_interactions_10M <- S_T_unite(MM_MI_10M)

Normal_TALL_MI_interactions_10M <- S_T_unite(Normal_TALL_MI_10M)
Normal_AML_MI_interactions_10M <- S_T_unite(Normal_AML_MI_10M)
Normal_BALL_MI_interactions_10M <- S_T_unite(Normal_BALL_MI_10M)
Normal_MM_MI_interactions_10M <- S_T_unite(Normal_MM_MI_10M)

####################### Get Sp Networks #################################################
TALL_Sp_10M <- read.table("Ready_Spearman_TALL_10M.sif", header = FALSE)
AML_Sp_10M <- read.table("Ready_Spearman_AML_10M.sif", header = FALSE)
BALL_Sp_10M <- read.table("Ready_Spearman_BALL_10M.sif", header = FALSE)
MM_Sp_10M <- read.table("Ready_Spearman_MM_10M.sif", header = FALSE)

Normal_TALL_Sp_10M <- read.table("Ready_Spearman_NormalBMvsTALL_10M.sif", header = FALSE)
Normal_AML_Sp_10M <- read.table("Ready_Spearman_NormalBMvsAML_10M.sif", header = FALSE)
Normal_BALL_Sp_10M <- read.table("Ready_Spearman_NormalBMvsBALL_10M.sif", header = FALSE)
Normal_MM_Sp_10M <- read.table("Ready_Spearman_NormalBMvsMM_10M.sif", header = FALSE)

colnames(TALL_Sp_10M) <- c("source", "sp", "target")
colnames(AML_Sp_10M) <- c("source", "sp", "target")
colnames(BALL_Sp_10M) <- c("source", "sp", "target")
colnames(MM_Sp_10M) <- c("source", "sp", "target")

colnames(Normal_TALL_Sp_10M) <- c("source", "sp", "target")
colnames(Normal_AML_Sp_10M) <- c("source", "sp", "target")
colnames(Normal_BALL_Sp_10M) <- c("source", "sp", "target")
colnames(Normal_MM_Sp_10M) <- c("source", "sp", "target")

TALL_Sp_interactions_10M <- S_T_unite(TALL_Sp_10M)
AML_Sp_interactions_10M <- S_T_unite(AML_Sp_10M)
BALL_Sp_interactions_10M <- S_T_unite(BALL_Sp_10M)
MM_Sp_interactions_10M <- S_T_unite(MM_Sp_10M)

Normal_TALL_Sp_interactions_10M <- S_T_unite(Normal_TALL_Sp_10M)
Normal_AML_Sp_interactions_10M <- S_T_unite(Normal_AML_Sp_10M)
Normal_BALL_Sp_interactions_10M <- S_T_unite(Normal_BALL_Sp_10M)
Normal_MM_Sp_interactions_10M <- S_T_unite(Normal_MM_Sp_10M)

####################### Get Similarity #################################################

#Get the similarity for each phenotype. Load the cancer and normal network for each case accordingly.

Spearman_Cancer <- AML_Sp_interactions_10M
MI_Cancer <- AML_MI_interactions_10M
Spearman_Normal <- Normal_AML_Sp_interactions_10M
MI_Normal <- Normal_AML_MI_interactions_10M

h <- 1
j <- 1

while(h*10^(j) <= 1e+07) {
  
  Intesection_SpMI_1 <- intersect(Spearman_Cancer[1:(h*10^(j)),]$S_T, MI_Cancer[1:(h*10^(j)),]$S_T) 
  Cancer_SpMI <- MI_Cancer[(MI_Cancer$S_T %in% Intesection_SpMI_1), ]
  
  Normal_Intesection_SpMI_1 <- intersect(Spearman_Normal[1:(h*10^(j)),]$S_T, MI_Normal[1:(h*10^(j)),]$S_T)
  Normal_SpMI <- MI_Normal[(MI_Normal$S_T %in% Normal_Intesection_SpMI_1), ]
  
  if(h*10^(j) == 10) {
    
    SpMI_1 <- data.frame(Condition="Cancer",total_interactions=nrow(Spearman_Cancer[1:(h*10^(j)),]), shared_edges=(nrow(Cancer_SpMI)/(nrow(Spearman_Cancer[1:(h*10^(j)),]))))
    SpMI_2 <- data.frame(Condition="Normal",total_interactions=nrow(Spearman_Normal[1:(h*10^(j)),]), shared_edges=(nrow(Normal_SpMI)/(nrow(Spearman_Normal[1:(h*10^(j)),]))))
    
    z <- rbind(SpMI_1, SpMI_2)
  } else {
    SpMI_1 <- data.frame(Condition="Cancer",total_interactions=nrow(Spearman_Cancer[1:(h*10^(j)),]), shared_edges=(nrow(Cancer_SpMI)/(nrow(Spearman_Cancer[1:(h*10^(j)),]))))
    SpMI_2 <- data.frame(Condition="Normal",total_interactions=nrow(Spearman_Normal[1:(h*10^(j)),]), shared_edges=(nrow(Normal_SpMI)/(nrow(Spearman_Normal[1:(h*10^(j)),]))))
    
    z <- rbind(z, SpMI_1, SpMI_2)
  }
  
  if(h != 9) {
    h = h + 1
  } else {
    h = 1
    j = j + 1
  }
}

SharedEdges_SpMI_AML_log  <- z
#SharedEdges_SpMI_MM_log <- z
#SharedEdges_SpMI_BALL_log  <- z
#SharedEdges_SpMI_TALL_log  <- z


ggplot(SharedEdges_SpMI_AML_log, aes(x=total_interactions, y=shared_edges, color=Condition)) + geom_line(aes(color=Condition)) + geom_point() + 
  ylim(0, 1) + scale_colour_manual(values = cols) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold")) +
  ggtitle("Shared interactions between Spearman correlation and Mutual information networks", subtitle= "Acute myeloid leukemia") +
  labs(y = "Fraction of shared interactions", x = "Total edges") +  scale_x_log10() + ylim(0.5, 1)

ggsave("AML_log_sharedEdges.tiff", units="in", width=10, height=5, dpi=300, compression = 'lzw')

