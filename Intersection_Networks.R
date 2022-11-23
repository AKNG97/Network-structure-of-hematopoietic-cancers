#Step 5.4: Get the intersection networks (shared interactions among phenotypes).

#This script uses the top 100,000 interactions in each each mutual information network.

library(tidyr)
library(dplyr)

##### Intersections of Cancers ####
#Load the top 100k interactions of each network.

TALL <- read.table("TALL_100k.sif", header = TRUE)
AML <- read.table("AML_100k.sif", header = TRUE)
BALL <- read.table("BALL_100k.sif", header = TRUE)
MM <- read.table("MM_100k.sif", header = TRUE)

# Getting shared interactions
S_T_unite <- function(x) {
  x %>% unite(S_T, sep = "-", c("source", "target"), remove = FALSE)
}

TALL_interactions <- S_T_unite(TALL)
AML_interactions <- S_T_unite(AML)
BALL_interactions <- S_T_unite(BALL)
MM_interactions <- S_T_unite(MM)

Cancer_intersections <- Reduce(intersect, list(TALL_interactions$S_T, AML_interactions$S_T, BALL_interactions$S_T, MM_interactions$S_T))
Cancer_intersections_network <- TALL_interactions[TALL_interactions$S_T %in% Cancer_intersections, ] %>% select("source", "mi", "target")
write.table(Cancer_intersections_network,file = "Cancer_network_intersections.sif", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)

##### Intersections of NormalBM ####
Normal_TALL <- read.table("NormalBMvsTALL.csv", header = TRUE)
Normal_AML <- read.table("NormalBMvsAML_100k.sif", header = TRUE)
Normal_BALL <- read.table("NormalBMvsBALL_100k.sif", header = TRUE)
Normal_MM <- read.table("NormalBMvsMM_100k.sif", header = TRUE)

Normal_TALL_interactions <- S_T_unite(Normal_TALL)
Normal_AML_interactions <- S_T_unite(Normal_AML)
Normal_BALL_interactions <- S_T_unite(Normal_BALL)
Normal_MM_interactions <- S_T_unite(Normal_MM)

Normal_intersections <- Reduce(intersect, list(Normal_TALL_interactions$S_T, Normal_AML_interactions$S_T, 
                                               Normal_BALL_interactions$S_T, Normal_MM_interactions$S_T))
Norma_intersections_network <- Normal_TALL_interactions[Normal_TALL_interactions$S_T %in% Normal_intersections, ] %>% select("source", "mi", "target")
write.table(Norma_intersections_network,file = "Normal_network_intersections.sif", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
