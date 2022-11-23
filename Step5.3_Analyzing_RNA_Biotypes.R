# Step 5.3: Analyze the interactions among different RNA biotypes in each network.

#This script was run using the top 10,000 interactions of each mutual information network, but it can be run using different cut-off points.

library(dplyr)
library(ggplot)

Network <- read.table("/Users/kenzuke/Documents/Paper_1/MM_MI_10k.sif", header = TRUE)
Annot <- read.csv("/Users/kenzuke/Documents/R_projects/Intersections_2/mart_export.txt")

#Detect the different subtypes of genes and group them into a bigger category.
table(as.factor(Annot$Gene.type))
Annot$Gene.type <- gsub(".*pseudogene", "pseudogene", Annot$Gene.type)
Annot$Gene.type <- gsub("IG.*", "IG", Annot$Gene.type)
Annot$Gene.type <- gsub("TR.*", "TR", Annot$Gene.type)
table(as.factor(Annot$Gene.type))

# Hay que revisar Annot para editar con Annot$Gene.type <- gsub("TR.*", "TR", Annot$Gene.type)
#Get the biotypes for each gene pair

Sources <- merge(Network[1:10000,], Annot, by.x = "source", by.y = "name")
Sources <- Sources %>% dplyr::select("source", "mi", "target", "Gene.type")
colnames(Sources) <- c("source", "mi", "target", "Gene.type.source")
Targets <- merge(Sources, Annot, by.x = "target", by.y = "name")
Targets <- Targets %>% dplyr::select("source", "mi", "target", "Gene.type.source", "Gene.type")
colnames(Targets) <- c("source", "mi", "target", "Gene.type.source", "Gene.type.target")

for(i in 1:nrow(Targets)) {
  
  if(Targets$Gene.type.source[i] == Targets$Gene.type.target[i]) {
    Targets$Homo_Hetero[i] <- "Homo"
    Targets$Genes.linked[i] <- paste(Targets$Gene.type.source[i])} else {
      Targets$Homo_Hetero[i] <- "Heter"
      Targets$Genes.linked[i] <- paste(pmin(Targets$Gene.type.source[i],Targets$Gene.type.target[i]), 
                                       pmax(Targets$Gene.type.source[i],Targets$Gene.type.target[i]),sep="-")
      
      Targets[order(abs(Targets$mi), decreasing = TRUE),]
    }
  
}

MM_10k_EdgesTypes <- as.data.frame((table(as.factor(Targets$Genes.linked))))
MM_10k_EdgesTypes <- MM_100k_EdgesTypes %>% mutate(fraction = Freq/100000) %>% arrange(desc(fraction))

#Repeat this steps for the four normal bone marrow files and construct a single object containing the top 5 interactions type in each network.

Global_Edges_Types <- rbind(AML_10k_EdgesTypes, MM_10k_EdgesTypes, TALL_10k_EdgesTypes, BALL_10k_EdgesTypes,
                 NormalAML_10k_EdgesTypes, NormalMM_10k_EdgesTypes, NormalTALL_10k_EdgesTypes, NormalBALL_10k_EdgesTypes)

Global_Edges_Types[1:20,4] <- "Cancer"
Global_Edges_Types[21:40,4] <- "Normal"

Global_Edges_Types$V4 = factor(Global_Edges_Types$V4, levels=c("Normal", "Cancer"))

Global_Edges_Types[1:5,5] <- "AML"
Global_Edges_Types[5:10,5] <- "MM"
Global_Edges_Types[11:15,5] <- "TALL"
Global_Edges_Types[16:20,5] <- "BALL"

Global_Edges_Types[21:25,5] <- "AML"
Global_Edges_Types[26:30,5] <- "MM"
Global_Edges_Types[31:35,5] <- "TALL"
Global_Edges_Types[36:40,5] <- "BALL"

ggplot(Global_Edges_Types, aes(fill=Var1, x=V5, y=fraction)) + geom_bar(position="stack", stat="identity") + ylim(0,1) +
  facet_grid(rows = vars(V4))  + scale_fill_brewer(palette = "Paired") +
  theme(plot.title = element_text(face = "bold")) +
  ggtitle("Fractions of Top 5 of interactions \nbetween RNA biotypes in the networks") +
  labs(y = "Fraction of initeractions by biotype", x = "Cancer and respective control")



