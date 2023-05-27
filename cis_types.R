#Step 5.2: Analyze the cis/trans-chromosonal fractions of the networks.

#The input fot his was the top 10,000,000 interactions of each network. 

library(tidyr)
library(dplyr)
library(ggplot2)
library(biomaRt)

#### Load 10M networks ####

TALL <- read.table("TALL_10M_MI.sif", header = TRUE)
NormalBMvsTALL <- read.table("NBM_TALL_10M_MI.sif", header = TRUE)
AML <- read.table("AML_10M_MI.sif", header = TRUE)
NormalBMvsAML <- read.table("NBM_AML_10M_MI.sif", header = TRUE)
MM <- read.table("MM_8th_MI.sif", header = TRUE)
NormalBMvsMM <- read.table("NBM_MM_10M_MI.sif", header = TRUE)
BALL <- read.table("BALL_10M_MI.sif", header = TRUE)
NormalBMvsBALL <- read.table("NBM_BALL_10M_MI.sif", header = TRUE)

#### Get an annotation file ####
#httr::set_config(httr::config(ssl_verifypeer = FALSE)) only on server 2
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "www")

features <- c("ensembl_gene_id", "chromosome_name", 
              "start_position", "end_position", "hgnc_symbol",	
              "percentage_gene_gc_content", "gene_biotype", "ensembl_gene_id_version", "hgnc_id")
chrs <- c(1:22, "X", "Y")

annot <- getBM(attributes = features,
               filters = "chromosome_name",
               values = chrs, 
               mart = ensembl)

colnames(annot)<-c("ensembl_gene_id", "Chr", "Start", "End", "HGNC_symbol", "GC", "Type", "Ensembl_ID_Version", "HGNC_ID")
annot$Length <- abs(annot$End - annot$Start)
dim(annot)
annot <- annot[!duplicated(annot$ensembl_gene_id), ]
dim(annot)
annot[annot == ""] <- NA
annot <- annot %>% drop_na
dim(annot)

# Add_Chromosomes_Types_SP <- function(x,y) {
#   
#   x <- x %>% dplyr::inner_join(y, by = c("Source"="Ensembl_ID_Version")) %>% 
#     dplyr::select("Source", "Sp_Coeff", "Target", "Chr", "Type") %>% dplyr::rename("Chr_Source"="Chr", "Type_Source"="Type") %>%
#     dplyr::inner_join(y, by = c("Target"="Ensembl_ID_Version")) %>% 
#     dplyr::select("Source", "Sp_Coeff", "Target", "Chr_Source", "Type_Source", "Chr", "Type") %>% dplyr::rename("Chr_Target"="Chr", "Type_Target"="Type") %>%
#     arrange(desc(Sp_Coeff))
#   
#   return(x)
# }

#### Prerp functions ####
Add_Chromosomes_Types_MI <- function(x,y) {
  
  x <- x %>% dplyr::inner_join(y, by = c("source"="ensembl_gene_id")) %>% 
    dplyr::select("source", "mi", "target", "Chr", "Type") %>% dplyr::rename("Chr_Source"="Chr", "Type_Source"="Type") %>%
    dplyr::inner_join(y, by = c("target"="ensembl_gene_id")) %>% 
    dplyr::select("source", "mi", "target", "Chr_Source", "Type_Source", "Chr", "Type") %>% dplyr::rename("Chr_Target"="Chr", "Type_Target"="Type") %>%
    arrange(desc(mi))
  
  return(x)
}

cis_counting <- function(x) {
  
  h <- 1
  j <- 1
  y <- 1
  cis_count <- 0
  
  while(h*10^(j) <= 1e+07) {
    for(i in y:nrow(x[1:(h*10^(j)),])) {
      if(x$Chr_Source[i] == x$Chr_Target[i]) { 
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

change_gene_types <- function(x) {
  print(table(as.factor(x$Type_Source)))
  x$Type_Source <- gsub(".*pseudogene", "pseudogene", x$Type_Source)
  x$Type_Source <- gsub("IG.*", "IG", x$Type_Source)
  x$Type_Source <- gsub("TR.*", "TR", x$Type_Source)
  print(table(as.factor(x$Type_Source)))
  
  print(table(as.factor(x$Type_Target)))
  x$Type_Target <- gsub(".*pseudogene", "pseudogene", x$Type_Target)
  x$Type_Target <- gsub("IG.*", "IG", x$Type_Target)
  x$Type_Target <- gsub("TR.*", "TR", x$Type_Target)
  print(table(as.factor(x$Type_Target)))
  
  return(x)
  
}

link_types <- function(x) {
  for(i in 1:nrow(x)) {
    x$link_type[i] <- paste(pmin(x$Type_Source[i],x$Type_Target[i]), 
                            pmax(x$Type_Source[i],x$Type_Target[i]),sep="-")
  }
  return(x)
}

#### Get gene chromosomes ####

TALL <- Add_Chromosomes_Types_MI(TALL, annot)
dim(TALL)
AML <- Add_Chromosomes_Types_MI(AML, annot)
dim(AML)
MM <- Add_Chromosomes_Types_MI(MM, annot)
dim(MM)
BALL <- Add_Chromosomes_Types_MI(BALL, annot)
dim(BALL)

NormalBMvsTALL <- Add_Chromosomes_Types_MI(NormalBMvsTALL, annot)
dim(NormalBMvsTALL)
NormalBMvsAML <- Add_Chromosomes_Types_MI(NormalBMvsAML, annot)
dim(NormalBMvsAML)
NormalBMvsMM <- Add_Chromosomes_Types_MI(NormalBMvsMM, annot)
dim(NormalBMvsMM)
NormalBMvsBALL <- Add_Chromosomes_Types_MI(NormalBMvsBALL, annot)
dim(NormalBMvsBALL)

#### Get cis counts ####

TALL_cis_count <- cis_counting(TALL)
AML_cis_count <- cis_counting(AML)
MM_cis_count <- cis_counting(MM)
BALL_cis_count <- cis_counting(BALL)

NormalBMvsTALL_cis_count <- cis_counting(NormalBMvsTALL)
NormalBMvsAML_cis_count <- cis_counting(NormalBMvsAML)
NormalBMvsMM_cis_count <- cis_counting(NormalBMvsMM)
NormalBMvsBALL_cis_count <- cis_counting(NormalBMvsBALL)

global_count <- data.frame(total_edges = NormalBMvsTALL_cis_count$total_edges, 
                           TALL_count = TALL_cis_count$cis_proportion,
                           AML_count = AML_cis_count$cis_proportion,
                           MM_count = MM_cis_count$cis_proportion,
                           BALL_count = BALL_cis_count$cis_proportion,
                           NBM_TALL_count = NormalBMvsTALL_cis_count$cis_proportion,
                           NBM_AML_count = NormalBMvsAML_cis_count$cis_proportion,
                           NBM_MM_count = NormalBMvsMM_cis_count$cis_proportion,
                           NBM_BALL_count = NormalBMvsBALL_cis_count$cis_proportion)
write.csv(global_count, file="global_cis_count_HC_8th_run.csv")

#### edit gene types tags, use 10k ####

TALL <- change_gene_types(TALL)
dim(TALL)
AML  <- change_gene_types(AML)
dim(AML)
MM  <- change_gene_types(MM)
dim(MM)
BALL  <- change_gene_types(BALL)
dim(BALL)

NormalBMvsTALL <- change_gene_types(NormalBMvsTALL)
dim(NormalBMvsTALL)
NormalBMvsAML  <- change_gene_types(NormalBMvsAML)
dim(NormalBMvsAML)
NormalBMvsMM  <- change_gene_types( NormalBMvsMM)
dim( NormalBMvsMM)
NormalBMvsBALL  <- change_gene_types( NormalBMvsBALL)
dim( NormalBMvsBALL)

#### get gene_type column ####
TALL_10k <- link_types(TALL[1:10000,])
dim(TALL_10k)
AML_10k  <- link_types(AML[1:10000,])
dim(AML_10k)
MM_10k  <- link_types(MM[1:10000,])
dim(MM_10k)
BALL_10k  <- link_types(BALL[1:10000,])
dim(BALL_10k)

NormalBMvsTALL_10k <- link_types(NormalBMvsTALL[1:10000,])
dim(NormalBMvsTALL_10k)
NormalBMvsAML_10k  <- link_types(NormalBMvsAML[1:10000,])
dim(NormalBMvsAML_10k)
NormalBMvsMM_10k  <- link_types( NormalBMvsMM[1:10000,])
dim(NormalBMvsMM_10k)
NormalBMvsBALL_10k  <- link_types( NormalBMvsBALL[1:10000,])
dim(NormalBMvsBALL_10k)

#### Get bar plots ####

TALL_EdgesTypes <- as.data.frame((table(as.factor(TALL_10k$link_type))))
TALL_EdgesTypes <- TALL_EdgesTypes %>% mutate(fraction = Freq/10000) %>% arrange(desc(fraction))

AML_EdgesTypes <- as.data.frame((table(as.factor(AML_10k$link_type))))
AML_EdgesTypes <- AML_EdgesTypes %>% mutate(fraction = Freq/10000) %>% arrange(desc(fraction))

MM_EdgesTypes <- as.data.frame((table(as.factor(MM_10k$link_type))))
MM_EdgesTypes <- MM_EdgesTypes %>% mutate(fraction = Freq/10000) %>% arrange(desc(fraction))

BALL_EdgesTypes <- as.data.frame((table(as.factor(BALL_10k$link_type))))
BALL_EdgesTypes <- BALL_EdgesTypes %>% mutate(fraction = Freq/10000) %>% arrange(desc(fraction))

NormalBMvsTALL_EdgesTypes <- as.data.frame((table(as.factor(NormalBMvsTALL_10k$link_type))))
NormalBMvsTALL_EdgesTypes <- NormalBMvsTALL_EdgesTypes %>% mutate(fraction = Freq/10000) %>% arrange(desc(fraction))

NormalBMvsAML_EdgesTypes <- as.data.frame((table(as.factor(NormalBMvsAML_10k$link_type))))
NormalBMvsAML_EdgesTypes <- NormalBMvsAML_EdgesTypes %>% mutate(fraction = Freq/10000) %>% arrange(desc(fraction))

NormalBMvsMM_EdgesTypes <- as.data.frame((table(as.factor(NormalBMvsMM_10k$link_type))))
NormalBMvsMM_EdgesTypes <- NormalBMvsMM_EdgesTypes %>% mutate(fraction = Freq/10000) %>% arrange(desc(fraction))

NormalBMvsBALL_EdgesTypes <- as.data.frame((table(as.factor(NormalBMvsBALL_10k$link_type))))
NormalBMvsBALL_EdgesTypes <- NormalBMvsBALL_EdgesTypes %>% mutate(fraction = Freq/10000) %>% arrange(desc(fraction))

Global_Edges_Types <- rbind(AML_EdgesTypes[1:5,], TALL_EdgesTypes[1:5,], BALL_EdgesTypes[1:5,], MM_EdgesTypes[1:5,], 
                            NormalBMvsAML_EdgesTypes[1:5,], NormalBMvsTALL_EdgesTypes[1:5,],
                            NormalBMvsBALL_EdgesTypes[1:5,], NormalBMvsMM_EdgesTypes[1:5,])

Global_Edges_Types[1:5,4] <- "AML"
Global_Edges_Types[6:10,4] <- "TALL"
Global_Edges_Types[11:15,4] <- "BALL"
Global_Edges_Types[16:20,4] <- "MM"
Global_Edges_Types[21:25,4] <- "AML"
Global_Edges_Types[21:25,4] <- "TALL"
Global_Edges_Types[21:25,4] <- "BALL"
Global_Edges_Types[21:25,4] <- "MM"

Global_Edges_Types[,4] <- c(rep("AML", 5), rep("TALL", 5), rep("BALL", 5), rep("MM", 5),
                            rep("AML", 5), rep("TALL", 5), rep("BALL", 5), rep("MM", 5))

Global_Edges_Types[,5] <- c(rep("Cancer", 20), rep("Normal", 20))
colnames(Global_Edges_Types)

Global_Edges_Types$V5 <- factor(Global_Edges_Types$V5, levels = c("Normal", "Cancer"))
Global_Edges_Types$V4 <- factor(Global_Edges_Types$V4, levels = c("AML", "BALL", "TALL", "MM"))

cols_barplot <- c("pseudogene-pseudogene" = "#33a02c", "protein_coding-pseudogene" = "#6a3d9a", "protein_coding-protein_coding"="#cab2d6",
                  "lncRNA-protein_coding" = "#ff7f00", "IG-IG" = "#fdbf6f", "lncRNA-lncRNA" = "#e31a1c", "snRNA-snRNA" = "#fb9a99",
                  "protein_coding-snoRNA" = "#ffff99", "protein_coding-scaRNA" = "#b2df8a", "misc_RNA-protein_coding" = "#1f78b4", 
                  "TR-TR" = "#a6cee3", "lncRNA-pseudogene" = "black")

tiff("RNA_biotype_MI_8th_3.tiff", units="in", width=12, height=10, res=300)
ggplot(Global_Edges_Types, aes(fill=Var1, x=V4, y=fraction)) + geom_bar(position="stack", stat="identity") + ylim(0,1) + 
  facet_grid(rows = vars(V5))  + scale_fill_manual(values = cols_barplot) + 
  theme(plot.title = element_text(face = "bold"), strip.text.x = element_text(size = 20)) + 
  ggtitle("Fractions of Top 5 interactions between interactions between RNA biotype") +
  labs(y = "Fraction of initeractions by biotype", x = "Cancer and respective control")

dev.off()

write.table(Global_Edges_Types, file = "Global_Edges_Types.txt", col.names = TRUE, row.names = FALSE)
