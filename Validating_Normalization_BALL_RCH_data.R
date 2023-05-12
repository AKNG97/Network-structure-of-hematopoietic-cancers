# The RCH dataset was retrieved form the article:
# https://ashpublications.org/bloodadvances/article/6/14/4093/485107/ALLSorts-an-RNA-Seq-subtype-classifier-for-B-cell

# This raw counts file contains mainly protein coding genes, hence we only used this RNA biotype in this analysis.

library(SummarizedExperiment)
library(TCGAbiolinks)
require(EDASeq)
require(dplyr)
require(NOISeq)
library(DESeq2)
library(biomaRt)

###################################################################################################################################
#Download annotation file from BioMart.

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

###################################################################################################################################
#Load expression matrices

############## Get Normal BM ###############
# As of the day that this script was run (14/12/2022), the TARGET-AML project on TCGA had more than 200 samples
# labaled as "Bone Marrow Normal". However, most of them correspond to patients with a cancer diagnosis. Here only included
# samples of patients without a diagnosis in order to assuare that the bone marrow samples was not affected by a malignancy 

qry.rna_Normal_BM <- GDCquery(project = "TARGET-AML",
                              data.category= "Transcriptome Profiling",
                              data.type = "Gene Expression Quantification",
                              workflow.type = "STAR - Counts",
                              sample.type = "Bone Marrow Normal")
GDCdownload(qry.rna_Normal_BM)
Normal_BoneMarrow <- GDCprepare(qry.rna_Normal_BM, summarizedExperiment = TRUE) 

############## Get BALL Bone Marrow ###############

qry.rna_ALL_BM <- GDCquery(project = "TARGET-ALL-P2",
                        data.category= "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts",
                        sample.type = c("Primary Blood Derived Cancer - Bone Marrow",
                                        "Recurrent Blood Derived Cancer - Bone Marrow"))
GDCdownload(qry.rna_ALL_BM)
ALL_BM <- GDCprepare(qry.rna_ALL_BM, summarizedExperiment = TRUE) 

BALL_BM <- ALL_BM[ , ALL_BM$primary_diagnosis == "Precursor B-cell lymphoblastic leukemia"]

############## Get BALL PM & RCH ###############

BALL_PM_RCH <- read.csv("combined_raw-counts.csv")
BALL_PM_RCH <- t(BALL_PM_RCH)

colnames(BALL_PM_RCH) <- BALL_PM_RCH["sample_id",]
BALL_PM_RCH <- BALL_PM_RCH[-1,]  %>% type.convert(as.is = TRUE) 

head(rownames(BALL_PM_RCH))
rownames(BALL_PM_RCH) <- gsub("\\.", "-" , rownames(BALL_PM_RCH))
rownames(BALL_PM_RCH) <- annot$Ensembl_ID_Version[match(rownames(BALL_PM_RCH), annot$HGNC_symbol)]
BALL_PM_RCH <- BALL_PM_RCH[!is.na(rownames(BALL_PM_RCH)),]
BALL_PM_RCH <- BALL_PM_RCH[!duplicated(rownames(BALL_PM_RCH)),]

BALL_RCH <- BALL_PM_RCH [,1:127]
BALL_PM <- BALL_PM_RCH[,128:195]

r_PM_RCH <- intersect(rownames(Normal_BoneMarrow), rownames(BALL_PM_RCH))

BALL_PM <- BALL_PM[rownames(BALL_PM) %in% r_PM_RCH, ]
BALL_RCH <- BALL_RCH[rownames(BALL_RCH) %in% r_PM_RCH, ]

#Get the factors and rnas objects for each normalization run  

rnas_BALL_BM_RCH <- cbind(BALL_RCH, assay(Normal_BoneMarrow)[rownames(Normal_BoneMarrow) %in% r_PM_RCH,])
rnas_BALL_BM_TCGA <- cbind(assay(BALL_BM), assay(Normal_BoneMarrow))

factors_BALL_RCH <- data.frame(Group = "BALL RCH", Sample =colnames(BALL_RCH))
factors_BALL_BM_TCGA <- data.frame(Group = "BALL BM TCGA", Sample =colnames(BALL_BM))
factors_NormalBM_TCGA <- data.frame(Group = "Normal BM TCGA", Sample =colnames(Normal_BoneMarrow))

factorsRCH_BM <- rbind(factors_BALL_RCH, factors_NormalBM_TCGA)
rownames(factorsRCH_BM) <- factorsRCH_BM$Sample
Ready_factorRCH_BM <- as.data.frame(factorsRCH_BM$Group)

factorsTCGA_BM <- rbind(factors_BALL_BM_TCGA, factors_NormalBM_TCGA)
rownames(factorsTCGA_BM) <- factorsTCGA_BM$Sample
Ready_factors_TCGA_BM<- as.data.frame(factorsTCGA_BM$Group)

#Filter

filter_TCGA <- function(x) {dataFilt <- TCGAanalyze_Filtering(tabDF = x,
                                                              method = "quantile",
                                                              qnt.cut = 0.25)
threshold <- round(dim(x)[2]/2)
ridx <- rowSums(dataFilt == 0) <= threshold
dataFilt <- dataFilt[ridx, ]
ridx <- rowMeans(dataFilt) >= 10
dataFilt <- dataFilt[ridx, ]
x <- x[rownames(x) %in% rownames(dataFilt), ]
print(dim(x))
return(x)
}

rnas_BM_RCH_filter <- filter_TCGA(rnas_BALL_BM_RCH)

rnas_BM_TCGA_filter <- filter_TCGA(rnas_BALL_BM_TCGA)


#Filter annot

get_annot <- function(x,y) {
  inter <- intersect(rownames(x), y$Ensembl_ID_Version)
  length(inter)
  x <- x[rownames(x) %in% inter,] #This is the raw expression matrix used in Step 2 as input for DESeq2
  print(dim(x))
  annot1 <- y[y$Ensembl_ID_Version  %in% inter,]
  print(dim(annot1))
  annot1 <- annot1[!duplicated(annot1$Ensembl_ID_Version),]
  print(dim(annot1))
  annot1[annot1 == ""] <- NA
  annot1 <- annot1[!is.na(annot1$HGNC_symbol),]
  annot1 <- annot1[annot1$Type=="protein_coding",]
  x <- x[rownames(x) %in% annot1$Ensembl_ID_Version,]
  annot1 <- annot1[match(rownames(x), annot1$Ensembl_ID_Version), ]
  print(dim(annot1))
  print(dim(x))
  return(list(x,annot1))
}



BM_RCH <- get_annot(rnas_BM_RCH_filter, annot)
BM_RCH_bf <- BM_RCH[[1]]
BM_RCH_annot <- BM_RCH[[2]]


BM_TCGA <- get_annot(rnas_BM_TCGA_filter, annot)
BM_TCGA_bf <- BM_TCGA[[1]]
BM_TCGA_annot <- BM_TCGA[[2]]


norm <- function(x, y, z) {
  ln.data <- withinLaneNormalization(x, y$Length, which = "full")
  gcn.data <- withinLaneNormalization(ln.data , y$GC, which = "full")
  Btwn.Norm <- betweenLaneNormalization(gcn.data, which = "full") 
  norm.counts <- tmm(Btwn.Norm, long = 1000, lc = 0, k = 0)
  noiseqData <- NOISeq::readData(norm.counts, factors = z)
  mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = FALSE)
  rnas2 <- exprs(mydata2corr1)
  return(rnas2)
}

BM_RCH_after <- norm(BM_RCH_bf, BM_RCH_annot, Ready_factorRCH_BM)

BM_TCGA_after <- norm(BM_TCGA_bf, BM_TCGA_annot, Ready_factors_TCGA_BM)


#RCH
BM_RCH_Norm <- BM_RCH_after[, factorsRCH_BM$Group=="BALL RCH"]
BM_RCH_Norm <- cbind(rownames(BM_RCH_Norm), BM_RCH_Norm)
colnames(BM_RCH_Norm)[1] <- "gene"
write.table(BM_RCH_Norm, file = "rnas_norm_BALL_RCH_ProtCoding.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)

Normal_BM_RCH_Norm <- BM_RCH_after[, factorsRCH_BM$Group=="Normal BM TCGA"]
Normal_BM_RCH_Norm <- cbind(rownames(Normal_BM_RCH_Norm), Normal_BM_RCH_Norm)
colnames(Normal_BM_RCH_Norm)[1] <- "gene"
write.table(Normal_BM_RCH_Norm, file = "rnas_norm_NormalBMvsBALL_RCH_ProtCoding.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)

#BM TCGA
BM_TCGA_Norm <- BM_TCGA_after[, factorsTCGA_BM$Group=="BALL BM TCGA"]
BM_TCGA_Norm <- cbind(rownames(BM_TCGA_Norm), BM_TCGA_Norm)
colnames(BM_TCGA_Norm)[1] <- "gene"
write.table(BM_TCGA_Norm, file = "rnas_norm_BALL_BM_TCGA_ProtCoding.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)

Normal_BM_TCGA_Norm <- BM_TCGA_after[, factorsTCGA_BM$Group=="Normal BM TCGA"]
Normal_BM_TCGA_Norm <- cbind(rownames(Normal_BM_TCGA_Norm), Normal_BM_TCGA_Norm)
colnames(Normal_BM_TCGA_Norm)[1] <- "gene"
write.table(Normal_BM_TCGA_Norm, file = "rnas_norm_NormalBMvsBALL_BM_TCGA_ProtCoding.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)

