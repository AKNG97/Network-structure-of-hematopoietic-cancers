#NOTE 1: The results presented in the paper were made using HT-Seq counts. Unfortunnaly, the HT-Seq pipeline was discontinued 
#from GDC and replaced by the STAR - Count pipeline.

#NOTE 2: All data used here was downloaded during autum of 2021. Since then there has been
# updates to the datasets and now includes more samples in the case of data from TARGET-AML.
# However, we decided to only include the original data that were processed using the HT-seq counts pipeline.

#### Load libraries ####

library(SummarizedExperiment)
library(TCGAbiolinks)
require(EDASeq)
require(dplyr)
require(NOISeq)
library(DESeq2)
library(biomaRt)

#### Prepare functions ####

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

get_annot <- function(x,y) {
  inter <- intersect(rownames(x), y$HGNC_symbol)
  length(inter)
  annot1 <- y[y$HGNC_symbol  %in% inter,]
  print(dim(annot1))
  annot1 <- annot1[!duplicated(annot1$HGNC_symbol),]
  print(dim(annot1))
  x <- x[rownames(x) %in% annot1$HGNC_symbol,]
  print(dim(x))
  x <- x[!duplicated(rownames(x)),] #
  print(dim(x))
  print(head(rownames(x)))
  print(head(annot1$HGNC_symbol))
  annot1 <- annot1[match(rownames(x), annot1$HGNC_symbol), ]
  print(dim(annot1))
  print(dim(x))
  print(head(rownames(x)))
  print(head(annot1$HGNC_symbol))
  return(list(x,annot1))
}

norm <- function(x, y, z) {
  ln.data <- withinLaneNormalization(x, y$Length, which = "full")
  gcn.data <- withinLaneNormalization(ln.data , y$GC, which = "full")
  Btwn.Norm <- betweenLaneNormalization(gcn.data, which = "full") 
  norm.counts <- tmm(Btwn.Norm, long = 1000, lc = 0, k = 0)
  noiseqData <- NOISeq::readData(norm.counts, factors = as.data.frame(z$Group))
  mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = FALSE)
  rnas2 <- exprs(mydata2corr1)
  return(rnas2)
}

Get_raw_matrix <- function(x,y) {
  z <- cbind(assay(x), assay(y))
  dim(z)
  head(rownames(z))
  rownames(z) <- rowData(x)$external_gene_name
  head(rownames(z))
  x <- x[!duplicated(rownames(z)),]
  dim(z)
  head(rownames(z))
  return(z)
  }

Get_factors_objects <- function(x,y){
  factors <- rbind(x, y)
  rownames(factors) <- factors$Sample
  return(factors)
}

get_norm_matrices <- function(x, a, y, z) { 
  rownames(x) <- a$ensembl_gene_id[match(rownames(x), a$HGNC_symbol)]
  x <- x[, y$Group==z]
  print(dim(x))
  x <- cbind(rownames(x), x)
  colnames(x)[1] <- "gene"
  print(dim(x))
  return(x)
}

#### Get annotation file ####
httr::set_config(httr::config(ssl_verifypeer = FALSE))
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
annot <- annot[!duplicated(annot$ensembl_gene_id),]

#### Download data ####

qry.AML <- GDCquery(project = "TARGET-AML",
                    data.category= "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts")
GDCdownload(qry.AML)
AML_Normal_BM <- GDCprepare(qry.AML, summarizedExperiment = TRUE) #For MM

qry.ALL <- GDCquery(project = "TARGET-ALL-P2",
                    data.category= "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts")
GDCdownload(qry.ALL)
ALL <- GDCprepare(qry.ALL, summarizedExperiment = TRUE) #For MM

qry.MM <- GDCquery(project = "MMRF-COMMPASS",
                    data.category= "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts")
GDCdownload(qry.MM)
MM_raw <- GDCprepare(qry.MM, summarizedExperiment = TRUE) #For MM

#AML_Normal_BM <- readRDS("AML_Normal_BM.rds")
#ALL <- readRDS("../../rnas_raw_TAP2.RDS")
#MM_raw <- readRDS("rnas_raw_MM.rds")

AML_BM <- AML_Normal_BM[, AML_Normal_BM$sample_type == "Primary Blood Derived Cancer - Bone Marrow" | AML_Normal_BM$sample_type =="Recurrent Blood Derived Cancer - Bone Marrow"] 
B_ALL <- ALL[ , (ALL$primary_diagnosis == "Precursor B-cell lymphoblastic leukemia" & ALL$sample_type == "Primary Blood Derived Cancer - Bone Marrow") |
                (ALL$primary_diagnosis == "Precursor B-cell lymphoblastic leukemia" & ALL$sample_type == "Recurrent Blood Derived Cancer - Bone Marrow")]
MM <- MM_raw[ , (MM_raw$sample_type == "Primary Blood Derived Cancer - Bone Marrow")|
                (MM_raw$sample_type == "Recurrent Blood Derived Cancer - Bone Marrow")]
T_ALL <- ALL[ , ALL$primary_diagnosis == "T lymphoblastic leukemia/lymphoma" & ALL$sample_type == "Primary Blood Derived Cancer - Bone Marrow"]
Normal_BoneMarrow <- AML_Normal_BM[ , AML_Normal_BM$sample_type == "Bone Marrow Normal"]

##### Get raw matrices ####

MM_NBM <- Get_raw_matrix(MM, Normal_BoneMarrow)
AML_NBM <- Get_raw_matrix(AML_BM, Normal_BoneMarrow)
BALL_NBM <- Get_raw_matrix(B_ALL, Normal_BoneMarrow)
TALL_NBM <- Get_raw_matrix(T_ALL, Normal_BoneMarrow)

#### Get factors objects ####

factorsAML <- data.frame(Group = "AML", Sample = colnames(AML_BM))
factorsTALL <- data.frame(Group = "TALL", Sample = colnames(T_ALL))
factors_BALL <- data.frame(Group = "BALL", Sample = colnames(B_ALL))
factorsMM <- data.frame(Group = "MM", Sample = colnames(MM))
factors_NBM <- data.frame(Group = "NormalBM", Sample = colnames(Normal_BoneMarrow))

factors_MM_NBM <- Get_factors_objects(factorsMM, factors_NBM)
factors_BALL_NBM <- Get_factors_objects(factors_BALL, factors_NBM)
factors_TALL_NBM <- Get_factors_objects(factorsTALL, factors_NBM)
factors_AML_NBM <- Get_factors_objects(factorsAML, factors_NBM)

factors_MM_NBM$Group <- as.factor(factors_MM_NBM$Group)
factors_BALL_NBM$Group <- as.factor(factors_BALL_NBM$Group)
factors_TALL_NBM$Group <- as.factor(factors_TALL_NBM$Group)
factors_AML_NBM$Group <- as.factor(factors_AML_NBM$Group)


#### Filter ####
MM_NBM <- filter_TCGA(MM_NBM)
AML_NBM <- filter_TCGA(AML_NBM)
BALL_NBM <- filter_TCGA(BALL_NBM)
TALL_NBM <- filter_TCGA(TALL_NBM)

#### Get annotation objects ####

MM_NBM <- get_annot(MM_NBM, annot)
TALL_NBM <- get_annot(TALL_NBM, annot)
BALL_NBM <- get_annot(BALL_NBM, annot)
AML_NBM <- get_annot(AML_NBM, annot)

#### normalization ####

MM_NBM_norm <- norm(MM_NBM[[1]], MM_NBM[[2]], factors_MM_NBM)
BALL_NBM_norm <- norm(BALL_NBM[[1]], BALL_NBM[[2]], factors_BALL_NBM)
TALL_NBM_norm <- norm(TALL_NBM[[1]], TALL_NBM[[2]], factors_TALL_NBM)
AML_NBM_norm <- norm(AML_NBM[[1]], AML_NBM[[2]], factors_AML_NBM)

#### QC ####

#QC bf
mydata_bf_MM <- NOISeq::readData(
  data = MM_NBM[[1]],
  factors = factors_MM_NBM,
  length = MM_NBM[[2]][,c("HGNC_symbol", "Length")],
  biotype = MM_NBM[[2]][,c("HGNC_symbol", "Type")],
  chromosome = MM_NBM[[2]][,c("Chr", "Start", "End")],
  gc = MM_NBM[[2]][, c("HGNC_symbol", "GC")])

QCreport(mydata_bf_MM, samples = NULL, factor = "Group", norm = FALSE)

mydata_bf_BALL <- NOISeq::readData(
  data = BALL_NBM[[1]],
  factors = factors_BALL_NBM,
  length = BALL_NBM[[2]][,c("HGNC_symbol", "Length")],
  biotype = BALL_NBM[[2]][,c("HGNC_symbol", "Type")],
  chromosome = BALL_NBM[[2]][,c("Chr", "Start", "End")],
  gc = BALL_NBM[[2]][, c("HGNC_symbol", "GC")])

QCreport(mydata_bf_BALL, samples = NULL, factor = "Group", norm = FALSE)

mydata_bf_TALL <- NOISeq::readData(
  data = TALL_NBM[[1]],
  factors = factors_TALL_NBM,
  length = TALL_NBM[[2]][,c("HGNC_symbol", "Length")],
  biotype = TALL_NBM[[2]][,c("HGNC_symbol", "Type")],
  chromosome = TALL_NBM[[2]][,c("Chr", "Start", "End")],
  gc = TALL_NBM[[2]][, c("HGNC_symbol", "GC")])

QCreport(mydata_bf_TALL, samples = NULL, factor = "Group", norm = FALSE)

mydata_bf_AML <- NOISeq::readData(
  data = AML_NBM[[1]],
  factors = factors_AML_NBM,
  length = AML_NBM[[2]][,c("HGNC_symbol", "Length")],
  biotype = AML_NBM[[2]][,c("HGNC_symbol", "Type")],
  chromosome = AML_NBM[[2]][,c("Chr", "Start", "End")],
  gc = AML_NBM[[2]][, c("HGNC_symbol", "GC")])

QCreport(mydata_bf_AML, samples = NULL, factor = "Group", norm = FALSE)

#QC after

mydata_after_MM <- NOISeq::readData(
  data = MM_NBM_norm,
  factors = factors_MM_NBM,
  length = MM_NBM[[2]][,c("HGNC_symbol", "Length")],
  biotype = MM_NBM[[2]][,c("HGNC_symbol", "Type")],
  chromosome = MM_NBM[[2]][,c("Chr", "Start", "End")],
  gc = MM_NBM[[2]][, c("HGNC_symbol", "GC")])

QCreport(mydata_after_MM, samples = NULL, factor = "Group", norm = TRUE)

mydata_after_BALL <- NOISeq::readData(
  data = BALL_NBM_norm,
  factors = factors_BALL_NBM,
  length = BALL_NBM[[2]][,c("HGNC_symbol", "Length")],
  biotype = BALL_NBM[[2]][,c("HGNC_symbol", "Type")],
  chromosome = BALL_NBM[[2]][,c("Chr", "Start", "End")],
  gc = BALL_NBM[[2]][, c("HGNC_symbol", "GC")])

QCreport(mydata_after_BALL, samples = NULL, factor = "Group", norm = TRUE)

mydata_after_TALL <- NOISeq::readData(
  data = TALL_NBM_norm,
  factors = factors_TALL_NBM,
  length = TALL_NBM[[2]][,c("HGNC_symbol", "Length")],
  biotype = TALL_NBM[[2]][,c("HGNC_symbol", "Type")],
  chromosome = TALL_NBM[[2]][,c("Chr", "Start", "End")],
  gc = TALL_NBM[[2]][, c("HGNC_symbol", "GC")])

QCreport(mydata_after_TALL, samples = NULL, factor = "Group", norm = TRUE)

mydata_after_AML <- NOISeq::readData(
  data = AML_NBM_norm,
  factors = factors_AML_NBM,
  length = AML_NBM[[2]][,c("HGNC_symbol", "Length")],
  biotype = AML_NBM[[2]][,c("HGNC_symbol", "Type")],
  chromosome = AML_NBM[[2]][,c("Chr", "Start", "End")],
  gc = AML_NBM[[2]][, c("HGNC_symbol", "GC")])

QCreport(mydata_after_AML, samples = NULL, factor = "Group", norm = TRUE)

#### Get norm expression matrices files ####
MM_expr <- get_norm_matrices(MM_NBM_norm, MM_NBM[[2]], factors_MM_NBM, "MM")
BALL_expr <- get_norm_matrices(BALL_NBM_norm, BALL_NBM[[2]], factors_BALL_NBM, "BALL")
TALL_expr <- get_norm_matrices(TALL_NBM_norm, TALL_NBM[[2]], factors_TALL_NBM, "TALL")
AML_expr <- get_norm_matrices(AML_NBM_norm, AML_NBM[[2]], factors_AML_NBM, "AML")

NBM_MM_expr <- get_norm_matrices(MM_NBM_norm, MM_NBM[[2]],factors_MM_NBM, "NormalBM")
NBM_BALL_expr <- get_norm_matrices(BALL_NBM_norm, BALL_NBM[[2]],factors_BALL_NBM, "NormalBM")
NBM_TALL_expr <- get_norm_matrices(TALL_NBM_norm, TALL_NBM[[2]],factors_TALL_NBM, "NormalBM")
NBM_AML_expr <- get_norm_matrices(AML_NBM_norm, AML_NBM[[2]],factors_AML_NBM, "NormalBM")

write.table(MM_expr, file = "rnas_norm_MM.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
write.table(BALL_expr, file = "rnas_norm_BALL.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
write.table(TALL_expr, file = "rnas_norm_TALL.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
write.table(AML_expr, file = "rnas_norm_AML.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)

write.table(NBM_MM_expr, file = "rnas_norm_NBM_MM.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
write.table(NBM_BALL_expr, file = "rnas_norm_NBM_BALL.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
write.table(NBM_TALL_expr, file = "rnas_norm_NBM_TALL.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
write.table(NBM_AML_expr, file = "rnas_norm_NBM_AML.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)

#### save bf ####
saveRDS(MM_NBM, file="MM_NBM_bf.RDS")
saveRDS(BALL_NBM, file="BALL_NBM_bf.RDS")
saveRDS(TALL_NBM, file="TALL_NBM_bf.RDS")
saveRDS(AML_NBM, file="AML_NBM_bf.RDS")
