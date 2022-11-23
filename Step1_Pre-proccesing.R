# Step 1: Pre-processing and normalization

#This script has to be run individually for each HC in order to obtain their normalized expression values. The normal bone marrow data has to be used
#each time the script is run. Hence, the final output for this are eight files. 

#NOTE: The results presented in the paper were made using HT-Seq counts. Unfortunnaly, the HT-Seq pipeline was discontinued from GDC and replaced 
#by the STAR - Count pipeline.

###################################################################################################################################
#Load packages.

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
              "percentage_gene_gc_content", "gene_biotype")
chrs <- c(1:22, "X", "Y")

annot <- getBM(attributes = features,
      filters = "chromosome_name",
      values = chrs, 
      mart = ensembl)

colnames(annot)<-c("ensembl_gene_id", "Chr", "Start", "End", "HGNC_symbol", "GC", "Type")
annot$Length <- abs(annot$End - annot$Start)

###################################################################################################################################
#Download RNA Seq data from TCGA, projects encompass "MMRF-COMMPASS" for MM, "TARGET-ALL-P2" for BALL and TALL,
#and "TARGET-AML" for AML and Normal bone marrow.

qry.rna <- GDCquery(project = "MMRF-COMMPASS",
                    data.category= "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts")

GDCdownload(qry.rna)

MM.raw <- GDCprepare(qry.rna, summarizedExperiment = TRUE) #For MM
#ALL.raw <- GDCprepare(qry.rna, summarizedExperiment = TRUE) #For BALL and TALL
#AML_Normal_BM.raw <- GDCprepare(qry.rna, summarizedExperiment = TRUE) #For AML and Normal Bone Marrow

###################################################################################################################################
#Filter the data to retain only relevant samples.

MM <- MM.raw[ , MM.raw$sample_type == "Primary Blood Derived Cancer - Bone Marrow"]
rMM <- MM.raw[ , MM.raw$sample_type == "Recurrent Blood Derived Cancer - Bone Marrow"]

#T_ALL <- ALL.raw[ , ALL.raw$primary_diagnosis == "T lymphoblastic leukemia/lymphoma" & 
#                    ALL.raw$sample_type == "Primary Blood Derived Cancer - Bone Marrow"]
                    
#B_ALL <- ALL.raw[ , ALL.raw$primary_diagnosis == "Precursor B-cell lymphoblastic leukemia" & 
#                    ALL.raw$sample_type == "Primary Blood Derived Cancer - Bone Marrow"]
#B_rALL <- ALL.raw[ , ALL.raw$primary_diagnosis == "Precursor B-cell lymphoblastic leukemia" & 
#                     ALL.raw$sample_type == "Recurrent Blood Derived Cancer - Bone Marrow"]

#AML_BM <- AML_Normal_BM.raw[, AML_Normal_BM.raw$sample_type == "Primary Blood Derived Cancer - Bone Marrow" | 
#                          AML_Normal_BM.raw$sample_type =="Recurrent Blood Derived Cancer - Bone Marrow"] 
 
Normal_BoneMarrow <- AML_Normal_BM.raw[ , AML_Normal_BM.raw$sample_type == "Bone Marrow Normal"]

###################################################################################################################################

#Bind raw expression matrices from the cancer phenotype and the normal phenotype into the object "rna"
rnas <- cbind(assay(MM), assay(rMM), assay(Normal_BoneMarrow))

#Construct and object containing the sample name and the group of which it belongs, the samples must be in the same order as in the expression matrix.
factorsMM <- data.frame(Group = "MM", Sample = c(colnames(MM), colnames(rMM)))
factorsNormalBM <- data.frame(Group = "NormalBM", Sample = colnames(Normal_BoneMarrow))
factors <- rbind(factorsMM, factorsNormalBM)
rownames(factors) <- factors$Sample
Ready_factors <- as.data.frame(factors$Group)

###################################################################################################################################
#Filter low expressed genes.

dataFilt <- TCGAanalyze_Filtering(tabDF = rnas,
                                  method = "quantile",
                                  qnt.cut = 0.25)
threshold <- round(dim(rnas)[2]/2)
ridx <- rowSums(dataFilt == 0) <= threshold
dataFilt <- dataFilt[ridx, ]
dim(dataFilt)
ridx <- rowMeans(dataFilt) >= 10
dataFilt <- dataFilt[ridx, ]
print(dim(dataFilt))
rnas <- rnas[rownames(rnas) %in% rownames(dataFilt), ]
dim(rnas)

###################################################################################################################################
#Filter the annotation file to get only the genes in the expression matrix. Check for duplicates and remove them if necessary.

inter <- intersect(rownames(rnas), annot$ensembl_gene_id)
length(inter)
rnas1 <- rnas[rownames(rnas) %in% inter,] #This is the raw expression matrix used in Step 2 as input for DESeq2
dim(rnas1)
annot1 <- annot[annot$ensembl_gene_id  %in% inter,]
dim(annot1)
annot1 <- annot1[!duplicated(annot1$ensembl_gene_id),]
dim(annot1)
annot1[annot1 == ""] <- NA  

###################################################################################################################################
#Normalization steps.

ln.data <- withinLaneNormalization(rnas1, annot1$Length, which = "full")
gcn.data <- withinLaneNormalization(ln.data , annot1$GC, which = "full")
Btwn.Norm <- betweenLaneNormalization(gcn.data, which = "full") 
norm.counts <- tmm(Btwn.Norm, long = 1000, lc = 0, k = 0)
noiseqData <- NOISeq::readData(norm.counts, factors = Ready_factors)
mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = FALSE)
rnas2 <- exprs(mydata2corr1)

###################################################################################################################################
#Quality control

library(ggbiplot)

before.pca <- prcomp(t(rnas1),center = TRUE,scale. = TRUE)
summary(before.pca)
ggbiplot(before.pca, var.axes=FALSE, ellipse=TRUE, groups=factors$Group)

after.pca <- prcomp(t(rnas2),center = TRUE,scale. = TRUE)
summary(after.pca)
ggbiplot(after.pca, var.axes=FALSE, ellipse=TRUE, groups=factors$Group)

#QC Report

mydata_bf <- NOISeq::readData(
  data = rnas1,
  factors = factors,
  length = annot[,c("ensembl_gene_id", "Length")],
  biotype = annot[,c("ensembl_gene_id", "Type")],
  chromosome = annot[,c("Chr", "Start", "End")],
  gc = annot[, c("ensembl_gene_id", "GC")])

QCreport(mydata_bf, samples = NULL, factor = "Group", norm = FALSE)

mydata_after <- NOISeq::readData(
  data = rnas_after,
  factors = factors,
  length = annot[,c("ensembl_gene_id", "Length")],
  biotype = annot[,c("ensembl_gene_id", "Type")],
  chromosome = annot[,c("Chr", "Start", "End")],
  gc = annot[, c("ensembl_gene_id", "GC")])

QCreport(mydata_after, samples = NULL, factor = "Group", norm = TRUE)

###################################################################################################################################
#Save normalizaed counts

MM_Norm <- rnas2[, factors$Group=="MM"]
NormalBM_Norm <- rnas2[, factors$Group=="NormalBM"] 

Aracne_MM_Norm <- cbind(rownames(MM_Norm), MM_Norm)
Aracne_NormalBM_Norm <- cbind(rownames(NormalBM_Norm), NormalBM_Norm)

colnames(Aracne_MM_Norm)[1] <- "gene"
colnames(Aracne_NormalBM_Norm)[1] <- "gene"

write.table(Aracne_MM_Norm, file = "rnas_norm_MM.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
write.table(Aracne_NormalBM_Norm, file = "rnas_norm_NormalBMvsMM.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)









