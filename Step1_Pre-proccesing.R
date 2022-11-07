library(SummarizedExperiment)
library(TCGAbiolinks)
require(EDASeq)
require(dplyr)
require(NOISeq)
library(DESeq2)

###################################################################################################################################
#Download RNA Seq data from TCGA, projects encompass "MMRF-COMMPASS" for MM, "TARGET-ALL-P2" for BALL and TALL,
#and "TARGET-AML" for AML and Normal bone marrow.

qry.rna <- GDCquery(project = "MMRF-COMMPASS",
                    data.category= "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - Counts",)

GDCdownload(qry.rna)

MM.raw <- GDCprepare(qry.rna, summarizedExperiment = TRUE) #For MM
ALL.raw <- GDCprepare(qry.rna, summarizedExperiment = TRUE) #For BALL and TALL
AML_Normal_BM.raw <- GDCprepare(qry.rna, summarizedExperiment = TRUE) #For AML and Normal Bone Marrow

###################################################################################################################################
#Filter the data to retain only relevant samples.

MM <- MM.raw[ , MM.raw$sample_type == "Primary Blood Derived Cancer - Bone Marrow"]
rMM <- MM.raw[ , MM.raw$sample_type == "Recurrent Blood Derived Cancer - Bone Marrow"]

T_ALL <- ALL.raw[ , ALL.raw$primary_diagnosis == "T lymphoblastic leukemia/lymphoma" & 
                    ALL.raw$sample_type == "Primary Blood Derived Cancer - Bone Marrow"]
                    
B_ALL <- ALL.raw[ , ALL.raw$primary_diagnosis == "Precursor B-cell lymphoblastic leukemia" & 
                    ALL.raw$sample_type == "Primary Blood Derived Cancer - Bone Marrow"]
B_rALL <- ALL.raw[ , ALL.raw$primary_diagnosis == "Precursor B-cell lymphoblastic leukemia" & 
                     ALL.raw$sample_type == "Recurrent Blood Derived Cancer - Bone Marrow"]

AML_BM <- AML_Normal_BM.raw[, AML_Normal_BM.raw$sample_type == "Primary Blood Derived Cancer - Bone Marrow" | 
                          AML_Normal_BM.raw$sample_type =="Recurrent Blood Derived Cancer - Bone Marrow"] 
 
Normal_BoneMarrow <- AML_Normal_BM.raw[ , AML_Normal_BM.raw$sample_type == "Bone Marrow Normal"]

###################################################################################################################################
#Arrange raw expression matrix and factor tables.

rownames(MM) <- rowData(MM)$external_gene_name
rownames(rMM) <- rowData(rMM)$external_gene_name
rownames(Normal_BoneMarrow) <- rowData(Normal_BoneMarrow)$external_gene_name

#Bind raw expression matrices from the cancer phenotype and the normal phenotype into the object "rna"
rnas <- cbind(assay(MM), assay(rMM), assay(Normal_BoneMarrow))
rnas <- rnas[!duplicated(rownames(rnas)),]
rnas_before <- rnas 

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
ridx <- rowMeans(dataFilt) >= 10
dataFilt <- dataFilt[ridx, ]
rnas <- rnas[rownames(rnas) %in% rownames(dataFilt), ] 

###################################################################################################################################
#Get an annotation file with the features of the genes on the expression matrix

annot<-read.delim(file="../../../mart_export.txt", sep="\t")
names(annot)<-c("Gene.name", "Chr", "Start", "End", "GC", "Type", "ensembl_gene_id")
annot$Length <- abs(annot$End - annot$Start)
inter <- intersect(rownames(rnas), annot$Gene.name)
rnas1 <- rnas[rownames(rnas) %in% inter,] #This is the raw expression matrix used in Step 2 as input for DESeq2
annot <- annot[annot$Gene.name %in% inter,]
annot <- annot[!duplicated(annot$Gene.name),]

###################################################################################################################################
#Normalization steps.

ln.data <- withinLaneNormalization(rnas1, annot$Length, which = "full")
gcn.data <- withinLaneNormalization(ln.data , annot$GC, which = "full")
Btwn.Norm <- betweenLaneNormalization(gcn.data, which = "full") 
norm.counts <- tmm(Btwn.Norm, long = 1000, lc = 0, k = 0)
noiseqData <- NOISeq::readData(norm.counts, factors = Ready_factors)
mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = FALSE)
rnas_after <- exprs(mydata2corr1)

###################################################################################################################################
#Get PCAs

library(ggbiplot)

before.pca <- prcomp(t( ),center = TRUE,scale. = TRUE)
summary(before.pca)
ggbiplot(before.pca, var.axes=FALSE, ellipse=TRUE, groups=factors$Group)

after.pca <- prcomp(t(rnas_after),center = TRUE,scale. = TRUE)
summary(after.pca)
ggbiplot(after.pca, var.axes=FALSE, ellipse=TRUE, groups=factors$Group)


###################################################################################################################################
#Save normalizaed counts

MM_Norm <- rnas_after[, factors$Group=="MM"]
NormalBM_Norm <- rnas_after[, factors$Group=="NormalBM"] 

Aracne_MM_Norm <- cbind(rownames(MM_Norm), MM_Norm)
Aracne_NormalBM_Norm <- cbind(rownames(NormalBM_Norm), NormalBM_Norm)

colnames(Aracne_MM_Norm)[1] <- "gene"
colnames(Aracne_NormalBM_Norm)[1] <- "gene"

write.table(Aracne_MM_Norm, file = "rnas_norm_MM.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
write.table(Aracne_NormalBM_Norm, file = "rnas_norm_NormalBMvsMM.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)









