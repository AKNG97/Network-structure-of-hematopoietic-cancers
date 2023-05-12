library(SummarizedExperiment)
library(TCGAbiolinks)
require(EDASeq)
require(dplyr)
require(NOISeq)
library(DESeq2)
library(biomaRt)

# The fist step is to obtain a file with the urls to download the data of each sample from St. Jude Cloud. This example shows the code for the normalization of AML samples

urls <- read.delim("file_containing_urls.txt",  header=FALSE)

for(url in urls) {
  
  download.file(url, destfile = basename(url))
}

#Get an expression matrix

files <- list.files(all.files=FALSE, full.names=FALSE)

for(i in 1:length(files)) {
  file_counts <- read.table(paste(files[i]), header = FALSE, row.names = 1)
  colnames(file_counts) <- paste(files[i]) %>% str_remove(regex(".RNA-Seq.feature-counts.txt"))
  
  if(i == 1) {
    expr_matrix <- file_counts
  } else {
    expr_matrix <- cbind(expr_matrix, file_counts)
  }
}

AML_raw_St <- expr_matrix

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
which(duplicated(annot$Ensembl_ID_Version))
annot1 <- annot[!duplicated(annot$Ensembl_ID_Version), ]
which(duplicated(annot$HGNC_symbol))
annot1 <- annot1[!duplicated(annot1$HGNC_symbol), ]
dim(annot1)
annot1[annot1 == ""] <- NA
annot1 <- annot1[!is.na(annot1$HGNC_symbol),]
dim(annot1)

###################################################################################################################################
##### Get normal BM #####

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

######## Change gene names to Ensembl on St Jude dataset

AML_raw_ensembl <- AML_raw

head(rownames(AML_raw_ensembl))
dim(AML_raw_ensembl)
which(duplicated(rownames(AML_raw_ensembl)))
genes_ST <- intersect(rownames(AML_raw_ensembl), annot1$HGNC_symbol)
AML_raw_ensembl <- AML_raw_ensembl[rownames(AML_raw_ensembl) %in% genes_ST,]
rownames(AML_raw_ensembl) <- annot1$Ensembl_ID_Version[match(rownames(AML_raw_ensembl), annot1$HGNC_symbol)]
dim(AML_raw_ensembl)
AML_raw_ensembl <- AML_raw_ensembl[!is.na(rownames(AML_raw_ensembl)),]
dim(AML_raw_ensembl)
AML_raw_ensembl <- AML_raw_ensembl[!duplicated(rownames(AML_raw_ensembl)),]
dim(AML_raw_ensembl)

###################################################################################################################################

NormalBM_AML <- Normal_BoneMarrow[rownames(Normal_BoneMarrow) %in% rownames(AML_raw_ensembl),]
dim(NormalBM_AML)
AML_raw_ensembl <- AML_raw_ensembl[rownames(AML_raw_ensembl) %in% rownames(Normal_BoneMarrow),]
dim(AML_raw_ensembl)
rnas_AML <- cbind(AML_raw_ensembl, assay(NormalBM_AML))

factors_AML_StJ <- data.frame(Group = "AML_StJ", Sample =colnames(AML_raw_ensembl))
factors_NormalBM <- data.frame(Group = "NormalBM", Sample =colnames(NormalBM_AML))

factors<- rbind(factors_AML_StJ, factors_NormalBM)
rownames(factors) <- factors$Sample
Ready_factor<- as.data.frame(factors$Group)

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

rnas_AML_filter <- filter_TCGA(rnas_AML)

get_annot <- function(x,y) {
  inter <- intersect(rownames(x), y$Ensembl_ID_Version)
  length(inter)
  x <- x[rownames(x) %in% inter,] #This is the raw expression matrix used in Step 2 as input for DESeq2
  print(dim(x))
  y <- y[y$Ensembl_ID_Version  %in% inter,]
  print(dim(y))
  y <- y[!duplicated(y$Ensembl_ID_Version),]
  print(dim(y))
  y[y == ""] <- NA
  y <- y[!is.na(y$HGNC_symbol),]
  x <- x[rownames(x) %in% y$Ensembl_ID_Version,]
  y <- y[match(rownames(x), y$Ensembl_ID_Version),]
  print(dim(y))
  print(dim(x))
  return(list(x,y))
}

AML <- get_annot(rnas_AML_filter, annot1)

AML_bf <- AML[[1]]
AML_annot <- AML[[2]]

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

AML_after <- norm(as.matrix(AML_bf), AML_annot, Ready_factor)

library(ggbiplot)

before.pca <- prcomp(t(AML_bf),center = TRUE,scale. = TRUE)
summary(before.pca)
ggbiplot(before.pca, var.axes=FALSE, ellipse=TRUE, groups=factors$Group)
ggsave("PCA_before.pdf")

after.pca <- prcomp(t(AML_after),center = TRUE,scale. = TRUE)
summary(after.pca)
ggbiplot(after.pca, var.axes=FALSE, ellipse=TRUE, groups=factors$Group)
ggsave("PCA_after.pdf")

################### Save normalized counts
AML_Norm <- AML_after[, factors$Group=="AML_StJ"]
AML_Norm <- cbind(rownames(AML_Norm), AML_Norm)
colnames(AML_Norm)[1] <- "gene"
write.table(AML_Norm, file = "rnas_norm_AML_StJ.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)

Normal_BM_TCGA_Norm <- AML_after[, factors$Group=="NormalBM"]
Normal_BM_TCGA_Norm <- cbind(rownames(Normal_BM_TCGA_Norm), Normal_BM_TCGA_Norm)
colnames(Normal_BM_TCGA_Norm)[1] <- "gene"
write.table(Normal_BM_TCGA_Norm, file = "rnas_norm_NormalBM_TCGAvs_AML_StJ.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)






