# Step2: Differential gene expression analysis
#The unnormalized expression matrix was used to detect differentially expressed genes using the DESeq2 package.

library(DESeq2)

#The input for this script is the unnormalized counts and the factors object created in the Step 1.

#### MM analysis ####
dds <- DESeqDataSetFromMatrix(countData = round(MM_NBM[[1]]),
                              colData = factors_MM_NBM,
                              design = ~ Group)

dds <- DESeq(dds)

#Set the NormalBM group as the reference in the analysis
dds$Group <- relevel(dds$Group, ref = "NormalBM")
dds <- DESeq(dds)

# Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="Group_MM_vs_NormalBM", type="apeglm")

write.table(resLFC, file = "resLFC_MM_vs_NormalBM.tsv", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

#### AML analysis ####
dds <- DESeqDataSetFromMatrix(countData = round(AML_NBM[[1]]),
                              colData = factors_AML_NBM,
                              design = ~ Group)

dds <- DESeq(dds)

#Set the NormalBM group as the reference in the analysis
dds$Group <- relevel(dds$Group, ref = "NormalBM")
dds <- DESeq(dds)

# Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="Group_AML_vs_NormalBM", type="apeglm")

write.table(resLFC, file = "resLFC_AML_vs_NormalBM.tsv", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

#### BALL analysis ####
dds <- DESeqDataSetFromMatrix(countData = round(BALL_NBM[[1]]),
                              colData = factors_BALL_NBM,
                              design = ~ Group)

dds <- DESeq(dds)

#Set the NormalBM group as the reference in the analysis
dds$Group <- relevel(dds$Group, ref = "NormalBM")
dds <- DESeq(dds)

# Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="Group_BALL_vs_NormalBM", type="apeglm")

write.table(resLFC, file = "resLFC_BALL_vs_NormalBM.tsv", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

#### TALL analysis ####
dds <- DESeqDataSetFromMatrix(countData = round(TALL_NBM[[1]]),
                              colData = factors_TALL_NBM,
                              design = ~ Group)

dds <- DESeq(dds)

#Set the NormalBM group as the reference in the analysis
dds$Group <- relevel(dds$Group, ref = "NormalBM")
dds <- DESeq(dds)

# Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="Group_TALL_vs_NormalBM", type="apeglm")

write.table(resLFC, file = "resLFC_TALL_vs_NormalBM.tsv", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

