#The co-expression matrix using ARACNe (mutual information) was obtained using the pipeline developed by...

#The co-expression matrix using the Spearman correlation was obtained using the rcorr function in R as follows.

library(dplyr)
library(Hmisc)
library(ggplot2)
library(reshape)

annot<-read.delim(file="mart_export.txt", sep="\t")
expr_matrix <- read.table(file = "rnas_norm_MM.tsv", sep = '\t', row.names = 1)

#Sort genes in expr_matrix according with chromosome and start point
names(annot)<-c("Gene.name", "Chr", "Start", "End", "GC", "Type", "ensembl_gene_id")
annot <- annot %>%
  mutate(Chr = factor(Chr, levels=c(as.character(1:22), "X", "Y"))) %>%
  arrange(Chr, Start)
annot <- annot[!duplicated(annot$Gene.name),]

#Filter the annot object to retain only the genes in expr_matrix
annot <- annot %>% filter(Gene.name %in% rownames(expr_matrix))

#Change the order of the genes in expr_matrix to correspond with the order in annot object
expr_matrix <- expr_matrix[annot$Gene.name, ]
expr_matrix <- t(expr_matrix)

#Get the correlation values between genes in expr_matrix. The output, coexpr_matrix, is list of two square matrices. The r matrix contains
#the correlation values and the P matrix contains the p-values associated with the correlation. We set the lower triangle of the matrices to NA
# in order to remove duplicated correlations
coexpr_matrix <- rcorr(expr_matrix, type="spearman")

coexpr_val <- coexpr_matrix$r
coexpr_val[lower.tri(coexpr_val, diag = TRUE)] <- NA

coexpr_P <- coexpr_matrix$P
coexpr_P[lower.tri(coexpr_P, diag = TRUE)] <- NA

#Converting the square matrix into a list of unique values
val_melt <- melt(coexpr_val)
val_melt <- val_melt[!is.na(val_melt$value), ]

P_melt <- melt(coexpr_P)
P_melt <- P_melt[!is.na(P_melt$value), ]

colnames(val_melt) <- c("Source", "Target", "Sp_Coeff")
val_melt <- val_melt %>% relocate(Sp_Coeff, .after=Source)
Sp_matrix <- cbind(val_melt, P_melt$value)
colnames(Sp_matrix) <- c("Source","Sp_Coeff", "Target", "P_value")

#We added a column with the absolute value of the correlation to sort them according to this
Sp_matrix$abs <- abs(Sp_matrix$Sp_Coeff)
sorted <- Sp_matrix[order(Sp_matrix$abs, decreasing = TRUE),]

write.table(sorted, file = "Spearman_MM.sif", row.names = FALSE, sep = "\t", quote = FALSE)


