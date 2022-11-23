# Step 3: Get the co-expression values for each gene pair

#We used the ARACNe algorithm to get the mutual information values for each pair of genes in the expression matrix. 
#This was performed as described in https://github.com/ddiannae/ARACNE-multicore

##This script calculates the Spearman correlation values. The co-expression matrix was obtained using the rcorr function in R as follows. 
#The input for this script is the normalized counts file from step 1.

library(dplyr)
library(Hmisc)
library(ggplot2)
library(reshape)

expr_matrix <- read.table(file = "rnas_norm_MM.tsv", sep = '\t', row.names = 1)

#Transpose the expression matrix
expr_matrix <- t(expr_matrix)

#Get the correlation values between genes in expr_matrix. The output, coexpr_matrix, is list of two square matrices. The r matrix contains
#the correlation values and the P matrix contains the p-values associated with the correlation. We set the lower triangle of the matrices to NA
#in order to remove duplicated correlations
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


