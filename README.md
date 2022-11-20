# Network-structure-of-hematopoietic-cancers



#Step 1: Pre-processing and normalization

The script is an example to perfom the pre-prosessing of the raw counts for the multiple myeloma dataset, and the next scripts presented here 
are written to analyze these sample set. In the script for analyzing the data, the inputs are files from every phenotype . The results presented
in the paper were made using HT-Seq counts. Unfortunnaly, the HT-Seq pipeline was discontinued from GDC, and replaced by the STAR - Count pipeline. 

This script has to be run individually for each HC in order to obtain their normalized expression values. The normal bone marrow data has to be used
each time the script is run. Hence, the final output for this are eight files.
 
#Step2: Differential gene expression analysis
The unnormalized expression matrix was used to detect differentially expressed genes using the DESeq2 package. This script is an example to detect DEGs in 
the multiple myeloma samples in comparison with the normal bone marrow data.

#Step 3: Get the co-expression for each gene pair

Following the normalization, we used the ARACNe algorithm to get the mutual information values for each pair of genes in the expression matrix. This was performed
as described in https://github.com/ddiannae/ARACNE-multicore

# Ste4: Functional enrichment

This scprit performs the functional enrichment of each community detected by HiDeF in Cytoscape. HiDeF was run using its default parameters ().





