# The network structure of hematopoietic cancers

This repository contains the code for the paper: "The network structure of hematopoietic cancers". We organized this code into 5 Steps.

- Step 1: Data download, pre-processing and normalization. 
- Step 2: Differential expression analysis.
- Step 3: Constructuing the co-expression matrix.
- Step 4: Communities functional enrichment.
- Step 5: Topological analysis of the networks.

The steps 1-4 yield results for a specific phenotype, in the scripts presented here we used the multiple myeloma dataset. The phenotype of interest has to be defined according to the instructions in each script. Step 5 contains five scripts to analyse the networks, communities and functional enrichments:

- Step 5.1: Get the fractions of differentially expressed genes in each community of a network.
- Step 5.2: Analyze the cis/trans-chromosonal fractions of the networks.
- Step 5.3: Analyze the interactions among different RNA biotypes in each network.
- Step 5.4: Get the intersection networks (shared interactions among phenotypes).
- Step 5.5: Analyze the common biological processes among the four hematopoietic cancer networks.

Additionally, the script `SpMI_Similarity.R` calculates de similarity between the mutual information networks and their Spearman correlation counterparts.



