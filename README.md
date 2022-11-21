# The network structure of hematopoietic cancers



# Step 1: Pre-processing and normalization

The script is an example to perfom the pre-prosessing of the raw counts for the multiple myeloma dataset, and the next scripts presented here 
are written to analyze these sample set. In the script for analyzing the data, the inputs are files from every phenotype . The results presented
in the paper were made using HT-Seq counts. Unfortunnaly, the HT-Seq pipeline was discontinued from GDC, and replaced by the STAR - Count pipeline. 

This script has to be run individually for each HC in order to obtain their normalized expression values. The normal bone marrow data has to be used
each time the script is run. Hence, the final output for this are eight files.
 
# Step2: Differential gene expression analysis
The unnormalized expression matrix was used to detect differentially expressed genes using the DESeq2 package. This script is an example to detect DEGs in 
the multiple myeloma samples in comparison with the normal bone marrow data.

# Step 3: Get the co-expression for each gene pair

Following the normalization, we used the ARACNe algorithm to get the mutual information values for each pair of genes in the expression matrix. This was performed
as described in https://github.com/ddiannae/ARACNE-multicore

The Spearman correlation values were calculated using the script Step3_Coexpression_matrix.R,

We used the top 10,000, 100,000, and 10,000,000 to analyze the networks in different aspects.

# Step 4: Functional enrichment

This scprit performs the functional enrichment of each community detected by HiDeF in Cytoscape. HiDeF was run using default parameters (Weight column: None, Max resolution parameter: 50, Consensus threshold: 75, Persistent threshold: 5, Algorith: Louvain). The files
of the communities in each network is found in /../Outs . 

# Network analysis

## Similarity between mutual information and Spearman networks.

The input for this was the top 10,000,000 interactions of each network.

## Analysis of the inter/intra-chromosomal regulation

The input fot his was, again, the top 10,000,000 interactions of each mutual information network.

## Get the shared interactions among the four HCs networks

This script uses the top 100,000 interactions in each each mutual information network.

## Get the fractions of over/under-expressed genes in each community.

This script calculated the fraction of under- or over-expressed genes in each community. The files
of the communities in each network is found in " ". This is represented in the bipartite networks of communities and BPs in the form of pie plots in each community.

## Shared processes among HCs

This script uses the results of the Step 4, it identifies the common BPs and the genes that are responsible for the each enrichment in the communities of the networks. This information is summarized in a heatmap. The namesof the different branches in the dendrograms were added in a image editor.

## Analyzing the RNA biotypes in the data

This script was run using the top 10,000 interactions of each mutual information network. To analyze addional cut-off points of the noetworks, the corresponding networks can be loaded in the Network object.





