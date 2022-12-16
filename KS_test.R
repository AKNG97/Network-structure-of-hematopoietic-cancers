######## Paper data ##########

install.packages("dgof")
library("dgof")
library(dplyr)

#The input for this script is a file with the count of cis- fractions of every phenotype in different cut-off points, like in the following example:

# total_edges AML_count NormalBMvsAML_count BALL_count NormalBMvsBALL_count TALL_count NormalBMvsTALL_count  MM_count NormalBMvsMM_count
#          10      0.00                0.10  0.1000000                  0.1       0.10            0.1000000 0.1000000          0.1000000
#          20      0.05                0.30  0.4500000                  0.3       0.50            0.3500000 0.5000000          0.3000000
#          30      0.10                0.30  0.5333333                  0.3       0.60            0.3333333 0.6666667          0.3666667
#          40      0.15                0.30  0.5750000                  0.3       0.70            0.3250000 0.7500000          0.3750000
#          50      0.12                0.28  0.6000000                  0.3       0.66            0.3000000 0.7800000          0.3400000

counts_paper <- read.csv("/Users/kenzuke/Documents/compilation_cis_Paper.csv")

#AML
count_AML <- counts_paper %>% dplyr::select("total_edges", "AML_count", "NormalBMvsAML_count")

for(i in 1:(nrow(count_AML)-1)){
  ks <- suppressWarnings(ks.test(count_AML[[2]][i:55], count_AML[[3]][i:55]))
  KS_df <- data.frame(total_edges=count_AML$total_edges[i],
                      ks = ks$statistic,
                      ks_pval= ks$p.value)
  if(i==1) {
    AML_KS <- KS_df
  } else {
    AML_KS <- rbind(AML_KS, KS_df)
  }
}
colnames(AML_KS) <- c("total_edges", "ks_AML", "ks_pval_AML")

#BALL   
count_BALL <- counts_paper %>% dplyr::select("total_edges", "BALL_count", "NormalBMvsBALL_count")

for(i in 1:(nrow(count_BALL)-1)){
  ks <- suppressWarnings(ks.test(count_BALL[[2]][i:55], count_BALL[[3]][i:55]))
  KS_df <- data.frame(total_edges=count_BALL$total_edges[i],
                      ks = ks$statistic,
                      ks_pval= ks$p.value)
  if(i==1) {
    BALL_KS <- KS_df
  } else {
    BALL_KS <- rbind(BALL_KS, KS_df)
  }
}
colnames(BALL_KS) <- c("total_edges", "ks_BALL", "ks_pval_BALL" )

#TALL
count_TALL <- counts_paper %>% dplyr::select("total_edges", "TALL_count", "NormalBMvsTALL_count")

for(i in 1:(nrow(count_TALL)-1)){
  ks <- suppressWarnings(ks.test(count_TALL[[2]][i:55], count_TALL[[3]][i:55]))
  KS_df <- data.frame(total_edges=count_TALL$total_edges[i],
                      ks = ks$statistic,
                      ks_pval= ks$p.value)
  if(i==1) {
    TALL_KS <- KS_df
  } else {
    TALL_KS <- rbind(TALL_KS, KS_df)
  }
}
colnames(TALL_KS) <- c("total_edges", "ks_TALL", "ks_pval_TALL" )

#MM
count_MM <- counts_paper %>% dplyr::select("total_edges", "MM_count", "NormalBMvsMM_count")

for(i in 1:(nrow(count_MM)-1)){
  ks <- suppressWarnings(ks.test(count_MM[[2]][i:55], count_MM[[3]][i:55]))
  KS_df <- data.frame(total_edges=count_MM$total_edges[i],
                      ks = ks$statistic,
                      ks_pval= ks$p.value)
  if(i==1) {
   MM_KS <- KS_df
  } else {
    MM_KS <- rbind(MM_KS, KS_df)
  }
}
colnames(MM_KS) <- c("total_edges", "ks_MM", "ks_pval_MM" )


KS_compile <- AML_KS %>% inner_join(BALL_KS, by = "total_edges") %>% inner_join(TALL_KS, by = "total_edges") %>% inner_join(MM_KS, by = "total_edges")

write.csv(KS_compile, file="KS_compile.csv")
