# total gene counts per sample ========================================

counts <- read_csv("stable_brains/raw_data/NaccCore_counts.csv")

library(tidyverse)

counts %>% 
  summarize_if(is.numeric,sum,na.rm = T) %>% 
  t() %>% 
  as.data.frame %>% 
  rename(genecounts = V1) %>% 
  rownames_to_column(var = 'sampleID') -> genecounts

library(ggplot2)
library(ggthemes)

genecounts %>% 
  ggplot(aes(genecounts)) +
  geom_histogram(bins = 40,alpha =0.5,color = 'grey') +
  #scale_x_continuous(breaks = pretty_breaks()) +
  #scale_y_continuous(breaks = pretty_breaks()) +
  theme_base() -> p_genecounts

print(p_genecounts)
# ggsave(p_genecounts,filename = 'results_figures/p_genecounts.png',width = 8, height = 5)
genecounts %>% 
  filter(genecounts < 800000) # S_22_P1F03(CA1)

genecounts %>% 
  filter(genecounts > 3000000) # V_31_P2A01 (MEA)



# filtering for each region separately 
rawcount_list







# DEseq2 & WGCNA ==============================================================

df <- counts 
filter_counts = 10

colnames(df)
coldata = data.frame(sampleID = as.character(colnames(df)[2:length(df)])) %>% 
  left_join (.,behav %>% 
               select(sampleID,region, status) %>% 
               mutate_if(is.character,factor), by = 'sampleID') %>% 
  column_to_rownames(var = 'sampleID')

# tidy count data format  
d.raw <- as.matrix(df[2:length(df)])
rownames(d.raw) <- df$ensgene

which(complete.cases(d.raw) == F) -> x

d.noNA <- d.raw[-x,]
countData <- d.noNA[rowSums(d.noNA > filter_counts) > 2, ]
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = coldata,
                              design = ~status+region)
dds_deseq <- DESeq(dds)

vsd <- vst(dds, blind=FALSE)

pca_region <- plotPCA(vsd,intgroup = c("region"))+
    theme_bw()

print(pca_region) # S_10 and S_28 (CA1 - ughhh they both are alpha... )
# ggsave(pca_region,filename = 'results_figures/pca_region.png',width = 8, height = 8) 


## remove all genes with counts < 15 in more than 75% of samples (33*0.75 -> 25)
## suggested by WGCNA on RNAseq FAQ
dds75 <- dds[ rowSums(counts(dds) >= 15) >= round((length(df)-1)*0.75), ]
nrow(dds75)   
## remove all genes with counts < 20 in more than 90% of samples (84*0.9=75.6)
## suggested by WGCNA on RNAseq FAQ
dds90 <- dds[ rowSums(counts(dds) >= 5) >= round((length(df)-1)*0.79), ]
nrow(dds90) 
## remove all genes with counts < 5in more than 80% of samples (84*0.8=67.2)
## suggested by another WGCNA pipeline
dds80 <- dds[ rowSums(counts(dds) >= 10) >= round((length(df)-1)*0.8), ]
nrow(dds80) 
 

 
## Perform the DESeq2 normalization, required before WGCNA analysis

dds2 <- DESeq(dds, betaPrior = FALSE, parallel = TRUE) #change here 

## Perform a variance-stabilizing transformation
vsd <- getVarianceStabilizedData(dds2)
## Many functions expect the matrix to be transposed
datExpr <- t(vsd) 
## check rows/cols
nrow(datExpr)
ncol(datExpr)
rownames(datExpr)
 
 
dds2x <- DESeq(dds90, betaPrior = FALSE, parallel = TRUE) #change here 
## Perform a variance-stabilizing transformation
vsd <- getVarianceStabilizedData(dds2x)
## Many functions expect the matrix to be transposed
datExpr <- t(vsd) 
## check rows/cols
nrow(datExpr)
ncol(datExpr)
rownames(datExpr)
 


Samples = rownames(datExpr);
collectGarbage()


# Re-cluster samples

pdf('results_figures/sampletree.pdf',width = 18, height = 5)
sampleTree2 = hclust(dist(datExpr), method = "average")
plot(sampleTree2)
dev.off()


# think getting rid of two alpha samples (S_10 and S_28) from CA1 region =========================

# 
# saveRDS(dds90, "results_RDS/dds90.RDS")
# ## Perform the DESeq2 normalization, required before WGCNA analysis
# 
# # dds2 <- DESeq(dds, betaPrior = FALSE, parallel = TRUE) #change here 
# # 
# # ## Perform a variance-stabilizing transformation
# # vsd <- getVarianceStabilizedData(dds2)
# # ## Many functions expect the matrix to be transposed
# # datExpr <- t(vsd) 
# # ## check rows/cols
# # nrow(datExpr)
# # ncol(datExpr)
# # rownames(datExpr)
# 
# 
# dds2x <- DESeq(dds90, betaPrior = FALSE, parallel = TRUE) #change here 
# ## Perform a variance-stabilizing transformation
# vsd <- getVarianceStabilizedData(dds2x)
# ## Many functions expect the matrix to be transposed
# datExpr <- t(vsd) 
# ## check rows/cols
# nrow(datExpr)
# ncol(datExpr)
# rownames(datExpr)
# 
# 
# 
# Samples = rownames(datExpr);
# collectGarbage()
# 
# 
# # Re-cluster samples
# 
# pdf('results_figures/sampletree_filtered.pdf',width = 18, height = 5)
# sampleTree2 = hclust(dist(datExpr), method = "average")
# plot(sampleTree2)
# dev.off()
# 
# ## save a copy 
# saveRDS(datExpr, "results_WGCNA/datExpr.RDS")
# 
# 
# 
