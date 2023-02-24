# total gene counts per sample ========================================

library(tidyverse)
library(ggthemes)
library(scales)
library(DESeq2)

# bring in raw counts
counts <- read_csv("stable_brains/raw_data/NaccCore_counts.csv")

# check for extremes
counts %>% 
  summarize_if(is.numeric,sum,na.rm = T) %>% 
  t() %>% 
  as.data.frame %>% 
  rename(genecounts = V1) %>% 
  rownames_to_column(var = 'sampleID') -> genecounts


genecounts %>% 
  ggplot(aes(genecounts)) +
  geom_histogram(bins = 40,alpha =0.5,color = 'grey') +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme_base() -> p_genecounts

p_genecounts





# filtering for each region separately 
rawcount_list
Nacc <- rawcount_list$NaccCore
Nacc






# DEseq2 & WGCNA ==============================================================

df <- counts 
filter_counts = 10

colnames(df)

# subjectID/sampleID/status.
coldata <-
full_join(
behavior %>% select(subjectID, glicko_rank, status),
sample_id %>% select(subjectID, sampleID,region)
) %>%
  filter(region=="NaccCore") %>% 
    column_to_rownames(var = 'sampleID')

# 
# coldata = data.frame(sampleID = as.character(colnames(df)[2:length(df)])) %>% 
#   left_join (.,behavior %>% 
#                select(sampleID,region, status) %>% 
#                mutate_if(is.character,factor), by = 'sampleID') %>% 
#   column_to_rownames(var = 'sampleID')




# tidy count data format  
d.raw <- as.matrix(df[2:length(df)])
rownames(d.raw) <- df$ensgene

# which(complete.cases(d.raw) == F) -> x
# d.noNA <- d.raw[-x,]

## Decide how you want to filter.
# e.g. could just do 10 across all subjects per gene (rowSums)

# the below only keeps genes where we have >90% of samples having at least 10 counts.

ncol(d.raw) # number of subjects
mincount <- 10
n_above_mincount <- apply(d.raw,1, function(x) sum(x>mincount,na.rm = T))
countData <- d.raw[n_above_mincount > (.9*ncol(d.raw)),]

countData

#make sure rownames and colnames of samples are in the same order.
all(row.names(coldata) %in% colnames(countData)) #if TRUE all samples are same in both but not necessarily in order
countData <- countData[,rownames(coldata)]
all(rownames(coldata) == colnames(countData)) #should be TRUE


## Put into DESeq matrix object
# warning about factors is fine.
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = coldata,
                              design = ~status)


#Normalize
dds_deseq <- DESeq(dds)
dds_deseq

#Variances stable transformation
#step before function for PCA plots
vsd <- vst(dds, blind=FALSE)

pca_region <- plotPCA(vsd,intgroup = c("status"))+
    theme_bw()

pca_region



## More filtering code (look at later)

## remove all genes with counts < 20 in more than 90% of samples (33*0.75 -> 25)
## suggested by WGCNA on RNAseq FAQ
dds90 <- dds[ rowSums(counts(dds) >= 20) >= round((length(df)-1)*0.90), ]
nrow(dds90)   

 

 

 
