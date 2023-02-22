## https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html

if (!require("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

#BiocManager::install("limma")
#BiocManager::install("Glimma")
#BiocManager::install("edgeR", force=TRUE)
#BiocManager::install("Mus.musculus")

library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
library(R.utils)


#####
# Import Example Data 
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb") 
utils::untar("GSE63310_RAW.tar", exdir = ".")
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
           "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
           "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")
for(i in paste(files, ".gz", sep=""))
  R.utils::gunzip(i, overwrite=TRUE)

#####
# Each of these text files contains the raw gene-level counts for a given sample. 
# Note analysis only includes the basal, LP and ML samples from this experiment 

files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", 
           "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt", 
           "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt", 
           "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", 
           "GSM1545545_JMS9-P8c.txt")
read.delim(files[1], nrow=5)



# Whilst each of the nine text files can be read into R separately and combined into a matrix
# of counts, edgeR offers a convenient way to do this in one step using the readDGE function.
# The resulting DGEList-object contains a matrix of counts with 27,179 rows associated 
# with unique Entrez gene identifiers (IDs) and nine columns associated with the individual 
# samples in the experiment.

x <- readDGE(files, columns=c(1,3))
class(x)
dim(x)


#' For downstream analysis, sample-level information related to the experimental design 
#' needs to be associated with the columns of the counts matrix. This should include
#'  experimental variables, both biological and technical, that could have an effect on 
#'  expression levels. Examples include cell type (basal, LP and ML in this experiment), 
#'  genotype (wild-type, knock-out), phenotype (disease status, sex, age), 
#'  sample treatment (drug, control) and batch information (date experiment was 
#'  performed if samples were collected and analysed at distinct time points) 
#'  to name just a few.


#' Our DGEList-object contains a samples data frame that stores both cell type 
#' (or group) and batch (sequencing lane) information, each of which consists of three
#'  distinct levels. Note that within x$samples, library sizes are automatically 
#'  calculated for each sample and normalisation factors are set to 1. 
#'  For simplicity, we remove the GEO sample IDs (GSM*) from the column names of our 
#'  DGEList-object x.




samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
samplenames

colnames(x) <- samplenames
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", 
                     "Basal", "ML", "LP"))
x$samples$group <- group
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
x$samples$lane <- lane
x$samples


##### Data Packaging

## Organising gene annotations


#' A second data frame named genes in the DGEList-object is used to store gene-level 
#' information associated with rows of the counts matrix. This information can be retrieved
#'  using organism specific packages such as Mus.musculus (Bioconductor Core Team 2016b) 
#'  for mouse (or Homo.sapiens (Bioconductor Core Team 2016a) for human) or the
#'   biomaRt package (Durinck et al. 2005, 2009) which interfaces the Ensembl genome
#'    databases in order to perform gene annotation.

#' The type of information that can be retrieved includes gene symbols, gene names, 
#' chromosome names and locations, Entrez gene IDs, Refseq gene IDs and Ensembl gene
#'  IDs to name just a few. biomaRt primarily works off Ensembl gene IDs, 
#'  whereas Mus.musculus packages information from various sources and allows users 
#'  to choose between many different gene IDs as the key.

#' The Entrez gene IDs available in our dataset were annotated using the Mus.musculus 
#' package to retrieve associated gene symbols and chromosome information.


geneid <- rownames(x)
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENTREZID")
head(genes)


# check for duplicated genes, e.g.if on multiple chromosomes,
# could keep eg first instance or combine into one entry

#keeping unduplicated
genes <- genes[!duplicated(genes$ENTREZID),]


#' The data frame of gene annotations is then added to the data object and 
#' neatly packaged in a DGEList-object containing raw count data with associated 
#' sample information and gene annotations.

# should use match() to do this probably, as unlikely to always have genes in same order
# as data object

x$genes <- genes
x


##### Data Pre-Processing

#' Transformations from the raw-scale
#' 
#' 
#' For differential expression and related analyses, gene expression is rarely considered 
#' at the level of raw counts since libraries sequenced at a greater depth will result in higher counts. 
#' Rather, it is common practice to transform raw counts onto a scale that accounts for such library size 
#' differences. Popular transformations include counts per million (CPM), log2-counts per million (log-CPM), 
#' reads per kilobase of transcript per million (RPKM), and fragments per kilobase of transcript per
#'  million (FPKM).
#'  
#'  In our analyses, CPM and log-CPM transformations are used regularly although 
#'  they do not account for gene length differences as RPKM and FPKM values do. Whilst RPKM and FPKM 
#'  values can just as well be used, CPM and log-CPM values can be calculated using 
#'  a counts matrix alone and will suffice for the type of comparisons we are interested in. 
#'  
#'  Assuming that there are no differences in isoform usage between conditions, differential expression analyses 
#'  look at gene expression changes between conditions rather than comparing expression across 
#'  multiple genes or drawing conclusions on absolute levels of expression. In other words, 
#'  gene lengths remain constant for comparisons of interest and any observed differences are
#'   a result of changes in condition rather than changes in gene length.

#' Here raw counts are converted to CPM and log-CPM values using the cpm function in edgeR. 
#' RPKM values are just as easily calculated as CPM values using the rpkm function in edgeR 
#' if gene lengths are available.

cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)


#' A CPM value of 1 for a gene equates to having 20 counts in the sample with the lowest sequencing depth
#'  (JMS0-P8c, library size approx. 20 million) or 76 counts in the sample with the greatest sequencing depth 
#'  (JMS8-3, library size approx. 76 million).

#' The log-CPM values will be used for exploratory plots. When log=TRUE, the cpm function adds an offset 
#' to the CPM values before converting to the log2-scale. By default, the offset is 2/L where 2 is the “prior count” and L 
#' is the average library size in millions, so the log-CPM values are related to the CPM values by log2(CPM + 2/L).
#'  This calculation ensures that any two read counts with identical CPM values will also have identical log-CPM values. 
#'  The prior count avoids taking the logarithm of zero, and also reduces spurious variability for genes with very low 
#'  counts by shrinking all the inter-sample log-fold-changes towards zero, something that is helpful for exploratory plotting. 
#'  For this dataset, the average library size is about 45.5 million, so L approx. 45.5 and the minimum log-CPM value for
#'   each sample becomes log2(2/45.5) = -4.51. In other words, a counr of zero for this data maps to a log-CPM value of -4.51 
#'   after adding the prior count or offset:


L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)

summary(lcpm)

#' Log-CPM values are also used in downstream linear modeling via limma’s voom function, although voom recomputes 
#' its own log-CPM values internally with a smaller prior count. 



### Removing genes that are lowly expressed

#' All datasets will include a mix of genes that are expressed and those that are not expressed. 
#' Whilst it is of interest to examine genes that are expressed in one condition but not in another,
#'  some genes are unexpressed throughout all samples. In fact, 19% of genes in this dataset have zero 
#'  counts across all nine samples.


table(rowSums(x$counts==0)==9)


#' Plotting the distribution log-CPM values shows that a sizeable proportion of genes within each sample 
#' are either unexpressed or lowly-expressed with log-CPM values that are small or negative (Figure 1A).

#' Genes that do not have a worthwhile number of reads in any sample should be filtered out of the downstream analyses.
#'  There are several reasons for this. From a biological point of view, genes that not expressed at a
#'   biologically meaningful level in any condition are not of interest and are therefore best ignored. 
#'   From a statistical point of view, removing low count genes allows the mean-variance relationship in the data
#'    to be estimated with greater reliability and also reduces the number of statistical tests that need to
#'     be carried out in downstream analyses looking at differential expression.


#' The filterByExpr function in the edgeR package provides an automatic way to filter genes, 
#' while keeping as many genes as possible with worthwhile counts.

dim(x) #[1] 27179     9

keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

#' By default, the function keeps genes with about 10 read counts or more in a minimum number of samples, 
#' where the number of samples is chosen according to the minimum group sample size. 
#' The actual filtering uses CPM values rather than counts in order to avoid giving preference to samples with 
#' large library sizes. For this dataset, the median library size is about 51 million and 10/51 approx. 0.2, so
#'  the filterByExpr function keeps genes that have a CPM of 0.2 or more in at least three samples. 
#'  A biologically interesting gene should be expressed in at least three samples because all the cell type 
#'  groups have three replicates. The cutoffs used depend on the sequencing depth and on the experimental design. 
#'  If the library sizes had been larger then a lower CPM cutoff would have been chosen, because larger library sizes 
#'  provide better resolution to explore more genes at lower expression levels. Alternatively, smaller library sizes
#'   decrease our ability to explore marginal genes and hence would have led to a higher CPM cutoff.

#' Using this criterion, the number of genes is reduced to 16,624, about 60% of the number that we started with 
#' (panel B of the next figure). Note that subsetting the entire DGEList-object removes both the counts and 
#' the associated gene information for the filtered genes. The filtered DGEList-object keeps the gene 
#' information and the counts for the retained genes correctly associated.



## Make Plot of Distributions

lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")




### Normalising gene expression distributions

#' During the sample preparation or sequencing process, external factors that are not of biological interest
#'  can affect the expression of individual samples. For example, samples processed in the first batch of an 
#'  experiment can have higher expression overall when compared to samples processed in a second batch. It is 
#'  assumed that all samples should have a similar range and distribution of expression values. Normalisation 
#'  is required to ensure that the expression distributions of each sample are similar across the entire experiment.

#' Any plot showing the per sample expression distributions, such as a density or boxplot, is useful in determining
#'  whether any samples are dissimilar to others. Distributions of log-CPM values are similar throughout all
#'   samples within this dataset (panel B of the figure above).

#' Nonetheless, normalisation by the method of trimmed mean of M-values (TMM) (Robinson and Oshlack 2010) is 
#' performed using the calcNormFactors function in edgeR. The normalisation factors calculated here are used 
#' as a scaling factor for the library sizes. When working with DGEList-objects, these normalisation factors 
#' are automatically stored in x$samples$norm.factors. For this dataset the effect of TMM-normalisation is mild,
#'  as evident in the magnitude of the scaling factors, which are all relatively close to 1.


x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

#' To give a better visual representation of the effects of normalisation, the data was duplicated then
#'  adjusted so that the counts of the first sample are reduced to 5% of 
#' their original values, and in the second sample they are inflated to be 5-times larger.
#' 
#' The figure below shows the expression distribution of samples for unnormalised and normalised data,
#'  where distributions are noticeably different pre-normalisation and are similar post-normalisation. 
#'  Here the first sample has a small TMM scaling factor of 0.06, whereas the second sample has a large 
#'  scaling factor of 6.08 – neither values are close to 1.



x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors


lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")



### Unsupervised clustering of samples

#' In our opinion, one of the most important exploratory plots to examine for gene expression analyses
#'  is the multi-dimensional scaling (MDS) plot, or similar. The plot shows similarities and dissimilarities
#'   between samples in an unsupervised manner so that one can have an idea of the extent to which 
#'   differential expression can be detected before carrying out formal tests. Ideally, samples would cluster
#'    well within the primary condition of interest, and any sample straying far from its group could be 
#'    identified and followed up for sources of error or extra variation. If present, technical replicates 
#'    should lie very close to one another.

#' Such a plot can be made in limma using the plotMDS function. The first dimension represents the leading-fold-change
#'  that best separates samples and explains the largest proportion of variation in the data, with subsequent 
#'  dimensions having a smaller effect and being orthogonal to the ones before it. When experimental design 
#'  involves multiple factors, it is recommended that each factor is examined over several dimensions. If samples cluster by
#'   a given factor in any of these dimensions, it suggests that the factor contributes to expression differences and is worth
#'    including in the linear modelling. On the other hand, factors that show little or no effect may be left
#'     out of downstream analysis.


#' In this dataset, samples can be seen to cluster well within experimental groups over dimension 1 and 2, 
#' and then separate by sequencing lane (sample batch) over dimension 3 (shown in the plot below).
#'  Keeping in mind that the first dimension explains the largest proportion of variation in the
#'   data, notice that the range of values over the dimensions become smaller as we move to higher dimensions.

#' Whilst all samples cluster by groups, the largest transcriptional difference is observed between basal 
#' and LP, and basal and ML over dimension 1. For this reason, it is expected that pairwise comparisons between
#'  cell populations will result in a greater number of DE genes for comparisons involving basal samples, and 
#'  relatively small numbers of DE genes when comparing ML to LP. Datasets where samples do not cluster by
#'   experimental group may show little or no evidence of differential expression in the downstream analysis.

#' To create the MDS plots, we assign different colours to the factors of interest. 
#' Dimensions 1 and 2 are examined using the color grouping defined by cell types.

#' Dimensions 3 and 4 are examined using the colour grouping defined by sequencing lanes (batch).

lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")


#' Alternatively, the Glimma package offers the convenience of an interactive MDS plot where multiple 
#' dimensions can be explored. The glMDSPlot function generates an html page (that is opened in a browser 
#' if launch=TRUE) with an MDS plot in the left panel and a barplot showing the proportion of variation explained
#'  by each dimension in the right panel. Clicking on the bars of the bar plot changes the pair of dimensions plotted
#'   in the MDS plot, and hovering over the individual points reveals the sample label. The colour scheme can 
#'   be changed as well to highlight cell population or sequencing lane (batch). An interactive MDS plot of 
#'   this dataset can be found at http://bioinf.wehi.edu.au/folders/limmaWorkflow/glimma-plots/MDS-Plot.html.


glMDSPlot(lcpm, labels=paste(group, lane, sep="_"), 
          groups=x$samples[,c(2,5)], launch=TRUE)




### Creating a design matrix and contrasts


# design matrix is set up with both the cell population and sequencing lane (batch) information
design <- model.matrix(~0+group+lane)
colnames(design) <- gsub("group", "", colnames(design))
design

#  Contrasts for pairwise comparisons between cell populations are set up in limma using the makeContrasts function.

contr.matrix <- makeContrasts(
  BasalvsLP = Basal-LP, 
  BasalvsML = Basal - ML, 
  LPvsML = LP - ML, 
  levels = colnames(design))
contr.matrix

 #  voom-plot provides a visual check on the level of filtering performed upstream. 

par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

# Examining the number of DE genes
summary(decideTests(efit))


#' The treat method (McCarthy and Smyth 2009) can be used to calculate p-values from 
#' empirical Bayes moderated t-statistics with a minimum log-FC requirement

# log fold change > 1
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

#' Genes that are DE in multiple comparisons can be extracted using the results from decideTests, 
#' where 0s represent genes that are not DE, 1s represent genes that are up-regulated, 
#' and -1s represent genes that are down-regulated. A total of 2,784 genes are DE in both 
#' basal versus LP and basal versus ML, twenty of which are listed below. 
#' 
#' The write.fit function can be used to extract and write results for all three comparisons
#'  to a single output file.

de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)


head(tfit$genes$SYMBOL[de.common], n=20)

dev.off()
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))


#' Examining individual DE genes from top to bottom
#' 

#' The top DE genes can be listed using topTreat for results using treat (or topTable for
#'  results using eBayes). By default topTreat arranges genes from smallest to largest
#'   adjusted p-value with associated gene information, log-FC, average log-CPM, moderated 
#'   t-statistic, raw and adjusted p-value for each gene. The number of top genes displayed
#'    can be specified, where n=Inf includes all genes. Genes Cldn7 and Rasef are amongst the
#'     top DE genes for both basal versus LP and basal versus ML.
#' 

basal.vs.lp <- topTreat(tfit, coef=1, n=Inf)
basal.vs.ml <- topTreat(tfit, coef=2, n=Inf)
head(basal.vs.lp)
head(basal.vs.ml)


## Graphical representations

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))


#' A heatmap is created for the top 100 DE genes (as ranked by adjusted p-value) 
#' from the basal versus LP contrast using the heatmap.2 function from the gplots package. 
#' The heatmap correctly clusters samples by cell type and reorders the genes into blocks
#'  with similar expression patterns. From the heatmap, we observe that the expression of ML
#'   and LP samples are very similar for the top 100 DE genes between basal and LP.

library(gplots)
basal.vs.lp.topgenes <- basal.vs.lp$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row",
          labRow=v$genes$SYMBOL[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")
