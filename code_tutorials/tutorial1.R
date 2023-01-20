## https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html

if (!require("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("Glimma")
BiocManager::install("edgeR", force=TRUE)
BiocManager::install("Mus.musculus")

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


