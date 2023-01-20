
# function to only install packages not installed yet. 
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}


packages <- c("BiocManager", 
              "readr",
              "splitstackshape",  
              "reshape2",
              "data.table",  
              "steepness",
              "PlayerRatings",
              "gridExtra",
              "ineq",
              "cluster",
              "factoextra",
              "Perc",
              "stringr",
              "brms",
              "coda",
              "broom",
              "multcomp",
              "bayesplot",
              "igraph",
              "ggthemes",
              "ggridges",
              "RColorBrewer",
              "viridis",
              "glue",
              "qwraps2",
              "lubridate",
              "cowplot",
              "tidybayes",
              "tidyverse",
              "mice",
              "PerformanceAnalytics",
              "ggsci",
              "msigdbr"
)
ipak(packages)
# Bioconductor packages require BioCManager 
BiocManager::install("enrichplot")
BiocManager::install("DESeq2")
BiocManager::install("VennDiagram")
BiocManager::install("genefilter")
BiocManager::install("AnnotationDbi")
BiocManager::install('EnhancedVolcano')
BiocManager::install("edgeR")
BiocManager::install("clusterProfiler")
BiocManager::install("DO.db")
BiocManager::install("STRINGdb")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")

# Packages not on CRAN or Bioconductor 
devtools::install_github("jalapic/compete")
devtools::install_github("stephenturner/annotables")
devtools::install_github("RRHO2/RRHO2")

