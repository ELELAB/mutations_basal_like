#Libaries -----------------------------------
library(tidyverse)

# Load Functions ----------------------------
source("/src/DMA.R")

# Load data ----------------------------------
load("data/BRCA_Basal_DEA_table.rda", verbose = TRUE) 

# Run Function --------------------------------
DMA(MafFile = "data/mutations.tsv",
    DEGs = DEA_table_extract,
    Drivers = "data/basal_specific_pra.RDS",
    coding_file = "data/css_coding.vcf.gz",
    noncoding_file = "data/css_noncoding.vcf.gz",
    results_folder = "results")

