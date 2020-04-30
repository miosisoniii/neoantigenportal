#Functions.R
library(tidyverse)
library(shinyFiles)
library(doMC)
registerDoMC(20)

#total_data <- read.csv("totaldata.csv", check.names = TRUE)

#Use TP53 Gene Frequency
genes <- c("TP53", "MYCN", "MLN", "BRAF", "PI3K")
gentab <- data.frame(matrix(ncol = 3, nrow=length(genes)))
colnames(gentab) <- c("gene","unitprot","seq")
gentab$gene <- genes

#Set Sequences
gene_seq_df <- read.csv("gene_seq.csv")

#Read in frequency file for HLA alleles
hla <- read.delim('HLAfreq.txt',header=T, stringsAsFactors=F)
hla_col <- read.delim('HLAfreq.txt',header=T, stringsAsFactors=F, check.names = F) 
hlanames <- colnames(hla_col[,-1])

