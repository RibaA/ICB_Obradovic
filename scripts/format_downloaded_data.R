## load libraries

library(GEOquery)
library(data.table)
library(readxl)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
#work_dir <- args[1]
work_dir <- "data"

# CLIN.txt
#gse <- getGEO("GSE195832", GSEMatrix = TRUE, AnnotGPL = TRUE)

# Extract clinical/phenotype (metadata / sample info)
clin <- read.csv(file.path(work_dir, 'metadata_TJcohort.csv'))
clin <- as.data.frame(clin[, -1])
rownames(clin) <- clin$Sample
clin <- clin[order(clin$Sample), ]

write.table(clin, file=file.path(work_dir, 'CLIN.txt'), sep = "\t" , quote = FALSE , row.names = FALSE)

# EXP_TPM.tsv 
expr <- read_tsv( file = file.path(work_dir, "GSE195832_TJU_featurecounts_symbol.tsv") ) 
expr <- as.data.frame(expr)
expr <- expr[order(expr$symbol), ]
rownames(expr) <- expr$symbol
expr <- expr[, -1]
expr <- expr[, order(colnames(expr))]

# compute TPM
counts_to_tpm <- function(counts, gene_length) {

  gene_length_kb <- gene_length / 1000
  rpk <- counts / gene_length_kb
  scaling_factors <- colSums(rpk, na.rm = TRUE)
  tpm <- t(t(rpk) / scaling_factors) * 1e6
  return(tpm)

}

annot <- read_tsv( file = file.path(work_dir, "Human.GRCh38.p13.annot.tsv") ) 
annot <- annot[order(annot$Symbol), ]
int <- intersect(rownames(expr), annot$Symbol)

expr <- expr[int, ]
expr <- expr[order(rownames(expr)), ]
annot <- annot[annot$Symbol %in% int, ]

tpm <- counts_to_tpm(expr, annot$Length)

write.table(tpm, file=file.path(work_dir, 'EXP_TPM.tsv'), sep = "\t" , quote = FALSE , row.names = TRUE, col.names=TRUE)


