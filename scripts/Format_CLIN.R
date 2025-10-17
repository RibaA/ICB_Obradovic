## load libraries
library(stringr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
#input_dir <- args[1]
#output_dir <- args[2]
#annot_dir <- args[3]

input_dir <- "data/input"
output_dir <- "data/output"
annot_dir <- "data/annot"

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_tissue.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_drug.R")

clin = read.csv( file.path(input_dir, "CLIN.txt"), stringsAsFactors=FALSE , sep="\t" , header=TRUE )
rownames(clin) <- clin$Sample
cols <- c('Sample', 'Gender', 'Age')

clin <- clin[, c(cols, colnames(clin)[!colnames(clin) %in% cols])]
colnames(clin)[colnames(clin) %in% cols] <- c('patient', 'sex', 'age')

new_cols <- c( "primary" , "pfs", "response", "response.other.info" ,  'histo', 'stage', "recist", "drug_type",
                                              'os', 't.os', 't.pfs', "dna" , "rna", "rna.info")
clin[new_cols] <- NA

clin$rna <- 'rnaseq'
clin$rna.info <- 'tpm'
clin$primary <- 'HNSCC'
clin$response.other.info <- ifelse(clin$'Clinical.Response.Primary_Site_TE.' %in% c('Responder', 'Complete Responder'), 'R', 'NR')
clin$drug_type <- ifelse(clin$Therapy == 'antiPD1', 'PD-1/PD-L1', 'PD-1+Tadalafil')
clin$response = Get_Response( data=clin )

# Tissue and drug annotation
annotation_tissue <- read.csv(file=file.path(annot_dir, 'curation_tissue.csv'))
clin <- annotate_tissue(clin=clin, study='Obradovic', annotation_tissue=annotation_tissue, check_histo=FALSE)

annotation_drug <- read.csv(file=file.path(annot_dir, 'curation_drug.csv'))
clin <- add_column(clin, treatmentid=clin$drug_type, .after='tissueid')

write.table( clin , file=file.path(output_dir, "CLIN.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

