# Assignment Script

# Luis Iglesias-Martinez PhD, 15/11/2023
# This script is to be used as a guide for your assignment
# It will help you to read and split the data.
# 

# Step 1. Download Data Into your Computer

# We will use the TCGA- Invasive Breast Carcinoma PanCancer Atlas from TCGA in cbioportal
# 

# We are using the tar file on your directory.

# Change working directory.

path  = path_wd 
# change this to your own directory

file_name = "brca_tcga_pan_can_atlas_2018.tar.gz"

# We will first extract the files into folders.

untar(file_name)

# change directory to the extracted folders

setwd(paste(path, "/brca_tcga_pan_can_atlas_2018", sep = ""))

# We will use the following files:

# data_clinical_patient.txt, data_mrna_seq_v2_rsem.txt, data_mutations.txt and data_cna.txt

clinical = read.delim("data_clinical_patient.txt")

rnaseq = read.delim("data_mrna_seq_v2_rsem.txt")

# in this assignment we will delete the genes for which there's more than one Hugo Symbol
# These are typically genes with no Hugo Symbol ("" as an entry) or pseudogenes.

# This is more for simplicity.If you keep your analysis would still be correct so no worries.

keep = !duplicated(rnaseq[,1])

rnaseq = rnaseq[keep,]

# set rownames of rnaseq to hugo symbols

rownames(rnaseq)  = rnaseq[,1]

# Read CNA Data

cna = read.delim('data_cna.txt')

# find ERBB2 in cna

erbb2_indx = which(cna[,1] == 'ERBB2')

# Plot histogram to visualize explore the data.

hist(as.numeric(cna[erbb2_indx,-c(1,2)]))

# match patients in rnaseq to patients in cna.

rna_cna_id = which(is.element(colnames(rnaseq[,-c(1,2)]), colnames(cna[,-c(1,2)])))

# select only the rna cases which have cna data.

rna_cna_sub = rnaseq[,2+rna_cna_id]

# check all patients in rna_can_sub are in cna



no_pats_in_rna_cna_sub_and_cna = sum(is.element(colnames(rnaseq[,2+rna_cna_id]), colnames(cna[,-c(1,2)]))) 



# sanity check.This will print an error if the result is not the same.

sanity_check = no_pats_in_rna_cna_sub_and_cna == dim(rna_cna_sub)[2]

# Pre-allocate memory for ERBB2

meta_erbb2 = matrix(0,length(rna_cna_id),1)

for (i in 1:length(rna_cna_id)){
  # access the colnames of i
  col_i = colnames(rna_cna_sub)[i]
  # get the index in cna for the same patient
  col_cna = which(colnames(cna)==col_i)
  # store if they're amplified.
  meta_erbb2[i,] = 1*(cna[erbb2_indx,col_cna]>0)
  
}

# This are some checks you can do to make sure your code worked.
# There's some more systematic checks you can do. See unit testing.


# simple checks to make sure. 

col_i = colnames(rna_cna_sub)[1]

col_cna = which(colnames(cna)==col_i)

# sanity check

(cna[erbb2_indx,col_cna]>0) == meta_erbb2[1,1]

# see now if a positive meta_erbb2 is amplified.

pos_example = which(meta_erbb2==1)[1]


col_i = colnames(rna_cna_sub)[pos_example]

col_cna = which(colnames(cna)==col_i)

# sanity check

(cna[erbb2_indx,col_cna]>0) == meta_erbb2[pos_example,1]

# botch checks should print true.

# We will add a title to the metadata.

colnames(meta_erbb2) = 'ERBB2Amp'

# transform into integers

rna_cna_sub = round(rna_cna_sub)

# Install DESeq2.
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("BiocManager", quietly = TRUE))
  BiocManager::install("DESeq2")