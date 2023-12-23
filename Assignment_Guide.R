clinical = data_clinical$data_clinical_patient
rnaseq = data_mrna$data_mrna_seq_v2_rsem
keep = !duplicated(rnaseq[, 1])
rnaseq = rnaseq[keep, ]
rownames(rnaseq)  = rnaseq[, 1]
cna = data_cna$data_cna
erbb2_indx = which(cna[, 1] == 'ERBB2')
hist(as.numeric(cna[erbb2_indx, -c(1, 2)]))
rna_cna_id = which(is.element(colnames(rnaseq[, -c(1, 2)]), colnames(cna[, -c(1, 2)])))
rna_cna_sub = rnaseq[, 2 + rna_cna_id]
no_pats_in_rna_cna_sub_and_cna = sum(is.element(colnames(rnaseq[, 2 + rna_cna_id]), colnames(cna[, -c(1, 2)])))
sanity_check = no_pats_in_rna_cna_sub_and_cna == dim(rna_cna_sub)[2]
meta_erbb2 = matrix(0, length(rna_cna_id), 1)

for (i in 1:length(rna_cna_id)) {
  col_i = colnames(rna_cna_sub)[i]
  col_cna = which(colnames(cna) == col_i)
  meta_erbb2[i, ] = 1 * (cna[erbb2_indx, col_cna] > 0)
}

col_i = colnames(rna_cna_sub)[1]
col_cna = which(colnames(cna) == col_i)
(cna[erbb2_indx, col_cna] > 0) == meta_erbb2[1, 1]
pos_example = which(meta_erbb2 == 1)[1]
col_i = colnames(rna_cna_sub)[pos_example]
col_cna = which(colnames(cna) == col_i)
(cna[erbb2_indx, col_cna] > 0) == meta_erbb2[pos_example, 1]
colnames(meta_erbb2) = 'ERBB2Amp'
rna_cna_sub = round(rna_cna_sub)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("BiocManager", quietly = TRUE))
  BiocManager::install("DESeq2")