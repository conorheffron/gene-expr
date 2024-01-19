library(DESeq2)            
library(survival)       
library(glmnet)         
library(ggplot2)      
if (!require("BiocManager", quietly = TRUE))
  BiocManager::install("clusterProfiler")
library(clusterProfiler)  
library(factoextra)       
library(AnnotationDbi)  
if (!require("BiocManager", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)    

clinical = data_clinical$data_clinical_patient
rnaseq = data_mrna$data_mrna_seq_v2_rsem
cna <- data_cna$data_cna                
keep <- !duplicated(rnaseq[,1])                     
rnaseq <- rnaseq[keep,]                           
rownames(rnaseq) <- rnaseq[,1]                     
erbb2_indx <- which(cna[,1] == 'ERBB2')            
rna_cna_id <- which(colnames(rnaseq[-c(1,2)]) %in% colnames(cna[-c(1,2)])) 
rna_cna_sub <- rnaseq[, 2 + rna_cna_id]            
countData <- round(as.matrix(rna_cna_sub))          
colnames(countData) <- colnames(rna_cna_sub)        
erbb2_cna_levels <- as.numeric(cna[erbb2_indx, colnames(countData)]) 
metadata <- data.frame(row.names = colnames(countData), 
                       erbb2_amplified = erbb2_cna_levels > 0)
metadata$erbb2_amplified <- as.factor(metadata$erbb2_amplified) 
dds <- DESeqDataSetFromMatrix(countData = countData, colData = DataFrame(metadata), design = ~ erbb2_amplified) 
dds <- DESeq(dds)                                
res <- results(dds, contrast = c("erbb2_amplified", "TRUE", "FALSE")) 
resOrdered <- res[order(res$padj), ]             
top10 <- head(resOrdered, 10)                      
resOrdered_up <- res[order(res$log2FoldChange, decreasing = TRUE),] 
top10_genes_up <- head(resOrdered_up, 10)                 
resOrdered_down <- res[order(res$log2FoldChange, decreasing = FALSE),] 
top10_genes_down <- head(resOrdered_down, 10)               
write.csv(top10, "./gene_mutation/mostsignificant.csv")
write.csv(top10_genes_up, "./gene_mutation/upregulated.csv")
write.csv(top10_genes_down, "./gene_mutation/downregulated.csv")
gene_symbols <- rownames(resOrdered)[which(resOrdered$padj < 0.05)] 
entrez_ids <- mapIds(org.Hs.eg.db,          
                     keys = gene_symbols,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")
entrez_ids <- entrez_ids[!is.na(entrez_ids)]
enrichResult <- enrichKEGG(gene = entrez_ids, organism = 'hsa', pvalueCutoff = 0.05)
write.csv(enrichResult, "./gene_mutation/PathwayEnrichmentResult.csv")
vsd <- vst(dds, blind = FALSE)                    
pcaResult <- prcomp(t(assay(vsd)))                
pcaData <- as.data.frame(pcaResult$x[, 1:2])      
sampleDists <- dist(t(assay(vsd)))                 
sampleClustering <- hclust(sampleDists, method = "average") 
clusters <- cutree(sampleClustering, k = 5)        
pcaData$cluster <- factor(clusters)                 
ggplot(pcaData, aes(x = PC1, y = PC2, color = cluster)) + geom_point() + theme_minimal() + ggtitle("PCA Plot with Clustering") 
dds$cluster <- factor(clusters)                 
design(dds) <- formula(~ cluster)                  
dds <- DESeq(dds)                                  
resClusterComparison <- results(dds, contrast=c("cluster", "1", "2")) 
top_genes_cluster <- head(resClusterComparison[order(resClusterComparison$pvalue), ])
write.csv(top_genes_cluster, "./gene_mutation/topgenescluster.csv")