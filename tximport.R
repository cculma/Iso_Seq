# tximport salmon TMP using ZhongmuNo1
rm(list = ls())

library(WGCNA)
library(tximport)
library(edgeR)
library(DESeq2)
library(flashClust)
library(csaw)
library(tidyr)
library(tidyverse)
library(data.table)

#################
list_2 <- c("LID50568_1001", "LID50568_1006", "LID50568_1020",
            "LID50569_1002", "LID50569_1003", "LID50569_1004",
            "LID50570_1005", "LID50570_1008", "LID50570_1012",
            "LID50571_1018", "LID50571_1019", "LID50571_1023",
            "LID50572_1001", "LID50572_1003", "LID50573_1004",
            "LID50573_1005", "LID50573_1006", "LID50574_1008",
            "LID50574_1012", "LID50574_1018", "LID50575_1019")
list_3 <- c('PI467895-CK-Stem','Saranac-DS-Stem','Saranac-SS-Leaf','PI467895-CK-Leaf',
            'PI467895-SS-Stem','PI467895-SS-Leaf','Saranac-DS-Leaf','Saranac-DS-Root',
            'Saranac-CK-Leaf','Saranac-CK-Root','Saranac-CK-Stem','Wilson-DS-Leaf',
            'Wilson-DS-Stem','PI467895-CK-Root','PI467895-SS-Root','Saranac-SS-Root',
            'Saranac-SS-Stem','Wilson-DS-Root','Wilson-CK-Leaf','Wilson-CK-Stem',
            'Wilson-CK-Root')

list_4 <- c('Saranac-DS-Leaf','Saranac-DS-Stem','Saranac-DS-Root',
            'Saranac-CK-Leaf','Saranac-CK-Stem','Saranac-CK-Root',
            'Wilson-DS-Leaf','Wilson-DS-Stem','Wilson-DS-Root',
            'Wilson-CK-Leaf','Wilson-CK-Stem','Wilson-CK-Root')

#################
setwd("/home/hawkins/Documents/Cesar/RNA/globus/lordec_reports/rna_seq/rna_seq_fastq/shen_sf")
data_rna_seq1 = list.files(pattern = "quant.sf", full.names = T)
data_rna_seq2 <- c("./LID50570_1005_quant.sf", "./LID50570_1005_quant.sf", "./LID50570_1008_quant.sf",
                   "./LID50570_1012_quant.sf", "./LID50571_1019_quant.sf", "./LID50571_1018_quant.sf",
                   "./LID50571_1023_quant.sf", "./LID50572_1001_quant.sf", "./LID50574_1008_quant.sf",
                   "./LID50574_1012_quant.sf", "./LID50574_1018_quant.sf", "./LID50575_1019_quant.sf")

txi.1 <- tximport(data_rna_seq1, type="salmon", txOut=TRUE, countsFromAbundance="scaledTPM")
txi.2 <- tximport(data_rna_seq2, type="salmon", txOut=TRUE, countsFromAbundance="scaledTPM")
cts <- txi.2$counts
cts[1:5,1:5]
normMat <- txi.1$length
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat
eff.lib <- calcNormFactors(normCts) * colSums(normCts)

normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)
class(normMat)
dim(normMat)
normMat[1:5,1:5]
# Creating a DGEList object for use in edgeR.
y <- DGEList(cts)
y <- scaleOffset(y, normMat)
# filtering
keep <- filterByExpr(y)
y <- y[keep, ]
class(y)
se <- SummarizedExperiment(assays = list(counts = y$counts, offset = y$offset))
se$totals <- y$samples$lib.size

cpms <- calculateCPM(se, use.offsets = TRUE, log = FALSE)
dim(cpms)
cpms[1:5,1:5]
colnames(cpms) <- list_3
thres <- cpms > 0.5
keep <- rowSums(thres) >= 2 



cpms.2 <- as.data.frame(cpms)
cpms.2 <- tibble::rownames_to_column(cpms.2, "lid")
cpms.2[1:5,1:5]
cpms.3 <- separate(cpms.2, 1, c('gene_id', 'rest'), sep = ";", remove = TRUE, convert = FALSE, extra = "warn")
cpms.3[1:5,1:5]
cpms.3 <- cpms.3[,-2]
head(cpms.3)
cpms.3.1 <- cpms.3 %>% remove_rownames() %>% column_to_rownames(var = 'gene_id')
cpms.3.1[1:5,1:5]
colnames(cpms.3.1)

cpms.4 <- cpms.3.1[,c(2,7:13, 18:21)]
colnames(cpms.4)
cpms.4 <- cpms.4[,c(2,1,3,4,6,5,7:12)]
cpms.4[1:5,1:5]
dim(cpms.4)


tf_expr <- inner_join(cpms.3, all_sqanti_3, by = "isoform")



write.table(cpms.3, '~/Documents/Cesar/RNA/globus/lordec_reports/wgcna/cpms.3.tsv', sep="\t", quote = FALSE, col.names = T, row.names = F)
cpms.1 <- t(cpms.3.1)
cpms.1[1:5,1:5]
write.table(cpms.1, 'cpms.3T.tsv', sep="\t", quote = FALSE, col.names = T, row.names = F)
#colnames(cpms) <- list_3

cpms.1 <- t(cpms)
cpms.1[1:4, 1:4]
dim(cpms.1)

# save(cpms.1, file = "cpms.1.RData")

######################
# DegSeq2

samples
samples_1 <- unite(samples, col = condition_1, variety, condition, sep = "_") 

sampleTable <- data.frame(condition = factor(rep(c("A", "B"), each = 3)))
rownames(sampleTable) <- colnames(txi$counts)

setwd("~/Documents/Cesar/RNA/globus/lordec_reports/rna_seq/salmon_sf")
data_rna_seq1 = list.files(pattern = "quant.sf", full.names = T)

txi.1 <- tximport(data_rna_seq1, type="salmon", txOut=TRUE,
                  countsFromAbundance="scaledTPM")
txi.1$counts

#all samples
samples <- read.csv('~/Documents/Cesar/RNA/globus/lordec_reports/wgcna/sample_metadata2.csv')
samples1 <- samples %>% remove_rownames() %>% column_to_rownames(var = 'sample_id')
head(samples1)
# head(datTraits)
# datTraits1 <- datTraits[,-c(1:9)]
# head(datTraits1)
str(samples1)
list_4 <- colnames(samples1)
samples1[,list_4] <- lapply(samples1[,list_4], factor)

dds <- DESeqDataSetFromTximport(txi.2, samples1, ~ variety + condition)

#pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
summary(res)

# save.image("~/Documents/Cesar/RNA/globus/lordec_reports/rna_seq/rna_seq_fastq/shen_sf/tximport.RData")

dds.1 <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")
dds.1@assays
resultsNames(dds)

# Remove all rows with less than n counts across all samples, where n=#samples
raw_counts <- txi.1$counts
dim(raw_counts)
raw_counts1 <- t(raw_counts)
low_count_mask <- rowSums(raw_counts) < ncol(raw_counts)
dim(low_count_mask)
sprintf("Removing %d low-count genes (%d remaining).", sum(low_count_mask), 
        sum(!low_count_mask))

log_counts <- log2(raw_counts + 1)
sft = pickSoftThreshold (cpms.1, powerVector = powers, verbose = 5)

log_counts[1:5, 1:21]
dim(log_counts)
gsg = goodSamplesGenes(log_counts, verbose = 3);
gsg$allOK
gsg$goodSamples

options(stringsAsFactors = FALSE);
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(log_counts, powerVector = powers, verbose = 5)

sizeGrWindow(9,5);
par(mfrow= c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab = "Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.90, cool="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

######################

head(cpms.5)
library(pheatmap)
pheatmap(cpms.5)



cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(cpms.5, 1, cal_z_score))

pheatmap(data_subset_norm, drop_levels = TRUE,
         show_rownames=F, cluster_cols=F, cluster_rows=T, scale="row", 
         cex=1, clustering_distance_rows="euclidean", cex=1,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE)

data_subset_norm[1:4, 1:4]

pathways <- read.csv("modules_pathways.csv", header = T)

head(pathways)
dim(pathways)
colnames(pathways)
pathways.1 <- separate(pathways, 12, c('Pathways', 'ID'), sep = "\\.", remove = TRUE, convert = FALSE, extra = "merge")
head(pathways.1)
colnames(pathways.1)
pathways.2 <- separate(pathways.1, 12, c('Pathways_1', 'Pathways_2', 'Pathways_3'), sep = ";", remove = TRUE, convert = FALSE, extra = "warn")
cc1 <- count(pathways.2, Pathways_1, Pathways_2)
cc <- count(pathways.2, module, Pathways_1, Pathways_2)
cc2 <- count(pathways.2, treatment, Pathways_1, Pathways_2)
