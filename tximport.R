# tximport salmon TMP using ZhongmuNo1
rm(list = ls())
BiocManager::install("apeglm")

library(WGCNA)
library(tximport)
library(edgeR)
library(DESeq2)
library(flashClust)
library(csaw)
library(tidyr)
library(tidyverse)
library(data.table)
library(apeglm)

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

# wilson
data_rna_seq3 <- c("./LID50571_1023_quant.sf", "./LID50572_1001_quant.sf", "./LID50574_1008_quant.sf",
                   "./LID50574_1012_quant.sf", "./LID50574_1018_quant.sf", "./LID50575_1019_quant.sf")

# sanarac DS
data_rna_seq4 <- c("./LID50570_1005_quant.sf", "./LID50570_1005_quant.sf", "./LID50570_1008_quant.sf",
                   "./LID50570_1012_quant.sf", "./LID50571_1019_quant.sf", "./LID50571_1018_quant.sf")

txi.1 <- tximport(data_rna_seq1, type="salmon", txOut=TRUE, countsFromAbundance="scaledTPM")
txi.2 <- tximport(data_rna_seq2, type="salmon", txOut=TRUE, countsFromAbundance="scaledTPM")
txi.3 <- tximport(data_rna_seq3, type="salmon", txOut=TRUE, countsFromAbundance="scaledTPM")
txi.4 <- tximport(data_rna_seq4, type="salmon", txOut=TRUE, countsFromAbundance="scaledTPM")

cts <- txi.1$counts
cts[1:5,1:5]
dim(cts)
#############
cpms.2 <- as.data.frame(cts)
cpms.2 <- tibble::rownames_to_column(cpms.2, "lid")
cpms.2[1:5,1:5]
cpms.3 <- separate(cpms.2, 1, c('gene_id', 'rest'), sep = ";", remove = TRUE, convert = FALSE, extra = "warn")
cpms.3[1:5,1:5]
cpms.3 <- cpms.3[,-2]
cpms.3.1 <- cpms.3 %>% remove_rownames() %>% column_to_rownames(var = 'gene_id')
cpms.3.1[1:5,1:5]
colnames(cpms.3.1) <- list_4
dim(cpms.3.1)

cpms.3.2 <- tibble::rownames_to_column(cpms.3.1, "lid") %>% separate(1, c('gene_id', 'isoform'), sep = "\\.", remove = FALSE, convert = FALSE, extra = "warn")
cpms.3.2[1:5,1:5]
dim(cpms.3.2)
save(cpms.3.2, file = "~/Documents/Cesar/RNA/globus/lordec_reports/rna_seq/rna_seq_fastq/shen_sf/cpms.3.2.RData")

#############

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

setwd("~/Documents/Cesar/RNA/globus/lordec_reports/rna_seq/salmon_sf")
data_rna_seq1 = list.files(pattern = "quant.sf", full.names = T)

txi.1 <- tximport(data_rna_seq1, type="salmon", txOut=TRUE,
                  countsFromAbundance="scaledTPM")
txi.1$counts

#all samples
samples <- read.csv('~/Documents/Cesar/RNA/globus/lordec_reports/wgcna/sample_metadata3.csv')
str(samples)
list_4 <- colnames(samples)
samples[,list_4] <- lapply(samples[,list_4], factor)

dds <- DESeqDataSetFromTximport(txi.3, samples, ~ condition)

#pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
plotCounts(dds, gene = which.min(res$padj), intgroup = "condition")

summary(res)
sum(res$padj < 0.1, na.rm = T)
resultsNames(dds)
resOrdered <- res[order(res$pvalue), ]
resSig <- subset(resOrdered, padj < 0.1)
resSig.1 <- as.data.frame(resSig)


# log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef = "condition_DS_vs_CK", type = "apeglm")
summary(resLFC)
plotMA(resLFC, ylim = c(-2,2))

dds.1 <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")
dds.1@assays
resultsNames(dds)

######################

### Extract normalized expression for significant genes from the OE and control samples (4:9), and set the gene column (1) to row names
norm_OEsig <- normalized_counts[,c(1,4:9)] %>% 
  filter(gene %in% sigOE$gene) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 


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

##################
# edgeR
library(edgeR)

samples <- read.csv('~/Documents/Cesar/RNA/globus/lordec_reports/wgcna/sample_metadata3.csv')
str(samples)
list_4 <- colnames(samples)
samples[,list_4] <- lapply(samples[,list_4], factor)


cts <- txi.3$counts
normMat <- txi.3$length

# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

# Computing effective library sizes from scaled counts, to account for
# composition biases between samples.
eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Creating a DGEList object for use in edgeR.
y <- DGEList(cts, group = samples$condition)
y <- scaleOffset(y, normMat)
# filtering
keep <- filterByExpr(y)
table(keep)
y <- y[keep, ]
# y is now ready for estimate dispersion functions see edgeR User's Guide
dim(y)

# CPMs
# cpms <- edgeR::cpm(y, offset = y$offset, log = FALSE)

y <- calcNormFactors(y)
head(y$samples)

design <- model.matrix(~samples$condition)
y <- estimateDisp(y, design, robust = T)
plotBCV(y)
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)
qlf <- glmQLFTest(fit)
topTags(qlf, n = 15)
summary(decideTests(qlf))

et <- exactTest(y, pair = c("CK", "DS"))

topTags(et)
?glmQLFTest






top <- rownames(topTags(qlf))
cpm(y)[top,]


plotMD(qlf)
abline(h=c(-1,1), col="blue")

FDR <- p.adjust(qlf$table$PValue, method="BH")
sum(FDR < 0.05)

y1 <- estimateGLMCommonDisp(y, method = "deviance", robust = TRUE, subset = NULL)

qlf <- glmQLFTest(fit, coef=2:3)

mod <- model.matrix(~0 + samples$condition)
contrast.mat <- makeContrasts(Diff = CK - DS, levels = mod)

