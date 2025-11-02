getwd()
setwd('C:/Users/chaud/Desktop/biosetup/preprocess')

#install.packages("BiocManager")
#BiocManager::install('DESeq2')
#BiocManager::install('pheatmap')
#BiocManager::install('RColorBrewer')

count <- as.matrix(read.csv("count_data.csv", sep = ",", row.names = "gene_id"))
head(count)
dim(count)
colnames(count)

meta <- read.csv("metadata.csv", row.names = 1)
head(meta)

all(rownames(meta)==colnames(count))

library('DESeq2')
dds <- DESeqDataSetFromMatrix(countData = count, colData = meta, design = ~condition)
dim(dds)

keep <- rowSums(counts(dds))>=10
dds <- dds[keep,]
dim(dds)

#DESeq2 
dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)
head(assay(vsd),3)

#heatmap
library('pheatmap')
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing = TRUE)[1:20]
df <- as.data.frame(colData(dds))

# 1. Define the colors for the 'condition' variable.
#    You must provide one color for each level (e.g., 'treated' and 'untreated').
ann_colors = list(
  condition = c(treated = "orange", untreated = "purple")
)

# 2. Run pheatmap with the new annotation_colors argument.
pheatmap(assay(vsd)[select,], 
         cluster_rows = FALSE,
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df,
         annotation_colors = ann_colors)

library("RColorBrewer")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix)
colnames(sampleDistMatrix)
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(225)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

#PCA plot
plotPCA(vsd, intgroup='condition')

#resukt of Differetial gene expression analysis
res <- results(dds, contrast = c("condition","treated","untreated"),alpha=0.05)
summary(res)
write.csv(res, "DEG.csv")

#create MA plot
sum(res$padj <0.05, na.rm=TRUE)
plotMA(res, ylim=c(-2,2))
