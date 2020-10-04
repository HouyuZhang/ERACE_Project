# Import & pre-process ----------------------------------------------------
workDir <- "E:/Lu_RNA-seq/DESeq"
setwd(workDir)
sampleFile <- read.table("E:/Lu_RNA-seq/read_count/featureCounts/featureCounts_MA_simpliy.txt",header = T,row.names=1)
#rRNA <- read.table("E:/Lu_RNA-seq/read_count/featureCounts/rRNA.txt",header = F)
head(sampleFile)
# Convert to matrix
countdata <- as.matrix(sampleFile)
head(countdata)
condition <- factor(c(rep(c("MA9_Con0","MA9_Con12","MA9_Con24","MA9_Con3","MA9_Con6","MA9_Dox12","MA9_Dox24","MA9_Dox3","MA9_Dox6"),2)))
#condition <- factor(c(rep(c("WT_Con0","WT_Con12","WT_Con24","WT_Con3","WT_Con6","WT_Dox12","WT_Dox24","WT_Dox3","WT_Dox6"),2)))
# Analysis with DESeq2 ----------------------------------------------------
library(DESeq2)
coldata <- data.frame(row.names=colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- dds[rowSums(counts(dds)) > 1,]
dds
# Run the DESeq pipeline
dds <- DESeq(dds)
# Get differential expression results

res <- results(dds, alpha=0.05,contrast = c("condition","MA9_Dox3","MA9_Con0"))
sum(res$padj < 0.05, na.rm=TRUE)
summary(res)
## Order by adjusted p-value
resOrdered <- res[order(res$padj),]
## Write results
write.csv(as.data.frame(resOrdered),file="re_nonmRNA_MA9_Dox3vsMA9_Con0_results.csv")

#Visulization------------------------------------------------------------------
plotMA(dds, ylim=c(-1,1), cex=1)

#plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition",
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  scale_y_log10(breaks=c(25,100,400))

# Plot dispersions
plotDispEsts(dds, main="Dispersion plot")
# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd.fast <- vst(dds, blind=FALSE)
head(assay(rld))
head(assay(vsd))
head(assay(vsd.fast))
hist(assay(rld))

# Colors for plots below
library(RColorBrewer)

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(vsd))))
library(ComplexHeatmap)
Heatmap(as.matrix(sampleDists),
                  row_names_side = "left")

# Principal components analysis (PCA)
data <- plotPCA(vsd, intgroup=c("condition", "strain"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=strain)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()


## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)



## Examine independent filtering
attr(res, "filterThreshold")
plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")


## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
volcanoplot(res, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))




