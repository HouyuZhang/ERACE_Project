library(DESeq2)
directory <- "E:/Dropbox/Luting/Lab/data analysis/RNA-seq data analysis/new analysis/MA9_count/"
setwd("E:/Dropbox/Luting/Lab/data analysis/RNA-seq data analysis/new analysis/MA9_count/")

#### can merge individual sample files (i.e. ctrl1.counts, ctrl2.counts, etc.)####
sampleFiles <- grep("MA9",list.files("."),value=T)
sampleNames <- sub("_gene","h_",sampleFiles)
treatment <- factor(c("Con","Con","Con","Con","Con","Dox","Dox","Dox","Dox","Con","Con","Con","Con","Con","Dox","Dox","Dox","Dox"))
time <- factor(c("0hour","12hour","24hour","3hour","6hour","12hour","24hour","3hour","6hour","0hour","12hour","24hour","3hour","6hour","12hour","24hour","3hour","6hour"))
group <- factor(c("A","A","A","A","A","A","A","A","A","B","B","B","B","B","B","B","B","B"))

sampleTable <- data.frame(sampleNames=sampleNames,fileName=sampleFiles,treatment=treatment,time=time,group=group)

####DESeqDataSet object and constructors####
dds <- DESeqDataSetFromHTSeqCount(sampleTable,".",design =~treatment+time)

#Pre-filtering the dataset
nrow(dds)
dds <- dds[rowSums(counts(dds)) > 1,]
nrow(dds)

#Accessor functions for the normalization factors in a DESeqDataSet object.
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
#dds <- nbinomWaldTest(dds)

# using rlog transformed data and PCA analysis
rld <- rlog(dds,fitType ="local")
dists <- dist(t(assay(rld)))
plotPCA(rld,intgroup=c("time","treatment","group"))

####differential expressed gene####
dds$condition <- factor(paste0(dds$treatment,dds$time))
design(dds) <- ~condition
dds <- DESeq(dds)
resultsNames(dds)

resDox6_0 <- results(dds,contrast = c("condition","Dox6hour","Con0hour"),independentFiltering = F)
 resDox6_0['ENSMUSG00000024406',]
write.csv(as.data.frame(resDox6_0),file = "t.csv")

#####
results(dds,contrast = c("condition","Dox3hour","Con0hour"))
resDox3_0 <- results(dds,contrast = c("condition","Dox3hour","Con0hour"))
write.csv(as.data.frame(resDox3_0),file = "MA9Dox3_0_wald.csv")

results(dds,contrast = c("condition","Dox6hour","Con0hour"))
resDox6_0 <- results(dds,contrast = c("condition","Dox6hour","Con0hour"))
write.csv(as.data.frame(resDox6_0),file = "MA9Dox6_0_wald.csv")

results(dds,contrast = c("condition","Dox12hour","Con0hour"))
resDox12_0 <- results(dds,contrast = c("condition","Dox12hour","Con0hour"))
write.csv(as.data.frame(resDox12_0),file = "MA9Dox12_0_wald.csv")

results(dds,contrast = c("condition","Dox24hour","Con0hour"))
resDox24_0 <- results(dds,contrast = c("condition","Dox24hour","Con0hour"))
write.csv(as.data.frame(resDox24_0),file = "MA9Dox24_0_wald.csv")

results(dds,contrast = c("condition","Con3hour","Con0hour"))
resCon3_0 <- results(dds,contrast = c("condition","Con3hour","Con0hour"))
write.csv(as.data.frame(resCon3_0),file = "MA9Con3_0_wald.csv")

results(dds,contrast = c("condition","Con6hour","Con0hour"))
resCon6_0 <- results(dds,contrast = c("condition","Con6hour","Con0hour"))
write.csv(as.data.frame(resCon6_0),file = "MA9Con6_0_wald.csv")

results(dds,contrast = c("condition","Con12hour","Con0hour"))
resCon12_0 <- results(dds,contrast = c("condition","Con12hour","Con0hour"))
write.csv(as.data.frame(resCon12_0),file = "MA9Con12_0_wald.csv")

results(dds,contrast = c("condition","Con24hour","Con0hour"))
resCon24_0 <- results(dds,contrast = c("condition","Con24hour","Con0hour"))
write.csv(as.data.frame(resCon24_0),file = "MA9Con24_0_wald.csv")

######
MA9 <- counts(dds, normalized=TRUE) + 0.01
write.csv(MA9,file = 'E:/Lu_ATAC-seq/2. remove_bl/MA9_normalized_number.csv')


MA9 <- as.data.frame(MA9)
MA9$MA9_A_Con3vsCon0 <- MA9$MA9_A_Con3h_count/MA9$MA9_A_Con0h_count
MA9$MA9_A_Con6vsCon0 <- MA9$MA9_A_Con6h_count/MA9$MA9_A_Con0h_count
MA9$MA9_A_Con12vsCon0 <- MA9$MA9_A_Con12h_count/MA9$MA9_A_Con0h_count
MA9$MA9_A_Con24vsCon0 <- MA9$MA9_A_Con24h_count/MA9$MA9_A_Con0h_count

MA9$MA9_A_Dox3vsCon0 <- MA9$MA9_A_Dox3h_count/MA9$MA9_A_Con0h_count
MA9$MA9_A_Dox6vsCon0 <- MA9$MA9_A_Dox6h_count/MA9$MA9_A_Con0h_count
MA9$MA9_A_Dox12vsCon0 <- MA9$MA9_A_Dox12h_count/MA9$MA9_A_Con0h_count
MA9$MA9_A_Dox24vsCon0 <- MA9$MA9_A_Dox24h_count/MA9$MA9_A_Con0h_count

MA9$MA9_B_Con3vsCon0 <- MA9$MA9_B_Con3h_count/MA9$MA9_B_Con0h_count
MA9$MA9_B_Con6vsCon0 <- MA9$MA9_B_Con6h_count/MA9$MA9_B_Con0h_count
MA9$MA9_B_Con12vsCon0 <- MA9$MA9_B_Con12h_count/MA9$MA9_B_Con0h_count
MA9$MA9_B_Con24vsCon0 <- MA9$MA9_B_Con24h_count/MA9$MA9_B_Con0h_count

MA9$MA9_B_Dox3vsCon0 <- MA9$MA9_B_Dox3h_count/MA9$MA9_B_Con0h_count
MA9$MA9_B_Dox6vsCon0 <- MA9$MA9_B_Dox6h_count/MA9$MA9_B_Con0h_count
MA9$MA9_B_Dox12vsCon0 <- MA9$MA9_B_Dox12h_count/MA9$MA9_B_Con0h_count
MA9$MA9_B_Dox24vsCon0 <- MA9$MA9_B_Dox24h_count/MA9$MA9_B_Con0h_count

write.csv(MA9[19:ncol(MA9)],file = 'E:/Lu_ATAC-seq/2. remove_bl/MA9_normalized_fc.csv')









WT <- counts(dds, normalized=TRUE) + 0.01
write.csv(MA9,file = 'E:/Lu_ATAC-seq/2. remove_bl//WT_normalized_number.csv')

WT <- as.data.frame(WT)
WT$WT_A_Con3vsCon0 <- WT$WT_A_Con3h_count/WT$WT_A_Con0h_count
WT$WT_A_Con6vsCon0 <- WT$WT_A_Con6h_count/WT$WT_A_Con0h_count
WT$WT_A_Con12vsCon0 <- WT$WT_A_Con12h_count/WT$WT_A_Con0h_count
WT$WT_A_Con24vsCon0 <- WT$WT_A_Con24h_count/WT$WT_A_Con0h_count

WT$WT_A_Dox3vsCon0 <- WT$WT_A_Dox3h_count/WT$WT_A_Con0h_count
WT$WT_A_Dox6vsCon0 <- WT$WT_A_Dox6h_count/WT$WT_A_Con0h_count
WT$WT_A_Dox12vsCon0 <- WT$WT_A_Dox12h_count/WT$WT_A_Con0h_count
WT$WT_A_Dox24vsCon0 <- WT$WT_A_Dox24h_count/WT$WT_A_Con0h_count

WT$WT_B_Con3vsCon0 <- WT$WT_B_Con3h_count/WT$WT_B_Con0h_count
WT$WT_B_Con6vsCon0 <- WT$WT_B_Con6h_count/WT$WT_B_Con0h_count
WT$WT_B_Con12vsCon0 <- WT$WT_B_Con12h_count/WT$WT_B_Con0h_count
WT$WT_B_Con24vsCon0 <- WT$WT_B_Con24h_count/WT$WT_B_Con0h_count

WT$WT_B_Dox3vsCon0 <- WT$WT_B_Dox3h_count/WT$WT_B_Con0h_count
WT$WT_B_Dox6vsCon0 <- WT$WT_B_Dox6h_count/WT$WT_B_Con0h_count
WT$WT_B_Dox12vsCon0 <- WT$WT_B_Dox12h_count/WT$WT_B_Con0h_count
WT$WT_B_Dox24vsCon0 <- WT$WT_B_Dox24h_count/WT$WT_B_Con0h_count

write.csv(WT[19:ncol(WT)],file = 'E:/Lu_ATAC-seq/2. remove_bl//WT_normalized_fc.csv')






