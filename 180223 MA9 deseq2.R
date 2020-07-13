library(DESeq2)
directory <- "F:/Dropbox/1_Yenlab/1_1_Documents/Lu, Zhang share/2018 Apoptosis project/ATAC-seq/190411 re-annotate peak TSS+-1kb/genecount/"
setwd("F:/Dropbox/1_Yenlab/1_1_Documents/Lu, Zhang share/2018 Apoptosis project/ATAC-seq/190411 re-annotate peak TSS+-1kb/genecount/")

#### can merge individual sample files (i.e. ctrl1.counts, ctrl2.counts, etc.)####
sampleFiles <- grep("MA9.*count",list.files(directory),value=T)

treatment <- factor(c("Con","Con","Con","Con","Con","Dox","Dox","Dox","Dox","Con","Con","Con","Con","Con","Dox","Dox","Dox","Dox"))
time <- factor(c("0hour","12hour","24hour","3hour","6hour","12hour","24hour","3hour","6hour","0hour","12hour","24hour","3hour","6hour","12hour","24hour","3hour","6hour"))
group <- factor(c("A","A","A","A","A","A","A","A","A","B","B","B","B","B","B","B","B","B"))

sampleTable <- data.frame(sampleNames=sampleFiles,fileName=sampleFiles,treatment=treatment,time=time,group=group)

####DESeqDataSet object and constructors####
dds <- DESeqDataSetFromHTSeqCount(sampleTable,directory,design =~treatment+time)

####differential expressed gene####
dds$condition <- factor(paste0(dds$treatment,dds$time))
design(dds) <- ~condition
dds <- DESeq(dds)
resultsNames(dds)

write.csv(counts(dds, normalized=TRUE),"MA9_ATAC_genescounts_DESeq_normalized.csv")


















results(dds,contrast = c("condition","Dox3hour","Con0hour"))
resDox3_0 <- results(dds,contrast = c("condition","Dox3hour","Con0hour"))
write.csv(as.data.frame(resDox3_0),file = "MA9Dox3_0_wald.csv")

results(dds,contrast = c("condition","Dox6hour","Dox3hour"))
resDox6_3 <- results(dds,contrast = c("condition","Dox6hour","Dox3hour"))
write.csv(as.data.frame(resDox6_3),file = "MA9Dox6_3_wald.csv")

results(dds,contrast = c("condition","Dox12hour","Dox6hour"))
resDox12_6 <- results(dds,contrast = c("condition","Dox12hour","Dox6hour"))
write.csv(as.data.frame(resDox12_6),file = "MA9Dox12_6_wald.csv")

results(dds,contrast = c("condition","Dox24hour","Dox12hour"))
resDox24_12 <- results(dds,contrast = c("condition","Dox24hour","Dox12hour"))
write.csv(as.data.frame(resDox24_12),file = "MA9Dox24_12_wald.csv")

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

results(dds,contrast = c("condition","Dox6hour","Con0hour"))
resDox6_0 <- results(dds,contrast = c("condition","Dox6hour","Con0hour"))
write.csv(as.data.frame(resDox6_0),file = "MA9Dox6_0_wald.csv")

results(dds,contrast = c("condition","Dox24hour","Con0hour"))
resDox24_0 <- results(dds,contrast = c("condition","Dox24hour","Con0hour"))
write.csv(as.data.frame(resDox24_0),file = "MA9Dox24_0_wald.csv")



