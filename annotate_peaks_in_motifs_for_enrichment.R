library(ChIPseeker)
library(GenomicRanges)
library(clusterProfiler)
require(TxDb.Mmusculus.UCSC.mm10.ensGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene

files <- list.files(path = ".", pattern = "*VS_Dox*")

for (i in files){
base_name <- tools::file_path_sans_ext(i)
a <- read.table(i,sep = "\t")
my.granges <- GRanges(seqnames = a$V1,
                      ranges = IRanges(start = a$V2,
                                       end = a$V3,
                                       names = seq(1,nrow(a))),
                      PeakID = a$V4)

b <- annotatePeak(my.granges, tssRegion=c(-1000, 1000),TxDb=txdb, annoDb="org.Mm.eg.db")
i <- as.data.frame(b)
name <- i[i$distanceToTSS<=1000 & i$distanceToTSS >=-1000,]$geneId
write.csv(name,paste(base_name,".csv",sep=""))
}

