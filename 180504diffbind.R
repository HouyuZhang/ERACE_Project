####set working directory####
library(DiffBind)
library(rtracklayer)
library(ChIPseeker)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(ComplexHeatmap)
library(circlize)
setwd("E:/Dropbox/Luting/Lab/data analysis/ATAC-seq/MA9/macs2")

####improt data####
MA9 <- dba(sampleSheet = "MA9_new.csv")

####check data is proproly loaded and correalation between dataset(using occupancy,peak caller score data)####
MA9
plot(MA9)

####export raw count matrix form diffbind if necessary####
raw_count <- dba.peakset(MA9,bRetrieve = TRUE,DataType=DBA_DATA_FRAME,writeFile = "MA9_raw_count.txt")

####overlap rates####
olap.rate <- dba.overlap(MA9,mode = DBA_OLAP_RATE)
olap.rate
plot(olap.rate,type='b',ylab='#peaks',xlab='Overlap at least this many peaksets')

####deriving consensus peakset####
dba.overlap(MA9,MA9$masks$Con)


####counting reads ####
MA9_count <- dba.count(MA9,
                       summits = 250,
                       minOverlap = 2)

####export count matrix form diffbind if necessary####
counts <- dba.peakset(MA9_count,bRetrieve = TRUE,DataType=DBA_DATA_FRAME,writeFile = "MA9_count_summit250.txt")

####correalation between dataset(using affinity scores,read count)####

MA9_count
plot(WT_count$peaks)
dba.plotPCA(MA9_count,label = DBA_TREATMENT)




MA9_count <- dba.load("MA9")
####establishing a contrast####
MA9_conrast <- dba.contrast(MA9_count,categories = DBA_CONDITION,minMembers = 2)
MA9_conrast <- dba.analyze(MA9_conrast)
MA9_conrast

MA9_Con0_Dox3.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=5,th=1)
write.csv(MA9_Con0_Dox3.DB,file = "MA9Con0_VS_Dox3_500.csv")
MA9_Dox24Anno <- annotatePeak(MA9_Con0_Dox3.DB, tssRegion=c(-500, 500),TxDb=txdb, annoDb="org.Mm.eg.db")
write.csv(MA9_Dox24Anno,file = "MA9Con0_VS_Dox3_500_anno.csv")

MA9_Con0_Dox6.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=6,th=1)
write.csv(MA9_Con0_Dox6.DB,file = "MA9Con0_VS_Dox6_500.csv")

MA9_Con0_Dox12.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=7,th=1)
write.csv(MA9_Con0_Dox12.DB,file = "MA9Con0_VS_Dox12_500.csv")


MA9_Con0_Dox24.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=8,th=1)
write.csv(MA9_Con0_Dox24.DB,file = "MA9Con0_VS_Dox24_500.csv")


MA9_Con0_Con3.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=1,th=1)
write.csv(MA9_Con0_Con3.DB,file = "MA9Con0_VS_Con3_500.csv")


MA9_Con0_Con6.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=2,th=1)
write.csv(MA9_Con0_Con6.DB,file = "MA9Con0_VS_Con6_500.csv")


MA9_Con0_Con12.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=3,th=1)
write.csv(MA9_Con0_Con12.DB,file = "MA9Con0_VS_Con12_500.csv")

MA9_Con0_Con24.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=4,th=1)
write.csv(MA9_Con0_Con24.DB,file = "MA9Con0_VS_Con24_500.csv")






WT_count <- dba.load("WT")
WT_conrast <- dba.contrast(WT_count,categories = DBA_CONDITION,minMembers = 2)
WT_conrast <- dba.analyze(WT_conrast)
WT_conrast

WT_Con0_Dox3.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=5,th=1)
write.csv(WT_Con0_Dox3.DB,file = "WTCon0_VS_Dox3_500.csv")
MA9_Dox24Anno <- annotatePeak(WT_Con0_Dox3.DB, tssRegion=c(-500, 500),TxDb=txdb, annoDb="org.Mm.eg.db")
write.csv(MA9_Dox24Anno,file = "WTCon0_VS_Dox3_500_anno.csv")

WT_Con0_Dox6.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=6,th=1)
write.csv(WT_Con0_Dox6.DB,file = "WTCon0_VS_Dox6_500.csv")


WT_Con0_Dox12.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=7,th=1)
write.csv(WT_Con0_Dox12.DB,file = "WTCon0_VS_Dox12_500.csv")


WT_Con0_Dox24.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=8,th=1)
write.csv(WT_Con0_Dox24.DB,file = "WTCon0_VS_Dox24_500.csv")


WT_Con0_Con3.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=1,th=1)
write.csv(WT_Con0_Con3.DB,file = "WTCon0_VS_Con3_500.csv")


WT_Con0_Con6.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=2,th=1)
write.csv(WT_Con0_Con6.DB,file = "WTCon0_VS_Con6_500.csv")


WT_Con0_Con12.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=3,th=1)
write.csv(WT_Con0_Con12.DB,file = "WTCon0_VS_Con12_500.csv")


WT_Con0_Con24.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=4,th=1)
write.csv(WT_Con0_Con24.DB,file = "WTCon0_VS_Con24_500.csv")


plotDistToTSS(MA9_Dox3Anno_10000)

MA9_Con0_Con3.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=1,th=0.05)
MA9_Con3Anno <- annotatePeak(MA9_Con0_Con3.DB, tssRegion=c(-3000, 3000),
                                 TxDb=txdb, annoDb="org.Mm.eg.db")
write.csv(MA9_Con3Anno_500,file = 'E:/Lu_ATAC-seq/diffbind motif calculated with occupancy/anno_gene/500/MA9_Con3Anno_500.csv')

MA9_Con0_Con6.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=2,th=0.05)
MA9_Con6Anno <- annotatePeak(MA9_Con0_Con6.DB, tssRegion=c(-3000, 3000),
                                 TxDb=txdb, annoDb="org.Mm.eg.db")
write.csv(MA9_Con6Anno_500,file = 'E:/Lu_ATAC-seq/diffbind motif calculated with occupancy/anno_gene/500/MA9_Con6Anno_500.csv')

MA9_Con0_Con12.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=3,th=0.05)
MA9_Con12Anno <- annotatePeak(MA9_Con0_Con12.DB, tssRegion=c(-3000, 3000),
                                  TxDb=txdb, annoDb="org.Mm.eg.db")
write.csv(MA9_Con12Anno_500,file = 'E:/Lu_ATAC-seq/diffbind motif calculated with occupancy/anno_gene/500/MA9_Con12Anno_500.csv')

MA9_Con0_Con24.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=4,th=0.05)
MA9_Con24Anno <- annotatePeak(MA9_Con0_Con24.DB,tssRegion=c(-3000, 3000),
                                  TxDb=txdb, annoDb="org.Mm.eg.db")
write.csv(MA9_Con24Anno_500,file = 'E:/Lu_ATAC-seq/diffbind motif calculated with occupancy/anno_gene/500/MA9_Con24Anno_500.csv')

MA9_Con0_Dox3.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=5,th=0.05)
MA9_Dox3Anno <- annotatePeak(MA9_Con0_Dox3.DB,tssRegion=c(-3000, 3000),
                                 TxDb=txdb, annoDb="org.Mm.eg.db")
write.csv(MA9_Dox3Anno_500,file = 'E:/Lu_ATAC-seq/diffbind motif calculated with occupancy/anno_gene/500/MA9_Dox3Anno_500.csv')

MA9_Con0_Dox6.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=6,th=0.05)
MA9_Dox6Anno <- annotatePeak(MA9_Con0_Dox6.DB,tssRegion=c(-3000, 3000),
                                 TxDb=txdb, annoDb="org.Mm.eg.db")
write.csv(MA9_Dox6Anno_500,file = 'E:/Lu_ATAC-seq/diffbind motif calculated with occupancy/anno_gene/500/MA9_Dox6Anno_500.csv')

MA9_Con0_Dox12.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=7,th=0.05)
MA9_Dox12Anno <- annotatePeak(MA9_Con0_Dox12.DB, tssRegion=c(-3000, 3000),
                                  TxDb=txdb, annoDb="org.Mm.eg.db")
write.csv(MA9_Dox12Anno_500,file = 'E:/Lu_ATAC-seq/diffbind motif calculated with occupancy/anno_gene/500/MA9_Dox12Anno_500.csv')

MA9_Con0_Dox24.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=8,th=0.05)
MA9_Dox24Anno <- annotatePeak(MA9_Con0_Dox24.DB,tssRegion=c(-3000, 3000),
                                  TxDb=txdb, annoDb="org.Mm.eg.db")
write.csv(MA9_Dox24Anno_500,file = 'E:/Lu_ATAC-seq/diffbind motif calculated with occupancy/anno_gene/500/MA9_Dox24Anno_500.csv')


WT_Con0_Con3.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=1,th=0.05)
WT_Con3Anno <- annotatePeak(WT_Con0_Con3.DB,tssRegion=c(-3000, 3000),
                                TxDb=txdb, annoDb="org.Mm.eg.db")
write.csv(WT_Con3Anno_500,file = 'E:/Lu_ATAC-seq/diffbind motif calculated with occupancy/anno_gene/500/WT_Con3Anno_500.csv')

WT_Con0_Con6.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=2,th=0.05)
WT_Con6Anno <- annotatePeak(WT_Con0_Con6.DB,tssRegion=c(-3000, 3000),
                                TxDb=txdb, annoDb="org.Mm.eg.db")
write.csv(WT_Con6Anno_500,file = 'E:/Lu_ATAC-seq/diffbind motif calculated with occupancy/anno_gene/500/WT_Con6Anno_500.csv')

WT_Con0_Con12.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=3,th=0.05)
WT_Con12Anno <- annotatePeak(WT_Con0_Con12.DB,tssRegion=c(-3000, 3000),
                                 TxDb=txdb, annoDb="org.Mm.eg.db")
write.csv(WT_Con12Anno_500,file = 'E:/Lu_ATAC-seq/diffbind motif calculated with occupancy/anno_gene/500/WT_Con12Anno_500.csv')

WT_Con0_Con24.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=4,th=0.05)
WT_Con24Anno <- annotatePeak(WT_Con0_Con24.DB,tssRegion=c(-3000, 3000),
                                 TxDb=txdb, annoDb="org.Mm.eg.db")
write.csv(WT_Con24Anno_500,file = 'E:/Lu_ATAC-seq/diffbind motif calculated with occupancy/anno_gene/500/WT_Con24Anno_500.csv')

WT_Con0_Dox3.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=5,th=0.05)
WT_Dox3Anno <- annotatePeak(WT_Con0_Dox3.DB,tssRegion=c(-3000, 3000),
                                TxDb=txdb, annoDb="org.Mm.eg.db")
write.csv(WT_Dox3Anno_500,file = 'E:/Lu_ATAC-seq/diffbind motif calculated with occupancy/anno_gene/500/WT_Dox3Anno_500.csv')

WT_Con0_Dox6.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=6,th=0.05)
WT_Dox6Anno <- annotatePeak(WT_Con0_Dox6.DB,tssRegion=c(-3000, 3000),
                                TxDb=txdb, annoDb="org.Mm.eg.db")
write.csv(WT_Dox6Anno_500,file = 'E:/Lu_ATAC-seq/diffbind motif calculated with occupancy/anno_gene/500/WT_Dox6Anno_500.csv')

WT_Con0_Dox12.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=7,th=0.05)
WT_Dox12Anno <- annotatePeak(WT_Con0_Dox12.DB,tssRegion=c(-3000, 3000),
                                 TxDb=txdb, annoDb="org.Mm.eg.db")
write.csv(WT_Dox12Anno_500,file = 'E:/Lu_ATAC-seq/diffbind motif calculated with occupancy/anno_gene/500/WT_Dox12Anno_500.csv')

WT_Con0_Dox24.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=8,th=0.05)
WT_Dox24Anno <- annotatePeak(WT_Con0_Dox24.DB,tssRegion=c(-3000, 3000),
                                 TxDb=txdb, annoDb="org.Mm.eg.db")
write.csv(WT_Dox24Anno_500,file = 'E:/Lu_ATAC-seq/diffbind motif calculated with occupancy/anno_gene/500/WT_Dox24Anno_500.csv')


plotAnnoPie(MA9_Con3Anno)
plotAnnoPie(MA9_Con6Anno)
plotAnnoPie(MA9_Con12Anno)
plotAnnoPie(MA9_Con24Anno)
plotAnnoPie(MA9_Dox3Anno)
plotAnnoPie(MA9_Dox6Anno)
plotAnnoPie(MA9_Dox12Anno)
plotAnnoPie(MA9_Dox24Anno)


####annotation percentage####



####annotation pathway####
#MA9_0_3_list <- lapply(MA9_Con0_Dox3.DB, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000),verbose=FALSE)
#plotAnnoBar(peakAnnoList)
#plotDistToTSS(peakAnnoList)

pathway0_3_1 <- enrichPathway(as.data.frame(Dox3Anno)$geneId,organism = "mouse")
head(pathway0_3,2)
#gene <- seq2gene(MA9_Con0_Dox3.DB,tssRegion = c(-1000,1000),flankDistance = 3000,TxDb = txdb)
#pathway0_3 <- enrichPathway(gene)
#head(pathway0_3,2)
dotplot(pathway0_3_1)

pathway3_6_1 <- enrichPathway(as.data.frame(Dox3_6Anno)$geneId,organism = "mouse")
head(pathway3_6_1,2)
dotplot(pathway3_6_1)


####seq2gene pathway entichment####
MA9_Dox0_3_gene <- seq2gene(MA9_Con0_Dox3.DB,tssRegion = c(-1000,1000),flankDistance = 5000,TxDb = txdb)
pathway_Con0_Dox3 <- enrichPathway(MA9_Dox0_3_gene,organism = "mouse")
dotplot(pathway_Con0_Dox3,showCategory=15)

MA9_Dox3_6_gene <- seq2gene(MA9_Dox3_Dox6.DB,tssRegion = c(-1000,1000),flankDistance = 5000,TxDb = txdb)
pathway_Dox3_Dox6 <- enrichPathway(MA9_Dox3_6_gene,organism = "mouse")
dotplot(pathway_Dox3_Dox6,showCategory=15)

MA9_Dox6_12_gene <- seq2gene(MA9_Dox6_Dox12.DB,tssRegion = c(-1000,1000),flankDistance = 5000,TxDb = txdb)
pathway_Dox6_Dox12 <- enrichPathway(MA9_Dox6_12_gene,organism = "mouse")
dotplot(pathway_Dox6_Dox12,showCategory=15)

MA9_Dox12_24_gene <- seq2gene(MA9_Dox12_Dox24.DB,tssRegion = c(-1000,1000),flankDistance = 5000,TxDb = txdb)
pathway_Dox12_Dox24 <- enrichPathway(MA9_Dox12_24_gene,organism = "mouse")
dotplot(pathway_Dox12_Dox24,showCategory=15)

####Go entichmen####
MA9_Dox0_3 <- enrichGO(gene   = as.data.frame(MA9_Dox12Anno_500)$geneId,
                       OrgDb         = org.Mm.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = T)
dotplot(MA9_Dox0_3,showCategory=30)

MA9_Dox3_6 <- enrichGO(gene   = MA9_Dox3_6_gene,
                       OrgDb         = org.Mm.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = T)
dotplot(MA9_Dox3_6,showCategory=10)

MA9_Dox6_12 <- enrichGO(gene   = MA9_Dox6_12_gene,
                        OrgDb         = org.Mm.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = T)
dotplot(MA9_Dox6_12,showCategory=10)

MA9_Dox12_24 <- enrichGO(gene   = MA9_Dox12_24_gene,
                         OrgDb         = org.Mm.eg.db,
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable      = T)
dotplot(MA9_Dox12_24,showCategory=10)

########
MA9_Dox0_3 <- enrichGO(gene   = MA9_Dox0_3_gene,
                       OrgDb         = org.Mm.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = T)
MA9_Dox0_3 <- dropGO(MA9_Dox0_3,level = 15)
dotplot(MA9_Dox0_3,showCategory=15)

mydf <- read.csv("diffpeak_gene.csv",header = T)
mydf <- data.frame(Entrez=mydf$gene,group=mydf$group)
formula_res <- compareCluster(Entrez~group,
                              data = mydf,
                              fun  ='enrichGO', 
                              ont = "BP",
                              OrgDb    ='org.Mm.eg.db',
                              readable = T)
head(as.data.frame(formula_res))
dotplot(formula_res,showCategory=2)

