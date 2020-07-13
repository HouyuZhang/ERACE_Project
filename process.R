MA9_Con0_Con3.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=1,th=1)
MA9_Con3Anno_3000 <- annotatePeak(MA9_Con0_Con3.DB, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")

MA9_Con0_Con6.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=2,th=1)
MA9_Con6Anno_3000 <- annotatePeak(MA9_Con0_Con6.DB, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")


MA9_Con0_Con12.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=3,th=1)
MA9_Con12Anno_3000 <- annotatePeak(MA9_Con0_Con12.DB, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")


MA9_Con0_Con24.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=4,th=1)
MA9_Con24Anno_3000 <- annotatePeak(MA9_Con0_Con24.DB, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")


MA9_Con0_Dox3.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=5,th=1)
MA9_Dox3Anno_3000 <- annotatePeak(MA9_Con0_Dox3.DB, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")


MA9_Con0_Dox6.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=6,th=1)
MA9_Dox6Anno_3000 <- annotatePeak(MA9_Con0_Dox6.DB, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")


MA9_Con0_Dox12.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=7,th=1)
MA9_Dox12Anno_3000 <- annotatePeak(MA9_Con0_Dox12.DB, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")


MA9_Con0_Dox24.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=8,th=1)
MA9_Dox24Anno_3000 <- annotatePeak(MA9_Con0_Dox24.DB, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")



WT_Con0_Con3.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=1,th=1)
WT_Con3Anno_3000 <- annotatePeak(WT_Con0_Con3.DB, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")


WT_Con0_Con6.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=2,th=1)
WT_Con6Anno_3000 <- annotatePeak(WT_Con0_Con6.DB, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")


WT_Con0_Con12.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=3,th=1)
WT_Con12Anno_3000 <- annotatePeak(WT_Con0_Con12.DB, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")


WT_Con0_Con24.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=4,th=1)
WT_Con24Anno_3000 <- annotatePeak(WT_Con0_Con24.DB, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")


WT_Con0_Dox3.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=5,th=1)
WT_Dox3Anno_3000 <- annotatePeak(WT_Con0_Dox3.DB, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")


WT_Con0_Dox6.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=6,th=1)
WT_Dox6Anno_3000 <- annotatePeak(WT_Con0_Dox6.DB, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")


WT_Con0_Dox12.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=7,th=1)
WT_Dox12Anno_3000 <- annotatePeak(WT_Con0_Dox12.DB, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")


WT_Con0_Dox24.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=8,th=1)
WT_Dox24Anno_3000 <- annotatePeak(WT_Con0_Dox24.DB, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")


peakAnnoList <- list(MA9_Dox3 = as.data.frame(MA9_Dox3Anno_3000),
                     MA9_Dox6 = as.data.frame(MA9_Dox6Anno_3000),
                     MA9_Dox12 = as.data.frame(MA9_Dox12Anno_3000),
                     MA9_Dox24 = as.data.frame(MA9_Dox24Anno_3000),
                     #WT_Dox3 = as.data.frame(WT_Dox3Anno_3000),
                     WT_Dox6 = as.data.frame(WT_Dox6Anno_3000),
                     WT_Dox12 = as.data.frame(WT_Dox12Anno_3000),
                     WT_Dox24 = as.data.frame(WT_Dox24Anno_3000)                    
                     )
genes <- lapply(peakAnnoList, function(i) i[i$distanceToTSS<=2000 & i$distanceToTSS >=-2000,]$geneId)

formula_res <- compareCluster(geneCluster   = genes,
                              fun      ='enrichGO', 
                              ont = "BP",
                              OrgDb    ='org.Mm.eg.db',
                              readable = T)

dotplot(formula_res,showCategory=10,title='TSS+-3000  2000 filtered')



MA9_anno <- read.csv("MA9Con0_VS_Dox3_500_anno.csv",header = T,row.names = 1)
MA9_Con0_VS_Dox3 <- read.csv("MA9Con0_VS_Dox3_500.csv",header = T,row.names = 1)
MA9_Con0_VS_Dox6 <- read.csv("MA9Con0_VS_Dox6_500.csv",header = T,row.names = 1)
MA9_Con0_VS_Dox12 <- read.csv("MA9Con0_VS_Dox12_500.csv",header = T,row.names = 1)
MA9_Con0_VS_Dox24 <- read.csv("MA9Con0_VS_Dox24_500.csv",header = T,row.names = 1)

MA9_Con0_VS_Con3 <- read.csv("MA9Con0_VS_Con3_500.csv",header = T,row.names = 1)
MA9_Con0_VS_Con6 <- read.csv("MA9Con0_VS_Con6_500.csv",header = T,row.names = 1)
MA9_Con0_VS_Con12 <- read.csv("MA9Con0_VS_Con12_500.csv",header = T,row.names = 1)
MA9_Con0_VS_Con24 <- read.csv("MA9Con0_VS_Con24_500.csv",header = T,row.names = 1)

#####
MA9_Con0_VS_Dox3$distanceToTSS <- MA9_anno$distanceToTSS
MA9_Con0_VS_Dox3$ENSEMBL <- MA9_anno$ENSEMBL

#把所有注视到一个基因上的peak内的count数字加起来
one <- with (MA9_Con0_VS_Dox3 , 
             tapply (MA9_Con0_VS_Dox3[MA9_Con0_VS_Dox3$distanceToTSS<1000 & MA9_Con0_VS_Dox3$distanceToTSS >= -1000,][[12]] , 
                     list (MA9_Con0_VS_Dox3[MA9_Con0_VS_Dox3$distanceToTSS<1000 & MA9_Con0_VS_Dox3$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
two <- with (MA9_Con0_VS_Dox3 , 
             tapply (MA9_Con0_VS_Dox3[MA9_Con0_VS_Dox3$distanceToTSS<1000 & MA9_Con0_VS_Dox3$distanceToTSS >= -1000,][[13]] , 
                     list (MA9_Con0_VS_Dox3[MA9_Con0_VS_Dox3$distanceToTSS<1000 & MA9_Con0_VS_Dox3$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
three <- with (MA9_Con0_VS_Dox3 , 
               tapply (MA9_Con0_VS_Dox3[MA9_Con0_VS_Dox3$distanceToTSS<1000 & MA9_Con0_VS_Dox3$distanceToTSS >= -1000,][[14]] , 
                       list (MA9_Con0_VS_Dox3[MA9_Con0_VS_Dox3$distanceToTSS<1000 & MA9_Con0_VS_Dox3$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
four <- with (MA9_Con0_VS_Dox3 , 
              tapply (MA9_Con0_VS_Dox3[MA9_Con0_VS_Dox3$distanceToTSS<1000 & MA9_Con0_VS_Dox3$distanceToTSS >= -1000,][[15]] , 
                      list (MA9_Con0_VS_Dox3[MA9_Con0_VS_Dox3$distanceToTSS<1000 & MA9_Con0_VS_Dox3$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )

MA9_Con0_VS_Dox3 <- as.data.frame(cbind(one,two,three,four))


MA9_Con0_VS_Dox6$distanceToTSS <- MA9_anno$distanceToTSS
MA9_Con0_VS_Dox6$ENSEMBL <- MA9_anno$ENSEMBL

one <- with (MA9_Con0_VS_Dox6 , 
             tapply (MA9_Con0_VS_Dox6[MA9_Con0_VS_Dox6$distanceToTSS<1000 & MA9_Con0_VS_Dox6$distanceToTSS >= -1000,][[12]] , 
                     list (MA9_Con0_VS_Dox6[MA9_Con0_VS_Dox6$distanceToTSS<1000 & MA9_Con0_VS_Dox6$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
two <- with (MA9_Con0_VS_Dox6 , 
             tapply (MA9_Con0_VS_Dox6[MA9_Con0_VS_Dox6$distanceToTSS<1000 & MA9_Con0_VS_Dox6$distanceToTSS >= -1000,][[13]] , 
                     list (MA9_Con0_VS_Dox6[MA9_Con0_VS_Dox6$distanceToTSS<1000 & MA9_Con0_VS_Dox6$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
three <- with (MA9_Con0_VS_Dox6 , 
               tapply (MA9_Con0_VS_Dox6[MA9_Con0_VS_Dox6$distanceToTSS<1000 & MA9_Con0_VS_Dox6$distanceToTSS >= -1000,][[14]] , 
                       list (MA9_Con0_VS_Dox6[MA9_Con0_VS_Dox6$distanceToTSS<1000 & MA9_Con0_VS_Dox6$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
four <- with (MA9_Con0_VS_Dox6 , 
              tapply (MA9_Con0_VS_Dox6[MA9_Con0_VS_Dox6$distanceToTSS<1000 & MA9_Con0_VS_Dox6$distanceToTSS >= -1000,][[15]] , 
                      list (MA9_Con0_VS_Dox6[MA9_Con0_VS_Dox6$distanceToTSS<1000 & MA9_Con0_VS_Dox6$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )

MA9_Con0_VS_Dox6 <- as.data.frame(cbind(one,two,three,four))

MA9_Con0_VS_Dox12$distanceToTSS <- MA9_anno$distanceToTSS
MA9_Con0_VS_Dox12$ENSEMBL <- MA9_anno$ENSEMBL

one <- with (MA9_Con0_VS_Dox12 , 
             tapply (MA9_Con0_VS_Dox12[MA9_Con0_VS_Dox12$distanceToTSS<1000 & MA9_Con0_VS_Dox12$distanceToTSS >= -1000,][[12]] , 
                     list (MA9_Con0_VS_Dox12[MA9_Con0_VS_Dox12$distanceToTSS<1000 & MA9_Con0_VS_Dox12$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
two <- with (MA9_Con0_VS_Dox12 , 
             tapply (MA9_Con0_VS_Dox12[MA9_Con0_VS_Dox12$distanceToTSS<1000 & MA9_Con0_VS_Dox12$distanceToTSS >= -1000,][[13]] , 
                     list (MA9_Con0_VS_Dox12[MA9_Con0_VS_Dox12$distanceToTSS<1000 & MA9_Con0_VS_Dox12$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
three <- with (MA9_Con0_VS_Dox12 , 
               tapply (MA9_Con0_VS_Dox12[MA9_Con0_VS_Dox12$distanceToTSS<1000 & MA9_Con0_VS_Dox12$distanceToTSS >= -1000,][[14]] , 
                       list (MA9_Con0_VS_Dox12[MA9_Con0_VS_Dox12$distanceToTSS<1000 & MA9_Con0_VS_Dox12$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
four <- with (MA9_Con0_VS_Dox12 , 
              tapply (MA9_Con0_VS_Dox12[MA9_Con0_VS_Dox12$distanceToTSS<1000 & MA9_Con0_VS_Dox12$distanceToTSS >= -1000,][[15]] , 
                      list (MA9_Con0_VS_Dox12[MA9_Con0_VS_Dox12$distanceToTSS<1000 & MA9_Con0_VS_Dox12$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )

MA9_Con0_VS_Dox12 <- as.data.frame(cbind(one,two,three,four))


MA9_Con0_VS_Dox24$distanceToTSS <- MA9_anno$distanceToTSS
MA9_Con0_VS_Dox24$ENSEMBL <- MA9_anno$ENSEMBL

one <- with (MA9_Con0_VS_Dox24 , 
             tapply (MA9_Con0_VS_Dox24[MA9_Con0_VS_Dox24$distanceToTSS<1000 & MA9_Con0_VS_Dox24$distanceToTSS >= -1000,][[12]] , 
                     list (MA9_Con0_VS_Dox24[MA9_Con0_VS_Dox24$distanceToTSS<1000 & MA9_Con0_VS_Dox24$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
two <- with (MA9_Con0_VS_Dox24 , 
             tapply (MA9_Con0_VS_Dox24[MA9_Con0_VS_Dox24$distanceToTSS<1000 & MA9_Con0_VS_Dox24$distanceToTSS >= -1000,][[13]] , 
                     list (MA9_Con0_VS_Dox24[MA9_Con0_VS_Dox24$distanceToTSS<1000 & MA9_Con0_VS_Dox24$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
three <- with (MA9_Con0_VS_Dox24 , 
               tapply (MA9_Con0_VS_Dox24[MA9_Con0_VS_Dox24$distanceToTSS<1000 & MA9_Con0_VS_Dox24$distanceToTSS >= -1000,][[14]] , 
                       list (MA9_Con0_VS_Dox24[MA9_Con0_VS_Dox24$distanceToTSS<1000 & MA9_Con0_VS_Dox24$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
four <- with (MA9_Con0_VS_Dox24 , 
              tapply (MA9_Con0_VS_Dox24[MA9_Con0_VS_Dox24$distanceToTSS<1000 & MA9_Con0_VS_Dox24$distanceToTSS >= -1000,][[15]] , 
                      list (MA9_Con0_VS_Dox24[MA9_Con0_VS_Dox24$distanceToTSS<1000 & MA9_Con0_VS_Dox24$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )

MA9_Con0_VS_Dox24 <- as.data.frame(cbind(one,two,three,four))

MA9_Con0_VS_Con3$distanceToTSS <- MA9_anno$distanceToTSS
MA9_Con0_VS_Con3$ENSEMBL <- MA9_anno$ENSEMBL

one <- with (MA9_Con0_VS_Con3 , 
             tapply (MA9_Con0_VS_Con3[MA9_Con0_VS_Con3$distanceToTSS<1000 & MA9_Con0_VS_Con3$distanceToTSS >= -1000,][[12]] , 
                     list (MA9_Con0_VS_Con3[MA9_Con0_VS_Con3$distanceToTSS<1000 & MA9_Con0_VS_Con3$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
two <- with (MA9_Con0_VS_Con3 , 
             tapply (MA9_Con0_VS_Con3[MA9_Con0_VS_Con3$distanceToTSS<1000 & MA9_Con0_VS_Con3$distanceToTSS >= -1000,][[13]] , 
                     list (MA9_Con0_VS_Con3[MA9_Con0_VS_Con3$distanceToTSS<1000 & MA9_Con0_VS_Con3$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
three <- with (MA9_Con0_VS_Con3 , 
               tapply (MA9_Con0_VS_Con3[MA9_Con0_VS_Con3$distanceToTSS<1000 & MA9_Con0_VS_Con3$distanceToTSS >= -1000,][[14]] , 
                       list (MA9_Con0_VS_Con3[MA9_Con0_VS_Con3$distanceToTSS<1000 & MA9_Con0_VS_Con3$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
four <- with (MA9_Con0_VS_Con3 , 
              tapply (MA9_Con0_VS_Con3[MA9_Con0_VS_Con3$distanceToTSS<1000 & MA9_Con0_VS_Con3$distanceToTSS >= -1000,][[15]] , 
                      list (MA9_Con0_VS_Con3[MA9_Con0_VS_Con3$distanceToTSS<1000 & MA9_Con0_VS_Con3$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )

MA9_Con0_VS_Con3 <- as.data.frame(cbind(one,two,three,four))

MA9_Con0_VS_Con6$distanceToTSS <- MA9_anno$distanceToTSS
MA9_Con0_VS_Con6$ENSEMBL <- MA9_anno$ENSEMBL

one <- with (MA9_Con0_VS_Con6 , 
             tapply (MA9_Con0_VS_Con6[MA9_Con0_VS_Con6$distanceToTSS<1000 & MA9_Con0_VS_Con6$distanceToTSS >= -1000,][[12]] , 
                     list (MA9_Con0_VS_Con6[MA9_Con0_VS_Con6$distanceToTSS<1000 & MA9_Con0_VS_Con6$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
two <- with (MA9_Con0_VS_Con6 , 
             tapply (MA9_Con0_VS_Con6[MA9_Con0_VS_Con6$distanceToTSS<1000 & MA9_Con0_VS_Con6$distanceToTSS >= -1000,][[13]] , 
                     list (MA9_Con0_VS_Con6[MA9_Con0_VS_Con6$distanceToTSS<1000 & MA9_Con0_VS_Con6$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
three <- with (MA9_Con0_VS_Con6 , 
               tapply (MA9_Con0_VS_Con6[MA9_Con0_VS_Con6$distanceToTSS<1000 & MA9_Con0_VS_Con6$distanceToTSS >= -1000,][[14]] , 
                       list (MA9_Con0_VS_Con6[MA9_Con0_VS_Con6$distanceToTSS<1000 & MA9_Con0_VS_Con6$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
four <- with (MA9_Con0_VS_Con6 , 
              tapply (MA9_Con0_VS_Con6[MA9_Con0_VS_Con6$distanceToTSS<1000 & MA9_Con0_VS_Con6$distanceToTSS >= -1000,][[15]] , 
                      list (MA9_Con0_VS_Con6[MA9_Con0_VS_Con6$distanceToTSS<1000 & MA9_Con0_VS_Con6$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )

MA9_Con0_VS_Con6 <- as.data.frame(cbind(one,two,three,four))


MA9_Con0_VS_Con12$distanceToTSS <- MA9_anno$distanceToTSS
MA9_Con0_VS_Con12$ENSEMBL <- MA9_anno$ENSEMBL

one <- with (MA9_Con0_VS_Con12 , 
             tapply (MA9_Con0_VS_Con12[MA9_Con0_VS_Con12$distanceToTSS<1000 & MA9_Con0_VS_Con12$distanceToTSS >= -1000,][[12]] , 
                     list (MA9_Con0_VS_Con12[MA9_Con0_VS_Con12$distanceToTSS<1000 & MA9_Con0_VS_Con12$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
two <- with (MA9_Con0_VS_Con12 , 
             tapply (MA9_Con0_VS_Con12[MA9_Con0_VS_Con12$distanceToTSS<1000 & MA9_Con0_VS_Con12$distanceToTSS >= -1000,][[13]] , 
                     list (MA9_Con0_VS_Con12[MA9_Con0_VS_Con12$distanceToTSS<1000 & MA9_Con0_VS_Con12$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
three <- with (MA9_Con0_VS_Con12 , 
               tapply (MA9_Con0_VS_Con12[MA9_Con0_VS_Con12$distanceToTSS<1000 & MA9_Con0_VS_Con12$distanceToTSS >= -1000,][[14]] , 
                       list (MA9_Con0_VS_Con12[MA9_Con0_VS_Con12$distanceToTSS<1000 & MA9_Con0_VS_Con12$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
four <- with (MA9_Con0_VS_Con12 , 
              tapply (MA9_Con0_VS_Con12[MA9_Con0_VS_Con12$distanceToTSS<1000 & MA9_Con0_VS_Con12$distanceToTSS >= -1000,][[15]] , 
                      list (MA9_Con0_VS_Con12[MA9_Con0_VS_Con12$distanceToTSS<1000 & MA9_Con0_VS_Con12$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )

MA9_Con0_VS_Con12 <- as.data.frame(cbind(one,two,three,four))

MA9_Con0_VS_Con24$distanceToTSS <- MA9_anno$distanceToTSS
MA9_Con0_VS_Con24$ENSEMBL <- MA9_anno$ENSEMBL

one <- with (MA9_Con0_VS_Con24 , 
             tapply (MA9_Con0_VS_Con24[MA9_Con0_VS_Con24$distanceToTSS<1000 & MA9_Con0_VS_Con24$distanceToTSS >= -1000,][[12]] , 
                     list (MA9_Con0_VS_Con24[MA9_Con0_VS_Con24$distanceToTSS<1000 & MA9_Con0_VS_Con24$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
two <- with (MA9_Con0_VS_Con24 , 
             tapply (MA9_Con0_VS_Con24[MA9_Con0_VS_Con24$distanceToTSS<1000 & MA9_Con0_VS_Con24$distanceToTSS >= -1000,][[13]] , 
                     list (MA9_Con0_VS_Con24[MA9_Con0_VS_Con24$distanceToTSS<1000 & MA9_Con0_VS_Con24$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
three <- with (MA9_Con0_VS_Con24 , 
               tapply (MA9_Con0_VS_Con24[MA9_Con0_VS_Con24$distanceToTSS<1000 & MA9_Con0_VS_Con24$distanceToTSS >= -1000,][[14]] , 
                       list (MA9_Con0_VS_Con24[MA9_Con0_VS_Con24$distanceToTSS<1000 & MA9_Con0_VS_Con24$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )
four <- with (MA9_Con0_VS_Con24 , 
              tapply (MA9_Con0_VS_Con24[MA9_Con0_VS_Con24$distanceToTSS<1000 & MA9_Con0_VS_Con24$distanceToTSS >= -1000,][[15]] , 
                      list (MA9_Con0_VS_Con24[MA9_Con0_VS_Con24$distanceToTSS<1000 & MA9_Con0_VS_Con24$distanceToTSS >= -1000,]$ENSEMBL) , sum ) )

MA9_Con0_VS_Con24 <- as.data.frame(cbind(one,two,three,four))
######
f_MA9_fc <- data.frame()
f_MA9_fc <- row.names(MA9_Con0_VS_Dox3)
f_MA9_fc <- as.data.frame(f_MA9_fc)

f_MA9_fc$MA9_A_Con3vsCon0 <-MA9_Con0_VS_Con3$three + 0.01 / MA9_Con0_VS_Con3$one + 0.01
f_MA9_fc$MA9_A_Con6vsCon0 <-MA9_Con0_VS_Con6$three + 0.01 / MA9_Con0_VS_Con6$one + 0.01
f_MA9_fc$MA9_A_Con12vsCon0 <-MA9_Con0_VS_Con12$three + 0.01 / MA9_Con0_VS_Con12$one + 0.01
f_MA9_fc$MA9_A_Con24vsCon0 <-MA9_Con0_VS_Con24$three + 0.01 / MA9_Con0_VS_Con24$one + 0.01
f_MA9_fc$MA9_A_Dox3vsCon0 <-MA9_Con0_VS_Dox3$three + 0.01 / MA9_Con0_VS_Dox3$one + 0.01
f_MA9_fc$MA9_A_Dox6vsCon0 <-MA9_Con0_VS_Dox6$three + 0.01 / MA9_Con0_VS_Dox6$one + 0.01
f_MA9_fc$MA9_A_Dox12vsCon0 <-MA9_Con0_VS_Dox12$three + 0.01 / MA9_Con0_VS_Dox12$one + 0.01
f_MA9_fc$MA9_A_Dox24vsCon0 <-MA9_Con0_VS_Dox24$three + 0.01 / MA9_Con0_VS_Dox24$one + 0.01

f_MA9_fc$MA9_B_Con3vsCon0 <-MA9_Con0_VS_Con3$three + 0.01 / MA9_Con0_VS_Con3$two + 0.01
f_MA9_fc$MA9_B_Con6vsCon0 <-MA9_Con0_VS_Con6$three + 0.01 / MA9_Con0_VS_Con6$two + 0.01
f_MA9_fc$MA9_B_Con12vsCon0 <-MA9_Con0_VS_Con12$three + 0.01 / MA9_Con0_VS_Con12$two + 0.01
f_MA9_fc$MA9_B_Con24vsCon0 <-MA9_Con0_VS_Con24$three + 0.01 / MA9_Con0_VS_Con24$two + 0.01
f_MA9_fc$MA9_B_Dox3vsCon0_ <-(MA9_Con0_VS_Dox3$three + 0.01) / (MA9_Con0_VS_Dox3$two + 0.01)
f_MA9_fc$MA9_B_Dox6vsCon0 <-MA9_Con0_VS_Dox6$three + 0.01 / MA9_Con0_VS_Dox6$two + 0.01
f_MA9_fc$MA9_B_Dox12vsCon0 <-MA9_Con0_VS_Dox12$three + 0.01 / MA9_Con0_VS_Dox12$two + 0.01
f_MA9_fc$MA9_B_Dox24vsCon0 <-MA9_Con0_VS_Dox24$three + 0.01 / MA9_Con0_VS_Dox24$two + 0.01


write.csv(f_MA9_fc,file = 'diffbind_nored_ATAC_MA9_fc.csv')























