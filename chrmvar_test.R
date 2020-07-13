###chromVAR test###
#directory <- "D:/Dropbox/pre-existing/Footprint/test_zhang_file/rawdata/"
setwd("E:/Dropbox/1_Yenlab/1_1_Documents/Lu, Zhang share/2018 Apoptosis project/ATAC-seq/Footprint analysis/chromVAR")

library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2016)
BiocManager::install("rtracklayer", version = "3.8")
register(SerialParam())

####get count table in ATAC-seq data####
peakfile <- c("Zhang_ATAC_sort_filter.narrowPeak")
peaks <- readNarrowpeaks(peakfile,width = 500,non_overlapping = TRUE)

B_bamfile <- c("B_ATAC_demit.bam")
T_bamfile <- c("CD4_ATAC_demit.bam")
fragment_counts <- getCounts(c(B_bamfile,T_bamfile), 
                             peaks, 
                             paired =  FALSE, 
                             by_rg = FALSE, 
                             format = "bam", 
                             colData = S4Vectors::DataFrame(celltype = c("B","T")))


####Adding GC content####
head(fragment_counts)
fragment_counts <- addGCBias(fragment_counts, 
                             genome = BSgenome.Mmusculus.UCSC.mm10)
head(fragment_counts)

####Filtering peaks####
counts_filtered <- filterPeaks(fragment_counts,non_overlapping = TRUE)

####get motifs####
motifs_mm <- getJasparMotifs(species = "Mus musculus")

library(universalmotif)
cisbp <- read_cisbp(system.file("extdata", "cisbp.txt",package = "universalmotif"))
cisbp <- read_cisbp("cisbp_mm_experiment_validated_502.txt")
motifs_mm <- convert_motifs(cisbp, "TFBSTools-PFMatrix")

motif_ix <- matchMotifs(motifs_mm,counts_filtered,
                        genome = BSgenome.Mmusculus.UCSC.mm10)

kmer_ix <- matchKmers(6, counts_filtered, 
                      genome = BSgenome.Mmusculus.UCSC.mm10)

####compute deviations####
dev <- computeDeviations(object = counts_filtered,annotations = motif_ix)

####background peaks####
bg <- getBackgroundPeaks(object = counts_filtered)

dev <- computeDeviations(object = counts_filtered,annotations = motif_ix,
                         background_peaks = bg)

####variability####
variability <- computeVariability(dev)
plotVariability(variability,use_plotly = FALSE,n=5)

####visualizing deviations####
#tsne_results <- deviationsTsne(dev,threshold = 1.5)
