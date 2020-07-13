MA9_neg_com <- read.table("F:/Dropbox/1_Yenlab/1_7_sz_download/MA9_neg_TSS.txt",sep = "\t",header = F)

colnames(MA9_neg_com) <- seq(-1000,999)
pdf("MA9_neg_com.pdf")
pheatmap::pheatmap(MA9_neg_com,
                   color = rainbow(12),
                   cluster_cols = F,
                   cluster_rows = T,
                   show_rownames = F,
                   show_colnames = F)
dev.off()

MA9_posi_com <- read.table("F:/Dropbox/1_Yenlab/1_7_sz_download/test.txt",sep = "\t",header = F)

pdf("MA9_posi_composit.pdf")
pheatmap::pheatmap(MA9_posi_com,
                   color = rainbow(12),
                   cluster_cols = F,
                   cluster_rows = T,
                   show_rownames = F,
                   show_colnames = F)
dev.off()
