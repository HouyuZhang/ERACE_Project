library(ComplexHeatmap)
library(circlize)

setwd("F:/Dropbox/1_Yenlab/1_1_Documents/Lu, Zhang share/2018 Apoptosis project/ATAC-seq/190411 re-annotate peak TSS+-1kb/")
### Caluculate correlation between RNA and ATAC (Z-score)----
MA9_Dox_Z <- read.csv("MA9_ATAC_genescounts_DESeq_normalized_Dox_z-score.csv",header = T, row.names = 1)
RNA_Z <- read.csv("../190401 replot ATAC enrichment heatmap/MA9/190401MA9 RNA diff zscore.csv",header = T, row.names = 1)
MA9_DoxZ_G <- MA9_DoxZ[rownames(RNA_Z),]
sum(rownames(RNA_Z) != rownames(MA9_DoxZ_G))

MA9_Dox_res <- as.data.frame(matrix(0, nrow = 5, ncol = 5))

for (i in 1:ncol(RNA_Z)){
  for (j in 1:ncol(MA9_DoxZ_G)){
    MA9_Dox_res[ncol(RNA_Z)+1-i,j] <- cor(RNA_Z[,i], MA9_DoxZ_G[,j],method = "pearson", use = "na.or.complete")
  }
}
rownames(MA9_Dox_res) <- rev(c("RNA_0","RNA_3","RNA_6","RNA_12","RNA_24"))
colnames(MA9_Dox_res) <- c("ATAC_0","ATAC_3","ATAC_6","ATAC_12","ATAC_24")

Heatmap(MA9_Dox_res,
        cluster_rows = F,cluster_columns = F,
        row_names_side = "left",
        column_names_side = "bottom",
        name = "MA9_Dox",
        col = circlize::colorRamp2(c(-0.2,0,0.2), c("#4575B4","white", "#D73027")),
        rect_gp = gpar(col = "grey", lwd = 2),
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(sprintf("%.2f", MA9_Dox_res[i, j]), x, y,gp = gpar(fontsize = 12))
        }
)

#WT
WT_DoxZ <- read.csv("WT_ATAC_genescounts_DESeq_normalized_Dox_z-score.csv",header = T, row.names = 1)
RNA_Z <- read.csv("../190401 replot ATAC enrichment heatmap/WT/190401WT RNA diff zscore.csv",header = T, row.names = 1)
WT_DoxZ_G <- WT_DoxZ[rownames(RNA_Z),]
sum(rownames(RNA_Z) != rownames(WT_DoxZ_G))

WT_Dox_res <- as.data.frame(matrix(0, nrow = 5, ncol = 5))

for (i in 1:ncol(RNA_Z)){
  for (j in 1:ncol(WT_DoxZ_G)){
    WT_Dox_res[ncol(RNA_Z)+1-i,j] <- cor(RNA_Z[,i], WT_DoxZ_G[,j],method = "pearson", use = "na.or.complete")
  }
}
rownames(WT_Dox_res) <- rev(c("RNA_0","RNA_3","RNA_6","RNA_12","RNA_24"))
colnames(WT_Dox_res) <- c("ATAC_0","ATAC_3","ATAC_6","ATAC_12","ATAC_24")

Heatmap(WT_Dox_res,
        cluster_rows = F,cluster_columns = F,
        row_names_side = "left",
        column_names_side = "bottom",
        name = "WT_Dox",
        col = circlize::colorRamp2(c(-0.2,0,0.2), c("#4575B4","white", "#D73027")),
        rect_gp = gpar(col = "grey", lwd = 2),
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(sprintf("%.2f", WT_Dox_res[i, j]), x, y,gp = gpar(fontsize = 12))
        }
)

####seperate ----

gene_list <- xlsx::read.xlsx("0414/4group gene.xlsx",1)
MA9.positive <- gene_list$MA9.positive.correlation
MA9.negative <- gene_list$MA9.negative.correlation
WT.positive <- gene_list$WT.positive.correlation
WT.negative <- gene_list$WT.negative.correlation

#MA9_posi
RNA_MA9_Z <- read.csv("../190401 replot ATAC enrichment heatmap/MA9/190401MA9 RNA diff zscore.csv",header = T, row.names = 1)
RR <- RNA_MA9_Z[as.vector(MA9.positive),]
a1 <- read.csv("MA9_ATAC_genescounts_DESeq_normalized_Dox_z-score.csv",header = T, row.names = 1)
AA <- a1[as.vector(MA9.positive),]
sum(rownames(RR) != rownames(AA))

res <- as.data.frame(matrix(0, nrow = 5, ncol = 5))

for (i in 1:ncol(RR)){
  for (j in 1:ncol(AA)){
    res[ncol(RR)+1-i,j] <- cor(RR[,i], AA[,j],method = "pearson", use = "na.or.complete")
  }
}
rownames(res) <- rev(c("RNA_0","RNA_3","RNA_6","RNA_12","RNA_24"))
colnames(res) <- c("ATAC_0","ATAC_3","ATAC_6","ATAC_12","ATAC_24")

Heatmap(res,
        cluster_rows = F,cluster_columns = F,
        row_names_side = "left",
        column_names_side = "bottom",
        name = "MA9_posi",
        col = circlize::colorRamp2(c(-0.2,0,0.2), c("#4575B4","white", "#D73027")),
        rect_gp = gpar(col = "grey", lwd = 2),
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(sprintf("%.2f", res[i, j]), x, y,gp = gpar(fontsize = 12))
        }
)

#MA9_neg
RNA_MA9_Z <- read.csv("../190401 replot ATAC enrichment heatmap/MA9/190401MA9 RNA diff zscore.csv",header = T, row.names = 1)
RR <- RNA_MA9_Z[as.vector(MA9.negative),]
a1 <- read.csv("MA9_ATAC_genescounts_DESeq_normalized_Dox_z-score.csv",header = T, row.names = 1)
AA <- a1[as.vector(MA9.negative),]
sum(rownames(RR) != rownames(AA))

res <- as.data.frame(matrix(0, nrow = 5, ncol = 5))

for (i in 1:ncol(RR)){
  for (j in 1:ncol(AA)){
    res[ncol(RR)+1-i,j] <- cor(RR[,i], AA[,j],method = "pearson", use = "na.or.complete")
  }
}
rownames(res) <- rev(c("RNA_0","RNA_3","RNA_6","RNA_12","RNA_24"))
colnames(res) <- c("ATAC_0","ATAC_3","ATAC_6","ATAC_12","ATAC_24")

Heatmap(res,
        cluster_rows = F,cluster_columns = F,
        row_names_side = "left",
        column_names_side = "bottom",
        name = "MA9_neg",
        col = circlize::colorRamp2(c(-0.2,0,0.2), c("#4575B4","white", "#D73027")),
        rect_gp = gpar(col = "grey", lwd = 2),
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(sprintf("%.2f", res[i, j]), x, y,gp = gpar(fontsize = 12))
        }
)

#WT_posi
RNA_WT_Z <- read.csv("../190401 replot ATAC enrichment heatmap/WT/190401WT RNA diff zscore.csv",header = T, row.names = 1)
RR <- RNA_WT_Z[as.vector(WT.positive),]
a1 <- read.csv("WT_ATAC_genescounts_DESeq_normalized_Dox_z-score.csv",header = T, row.names = 1)
AA <- a1[as.vector(WT.positive),]
sum(rownames(RR) != rownames(AA))

res <- as.data.frame(matrix(0, nrow = 5, ncol = 5))

for (i in 1:ncol(RR)){
  for (j in 1:ncol(AA)){
    res[ncol(RR)+1-i,j] <- cor(RR[,i], AA[,j],method = "pearson", use = "na.or.complete")
  }
}
rownames(res) <- rev(c("RNA_0","RNA_3","RNA_6","RNA_12","RNA_24"))
colnames(res) <- c("ATAC_0","ATAC_3","ATAC_6","ATAC_12","ATAC_24")

Heatmap(res,
        cluster_rows = F,cluster_columns = F,
        row_names_side = "left",
        column_names_side = "bottom",
        name = "WT_posi",
        col = circlize::colorRamp2(c(-0.2,0,0.2), c("#4575B4","white", "#D73027")),
        rect_gp = gpar(col = "grey", lwd = 2),
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(sprintf("%.2f", res[i, j]), x, y,gp = gpar(fontsize = 12))
        }
)

#WT_neg
RNA_WT_Z <- read.csv("../190401 replot ATAC enrichment heatmap/WT/190401WT RNA diff zscore.csv",header = T, row.names = 1)
RR <- RNA_WT_Z[as.vector(WT.negative),]
a1 <- read.csv("WT_ATAC_genescounts_DESeq_normalized_Dox_z-score.csv",header = T, row.names = 1)
AA <- a1[as.vector(WT.negative),]
sum(rownames(RR) != rownames(AA))

res <- as.data.frame(matrix(0, nrow = 5, ncol = 5))

for (i in 1:ncol(RR)){
  for (j in 1:ncol(AA)){
    res[ncol(RR)+1-i,j] <- cor(RR[,i], AA[,j],method = "pearson", use = "na.or.complete")
  }
}
rownames(res) <- rev(c("RNA_0","RNA_3","RNA_6","RNA_12","RNA_24"))
colnames(res) <- c("ATAC_0","ATAC_3","ATAC_6","ATAC_12","ATAC_24")

Heatmap(res,
        cluster_rows = F,cluster_columns = F,
        row_names_side = "left",
        column_names_side = "bottom",
        name = "WT_neg",
        col = circlize::colorRamp2(c(-0.2,0,0.2), c("#4575B4","white", "#D73027")),
        rect_gp = gpar(col = "grey", lwd = 2),
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(sprintf("%.2f", res[i, j]), x, y,gp = gpar(fontsize = 12))
        }
)

