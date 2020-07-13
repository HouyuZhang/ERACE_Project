setwd("F:/Dropbox/1_Yenlab/1_1_Documents/Lu, Zhang share/2018 Apoptosis project/ATAC-seq/190403 ATAC frag nor count")

####1. Read Diffbind exported results (FDR,normalized counts...) ----
MA9_anno <- read.csv("Diffbind/Diffbind_results_all/MA9Con0_VS_Dox3_500_anno.csv",header = T,row.names = 1)
MA9_Con0_VS_Dox3 <- read.csv("Diffbind/Diffbind_results_all/MA9Con0_VS_Dox3_500.csv",header = T,row.names = 1)
MA9_Con0_VS_Dox6 <- read.csv("Diffbind/Diffbind_results_all/MA9Con0_VS_Dox6_500.csv",header = T,row.names = 1)
MA9_Con0_VS_Dox12 <- read.csv("Diffbind/Diffbind_results_all/MA9Con0_VS_Dox12_500.csv",header = T,row.names = 1)
MA9_Con0_VS_Dox24 <- read.csv("Diffbind/Diffbind_results_all/MA9Con0_VS_Dox24_500.csv",header = T,row.names = 1)

MA9_Con0_VS_Con3 <- read.csv("../../Diffbind/Diffbind_results_all/MA9Con0_VS_Con3_500.csv",header = T,row.names = 1)
MA9_Con0_VS_Con6 <- read.csv("../../Diffbind/Diffbind_results_all/MA9Con0_VS_Con6_500.csv",header = T,row.names = 1)
MA9_Con0_VS_Con12 <- read.csv("../../Diffbind/Diffbind_results_all/MA9Con0_VS_Con12_500.csv",header = T,row.names = 1)
MA9_Con0_VS_Con24 <- read.csv("../../Diffbind/Diffbind_results_all/MA9Con0_VS_Con24_500.csv",header = T,row.names = 1)

filenames <- c("MA9_Con0_VS_Dox3","MA9_Con0_VS_Dox6","MA9_Con0_VS_Dox12","MA9_Con0_VS_Dox24",
               "MA9_Con0_VS_Con3","MA9_Con0_VS_Con6","MA9_Con0_VS_Con12","MA9_Con0_VS_Con24")
###2. Make normalized counts matrix----

#assign basic information to head cols
MA9_matrix <- MA9_Con0_VS_Dox3[,c(1:3)]
MA9_matrix$distanceToTSS <- MA9_anno$distanceToTSS
MA9_matrix$ENSEMBL <- MA9_anno$ENSEMBL
rownames(MA9_matrix) <- paste("Peak_",seq(1,nrow(MA9_matrix)),sep = "")

#assign Con0 data through average
A_con0 <- MA9_Con0_VS_Dox3[,c(1:3)]
B_con0 <- MA9_Con0_VS_Dox3[,c(1:3)]
for (i in filenames){
  tmp <- eval(parse(text = i))
  A_con0 <- cbind(A_con0,tmp$MA9_A_Con0)
  B_con0 <- cbind(B_con0,tmp$MA9_B_Con0)
}

MA9_matrix$MA9_A_Con0 <- apply(A_con0[,4:ncol(A_con0)],1,mean)
MA9_matrix$MA9_B_Con0 <- apply(B_con0[,4:ncol(B_con0)],1,mean)

#assign normalized counts
for (i in filenames){
  tmp <- eval(parse(text = i))
  MA9_matrix <- cbind(MA9_matrix,tmp[,c(14:15)])
}
#assign FDR 
for (i in filenames){
  tmp <- eval(parse(text = i))
  MA9_matrix <- cbind(MA9_matrix,tmp[,11])
}

ncol(MA9_matrix)
write.csv(MA9_matrix[,-c((ncol(MA9_matrix)-7):ncol(MA9_matrix))],
          "190401 replot ATAC enrichment heatmap/MA9/ATAC_MA9_normalized_countMatrix.csv")

###3. Pick out significant peaks (FDR<0.05)----
FDR_sig <- MA9_matrix[,c((ncol(MA9_matrix)-7):ncol(MA9_matrix))] < 0.05
MA9_matrix_sig <- MA9_matrix[rowSums(FDR_sig)!=0,]
nrow(MA9_matrix_sig)
head(MA9_matrix_sig)

write.csv(MA9_matrix_sig[,-c((ncol(MA9_matrix_sig)-7):ncol(MA9_matrix_sig))],
          "190401 replot ATAC enrichment heatmap/MA9/ATAC_MA9_normalized_countMatrix_FDR0.05.csv")

###4. output genes after TSS +-1kb filter (discarded, using before calculated genes(gene body +-2kb) instead)----
if(F){
# with FDR
MA9_matrix_sig_TSS <- MA9_matrix_sig[MA9_matrix_sig$distanceToTSS<1000 & MA9_matrix_sig$distanceToTSS >= -1000,]
write.csv(MA9_matrix_sig_TSS,"190401 replot ATAC enrichment heatmap/MA9/ATAC_MA9_normalized_countMatrix_FDR0.05_TSS.csv")

gene_counts <- with(MA9_matrix_sig_TSS, tapply(MA9_matrix_sig_TSS[,6], list(MA9_matrix_sig_TSS$ENSEMBL), sum))
gene_counts <- as.data.frame(gene_counts)
for (i in seq(7,23)){
  y <- with(MA9_matrix_sig_TSS, tapply(MA9_matrix_sig_TSS[,i], list(MA9_matrix_sig_TSS$ENSEMBL), sum))
  y <- as.data.frame(y)
  gene_counts <- cbind(gene_counts,y)
}
colnames(gene_counts) <- colnames(MA9_matrix_sig_TSS[,c(6:23)])
gene_counts_FDR_TSS <- gene_counts[!is.na(rowSums(gene_counts)),]
write.csv(gene_counts_FDR_TSS ,"190401 replot ATAC enrichment heatmap/MA9/ATAC_MA9_normalized_peak2genes_FDR0.05_TSS.csv")
}

###5. map ATAC to RNA for heatmap----
setwd("F:/Dropbox/1_Yenlab/1_1_Documents/Lu, Zhang share/2018 Apoptosis project/ATAC-seq/190411 re-annotate peak TSS+-1kb/")

plot_correlation <- function(subset="diff",sample="MA9",condition="Dox",removeNA=T,outlist=F,match2=F,kmeans=5){
  #diff/total
  RNA_path<- paste("../190401 replot ATAC enrichment heatmap\\",sample,"\\190401",sample," RNA ",subset," zscore.csv",sep="")
  RNA_genes <- read.csv(RNA_path, header = T,row.names = 1)
  
  ATAC_path <- paste(sample,"_ATAC_genescounts_DESeq_normalized.csv",sep="")
  ATAC_genes <- read.csv(ATAC_path,header = T,row.names = 1)
  
  if(match2){
    out <- pheatmap::pheatmap(RNA_genes,
                                           cluster_cols = F,
                                           show_rownames = F,
                                           clustering_method = "ward.D2",
                                           silent = T
    )    
    genelist <- rownames(RNA_genes)[out$tree_row$order]
  }
  else{genelist <- rownames(RNA_genes)}
  
  ATAC_genes2RNA_genes <- ATAC_genes[as.vector(genelist),]
  #nrow(ATAC_genes2RNA_genes) == length(genelist)
  rownames(ATAC_genes2RNA_genes) <- genelist
  ATAC_genes2RNA_genes[is.na(ATAC_genes2RNA_genes)] <- 0
  #ATAC_genes2RNA_genes <- ATAC_genes
  
  if(condition=="Dox"){
  ATAC_genes2RNA_genes$MA9_Con0 <- (ATAC_genes2RNA_genes[,1] + ATAC_genes2RNA_genes[,10])/2
  ATAC_genes2RNA_genes$MA9_Dox3 <- (ATAC_genes2RNA_genes[,8] + ATAC_genes2RNA_genes[,17])/2
  ATAC_genes2RNA_genes$MA9_Dox6 <- (ATAC_genes2RNA_genes[,9] + ATAC_genes2RNA_genes[,18])/2
  ATAC_genes2RNA_genes$MA9_Dox12 <- (ATAC_genes2RNA_genes[,6] + ATAC_genes2RNA_genes[,15])/2
  ATAC_genes2RNA_genes$MA9_Dox24 <- (ATAC_genes2RNA_genes[,7] + ATAC_genes2RNA_genes[,16])/2
  ATAC_genes2RNA_genes_scaled <- as.data.frame(t(scale(t(ATAC_genes2RNA_genes[,c(19:23)]))))
  }
  else{
  ATAC_genes2RNA_genes$MA9_Con0 <- (ATAC_genes2RNA_genes[,1] + ATAC_genes2RNA_genes[,10])/2
  ATAC_genes2RNA_genes$MA9_Con3 <- (ATAC_genes2RNA_genes[,4] + ATAC_genes2RNA_genes[,13])/2
  ATAC_genes2RNA_genes$MA9_Con6 <- (ATAC_genes2RNA_genes[,5] + ATAC_genes2RNA_genes[,14])/2
  ATAC_genes2RNA_genes$MA9_Con12 <- (ATAC_genes2RNA_genes[,2] + ATAC_genes2RNA_genes[,11])/2
  ATAC_genes2RNA_genes$MA9_Con24 <- (ATAC_genes2RNA_genes[,3] + ATAC_genes2RNA_genes[,12])/2
  ATAC_genes2RNA_genes_scaled <- as.data.frame(t(scale(t(ATAC_genes2RNA_genes[,c(19:23)]))))    
  }
  
  ATAC_genes2RNA_genes_scaled[is.na(ATAC_genes2RNA_genes_scaled)] <- 0
  
  RNA_ATAC <- cbind(RNA_genes[as.vector(genelist),],ATAC_genes2RNA_genes_scaled)
  
  colnames(RNA_ATAC) <- c("RNA_0h","RNA_3h","RNA_6h","RNA_12h","RNA_24h","ATAC_0h","ATAC_3h","ATAC_6h","ATAC_12h","ATAC_24h")
  a <- nrow(RNA_ATAC)
  
  if(removeNA){
    RNA_ATAC <- RNA_ATAC[!rowSums(RNA_ATAC[,6:10])==0,]
    RNA_ATAC <- RNA_ATAC[!rowSums(RNA_ATAC[,1:5])==0,]
  }

  b <- nrow(RNA_ATAC)
  cat(paste(a-b," genes removed because of ATAC-seq unmatched and ",nrow(RNA_ATAC)," genes for plotting...\n",sep = ""))
  
  #colors <- colorRampPalette(("RdYlBu")(100))
  #colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlGn")))(100)
  
  RNA_ATAC <- round(RNA_ATAC,1)
  get_quan <- function(df){
    x <- (abs(as.matrix(df)))
    return(quantile(x,c(.99)[[1]]))
  }
  library(ComplexHeatmap)
  p <- Heatmap(RNA_ATAC,
                          show_row_names = F,
                          cluster_columns = F,
                          clustering_method_rows = "ward.D2",
                          km = kmeans,
                          column_title = paste(sample,"_",condition,"_",subset,sep = ""),
                          name = paste(b,"genes"),
                          col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100)
                          
  )
  ht <- draw(p)
  pdf(paste("0414/Heatmap_ATAC_RNA_",sample,"_",condition,"_",subset,"_",removeNA,"_",kmeans,"clusters.pdf",sep = ""))
  print(ht)
  dev.off()
  
  #get genes for each clusters 
  if (outlist){
  rank <- row_order(ht)
  all_groups <- rowr::cbind.fill(rownames(RNA_ATAC)[rank[[1]]],rownames(RNA_ATAC)[rank[[2]]],fill = NA)
  for (i in 3:kmeans){
    all_groups <- rowr::cbind.fill(all_groups,rownames(RNA_ATAC)[rank[[i]]],fill = NA)
  }
  colnames(all_groups) <- c(paste("cluster",seq(1,kmeans),sep = ""))
  write.csv(all_groups,paste("0414/Heatmap_ATAC_RNA_",sample,"_",condition,"_",subset,"_",removeNA,"_genelist_for_",kmeans,"clusters.csv",sep = ""), na = "",row.names = F)
  }
}
for(i in c("MA9","WT")){
  for (j in c("Dox","Con")){
    for (k in seq(3,10)){
    plot_correlation(subset = "diff",sample = i,condition=j,removeNA=T,outlist=T,match2 = F,kmeans=k)
    }
  }
}
plot_correlation(subset = "diff",sample = "MA9",condition="Dox",removeNA=T,outlist=T,match2 = F,kmeans=5)
#plot_correlation(subset = "diff",sample = "WT",condition="Dox",removeNA=T,outlist=F,match2 = F,kmeans=5)

fun1 <- function(sample="MA9",condition="Dox"){
  #diff/total

  ATAC_path <- paste(sample,"_ATAC_genescounts_DESeq_normalized.csv",sep="")
  ATAC_genes <- read.csv(ATAC_path,header = T,row.names = 1)

  #ATAC_genes2RNA_genes <- ATAC_genes
  
  if(condition=="Dox"){
    ATAC_genes2RNA_genes$MA9_Con0 <- (ATAC_genes2RNA_genes[,1] + ATAC_genes2RNA_genes[,10])/2
    ATAC_genes2RNA_genes$MA9_Dox3 <- (ATAC_genes2RNA_genes[,8] + ATAC_genes2RNA_genes[,17])/2
    ATAC_genes2RNA_genes$MA9_Dox6 <- (ATAC_genes2RNA_genes[,9] + ATAC_genes2RNA_genes[,18])/2
    ATAC_genes2RNA_genes$MA9_Dox12 <- (ATAC_genes2RNA_genes[,6] + ATAC_genes2RNA_genes[,15])/2
    ATAC_genes2RNA_genes$MA9_Dox24 <- (ATAC_genes2RNA_genes[,7] + ATAC_genes2RNA_genes[,16])/2
    ATAC_genes2RNA_genes_scaled <- as.data.frame(t(scale(t(ATAC_genes2RNA_genes[,c(19:23)]))))
  }
  else{
    ATAC_genes2RNA_genes$MA9_Con0 <- (ATAC_genes2RNA_genes[,1] + ATAC_genes2RNA_genes[,10])/2
    ATAC_genes2RNA_genes$MA9_Con3 <- (ATAC_genes2RNA_genes[,4] + ATAC_genes2RNA_genes[,13])/2
    ATAC_genes2RNA_genes$MA9_Con6 <- (ATAC_genes2RNA_genes[,5] + ATAC_genes2RNA_genes[,14])/2
    ATAC_genes2RNA_genes$MA9_Con12 <- (ATAC_genes2RNA_genes[,2] + ATAC_genes2RNA_genes[,11])/2
    ATAC_genes2RNA_genes$MA9_Con24 <- (ATAC_genes2RNA_genes[,3] + ATAC_genes2RNA_genes[,12])/2
    ATAC_genes2RNA_genes_scaled <- as.data.frame(t(scale(t(ATAC_genes2RNA_genes[,c(19:23)]))))    
  }
  
  ATAC_genes2RNA_genes_scaled[is.na(ATAC_genes2RNA_genes_scaled)] <- 0
  write.csv(ATAC_genes2RNA_genes_scaled,paste(sample,"_ATAC_genescounts_DESeq_normalized_",condition,"_z-score.csv",sep = ""), na = "")
}
for(i in c("MA9","WT")){
  for (j in c("Dox","Con")){
      fun1(sample=i,condition=j)
  }
}









