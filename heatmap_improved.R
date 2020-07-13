library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(xlsx)

setwd("./")
######## 0. Pathway code (from "180415go_list.xlsx")########
"
1: methylation+chromatin organization
2: methylation
3: chromatin organization
4: mitotic cell cycle
5: autophagy
6: apoptosis main
7: apoptosis
8: pluripotency markers
"
"new code (from 180601go_list.xlsx)
1: methylation
2: apoptosis
3: mitotic cell cycle
4: pluripotency markers
5: cell death
6: activation of innate immune response + T cell receptor signaling pathway
7: apoptosis + cell death
8: T cell receptor signaling pathway
"

######## 1. define functions ########

#### 1.1 sort your data frame ####
mysort <- function(df,long = T){
  if (isTRUE(long)){
    return(df[order(-df[,1],-df[,2],-df[,3],-df[,4],-df[,5],-df[,6],-df[,7],-df[,8]),])
  }
  else{
    return(df[order(-df[,1],-df[,2],-df[,3],-df[,4]),])
  }
}
rev_sort <- function(df){
  return(df[order(-df[,5],-df[,6],-df[,7],-df[,8]),])
}

##### 1.2 get upper quantile ####
get_quan <- function(df){
  x <- (abs(as.matrix(df)))
  x <- x[which(x!= 0)]
  return(quantile(x,c(.90)[[1]]))
}

##### 1.3 substract method ####
substract <- function(df,long = T){
  if (isTRUE(long)){
    df$MA9_3h <- df[,1]-df[,5]
    df$MA9_6h <- df[,2]-df[,6]
    df$MA9_12h <- df[,3]-df[,7]
    df$MA9_24h <- df[,4]-df[,8]
    df$WT_3h <- df[,9]-df[,13]
    df$WT_6h <- df[,10]-df[,14]
    df$WT_12h <- df[,11]-df[,15]
    df$WT_24h <- df[,12]-df[,16]  
    return(df[,17:ncol(df)])
  }
  else{
    df$h3 <- df[,1]-df[,5]
    df$h6 <- df[,2]-df[,6]
    df$h12 <- df[,3]-df[,7]
    df$h24 <- df[,4]-df[,8]
    return(df[,9:ncol(df)])
  }
}
##### 1.4 remove zeros (discarded)####
remove_zero <- function(df){
  zero_line <- c()
  num <- 0
  for (line in 1:nrow(df))
    if (sum(df[line,])==0){
      num <- num +1
      zero_line[num] <- line
    }
  return(df[-zero_line,])
}



##### 1.5 correlation####

cal_correlation <- function(df){
  h3 <- cor(df[[1]],df[[5]],method="pearson")
  h6 <- cor(df[[2]],df[[6]],method="pearson")
  h12 <- cor(df[[3]],df[[7]],method="pearson")
  h24 <- cor(df[[4]],df[[8]],method="pearson")
  tmp <- c(h3,h6,h12,h24)
  return(tmp)
}

######## 2. Run your own data ########

#### 2.1 Run your RNA-seq data ####
RNA_run <- function(sheet = 1, MA9 = T){
  #read your RNA-seq FC data (i.e. a matrix include all MA9 and WT samples)
  if (isTRUE(MA9)){
  RNA <- read.table("RNA_MA_WT_FINALLIST.txt",header = T,row.names = 1)[1:8]
  }
  else{
  RNA <- read.table("RNA_MA_WT_FINALLIST.txt",header = T,row.names = 1)[9:16]
  }
  #read your pathway data
  pathway_list <- read.xlsx("180601go_list_motif.xlsx",sheet)
  #get values from RNA of your current pathway genelist
  df <- na.omit(RNA[as.vector(pathway_list[,1]),])
  
  #df <- substract(df,T)
  
  #if long = T is assigned, it will sort all 8 columns,else sort head 4 columns. IF you want sort the tail 4 columns, please use :
  #df  <- rev_sort(df)
  df <- mysort(df, long = T)
  
  cat("There are",nrow(df),"genes plotted\nThe upper quantile is",get_quan(df),"\n")
  
  Heatmap(df,cluster_rows = F,
          cluster_columns = F,
          row_names_side = "left",
          column_names_side = "top",
          col = colorRamp2(c(-2,0,2), c("blue","white", "red"))
          #col = colorRamp2(c(-get_quan(df),0,get_quan(df)), c("blue","white", "red"))
          )
  }

RNA_run(5)


#### 2.2 Run your ATAC-seq data ####

ATAC_run <- function(sheet = 1, map_RNA = T,MA9 = T){
  if (isTRUE(MA9)){
  #read your ATAC-seq FC data (i.e. a matrix include all MA9 and WT samples)
  ATAC <- read.table("ATAC_frag_MA9_WT_FINAL.txt",header = T,row.names = 1)[1:8]}
  else{
  ATAC <- read.table("ATAC_frag_MA9_WT_FINAL.txt",header = T,row.names = 1)[9:16]
  }
  #read your pathway data
  pathway_list <- read.xlsx("180601go_list_motif.xlsx",sheet)
  
  if (isTRUE(map_RNA)){
    if (isTRUE(MA9)){
      RNA <- read.table("RNA_MA_WT_FINALLIST.txt",header = T,row.names = 1)[1:8]
    }
    else{
      RNA <- read.table("RNA_MA_WT_FINALLIST.txt",header = T,row.names = 1)[9:16]
    }
    df <- na.omit(RNA[as.vector(pathway_list[,1]),])
    #df <- substract(df,long = T)
    #df  <- rev_sort(df)
    df <- mysort(df, long = T)
    df <- na.omit(ATAC[as.vector(rownames(df)),])
    #df <- substract(df ,long = T)
    cat("There are",nrow(df),"genes plotted\nThe upper quantile is",get_quan(df),"\n")
    
    Heatmap(df,cluster_rows = F,
            cluster_columns = F,
            row_names_side = "left",
            column_names_side = "top",
            col = colorRamp2(c(-4,0,4), c("blue","white", "red"))
            #col = colorRamp2(c(-get_quan(df),0,get_quan(df)), c("blue","white", "red"))
    )
  }
  else{
  #get values from ATAC of your current pathway genelist
  df <- na.omit(ATAC[as.vector(pathway_list[,1]),])
  #df <- substract(df,T)

  #if long = T is assigned, it will sort all 8 columns,else sort head 4 columns. If you want sort the tail 4 columns, please use:
  #df  <- rev_sort(df)
  df <- mysort(df, long = F)
  cat("There are",nrow(df),"genes plotted\nThe upper quantile is",get_quan(df),"\n")
  
  Heatmap(df,cluster_rows = F,
          cluster_columns = F,
          row_names_side = "left",
          column_names_side = "top",
          #col = colorRamp2(c(-10,0,10), c("deepskyblue2","black", "yellow"))
          col = colorRamp2(c(-get_quan(df),0,get_quan(df)), c("blue","white", "red"))
  )
  }
}
# the genes for map_RNA = F are usually more than map_RNA = T, because the input file include all genes, which many are zero,
# so you may PAY ATTENTION!
ATAC_run(2,map_RNA = F) 

ATAC_run(2)

#### 2.3 Run your microarray data ####

micro_run <- function(sheet = 1, map_RNA = T){
  #read your micro FC data (i.e. a matrix include all MA9 and WT samples)
  micro <- read.table("MICRO_LIST.txt",header = T,row.names = 1)
  #read your pathway data
  pathway_list <- read.xlsx("180415go_list.xlsx",sheet)
  
  if (isTRUE(map_RNA)){
    RNA <- read.table("RNA_MA_WT_FINALLIST.txt",header = T,row.names = 1)
    RNA_raw_df <- na.omit(RNA[as.vector(pathway_list[,1]),])
    RNA_substracted_df <- substract(RNA_raw_df,long = T)
    #RNA_sorted_substracted_df  <- rev_sort(RNA_substracted_df)
    RNA_sorted_substracted_df <- mysort(RNA_substracted_df, long = F)
    micro_mapped_raw_df <- na.omit(micro[as.vector(rownames(RNA_sorted_substracted_df)),])
    micro_sorted_mapped_raw_df <- micro_mapped_raw_df[order(- micro_mapped_raw_df[,1],- micro_mapped_raw_df[,2]),]
    cat("There are",nrow(micro_sorted_mapped_raw_df),"genes plotted\nThe upper quantile is",get_quan(micro_sorted_mapped_raw_df),"\n")
    
    Heatmap(micro_sorted_mapped_raw_df,cluster_rows = F,
            cluster_columns = F,
            row_names_side = "left",
            column_names_side = "top",
            #col = colorRamp2(c(-10,0,10), c("deepskyblue2","black", "yellow"))
            col = colorRamp2(c(-get_quan(micro_sorted_mapped_raw_df),0,get_quan(micro_sorted_mapped_raw_df)), c("blue","white", "red"))
    )
  }
  else{
    #get values from microarray of your current pathway genelist
    micro_raw_df <- na.omit(micro[as.vector(pathway_list[,1]),])
    micro_sorted_mapped_raw_df <- micro_raw_df[order(- micro_raw_df[,1],- micro_raw_df[,2]),]
    #micro_sorted_mapped_raw_df <- round(micro_sorted_mapped_raw_df,4)
    cat("There are",nrow(micro_sorted_mapped_raw_df),"genes plotted\nThe upper quantile is",get_quan(micro_sorted_mapped_raw_df),"\n")
      
    Heatmap(micro_sorted_mapped_raw_df,cluster_rows = F,
            cluster_columns = F,
            row_names_side = "left",
            column_names_side = "top",
            #col = colorRamp2(c(-0.4,0,0.4), c("blue","white", "red"))
            col = colorRamp2(c(-get_quan(micro_sorted_mapped_raw_df),0,get_quan(micro_sorted_mapped_raw_df)), c("blue","white", "red"))
    )
  }
}

micro_run(1,map_RNA = T)
micro_run(2,map_RNA = F)


#### 2.4 confidece interval heatmap ####
#RNA <- read.table("RNA_MA_WT_FINALLIST.txt",header = T,row.names = 1)
RNA_filter_run <- function(sheet = 1,log = T, MA9 = T){
  
  pathway_list <- read.xlsx("E:/Lu_ATAC-seq/fragment_genecount/180601go_list_motif.xlsx",sheet)
  #df <- na.omit(RNA[as.vector(pathway_list[,1]),])
  if (isTRUE(MA9)){
    RNA_MA9_fc <-read.csv('E:/Lu_RNA-seq/corr/DESeq normalized count/MA9_normalized_fc.csv',header = T,row.names = 1)
    RNA_MA9_df <- na.omit(RNA_MA9_fc[as.vector(pathway_list[,1]),])
    RNA_MA9_df$MA9_h3 <- 0
    RNA_MA9_df$MA9_h6 <- 0
    RNA_MA9_df$MA9_h12 <- 0
    RNA_MA9_df$MA9_h24 <- 0
    
    for (sample in 1:4){
      for (irow in 1:nrow(RNA_MA9_df)){
        rep1 <- RNA_MA9_df[[4 + sample]][irow]/RNA_MA9_df[[sample]][irow]
        rep2 <- RNA_MA9_df[[12 + sample]][irow]/RNA_MA9_df[[8 + sample]][irow]
        change <- 0 
        if (rep1 != rep2){
          if (t.test(c(rep1,rep2))$p.value < 0.05){
            if (!data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
              if(isTRUE(log)){change <- log10((rep1 + rep2)/2)}
              else{change <- (rep1 + rep2)/2}
              
            }
          }
        }
        #print(change)
        RNA_MA9_df[row.names(RNA_MA9_df)[irow],][[16 + sample]] <- change
      }
    }
    RNA_MA9_df <- RNA_MA9_df[17:20]
    RNA_MA9_df <- mysort(RNA_MA9_df,long = F)
    cat("There are",nrow(RNA_MA9_df),"genes plotted\nThe upper quantile is",get_quan(RNA_MA9_df),"\n")
    
    if(isTRUE(log)){
      Heatmap(RNA_MA9_df,cluster_rows = F,
              cluster_columns = F,
              row_names_side = "left",
              column_names_side = "top",
              col = colorRamp2(c(-get_quan(RNA_MA9_df),0,get_quan(RNA_MA9_df)), c("blue","white", "red"))
              #col = colorRamp2(c(-0.4,0,0.4), c("blue","white", "red"))
      )
    }
    else{
      Heatmap(RNA_MA9_df,cluster_rows = F,
              cluster_columns = F,
              row_names_side = "left",
              column_names_side = "top",
              #col = colorRamp2(c(0,get_quan(RNA_MA9_df)), c("white", "red"))
              col = colorRamp2(c(0,2), c("white", "red"))
      )
    }
  }


  else{
    RNA_WT_fc <-read.csv('E:/Lu_RNA-seq/corr/DESeq normalized count/WT_normalized_fc.csv',header = T,row.names = 1) 
    RNA_WT_df <- na.omit(RNA_WT_fc[as.vector(pathway_list[,1]),])
    RNA_WT_df$WT_h3 <- 0
    RNA_WT_df$WT_h6 <- 0
    RNA_WT_df$WT_h12 <- 0
    RNA_WT_df$WT_h24 <- 0
    
    for (sample in 1:4){
      for (irow in 1:nrow(RNA_WT_df)){
        rep1 <- RNA_WT_df[[4 + sample]][irow]/RNA_WT_df[[sample]][irow]
        rep2 <- RNA_WT_df[[12 + sample]][irow]/RNA_WT_df[[8 + sample]][irow]
        change <- 0 
        if (rep1 != rep2){
          if (t.test(c(rep1,rep2))$p.value < 0.05){
            if (!data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
              if(isTRUE(log)){change <- log10((rep1 + rep2)/2)}
              else{change <- (rep1 + rep2)/2}
            }
          }
        }
        RNA_WT_df[row.names(RNA_WT_df)[irow],][[16 + sample]] <- change
      }
    }
    RNA_WT_df <- RNA_WT_df[17:20]
    RNA_WT_df <- mysort(RNA_WT_df,long = F)
    cat("There are",nrow(RNA_WT_df),"genes plotted\nThe upper quantile is",get_quan(RNA_WT_df),"\n")
    
    if(isTRUE(log)){
      Heatmap(RNA_WT_df,cluster_rows = F,
              cluster_columns = F,
              row_names_side = "left",
              column_names_side = "top",
              col = colorRamp2(c(-get_quan(RNA_WT_df),0,get_quan(RNA_WT_df)), c("blue","white", "red"))
              #col = colorRamp2(c(-0.4,0,0.4), c("blue","white", "red"))
      )
    }
    else{
      Heatmap(RNA_WT_df,cluster_rows = F,
              cluster_columns = F,
              row_names_side = "left",
              column_names_side = "top",
              #col = colorRamp2(c(0,get_quan(RNA_WT_df)), c("white", "red"))
              col = colorRamp2(c(0,2), c("white", "red"))
      )
  }
  #RNA_df <- cbind(RNA_MA9_df[17:20],RNA_WT_df[17:20])
  #RNA_df <- mysort(RNA_df,long = T)
}
}



WT_fc <- WT[19:ncol(WT)]
MA9_fc <- MA9[19:ncol(MA9)]

#***************************************
ATAC_filter_run <- function(sheet = 1,log = T, MA9 = T){
pathway_list <- read.xlsx("E:/Lu_ATAC-seq/fragment_genecount/180601go_list_motif.xlsx",sheet)
#df <- na.omit(RNA[as.vector(pathway_list[,1]),])

if (isTRUE(MA9)){
  RNA_MA9_fc <-read.csv('E:/Lu_RNA-seq/corr/DESeq normalized count/MA9_normalized_fc.csv',header = T,row.names = 1)
  RNA_MA9_df <- na.omit(RNA_MA9_fc[as.vector(pathway_list[,1]),])
  RNA_MA9_df$MA9_h3 <- 0
  RNA_MA9_df$MA9_h6 <- 0
  RNA_MA9_df$MA9_h12 <- 0
  RNA_MA9_df$MA9_h24 <- 0
  
  for (sample in 1:4){
    for (irow in 1:nrow(RNA_MA9_df)){
      rep1 <- RNA_MA9_df[[4 + sample]][irow]/RNA_MA9_df[[sample]][irow]
      rep2 <- RNA_MA9_df[[12 + sample]][irow]/RNA_MA9_df[[8 + sample]][irow]
      change <- 0 
      if (rep1 != rep2){
        if (t.test(c(rep1,rep2))$p.value < 0.05){
          if (!data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
            if(isTRUE(log)){change <- log10((rep1 + rep2)/2)}
            else{change <- (rep1 + rep2)/2}
            
          }
        }
      }
      #print(change)
      RNA_MA9_df[row.names(RNA_MA9_df)[irow],][[16 + sample]] <- change
    }
  }
  
  RNA_MA9_df <- RNA_MA9_df[17:20]
  RNA_MA9_df <- mysort(RNA_MA9_df,long = F)
  MA9_df <- na.omit(MA9_fc[as.vector(rownames(RNA_MA9_df)),])
  
  MA9_df$MA9_h3 <- 0
  MA9_df$MA9_h6 <- 0
  MA9_df$MA9_h12 <- 0
  MA9_df$MA9_h24 <- 0
  
  for (sample in 1:4){
    for (irow in 1:nrow(MA9_df)){
      rep1 <- MA9_df[[4 + sample]][irow]/MA9_df[[sample]][irow]
      rep2 <- MA9_df[[12 + sample]][irow]/MA9_df[[8 + sample]][irow]
      change <- 0 
      if (rep1 != rep2){
        if (t.test(c(rep1,rep2))$p.value < 0.05){
          if (!data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
            if(isTRUE(log)){change <- log10((rep1 + rep2)/2)}
            else{change <- (rep1 + rep2)/2}
          }
        }
      }
      #print(change)
      MA9_df[row.names(MA9_df)[irow],][[16 + sample]] <- change
    }
  }
  
  MA9_df <- MA9_df[17:20]
  cat("There are",nrow(MA9_df),"genes plotted\nThe upper quantile is",get_quan(MA9_df),"\n")
  if(isTRUE(log)){
    Heatmap(MA9_df,cluster_rows = F,
            cluster_columns = F,
            row_names_side = "left",
            column_names_side = "top",
            #col = colorRamp2(c(-get_quan(MA9_df),0,get_quan(MA9_df)), c("blue","white", "red"))
            col = colorRamp2(c(-1,0,1), c("blue","white", "red"))
    )
  }
  else{
    Heatmap(MA9_df,cluster_rows = F,
            cluster_columns = F,
            row_names_side = "left",
            column_names_side = "top",
            #col = colorRamp2(c(0,get_quan(MA9_df)), c("white", "red"))
            col = colorRamp2(c(0,2), c("white", "red"))
    )
  }
}
#####*************************************************************
else{
  RNA_WT_fc <-read.csv('E:/Lu_RNA-seq/corr/DESeq normalized count/WT_normalized_fc.csv',header = T,row.names = 1) 
  RNA_WT_df <- na.omit(RNA_WT_fc[as.vector(pathway_list[,1]),])
  RNA_WT_df$WT_h3 <- 0
  RNA_WT_df$WT_h6 <- 0
  RNA_WT_df$WT_h12 <- 0
  RNA_WT_df$WT_h24 <- 0

  for (sample in 1:4){
    for (irow in 1:nrow(RNA_WT_df)){
      rep1 <- RNA_WT_df[[4 + sample]][irow]/RNA_WT_df[[sample]][irow]
      rep2 <- RNA_WT_df[[12 + sample]][irow]/RNA_WT_df[[8 + sample]][irow]
      change <- 0 
      if (rep1 != rep2){
        if (t.test(c(rep1,rep2))$p.value < 0.05){
          if (!data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
            if(isTRUE(log)){change <- log10((rep1 + rep2)/2)}
            else{change <- (rep1 + rep2)/2}
          }
        }
      }
      RNA_WT_df[row.names(RNA_WT_df)[irow],][[16 + sample]] <- change
    }
  }
  RNA_WT_df <- RNA_WT_df[17:20]
  RNA_WT_df <- mysort(RNA_WT_df,long = F)
  
  WT_df <- na.omit(WT_fc[as.vector(rownames(RNA_WT_df)),])
  WT_df$WT_h3 <- 0
  WT_df$WT_h6 <- 0
  WT_df$WT_h12 <- 0
  WT_df$WT_h24 <- 0
  for (sample in 1:4){
    for (irow in 1:nrow(WT_df)){
      rep1 <- WT_df[[4 + sample]][irow]/WT_df[[sample]][irow]
      rep2 <- WT_df[[12 + sample]][irow]/WT_df[[8 + sample]][irow]
      change <- 0 
      if (rep1 != rep2){
        if (t.test(c(rep1,rep2))$p.value < 0.05){
          if (!data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
            if(isTRUE(log)){change <- log10((rep1 + rep2)/2)}
            else{change <- (rep1 + rep2)/2}
          }
        }
      }
      WT_df[row.names(WT_df)[irow],][[16 + sample]] <- change
    }
  }
  
  WT_df <- WT_df[17:20]
  cat("There are",nrow(WT_df),"genes plotted\nThe upper quantile is",get_quan(WT_df),"\n")
  if(isTRUE(log)){
    Heatmap(WT_df,cluster_rows = F,
            cluster_columns = F,
            row_names_side = "left",
            column_names_side = "top",
            #col = colorRamp2(c(-get_quan(WT_df),0,get_quan(WT_df)), c("blue","white", "red"))
            col = colorRamp2(c(-1,0,1), c("blue","white", "red"))
    )
  }
  else{
    Heatmap(WT_df,cluster_rows = F,
            cluster_columns = F,
            row_names_side = "left",
            column_names_side = "top",
            #col = colorRamp2(c(0,get_quan(WT_df)), c("white", "red"))
            col = colorRamp2(c(0,2), c("white", "red"))
    )
  }
  
}
#####
#RNA_df <- cbind(RNA_MA9_df[17:20],RNA_WT_df[17:20])
#df <- rev_sort(df)
#df <- remove_zero(df)
}

RNA <- read.table("E:/Lu_ATAC-seq/fragment_genecount/ALL_RNA_MA_WT_FINALLIST.txt",header = T, row.names = 1)
ATAC <- read.table("E:/Lu_ATAC-seq/2. remove_bl/ATAC_MA9_WT_FINALLIST_all.txt",header = T, row.names = 1)

pathway_list <- read.xlsx("E:/Lu_ATAC-seq/fragment_genecount/180601go_list_motif.xlsx",7)
RNA_df <- na.omit(RNA[as.vector(pathway_list[,1]),])

#####
RNA_MA9_fc <-read.csv('E:/Lu_RNA-seq/corr/DESeq normalized count/MA9_normalized_fc.csv',header = T,row.names = 1)
RNA_MA9_df <- na.omit(RNA_MA9_fc[as.vector(pathway_list[,1]),])
RNA_MA9_df$MA9_h3 <- 0
RNA_MA9_df$MA9_h6 <- 0
RNA_MA9_df$MA9_h12 <- 0
RNA_MA9_df$MA9_h24 <- 0

for (sample in 1:4){
  for (irow in 1:nrow(RNA_MA9_df)){
    rep1 <- RNA_MA9_df[[4 + sample]][irow]/RNA_MA9_df[[sample]][irow]
    rep2 <- RNA_MA9_df[[12 + sample]][irow]/RNA_MA9_df[[8 + sample]][irow]
    change <- 0 
    if (rep1 != rep2){
      if (t.test(c(rep1,rep2))$p.value < 0.05){
        if (data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
          RNA_df[rownames(RNA_MA9_df)[irow],][[sample]] <- 0
          RNA_df[rownames(RNA_MA9_df)[irow],][[4 + sample]] <- 0
        }
      }
      else{
        RNA_df[rownames(RNA_MA9_df)[irow],][[sample]] <- 0
        RNA_df[rownames(RNA_MA9_df)[irow],][[4 + sample]] <- 0
      }
    }
    else{
      RNA_df[rownames(RNA_MA9_df)[irow],][[sample]] <- 0
      RNA_df[rownames(RNA_MA9_df)[irow],][[4 + sample]] <- 0
    }
  }
}



RNA_WT_fc <-read.csv('E:/Lu_RNA-seq/corr/DESeq normalized count/WT_normalized_fc.csv',header = T,row.names = 1)
RNA_WT_df <- na.omit(RNA_WT_fc[as.vector(pathway_list[,1]),])
RNA_WT_df$WT_h3 <- 0
RNA_WT_df$WT_h6 <- 0
RNA_WT_df$WT_h12 <- 0
RNA_WT_df$WT_h24 <- 0

for (sample in 1:4){
  for (irow in 1:nrow(RNA_WT_df)){
    rep1 <- RNA_WT_df[[4 + sample]][irow]/RNA_WT_df[[sample]][irow]
    rep2 <- RNA_WT_df[[12 + sample]][irow]/RNA_WT_df[[8 + sample]][irow]
    change <- 0 
    if (rep1 != rep2){
      if (t.test(c(rep1,rep2))$p.value < 0.05){
        if (data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
          RNA_df[rownames(RNA_WT_df)[irow],][[8 + sample]] <- 0
          RNA_df[rownames(RNA_WT_df)[irow],][[12 + sample]] <- 0
        }
      }
      else{
        RNA_df[rownames(RNA_WT_df)[irow],][[8 + sample]] <- 0
        RNA_df[rownames(RNA_WT_df)[irow],][[12 + sample]] <- 0
      }
    }
    else{
      RNA_df[rownames(RNA_WT_df)[irow],][[8 + sample]] <- 0
      RNA_df[rownames(RNA_WT_df)[irow],][[12 + sample]] <- 0
    }
}
}
#************************************************
res_MA9_df <- RNA_df[1:4]
res_WT_df <- RNA_df[9:12]
res_MA9_df <- mysort(res_MA9_df,long = F)
res_WT_df <- mysort(res_WT_df,long = F)




ATAC_MA9_reads_df <- na.omit(ATAC[as.vector(rownames(res_MA9_df)),])[1:8]

ATAC_MA9_fc <-read.csv('E:/Lu_ATAC-seq/2. remove_bl/MA9_normalized_fc.csv',header = T,row.names = 1)
ATAC_MA9_df <- na.omit(ATAC_MA9_fc[as.vector(pathway_list[,1]),])
ATAC_MA9_df$MA9_h3 <- 0
ATAC_MA9_df$MA9_h6 <- 0
ATAC_MA9_df$MA9_h12 <- 0
ATAC_MA9_df$MA9_h24 <- 0

for (sample in 1:4){
  for (irow in 1:nrow(ATAC_MA9_df)){
    rep1 <- ATAC_MA9_df[[4 + sample]][irow]/ATAC_MA9_df[[sample]][irow]
    rep2 <- ATAC_MA9_df[[12 + sample]][irow]/ATAC_MA9_df[[8 + sample]][irow]

    if (rep1 != rep2){
      if (t.test(c(rep1,rep2))$p.value < 0.05){
        if (data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
          ATAC_MA9_reads_df[rownames(ATAC_MA9_df)[irow],][[sample]] <- 0
          ATAC_MA9_reads_df[rownames(ATAC_MA9_df)[irow],][[4 + sample]] <- 0
        }
      }
      else{
        ATAC_MA9_reads_df[rownames(ATAC_MA9_df)[irow],][[sample]] <- 0
        ATAC_MA9_reads_df[rownames(ATAC_MA9_df)[irow],][[4 + sample]] <- 0
      }
    }
    else{
      ATAC_MA9_reads_df[rownames(ATAC_MA9_df)[irow],][[sample]] <- 0
      ATAC_MA9_reads_df[rownames(ATAC_MA9_df)[irow],][[4 + sample]] <- 0
    }
  }
}

ATAC_WT_reads_df <- na.omit(ATAC[as.vector(rownames(res_WT_df)),])[9:16]

ATAC_WT_fc <-read.csv('E:/Lu_ATAC-seq/2. remove_bl/WT_normalized_fc.csv',header = T,row.names = 1)
ATAC_WT_df <- na.omit(ATAC_WT_fc[as.vector(pathway_list[,1]),])
ATAC_WT_df$WT_h3 <- 0
ATAC_WT_df$WT_h6 <- 0
ATAC_WT_df$WT_h12 <- 0
ATAC_WT_df$WT_h24 <- 0

for (sample in 1:4){
  for (irow in 1:nrow(ATAC_WT_df)){
    rep1 <- ATAC_WT_df[[4 + sample]][irow]/ATAC_WT_df[[sample]][irow]
    rep2 <- ATAC_WT_df[[12 + sample]][irow]/ATAC_WT_df[[8 + sample]][irow]
    change <- 0 
    if (rep1 != rep2){
      if (t.test(c(rep1,rep2))$p.value < 0.05){
        if (data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
          ATAC_WT_reads_df[rownames(ATAC_WT_df)[irow],][[sample]] <- 0
          ATAC_WT_reads_df[rownames(ATAC_WT_df)[irow],][[4 + sample]] <- 0
        }
      }
      else{
        ATAC_WT_reads_df[rownames(ATAC_WT_df)[irow],][[sample]] <- 0
        ATAC_WT_reads_df[rownames(ATAC_WT_df)[irow],][[4 + sample]] <- 0
      }
    }
    else{
      ATAC_WT_reads_df[rownames(ATAC_WT_df)[irow],][[sample]] <- 0
      ATAC_WT_reads_df[rownames(ATAC_WT_df)[irow],][[4 + sample]] <- 0
    }
  }
}

#####
ATAC_MA9_reads_df <- ATAC_MA9_reads_df[1:4]

ATAC_WT_reads_df <- ATAC_WT_reads_df[1:4]

cat("There are",nrow(ATAC_MA9_reads_df),"genes plotted\nThe upper quantile is",get_quan(ATAC_MA9_reads_df),"\n")

Heatmap(ATAC_MA9_reads_df,cluster_rows = F,
        cluster_columns = F,
        row_names_side = "left",
        column_names_side = "top",
        na_col = 'ivory',
        #col = colorRamp2(c(-get_quan(ATAC_MA9_reads_df),0,get_quan(ATAC_MA9_reads_df)), c("blue","ivory", "red"))
        col = colorRamp2(c(-4,0,4), c("blue","ivory", "red"))
)



cat("There are",nrow(ATAC_WT_reads_df),"genes plotted\nThe upper quantile is",get_quan(ATAC_WT_reads_df),"\n")
Heatmap(ATAC_WT_reads_df,cluster_rows = F,
        cluster_columns = F,
        row_names_side = "left",
        column_names_side = "top",
        na_col = 'ivory'
        #col = colorRamp2(c(-get_quan(ATAC_WT_reads_df),0,get_quan(ATAC_WT_reads_df)), c("blue","ivory", "red"))
        col = colorRamp2(c(-4,0,4), c("blue","ivory", "red"))
)














