MA9_fc <- re_df[19:ncol(MA9)]
WT_fc <- WT[19:ncol(WT)]

#pathway_list <- read.xlsx("180521go_list_motif.xlsx",4)
pathway_list <- read.xlsx("all_diffgene_OSKM.xlsx",1)

MA9_df <- na.omit(MA9_fc[as.vector(pathway_list[,1]),])
WT_df <- na.omit(WT_fc[as.vector(pathway_list[,1]),])


k = 0

for (irow in 1:nrow(MA9_df)){
  rep1 <- MA9_df$MA9_A_Dox3vsCon0[irow]/MA9_df$MA9_A_Con3vsCon0[irow]
  rep2 <- MA9_df$MA9_B_Dox3vsCon0[irow]/MA9_df$MA9_B_Con3vsCon0[irow]
  if (rep1 != rep2){
    if (t.test(c(rep1,rep2))$p.value < 0.05){
      if (!data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
        k = k + 1
      }
    }
  }}
print(k)


t = 0
for (irow in 1:nrow(WT_df)){
  rep1 <- WT_df$WT_A_Dox3vsCon0[irow]/WT_df$WT_A_Con3vsCon0[irow]
  rep2 <- WT_df$WT_B_Dox3vsCon0[irow]/WT_df$WT_B_Con3vsCon0[irow]
  if (rep1 != rep2){
    if (t.test(c(rep1,rep2))$p.value < 0.05){
      if (!data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
        t = t + 1
      }
    }
  }}
print(t)



k = 0
for (irow in 1:nrow(MA9_df)){
  rep1 <- MA9_df$MA9_A_Dox6vsCon0[irow]/MA9_df$MA9_A_Con6vsCon0[irow]
  rep2 <- MA9_df$MA9_B_Dox6vsCon0[irow]/MA9_df$MA9_B_Con6vsCon0[irow]
  if (rep1 != rep2){
    if (t.test(c(rep1,rep2))$p.value < 0.05){
      if (!data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
        k = k + 1
      }
    }
  }}
print(k)


t = 0
for (irow in 1:nrow(WT_df)){
  rep1 <- WT_df$WT_A_Dox6vsCon0[irow]/WT_df$WT_A_Con6vsCon0[irow]
  rep2 <- WT_df$WT_B_Dox6vsCon0[irow]/WT_df$WT_B_Con6vsCon0[irow]
  if (rep1 != rep2){
    if (t.test(c(rep1,rep2))$p.value < 0.05){
      if (!data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
        t = t + 1
      }
    }
  }}
print(t)


k = 0
for (irow in 1:nrow(MA9_df)){
  rep1 <- MA9_df$MA9_A_Dox12vsCon0[irow]/MA9_df$MA9_A_Con12vsCon0[irow]
  rep2 <- MA9_df$MA9_B_Dox12vsCon0[irow]/MA9_df$MA9_B_Con12vsCon0[irow]
  if (rep1 != rep2){
    if (t.test(c(rep1,rep2))$p.value < 0.05){
      if (!data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
        k = k + 1
      }
    }
  }}
print(k)


t = 0
for (irow in 1:nrow(WT_df)){
  rep1 <- WT_df$WT_A_Dox12vsCon0[irow]/WT_df$WT_A_Con12vsCon0[irow]
  rep2 <- WT_df$WT_B_Dox12vsCon0[irow]/WT_df$WT_B_Con12vsCon0[irow]
  if (rep1 != rep2){
    if (t.test(c(rep1,rep2))$p.value < 0.05){
      if (!data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
        t = t + 1
      }
    }
  }}
print(t)

MA9_df <- read.csv("diffbind_nored_ATAC_MA9_fc.csv",header = T,row.names = 1)


t = 0
for (irow in 1:nrow(MA9_df)){
  rep1 <- log2(MA9_df$MA9_A_Dox3vsCon0[irow])/log2(MA9_df$MA9_A_Con3vsCon0[irow])
  rep2 <- log2(MA9_df$MA9_B_Dox3vsCon0[irow])/log2(MA9_df$MA9_B_Con3vsCon0[irow])
  if (rep1 != rep2){
    if (t.test(c(rep1,rep2))$p.value < 0.05){
      if (!data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
        MA9_df[irow,]$MA9_A_Dox3vsCon0 <- 666
        t = t + 1
      }
    }
    }
}
print(t)


write.csv(MA9_df,file = 'MA9_df.csv')

WT_df <- read.csv("diffbind_nored_ATAC_WT_fc.csv",header = T,row.names = 1)

t = 0

for (irow in 1:nrow(WT_df)){
  rep1 <- WT_df$WT_A_Dox24vsCon0[irow]/WT_df$WT_A_Con24vsCon0[irow]
  rep2 <- WT_df$WT_B_Dox24vsCon0[irow]/WT_df$WT_B_Con24vsCon0[irow]
  if (rep1 != rep2){
    if (t.test(c(rep1,rep2))$p.value < 0.05){
    if (!data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
      WT_df[irow,]$WT_A_Dox24vsCon0 <- 666
      t = t + 1
    }
  }
}}
print(t)


write.csv(WT_df,file = 'WT_df.csv')