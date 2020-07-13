MA9_fc <- read.csv("DESeq normalized count/MA9_normalized_fc.csv",header = T,row.names = 1)
WT_fc <- read.csv("DESeq normalized count/WT_normalized_fc.csv",header = T,row.names = 1)


k = 0

for (irow in 1:nrow(MA9_fc)){
  rep1 <- MA9_fc$MA9_A_Dox3vsCon0[irow]/MA9_fc$MA9_A_Con3vsCon0[irow]
  rep2 <- MA9_fc$MA9_B_Dox3vsCon0[irow]/MA9_fc$MA9_B_Con3vsCon0[irow]
  if (rep1 != rep2){
    if (t.test(c(rep1,rep2))$p.value < 0.05){
      if (!data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
        k = k + 1
      }
    }
  }}
print(k)


t = 0
for (irow in 1:nrow(WT_fc)){
  rep1 <- WT_fc$WT_A_Dox3vsCon0[irow]/WT_fc$WT_A_Con3vsCon0[irow]
  rep2 <- WT_fc$WT_B_Dox3vsCon0[irow]/WT_fc$WT_B_Con3vsCon0[irow]
  if (rep1 != rep2){
    if (t.test(c(rep1,rep2))$p.value < 0.05){
      if (!data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
        t = t + 1
      }
    }
  }}
print(t)



#20190130 测试ATAC-seq 用置信区间筛选过后的peak number

MA9_fc <- read.csv("diffbind_nored_ATAC_MA9_fc.csv",header = T,row.names = 1)


t = 0
for (irow in 1:nrow(MA9_fc)){
  rep1 <- log2(MA9_fc$MA9_A_Dox3vsCon0[irow])/log2(MA9_fc$MA9_A_Con3vsCon0[irow])
  rep2 <- log2(MA9_fc$MA9_B_Dox3vsCon0[irow])/log2(MA9_fc$MA9_B_Con3vsCon0[irow])
  if (rep1 != rep2){
    if (t.test(c(rep1,rep2))$p.value < 0.05){
      if (!data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
        MA9_fc[irow,]$MA9_A_Dox3vsCon0 <- 666
        t = t + 1
      }
    }
  }
}
print(t)


WT_fc <- read.csv("diffbind_nored_ATAC_WT_fc.csv",header = T,row.names = 1)

t = 0

for (irow in 1:nrow(WT_fc)){
  rep1 <- WT_fc$WT_A_Dox24vsCon0[irow]/WT_fc$WT_A_Con24vsCon0[irow]
  rep2 <- WT_fc$WT_B_Dox24vsCon0[irow]/WT_fc$WT_B_Con24vsCon0[irow]
  if (rep1 != rep2){
    if (t.test(c(rep1,rep2))$p.value < 0.05){
      if (!data.table::between(1, t.test(c(rep1,rep2))$conf.int[1],t.test(c(rep1,rep2))$conf.int[2])){
        WT_fc[irow,]$WT_A_Dox24vsCon0 <- 666
        t = t + 1
      }
    }
  }}
print(t)

