library(DiffBind)
MA9_count <- dba.load('MA9')
WT_count <- dba.load('WT')

for (i in 1:18){
  cat(MA9_count$samples[[1]][i],'\t',sum(MA9_count$peaks[[i]]$Reads)/sum(MA9_count$peaks[[1]]$Reads),'\n')
  cat(WT_count$samples[[1]][i],'\t',sum(WT_count$peaks[[i]]$Reads)/sum(WT_count$peaks[[1]]$Reads),'\n')
}
for (i in 1:18){density(length(x), bw = bandwidth)
  MA9_count$peaks[[i]]$peak_name <- paste('peak_',rownames(MA9_count$peaks[[i]]),sep = '')
  WT_count$peaks[[i]]$peak_name <- paste('peak_',rownames(WT_count$peaks[[i]]),sep = '')
}




WT_conrast <- dba.contrast(WT_count,categories = DBA_CONDITION,minMembers = 2)
WT_conrast <- dba.analyze(WT_conrast)
WT_conrast

WT_Con0_Con3.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=1,th=1)
write.table(WT_Con0_Con3.DB,file = "WTCon0_VS_Con3_500.txt",sep = '\t')

WT_Con3_Con6.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=9,th=1)
write.table(WT_Con3_Con6.DB,file = "WTCon3_VS_Con6_500.txt",sep = '\t')

WT_Con6_Con12.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=16,th=1)
write.table(WT_Con6_Con12.DB,file = "WTCon6_VS_Con12_500.txt",sep = '\t')

WT_Con12_Con24.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=22,th=1)
write.table(WT_Con12_Con24.DB,file = "WTCon12_VS_Con24_500.txt",sep = '\t')

WT_Con0_Dox3.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=5,th=1)
write.table(WT_Con0_Dox3.DB,file = "WTCon0_VS_Dox3_500.txt",sep = '\t')

WT_Dox3_Dox6.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=31,th=1)
write.table(WT_Dox3_Dox6.DB,file = "WTDox3_VS_Dox6_500.txt",sep = '\t')

WT_Dox6_Dox12.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=34,th=1)
write.table(WT_Dox6_Dox12.DB,file = "WTDox6_VS_Dox12_500.txt",sep = '\t')

WT_Dox12_Dox24.DB <- dba.report(WT_conrast,bCounts=TRUE,contrast=36,th=1)
write.table(WT_Dox12_Dox24.DB,file = "WTDox12_VS_Dox24_500.txt",sep = '\t')












MA9_conrast <- dba.contrast(MA9_count,categories = DBA_CONDITION,minMembers = 2)
MA9_conrast <- dba.analyze(MA9_conrast)
MA9_conrast


MA9_Con0_Con3.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=1,th=1)
write.table(MA9_Con0_Con3.DB,file = "MA9Con0_VS_Con3_500.txt",sep = '\t')

MA9_Con3_Con6.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=9,th=1)
write.table(MA9_Con3_Con6.DB,file = "MA9Con3_VS_Con6_500.txt",sep = '\t')

MA9_Con6_Con12.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=16,th=1)
write.table(MA9_Con6_Con12.DB,file = "MA9Con6VS_Con12_500.txt",sep = '\t')

MA9_Con12_Con24.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=22,th=1)
write.table(MA9_Con12_Con24.DB,file = "MA9Con12_VS_Con24_500.txt",sep = '\t')

MA9_Con0_Dox3.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=5,th=1)
write.table(MA9_Con0_Dox3.DB,file = "MA9Con0_VS_Dox3_500.txt",sep = '\t')

MA9_Dox3_Dox6.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=31,th=1)
write.table(MA9_Dox3_Dox6.DB,file = "MA9Dox3_VS_Dox6_500.txt",sep = '\t')

MA9_Dox6_Dox12.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=34,th=1)
write.table(MA9_Dox6_Dox12.DB,file = "MA9Dox6_VS_Dox12_500.txt",sep = '\t')

MA9_Dox12_Dox24.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=36,th=1)
write.table(MA9_Dox12_Dox24.DB,file = "MA9Dox12_VS_Dox24_500.txt",sep = '\t')






MA9_Con0_Con3.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=1,th=1)
write.table(MA9_Con0_Con3.DB,file = "MA9Con0_VS_Con3_500.txt",sep = '\t')

MA9_Con0_Con6.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=9,th=1)
write.table(MA9_Con0_Con6.DB,file = "MA9Con0_VS_Con6_500.txt",sep = '\t')

MA9_Con0_Con12.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=16,th=1)
write.table(MA9_Con0_Con12.DB,file = "MA9Con0_VS_Con12_500.txt",sep = '\t')

MA9_Con0_Con24.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=22,th=1)
write.table(MA9_Con0_Con24.DB,file = "MA9Con0_VS_Con24_500.txt",sep = '\t')

MA9_Con0_Dox3.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=5,th=1)
write.table(MA9_Con0_Dox3.DB,file = "MA9Con0_VS_Dox3_500.txt",sep = '\t')

MA9_Con0_Dox6.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=31,th=1)
write.table(MA9_Con0_Dox6.DB,file = "MA9Con0_VS_Dox6_500.txt",sep = '\t')

MA9_Con0_Dox12.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=34,th=1)
write.table(MA9_Con0_Dox12.DB,file = "MA9Con0_VS_Dox12_500.txt",sep = '\t')

MA9_Con0_Dox24.DB <- dba.report(MA9_conrast,bCounts=TRUE,contrast=36,th=1)
write.table(MA9_Con0_Dox24.DB,file = "MA9DCon0_VS_Dox24_500.txt",sep = '\t')
