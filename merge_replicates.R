suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-A", "--rep_A"), type="character", default=NULL, help="rep_A directory", metavar="character"),
  make_option(c("-B", "--rep_B"), type="character", default=NULL, help="rep_B directory", metavar="character")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

A <- list.files(opt$rep_A,"footprint_depth_table.csv")
B <- list.files(opt$rep_B,"footprint_depth_table.csv")
dat_A = read.table(paste(opt$rep_A,"/",A,sep = ""))
dat_B = read.table(paste(opt$rep_B,"/",B,sep = ""))
dat <- dat_A

for (i in seq(ncol(dat_A)-8)){
  dat[i+8] <- (dat_A[i+8]+dat_B[i+8])/2
}
dataname <- strsplit(opt$rep_A,"vs")
dataname1 <- sub("_[AB]_","_merge_",dataname[[1]][1])
dataname2 <- sub("_[AB]_","_merge_",dataname[[1]][2])

source("gen_bagplot.R")
my_gen_bagplot_chisq(dat, dataname1=dataname1, dataname2=dataname2, factor=3)