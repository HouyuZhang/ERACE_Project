
my_gen_bagplot<-function(dat, dataname1=dataname1, dataname2=dataname2, factor=3, pdf=T,  noBag=F, noFence=F, noMedianPoint=F) {
  x=dat$lrtotal 
  y=dat$fdepth1-dat$fdepth2
  
  trimbp <- function(sss) {
    doit <- function(sss) {
      ttt= gregexpr('_[0-9]+bp', sss)[[1]][1]
      ifelse(ttt<1, sss, substring(sss, 1,ttt-1))
    }
    sapply(sss, doit)
  }
  
  name = trimbp(as.character(dat$motif))
  
  datan = data.frame(x=x, y=y, name=as.character(name),stringsAsFactors=F)
  ss <- aplpack::compute.bagplot(datan$x, datan$y, factor=factor, approx.limit = nrow(datan))
  
  outlier = ss$pxy.outlier
  outer = ss$pxy.outer
  bag = ss$pxy.bag
  hdepths = ss$hdepths
  dat = datan
  
  findOutliersIndex<-function(data, outlier) {
    cols <- colnames(outlier)
    dat2 <- data[,cols]
    nr = nrow(outlier)
    rg = 1:nrow(dat2)
    matchingIdx = unique(sort(unlist(sapply(1:nr, function(ii) rg[(outlier[ii,2] == dat2[,2]) & (outlier[ii,1]==dat2[,1])]))))
    if (length(matchingIdx) != nr) {
      stop('cannot locate outliers in the data.')
    }
    matchingIdx
  }
  
  nametoshow = rep("", nrow(datan))
  nametoshowInGray = rep("", nrow(datan))
  if (!is.null(outlier)) {
    outlieridx = findOutliersIndex(datan, outlier)
    nametoshow[outlieridx] = as.character(datan$name[outlieridx])
  }
  
  outeridx = findOutliersIndex(datan, outer)
  inneridx = findOutliersIndex(datan, bag)
  nametoshowInGray[outeridx] = as.character(datan$name[outeridx])
  
  category = vector("character", length=nrow(datan))
  category[inneridx] = "bag"
  category[outeridx] = "fence"
  category[outlieridx] = "outlier"
  
  mx = mean(datan$x)
  my = mean(datan$y)
  
  udist = sqrt((datan$x-mx)^2 + (datan$y-my)^2)
  bagplottable= data.frame(name=datan$name, deltahyp=datan$x, deltafd = datan$y, category=category, hdepth= hdepths, udist=udist)
  
  outputfilename=sprintf("bagplot_cutcount_diff_total_footprinting_depth_%s-%s",dataname2,dataname1)
  write.csv(bagplottable, file=sprintf('%s_bagplot_output.csv', outputfilename))
  
  labx=sprintf("Normalized Cutcount Diff. (%s-%s)",dataname2,dataname1)
  laby=sprintf("-Footprinting Depth Diff. (%s-%s)",dataname2,dataname1)
  
  do_draw<-function() {
    ss <- aplpack::bagplot(datan$x, datan$y, factor=3,
                           show.looppoints = F,show.bagpoints=F, 
                           show.baghull = !noBag,  show.loophull= !noFence,  
                           show.whiskers=F,xlab=labx,ylab=laby, 
                           col.baghull="turquoise3",col.loophull="turquoise",
                           cex=1.5,pch=20, cex.lab=1.5,cex.axis=1.5)
    
    points(datan$x[outeridx], datan$y[outeridx], pch=16, col= "black",cex=0.5)      # outer
    points(datan$x[inneridx], datan$y[inneridx], pch=16, col= "black",cex=0.5)   # bag
    text(datan$x, datan$y, labels=nametoshow,adj=1.1, cex=1,font = 2)
    text(datan$x, datan$y, labels=nametoshowInGray,col="#000000A0", adj=1.1, cex=0.5)
    plotrange <- par("usr")
    xmin = plotrange[1]
    xmax = plotrange[2]
    ymin = plotrange[3]
    ymax = plotrange[4]
    lines(c(xmin, xmax),c(0,0),lty = 2,col="grey",lwd=2)
    lines(c(0,0),c(ymin, ymax),lty = 2,col="grey",lwd=2)
  }
  do_draw()
  
  outputfilename=sprintf("bagplot_cutcount_diff_total_footprinting_depth_%s-%s",dataname2,dataname1)
  tiff(filename= sprintf("%s.tiff",outputfilename),compression="zip", width=2000,height=1500, units="px", pointsize=20)
  do_draw()
  dev.off()
  
  if (pdf) {
    pdf(file= sprintf("%s.pdf",outputfilename), width=16,height=12, pointsize=10,useDingbats=FALSE)
    do_draw()
    dev.off()
  }
  
}

my_gen_bagplot_chisq<-function(dat, dataname1=dataname1, dataname2=dataname2, qvaluethreshold=0.05,factor=3, pdf=T,  
                            noBag=F, noFence=F, noMedianPoint=F,pdfwidth=16, pdfheight=12, textsizefactor=1.0, printmode=F) {
  x=dat$lrtotal 
  y=dat$fdepth1-dat$fdepth2 
  
  trimbp <- function(sss) {
    doit <- function(sss) {
      ttt= gregexpr('_[0-9]+bp', sss)[[1]][1]
      ifelse(ttt<1, sss, substring(sss, 1,ttt-1))
    }
    sapply(sss, doit)
  }
  
  mat <- data.frame(dhp=dat$lrtotal, dfd=dat$fdepth1-dat$fdepth2)
  name = trimbp(as.character(dat$motif))
  
  datan = data.frame(x=x, y=y, name=as.character(name),stringsAsFactors=F)
  ss <- aplpack::compute.bagplot(datan$x, datan$y, factor=factor, approx.limit = nrow(datan))
  
  outlier = ss$pxy.outlier
  outer = ss$pxy.outer
  bag = ss$pxy.bag
  hdepths = ss$hdepths
  dat = datan
  
  findOutliersIndex<-function(data, outlier) {
    cols <- colnames(outlier)
    dat2 <- data[,cols]
    nr = nrow(outlier)
    rg = 1:nrow(dat2)
    matchingIdx = unique(sort(unlist(sapply(1:nr, function(ii) rg[(outlier[ii,2] == dat2[,2]) & (outlier[ii,1]==dat2[,1])]))))
    if (length(matchingIdx) != nr) {
      stop('cannot locate outliers in the data.')
    }
    matchingIdx
  }
  
  nametoshow = rep("", nrow(datan))
  nametoshowInGray = rep("", nrow(datan))
  if (!is.null(outlier)) {
    outlieridx = findOutliersIndex(datan, outlier)
    nametoshow[outlieridx] = as.character(datan$name[outlieridx])
  }
  
  outeridx = findOutliersIndex(datan, outer)
  inneridx = findOutliersIndex(datan, bag)
  nametoshowInGray[outeridx] = as.character(datan$name[outeridx])
  
  category = vector("character", length=nrow(datan))
  category[inneridx] = "bag"
  category[outeridx] = "fence"
  category[outlieridx ] = "outlier"
  
  mx = mean(datan$x)
  my = mean(datan$y)
  
  Sx <- cov(mat)
  #D2 <- mahalanobis(mat, colMeans(mat), Sx)
  D2 <- mahalanobis(mat,apply(mat, 2, median), Sx)
  pvalue = 1-pchisq(D2, ncol(mat))                
  qvalue = p.adjust(pvalue, method="BH")  # FDR adjusted p-value (BH)
  #browser()
  
  udist = sqrt((datan$x-mx)^2 + (datan$y-my)^2)
  bagplottable= data.frame(name=datan$name, deltahyp=datan$x, deltafd = datan$y, category=category, hdepth= hdepths, udist=udist,pvalue=pvalue,qvalue=qvalue)
  
  outputfilename=sprintf("bagplot_cutcount_diff_total_footprinting_depth_%s-%s_qvalue",dataname2,dataname1)
  write.csv(bagplottable, file=sprintf('%s_bagplot_output.csv', outputfilename))
  
  labx=sprintf("Normalized Cutcount Diff. (%s-%s)",dataname2,dataname1)
  laby=sprintf("-Footprinting Depth Diff. (%s-%s)",dataname2,dataname1)
  
  do_draw<-function() {
    ss <- aplpack::bagplot(datan$x, datan$y, factor=3,
                           show.looppoints = F,show.bagpoints=F, 
                           show.baghull = !noBag,  show.loophull= !noFence,  
                           show.whiskers=F,xlab=labx,ylab=laby, 
                           #col.baghull="turquoise3",col.loophull="turquoise",
                           cex=1,pch=1, cex.lab=1.5,cex.axis=1.5)
    
    points(datan$x[outeridx], datan$y[outeridx], pch=16, col= "black",cex=0.5)   # outer
    points(datan$x[inneridx], datan$y[inneridx], pch=16, col= "black",cex=0.5)   # bag
    
    #PVALUETHRESHOLD = 0.05
    points(datan$x[qvalue < qvaluethreshold], datan$y[qvalue < qvaluethreshold], cex=1.5, pch=19, col="red")   #  p-value < 0.05q
    
    text(datan$x, datan$y, labels=nametoshow,adj=1.1, cex=1,font = 2)
    text(datan$x, datan$y, labels=nametoshowInGray,col="#000000A0", adj=1.1, cex=0.5)
    plotrange <- par("usr")
    xmin = plotrange[1]
    xmax = plotrange[2]
    ymin = plotrange[3]
    ymax = plotrange[4]
    lines(c(xmin, xmax),c(0,0),lty = 2,col="grey",lwd=2)
    lines(c(0,0),c(ymin, ymax),lty = 2,col="grey",lwd=2)

  }
  #do_draw()
  
  #outputfilename=sprintf("bagplot_cutcount_diff_total_footprinting_depth_%s-%s_qvalue%g",dataname2,dataname1,qvaluethreshold)
  #tiff(filename= sprintf("%s.tiff",outputfilename),compression="zip", width=2000,height=1500, units="px", pointsize=20)
  #do_draw()
  #dev.off()
  
  if (pdf) {
    pdf(file= sprintf("%s.pdf",outputfilename), width=16,height=12, pointsize=10,useDingbats=FALSE)
    do_draw()
    dev.off()
  }
  
}
