#setwd('/mnt/Data1/baeks/motifs/progs');

library(parallel);
library(gdata)   
#source('/mnt/Data1/MA/sjbaek/project/footprinting_NatureCorres/Aggregation_stam_bias_correction/calcExpectedCuts.R');


###################################################  Util Functions ###############
dna=c("A","C","G","T");
str_trim <- function (x) gsub("^\\s+|\\s+$", "", x);

readMeme<- function(memefile) {
	meme = readLines(memefile);
	motifstart= grep('MOTIF', meme);
	motifnames=str_trim(sub('MOTIF ','',meme[motifstart]));
	
	motifnames = sapply(motifnames, function(x) strsplit(x," ")[[1]][1]);
	
	#browser();

	inforow= grep("letter-probability matrix: alength=", meme);
	rr=meme[inforow];
	nmotif= length(inforow);
	rr=gsub('letter-probability matrix: alength= ','',rr);
	rr=gsub('w=','',rr);
	rr=gsub('nsites=','',rr);
	rr=gsub('E=','',rr);
	tt=do.call("rbind", lapply(strsplit( rr,"  "), as.numeric));
	alength= tt[,1];
	wmeme = tt[,2];
	nsites=tt[,3];
	Ememe = tt[,4];

	pwmstart= inforow+1;
	pwmend = inforow+wmeme;


	string2double <- function(s) {
		ss = gsub('  ',' ',str_trim(s));
		ss = gsub('\t',' ',ss);
		ss = gsub('  ',' ',ss);
		t<-read.csv(text=ss, sep=" ", header=F);
		names(t) <- c('A','C','G','T');
		t;
	}
	pwm=lapply(1:nmotif,	function(x) { string2double(meme[pwmstart[x]:pwmend[x]]);});

	list(name=motifnames, alength=alength,w=wmeme, nsites=nsites, E=Ememe, pwm=pwm);
}

revComplement <- function(M) {
	nr = nrow(M);
	RC <- cbind(A=M[nr:1,4],C=M[nr:1,3],G=M[nr:1,2], T=M[nr:1,1]);
	RC;
}

Seqs2PWM <- function(TT) {
	seqarray= do.call(rbind, lapply(TT, function(s) { sapply(1:nchar(s), function(x) { substring(s,x,x); }); } ));
	freqmatrix=t(sapply(1:ncol(seqarray), function(y) { sapply(dna, function(x) { sum(seqarray[,y]==x);}); })/nrow(seqarray));
	freqmatrix;
}

distanceFMatrix<-function(freq1, freq2)  {
# returns 10000 if dimension of freq1 and freq2 differ. 
	nr1 = nrow(freq1);
	nr2 = nrow(freq2);
	if (nr1 != nr2) {
		ss= 10000;
	} else {
		diff= abs(freq1 - freq2);
		ss =  sum(diff);
	}
	ss;

}


palindromicMetric<-function(freqmatrix,center) {
	aa=t(freqmatrix);
	bb=aa[4:1, ncol(aa):1];
	#browser();
	right = ncol(aa);
	#center = 6;
	shift = right - 2*center + 1;
	ap = array(0, dim=c(4, right+abs(shift)));
	bp = ap;

	if (shift >0 ) {
		ap[,-(1:shift)] = aa;
		bp[,1:right] = bb;
	} else if (shift <0){
		ap[,1:right] = aa;
		bp[,(1-shift):(right-shift)] = bb;
	} else {
		ap = aa; bp =bb;
	}
	sum(sapply(1:ncol(ap), function(ii) { abs(ap[1,ii]-bp[4,ii])+abs(ap[2,ii]-bp[3,ii])+abs(bp[1,ii]-ap[4,ii])+abs(bp[2,ii]-ap[3,ii]); }));
}

filenameWithNoExt<-function(filename) {
	str=basename(filename);
	dir = dirname(filename);
	idxperiod= gregexpr("\\.",str)[[1]][1];
	if (idxperiod>0) {
		filename = substring(str, 1, idxperiod-1);
	};
	filename;
}




makeDirectory <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir);
  }
}

###################################################  Util Functions ###############


FIMOandMakeMotifList <- function(memefile = '',
	memeindexfile = '',
	genomefile = '',
	outputfilename= '',
	fimooutpath = '',
	minsize = 10000,
	maxstoredscores = 1500000,
	thresh = 0.0001,
	MCCORES = NA)
 {



	if (is.na(MCCORES)) {
		MCCORES =  detectCores(all.tests = FALSE, logical = F);
	}
	
	fimoprog = 'fimo';
	gg=system(sprintf('which %s',fimoprog),intern=T);

	if (length(gg) < 1) stop(paste(fimoprog,'is not found.'));

	if (!file.exists(memefile)) stop(paste('The meme output-format',memefile,'does not exist.'));
	if (!file.exists(genomefile)) stop(paste('The reference genome',genomefile,'does not exist.'));

	meme <- readMeme(memefile);

	memename = str_trim(meme$name);
	if (memeindexfile !='' || !file.exists(memeindexfile)) {
		memeindex = read.csv(memeindexfile,stringsAsFactors =F);
		motifname = str_trim(as.character(memeindex$Name));
	} else {
		motifname = memename;
	}

	shortnames = motifname;
	lengths= meme$w;
	motiffilenames= vector("character",length=length(shortnames));       
	motifcenter=vector("numeric",length=length(shortnames));       
	motifseqs=vector("character",length=length(shortnames));       
	motiflength=vector("numeric",length=length(shortnames));       
	motifs=vector("character",length=length(shortnames));       
	memeNo=vector("numeric",length=length(shortnames));       
	memeEntry=vector("character",length=length(shortnames));       

	makeDirectory(fimooutpath);

	runFIMO <- function(ii) {
		motiffilename = file.path(fimooutpath,make.names(sprintf('%s_%g_fimo_hit.txt', make.names(meme$name[ii]),thresh,filenameWithNoExt(outputfilename))));
		if (!file.exists(motiffilename)) {
			fimocmd = sprintf('fimo --text --motif %s --max-stored-scores %f --thresh %f   --verbosity 1 %s %s >%s', memename[ii], maxstoredscores,
			thresh,memefile, genomefile, motiffilename);
			cat(sprintf("running:%s\n", fimocmd));
			system(fimocmd);
		}
	}

	r <- mclapply(1:length(shortnames), runFIMO,  mc.cores = MCCORES);
	#r <- lapply(1:length(shortnames), runFIMO);

	summerizeMotif<-function(ii) {
		groupname = shortnames[ii];
		fq = meme$pwm[[ii]];

		ll=sapply(1:nrow(fq),function(x) { cc=dna[fq[x,]>0.5]; ifelse(length(cc)!=1, '_',cc);});
		seqpattern=paste(unlist(ll),sep="",collapse="");

		motiffilename = file.path(make.names(sprintf('%s_%g_fimo_hit.txt', make.names(meme$name[ii]),thresh,filenameWithNoExt(outputfilename))));

		mcenter=which.min(sapply(1:nrow(fq), function(x) { palindromicMetric(fq, x); } ));
		motifwidth = lengths[ii];
		list(filename= motiffilename,    
		center= mcenter,
		logo = seqpattern,
		motiflength = motifwidth,
		motif = sprintf("%s_%dbp",groupname, motifwidth),
		memeNo = ii,
		memeEntry = motifname[ii]);
	}

	ll <- mclapply(1:length(shortnames), summerizeMotif, mc.cores=MCCORES);

	motiftable = do.call(rbind, lapply(ll, function(x) data.frame(x, stringsAsFactors=FALSE)));
	uniqueIDs = cumsum(rle(as.character(motiftable$motif))$lengths);
	motiftable = motiftable[uniqueIDs,];     # remove non-unique tags   #724

	filepaths= file.path(normalizePath(fimooutpath),motiftable$filename);
	isFileBigEnough <- function(filename, threshold= minsize) {
		test = file.info(filename)$size > threshold;
		test;
	}

	bigenough = sapply(filepaths, isFileBigEnough)
	motiftable = subset(motiftable,bigenough);  
	motiflistfile = make.names(sprintf('motiflist_%s.txt', filenameWithNoExt(outputfilename)));
	cat(sprintf('\n\n\n%s is now created.\n', motiflistfile));
	write.table(motiftable, file=motiflistfile);
}
