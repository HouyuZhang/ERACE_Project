source('FimoAndMakeMotifList.R');


#~ ##### hg19 ####################################
#~ FIMOandMakeMotifList(memefile = 'hg_transfac_jaspar_uniprobe_meme.txt',
#~ 	memeindexfile = 'hg_transfac_jaspar_uniprobe_meme_index.csv',
#~ 	genomefile = '/home/baeks/Data/hg19/hg19.fa',
#~ 	outputfilename= 'hg19motifs715',
#~ 	fimooutpath = 'hg19_fimo',
#~ 	minsize = 10000,
#~ 	maxstoredscores = 1500000,
#~ 	thresh = 0.0001,
#~ 	MCCORES = NA)


#~ ##### mm10 #####################################

#~ ref = 'mm10';
#~ motiffile = 'mm10motifs634';
#~ minsize = 1000;
#~ fimooutpath = '/mnt/Data1/baeks/motifs/mm10_fimo';
#~ memefile = '/mnt/Data1/baeks/project/footprinting_casellas/dhs/fimo/All_JasparHomerUniLeoLitNFKB_pretty.meme';
#~ memeindexfile = '/mnt/Data1/baeks/project/footprinting_casellas/dhs/fimo/All_JasparHomerUniLeoLit_index.txt';

#~ FIMOandMakeMotifList(memefile = 'transfac_jaspar_uniprobe_meme.txt',
#~ 	memeindexfile = 'transfac_jaspar_uniprobe_meme_index.csv',
#~ 	genomefile = '/home/baeks/Data/hg19/hg19.fa',
#~ 	outputfilename= 'hg19motifs715',
#~ 	fimooutpath = 'hg19_fimo',
#~ 	minsize = 10000,
#~ 	maxstoredscores = 1500000,
#~ 	thresh = 0.0001,
#~ 	MCCORES = NA)


##### hg19 CTCF ####################################
 FIMOandMakeMotifList(memefile = 'hg_ctcf_meme.txt',
 	memeindexfile = 'hg_ctcf_meme_index.csv',
 	genomefile = '/home/baeks/Data/hg19/hg19.fa',
 	outputfilename= 'hg19_CTCF',
 	fimooutpath = 'hg19_CTCF_fimo',
 	minsize = 10000,
 	maxstoredscores = 1500000,
 	thresh = 0.0001,
 	MCCORES = NA)
