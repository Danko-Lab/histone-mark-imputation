source("../scripts/hist.svm.com.R")

tb <- read.table("./GBM20.file.tab", stringsAsFactors=F);
files.bigwig.plus <- tb[,1]
gen.bed <- read.table("gen3-pred.bed", stringsAsFactors=F)

for(i in 1:NROW(files.bigwig.plus))
{
	file.H3k27ac.bed  <- paste( "./GBM20.H3k27ac.",  basename(files.bigwig.plus[i]), ".pred.bed.gz", sep="");
	file.H3k27me3.bed <- paste( "./GBM20.H3k27me3.", basename(files.bigwig.plus[i]), ".pred.bed.gz", sep="");
	file.H3k36me3.bed <- paste( "./GBM20.H3k36me3.", basename(files.bigwig.plus[i]), ".pred.bed.gz", sep="");
	file.H3k4me1.bed  <- paste( "./GBM20.H3k4me1.",  basename(files.bigwig.plus[i]), ".pred.bed.gz", sep="");
	file.H3k4me3.bed  <- paste( "./GBM20.H3k4me3.",  basename(files.bigwig.plus[i]), ".pred.bed.gz", sep="");
	file.H3k9me3.bed  <- paste( "./GBM20.H3k9me3.",  basename(files.bigwig.plus[i]), ".pred.bed.gz", sep="");

    for(file.hist in c(file.H3k27ac.bed, file.H3k27me3.bed, file.H3k36me3.bed, file.H3k4me1.bed, file.H3k4me3.bed, file.H3k9me3.bed))
    {
	    tb <- read.table(file.hist, stringsAsFactors=F);
	    for(k in 1:NROW(gen.bed))
	    {
	    	idx.rem <- which( tb[,1]==as.character(gen.bed[k,1]) & tb[,2]>=gen.bed[k,2] & tb[,2]<=gen.bed[k,3])
	    	if(NROW(idx.rem)>0) tb <- tb[-idx.rem,];
	    }	
	    
	    new.bed <- paste("GBM20G3", paste(strsplit(file.hist, "[.]")[[1]][-c(1:2)], collapse="."), sep=".");
	    tb2 <- read.table(new.bed, stringsAsFactors=F);
	    
	    file.tmp <- tempfile(fileext=".bed.gz")
	    write.bed(as.data.frame(rbind(tb,tb2)), file=file.hist, compress=TRUE); 
    }
}
