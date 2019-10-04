(Rgtsvm)
library(bigWig)
library(dREG)
library(data.table);
library(snowfall);

#source("../script/hist.svm.main.R");
#source("../script/hist.svm.com.R");
source("../script/hist.param.R");

svm_model_scatter <- function( file.rdata, file.png, title="", ncores=15 )
{
	load(file.rdata);

	if(is.null(model$y_pred))
	{
		mloaded <- Rgtsvm::predict.load(model$svm);

		y_test <- y_pred <- c();
		for(i in 1:length(model$test))
		{
			x1 <- extract_feature_matrix( model$test[[i]]$pos, paste(model$src$path.proseq, model$src$file.bws.plus[i], sep="/"),
							paste(model$src$path.proseq, model$src$file.bws.minus[i], sep="/"), model$gdm, linear_scale=F, ncores=ncores )

			y_test <- c(y_test, model$test[[i]]$y);
			y_pred <- c(y_pred, Rgtsvm::predict.run( mloaded, x1$mat) );
		}

		model$y_pred <- y_pred;
		model$y_test <- y_test;
		model$pred.cor <- cor( y_test, y_pred );
		save(model, file=file.rdata);

		Rgtsvm::predict.unload( mloaded );
	}

	cat("PRED COR=", model$pred.cor, "\n");

	png(file.png, width=900, height=900);
	densScatterplot( model$y_test, model$y_pred, main=paste(title, "(r=",model$pred.cor, ")")  );
	dev.off();
	return;
}

if(0)
{
	svm_model_scatter( "H3k27ac.train.rdata", "H3k27ac.test.png",  title="H3k27ac");
	svm_model_scatter("H3k27me3.train.rdata", "H3k27me3.test.png", title="H3k27me3");
	svm_model_scatter("H3k36me3.train.rdata", "H3k36me3.test.png", title="H3k36me3");
	svm_model_scatter( "H3k4me1.train.rdata", "H3k4me1.test.png",  title="H3k4me1");
	svm_model_scatter( "H3k4me2.train.rdata", "H3k4me2.test.png",  title="H3k4me2");
	svm_model_scatter( "H3k4me3.train.rdata", "H3k4me3.test.png",  title="H3k4me3");
	svm_model_scatter( "H3K9me3.train.rdata", "H3k9me3.test.png",  title="H3k9me3");
	svm_model_scatter(  "H3k9ac.train.rdata", "H3k9ac.test.png",   title="H3k9ac");
	svm_model_scatter("H4k20me1.train.rdata", "H4k20me1.test.png", title="H4k20me1");

}


draw_scatter <- function( file.rdata, file.png, title="", ncores=15 )
{
	load(file.rdata);
	cat("PRED COR=", model$pred.cor, "\n");

	png(file.png, width=900, height=900);
	lim <- min( max(model$y_test), max(model$y_pred) )*2;
	densScatterplot( model$y_test, model$y_pred, main=paste(title, "(cor=",model$pred.cor, ")"), xlim=c(0, lim), ylim=c(0, lim), xlab="Histone data", ylab="Prediction" );
	segments(0,0, lim, lim);
	dev.off();
	return;
}


if(0)
{
	draw_scatter( "H3k27ac.train.rdata", "H3k27ac.test.opt.png",  title="H3k27ac");
	draw_scatter("H3k27me3.train.rdata", "H3k27me3.test.opt.png", title="H3k27me3");
	draw_scatter("H3k36me3.train.rdata", "H3k36me3.test.opt.png", title="H3k36me3");
	draw_scatter( "H3k4me1.train.rdata", "H3k4me1.test.opt.png",  title="H3k4me1");
	draw_scatter( "H3k4me2.train.rdata", "H3k4me2.test.opt.png",  title="H3k4me2");
	draw_scatter( "H3k4me3.train.rdata", "H3k4me3.test.opt.png",  title="H3k4me3");
	draw_scatter( "H3k9me3.train.rdata", "H3k9me3.test.opt.png",  title="H3k9me3");
	draw_scatter(  "H3k9ac.train.rdata", "H3k9ac.test.opt.png",   title="H3k9ac");
	draw_scatter("H4k20me1.train.rdata", "H4k20me1.test.opt.png", title="H4k20me1");
}


genomewide_cor <- function(file.plus, file.minus, file.histone, file.imputed, str.chr=NULL)
{
   chrom.info.table <- get.chromosome.info( file.plus, file.minus );

   dnase_bed <- as.data.frame(rbindlist( lapply(1:NROW(chrom.info.table), function(i){ data.frame(chrom.info.table[i,1], seq(1, (chrom.info.table[i,2]-1), 10), seq(1, (chrom.info.table[i,2]-1), 10)+1) }) ));
   colnames(dnase_bed)<-NULL

   if(!is.null(str.chr))
       dnase_bed <- dnase_bed[ as.character(dnase_bed[,1]) %in% c( str.chr ),];

   pred <- get_histone_read(dnase_bed, file.imputed);
   label <- get_histone_read(dnase_bed, file.histone);
   r.cor1 <- cor(pred, label, method ="pearson");
   r.cor2 <- cor(pred, label, method ="spearman");
   cat("COR=", r.cor1, r.cor2, "\n" );

   return(r.cor2);
}

compare_correlation<-function(file.bed, bigwig_histone)
{
   r.cor <- NA;
   if( file.exists(file.bed) && !is.null(bigwig_histone) )
   {
       ret <- read.table(file.bed);
       label <-  get_histone_read( ret[,c(1:3)], bigwig_histone);
       r.cor <- cor(ret[,4], label, );
       cat("COR=", cor(ret[,4], label), "\n" )
    }

	return(r.cor);
}



compare_predict0<-function(file.bed, bigwig_histone)
{
   if( !file.exists(file.bed))
   	   stop("file.bed is not existing.");
   if( is.null(bigwig_histone) )
   	   stop("file.bigwig is not existing.");

   tb.pred <- read.table(file.bed);
   y.pred <- y.bigwig <- c();
   y.chr <- y.cor1 <- y.cor2 <- y.cor3 <- y.mad <- c();

   for(chr in unique(tb.pred[,1]))
   {
      label <-  get_histone_read( tb.pred[ tb.pred[,1]==chr, c(1:3)], bigwig_histone);
      y.pred <- c( y.pred, tb.pred[ tb.pred[,1]==chr ,4] );
      y.bigwig <- c(y.bigwig, label);

	  r.cor1 <- cor(tb.pred[ tb.pred[,1]==chr ,4], label, method = "pearson" );
      r.cor3 <- cor(tb.pred[ tb.pred[,1]==chr ,4], label, method = "spearman" );
      r.mad  <- mad( tb.pred[ tb.pred[,1]==chr ,4] - label);
      cat(chr, "COR=", r.cor1, r.cor3, r.mad, "\n" );

      y.chr <- c(y.chr, chr);
      y.cor1 <- c(y.cor1, r.cor1);
      y.cor3 <- c(y.cor3, r.cor3);
      y.mad  <- c(y.mad, r.mad);
   }

   r.cor1 <- cor(tb.pred[ tb.pred[,1]==chr ,4], label, method = "pearson" );
   r.cor3 <- cor(tb.pred[ tb.pred[,1]==chr ,4], label, method = "spearman" );
   r.mad  <- mad( tb.pred[ tb.pred[,1]==chr ,4] - label);
   y.chr <- c(y.chr, "all");
   cat("Final COR=", r.cor1, r.cor3, r.mad, "\n" );

   y.cor1 <- c(y.cor1, r.cor1);
   y.cor3 <- c(y.cor3, r.cor3);
   y.mad  <- c(y.mad, r.mad);

   tb <- cbind(pearson=y.cor1, kendall=y.cor2, spearman=y.cor3, MAD=y.mad);
   rownames(tb) <- y.chr

   return(tb);
}

compare_predict1<-function( file.bed, bigwig_histone)
{
cat("BED file", file.bed, "\n");
cat("histone", bigwig_histone, "\n");

   tb.pred <- read.table(file.bed);
   y.pred <- y.bigwig <- c();
   y.chr <- y.cor1 <- y.cor3 <- y.mad <- c();

   for(chr in unique(tb.pred[,1]))
   {
	  file.input <- tempfile(fileext=".bed");
	  file.tmp.bed <- tempfile(fileext=".bed");

      y.label <-  get_histone_read( tb.pred[ tb.pred[,1]==chr, c(1:3)], bigwig_histone);
	  tb.chr <- cbind( tb.pred[ tb.pred[,1]==chr ,1:4], y.label );

      cutoff <- quantile(y.label, 0.95)

	  write.table(tb.chr, file=file.input, quote=F, row.names=F, col.names=F, sep="\t");
	  system( paste("macs2 bdgpeakcall -i ", file.input, "-c ", cutoff, " -o ", file.tmp.bed));
	  tbp.pred <- read.table(file.tmp.bed, skip=1);
      unlink(c(file.tmp.bed, file.input));

	  system( paste("bigWigToBedGraph ", bigwig_histone, " ", file.input, " -chrom=", chr, sep=""));
	  system( paste("macs2 bdgpeakcall -i ", file.input, "-c ", cutoff, " -o ", file.tmp.bed));
	  tbp.org <- read.table(file.tmp.bed, skip=1);
      unlink(c(file.tmp.bed, file.input));

      cat(chr, "cutoff=", cutoff, "NROW.pred=", NROW(tbp.pred), "NROW.org=", NROW(tbp.org), "\n" );

	  ## get overlapped region
	  TPR <- NROW(bed.intersect( tbp.pred, tbp.org ))/NROW(tbp.org);
	  FDR <- 1- NROW(bed.intersect( tbp.pred, tbp.org ))/NROW(tbp.pred);
      cat(chr, "TPR=", TPR,  "FDR=", FDR, "\n" );

	  # get regions overlapped with peaks
      tb.over0 <- bed.intersect( tb.chr, tbp.org );
	  r.cor1 <- cor( tb.over0[,4], tb.over0[,5], method = "pearson" );
      r.cor3 <- cor( tb.over0[,4], tb.over0[,5], method = "spearman" );
      r.mad  <- mad( tb.over0[,4] - tb.over0[,5] );
      cat(chr, "COR=", r.cor1,  r.cor3, r.mad, "\n" );

	  # get regions overlapped with peaks
      tb.nonover <- bed.unintersect( tb.chr, tbp.org );
	  r.cor1 <- cor( tb.nonover[,4], tb.nonover[,5], method = "pearson" );
      r.cor3 <- cor( tb.nonover[,4], tb.nonover[,5], method = "spearman" );
      r.mad  <- mad( tb.nonover[,4] - tb.nonover[,5] );
      cat(chr, "COR=", r.cor1,  r.cor3, r.mad, "\n" );
   }

   return();
}

compare_predict2<-function( file.bed, bigwig_histone, bigwig_histone2=NULL, bigwig_peak=NULL, file.prefix)
{
cat("BED file=", file.bed, "\n");
cat("histone=", bigwig_histone, "\n");
cat("peak=", bigwig_peak, "\n");

   tb.pred <- read.table(file.bed);
   y.pred <- y.bigwig <- c();
   y.chr <- y.cor1 <- y.cor3 <- y.mad <- c();

   tbp.org <- NULL;
   if(!is.null(bigwig_peak))
      tbp.org <- read.table(bigwig_peak);

   for(chr in unique(tb.pred[,1]))
   {
	  file.input <- tempfile(fileext=".bed");
	  file.tmp.bed <- tempfile(fileext=".bed");

      y.label <-  get_histone_read( tb.pred[ tb.pred[,1]==chr, c(1:3)], bigwig_histone);
	  tb.chr <- cbind( tb.pred[ tb.pred[,1]==chr ,1:4], y.label);
      if(!is.null(bigwig_histone2))
      {
		 y.label2 <-  get_histone_read( tb.pred[ tb.pred[,1]==chr, c(1:3)], bigwig_histone2);
	     tb.chr <- cbind(tb.chr, y.label2 );
      }

	  r.cor1 <- cor( tb.chr[,4], tb.chr[,5], method = "pearson" );
      r.cor3 <- cor( tb.chr[,4], tb.chr[,5], method = "spearman" );
      r.mad  <- mad( tb.chr[,4] - tb.chr[,5] );
      cat(chr, "COR=", r.cor1,  r.cor3, r.mad, "\n" );


	png(paste(file.prefix, "all", "png", sep="."), width=900, height=900);
	densScatterplot( tb.chr[,5], tb.chr[,4], main=paste(file.prefix, "(r=", r.cor1, ")"), xlab="histone", ylab="predict");
	dev.off();

	png(paste(file.prefix, "all", "log", "png", sep="."), width=900, height=900);
	densScatterplot( tb.chr[,5], tb.chr[,4], uselog=TRUE, main=paste(file.prefix, "(r=", r.cor1, ")"), xlab="log(histone)", ylab="log(predict)");
	dev.off();

if(!is.null(bigwig_histone2))
{
	png(paste(file.prefix, "all", "uw", "png", sep="."), width=900, height=900);
	densScatterplot( tb.chr[,5], tb.chr[,6], main=paste(file.prefix, "(r=", round(cor(tb.chr[,5], tb.chr[,6]),2), ")"), xlab="Broad", ylab="UW");
	dev.off();

	png(paste(file.prefix, "all", "uw", "log", "png", sep="."), width=900, height=900);
	densScatterplot( tb.chr[,5], tb.chr[,6], main=paste(file.prefix, "(r=", round(cor(tb.chr[,5], tb.chr[,6]),2), ")"), xlab="log(Broad)", ylab="log(UW)");
	dev.off();
}


if(!is.null(tbp.org))
{
	  # get regions overlapped with peaks
      tb.over0 <- bed.intersect( tb.chr, tbp.org );
	  r.cor1 <- cor( tb.over0[,4], tb.over0[,5], method = "pearson" );
      r.cor3 <- cor( tb.over0[,4], tb.over0[,5], method = "spearman" );
      r.mad  <- mad( tb.over0[,4] - tb.over0[,5] );
      cat(chr, "COR=", r.cor1,  r.cor3, r.mad, "\n" );

	png(paste(file.prefix, "peak", "png", sep="."), width=900, height=900);
	densScatterplot( tb.over0[,5], tb.over0[,4], main=paste(file.prefix, "peak(r=", r.cor1, ")"), xlab="histone", ylab="predict");
	dev.off();

	png(paste(file.prefix, "peak", "log", "png", sep="."), width=900, height=900);
	densScatterplot( tb.over0[,5], tb.over0[,4], uselog=TRUE, main=paste(file.prefix, "peak(r=", r.cor1, ")"), xlab="log(histone)", ylab="log(predict)");
	dev.off();


if(!is.null(bigwig_histone2))
{
	png(paste(file.prefix, "peak", "uw", "png", sep="."), width=900, height=900);
	densScatterplot( tb.over0[,5], tb.over0[,6], main=paste(file.prefix, "peak(r=", round(cor(tb.over0[,5], tb.over0[,6]),2), ")"), xlab="Broad", ylab="UW");
	dev.off();

	png(paste(file.prefix, "peak", "uw", "log", "png", sep="."), width=900, height=900);
	densScatterplot( tb.over0[,5], tb.over0[,6], uselog=TRUE, main=paste(file.prefix, "peak(r=", round(cor(tb.over0[,5], tb.over0[,6]),2), ")"), xlab="log(Broad)", ylab="log(UW)");
	dev.off();
}
}


if(!is.null(tbp.org))
{

	  # get regions un-overlapped with peaks
      tb.nonover <- bed.unintersect( tb.chr, tbp.org );
	  r.cor1 <- cor( tb.nonover[,4], tb.nonover[,5], method = "pearson" );
      r.cor3 <- cor( tb.nonover[,4], tb.nonover[,5], method = "spearman" );
      r.mad  <- mad( tb.nonover[,4] - tb.nonover[,5] );
      cat(chr, "COR=", r.cor1,  r.cor3, r.mad, "\n" );

	png(paste(file.prefix, "nonpeak", "png", sep="."), width=900, height=900);
	densScatterplot( tb.nonover[,5], tb.nonover[,4], main=paste(file.prefix, "nonpeak(r=", r.cor1, ")"), xlab="histone", ylab="predict");
	dev.off();

	png(paste(file.prefix, "nonpeak", "log", "png", sep="."), width=900, height=900);
	densScatterplot( tb.nonover[,5], tb.nonover[,4], uselog=TRUE, main=paste(file.prefix, "nonpeak(r=", r.cor1, ")"), xlab="log(histone)", ylab="log(predict)");
	dev.off();

if(!is.null(bigwig_histone2))
{
	png(paste(file.prefix, "nonpeak", "uw", "png", sep="."), width=900, height=900);
	densScatterplot( tb.nonover[,5], tb.nonover[,6], main=paste(file.prefix, "nonpeak(r=", round(cor(tb.nonover[,5], tb.nonover[,6]),2), ")"), xlab="log(Broad)", ylab="log(UW)");
	dev.off();

	png(paste(file.prefix, "nonpeak", "uw", "log", "png", sep="."), width=900, height=900);
	densScatterplot( tb.nonover[,5], tb.nonover[,6], uselog=TRUE, main=paste(file.prefix, "nonpeak(r=", round(cor(tb.nonover[,5], tb.nonover[,6]),2), ")"), xlab="log(Broad)", ylab="log(UW)");
	dev.off();
}
}

   }

   return();
}


source("hist.svm.com.R")
#compare_predict1("/local/ftp/pub/hub/dreg.test/G1-hist/H3k27ac.G1.chr22.bed.gz", "/local/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig" );
#compare_predict1("/local/ftp/pub/hub/dreg.test/G1-hist/K562-H3K27ac9-chr22-imputed.bed", "/local/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig" );
#compare_predict1("./H3K27ac.S1.G1.chr21.bed", "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig" );
#compare_predict2("./H3K27ac.S1.G1.chr21.bed",
#                 "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig",
#                 "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz");

#compare_predict0("H3k27ac.S1.exc21_chr21.bed",
#                 "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig");
#compare_predict2("./H3k27ac.S1.exc21_chr21.bed",
#                 "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig",
#                 "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz");


#compare_predict2("./H3k27ac.S1.exc21_chr21.bed",
#                 "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig",
#                 "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz");


## 3/27/2018
histpath<-function( file.hist )
{
	paste(path.histone, file.hist, sep="/");
}

if(0)
{
	compare_predict2("../predict/H3k122ac.S1.V3.G1_chr22.bed",   histpath(file.H3k122ac.bw), NULL, histpath(file.H3k122ac.peak), "H3k122ac.S1.V3.G1_chr22");
	compare_predict2("../predict/H3k27ac.S1.V3.G1_chr22.bed",    histpath(file.H3k27ac.bw),  NULL, histpath(file.H3k27ac.peak),  "H3k27ac.S1.V3.G1_chr22" );
	compare_predict2("../predict/H3K27ac.S1.exl21.G1_chr22.bed", histpath(file.H3k27ac.bw),  NULL, histpath(file.H3k27ac.peak),  "H3K27ac.S1.exl21.G1_chr22" );
	compare_predict2("../predict/H3K27ac.S1.G1_chr22.bed",       histpath(file.H3k27ac.bw),  NULL, histpath(file.H3k27ac.peak),  "H3K27ac.S1.G1_chr22" );
	compare_predict2("../predict/H3k4me1.S1.V2.G1_chr22.bed",    histpath(file.H3k4me1.bw),  NULL, histpath(file.H3k4me1.peak),  "H3k4me1.S1.V2.G1_chr22" );
	compare_predict2("../predict/H3k4me2.S1.V2.G1_chr22.bed",    histpath(file.H3k4me2.bw),  NULL, histpath(file.H3k4me2.peak),  "H3k4me2.S1.V2.G1_chr22" );
	compare_predict2("../predict/H3k9ac.S1.V2.G1_chr22.bed",     histpath(file.H3k9ac.bw),   NULL, histpath(file.H3k9ac.peak) ,  "H3k9ac.S1.V2.G1_chr22" );
	compare_predict2("../predict/H4k20me1.S1.V3.G1_chr22.bed",   histpath(file.H4k20me1.bw), NULL, histpath(file.H4k20me1.peak), "H4k20me1.S1.V3.G1_chr22");
	compare_predict2("../predict/H3k9me3.S1.V2.G1_chr22.bed",    histpath(file.H3k9me3.bw),  NULL, histpath(file.H3k9me3.peak),  "H3k9me3.S1.V2.G1_chr22" );

	compare_predict2("../predict/H3k4me3.S1.V3.G1_chr22.bed",    histpath(file.H3k4me3.bw),  histpath(file.H3k4me3.bw2),  histpath(file.H3k4me3.peak),  "H3k4me3.S1.V3.G1_chr22" );
	compare_predict2("../predict/H3k27me3.S1.V3.G1_chr22.bed",   histpath(file.H3k27me3.bw), histpath(file.H3k27me3.bw2), histpath(file.H3k27me3.peak), "H3k27me3.S1.V3.G1_chr22");
	compare_predict2("../predict/H3k36me3.S1.V3.G1_chr22.bed",   histpath(file.H3k36me3.bw), histpath(file.H3k36me3.bw2), histpath(file.H3k36me3.peak), "H3k36me3.S1.V3.G1_chr22");
}



compare_peak_cor<-function( file.bed, bigwig_histone, bigwig_histone2=NULL, bigwig_peak, file.prefix)
{
cat("BED file=", file.bed, "\n");
cat("histone=", bigwig_histone, "\n");
cat("peak=", bigwig_peak, "\n");

   if (is.null(bigwig_peak))
   	   return();

   tb.pred <- read.table(file.bed);
   y.pred <- y.bigwig <- c();
   y.chr <- y.cor1 <- y.cor3 <- y.mad <- c();

   tbp.org <- read.table(bigwig_peak);

   for(chr in unique(tb.pred[,1]))
   {
	  file.input <- tempfile(fileext=".bed");
	  file.tmp.bed <- tempfile(fileext=".bed");

      y.label <-  get_histone_read( tb.pred[ tb.pred[,1]==chr, c(1:3)], bigwig_histone);
      y.label2 <- c();
      if(!is.null(bigwig_histone2))
         y.label2 <-  get_histone_read( tb.pred[ tb.pred[,1]==chr, c(1:3)], bigwig_histone2);

	  tb.chr <- cbind( tb.pred[ tb.pred[,1]==chr ,1:4], y.label);
      temp.file <- tempfile(fileext=".bed");
   	  write.table(tb.chr, file=temp.file, quote=F, row.names=F, col.names=F, sep="\t");

      temp.file2 <- tempfile(fileext=".bed");
   	  write.table(tbp.org[tbp.org[,1]==chr, 1:3], file=temp.file2, quote=F, row.names=F, col.names=F, sep="\t");

	  tb.peak1 <- read.table(pipe(paste("bedtools map -a ", temp.file2, " -b ", temp.file, " -c 4 -o mean ")));
	  tb.peak2 <- read.table(pipe(paste("bedtools map -a ", temp.file2, " -b ", temp.file, " -c 5 -o mean ")));

	  if(NROW(tb.peak1)!=NROW(tb.peak1))
	  	  error("Not same size!")

	  r.cor1 <- cor( tb.peak1[,4], tb.peak2[,4], method = "pearson" );
      r.cor3 <- cor( tb.peak1[,4], tb.peak2[,4], method = "spearman" );
      r.mad  <- mad( tb.peak1[,4] - tb.peak2[,4] );
      cat(chr, "COR=", r.cor1,  r.cor3, r.mad, "\n" );

      unlink(c(temp.file2,temp.file));
   }

   return();
}


if(0)
{
#chr 22
compare_peak_cor("../predict/H3k122ac.S1.V3.G1_chr22.bed.gz",   histpath(file.H3k122ac.bw), NULL, histpath(file.H3k122ac.peak), "H3k122ac.S1.V3.G1_chr22");
compare_peak_cor("../predict/H3k27ac.S1.V3.G1_chr22.bed.gz",    histpath(file.H3k27ac.bw),  NULL, histpath(file.H3k27ac.peak),  "H3k27ac.S1.V3.G1_chr22" );
compare_peak_cor("../predict/H3K27ac.S1.exl21.G1_chr22.bed.gz", histpath(file.H3k27ac.bw),  NULL, histpath(file.H3k27ac.peak),  "H3K27ac.S1.exl21.G1_chr22" );
compare_peak_cor("../predict/H3K27ac.S1.G1_chr22.bed.gz",       histpath(file.H3k27ac.bw),  NULL, histpath(file.H3k27ac.peak),  "H3K27ac.S1.G1_chr22" );
compare_peak_cor("../predict/H3k4me1.S1.V2.G1_chr22.bed.gz",    histpath(file.H3k4me1.bw),  NULL, histpath(file.H3k4me1.peak),  "H3k4me1.S1.V2.G1_chr22" );
compare_peak_cor("../predict/H3k4me2.S1.V2.G1_chr22.bed.gz",    histpath(file.H3k4me2.bw),  NULL, histpath(file.H3k4me2.peak),  "H3k4me2.S1.V2.G1_chr22" );
compare_peak_cor("../predict/H3k9ac.S1.V2.G1_chr22.bed.gz",     histpath(file.H3k9ac.bw),   NULL, histpath(file.H3k9ac.peak) ,  "H3k9ac.S1.V2.G1_chr22" );
compare_peak_cor("../predict/H4k20me1.S1.V3.G1_chr22.bed.gz",   histpath(file.H4k20me1.bw), NULL, histpath(file.H4k20me1.peak), "H4k20me1.S1.V3.G1_chr22");
compare_peak_cor("../predict/H3k9me3.S1.V2.G1_chr22.bed.gz",    histpath(file.H3k9me3.bw),  NULL, histpath(file.H3k9me3.peak),  "H3k9me3.S1.V2.G1_chr22" );
compare_peak_cor("../predict/H3k4me3.S1.V3.G1_chr22.bed.gz",    histpath(file.H3k4me3.bw),  NULL, histpath(file.H3k4me3.peak),  "H3k4me3.S1.V3.G1_chr22" );
compare_peak_cor("../predict/H3k27me3.S1.V3.G1_chr22.bed.gz",   histpath(file.H3k27me3.bw), NULL, histpath(file.H3k27me3.peak), "H3k27me3.S1.V3.G1_chr22");
compare_peak_cor("../predict/H3k36me3.S1.V3.G1_chr22.bed.gz",   histpath(file.H3k36me3.bw), NULL, histpath(file.H3k36me3.peak), "H3k36me3.S1.V3.G1_chr22");

#chr 1
compare_peak_cor("../predict/H3k122ac.S1.V3.G1_chr1.bed.gz",   histpath(file.H3k122ac.bw), NULL, histpath(file.H3k122ac.peak), "H3k122ac.S1.V3.G1_chr22");
compare_peak_cor("../predict/H3k27ac.S1.V3.G1_chr1.bed.gz",    histpath(file.H3k27ac.bw),  NULL, histpath(file.H3k27ac.peak),  "H3k27ac.S1.V3.G1_chr22" );
compare_peak_cor("../predict/H3k27me3.S1.V3.G1_chr1.bed.gz",   histpath(file.H3k27me3.bw), NULL, histpath(file.H3k27me3.peak), "H3k27me3.S1.V3.G1_chr22");
compare_peak_cor("../predict/H3k4me1.S1.V2.G1_chr1.bed.gz",    histpath(file.H3k4me1.bw),  NULL, histpath(file.H3k4me1.peak),  "H3k4me1.S1.V2.G1_chr22" );
compare_peak_cor("../predict/H3k4me2.S1.V2.G1_chr1.bed.gz",    histpath(file.H3k4me2.bw),  NULL, histpath(file.H3k4me2.peak),  "H3k4me2.S1.V2.G1_chr22" );
compare_peak_cor("../predict/H3k4me3.S1.V3.G1_chr1.bed.gz",    histpath(file.H3k4me3.bw),  NULL, histpath(file.H3k4me3.peak),  "H3k4me3.S1.V3.G1_chr22" );
compare_peak_cor("../predict/H3k9ac.S1.V2.G1_chr1.bed.gz",     histpath(file.H3k9ac.bw),   NULL, histpath(file.H3k9ac.peak) ,  "H3k9ac.S1.V2.G1_chr22" );
compare_peak_cor("../predict/H4k20me1.S1.V3.G1_chr1.bed.gz",   histpath(file.H4k20me1.bw), NULL, histpath(file.H4k20me1.peak), "H4k20me1.S1.V3.G1_chr22");
compare_peak_cor("../predict/H3k9me3.S1.V2.G1_chr1.bed.gz",    histpath(file.H3k9me3.bw),  NULL, histpath(file.H3k9me3.peak),  "H3k9me3.S1.V2.G1_chr22" );
compare_peak_cor("../predict/H3k36me3.S1.V3.G1_chr1.bed.gz",   histpath(file.H3k36me3.bw), NULL, histpath(file.H3k36me3.peak), "H3k36me3.S1.V3.G1_chr22");


compare_predict2("../predict/H3k122ac.S1.V3.G1_chr1.bed.gz",   histpath(file.H3k122ac.bw), NULL, histpath(file.H3k122ac.peak), "H3k122ac.S1.V3.G1_chr22");
compare_predict2("../predict/H3k27ac.S1.V3.G1_chr1.bed.gz",    histpath(file.H3k27ac.bw),  NULL, histpath(file.H3k27ac.peak),  "H3k27ac.S1.V3.G1_chr22" );
compare_predict2("../predict/H3k27me3.S1.V3.G1_chr1.bed.gz",   histpath(file.H3k27me3.bw), NULL, histpath(file.H3k27me3.peak), "H3k27me3.S1.V3.G1_chr22");
compare_predict2("../predict/H3k4me1.S1.V2.G1_chr1.bed.gz",    histpath(file.H3k4me1.bw),  NULL, histpath(file.H3k4me1.peak),  "H3k4me1.S1.V2.G1_chr22" );
compare_predict2("../predict/H3k4me2.S1.V2.G1_chr1.bed.gz",    histpath(file.H3k4me2.bw),  NULL, histpath(file.H3k4me2.peak),  "H3k4me2.S1.V2.G1_chr22" );
compare_predict2("../predict/H3k4me3.S1.V3.G1_chr1.bed.gz",    histpath(file.H3k4me3.bw),  NULL, histpath(file.H3k4me3.peak),  "H3k4me3.S1.V3.G1_chr22" );
compare_predict2("../predict/H3k9ac.S1.V2.G1_chr1.bed.gz",     histpath(file.H3k9ac.bw),   NULL, histpath(file.H3k9ac.peak) ,  "H3k9ac.S1.V2.G1_chr22" );
compare_predict2("../predict/H4k20me1.S1.V3.G1_chr1.bed.gz",   histpath(file.H4k20me1.bw), NULL, histpath(file.H4k20me1.peak), "H4k20me1.S1.V3.G1_chr22");
compare_predict2("../predict/H3k9me3.S1.V2.G1_chr1.bed.gz",    histpath(file.H3k9me3.bw),  NULL, histpath(file.H3k9me3.peak),  "H3k9me3.S1.V2.G1_chr22" );
compare_predict2("../predict/H3k36me3.S1.V3.G1_chr1.bed.gz",   histpath(file.H3k36me3.bw), NULL, histpath(file.H3k36me3.peak), "H3k36me3.S1.V3.G1_chr22");


}

## CD4
if(0)
{

file.H3k27ac.bigw  <- "/fs/cbsudanko/storage/data/hg19/cd4/epiRoadmap_histone/H3K27ac.bw"
file.H3k27ac.peak  <- NULL
file.H3k27me3.bigw <- "/fs/cbsudanko/storage/data/hg19/cd4/epiRoadmap_histone/H3K27me3.extend.bw"
file.H3k27me3.peak <- NULL
file.H3k36me3.bigw <- "/fs/cbsudanko/storage/data/hg19/cd4/h3k36me3/h3k36me3.extend.bw"
file.H3k36me3.peak <- NULL
file.H3k4me1.bigw  <- "/fs/cbsudanko/storage/data/hg19/cd4/epiRoadmap_histone/H3K4me1.extend.bw"
file.H3k4me1.peak  <- "/fs/cbsudanko/storage/data/hg19/cd4/epiRoadmap_histone/h3k4me1.peaks_peaks.narrowPeak"
file.H3k4me3.bigw  <- "/fs/cbsudanko/storage/data/hg19/cd4/epiRoadmap_histone/H3K4me3.extend.bw"
file.H3k4me3.peak  <- "/fs/cbsudanko/storage/data/hg19/cd4/epiRoadmap_histone/h3k4me3.peaks_peaks.narrowPeak"
file.H3k9ac.bigw   <- "/fs/cbsudanko/storage/data/hg19/cd4/h3k9ac/h3k9ac.extend.bw"
file.H3k9ac.peak   <- "/fs/cbsudanko/storage/data/hg19/cd4/h3k9ac/h3k9ac.peaks_peaks.narrowPeak"
file.H3k9me3.bigw  <- "/fs/cbsudanko/storage/data/hg19/cd4/epiRoadmap_histone/H3K9me3.extend.bw"
file.H3k9me3.peak  <- NULL


compare_peak_cor("../cd4/H3k27me3.S1.V3.CD4_chr22.bed.gz",  file.H3k27me3.bigw, NULL,  file.H3k27me3.peak,  "H3k27me3.S1.V3.CD4.chr22" );
compare_peak_cor("../cd4/H3k36me3.S1.V3.CD4_chr22.bed.gz",  file.H3k36me3.bigw, NULL,  file.H3k36me3.peak,  "H3k36me3.S1.V3.CD4.chr22" );
compare_peak_cor("../cd4/H3k4me1.S1.V2.CD4_chr22.bed.gz",   file.H3k4me1.bigw,  NULL,  file.H3k4me1.peak,   "H3k4me1.S1.V2.CD4.chr22" );
compare_peak_cor("../cd4/H3k4me3.S1.V3.CD4_chr22.bed.gz",   file.H3k4me3.bigw,  NULL,  file.H3k4me3.peak,   "H3k4me3.S1.V3.CD4.chr22" );
compare_peak_cor("../cd4/H3k9ac.S1.V2.CD4_chr22.bed.gz",    file.H3k9ac.bigw,   NULL,  file.H3k9ac.peak,    "H3k9ac.S1.V2.CD4.chr22" );
compare_peak_cor("../cd4/H3k9me3.S1.V2.CD4_chr22.bed.gz",   file.H3k9me3.bigw,  NULL,  file.H3k9me3.peak,   "H3k9me3.S1.V2.CD4.chr22" );


compare_predict2("../cd4/H3k27me3.S1.V3.CD4_chr22.bed.gz",  file.H3k27me3.bigw, NULL,  file.H3k27me3.peak,  "H3k27me3.S1.V3.CD4.chr22" );
compare_predict2("../cd4/H3k36me3.S1.V3.CD4_chr22.bed.gz",  file.H3k36me3.bigw, NULL,  file.H3k36me3.peak,  "H3k36me3.S1.V3.CD4.chr22" );
compare_predict2("../cd4/H3k4me1.S1.V2.CD4_chr22.bed.gz",   file.H3k4me1.bigw,  NULL,  file.H3k4me1.peak,   "H3k4me1.S1.V2.CD4.chr22" );
compare_predict2("../cd4/H3k4me3.S1.V3.CD4_chr22.bed.gz",   file.H3k4me3.bigw,  NULL,  file.H3k4me3.peak,   "H3k4me3.S1.V3.CD4.chr22" );
compare_predict2("../cd4/H3k9ac.S1.V2.CD4_chr22.bed.gz",    file.H3k9ac.bigw,   NULL,  file.H3k9ac.peak,    "H3k9ac.S1.V2.CD4.chr22" );
compare_predict2("../cd4/H3k9me3.S1.V2.CD4_chr22.bed.gz",   file.H3k9me3.bigw,  NULL,  file.H3k9me3.peak,   "H3k9me3.S1.V2.CD4.chr22" );

}

## GM12878
if(0)
{

file.H3k27ac.bigw <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k27acStdSig.bigWig"
file.H3k27ac.peak <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k27acStdPk.broadPeak.gz"
file.H4k20me1.bigw<- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H4k20me1StdSig.bigWig"
file.H4k20me1.peak<- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H4k20me1StdPk.broadPeak.gz"
file.H3k4me1.bigw <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k4me1StdSig.bigWig"
file.H3k4me1.peak <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k4me1StdPk.broadPeak.gz"
file.H3k4me2.bigw <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k4me2StdSig.bigWig"
file.H3k4me2.peak <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k4me2StdPk.broadPeak.gz"
file.H3k4me3.bigw <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k4me3StdSig.bigWig"
file.H3k4me3.peak <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k4me3StdPk.broadPeak.gz"
file.H3k27me3.bigw<- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k27me3StdSig.bigWig"
file.H3k27me3.peak<- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k27me3StdPk.broadPeak.gz"
file.H3k36me3.bigw<- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k36me3StdSig.bigWig"
file.H3k36me3.peak<- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k36me3StdPk.broadPeak.gz"
file.H3k9ac.bigw  <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k9acStdSig.bigWig"
file.H3k9ac.peak  <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k9acStdPk.broadPeak.gz"
file.H3k9me3.bigw <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k9me3StdSig.bigWig"
file.H3k9me3.peak <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k9me3StdPk.broadPeak.gz"

compare_peak_cor("../gm12878/H3k27ac.S1.V3.GM12878_chr22.bed.gz",  file.H3k27ac.bigw,  NULL, file.H3k27ac.peak,  "H3k27ac.S1.V3.GM12878.chr22" );
compare_peak_cor("../gm12878/H3k27me3.S1.V3.GM12878_chr22.bed.gz", file.H3k27me3.bigw, NULL, file.H3k27me3.peak, "H3k27me3.S1.V3.GM12878.chr22" );
compare_peak_cor("../gm12878/H3k9ac.S1.V2.GM12878_chr22.bed.gz",   file.H3k9ac.bigw,   NULL, file.H3k9ac.peak,   "H3k9ac.S1.V2.GM12878.chr22" );
compare_peak_cor("../gm12878/H3k9me3.S1.V2.GM12878_chr22.bed.gz",  file.H3k9me3.bigw,  NULL, file.H3k9me3.peak,  "H3k9me3.S1.V2.GM12878.chr22" );
compare_peak_cor("../gm12878/H3k4me1.S1.V2.GM12878_chr22.bed.gz",  file.H3k4me1.bigw,  NULL, file.H3k4me1.peak,  "H3k4me1.S1.V2.GM12878.chr22" );
compare_peak_cor("../gm12878/H3k4me2.S1.V2.GM12878_chr22.bed.gz",  file.H3k4me2.bigw,  NULL, file.H3k4me2.peak,  "H3k4me2.S1.V2.GM12878.chr22" );
compare_peak_cor("../gm12878/H3k4me3.S1.V3.GM12878_chr22.bed.gz",  file.H3k4me3.bigw,  NULL, file.H3k4me3.peak,  "H3k4me3.S1.V3.GM12878.chr22" );
compare_peak_cor("../gm12878/H4k20me1.S1.V3.GM12878_chr22.bed.gz", file.H4k20me1.bigw, NULL, file.H4k20me1.peak, "H4k20me1.S1.V3.GM12878.chr22" );
compare_peak_cor("../gm12878/H3k36me3.S1.V3.GM12878_chr22.bed.gz", file.H3k36me3.bigw, NULL, file.H3k36me3.peak, "H3k36me3.S1.V3.GM12878.chr22" );

compare_predict2("../gm12878/H3k27ac.S1.V3.GM12878_chr22.bed.gz",  file.H3k27ac.bigw,  NULL, file.H3k27ac.peak,  "H3k27ac.S1.V3.GM12878.chr22" );
compare_predict2("../gm12878/H3k27me3.S1.V3.GM12878_chr22.bed.gz", file.H3k27me3.bigw, NULL, file.H3k27me3.peak, "H3k27me3.S1.V3.GM12878.chr22" );
compare_predict2("../gm12878/H3k9ac.S1.V2.GM12878_chr22.bed.gz",   file.H3k9ac.bigw,   NULL, file.H3k9ac.peak,   "H3k9ac.S1.V2.GM12878.chr22" );
compare_predict2("../gm12878/H3k9me3.S1.V2.GM12878_chr22.bed.gz",  file.H3k9me3.bigw,  NULL, file.H3k9me3.peak,  "H3k9me3.S1.V2.GM12878.chr22" );
compare_predict2("../gm12878/H3k4me1.S1.V2.GM12878_chr22.bed.gz",  file.H3k4me1.bigw,  NULL, file.H3k4me1.peak,  "H3k4me1.S1.V2.GM12878.chr22" );
compare_predict2("../gm12878/H3k4me2.S1.V2.GM12878_chr22.bed.gz",  file.H3k4me2.bigw,  NULL, file.H3k4me2.peak,  "H3k4me2.S1.V2.GM12878.chr22" );
compare_predict2("../gm12878/H3k4me3.S1.V3.GM12878_chr22.bed.gz",  file.H3k4me3.bigw,  NULL, file.H3k4me3.peak,  "H3k4me3.S1.V3.GM12878.chr22" );
compare_predict2("../gm12878/H4k20me1.S1.V3.GM12878_chr22.bed.gz", file.H4k20me1.bigw, NULL, file.H4k20me1.peak, "H4k20me1.S1.V3.GM12878.chr22" );
compare_predict2("../gm12878/H3k36me3.S1.V3.GM12878_chr22.bed.gz", file.H3k36me3.bigw, NULL, file.H3k36me3.peak, "H3k36me3.S1.V3.GM12878.chr22" );

}


calculate_cor <- function (y.org, y.pred , range=2000 )
{
   if(range!=10)
   {
	   n.per.range <- round( range/10 );

	   y.range.org <-  unlist(lapply(1:floor(NROW(y.org)/n.per.range),  function(i) {sum(y.org [(i-1)*n.per.range+c(1:n.per.range)])}));
	   y.range.pred <- unlist(lapply(1:floor(NROW(y.pred)/n.per.range), function(i) {sum(y.pred[(i-1)*n.per.range+c(1:n.per.range)])}));
   }
   else
   {
		y.range.org <- y.org;
		y.range.pred <- y.pred;
   }

   return(list(r.pearson = cor( y.range.org, y.range.pred, method = "pearson" ),
      r.spearman = cor( y.range.org, y.range.pred, method = "spearman" ),
      r.mad = mad( y.range.org -  y.range.pred ) )  );
}

compare_predict_range<-function( file.bed, bigwig_histone, file.prefix, range=2000, est.range=c(10, 20, 50, 100, seq(200, 10000, 200)) )
{
	cat("BED file=", file.bed, "\n");
	cat("histone=", bigwig_histone, "\n");

   tb.pred <- read.table(file.bed);
   y.pred <- y.bigwig <- c();
   y.chr <- y.cor1 <- y.cor3 <- y.mad <- c();

   tbp.org <- NULL;


   for(chr in unique(tb.pred[,1]))
   {
	  file.input <- tempfile(fileext=".bed");
	  file.tmp.bed <- tempfile(fileext=".bed");

      y.org <-  get_histone_read( tb.pred[ tb.pred[,1]==chr, c(1:3)], bigwig_histone);
	  tb.chr <- cbind( tb.pred[ tb.pred[,1]==chr ,1:4], y.org);

     r <- rbindlist( lapply(est.range, function(x) { calculate_cor(y.org, tb.pred[,4], x) }  ) )
show(r);

	pdf(paste(file.prefix, "curve", "pdf", sep="."), width=6, height=6);
	plot(est.range, unlist(r[,1]), type="l", col="red", xlim=range(est.range), ylim=c(-0.2,1) )
	lines(est.range, unlist(r[,2]), col="blue")
	lines(est.range, unlist(r[,3]/max(r[,3])), col="black", lty="22");
	dev.off();

	n.per.range <- round( range/10 );
	y.range.org <-  unlist(lapply(1:floor(NROW(y.org)/n.per.range),  function(i) {sum(y.org [(i-1)*n.per.range+c(1:n.per.range)])}));
	y.range.pred <- unlist(lapply(1:floor(NROW(tb.pred)/n.per.range), function(i) {sum(tb.pred[ (i-1)*n.per.range+c(1:n.per.range),4 ])}));
	r.cor1 = cor( y.range.org, y.range.pred, method = "pearson" );
	r.cor2 = cor( y.range.org, y.range.pred, method = "spearman" );
    r.mad = mad( y.range.org -  y.range.pred );

    cat(chr, "COR=", round(r.cor1,3),  round(r.cor2,3), round(r.mad,3), "\n" );

source("/home/zw355/src/Rplot/denScatter.R");
	png(paste(file.prefix, "all", "png", sep="."), width=900, height=900);
	densScatterplot( y.range.org, y.range.pred, main=paste(file.prefix, "(r=", round(r.cor1,2),"/", round(r.cor2,2), ")"), xlab="histone", ylab="predict");
	dev.off();

	png(paste(file.prefix, "all", "log", "png", sep="."), width=900, height=900);
	densScatterplot( y.range.org, y.range.pred, uselog=TRUE, main=paste(file.prefix, "(r=", round(r.cor1,2),"/", round(r.cor2,2), ")"), xlab="log(histone)", ylab="log(predict)");
	dev.off();
}

}

compare_predict_range("../k562/H3k27ac.S1.V3.G1_chr22.bed.gz",    histpath(file.H3k27ac.bw),  "H3k27ac.S1.V3.G1_chr22", range=2000  );
#chr22 COR= 0.69292 1 352.5257
compare_predict_range("../k562/H3k122ac.S1.V3.G1_chr22.bed.gz",   histpath(file.H3k122ac.bw), "H3k122ac.S1.V3.G1_chr22", range=2000  );
compare_predict_range("../k562/H3K27ac.S1.exl21.G1_chr22.bed.gz", histpath(file.H3k27ac.bw),  "H3K27ac.S1.exl21.G1_chr22", range=2000   );
compare_predict_range("../k562/H3K27ac.S1.G1_chr22.bed.gz",       histpath(file.H3k27ac.bw),  "H3K27ac.S1.G1_chr22" , range=2000  );
compare_predict_range("../k562/H3k4me1.S1.V2.G1_chr22.bed.gz",    histpath(file.H3k4me1.bw),  "H3k4me1.S1.V2.G1_chr22", range=2000   );
compare_predict_range("../k562/H3k4me2.S1.V2.G1_chr22.bed.gz",    histpath(file.H3k4me2.bw),  "H3k4me2.S1.V2.G1_chr22", range=2000   );
compare_predict_range("../k562/H3k9ac.S1.V2.G1_chr22.bed.gz",     histpath(file.H3k9ac.bw),   "H3k9ac.S1.V2.G1_chr22", range=2000   );
compare_predict_range("../k562/H4k20me1.S1.V3.G1_chr22.bed.gz",   histpath(file.H4k20me1.bw), "H4k20me1.S1.V3.G1_chr22", range=2000  );
compare_predict_range("../k562/H3k9me3.S1.V2.G1_chr22.bed.gz",    histpath(file.H3k9me3.bw),  "H3k9me3.S1.V2.G1_chr22", range=2000   );
compare_predict_range("../k562/H3k4me3.S1.V3.G1_chr22.bed.gz",    histpath(file.H3k4me3.bw),  "H3k4me3.S1.V3.G1_chr22", range=2000   );
compare_predict_range("../k562/H3k27me3.S1.V3.G1_chr22.bed.gz",   histpath(file.H3k27me3.bw), "H3k27me3.S1.V3.G1_chr22", range=2000  );
compare_predict_range("../k562/H3k36me3.S1.V3.G1_chr22.bed.gz",   histpath(file.H3k36me3.bw), "H3k36me3.S1.V3.G1_chr22", range=2000  );

