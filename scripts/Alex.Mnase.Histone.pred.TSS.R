source("hist.svm.com.R")
source("hist.svm.pred.R")
source("hist.svm.main.R")
source("hist.param.R");


## H3k4me1
if(0)
{
file.bw.histone  <- file.Alex.Mnase.H3k4me1.bw;
file.peak.histone<- file.Alex.Mnase.H3k4me1.peak;
file.rdata.model <- "../models/Alex2.Mnase.H3k4me1.S1.train.rdata";
file.TSS.pred.gz <- "../pred-alex/Alex-Mnase-H3k4me1.pred.TSS.bed.gz"
}

## H3k4me2
if(0)
{
file.bw.histone  <- file.Alex.Mnase.H3k4me2.bw;
file.peak.histone<- file.Alex.Mnase.H3k4me2.peak;
file.rdata.model <- "../models/Alex2.Mnase.H3k4me2.S1.train.rdata";
file.TSS.pred.gz <- "../pred-alex/Alex-Mnase-H3k4me2.pred.TSS.bed.gz"
}

## H3k4me3
if(0)
{
file.bw.histone  <- file.Alex.Mnase.H3k4me3.bw;
file.peak.histone<- file.Alex.Mnase.H3k4me3.peak;
file.rdata.model <- "../models/Alex2.Mnase.H3k4me3.S1.train.rdata";
file.TSS.pred.gz <- "../pred-alex/Alex-Mnase-H3k4me3.pred.TSS.bed.gz"
}

## H3k27ac
if(0)
{
file.bw.histone  <- "/workdir/zw355/proj/prj15-histone/pred-alex/Zhg19_K27ac_merged.bw"
file.peak.histone<- "/workdir/zw355/proj/prj15-histone/narrowpeaks/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz";
file.rdata.model <- "../models/Alex2.Mnase.H3k27ac.S1.train.rdata";
file.TSS.pred.gz <- "../pred-alex/Alex-Mnase-H3k27ac.pred.TSS.bed.gz"
}

## H3k36me3
if(1)
{
file.bw.histone  <- file.Alex.Mnase.H3k36me3.bw;
file.peak.histone<- file.Alex.Mnase.H3k36me3.peak;
file.rdata.model <- "../models/Alex.Mnase.H3k36me3.S1.train.rdata";
file.TSS.pred.gz <- "../pred-alex/Alex-Mnase-H3k36me3.pred.TSS.sup.bed.gz"
}

## H3k79me3
if(1)
{
file.bw.histone  <- file.Alex.Mnase.H3k79me3.bw;
file.peak.histone<- file.Alex.Mnase.H3k79me3.peak;
file.rdata.model <- "../models/Alex.Mnase.H3k79me3.S1.train.rdata";
file.TSS.pred.gz <- "../pred-alex/Alex-Mnase-H3k79me3.pred.TSS.sup.bed.gz"
}


file.TSS.bed <- "/workdir/zw355/proj/prj15-histone/pred-k562/K562_DivergentPairs.bed"
file.bws.plus <- c(file.bw.G1[1], file.bw.G2[1],file.bw.G3[1], file.bw.G5[1], file.bw.G6[1]);
file.bws.minus <-c(file.bw.G1[2], file.bw.G2[2],file.bw.G3[2], file.bw.G5[2], file.bw.G6[2]);
path.histone <- ""
path.proseq  <- ""
ncores <- 5

selectGPUdevice(1);

predict_TSS_bed<-function()
{
	tb <- read.table(file.TSS.bed, skip=1)
	bed.Tss <- c(); 
	bed.region.temp<-c();
	
	for(i in 1:((NROW(tb)/2-1)))
	{
		bed.region.temp <- c(bed.region.temp, i);

		r0 <- tb[i*2+c(-1:0),];
		r0.max <- r0[ which.max(r0[,5]) , ,drop=F];
		r0.min <- r0[ which.min(r0[,5]) , ,drop=F];
		r1 <- tb[i*2+c(1:2),];
		r1.max <- r1[ which.max(r1[,5]) , ,drop=F];
		r1.min <- r1[ which.min(r1[,5]) , ,drop=F];

		gap <- max(r0.max[,2:3], r0.min[,2:3] ) - min(r1.max[,2:3], r1.min[,2:3] );
		if(abs(gap)<=600)
		{   
			cat("Remove small region", i, gap, "\n");
			next;
		}
		else
		{
			r0 <- tb[bed.region.temp*2+c(-1:0),];
			r0.max <- r0[ which.max(r0[,5]) , ,drop=F];
			r0.max.center <- (r0.max [1, 3] + r0.max [1, 2])/2

			r0 <- r0[ r0[,6] != r0.max[1,6],, drop=F];
			r0.min <- r0[ which.min(r0[,5]), ,drop=F];
			r0.min.center <- (r0.min[1, 2]+r0.min[1, 3])/2

			r0.dist <- r0.max.center - r0.min.center;

			r0 <- rbind(r0.max, r0.min);
			r0.plus <- r0[ r0[,6]=="+", ,drop=F];
			r0.minus <- r0[ r0[,6]=="-", ,drop=F];

			if(r0.plus[1,3] > r0.minus[1,3])
			   r.type = 1 #divergent
			else
			   r.type = -1  #divergent 

			bed.Tss<-rbind(bed.Tss, data.frame(chr=r0.max[1,1], start=r0.max.center-1000, 
				  end=r0.max.center+1000, 
				  dist=r0.dist, 
				  score=r0.max[1,5], 
				  strand=r0.max[1,6],
				  type=r.type));

		   bed.region.temp <-c();
	  }
	}

	bed.Tss$dist <- abs(bed.Tss$dist)
	bed.Tss <- bed.Tss[order(bed.Tss$chr, bed.Tss$start),]
	bed.Tss[,2] <- bed.Tss[,2]-100
	bed.Tss[,3] <- bed.Tss[,3]+100

cat("NROW(bed.Tss)", NROW(bed.Tss), "\n");

	bed.pred <- do.call("rbind", lapply(1:NROW(bed.Tss), function(i)
	{
		df <- data.frame(chr=bed.Tss[i,1], start=seq(round(bed.Tss[i,2]/10 -1 )*10, round(bed.Tss[i,3]/10 + 1)*10, 10));
		df$stop <- df$start+1;
		return(df);
	}));

	bed.pred <- unique(bed.pred);

cat("NROW(bed.pred)", NROW(bed.pred), "\n");
cat("Model:", file.rdata.model, "\n");

	load(file.rdata.model);

	x1 <- extract_feature_matrix( bed.pred, 
	          paste(model$src$path.proseq, model$src$file.bws.plus[1], sep="/"), 
	          paste(model$src$path.proseq, model$src$file.bws.minus[1], sep="/"), 
	          model$gdm, linear_scale=F, ncores=ncores );

	y_pred <- Rgtsvm::predict.gtsvm(model$svm, x1$mat);

cat("Bed file:", file.TSS.pred.gz, "\n");
	write.bed( data.frame(bed.pred, y_pred ), file.TSS.pred.gz, compress=TRUE);
	return(file.TSS.pred.gz);
}

predict_TSS_bed();

