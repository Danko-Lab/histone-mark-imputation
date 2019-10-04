source("/local/workdir/zw355/proj/prj15-histone/scripts/hist.param.R");
source("/local/workdir/zw355/proj/prj15-histone/scripts/hist.svm.com.R")

getProfile <- function( bed1, hMarkFile, dist= 5000, step=25 )
{
	## Load mark.
	hMark <- load.bigWig(hMarkFile)
	mat1 <- NULL;

	if(NROW(bed1)!=0)
	{
		if(NROW(bed1)<100)
		{
			hCountMatrix <- bed.step.bpQuery.bigWig(hMark, center.bed(bed1[,c(1,2,3)], dist, dist), step=step, abs.value=TRUE)
			hmat1 <- matrix(unlist(hCountMatrix), nrow= NROW(bed1), byrow=TRUE);
			mat1 <- colMeans(hmat1);
		}
		else
		{
			hmat1 <- metaprofile.bigWig(center.bed(bed1[,c(1,2,3)], dist, dist), hMark, step=step)
			mat1 <- abs(hmat1$middle)
		}
	}

    unload.bigWig(hMark);

    return(mat1)
}


writeMatplot<- function(file.pdf, file.bed1, hMarkFile, hPredFile, title, subs= NULL, breaks= NULL,  cols= NULL, dist= 5000, step=50, chrs=NULL)
{
    dregX <- read.table(file.bed1)
    #dregX <- dregX[dregX$V4>0.9 &dregX$V4<1.2, ];
    dregX <- dregX[,c(1,6,6)];
	dregX <- dregX[grep("_|chrY|chrX", dregX[,1], invert=TRUE),]
	dregX <- dregX[ dregX[,2]> dist,]

	if(!is.null(chrs))
	   dregX <- dregX[dregX$V1 %in% chrs,]

	mat1 <- getProfile( dregX, hMarkFile, dist, step )
	mat2 <- getProfile( dregX, hPredFile, dist, step )*10
    
    pdf(file.pdf)
	par( mar=c(5,4,2,2), plt=c(0.15, 0.99, 0.2, 0.8), mgp=c(1.5,0.4,0) );

	if(!is.null(mat1) || !is.null(mat2))
		y.max <- max(c(mat1, mat2))
	else
		y.max <- 10;

	plot(NA, NA, type="n", col="gray", xlim=c(1,ifelse(is.null(mat1), NROW(mat2), NROW(mat1))), ylim=c(0, y.max)/1000,  xlab="Distance (Kbp) ", ylab="", cex.axis=1/3, cex.lab=1/3,  cex.main=1/3, xaxt = "n", yaxt="n", main=title)

	if(!is.null(mat1)) lines(1:NROW(mat1), mat1/1000, col="blue", lwd=1/8);
	if(!is.null(mat2)) lines(1:NROW(mat2), mat2/1000, col="black", lwd=1/8);

	axis(1, c(0, NROW(mat1)/4, NROW(mat1)/2, NROW(mat1)*3/4, NROW(mat1) ), c(-dist/1000, -dist/2000, 0, dist/2000, dist/1000), cex.axis=1/3, cex=1/3 )
	if( y.max>10 )
		axis(2, c(0, y.max/2000, y.max/1000), c(0, round(y.max/2), round(y.max)), cex.axis=1/3, cex=1/3 )
	else if( y.max>=1 )
		axis(2, c(0, y.max/2000, y.max/1000), c(0, round(y.max/2, digits=1), round(y.max,digits=1)), cex.axis=1/3, cex=1/3 )
	else
		axis(2, c(0, y.max/2000, y.max/1000), c(0, signif(y.max/2,digits=2), signif(y.max,digits=2)), cex.axis=1/3, cex=1/3 )

    dev.off();
	return()
}

if(0)
{
# GM 
file.GM.bed <- "/local/workdir/zw355/proj/prj10-dreg/new-rf-201803/GM/GM.dREG.peak.full.bed.gz"
file.TSS.bed <- file.gm.H3k27ac.peak;
file.TSS.bed <- file.GM.bed;

file.pdf <- "metaplot-GM-H3K27ac.pdf"
writeMatplot( file.pdf,  file.TSS.bed, file.gm.H3k27ac.bw,  file.gm.H3k27ac.pred.bw, "H3k27ac" );
file.pdf <- "metaplot-GM-H3k27me3.pdf"
writeMatplot( file.pdf,  file.TSS.bed, file.gm.H3k27me3.bw, file.gm.H3k27me3.pred.bw, "H3k27me3" );
file.pdf <- "metaplot-GM-H3k36me3.pdf"
writeMatplot( file.pdf,  file.TSS.bed, file.gm.H3k36me3.bw, file.gm.H3k36me3.pred.bw, "H3k36me3" );
file.pdf <- "metaplot-GM-H3k4me1.pdf"
writeMatplot( file.pdf,  file.TSS.bed, file.gm.H3k4me1.bw,  file.gm.H3k4me1.pred.bw, "H3k4me1" );
file.pdf <- "metaplot-GM-H3k4me2.pdf"
writeMatplot( file.pdf,  file.TSS.bed, file.gm.H3k4me2.bw,  file.gm.H3k4me2.pred.bw, "H3k4me2" );
file.pdf <- "metaplot-GM-H3k4me3.pdf"
writeMatplot( file.pdf,  file.TSS.bed, file.gm.H3k4me3.bw,  file.gm.H3k4me3.pred.bw, "H3k4me3" );
file.pdf <- "metaplot-GM-H3k9ac.pdf"
writeMatplot( file.pdf,  file.TSS.bed, file.gm.H3k9ac.bw,  file.gm.H3k9ac.pred.bw, "H3k9ac" );
file.pdf <- "metaplot-GM-H3k9me3.pdf"
writeMatplot( file.pdf,  file.TSS.bed, file.gm.H3k9me3.bw,  file.gm.H3k9me3.pred.bw, "H3k9me3" );
file.pdf <- "metaplot-GM-H4k20me1.pdf"
writeMatplot( file.pdf,  file.TSS.bed, file.gm.H4k20me1.bw,  file.gm.H4k20me1.pred.bw, "H4k20me1" );

file.dreg.bed <- "/local/workdir/zw355/proj/prj10-dreg/new-rf-201803/G1/G1.dREG.peak.full.bed.gz"
file.TSS.bed <- file.dreg.bed;
file.pdf <- "metaplot-K562-Alex-Mnase-H3K27ac.pdf"
writeMatplot( file.pdf,  file.TSS.bed, file.Alex.Mnase.H3k27ac.bw,  file.Alex.Mnase.H3k27ac.pred.bw, "H3k27ac", chrs="chr22" );
file.pdf <- "metaplot-K562-Alex-Mnase-H3k4me1.pdf"
writeMatplot( file.pdf,  file.TSS.bed, file.Alex.Mnase.H3k4me1.bw,  file.Alex.Mnase.H3k4me1.pred.bw, "H3k4me1", chrs="chr22" );
file.pdf <- "metaplot-K562-Alex-Mnase-H3k4me2.pdf"
writeMatplot( file.pdf,  file.TSS.bed, file.Alex.Mnase.H3k4me2.bw,  file.Alex.Mnase.H3k4me2.pred.bw, "H3k4me2", chrs="chr22" );
file.pdf <- "metaplot-K562-Alex-Mnase-H3k4me3.pdf"
writeMatplot( file.pdf,  file.TSS.bed, file.Alex.Mnase.H3k4me3.bw,  file.Alex.Mnase.H3k4me3.pred.bw, "H3k4me3", chrs="chr22" );
}


writeMultiInMatplot<- function(file.pdf, bed.TSS, hMarkFiles, hPredFiles, marker.names, cols, subs= NULL, breaks= NULL,  dist= 5000, step=50, chrs=NULL)
{
    ##dregX <- dregX[dregX$V4>0.9 &dregX$V4<1.2, ];
    #dregX <- dregX[,c(1,6,6)];
	bed.TSS <- bed.TSS[grep("_|chrY|chrX|chrM", bed.TSS[,1], invert=TRUE),]
	bed.TSS <- bed.TSS[ bed.TSS[,2]> dist,]

	if(!is.null(chrs))
	   bed.TSS <- bed.TSS[bed.TSS$V1 %in% chrs,]

	y.max <- 10;
    for(i in 1:NROW(hMarkFiles))
    {
		mat1 <- getProfile( bed.TSS, hMarkFiles[i], dist, step )
		mat2 <- getProfile( bed.TSS, hPredFiles[i], dist, step )*10

        x.max <- ifelse(is.null(mat1), NROW(mat2), NROW(mat1));
        
		if(!is.null(mat1) || !is.null(mat2))
			y.max <- max(y.max, c(mat1), c(mat2))
cat("i=", i, "x.max=", x.max, "y.max=", y.max, "\n");	
	}
	

    pdf(file.pdf);

	plot(NA, NA, type="n", col="gray", xlim=c(0,x.max), ylim=c(0, y.max*1.1)/1000,  xlab="Distance (Kbp) ", ylab="Signal or Read", cex.axis=1, cex.lab=1,  cex.main=1, xaxt = "n", main="")
	par( mar=c(5,4,2,2), mgp=c(1.5,0.4,0) ); #plt=c(0.15, 0.99, 0.2, 0.8), 

    for(i in 1:NROW(hMarkFiles))
    {
		mat11 <- getProfile( bed.TSS[bed.TSS[,4]=="+",], hMarkFiles[i], dist, step )
		mat10 <- getProfile( bed.TSS[bed.TSS[,4]=="-",], hMarkFiles[i], dist, step )
		mat21 <- getProfile( bed.TSS[bed.TSS[,4]=="+",], hPredFiles[i], dist, step )*10
		mat20 <- getProfile( bed.TSS[bed.TSS[,4]=="-",], hPredFiles[i], dist, step )*10
        
        
        mat1 <- mat11*sum(bed.TSS[,4]=="+")/NROW(bed.TSS)+rev(mat10)*sum(bed.TSS[,4]=="-")/NROW(bed.TSS)
        mat2 <- mat21*sum(bed.TSS[,4]=="+")/NROW(bed.TSS)+rev(mat20)*sum(bed.TSS[,4]=="-")/NROW(bed.TSS)
    
		if(!is.null(mat1)) lines(1:NROW(mat1), mat1/1000, col=cols[i], lwd=1.2);
		if(!is.null(mat2)) lines(1:NROW(mat2), mat2/1000, col=cols[i], lwd=1.5, lty="22");
     }
  
	axis(1, c(0, NROW(mat1)/4, NROW(mat1)/2, NROW(mat1)*3/4, NROW(mat1) ), c(-dist/1000, -dist/2000, 0, dist/2000, dist/1000), cex.axis=1, cex=1 )
	#if( y.max>10 )
	#	axis(2, c(0, y.max/2000, y.max/1000), c(0, round(y.max/2), round(y.max)), cex.axis=1, cex=1 )
	#else if( y.max>=1 )
	#	axis(2, c(0, y.max/2000, y.max/1000), c(0, round(y.max/2, digits=1), round(y.max,digits=1)), cex.axis=1, cex=1 )
	#else
	#	axis(2, c(0, y.max/2000, y.max/1000), c(0, signif(y.max/2,digits=2), signif(y.max,digits=2)), cex.axis=1, cex=1 )

    legend("topright", c(marker.names), cex=1.0,  bty="n", fill=cols);

    dev.off();

	return()
}

getTssWithStrand <- function(file.TSS.bed)
{
  library(sqldf)
  tb <- read.table(file.TSS.bed);
  tb[,2] <- round((tb[,2]+tb[,3])/2);
  tb[,3] <- tb[,2]+1

  bed.TSS <- c();
  for(i in 1:200)
  {
    file.GM.temp.tss <- write.temp.bed(tb);
    err <- try( tmpTSS0 <- read.table(pipe(paste("bedtools intersect -a", file.GM.temp.tss, "-b", file.bed.strand, "-wa -wb" ))) , silent=TRUE);
    if(class(err)=="try-error")
    {
      tb[,2] <- tb[,2]-10;
      tb[,3] <- tb[,3]+10;
      next;
    }
    
    tmpTSS <- unique(tmpTSS0[,c(1,2,3,7)])
    bedsum <- sqldf("select V1, V2, V3, count(V7) from tmpTSS group by V1, V2, V3")
    tb <- tb[is.na(match(paste(tb$V1, tb$V2, sep=":"), paste(tmpTSS$V1, tmpTSS$V2, sep=":"))),]
    tb[,2] <- tb[,2]-10;
    tb[,3] <- tb[,3]+10;

    bedUnique <- bedsum[bedsum[,4]==1,]
    bed.TSS <- rbind(bed.TSS, tmpTSS[!is.na(match(paste(tmpTSS$V1, tmpTSS$V2, sep=":"), paste(bedUnique$V1, bedUnique$V2, sep=":"))),]); 

    if(NROW(tb)<10)
       break;
  }

  return(read.table(write.temp.bed(bed.TSS)));
}


if(0)
{
bed.strand <- read.table(pipe("zcat /fs/cbsudanko/storage/data/hg19/all/gencode/gencode.v19.annotation.gtf.gz | awk 'BEGIN{OFS=\"\\t\"}{if ($3==\"transcript\") print $0}' - | grep  \"transcript_type \\\"protein_coding\\\"\" - | awk 'BEGIN{OFS=\"\\t\"}{print $1, $4, $5, $7}' - "))
file.bed.strand <- write.temp.bed(bed.strand);

#file.GM.bed <- "/local/workdir/zw355/proj/prj10-dreg/new-rf-201803/GM/GM.dREG.peak.full.bed.gz"
#bed.TSS <- read.table(file.GM.bed);
file.GM.Tss.bed <- "/local/workdir/zw355/proj/prj10-dreg/GM12878/hg19.GM12878.grocap.pair.bed"
bed.TSS <- getTssWithStrand(file.GM.Tss.bed)
print(NROW(bed.TSS)/NROW(read.table(file.GM.Tss.bed)));
#0.80

file.pdf <- "metaplot-GM-H3K27ac-3k4me-3k9ac.pdf"
writeMultiInMatplot( file.pdf,  bed.TSS, 
                c(file.gm.H3k27ac.bw,  file.gm.H3k4me1.bw, file.gm.H3k4me2.bw, file.gm.H3k4me3.bw, file.gm.H3k9ac.bw),
                c(file.gm.H3k27ac.pred.bw,  file.gm.H3k4me1.pred.bw, file.gm.H3k4me2.pred.bw, file.gm.H3k4me3.pred.bw, file.gm.H3k9ac.pred.bw),
                c("H3k27ac",  "H3k4me1", "H3k4me2", "H3k4me3", "H3k9ac"),
                c("black",  "#cb5b42", "#6985cd", "#cbf21f",  "purple") );


file.G1.Tss.bed <- "/local/workdir/zw355/proj/prj10-dreg/k562/hg19.k562.grocap.pair.bed"
bed.TSS <- getTssWithStrand(file.G1.Tss.bed)
print(NROW(bed.TSS)/NROW(read.table(file.G1.Tss.bed)));
#0.75

file.pdf <- "metaplot-Alex-Mnase-H3k4me.pdf"
writeMultiInMatplot( file.pdf,  bed.TSS, 
                c(file.Alex.Mnase.H3k4me1.bw, file.Alex.Mnase.H3k4me2.bw, file.Alex.Mnase.H3k4me3.bw),
                c(file.Alex.Mnase.H3k4me1.pred.bw, file.Alex.Mnase.H3k4me2.pred.bw, file.Alex.Mnase.H3k4me3.pred.bw),
                c("H3k4me1", "H3k4me2", "H3k4me3" ),
                c("#cb5b42", "#6985cd", "#cbf21f" ), chrs="chr22" );
}


getStrechProfile <- function( bed1, hMarkFile, step=100, points=100 )
{
	## Load mark.
	hMark <- load.bigWig(hMarkFile)
	mat1 <- NULL;
	
	bed1[,2] <- round(bed1[,2]/10)*10 - 10;
	bed1[,3] <- round(bed1[,3]/10)*10 + 50;

	r <- lapply(1:NROW(bed1), function(i){

		hmat1 <- bed.step.bpQuery.bigWig(hMark, bed1[i,c(1,2,3)], step=step, abs.value=TRUE)
		meanvec <- abs(hmat1[[1]])
		if(bed1[i,4]=="-")
		  meanvec <- rev(meanvec);
		meanvec <- approx(1:NROW(meanvec), meanvec, n=points)
		return(meanvec$y);
	});

    mat1 <- do.call("rbind",r);

    unload.bigWig(hMark);

    return(colMeans(mat1))
}



writeStrectMatplot<- function(file.pdf, bed.strand, hMarkFiles, hPredFiles, marker.names, cols, subs= NULL, breaks= NULL,  chrs=NULL)
{
    ##dregX <- dregX[dregX$V4>0.9 &dregX$V4<1.2, ];
    #dregX <- dregX[,c(1,6,6)];
	bed.strand <- bed.strand[grep("_|chrY|chrX|chrM", bed.strand[,1], invert=TRUE),]
	bed.strand <- bed.strand[ bed.strand[,3] - bed.strand[,2]>300,]

	if(!is.null(chrs))
	   bed.strand <- bed.strand[bed.strand$V1 %in% chrs,]

    xsize = 100;
	y.max <- 10;
    for(i in 1:NROW(hMarkFiles))
    {
		mat1 <- getStrechProfile( bed.strand, hMarkFiles[i], points=xsize)
		mat2 <- getStrechProfile( bed.strand, hPredFiles[i], points=xsize)*10

        x.max <- ifelse(is.null(mat1), NROW(mat2), NROW(mat1));
        
		if(!is.null(mat1) || !is.null(mat2))
			y.max <- max(y.max, c(mat1), c(mat2))
cat("i=", i, "x.max=", x.max, "y.max=", y.max, "\n");	
	}
	

    pdf(file.pdf);

	plot(NA, NA, type="n", col="gray", xlim=c(0,xsize), ylim=c(0, y.max*1.1)/1000,  xlab="Gene Length ", ylab="Signal", cex.axis=1, cex.lab=1,  cex.main=1, xaxt = "n", main="")
	par( mar=c(5,4,2,2), mgp=c(1.5,0.4,0) ); #plt=c(0.15, 0.99, 0.2, 0.8), 

    for(i in 1:NROW(hMarkFiles))
    {
		mat1 <- getStrechProfile( bed.strand, hMarkFiles[i], points=xsize)
		mat2 <- getStrechProfile( bed.strand, hPredFiles[i], points=xsize)*10
        
		if(!is.null(mat1)) lines(1:NROW(mat1), mat1/1000, col=cols[i], lwd=1.2);
		if(!is.null(mat2)) lines(1:NROW(mat2), mat2/1000, col=cols[i], lwd=1.5, lty="22");
     }
  
	axis(1, c(0, NROW(mat1)/4, NROW(mat1)/2, NROW(mat1)*3/4, NROW(mat1) ), c(0, 0.25, 0.5, 0.75, 1), cex.axis=1, cex=1 )
	#if( y.max>10 )
	#	axis(2, c(0, y.max/2000, y.max/1000), c(0, round(y.max/2), round(y.max)), cex.axis=1, cex=1 )
	#else if( y.max>=1 )
	#	axis(2, c(0, y.max/2000, y.max/1000), c(0, round(y.max/2, digits=1), round(y.max,digits=1)), cex.axis=1, cex=1 )
	#else
	#	axis(2, c(0, y.max/2000, y.max/1000), c(0, signif(y.max/2,digits=2), signif(y.max,digits=2)), cex.axis=1, cex=1 )

    legend("topright", c(marker.names), cex=1.0,  bty="n", fill=cols);

    dev.off();

	return()
}

if(0)
{
bed.strand <- read.table(pipe("zcat /fs/cbsudanko/storage/data/hg19/all/gencode/gencode.v19.annotation.gtf.gz | awk 'BEGIN{OFS=\"\\t\"}{if ($3==\"gene\") print $0}' - | grep  \"gene_type \\\"protein_coding\\\"\" - | awk 'BEGIN{OFS=\"\\t\"}{print $1, $4, $5, $7}' - "))
file.bed.strand <- write.temp.bed(bed.strand);


file.pdf <- "stretch-metaplot-GM-H3k4me1-4k20me1-3k36me3.pdf"
writeStrectMatplot( file.pdf,  bed.strand, 
                c(file.gm.H3k4me1.bw,  file.gm.H4k20me1.bw, file.gm.H3k36me3.bw ),
                c(file.gm.H3k4me1.pred.bw,  file.gm.H4k20me1.pred.bw, file.gm.H3k36me3.pred.bw ),
                c("H3k4me1",  "H4k20me1", "H3k36me3"),
                c("black",  "#cb5b42", "#6985cd") );


file.pdf <- "stretch-metaplot-Alex-Mnase-H3k4me.pdf"
writeStrectMatplot( file.pdf,  bed.strand, 
                c(file.Alex.Mnase.H3k4me1.bw, file.Alex.Mnase.H3k4me2.bw, file.Alex.Mnase.H3k4me3.bw),
                c(file.Alex.Mnase.H3k4me1.pred.bw, file.Alex.Mnase.H3k4me2.pred.bw, file.Alex.Mnase.H3k4me3.pred.bw),
                c("H3k4me1", "H3k4me2", "H3k4me3" ),
                c("#cb5b42", "#6985cd", "#cbf21f" ), chrs="chr22" );
}
