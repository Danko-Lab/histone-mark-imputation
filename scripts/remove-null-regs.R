library(parallel);
library(data.table);
library(bigWig);

source("../scripts/hist.svm.com.R");
source("../scripts/hist.param.R");

remove_null_reg <- function(file.bed.gz, rate=1)
{
	tb <- read.table(file.bed.gz);
	tbx <- tb[,4];
	
	#tbx[ round(tbx) <= 3 ] <- 1;

	tb_10<- c( rep(tbx[1], 10), tbx[-NROW(tbx)+c(0:9)]);
	tb_9 <- c( rep(tbx[1], 9), tbx[-NROW(tbx)+c(0:8)]);
	tb_8 <- c( rep(tbx[1], 8), tbx[-NROW(tbx)+c(0:7)]);
	tb_7 <- c( rep(tbx[1], 7), tbx[-NROW(tbx)+c(0:6)]);
	tb_6 <- c( rep(tbx[1], 6), tbx[-NROW(tbx)+c(0:5)]);
	tb_5 <- c( rep(tbx[1], 5), tbx[-NROW(tbx)+c(0:4)]);
	tb_4 <- c( rep(tbx[1], 4), tbx[-NROW(tbx)+c(0:3)]);
	tb_3 <- c( rep(tbx[1], 3), tbx[-NROW(tbx)+c(0:2)]);
	tb_2 <- c( rep(tbx[1], 2), tbx[-NROW(tbx)+c(0:1)]);
	tb_1 <- c( tbx[1], tbx[-NROW(tbx)]);
	tb0 <- tbx;
	tb1 <- c( tbx[-1], tbx[NROW(tbx)]);
	tb2 <- c( tbx[-c(1:2)], rep( tbx[NROW(tbx)],2) );
	tb3 <- c( tbx[-c(1:3)], rep( tbx[NROW(tbx)],3) );
	tb4 <- c( tbx[-c(1:4)], rep( tbx[NROW(tbx)],4) );
	tb5 <- c( tbx[-c(1:5)], rep( tbx[NROW(tbx)],5) );
	tb6 <- c( tbx[-c(1:6)], rep( tbx[NROW(tbx)],6) );
	tb7 <- c( tbx[-c(1:7)], rep( tbx[NROW(tbx)],7) );
	tb8 <- c( tbx[-c(1:8)], rep( tbx[NROW(tbx)],8) );
	tb9 <- c( tbx[-c(1:9)], rep( tbx[NROW(tbx)],9) );
	tb10<- c( tbx[-c(1:10)],rep( tbx[NROW(tbx)],10) );

	b0 <- ( tb0==tb1 & tb0==tb2  & tb0==tb3  & tb0==tb4  & tb0==tb5  & tb0==tb6) 
	b1 <- (tb_1==tb0 & tb0==tb1  & tb0==tb2  & tb0==tb3  & tb0==tb4  & tb0==tb5) 
	b2 <- (tb_2==tb0 & tb_1==tb0 & tb0==tb1  & tb0==tb2  & tb0==tb3  & tb0==tb4)  
	b3 <- (tb_3==tb0 & tb_2==tb0 & tb_1==tb0 & tb0==tb1  & tb0==tb2  & tb0==tb3)  
	b4 <- (tb_4==tb0 & tb_3==tb0 & tb_2==tb0 & tb_1==tb0 & tb0==tb1  & tb0==tb2)  
	b5 <- (tb_5==tb0 & tb_4==tb0 & tb_3==tb0 & tb_2==tb0 & tb_1==tb0 & tb0==tb1)  
	b6 <- (tb_6==tb0 & tb_5==tb0 & tb_4==tb0 & tb_3==tb0 & tb_2==tb0 & tb_1==tb0)  
; 
	idx.rem <- which( b0+b1+b2+b3+b4+b5+b6>0 & tb0 < 4.0/rate);
	if (NROW(idx.rem)>0) 
	   tb <- tb[ -idx.rem,,drop=F]

    idx <- which(tb[,4] <5/rate);
    tb[idx,4] <- (tb[idx,4])^2/5

    idx <- which(tb[,4] <4/rate);
    tb[idx,4] <- (tb[idx,4])^2/4
   
    idx.rem <- which(tb[,4] < 0.01/rate);
	if (NROW(idx.rem)>0) 
	   tb <- tb[ -idx.rem,,drop=F]
    
    #mapped <- get_histone_read(tb, file.unmap.bw)
    #tb <- tb[mapped==0,] ;

    file.tmp <- write.temp.bed(tb);

    tb <- read.table(pipe(paste("bedmap --echo --delim '\t' --indicator ", file.tmp, file.removed.reg )));
    if (sum(tb$V5==0)>0) 
        tb <- tb[tb$V5==0,,drop=F];
    tb <- tb[,c(1:4),drop=F] 

    unlink(file.tmp);
    file.tmp <- NULL;
    if(NROW(tb)>0)
        file.tmp <- write.temp.bed(tb, compress=TRUE);

    gc();
	
	return( file.tmp );
}

remove_null_reg10 <- function(file.bed.gz)
{
	tb <- read.table(file.bed.gz);
	tbx <- tb[,4];

	tb_10<- c( rep(tbx[1], 10), tbx[-NROW(tbx)+c(0:9)]);
	tb_9 <- c( rep(tbx[1], 9), tbx[-NROW(tbx)+c(0:8)]);
	tb_8 <- c( rep(tbx[1], 8), tbx[-NROW(tbx)+c(0:7)]);
	tb_7 <- c( rep(tbx[1], 7), tbx[-NROW(tbx)+c(0:6)]);
	tb_6 <- c( rep(tbx[1], 6), tbx[-NROW(tbx)+c(0:5)]);
	tb_5 <- c( rep(tbx[1], 5), tbx[-NROW(tbx)+c(0:4)]);
	tb_4 <- c( rep(tbx[1], 4), tbx[-NROW(tbx)+c(0:3)]);
	tb_3 <- c( rep(tbx[1], 3), tbx[-NROW(tbx)+c(0:2)]);
	tb_2 <- c( rep(tbx[1], 2), tbx[-NROW(tbx)+c(0:1)]);
	tb_1 <- c( tbx[1], tbx[-NROW(tbx)]);
	tb0 <- tbx;
	tb1 <- c( tbx[-1], tbx[NROW(tbx)]);
	tb2 <- c( tbx[-c(1:2)], rep( tbx[NROW(tbx)],2) );
	tb3 <- c( tbx[-c(1:3)], rep( tbx[NROW(tbx)],3) );
	tb4 <- c( tbx[-c(1:4)], rep( tbx[NROW(tbx)],4) );
	tb5 <- c( tbx[-c(1:5)], rep( tbx[NROW(tbx)],5) );
	tb6 <- c( tbx[-c(1:6)], rep( tbx[NROW(tbx)],6) );
	tb7 <- c( tbx[-c(1:7)], rep( tbx[NROW(tbx)],7) );
	tb8 <- c( tbx[-c(1:8)], rep( tbx[NROW(tbx)],8) );
	tb9 <- c( tbx[-c(1:9)], rep( tbx[NROW(tbx)],9) );
	tb10<- c( tbx[-c(1:10)],rep( tbx[NROW(tbx)],10) );

	b0 <- ( tb0==tb1 & tb0==tb2  & tb0==tb3  & tb0==tb4  & tb0==tb5  & tb0==tb6  & tb0==tb7  & tb0==tb8  & tb0==tb9) 
	b1 <- (tb_1==tb0 & tb0==tb1  & tb0==tb2  & tb0==tb3  & tb0==tb4  & tb0==tb5  & tb0==tb6  & tb0==tb7  & tb0==tb8) 
	b2 <- (tb_2==tb0 & tb_1==tb0 & tb0==tb1  & tb0==tb2  & tb0==tb3  & tb0==tb4  & tb0==tb5  & tb0==tb6  & tb0==tb7)  
	b3 <- (tb_3==tb0 & tb_2==tb0 & tb_1==tb0 & tb0==tb1  & tb0==tb2  & tb0==tb3  & tb0==tb4  & tb0==tb5  & tb0==tb6)  
	b4 <- (tb_4==tb0 & tb_3==tb0 & tb_2==tb0 & tb_1==tb0 & tb0==tb1  & tb0==tb2  & tb0==tb3  & tb0==tb4  & tb0==tb5)  
	b5 <- (tb_5==tb0 & tb_4==tb0 & tb_3==tb0 & tb_2==tb0 & tb_1==tb0 & tb0==tb1  & tb0==tb2  & tb0==tb3  & tb0==tb4)  
	b6 <- (tb_6==tb0 & tb_5==tb0 & tb_4==tb0 & tb_3==tb0 & tb_2==tb0 & tb_1==tb0 & tb0==tb1  & tb0==tb2  & tb0==tb3)  
	b7 <- (tb_7==tb0 & tb_6==tb0 & tb_5==tb0 & tb_4==tb0 & tb_3==tb0 & tb_2==tb0 & tb_1==tb0 & tb0==tb1  & tb0==tb2)  
	b8 <- (tb_8==tb0 & tb_7==tb0 & tb_6==tb0 & tb_5==tb0 & tb_4==tb0 & tb_3==tb0 & tb_2==tb0 & tb_1==tb0 & tb0==tb1)  
	b9 <- (tb_9==tb0 & tb_8==tb0 & tb_7==tb0 & tb_6==tb0 & tb_5==tb0 & tb_4==tb0 & tb_3==tb0 & tb_2==tb0 & tb_1==tb0)  

	idx.rem <- which( b0+b1+b2+b3+b4+b5+b6+b7+b8+b9 >0 & tb0 < 3.0);
	if (NROW(idx.rem)>0) 
	   tb <- tb[ -idx.rem,,drop=F]

    #mapped <- get_histone_read(tb, file.unmap.bw)
    #tb <- tb[mapped==0,] ;
   
    file.tmp <- write.temp.bed(tb);

    tb <- read.table(pipe(paste("bedmap --echo --delim '\t' --indicator ", file.tmp, file.removed.reg )));
    if (sum(tb$V5==0)>0) 
        tb <- tb[tb$V5==0,,drop=F];
    tb <- tb[,c(1:4),drop=F] 

    unlink(file.tmp);
    file.tmp <- NULL;
    if(NROW(tb)>0)
        file.tmp <- write.temp.bed(tb, compress=TRUE);

    gc();
	
	return( file.tmp );
}

comvert_bed_to_bigwig<-function( input.path, prefix, output.path, rate=1 )
{
	file.beds <- list.files(path=input.path, pattern=glob2rx(paste0(prefix, "*.bed.gz")), full.names=T)

	file.tmp.beds <- unlist( mclapply(file.beds, function(file.bed)
	{
	   file.tmp <- remove_null_reg(file.bed, rate);
cat(file.bed, "=>", file.tmp, "\n");
	   return(file.tmp);
	}, mc.cores=4));

	file.tmp <- tempfile(fileext=".bed");
	system( paste("zcat", paste(file.tmp.beds, collapse="  " ), " > ", file.tmp ) );
	system( paste("bedGraphToBigWig", file.tmp, file.chromo.info, paste0(output.path, "/", prefix, ".bw")));
	unlink(file.tmp);
	unlink(file.tmp.beds);
	
}


file.temp.black  <- write.temp.bed(read.table(file.blacklist)[,c(1:3)], compress=FALSE )

tb.unmap.bed <- read.table(file.unmap.bed)[,c(1:3)];
tb.unmap.bed <- tb.unmap.bed[ tb.unmap.bed[,3] - tb.unmap.bed[,2]>100,,drop=F ]
file.temp.unmap  <- write.temp.bed(tb.unmap.bed, compress=FALSE )
rm(tb.unmap.bed);

file.removed.reg <- write.temp.bed(read.table(pipe(paste("cat ", file.temp.black, file.temp.unmap , " | sort-bed - | bedtools merge -i - "))));

if(0)
{
#comvert_bed_to_bigwig("../pred-gm/beds/", "H3k27ac.S1.V3.GM", "../pred-gm/");
comvert_bed_to_bigwig("../pred-gm/beds/", "H3k27me3.S1.V3.GM", "../pred-gm/");
comvert_bed_to_bigwig("../pred-gm/beds/", "H3k36me3.S1.V3.GM", "../pred-gm/");
comvert_bed_to_bigwig("../pred-gm/beds/", "H3k4me1.S1.V2.GM",  "../pred-gm/");
comvert_bed_to_bigwig("../pred-gm/beds/", "H3k4me2.S1.V2.GM",  "../pred-gm/");
comvert_bed_to_bigwig("../pred-gm/beds/", "H3k4me3.S1.V3.GM",  "../pred-gm/");
comvert_bed_to_bigwig("../pred-gm/beds/", "H3k9ac.S1.V2.GM",   "../pred-gm/");
comvert_bed_to_bigwig("../pred-gm/beds/", "H3k9me3.S1.V2.GM",  "../pred-gm/");
comvert_bed_to_bigwig("../pred-gm/beds/", "H4k20me1.S1.V3.GM", "../pred-gm/");
}

if(0)
{
#comvert_bed_to_bigwig("../pred-k562/beds/", "H3k27ac.S1.V3.G1", "../pred-k562/");
comvert_bed_to_bigwig("../pred-k562/beds/", "H3k122ac.S1.V3.G1", "../pred-k562/");
comvert_bed_to_bigwig("../pred-k562/beds/", "H3k27me3.S1.V3.G1", "../pred-k562/");
comvert_bed_to_bigwig("../pred-k562/beds/", "H3k36me3.S1.V3.G1", "../pred-k562/");
comvert_bed_to_bigwig("../pred-k562/beds/", "H3k4me1.S1.V2.G1",  "../pred-k562/");
comvert_bed_to_bigwig("../pred-k562/beds/", "H3k4me2.S1.V2.G1",  "../pred-k562/");
comvert_bed_to_bigwig("../pred-k562/beds/", "H3k4me3.S1.V3.G1",  "../pred-k562/");
comvert_bed_to_bigwig("../pred-k562/beds/", "H3k9ac.S1.V2.G1",   "../pred-k562/");
comvert_bed_to_bigwig("../pred-k562/beds/", "H3k9me3.S1.V2.G1",  "../pred-k562/");
comvert_bed_to_bigwig("../pred-k562/beds/", "H4k20me1.S1.V3.G1", "../pred-k562/");
}

if(0)
{
comvert_bed_to_bigwig("../pred-hct/beds/", "H3k122ac.S1.V3.HCT", "../pred-hct/");
comvert_bed_to_bigwig("../pred-hct/beds/", "H3k27ac.S1.V3.HCT",  "../pred-hct/");
comvert_bed_to_bigwig("../pred-hct/beds/", "H3k27me3.S1.V3.HCT", "../pred-hct/");
comvert_bed_to_bigwig("../pred-hct/beds/", "H3k36me3.S1.V3.HCT", "../pred-hct/");
comvert_bed_to_bigwig("../pred-hct/beds/", "H3k4me1.S1.V2.HCT",  "../pred-hct/");
comvert_bed_to_bigwig("../pred-hct/beds/", "H3k4me2.S1.V2.HCT",  "../pred-hct/");
comvert_bed_to_bigwig("../pred-hct/beds/", "H3k4me3.S1.V3.HCT",  "../pred-hct/");
comvert_bed_to_bigwig("../pred-hct/beds/", "H3k9ac.S1.V2.HCT",   "../pred-hct/");
comvert_bed_to_bigwig("../pred-hct/beds/", "H3k9me3.S1.V2.HCT",  "../pred-hct/");
comvert_bed_to_bigwig("../pred-hct/beds/", "H4k20me1.S1.V3.HCT", "../pred-hct/");
}


if(0)
{ 
comvert_bed_to_bigwig("../pred-cd4/beds/", "H3k122ac.S1.V3.CD4", "../pred-cd4/");
comvert_bed_to_bigwig("../pred-cd4/beds/", "H3k27ac.S1.V3.CD4",  "../pred-cd4/");
comvert_bed_to_bigwig("../pred-cd4/beds/", "H3k27me3.S1.V3.CD4", "../pred-cd4/");
comvert_bed_to_bigwig("../pred-cd4/beds/", "H3k36me3.S1.V3.CD4", "../pred-cd4/");
comvert_bed_to_bigwig("../pred-cd4/beds/", "H3k4me1.S1.V2.CD4",  "../pred-cd4/");
comvert_bed_to_bigwig("../pred-cd4/beds/", "H3k4me2.S1.V2.CD4",  "../pred-cd4/");
comvert_bed_to_bigwig("../pred-cd4/beds/", "H3k4me3.S1.V3.CD4",  "../pred-cd4/");
comvert_bed_to_bigwig("../pred-cd4/beds/", "H3k9ac.S1.V2.CD4",   "../pred-cd4/");
comvert_bed_to_bigwig("../pred-cd4/beds/", "H3k9me3.S1.V2.CD4",  "../pred-cd4/");
comvert_bed_to_bigwig("../pred-cd4/beds/", "H4k20me1.S1.V3.CD4", "../pred-cd4/");
}


if(0)
{
#comvert_bed_to_bigwig("../pred-hela/beds/", "H3k122ac.S1.V3.HELA", "../pred-hela/");
comvert_bed_to_bigwig("../pred-hela/beds/", "H3k27ac.S1.V3.HELA",  "../pred-hela/");
comvert_bed_to_bigwig("../pred-hela/beds/", "H3k27me3.S1.V3.HELA", "../pred-hela/");
comvert_bed_to_bigwig("../pred-hela/beds/", "H3k36me3.S1.V3.HELA", "../pred-hela/");
comvert_bed_to_bigwig("../pred-hela/beds/", "H3k4me1.S1.V2.HELA",  "../pred-hela/");
comvert_bed_to_bigwig("../pred-hela/beds/", "H3k4me2.S1.V2.HELA",  "../pred-hela/");
comvert_bed_to_bigwig("../pred-hela/beds/", "H3k4me3.S1.V3.HELA",  "../pred-hela/");
comvert_bed_to_bigwig("../pred-hela/beds/", "H3k9ac.S1.V2.HELA",   "../pred-hela/");
comvert_bed_to_bigwig("../pred-hela/beds/", "H3k9me3.S1.V2.HELA",  "../pred-hela/");
comvert_bed_to_bigwig("../pred-hela/beds/", "H4k20me1.S1.V3.HELA", "../pred-hela/");
}


if(0)
{
comvert_bed_to_bigwig("/workdir/zw355/proj/prj15-histone/horse_don/pred_by_k562_model", "H3k27ac.S1.V3.G1", "/workdir/zw355/proj/prj15-histone/horse_don/");
comvert_bed_to_bigwig("/workdir/zw355/proj/prj15-histone/horse_don/pred_by_k562_model", "H3k27me3.S1.V3.G1", "/workdir/zw355/proj/prj15-histone/horse_don/");
comvert_bed_to_bigwig("/workdir/zw355/proj/prj15-histone/horse_don/pred_by_k562_model", "H3k4me1.S1.V2.G1", "/workdir/zw355/proj/prj15-histone/horse_don/");
comvert_bed_to_bigwig("/workdir/zw355/proj/prj15-histone/horse_don/pred_by_k562_model", "H3k4me3.S1.V3.G1", "/workdir/zw355/proj/prj15-histone/horse_don/");
}


if(1)
{
#comvert_bed_to_bigwig("../pred-alex/", "Alex.H3K27ac.S1", "../pred-alex/");
#comvert_bed_to_bigwig("../pred-alex/", "Alex.H3k4me1.S1", "../pred-alex/");
#comvert_bed_to_bigwig("../pred-alex/", "Alex.H3k4me2.S1", "../pred-alex/");
#comvert_bed_to_bigwig("../pred-alex/", "Alex.H3k4me3.S1", "../pred-alex/");

#comvert_bed_to_bigwig("../pred-alex/", "Alex.Mnase.H3k4me1.S1", "../pred-alex/");
#comvert_bed_to_bigwig("../pred-alex/", "Alex.Mnase.H3k4me2.S1", "../pred-alex/");
#comvert_bed_to_bigwig("../pred-alex/", "Alex.Mnase.H3k4me3.S1", "../pred-alex/");
comvert_bed_to_bigwig("../pred-alex/", "Alex2.Mnase.H3k27ac.S1", "../pred-alex/");
}

if(0)
{
comvert_bed_to_bigwig("../atac/", "ATAC.K562.S1", "../atac/");
comvert_bed_to_bigwig("../atac/", "ATAC.exlchr22.K562.S1", "../atac/");
comvert_bed_to_bigwig("../dnase/", "dNase.exlchr22.K562.S1", "../dnase/", rate=2000);
}

unlink(file.removed.reg)
unlink(file.temp.black)