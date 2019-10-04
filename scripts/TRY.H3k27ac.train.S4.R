source("hist.svm.com.R")
source("hist.svm.pred.R")
source("hist.svm.main.R")
source("hist.param.R");

#args = commandArgs(trailingOnly=TRUE)
#file.model <- args[1];
#file.bw.plus <- args[2];
#file.bw.minus <- args[3];
#str.FID <- args[4];
#str.chr <- args[5];
#file.bed <- args[6];

#path.histone defined in "hist.param.R"
#path.proseq  defined in "hist.param.R"
#file.rdata.negative defined in "hist.param.R"
#file.rdata.positive defined in "hist.param.R"


detect_spike_reg<-function(file.org.peak, file.org.bw)
{
   peak_max <- get_histone_peak( file.org.peak, file.org.bw );
   y_max <- quantile(peak_max[,4],0.98);
   peak_spike <- peak_max[peak_max[,4] > y_max,] ;

   bw.hist <- load.bigWig(file.org.bw);
   max_count <- unlist(lapply(1:NROW(peak_spike), function(i) {
   	  bedv <- step.bpQuery.bigWig( bw.hist, as.character(peak_spike[i,1]),peak_spike[i,2], peak_spike[i,3],1)
   	  bedv_max <- sort(abs(bedv[-NROW(bedv)] - bedv[-1]), decreasing=T);
   	  return( (bedv_max[1]+bedv_max[2])/max(bedv)  );
   }));
   
   unload.bigWig( bw.hist)
   
   return(data.frame(peak_spike, max_count));
}


file.temp.black  <- write.temp.bed(read.table(file.blacklist)[,c(1:3)], compress=FALSE )
tb.unmap.bed <- read.table(file.unmap.bed)[,c(1:3)];
tb.unmap.bed <- tb.unmap.bed[ tb.unmap.bed[,3] - tb.unmap.bed[,2]>100,,drop=F ]
file.temp.unmap  <- write.temp.bed(tb.unmap.bed, compress=FALSE )
rm(tb.unmap.bed);

bed_spike <- detect_spike_reg(file.k562.H3k27ac.peak, file.k562.H3k27ac.bw)
file.temp.spike  <- write.temp.bed( bed_spike[bed_spike[,5]>1,1:3], compress=FALSE )
file.removed.reg <<- write.temp.bed(read.table(pipe(paste("cat ", file.temp.black, file.temp.unmap, file.temp.spike, " | sort-bed - | bedtools merge -i - "))));

file.bws.plus <- c(file.bw.G1[1], file.bw.G2[1],file.bw.G3[1], file.bw.G5[1], file.bw.G6[1]);
file.bws.minus <-c(file.bw.G1[2], file.bw.G2[2],file.bw.G3[2], file.bw.G5[2], file.bw.G6[2]);

file.all.histone <- c();
#1	Broad	K562	H3k27ac
file.all.histone[1] <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig"
#0	SYDH	HCT116	H3k27ac	 
file.all.histone[2] <- "/workdir/zw355/proj/prj15-histone/SydhHistone/wgEncodeSydhHistoneHct116H3k27acUcdSig.bigWig"
file.all.histone[3] <- "/local/workdir/zw355/proj/prj15-histone/VerifyHistone/GSM646435_K562_H3K27ac_rep2.bigWig"
file.all.histone[4] <- "/local/workdir/zw355/proj/prj15-histone/VerifyHistone/GSM646434_K562_H3K27ac_rep1.bigWig"
file.all.histone[5] <- "/local/workdir/zw355/proj/prj15-histone/VerifyHistone/GSM1782704_K562_10M_ChIP-seq_H3K27AC_nan_nan_1_0_hg19.bigWig"
file.all.histone[6] <- "/local/workdir/zw355/proj/prj15-histone/VerifyHistone/GSM2054696_K562_H3K27ac.10bpres.bigWig"
file.all.histone[7] <- "/local/workdir/zw355/proj/prj15-histone/VerifyHistone/GSM2309710_H3K27ac_K562_shNT.bigWig"
file.all.histone[8] <- "/local/workdir/zw355/proj/prj15-histone/VerifyHistone/GSM2877103_ChIP-seq_K562_H3K27ac_rep1.bw"
file.all.histone[9] <- "/local/workdir/zw355/proj/prj15-histone/VerifyHistone/GSM2877104_ChIP-seq_K562_H3K27ac_rep2.bw"
file.all.histone[10] <- "/local/workdir/zw355/proj/prj15-histone/VerifyHistone/GSM646434_K562_H3K27ac_rep1.bigWig"
file.all.histone[11] <- "/local/workdir/zw355/proj/prj15-histone/VerifyHistone/GSM646435_K562_H3K27ac_rep2.bigWig"

library(Rgtsvm);
selectGPUdevice(1);

file.model.rdata <- "../models/H3K27ac.S4.w10.V1.train.rdata";
if(!file.exists(file.model.rdata))
{
   model <- create_train_model( "/", file.bws.plus, file.bws.minus, file.rdata.negative, file.rdata.positive, 
        "/", file.k562.H3k27ac.bw, file.k562.H3k27ac.peak, 
        bed.blacklist=file.removed.reg,
        file.all.histone=file.all.histone,
        ratio = 0.1, samples=600000, 
        strategy=4, window.size=10,
        exclude=c("chr22", "chrY", "chrM"))
        
   model <- build_train_model( model )
   model <- svm_train_model(model, gdm, file.model.rdata, ncores=5);
   save(model, file=file.model.rdata);
}
cat("Predicting chr.22\n");
pred <- svm_predict_chr_parallel( "H3K27ac.S4.w10.V1", file.bws.plus[1], file.bws.minus[1], file.model.rdata, chr="chr22", bigwig_compare=file.all.histone[1], ncores=5, linear_scale=F, gpu.idx=0  )



file.model.rdata <- "../models/H3K27ac.S4.w50.V1.train.rdata";
if(!file.exists(file.model.rdata))
{
   model <- create_train_model( "/", file.bws.plus, file.bws.minus, file.rdata.negative, file.rdata.positive, 
        "/", file.k562.H3k27ac.bw, file.k562.H3k27ac.peak, 
        bed.blacklist=file.removed.reg,
        file.all.histone=file.all.histone,
        ratio = 0.1, samples=600000, 
        strategy=4, window.size=50,
        exclude=c("chr22", "chrY", "chrM"))
        
   model <- build_train_model( model )
   model <- svm_train_model(model, gdm, file.model.rdata, ncores=5);
   save(model, file=file.model.rdata);
}
cat("Predicting chr.22\n");
pred <- svm_predict_chr_parallel( "H3K27ac.S4.w50.V1", file.bws.plus[1], file.bws.minus[1], file.model.rdata, chr="chr22", bigwig_compare=file.all.histone[1], ncores=5, linear_scale=F, gpu.idx=0  )

file.model.rdata <- "../models/H3K27ac.S4.w100.V1.train.rdata";
if(!file.exists(file.model.rdata))
{
   model <- create_train_model( "/", file.bws.plus, file.bws.minus, file.rdata.negative, file.rdata.positive, 
        "/", file.k562.H3k27ac.bw, file.k562.H3k27ac.peak, 
        bed.blacklist=file.removed.reg,
        file.all.histone=file.all.histone,
        ratio = 0.1, samples=600000, 
        strategy=4, window.size=50*2,
        exclude=c("chr22", "chrY", "chrM"))
        
   model <- build_train_model( model )
   model <- svm_train_model(model, gdm, file.model.rdata, ncores=5);
   save(model, file=file.model.rdata);
}
cat("Predicting chr.22\n");
pred <- svm_predict_chr_parallel( "H3K27ac.S4.w100.V1", file.bws.plus[1], file.bws.minus[1], file.model.rdata, chr="chr22", bigwig_compare=file.all.histone[1], ncores=5, linear_scale=F, gpu.idx=0  )

