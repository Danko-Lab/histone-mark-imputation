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

file.CFP1.H3k4me3.peak  <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/ChIP_alignments/k4me3.merge.all.narrowPeak.bed.gz"

file.CFP1.wt.plus.bw <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/NRO.wt_QC_plus.bw"
file.CFP1.wt.minus.bw <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/NRO.wt_QC_minus.bw"
file.CFP1.C169A.plus.bw <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/NRO.C169A_QC_plus.bw"
file.CFP1.C169A.minus.bw <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/NRO.C169A_QC_minus.bw"
file.CFP1.ko.plus.bw <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/NRO.ko_QC_plus.bw"
file.CFP1.ko.minus.bw <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/NRO.ko_QC_minus.bw"
file.CFP1.wt_resuce.plus.bw <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/NRO.wt_rescue_QC_plus.bw"
file.CFP1.wt_resuce.minus.bw <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/NRO.wt_rescue_QC_minus.bw"

file.CFP1.wt.H3k4me3.bw1 <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/ChIP_alignments/k4me3.wt.1.bw"
file.CFP1.wt.H3k4me3.bw2 <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/ChIP_alignments/k4me3.wt.2.bw"

file.CFP1.wt.H3k4me3.bw  <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/ChIP_alignments/k4me3.wt.merge.bw"
file.CFP1.ko.H3k4me3.bw  <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/ChIP_alignments/k4me3.ko.1.bw"
file.CFP1.C169A.H3k4me3.bw  <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/ChIP_alignments/k4me3.C169A.merge.bw"
file.CFP1.wt_resuce.H3k4me3.bw  <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/ChIP_alignments/k4me3.wt_rescue.merge.bw"

file.CFP1.wt.H3k4me3.peak <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/ChIP_alignments/k4me3.wt.merge.peak";
file.dreg.bed.peak <- "/fs/cbsudanko/storage/data/mm10/esc/cfp1ko/NRO.wt_QC_dreg.peak.bed.gz";

file.bw.histone <- file.CFP1.wt.H3k4me3.bw;
file.peak.histone <- file.CFP1.wt.H3k4me3.peak;
path.histone="";
path.proseq="";

file.bw.plus <- file.CFP1.wt.plus.bw;
file.bw.minus <- file.CFP1.wt.minus.bw;
file.rdata.model <- "../models/CFP1.wt.H3k4me3.S1.train.rdata";


if(0)
{
    if(!file.exists(file.rdata.model))
    {
       tmp.bed <- tempfile(fileext=".bed")
       tb <- read.table(file.dreg.bed.peak);
       tb[,2] <- tb[,2] - 100000
       tb[,3] <- tb[,3] + 100000
       write.bed(tb, file=tmp.bed,compress=FALSE);

       tb <- read.table(pipe(paste("bedtools complement -i ",  tmp.bed, " -g /fs/cbsudanko/storage/data/mm10/mm10.chromInfo | sort-bed -")));
       tmp.neg.bed <- tempfile(fileext=".bed.gz")
       write.bed(tb, file=tmp.neg.bed,compress=TRUE);

       model <- create_train_model( path.proseq, 
          file.bw.plus, file.bw.minus, 
          tmp.neg.bed, file.dreg.bed.peak, 
          path.histone, 
          file.bw.histone, 
          file.peak.histone, 
          ratio = 0.1, 
          samples=600000*5, 
          strategy=1, 
          exclude=c("chr19", "chrY", "chrM"))

       model <- build_train_model( model )
       gc();

       model <- svm_train_model(model, gdm, file.rdata.model, ncores=15);
       save(model, file=file.rdata.model);
    }
    gc(reset=TRUE); 
    pred <- svm_predict_chr_parallel( "CFP1.wt.H3k4me3.S1", 
        file.bw.plus, file.bw.minus, 
        file.rdata.model, 
        chr="chr19", 
        bigwig_compare=file.bw.histone, 
        ncores=5, linear_scale=F, gpu.idx=1  )
}

predict_TSS_bed<-function(file.rdata.model, file.peak.bed, file.bw.histone, file.bw.plus, file.bw.minus, file.peak.pred.gz)
{
    tb <- read.table(file.peak.bed, stringsAsFactors=FALSE)
    tb <- tb[grep("_|chrM|chrY|chrX", tb[,1], invert=TRUE),]
    #tb <- tb[tb[,1]=="chr19",]
    tb [,2] <- tb [,2]-10000
    tb [,3] <- tb [,3]+10000
    file.tmp.bed  <- write.temp.bed(tb, compress=FALSE)

    bed.peak <- read.table(pipe(paste("bedtools merge -i", file.tmp.bed )))
    bed.pred <- as.data.frame(rbindlist( mclapply(1:NROW(bed.peak ), function(i)
    {
        df <- data.frame(chr=bed.peak [i,1], start=seq(round(bed.peak [i,2]/10 -1 )*10, round(bed.peak [i,3]/10 + 1)*10, 10));
        df$stop <- df$start+1;
        return(df);
    }, mc.cores=10)));
    gc(reset=TRUE);
   
    #bed.pred <- unique(bed.pred);

    svm_predict_bed2( bed.pred, file.rdata.model, file.bw.histone, file.bw.plus, file.bw.minus, file.peak.pred.gz, ncores=10)
    gc(reset=TRUE);
    
    return(file.peak.pred.gz);
}

#predict_TSS_bed(file.rdata.model, file.CFP1.H3k4me3.peak, file.CFP1.wt.H3k4me3.bw, file.CFP1.wt.plus.bw, file.CFP1.wt.minus.bw, "../pred-cfp1/CFP1.wt.H3k4me3.S1.peak.bed.gz")
#predict_TSS_bed(file.rdata.model, file.CFP1.H3k4me3.peak, file.CFP1.C169A.H3k4me3.bw, file.CFP1.C169A.plus.bw, file.CFP1.C169A.minus.bw, "../pred-cfp1/CFP1.C169A.H3k4me3.S1.peak.bed.gz")
#selectGPUdevice(0);predict_TSS_bed(file.rdata.model, file.CFP1.H3k4me3.peak, file.CFP1.ko.H3k4me3.bw, file.CFP1.ko.plus.bw, file.CFP1.ko.minus.bw, "../pred-cfp1/CFP1.ko.H3k4me3.S1.peak.bed.gz")
#selectGPUdevice(0);predict_TSS_bed(file.rdata.model, file.CFP1.H3k4me3.peak, file.CFP1.wt_resuce.H3k4me3.bw, file.CFP1.wt_resuce.plus.bw, file.CFP1.wt_resuce.minus.bw, "../pred-cfp1/CFP1.wt_resuce.H3k4me3.S1.peak.bed.gz")



