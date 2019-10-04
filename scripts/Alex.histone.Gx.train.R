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


file.Alex.H3k4me1.bw <- "/local/workdir/agc94/ChIP_sonication_Nov_2018/bam_files/10069_11142_82835_HTNJHBGX7_me_1_CGATGT.bw"
file.Alex.H3k4me2.bw <- "/local/workdir/agc94/ChIP_sonication_Nov_2018/bam_files/10069_11142_82836_HTNJHBGX7_me_2_TTAGGC.bw"
file.Alex.H3k4me3.bw <- "/local/workdir/agc94/ChIP_sonication_Nov_2018/bam_files/10069_11142_82837_HTNJHBGX7_me_3_TGACCA.bw"
file.Alex.H3k27ac.bw <- "/local/workdir/agc94/ChIP_sonication_Nov_2018/bam_files/10069_11142_82838_HTNJHBGX7_H3k27ac_ACAGTG.bw"

file.Alex.H3k4me1.peak <- "/local/workdir/zw355/proj/prj15-histone/scripts/Alex/Alex.H3k4me1_c12.0_l200_g30_peaks.narrowPeak"
file.Alex.H3k4me2.peak <- "/local/workdir/zw355/proj/prj15-histone/scripts/Alex/Alex.H3k4me2_c12.0_l200_g30_peaks.narrowPeak"
file.Alex.H3k4me3.peak <- "/local/workdir/zw355/proj/prj15-histone/scripts/Alex/Alex.H3k4me3_c12.0_l200_g30_peaks.narrowPeak"
file.Alex.H3k27ac.peak <- "/local/workdir/zw355/proj/prj15-histone/scripts/Alex/Alex.H3k27ac_c12.0_l200_g30_peaks.narrowPeak"


if(0)
{
   file.bw.histone <- file.Alex.H3k27ac.bw;
   file.peak.histone <- file.Alex.H3k27ac.peak;
   file.rdata.model <- "../models/Alex.H3k27ac.Gx.S1.train.rdata";
   str.prefix <- "Alex.H3k27ac.Gx.S1"
}

if(0)
{
    file.bw.histone <- file.Alex.H3k4me1.bw;
    file.peak.histone <- file.Alex.H3k4me1.peak;
    file.rdata.model <- "../models/Alex.H3k4me1.Gx.S1.train.rdata";
    str.prefix <- "Alex.H3k4me1.Gx.S1"
}

if(0)
{
    file.bw.histone <- file.Alex.H3k4me2.bw;
    file.peak.histone <- file.Alex.H3k4me2.peak;
    file.rdata.model <- "../models/Alex.H3k4me2.Gx.S1.train.rdata";
    str.prefix <- "Alex.H3k4me2.Gx.S1"
}

if(0)
{
    file.bw.histone <- file.Alex.H3k4me3.bw;
    file.peak.histone <- file.Alex.H3k4me3.peak;
    file.rdata.model <- "../models/Alex.H3k4me3.Gx.S1.train.rdata";
    str.prefix <- "Alex.H3k4me3.Gx.S1"
}

file.bws.plus <- c(file.bw.G1[1], file.bw.G2[1],file.bw.G3[1], file.bw.G5[1], file.bw.G6[1]);
file.bws.minus <-c(file.bw.G1[2], file.bw.G2[2],file.bw.G3[2], file.bw.G5[2], file.bw.G6[2]);
path.histone="";

if(!file.exists(file.rdata.model))
{
  model <- create_train_model( path.proseq, 
     file.bws.plus, file.bws.minus, 
     file.rdata.negative, file.rdata.positive, 
     path.histone, file.bw.histone, file.peak.histone, 
     ratio = 0.1, samples=600000, strategy=1, 
     exclude=c("chr22", "chrY", "chrM") )
  model <- build_train_model( model )
  model <- svm_train_model(model, gdm, file.rdata.model, ncores=15);
  save(model, file=file.rdata.model);
  rm(model);
  gc(reset=TRUE);
}

pred <- svm_predict_chr_parallel( str.prefix, 
   file.bws.plus[1], file.bws.minus[1], 
   file.rdata.model, 
   chr="chr22", 
   bigwig_compare=file.bw.histone, 
   ncores=5, linear_scale=F, gpu.idx=1  )
