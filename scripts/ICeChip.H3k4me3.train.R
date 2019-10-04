source("hist.svm.com.R")
source("hist.svm.pred.R")
source("hist.svm.main.R")
source("hist.param.R");

#path.histone defined in "hist.param.R"
#path.proseq  defined in "hist.param.R"
#file.rdata.negative defined in "hist.param.R"
#file.rdata.positive defined in "hist.param.R"

file.H3k4me3.bw <- "../pred-iceChip/GSK562_AB-12209_H3K4me3_ICeChIP.hg19.bw"
file.H3k4me3.peak <- "../pred-iceChip/GSK562_AB-12209_H3K4me3_ICeChIP.hg19.macs2.peak.bed";
file.bws.plus <- c(file.bw.G1[1], file.bw.G2[1],file.bw.G3[1], file.bw.G5[1], file.bw.G6[1]);
file.bws.minus <-c(file.bw.G1[2], file.bw.G2[2],file.bw.G3[2], file.bw.G5[2], file.bw.G6[2]);

file.bw.histone <- file.H3k4me3.bw
file.peak.histone <- file.H3k4me3.peak;
path.histone="";
path.proseq="";

file.rdata.model <- "../models/IceChip.H3k4me3.S1.train.rdata";

selectGPUdevice(0);

if(!file.exists(file.rdata.model))
{
 model <- create_train_model( path.proseq, file.bws.plus, file.bws.minus, file.rdata.negative, file.rdata.positive, path.histone, file.bw.histone, file.peak.histone, ratio = 0.1, samples=600000, strategy=1, exclude=c("chr22", "chrY", "chrM"))
 model <- build_train_model( model )
 model <- svm_train_model(model, gdm, file.rdata.model, ncores=15);
 save(model, file=file.rdata.model);
}
 
pred <- svm_predict_chr_parallel( "IceChip.H3k4me3.S1", file.bws.plus, file.bws.minus, file.rdata.model, chr="chr22", bigwig_compare=file.bw.histone, ncores=5, linear_scale=F, gpu.idx=0  )
