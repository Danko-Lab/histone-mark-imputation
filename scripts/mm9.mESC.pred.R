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


file.mm9.mESC.proseq.minus <- "/fs/cbsudanko/storage/data/mm9/esc/groseq/GSE48895_V6.5_untreated_Minus.bw"
file.mm9.mESC.proseq.plus <- "/fs/cbsudanko/storage/data/mm9/esc/groseq/GSE48895_V6.5_untreated_Plus.bw"

file.mm9.mESC.H3k4me3 <-  "../pred-mm9-mESC/mm9.GSM307618_ES.H3K4me3.bw"
file.mm9.mESC.H3k27me3 <- "../pred-mm9-mESC/mm9.GSM307619_ES.H3K27me3.bw"

path.histone <- NA;
path.proseq <- NA

library(Rgtsvm);
selectGPUdevice(0);

#pred <- svm_predict_chr_parallel( "MM9.mEsc.H3k4me3.S1", file.mm9.mESC.proseq.plus, file.mm9.mESC.proseq.minus, file.H3k4me3.model, chr="chr1", bigwig_compare=file.mm9.mESC.H3k4me3, ncores=5, linear_scale=F, gpu.idx=1)
#pred <- svm_predict_chr_parallel( "MM9.mEsc.H3k27me3.S1", file.mm9.mESC.proseq.plus, file.mm9.mESC.proseq.minus, file.H3k27me3.model, chr="chr1", bigwig_compare=file.mm9.mESC.H3k27me3, ncores=5, linear_scale=F, gpu.idx=1)
#pred <- svm_predict_chr_parallel( "MM9.mEsc.H3k4me3.S1", file.mm9.mESC.proseq.plus, file.mm9.mESC.proseq.minus, file.H3k4me3.model, chr="chr2", bigwig_compare=file.mm9.mESC.H3k4me3, ncores=5, linear_scale=F, gpu.idx=1)
#pred <- svm_predict_chr_parallel( "MM9.mEsc.H3k27me3.S1", file.mm9.mESC.proseq.plus, file.mm9.mESC.proseq.minus, file.H3k27me3.model, chr="chr2", bigwig_compare=file.mm9.mESC.H3k27me3, ncores=5, linear_scale=F, gpu.idx=1)
gc(reset=TRUE);
