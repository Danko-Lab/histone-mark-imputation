source("hist.svm.com.R")
source("hist.svm.pred.R")
source("hist.svm.main.R")

args = commandArgs(trailingOnly=TRUE)
file.model <- args[1];
file.bw.plus <- args[2];
file.bw.minus <- args[3];
str.prefix <- args[4];
str.chr <- args[5];
file.histone <- if(is.na(args[6])) NULL else args[6];
gpu.count <- if(is.na(args[7])) 1 else as.numeric(args[7])

if (file.histone=="NA") file.histone <- NULL
 
cat("MODEL file = ", file.model, "\n");
cat("PLUS bigwig = ", file.bw.plus, "\n");
cat("MINUS bigwig = ", file.bw.minus, "\n");
cat("Prefix = ", str.prefix, "\n");
cat("CHR = ", str.chr, "\n");
cat("histone file = ", file.histone, "\n");
cat("gpu count = ", gpu.count, "\n");

if(!file.exists(file.model)) stop("No model file!\n");
if(!file.exists(file.bw.plus)) stop("No plus file!\n");
if(!file.exists(file.bw.minus)) stop("No minus file!\n");
if(!is.null(file.histone) && !file.exists(file.histone)) stop("No histone file!\n");

library(Rgtsvm)
selectGPUdevice(0);

if(str.chr=="all" || str.chr=="ALL")
   str.chr <- NULL;
   
svm_predict_chr_parallel( str.prefix, file.bw.plus, file.bw.minus, file.model, chr=str.chr, bigwig_compare=file.histone, ncores=5, linear_scale=F, gpu.idx=unique(c(1:gpu.count)-1)  )

