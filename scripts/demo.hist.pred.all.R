source("hist.svm.com.R")
source("hist.svm.pred.R")
source("hist.svm.main.R")

args = commandArgs(trailingOnly=TRUE)
file.model <- args[1];
file.bw.plus <- args[2];
file.bw.minus <- args[3];
str.prefix <- args[4];
str.chr <- args[5];
file.histone <- if(args[6]=="NA") NULL else args[6];
gpu.count <- as.numeric(args[7])

cat("MODEL file = ", file.model, "\n");
cat("PLUS bigwig = ", file.bw.plus, "\n");
cat("MINUS bigwig = ", file.bw.minus, "\n");
cat("Prefix = ", str.prefix, "\n");
cat("histone file = ", file.histone, "\n");
cat("gpu count = ", gpu.count, "\n");

if(!file.exists(file.model)) stop("No model file!\n");
if(!file.exists(file.bw.plus)) stop("No plus file!\n");
if(!file.exists(file.bw.minus)) stop("No minus file!\n");
if(!is.null(file.histone) && !file.exists(file.histone)) stop("No histone file!\n");

pred <- svm_predict_all(str.prefix, file.bw.plus, file.bw.minus, file.model, bigwig_compare=file.histone, linear_scale=F, ncores=8, gpu.idx=unique(c(1:gpu.count)-1))
pred
