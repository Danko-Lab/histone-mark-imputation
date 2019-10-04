source("hist.svm.com.R")
source("hist.svm.pred.R")
source("hist.svm.main.R")

args = commandArgs(trailingOnly=TRUE)
file.model <- args[1];
file.bw.plus <- args[2];
file.bw.minus <- args[3];
str.prefix <- args[4];
file.histone <- args[5];

cat("MIDEL file = ", file.model, "\n");
cat("PLUS bigwig = ", file.bw.plus, "\n");
cat("MINUS bigwig = ", file.bw.minus, "\n");
cat("Prefix = ", str.prefix, "\n");
cat("histone file = ", file.histone, "\n");

if(!file.exists(file.model)) stop("No model file!\n");
if(!file.exists(file.bw.plus)) stop("No plus file!\n");
if(!file.exists(file.bw.minus)) stop("No minus file!\n");
if(!is.null(file.histone) && !file.exists(file.histone)) stop("No histone file!\n");


library(Rgtsvm);

selectGPUdevice(0);

#load model object
load( file.model );
t0 <- proc.time();
pred <- svm_predict( str.prefix, file.bw.plus, file.bw.minus, model$svm, model$gdm, chrs=NULL, file.bed=NULL, bigwig_compare=file.histone, ncores=15, linear_scale=F  )
t1 <- proc.time() - t0;
show(t1);
cat("Output File=", pred$file.bed, "\n");
show(pred$r.cors);
