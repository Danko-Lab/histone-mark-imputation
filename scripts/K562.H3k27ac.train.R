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

file.bws.plus <- c(file.bw.G1[1], file.bw.G2[1],file.bw.G3[1], file.bw.G5[1], file.bw.G6[1]);
file.bws.minus <-c(file.bw.G1[2], file.bw.G2[2],file.bw.G3[2], file.bw.G5[2], file.bw.G6[2]);

file.bw.histone <- file.H3k27ac.bw;
file.peak.histone <- file.H3k27ac.peak;

#model <- create_train_model( path.proseq, file.bws.plus, file.bws.minus, file.rdata.negative, file.rdata.positive, path.histone, file.bw.histone, file.peak.histone, ratio = 0.1, samples=600000, strategy=1, exclude=c("chr21", "chrY", "chrM"))
#model <- build_train_model( model )
#model <- svm_train_model(model, gdm, "../models/H3K27ac.S1.exc21.train.rdata", ncores=15);
#save(model, file="../models/H3K27ac.S1.exc21.train.rdata");

model <- create_train_model( path.proseq, file.bws.plus, file.bws.minus, file.rdata.negative, file.rdata.positive, path.histone, file.bw.histone, file.peak.histone, ratio = 0.1, samples=600000,strategy=1)
model <- build_train_model( model )
model <- svm_train_model(model, gdm, "../models/H3K27ac.S1.train.rdata", ncores=15);
save(model, file="../models/H3K27ac.S1.train.rdata");

