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


file.bw.histone <- file.k562.H3k27me3.bw;
file.peak.histone <- file.k562.H3k27me3.peak;

df.bed <- rbind(data.frame(chr="chr17", start=25795759, stop=48333103),
                data.frame(chr="chr20", start=29520574, stop=61938870),   
                data.frame(chr="chr18", start=19651418, stop=43559135),
                data.frame(chr="chr21", start=17734911, stop=49415365));

gdm <- genomic_data_model(window_sizes= c(10, 25, 50, 500, 5000), half_nWindows= c(5, 5, 15, 10, 10) )


file.rdata.model <- "../models/H3k27me3.S5.V2.train.rdata";

selectGPUdevice(1);

if(!file.exists(file.rdata.model))
{
  model <- create_train_model( "/", file.bws.plus, file.bws.minus, file.rdata.negative, file.rdata.positive, "/", file.bw.histone, file.peak.histone, ratio = 0.1, samples=600000, strategy=5, include.bed=df.bed, exclude=c("chr22", "chrY", "chrM", "chrUn")  )
  model <- build_train_model( model )
  model <- svm_train_model(model, gdm, file.rdata.model, ncores=15);
  save(model, file=file.rdata.model);
}

red <- svm_predict_chr_parallel( "K562.H3k27me3.S5.V2", file.bws.plus, file.bws.minus, file.rdata.model, chr="chr22", bigwig_compare=file.bw.histone, ncores=5, linear_scale=F, gpu.idx=1  )

