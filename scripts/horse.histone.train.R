source("hist.svm.com.R")
source("hist.svm.pred.R")
source("hist.svm.main.R")
source("horse.param.R");

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

file.bws.plus <- c(file.proseq.horse[1]);
file.bws.minus <-c(file.proseq.horse[2]);

library(Rgtsvm);
selectGPUdevice(1);

if(0)
{
   file.model.rdata <- "../models/horse.exchr22.H3k27ac.S1.train.rdata"

   model <- create_train_model( path.proseq, file.bws.plus, file.bws.minus, file.H3k27ac.negative, file.H3k27ac.positive, path.histone, file.horse.H3k27ac.bw, file.horse.H3k27ac.peak, ratio = 0.1, samples=600000*4, strategy=1, exclude=c("chr22", "chrY", "chrM", "chrUn") )
   save(model, gdm, file=file.model.rdata);
   model <- build_train_model( model )
   save(model, gdm, file=file.model.rdata);
   model <- svm_train_model(model, gdm, file.model.rdata, ncores=1, skip.svm=TRUE);
   save(model, gdm, file=file.model.rdata);
}

if(0)
{
   file.model.rdata <- "../models/horse.exchr22.H3k27me3.S1.train.rdata"

   model <- create_train_model( path.proseq, file.bws.plus, file.bws.minus, file.H3k27me3.negative, file.H3k27me3.positive, path.histone,  file.horse.H3k27me3.bw, file.horse.H3k27me3.peak, ratio = 0.1, samples=600000*4, strategy=1, exclude=c("chr22", "chrY", "chrM", "chrUn") )
   save(model, file=file.model.rdata);
   model <- build_train_model( model )
   save(model, gdm, file=file.model.rdata);
   model <- svm_train_model(model, gdm, file.model.rdata, ncores=15, skip.svm=TRUE);
   save(model, gdm, file=file.model.rdata);

}

if(0)
{
   file.model.rdata <- "../models/horse.exchr22.H3k4me1.S1.train.rdata"

   model <- create_train_model( path.proseq, file.bws.plus, file.bws.minus, file.H3k4me3.negative, file.H3k4me1.positive, path.histone, file.horse.H3k4me1.bw, file.horse.H3k4me1.peak, ratio = 0.1, samples=600000*4, strategy=1, exclude=c("chr22", "chrY", "chrM", "chrUn") )
   save(model, gdm, file=file.model.rdata);
   model <- build_train_model( model )
   save(model, gdm, file=file.model.rdata);
   model <- svm_train_model(model, gdm, file.model.rdata, ncores=15, skip.svm=TRUE);
   save(model, gdm, file=file.model.rdata);
}

if(0)
{
    
   file.model.rdata <- "../models/horse.exchr22.H3k4me3.S1.train.rdata"
    
   model <- create_train_model( path.proseq, file.bws.plus, file.bws.minus, file.H3k4me3.negative, file.H3k4me3.positive, path.histone, file.horse.H3k4me3.bw, file.horse.H3k4me3.peak, ratio = 0.1, samples=600000*4, strategy=1, exclude=c("chr22", "chrY", "chrM", "chrUn") )
   save(model, gdm, file=file.model.rdata);
   model <- build_train_model( model );
   save(model, gdm, file=file.model.rdata);
   model <- svm_train_model(model, gdm, file.model.rdata, ncores=15, skip.svm=TRUE);
   save(model, gdm, file=file.model.rdata);
}
