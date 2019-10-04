source("hist.svm.com.R")
source("hist.svm.pred.R")
source("hist.svm.main.R")
source("hist.param.R");

#path.histone defined in "hist.param.R"
#path.proseq  defined in "hist.param.R"
#file.rdata.negative defined in "hist.param.R"
#file.rdata.positive defined in "hist.param.R"

file.mm9.mESC.proseq.minus <- "/fs/cbsudanko/storage/data/mm9/esc/groseq/GSE48895_V6.5_untreated_Minus.bw"
file.mm9.mESC.proseq.plus  <- "/fs/cbsudanko/storage/data/mm9/esc/groseq/GSE48895_V6.5_untreated_Plus.bw"

file.mm9.mESC.H3k4me3  <- "../pred-mm9-mESC/mm9.GSM307618_ES.H3K4me3.int.bw"
file.mm9.mESC.H3k27me3 <- "../pred-mm9-mESC/mm9.GSM307619_ES.H3K27me3.bw"

file.mm9.mESC.H3k4me3.peak  <- "../pred-mm9-mESC/mm9.GSM307618_ES.H3K4me3-dREG.peak"
file.mm9.mESC.H3k27me3.peak <- "../pred-mm9-mESC/mm9.GSM307619_ES.H3K27me3-dREG.peak"

file.H3k4me3.model  <- "../models/mESC.H3k4me3.S1.train.rdata";
file.H3k27me3.model <- "../models/mESC.H3k27me3.S1.train.rdata";

path.histone="";
path.proseq="";


selectGPUdevice(0);

#training H3k4me3 Model
if(1)
{
if(!file.exists(file.H3k4me3.model))
{
  tb <- read.table(file.mm9.mESC.H3k4me3.peak);
  tb[,2] <- tb[,2] #- 10000
  tb[,3] <- tb[,3] #+ 10000
  file.tmp.pos.bed <- write.temp.bed(tb, compress=FALSE);

  tb <- read.table(pipe(paste("bedtools complement -i ",  file.tmp.pos.bed, " -g /fs/cbsudanko/storage/data/mm9/mm9.chromInfo.sorted | sort-bed -")));
  file.tmp.neg.bed <- write.temp.bed(tb, compress=FALSE);

  model <- create_train_model( path.proseq, file.mm9.mESC.proseq.plus, file.mm9.mESC.proseq.minus, 
                file.tmp.neg.bed, file.mm9.mESC.H3k4me3.peak, 
                path.histone, file.mm9.mESC.H3k4me3, 
                file.mm9.mESC.H3k4me3.peak, ratio = 0.1, 
                samples=400000*5, strategy=1, 
                exclude=c("chr1", "chr2", "chrX", "chrY", "chrM"))
  model <- build_train_model( model )
  model <- svm_train_model(model, gdm, file.H3k4me3.model, ncores=8);
  save(model, file = file.H3k4me3.model);
}
 
pred1 <- svm_predict_chr_parallel( "mESC.H3k4me3.trainbyself.S1", 
                file.mm9.mESC.proseq.plus, 
                file.mm9.mESC.proseq.minus, 
                file.H3k4me3.model, 
                chr=c("chr1"), 
                bigwig_compare=file.mm9.mESC.H3k4me3, 
                ncores=5, linear_scale=F, gpu.idx=0  );
                
pred2 <- svm_predict_chr_parallel( "mESC.H3k4me3.trainbyself.S1", 
                file.mm9.mESC.proseq.plus, 
                file.mm9.mESC.proseq.minus, 
                file.H3k4me3.model, 
                chr=c("chr2"), 
                bigwig_compare=file.mm9.mESC.H3k4me3, 
                ncores=5, linear_scale=F, gpu.idx=0  );
                
}

#training H3k27me3 Model
if(1) {
if(!file.exists(file.H3k27me3.model))
{
  tb <- read.table(file.mm9.mESC.H3k27me3.peak);
  tb[,2] <- tb[,2] #- 10000
  tb[,3] <- tb[,3] #+ 10000
  file.tmp.pos.bed <- write.temp.bed(tb, compress=FALSE);

  tb <- read.table(pipe(paste("bedtools complement -i ",  file.tmp.pos.bed, " -g /fs/cbsudanko/storage/data/mm9/mm9.chromInfo.sorted | sort-bed -")));
  file.tmp.neg.bed <- write.temp.bed(tb, compress=FALSE);

  model <- create_train_model( path.proseq, file.mm9.mESC.proseq.plus, file.mm9.mESC.proseq.minus, 
                file.tmp.neg.bed, file.mm9.mESC.H3k27me3.peak, 
                path.histone, file.mm9.mESC.H3k27me3, 
                file.mm9.mESC.H3k27me3.peak, ratio = 0.1, 
                samples=400000*5, strategy=1, 
                exclude=c("chr1", "chr2", "chrX", "chrY", "chrM"))
  model <- build_train_model( model )
  model <- svm_train_model(model, gdm, file.H3k27me3.model, ncores=8);
  save(model, file = file.H3k27me3.model);
}
 
pred1 <- svm_predict_chr_parallel( "mESC.H3k27me3.trainbyself.S1", 
                file.mm9.mESC.proseq.plus, 
                file.mm9.mESC.proseq.minus, 
                file.H3k27me3.model, 
                chr=c("chr1"), 
                bigwig_compare=file.mm9.mESC.H3k27me3, 
                ncores=5, linear_scale=F, gpu.idx=0  );
                
pred2 <- svm_predict_chr_parallel( "mESC.H3k27me3.trainbyself.S1", 
                file.mm9.mESC.proseq.plus, 
                file.mm9.mESC.proseq.minus, 
                file.H3k27me3.model, 
                chr=c("chr2"), 
                bigwig_compare=file.mm9.mESC.H3k27me3, 
                ncores=5, linear_scale=F, gpu.idx=0  );
}                
