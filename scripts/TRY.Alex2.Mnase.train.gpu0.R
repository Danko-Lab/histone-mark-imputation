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


do.gpu<-function(gpu.idx)
{
    path.histone="";
    path.proseq="";

    file.bws.plus <- file.Alex.proseq.plus.m.bw
    file.bws.minus <- file.Alex.proseq.minus.m.bw
    selectGPUdevice(gpu.idx);

    if(!file.exists(file.rdata.model))
    {
     model <- create_train_model( path.proseq, 
		    file.bws.plus, file.bws.minus, 
		    file.rdata.negative, file.rdata.positive, 
		    path.histone, file.bw.histone, file.peak.histone, 
		    ratio = 0.1, samples=600000*5, strategy=1, 
		    exclude=c("chr22", "chrY", "chrX","chrM"))

     model <- build_train_model( model )
     gc();
     model <- svm_train_model(model, gdm, file.rdata.model, ncores=5);
     gc();
     save(model, file=file.rdata.model);
     gc(reset=TRUE);
    }

    pred <- svm_predict_chr_parallel( str.prefix, 
		    file.bws.plus, file.bws.minus, 
		    file.rdata.model, chr="chr22", 
		    bigwig_compare=file.bw.histone, 
		    ncores=5, linear_scale=F, gpu.idx=gpu.idx  )
     gc(reset=TRUE);
		 
}

##H3k4me1
if(1)
{
file.bw.histone <- file.Alex.Mnase.H3k4me1.bw;
file.peak.histone <- file.Alex.Mnase.H3k4me1.peak;
file.rdata.model <- file.Alex.Mnase.H3k4me1.AL.model;
str.prefix <- "Alex.Mnase.H3k4me1.AL.S1"
}

do.gpu(0);

##H3k4me2
if(1)
{
file.bw.histone <- file.Alex.Mnase.H3k4me2.bw;
file.peak.histone <- file.Alex.Mnase.H3k4me2.peak;
file.rdata.model <- file.Alex.Mnase.H3k4me2.AL.model;
str.prefix <- "Alex.Mnase.H3k4me2.AL.S1"
}

do.gpu(0);

##H3k4me3
if(1)
{
file.bw.histone <- file.Alex.Mnase.H3k4me3.bw;
file.peak.histone <- file.Alex.Mnase.H3k4me3.peak;
file.rdata.model <- file.Alex.Mnase.H3k4me3.AL.model;
str.prefix <- "Alex.Mnase.H3k4me3.AL.S1"
}

do.gpu(0);


		
