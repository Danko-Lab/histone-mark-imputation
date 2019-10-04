source("hist.svm.com.R")
source("hist.param.R");
source("h2h.svm.main.R")

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

gdm <- genomic_data_model(window_sizes= c(10, 50, 500, 5000), half_nWindows= c(10, 20, 20, 20) )


file.df.hist <- data.frame( rbind(
    c(mark="H3k122ac", file=file.k562.H3k122ac.bw),
    c(mark="H3k27ac",  file=file.k562.H3k27ac.bw ),
    c(mark="H3k27me3", file=file.k562.H3k27me3.bw),
    c(mark="H3k36me3", file=file.k562.H3k36me3.bw),
    c(mark="H3k4me1",  file=file.k562.H3k4me1.bw ),
    c(mark="H3k4me2",  file=file.k562.H3k4me2.bw ),
    c(mark="H3k4me3",  file=file.k562.H3k4me3.bw ),
    c(mark="H3k9ac",   file=file.k562.H3k9ac.bw  ),
    c(mark="H3k9me3",  file=file.k562.H3k9me3.bw ),
    c(mark="H4k20me1", file=file.k562.H4k20me1.bw)), stringsAsFactors=FALSE );


h2toh_parallel <- function( hist.train.set, chr.train="chr1", chr.pred="chr22", ncores=6, linear_scale=F, gpu.count=1 )
{
    gpu.fun<-function(gpu.idx)
    {
        source("hist.svm.com.R")
        source("h2h.svm.main.R")
        
        library(dREG);
        library(bigWig);
        library(Rgtsvm);

cat("gpu.index=", gpu.idx, "\n", file="outfile.txt",append=TRUE);

        #selectGPUdevice( gpu.idx )
        selectGPUdevice( 1 )
        
        for(task.idx in 1:NROW(file.df.predict))
        {
            if ( task.idx %% gpu.count != gpu.idx) next;
             
            trainset <- paste(file.df.train[,1],collapse="-");  
            cat("Train", trainset, "==>", file.df.predict[task.idx,1], "on", chr.train, "\n", file="outfile.txt",append=TRUE)
            file.rdata.model <- paste("h2h", trainset, "to", file.df.predict[task.idx,1], chr.train, "rdata", sep=".");
            if(!file.exists(file.rdata.model))
            {
                model <- create_h2h_model( file.df.train[,2], file.df.predict[task.idx,2], ratio = 0.1, samples=500000*2, strategy=1, include=chr.train )
                model <- build_h2h_model( model)
                model <- h2h_train_model(model, gdm, file.rdata.model, ncores=ncores );
                cat("model cor=", model$svm$fitted.r2, " pred  cor=", model$pred.cor, "\n");
    
                save(model, file=file.rdata.model);
                rm(model)
                gc();
            }

            file.predict.bed <-  paste("h2h", trainset, "to", file.df.predict[task.idx,1], chr.pred, "pred.bed.gz", sep=".");
            cat("Predict", trainset, "on", chr.pred, "==>", file.predict.bed, "\n", file="outfile.txt",append=TRUE)
            if(!file.exists(file.predict.bed))
            {
                bed.pred <- h2h_predict_chr_parallel( file.rdata.model, file.df.train[,2], file.df.predict[task.idx,2], chr=chr.pred, ncores=ncores, linear_scale=linear_scale )
                write.bed(bed.pred, file.predict.bed, compress=TRUE);
                rm(bed.pred)
                gc(reset=TRUE);
            }
        }

        gc(reset=TRUE);
        return(1)
    }

    idx.rem <- match(hist.train.set, file.df.hist[,1])
    file.df.train <- file.df.hist[idx.rem,]
    file.df.predict <- file.df.hist[-idx.rem,]

    if(gpu.count>1)
    {
        library(snowfall);
        sfInit(parallel = TRUE, cpus = gpu.count, type = "SOCK" )
        sfExport("file.df.predict", "file.df.train", "chr.train", "chr.pred", "ncores", "gpu.count", "gdm", "linear_scale" );

        fun <- as.function(gpu.fun);
        environment(fun)<-globalenv();

        sfLapply((1:gpu.count)-1, fun)
        sfStop();
    }
    else
       for(i in (1:gpu.count)-1) gpu.fun( i );
    
    gc(reset=TRUE); 
    return(0);
}

if(0)
{
# get the train results (Histione*Histone --> Histone);
#h2toh_parallel(c("H3k4me2", "H3k36me3"), gpu.count=2);
#h2toh_parallel(c("H3k4me2", "H3k36me3", "H3k27ac"), gpu.count=1);
#h2toh_parallel(c("H3k4me2", "H3k36me3", "H3k27ac", "H3k27me3"), gpu.count=1);
#h2toh_parallel(c("H3k4me2", "H3k36me3", "H3k27ac", "H3k27me3", "H3k9me3"), gpu.count=1);
}

library(lmtest);
fun_granger <- function(x, y)
{
    r <- grangertest(y, x, 1 );
    return(r$F[2]);
}

fun_L1 <- function(x, y)
{
    mean( abs((x- median(x))/sd(x) - (y- median(y))/sd(y) ) )
}

fun_RMSE <- function(x, y)
{
    sqrt(sum( ((x- median(x))/sd(x) - (y- median(y))/sd(y) )^2 )/NROW(x))
}

calc_dist<-function(file.predict.bed)
{
    tb <- read.table(file.predict.bed);
    return(c(fun.dist(tb[,5], tb[,4] ), cor(tb[,4], tb[,5])));
}

if(0) 
{ 
    setwd("../h2h-model");
    
    fun.dist <- fun_L1;
    file.rdata <- "summary.hx2h.train.L1.rdata"
    file.pdf <- "summary.hx2h.train.L1.pdf"

    #fun.dist <- fun_RMSE;
    #file.rdata <- "summary.h2h.train.RMSE.rdata"
    #file.pdf <- "summary.h2h.train.RMSE.pdf"

    df <- data.frame(rbind(
    c("H3k4me2-H3k36me3", "H3k122ac", "h2h.H3k4me2-H3k36me3.to.H3k122ac.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3", "H3k27me3", "h2h.H3k4me2-H3k36me3.to.H3k27me3.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3", "H3k4me1",  "h2h.H3k4me2-H3k36me3.to.H3k4me1.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3", "H3k4me3",  "h2h.H3k4me2-H3k36me3.to.H3k4me3.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3", "H3k9ac",   "h2h.H3k4me2-H3k36me3.to.H3k9ac.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3", "H3k9me3",  "h2h.H3k4me2-H3k36me3.to.H3k9me3.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3", "H4k20me1", "h2h.H3k4me2-H3k36me3.to.H4k20me1.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3", "H3k27ac",  "h2h.H3k4me2-H3k36me3.to.H3k27ac.chr22.pred.bed.gz"),
    
    c("H3k4me2-H3k36me3-H3k27ac", "H3k122ac",  "h2h.H3k4me2-H3k36me3-H3k27ac.to.H3k122ac.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3-H3k27ac", "H3k27me3",  "h2h.H3k4me2-H3k36me3-H3k27ac.to.H3k27me3.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3-H3k27ac", "H3k4me1",   "h2h.H3k4me2-H3k36me3-H3k27ac.to.H3k4me1.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3-H3k27ac", "H3k4me3",   "h2h.H3k4me2-H3k36me3-H3k27ac.to.H3k4me3.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3-H3k27ac", "H3k9ac",    "h2h.H3k4me2-H3k36me3-H3k27ac.to.H3k9ac.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3-H3k27ac", "H3k9me3",   "h2h.H3k4me2-H3k36me3-H3k27ac.to.H3k9me3.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3-H3k27ac", "H4k20me1",  "h2h.H3k4me2-H3k36me3-H3k27ac.to.H4k20me1.chr22.pred.bed.gz"),
    
    c("H3k4me2-H3k36me3-H3k27ac-H3k27me3", "H3k122ac",  "h2h.H3k4me2-H3k36me3-H3k27ac-H3k27me3.to.H3k122ac.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3-H3k27ac-H3k27me3", "H3k4me1",   "h2h.H3k4me2-H3k36me3-H3k27ac-H3k27me3.to.H3k4me1.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3-H3k27ac-H3k27me3", "H3k4me3",   "h2h.H3k4me2-H3k36me3-H3k27ac-H3k27me3.to.H3k4me3.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3-H3k27ac-H3k27me3", "H3k9ac",    "h2h.H3k4me2-H3k36me3-H3k27ac-H3k27me3.to.H3k9ac.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3-H3k27ac-H3k27me3", "H3k9me3",   "h2h.H3k4me2-H3k36me3-H3k27ac-H3k27me3.to.H3k9me3.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3-H3k27ac-H3k27me3", "H4k20me1",  "h2h.H3k4me2-H3k36me3-H3k27ac-H3k27me3.to.H4k20me1.chr22.pred.bed.gz"),
    
    c("H3k4me2-H3k36me3-H3k27ac-H3k27me3-H3k9me3", "H3k122ac",  "h2h.H3k4me2-H3k36me3-H3k27ac-H3k27me3-H3k9me3.to.H3k122ac.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3-H3k27ac-H3k27me3-H3k9me3", "H3k4me1",   "h2h.H3k4me2-H3k36me3-H3k27ac-H3k27me3-H3k9me3.to.H3k4me1.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3-H3k27ac-H3k27me3-H3k9me3", "H3k4me3",   "h2h.H3k4me2-H3k36me3-H3k27ac-H3k27me3-H3k9me3.to.H3k4me3.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3-H3k27ac-H3k27me3-H3k9me3", "H3k9ac",    "h2h.H3k4me2-H3k36me3-H3k27ac-H3k27me3-H3k9me3.to.H3k9ac.chr22.pred.bed.gz"),
    c("H3k4me2-H3k36me3-H3k27ac-H3k27me3-H3k9me3", "H4k20me1",  "h2h.H3k4me2-H3k36me3-H3k27ac-H3k27me3-H3k9me3.to.H4k20me1.chr22.pred.bed.gz")))

    r <- c();
    for(i in 1:NROW(df))
    {
        rD <- calc_dist(as.character(df[i,3]))
        r <- rbind(r, data.frame(train_mark=df[i,1], pred_mark=df[i,2], Dist=rD[1], cor=rD[2] ));
    }

    rDist <- matrix(0, ncol=NROW(unique(r$pred_mark)), nrow=NROW(unique(r$train_mark)))
    rownames(rDist) <- unique(r$train_mark)
    colnames(rDist) <- unique(r$pred_mark)
    for(i in 1:NROW(r))
       rDist[as.character(r[i, "train_mark"]), as.character(r[i, "pred_mark"])] <- r[i, "Dist"]

    rCor <- matrix(1, ncol=NROW(unique(r$pred_mark)), nrow=NROW(unique(r$train_mark)))
    rownames(rCor) <- unique(r$train_mark)
    colnames(rCor) <- unique(r$pred_mark)
    for(i in 1:NROW(r))
       rCor[as.character(r[i, "train_mark"]), as.character(r[i, "pred_mark"])] <- r[i, "cor"]

    save(r, rCor, rDist, file=file.rdata);
rownames(rDist)<-c("4me2+36me3", "++3k27ac", "++3k27me3", "++3k9me3"); 

    library(gplots);
    library(corrplot);

    pdf(file.pdf)
    hM <- format(round(rDist, 2))
    heatmap.2(rDist, Rowv=F, Colv=T, symm=FALSE, col=rev(heat.colors(16)), trace="none",cellnote=hM, cexRow=0.8, cexCol=0.8, notecex=0.8, notecol="black" )
    dev.off()
}


if(0)
{
    file.pdf <- "summary.hx2h.train.L1.pdf"

    load("summary.hx2h.train.L1.rdata")
    rDist4 <- rDist;
    load("summary.h2h.train.L1.rdata")
    rDist4 <- rbind(rDist4, c(rDist["G1","H3k122ac"], rDist["G1","H3k27me3"], rDist["G1","H3k4me1"], rDist["G1","H3k4me3"], rDist["G1","H3k9ac"], rDist["G1","H3k9me3"], rDist["G1","H4k20me1"], rDist["G1","H3k27ac"]))
    rownames(rDist4) <- c("4me2+36me3", "++3k27ac", "++3k27me3", "++3k9me3", "G1")
    rDist <- rDist4

    library(gplots);
    library(corrplot);

    pdf(file.pdf)
    hM <- format(round(rDist, 2))
    heatmap.2(rDist, Rowv=F, Colv=T, symm=FALSE, col=rev(heat.colors(16)), trace="none",cellnote=hM, cexRow=0.8, cexCol=0.8, notecex=0.8, notecol="black" )
    dev.off()
}