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


proseq_hist_parallel <- function( chr.train="chr1", chr.pred="chr22", ncores=6, linear_scale=F, gpu.count=1 )
{
    cpu.fun<-function(gpu.idx)
    {
        source("hist.svm.com.R")
        source("h2h.svm.main.R")
        source("hist.param.R")
        
        library(dREG);
        library(bigWig);
        library(Rgtsvm);

cat("gpu.index=", gpu.idx, "\n", file="outfile.txt",append=TRUE);

        selectGPUdevice( gpu.idx )
        
        for(task.idx in 1:NROW(file.df.hist))
        {
            if ( task.idx %% gpu.count != gpu.idx) next;

            cat("Predict PROSEQ ==>", file.df.hist[task.idx,1], "\n", file="outfile.txt",append=TRUE)
            file.rdata.model <- paste("h2h", "G1", "to", basename(file.df.hist[task.idx,1]), chr.train, "rdata", sep=".");
            if(!file.exists(file.rdata.model))
            {
                model <- create_h2h_model( file.bw.G1, file.df.hist[task.idx,2], ratio = 0.1, samples=500000*2, strategy=1, include=chr.train )
                model <- build_h2h_model( model)
                model <- h2h_train_model(model, gdm, file.rdata.model, ncores=ncores );
                cat("model cor=", model$svm$fitted.r2, " pred  cor=", model$pred.cor, "\n");
    
                save(model, file=file.rdata.model);
                rm(model)
                gc();
            }

            file.predict.bed <-  paste("h2h", "G1", "to", file.df.hist[task.idx,1], chr.pred, "pred.bed.gz", sep=".");
            if(!file.exists(file.predict.bed))
            {
                bed.pred <- h2h_predict_chr_parallel( file.rdata.model, file.bw.G1, file.df.hist[task.idx,2], chr=chr.pred, ncores=ncores, linear_scale=linear_scale )
                write.bed(bed.pred, file.predict.bed, compress=TRUE);
                rm(bed.pred)
                gc(reset=TRUE);
            }
        }

        return(1)
    }

    if(gpu.count>1)
    {
        library(snowfall);
        sfInit(parallel = TRUE, cpus = gpu.count, type = "SOCK" )
        sfExport("file.df.hist", "chr.train", "chr.pred", "ncores", "gpu.count", "gdm", "linear_scale" );

        fun <- as.function(cpu.fun);
        environment(fun)<-globalenv();

        sfLapply((1:gpu.count)-1, fun)
        sfStop();
    }
    else
       for(i in (1:gpu.count)-1) cpu.fun( i );

    return(0);
}

# get the train results (G1-->Histone);
#proseq_hist_parallel(gpu.count=1);

h2h_parallel <- function( chr.train="chr1", chr.pred="chr22", ncores=6, linear_scale=F, gpu.count=1 )
{
    cpu.fun<-function(gpu.idx)
    {
        source("hist.svm.com.R")
        source("h2h.svm.main.R")
        
        library(dREG);
        library(bigWig);
        library(Rgtsvm);

cat("gpu.index=", gpu.idx, "\n", file="outfile.txt",append=TRUE);

        selectGPUdevice( gpu.idx )
        
        for(task.idx in 1:NROW(df.comb))
        {
            if ( task.idx %% gpu.count != gpu.idx) next;

            cat("Predict", df.comb[task.idx,1], "==>", df.comb[task.idx,2], "\n", file="outfile.txt",append=TRUE)
            file.rdata.model <- paste("h2h", df.comb[task.idx,1], "to", df.comb[task.idx,2], chr.train, "rdata", sep=".");
            if(!file.exists(file.rdata.model))
            {
                model <- create_h2h_model( as.character(df.comb[task.idx,3]), as.character(df.comb[task.idx,4]), ratio = 0.1, samples=500000*2, strategy=1, include=chr.train )
                model <- build_h2h_model( model)
                model <- h2h_train_model(model, gdm, file.rdata.model, ncores=ncores );
                cat("model cor=", model$svm$fitted.r2, " pred  cor=", model$pred.cor, "\n");
    
                save(model, file=file.rdata.model);
                rm(model)
                gc();
            }

            file.predict.bed <-  paste("h2h", df.comb[task.idx,1], "to", df.comb[task.idx,2], chr.pred, "pred.bed.gz", sep=".");
            if(!file.exists(file.predict.bed))
            {
                bed.pred <- h2h_predict_chr_parallel( file.rdata.model, as.character(df.comb[task.idx,3]), as.character(df.comb[task.idx,4]), chr=chr.pred, ncores=ncores, linear_scale=linear_scale )
                write.bed(bed.pred, file.predict.bed, compress=TRUE);
                rm(bed.pred)
                gc(reset=TRUE);
            }
        }

        return(1)
    }

    df.comb0 <- combn( 1:10, 2) 
    df.comb <- t(cbind(df.comb0, df.comb0[c(2,1),]))
    df.comb <- df.comb[order(df.comb[,1], df.comb[,2]),]

    df.comb <- data.frame(train_mark=as.character(file.df.hist[df.comb[,1],1]), pred_mark=as.character(file.df.hist[df.comb[,2],1]), train_hist=as.character(file.df.hist[df.comb[,1],2]), pred_hist=as.character(file.df.hist[df.comb[,2],2]), stringsAsFactors=FALSE)

    if(gpu.count>1)
    {
        library(snowfall);
        sfInit(parallel = TRUE, cpus = gpu.count, type = "SOCK" )
        sfExport("df.comb", "chr.train", "chr.pred", "ncores", "gpu.count", "gdm", "linear_scale" );

        fun <- as.function(cpu.fun);
        environment(fun)<-globalenv();

        sfLapply((1:gpu.count)-1, fun)
        sfStop();
    }
    else
       for(i in (1:gpu.count)-1) cpu.fun( i );

    return(0);
}

# get the train results (Histione --> Histone);
#h2h_parallel(gpu.count=4);


dKL <- function(v.pred, v.exp, win.size=1000, filter=FALSE, breaks=10000)
{
    win.blk <- unique(seq( 1, NROW(v.exp), win.size), NROW(v.exp));
    s.exp<- unlist(lapply(1:(NROW(win.blk)-1), function(i){ sum(v.exp[win.blk[i]:(win.blk[i+1]-1)]);}));
    s.pred<- unlist(lapply(1:(NROW(win.blk)-1), function(i){ sum(v.pred[win.blk[i]:(win.blk[i+1]-1)]);}));

    if(filter==TRUE)
    {
        h <- hist(c(s.exp,s.pred), breaks=breaks, plot=FALSE);
        hp <- hist(s.pred, breaks=h$breaks, plot=FALSE);

        s.pred.offset <- hp$breaks[which.max(hp$density)];
        s.pred <- s.pred - s.pred.offset
        s.pred[s.pred<0] <- 0;
    }

    #if breaks=100K, sum(h$density) = 1
    #if breaks=1K, sum(h$density) = 0.01
    h <- hist(c(s.exp,s.pred), breaks=breaks, plot=FALSE);
    hp <- hist(s.pred, breaks=h$breaks, plot=FALSE);
    hq <- hist(s.exp, breaks=h$breaks, plot=FALSE);
    hpd <- hp$density[h$density!=0.0]/sum(hp$density)
    hqd <- hq$density[h$density!=0.0]/sum(hq$density) 
    
    nozero <- .Machine$double.xmin
    r.KL <- sum(hpd*log((hpd+nozero)/(hqd+nozero)))
    r.cor <- cor( s.exp, s.pred )

    return( c( r.KL, r.cor)  );
}

    
calc_dKL<-function(df.pred, file.bw.histion, chr="chr22", filter=FALSE)
{
    df.comb0 <- combn( 1:10, 2) 
    df.comb <- t(cbind(df.comb0, df.comb0[c(2,1),]))
    df.comb <- df.comb[order(df.comb[,1], df.comb[,2]),]
    
    KLs <- win.cors <- cors <- c()
    df.comb <- data.frame(train_mark=as.character(file.df.hist[df.comb[,1],1]), pred_mark=as.character(file.df.hist[df.comb[,2],1]), train_hist=as.character(file.df.hist[df.comb[,1],2]), pred_hist=as.character(file.df.hist[df.comb[,2],2]), stringsAsFactors=FALSE)
    df.ret <- do.call("rbind", mclapply(1:NROW(df.comb), function(task.idx){
        file.predict.bed <-  paste("h2h", df.comb[task.idx,1], "to", df.comb[task.idx,2], chr, "pred.bed.gz", sep=".");
        tb <- read.table(file.predict.bed);
        r <- dKL(tb[,4], tb[,5], filter=filter);
        return(c(KL=r[1], win.cor=r[2], cor=cor(tb[,4], tb[,5]) ) )
    }, mc.cores=18));
    
    #show(df.comb);
    
    return(data.frame(df.comb,df.ret) );
}

if(0) 
{
    setwd("../h2h-model");

    KLs<- win.cors <- cors <- c();
    for(i in 1:NROW(file.df.hist))
     {
         file.predict.bed <-  paste("h2h", "G1", "to", file.df.hist[i,1], "chr22", "pred.bed.gz", sep=".");
        tb <- read.table(file.predict.bed);
        r.vec <- dKL(tb[,4], tb[,5], filter=F)
        KLs <- c(KLs, r.vec[1] )
        win.cors <- c(win.cors, r.vec[2] )
        cors <- c(cors, cor(tb[,4], tb[,5]));
    }
    r1 <- data.frame(train_mark="G1", pred_mark=file.df.hist[,1], train_hist="G1", pred_hist=file.df.hist[,2], KL=KLs, win.cor=win.cors, cor=cors);
    r0 <- calc_dKL(chr="chr22", filter=F);
    r <- rbind(r0, r1);

    rKL <- matrix(0, ncol=NROW(unique(r$pred_mark)), nrow=NROW(unique(r$train_mark)))
    rownames(rKL) <- unique(r$train_mark)
    colnames(rKL) <- unique(r$pred_mark)
    for(i in 1:NROW(r))
       rKL[r[i, "train_mark"], r[i, "pred_mark"]] <- r[i, "KL"]

    rWinCor <- matrix(1, ncol=NROW(unique(r$pred_mark)), nrow=NROW(unique(r$train_mark)))
    rownames(rWinCor) <- unique(r$train_mark)
    colnames(rWinCor) <- unique(r$pred_mark)
    for(i in 1:NROW(r))
       rWinCor[r[i, "train_mark"], r[i, "pred_mark"]] <- r[i, "win.cor"]

    rCor <- matrix(1, ncol=NROW(unique(r$pred_mark)), nrow=NROW(unique(r$train_mark)))
    rownames(rCor) <- unique(r$train_mark)
    colnames(rCor) <- unique(r$pred_mark)
    for(i in 1:NROW(r))
       rCor[r[i, "train_mark"], r[i, "pred_mark"]] <- r[i, "cor"]

    #save(r, rWinCor, rCor, rKL, file="DKL-h2h.train.filter.1k.10k.rdata");
    save(r, rWinCor, rCor, rKL, file="DKL-h2h.train.nofilter.1k.10k.rdata");

    library(gplots);
    library(corrplot);

    pdf("DKL-h2h.rkl.pdf")
    hM <- format(round(rKL, 1))
    heatmap.2(rKL, Rowv=F, Colv=T, symm=FALSE, col=rev(heat.colors(16)), trace="none",cellnote=hM, cexRow=0.8, cexCol=0.8, notecex=0.8, notecol="black" )
    dev.off()

    pdf("DKL-h2h.rWinCor.pdf")
    corrplot(as.matrix(rWinCor), number.cex = 0.7, addCoef.col = "black", tl.col = "black", tl.srt = 90 )
    dev.off()

    rCor["H3k4me2", "H3k27me3"] <-0.05
    pdf("DKL-h2h.rCor.pdf")
    corrplot(as.matrix(rCor), number.cex = 0.7, addCoef.col = "black", tl.col = "black", tl.srt = 90 ) #method = "square", order = "hclust", 
    dev.off()

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

calc_dist<-function(chr="chr22", filter=FALSE, fun_dist)
{
    df.comb0 <- combn( 1:10, 2) 
    df.comb <- t(cbind(df.comb0, df.comb0[c(2,1),]))
    df.comb <- df.comb[order(df.comb[,1], df.comb[,2]),]
    
    df.comb <- data.frame(train_mark=as.character(file.df.hist[df.comb[,1],1]), pred_mark=as.character(file.df.hist[df.comb[,2],1]), train_hist=as.character(file.df.hist[df.comb[,1],2]), pred_hist=as.character(file.df.hist[df.comb[,2],2]), stringsAsFactors=FALSE)
    df.ret <- do.call("rbind", mclapply(1:NROW(df.comb), function(task.idx){
        file.predict.bed <-  paste("h2h", df.comb[task.idx,1], "to", df.comb[task.idx,2], chr, "pred.bed.gz", sep=".");
        tb <- read.table(file.predict.bed);
        r <- fun_dist( tb[,5], tb[,4]);
        return(c(Dist=r, cor=cor(tb[,4], tb[,5]) ) )
    }, mc.cores=18));
    
    #show(df.comb);
    
    return(data.frame(df.comb,df.ret) );
}


if(0) 
{ 
    setwd("../h2h-model");
    
    fun.dist <- fun_L1;
    file.rdata <- "summary.h2h.train.L1.rdata"
    file.pdf <- "summary.h2h.train.L1.pdf"

    fun.dist <- fun_RMSE;
    file.rdata <- "summary.h2h.train.RMSE.rdata"
    file.pdf <- "summary.h2h.train.RMSE.pdf"
    
    Dists<- cors <- c();
    for(i in 1:NROW(file.df.hist))
     {
        file.predict.bed <-  paste("h2h", "G1", "to", file.df.hist[i,1], "chr22", "pred.bed.gz", sep=".");
        tb <- read.table(file.predict.bed);
        r.dist <- fun.dist(tb[,5], tb[,4] );
        Dists <- c(Dists, r.dist )
        cors <- c(cors, cor(tb[,4], tb[,5]));
    }
    r1 <- data.frame(train_mark="G1", pred_mark=file.df.hist[,1], train_hist="G1", pred_hist=file.df.hist[,2], Dist=Dists, cor=cors);
    r0 <- calc_dist(chr="chr22", filter=F, fun.dist);
    r <- rbind(r0, r1);

    rDist <- matrix(0, ncol=NROW(unique(r$pred_mark)), nrow=NROW(unique(r$train_mark)))
    rownames(rDist) <- unique(r$train_mark)
    colnames(rDist) <- unique(r$pred_mark)
    for(i in 1:NROW(r))
       rDist[r[i, "train_mark"], r[i, "pred_mark"]] <- r[i, "Dist"]

    rCor <- matrix(1, ncol=NROW(unique(r$pred_mark)), nrow=NROW(unique(r$train_mark)))
    rownames(rCor) <- unique(r$train_mark)
    colnames(rCor) <- unique(r$pred_mark)
    for(i in 1:NROW(r))
       rCor[r[i, "train_mark"], r[i, "pred_mark"]] <- r[i, "cor"]

    save(r, rCor, rDist, file=file.rdata);

    library(gplots);
    library(corrplot);

    pdf(file.pdf)
    hM <- format(round(rDist, 2))
    heatmap.2(rDist, Rowv=F, Colv=T, symm=FALSE, col=rev(heat.colors(16)), trace="none",cellnote=hM, cexRow=0.8, cexCol=0.8, notecex=0.8, notecol="black" )
    dev.off()
}






