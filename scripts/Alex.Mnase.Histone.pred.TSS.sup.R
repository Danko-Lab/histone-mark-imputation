source("hist.svm.com.R")
source("hist.svm.pred.R")
source("hist.svm.main.R")
source("hist.param.R");


file.TSS.bed <- "/local/workdir/zw355/proj/prj15-histone/figures/Alex.Mnase.Tss.sup.bed"
file.bws.plus <- c(file.bw.G1[1], file.bw.G2[1],file.bw.G3[1], file.bw.G5[1], file.bw.G6[1]);
file.bws.minus <-c(file.bw.G1[2], file.bw.G2[2],file.bw.G3[2], file.bw.G5[2], file.bw.G6[2]);
path.histone <- ""
path.proseq  <- ""
ncores <- 5


predict_TSS_bed<-function()
{
cat("Model:", file.rdata.model, "\n");
    bed.Tss <- read.table(file.TSS.bed)
cat("NROW(bed.Tss)", NROW(bed.Tss), "\n");

    bed.pred <- do.call("rbind", mclapply(1:NROW(bed.Tss), function(i) {
        df <- data.frame(chr=bed.Tss[i,1], start=seq(round(bed.Tss[i,2]/10 -1 )*10, round(bed.Tss[i,3]/10 + 1)*10, 10));
        df$stop <- df$start+1;
        return(df);
    }, mc.cores=8));

    bed.pred <- unique(bed.pred);

cat("NROW(bed.pred)", NROW(bed.pred), "\n");

    y_pred <- c();
    task.seq <- unique(c(seq(1, NROW(bed.pred)+1, 500000), NROW(bed.pred)+1))

    load(file.rdata.model);
    model$svm <- Rgtsvm::predict.load(model$svm);
    gc(reset=TRUE);

    for(i in 1:(NROW(task.seq)-1))
    {
       idx.start <- task.seq[i];
       idx.stop <- task.seq[i+1]-1;
       x1 <- extract_feature_matrix( bed.pred[idx.start:idx.stop,], 
              paste(model$src$path.proseq, model$src$file.bws.plus[1], sep="/"), 
              paste(model$src$path.proseq, model$src$file.bws.minus[1], sep="/"), 
              model$gdm, linear_scale=F, ncores=ncores );

       y_pred0 <- Rgtsvm::predict.run(model$svm, x1$mat);
       y_pred <- c(y_pred, y_pred0)
       gc(reset=TRUE);
    }
 
cat("Bed file:", file.TSS.pred.gz, "\n");
    write.bed( data.frame(bed.pred, y_pred ), file.TSS.pred.gz, compress=TRUE);

    Rgtsvm::predict.unload(model$svm);

    return(file.TSS.pred.gz);
}



## H3k36me3
if(0)
{
file.bw.histone  <- file.Alex.Mnase.H3k36me3.bw;
file.peak.histone<- file.Alex.Mnase.H3k36me3.peak;
file.rdata.model <- "../models/Alex.Mnase.H3k36me3.Gx.S1.train.rdata";
file.TSS.pred.gz <- "../pred-alex/Alex-Mnase-H3k36me3.pred.TSS.sup.bed.gz"
selectGPUdevice(1);
predict_TSS_bed();
gc(reset=TRUE);
}



## H3k79me3
if(1)
{
file.bw.histone  <- file.Alex.Mnase.H3k79me3.bw;
file.peak.histone<- file.Alex.Mnase.H3k79me3.peak;
file.rdata.model <- "../models/Alex.Mnase.H3k79me3.Gx.S1.train.rdata";
file.TSS.pred.gz <- "../pred-alex/Alex-Mnase-H3k79me3.pred.TSS.sup.bed.gz"
selectGPUdevice(0);
predict_TSS_bed();
gc(reset=TRUE);
}

