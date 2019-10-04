svm_predict_multfile_parallel <- function( str.prefix, file.model.rdata,  file.bigwig.plus.list, file.bigwig.minus.list, bed.region, ncores=4, linear_scale=F, gpu.count=1, str.wd="." )
{
    cpu.fun<-function(gpu.idx)
    {
        source("hist.svm.com.R")
        source("hist.svm.pred.R")
        source("hist.svm.main.R")
        
        library(dREG);
        library(bigWig);
        library(Rgtsvm);

        selectGPUdevice(gpu.idx);
        setwd(str.wd)

        model <- svm_load_model(file.model.rdata);   
        
        files.bigwig.plus <- file.bigwig.plus.list[seq(gpu.idx,NROW(file.bigwig.plus.list),gpu.count)]
        files.bigwig.minus <- file.bigwig.minus.list[seq(gpu.idx,NROW(file.bigwig.minus.list),gpu.count)]
                
        for(i in 1:NROW(files.bigwig.plus) )
        { 
            file.imputed.bed <- paste( str.prefix, ".", basename(files.bigwig.plus[i]), ".pred.bed.gz", sep="");
            if(file.exists(file.imputed.bed)) next;
cat("predicting =>", file.imputed.bed, "\n");

            task.step <- 50*1000;
            task.blks <- unique(c(seq(1, NROW(bed.region), task.step), NROW(bed.region)));
            ret.list <- list();
            for(k in 1:(NROW(task.blks)-1))
            {
                idx.start <- task.blks [ k ] 
                idx.end   <- task.blks [ k + 1 ] - 1;
                cat(task.blks [ k ], task.blks [ k + 1 ] - 1, "\n");

                ret.list[[k]] = run_svm_pred( model$gdm, bed.region[idx.start:idx.end, ], files.bigwig.plus[i], files.bigwig.minus[i], model$model, NULL, ncores=ncores, linear_scale=linear_scale) 
                gc(reset=TRUE);
            }

            pred.bed <- do.call("rbind", ret.list)
            write.bed(pred.bed, file=file.imputed.bed, compress=TRUE);
            rm(pred.bed);
            gc(reset=TRUE);
        }  

        Rgtsvm::predict.unload( model$model );
        
        rm(model);
        gc();
        
        return(1)
    }

    library(data.table)
    bed.region <- as.data.frame( rbindlist(mclapply(1:NROW(bed.region), function(i){ return(data.frame(bed.region[i,1], seq(bed.region[i,2], bed.region[i,3], 10)))}, mc.cores=16)) );
    bed.region$V3 <- bed.region[,2]+1;
    colnames(bed.region)<-c("chr", "start", "stop");
    cat("Bed.reegion=", NROW(bed.region), "\n")
    
    if(gpu.count>1)
    {
        library(snowfall);
        sfInit(parallel = TRUE, cpus = gpu.count, type = "SOCK" )
        sfExport("str.wd", "str.prefix", "file.model.rdata", "file.bigwig.plus.list", "file.bigwig.minus.list", "bed.region", "ncores", "gpu.count", "linear_scale" );

        fun <- as.function(cpu.fun);
        environment(fun)<-globalenv();

        do.call("rbindlist", sfLapply( c(1:gpu.count)-1, fun));
        sfStop();
    }
    else
       cpu.fun( gpu.count-1 );

    return(0);
}
