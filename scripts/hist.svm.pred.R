options(scipen =99)
options(error=traceback);

bedtoBigWig<-function(file.bw.plus, file.bw.minus, file.beds, file.bigWig)
{
    temp.bed <- do.call("rbind", lapply(file.beds, function(x){return(read.table(x, header=F))}));

    temp.file <- tempfile(fileext=".bed");
    write.table(temp.bed, file=temp.file, row.names=F, col.names=F, quote=F, sep="\t");

    chromInfo <- get.chromosome.info( file.bw.plus, file.bw.minus );
    tobigwig(filename=file.bigWig, temp.bg = temp.file, chromInfo= chromInfo)

    unlink(temp.file);

    return(file.bigWig);
}

compare_correlation<-function( file.bed, bigwig_compare )
{
   r.cor1 <- r.cor2 <- NA;
   if( file.exists(file.bed) && !is.null(bigwig_compare) )
   {
       ret <- read.table(file.bed);
       y.label <-  get_histone_read( ret[,c(1:3)], bigwig_compare);
       if(!is.null(y.label))
       {
           r.cor1 <- cor(ret[,4], y.label, method="pearson");
           r.cor2 <- cor(ret[,4], y.label, method="spearman");
           cat("COR=", r.cor1, r.cor2, "\n" )
       }
    }
    return( c(r.cor1,  r.cor2) );
}

svm_load_model<-function(file.model.rdata)
{    
    load(file.model.rdata);       

    gdm <- model$gdm;
    model <- model$svm;
    class(model)<-"gtsvm";
    
    model$decision.values <- NULL;
    model$fitted <- NULL;
    model$residuals <- NULL;

    model <- Rgtsvm::predict.load( model );
    gc(reset=TRUE);
    
    return(list(model=model, gdm=gdm));
}

get_bedgraph_read <- function (dnase_bed, file.imputed)
{
    file.tmp <- tempfile(fileext=".bed")
    write.bed(dnase_bed, file=file.tmp);
    
    tb <- read.table(pipe(paste("zcat ", file.imputed, " | bedtools intersect -a ", file.tmp, "-b -  -loj" )) );
    unlink(file.tmp);
    
    return(c(as.numeric(tb[,7])));
}

genomewide_cor <- function(file.plus, file.minus, file.histone, file.imputed, chr.include=NULL, chr.exclude=c("chrX", "chrY", "chrM"), sample_rate = 0.1 )
{
   if( is.na(file.histone) || is.na(file.imputed) )
       return(c( pearson=NA, spearman=NA ));
   
   if( !all(file.exists(c(file.histone,file.imputed))))
       return(c( pearson=NA, spearman=NA ));
   
   chrom.info.table <- get.chromosome.info( file.plus, file.minus );

   dnase_bed <- as.data.frame(rbindlist( lapply(1:NROW(chrom.info.table), function(i){ data.frame(chrom.info.table[i,1], seq(1, (chrom.info.table[i,2]-1), 10), seq(1, (chrom.info.table[i,2]-1), 10)+1) }) ));
   colnames(dnase_bed)<-NULL

   if(!is.null(chr.include))
       dnase_bed <- dnase_bed[ as.character(dnase_bed[,1]) %in% c( chr.include ),];

   if(!is.null(chr.exclude))
       dnase_bed <- dnase_bed[ !as.character(dnase_bed[,1]) %in% c( chr.exclude ),];

   if(sample_rate!=1)
   {
       idx <- sample( NROW(dnase_bed) )[ 1:round(NROW(dnase_bed)*sample_rate)];
       dnase_bed <- dnase_bed[sort(idx),];
   }

   err <- try(load.bigWig(file.imputed), silent=TRUE);
   if(class(err)=="try-error")
     pred <- get_bedgraph_read(dnase_bed, file.imputed)
   else
     pred <- get_histone_read(dnase_bed, file.imputed);
   
   label <- get_histone_read(dnase_bed, file.histone);
   
   r.cor1 <- cor(pred, label, method ="pearson");
   r.cor2 <- cor(pred, label, method ="spearman");
   cat("COR=", r.cor1, r.cor2, "\n" );

   return(c(pearson=r.cor1,  spearman=r.cor2 ));
}

run_svm_pred<-function(gdm, bed3, file.bw.plus, file.bw.minus, model, temp.bg, ncores, linear_scale) {

    stopifnot(NROW(gdm@window_sizes) == NROW(gdm@half_nWindows))
    zoom<- list(as.integer(gdm@window_sizes), as.integer(gdm@half_nWindows))

cat(".....BED DIM=", NROW(bed3), "\n");

    fp <- extract_feature_matrix( bed3, file.bw.plus, file.bw.minus, gdm, linear_scale, ncores=ncores  )
    pred <- Rgtsvm::predict.run(model, fp$mat )
    ret <- cbind( fp$pos, pred );

    stopifnot(nrow(ret)==NROW(bed3))
    
    rm(fp);
    gc(verbose=TRUE, reset=TRUE)
    if(!is.null(temp.bg))
    {
        options(scipen =99)
        write.table( ret, file=temp.bg, quote=F,sep="\t",col.names=F,row.names=F,append = TRUE)
        rm(ret);
    }
    else
    {
        colnames(ret) <- c("chr", "start", "end", "pred");
        return(ret);
    }    
}

svm_predict_bed <- function( file.imputed.bed, file.bw.plus, file.bw.minus, model, gdm, loci_bed, ncores=1, linear_scale=F  )
{
    if( file.exists(file.imputed.bed) )
    {
        tb.size <-  read.table(pipe(paste("wc -l ",file.imputed.bed)))$V1;
        if(tb.size < NROW(loci_bed))
        {
            cat(file.imputed.bed, "is existing\n");
            loci_bed <- loci_bed[ -c(1:tb.size),]
            gc(reset=TRUE);
        }
        else
        {
            cat(file.imputed.bed, "is done\n");
            return(TRUE);
        }
    }

    line.cutoff <- seq(1, NROW(loci_bed), 1000*100);
    if ( line.cutoff[NROW(line.cutoff)] < NROW(loci_bed)+1 ) 
        line.cutoff <- c(line.cutoff,NROW(loci_bed)+1)

    cat("number of bed =",NROW(loci_bed),"\n")
    cat("number of total blocks=",NROW(line.cutoff)-1,"\n")

    cpu.fun<-function(idx){
        print(idx);
        loci_bed_part <- as.data.frame( loci_bed[c(line.cutoff[idx]:(line.cutoff[idx+1]-1)),,drop=F] );
        run_svm_pred(gdm=gdm, bed3=loci_bed_part,
                    file.bw.plus= file.bw.plus, file.bw.minus= file.bw.minus,
                    model= model, temp.bg= file.imputed.bed,
                    ncores= ncores, linear_scale=linear_scale);
        rm( loci_bed_part );
        gc(verbose=TRUE, reset=TRUE);
        return(NULL);
      }

    lapply(c(1:(NROW(line.cutoff)-1)),FUN= cpu.fun)

    return( TRUE );
}

svm_predict_chr <- function( file.imputed.bed, file.bw.plus, file.bw.minus, model, gdm, chr, ncores=1, linear_scale=F  )
{
    chrom.info.table <- get.chromosome.info( file.bw.plus, file.bw.minus );
    if( NROW(which(chr==as.character(chrom.info.table[,1])))==0 )
       return(FALSE);

    idx <- which(chr==as.character(chrom.info.table[,1]));
cat(idx, file="debug.out", append=TRUE);    
    loci<- seq(1, chrom.info.table[idx,2], 10);

    if( file.exists(file.imputed.bed) )
    {
        tb.size <-  read.table(pipe(paste("wc -l ",file.imputed.bed)))$V1;
        if( tb.size < NROW(loci))
        {
            cat(file.imputed.bed, "is existing\n");
            loci <- loci[ -c(1:tb.size)]
            gc(reset=TRUE);
        }
        else
        {
            cat(file.imputed.bed, "is done\n");
            return(TRUE);
        }
    }

    if (NROW(loci)<=1) return(TRUE);
       
    line.cutoff <- seq(1, NROW(loci), 1000*100);
    if (line.cutoff[NROW(line.cutoff)] <= NROW(loci) ) 
       line.cutoff <- c(line.cutoff, NROW(loci)+1 );

    cat("number of bed =",NROW(loci), file="debug.out", append=TRUE, "\n")
    cat("number of total blocks=",NROW(line.cutoff)-1,file="debug.out", append=TRUE,"\n")

    cpu.fun<-function(idx ){
        print(idx);
        start.p <- unique(loci[c(line.cutoff[idx]:(line.cutoff[idx+1]-1))]);
        loci_bed_part <- data.frame( chr, start.p, start.p + 1 );
show(head(loci_bed_part));

        run_svm_pred(gdm=gdm, bed3=loci_bed_part,
                    file.bw.plus= file.bw.plus, file.bw.minus= file.bw.minus,
                    model= model,temp.bg= file.imputed.bed,
                    ncores= ncores, linear_scale=linear_scale);
          rm(loci_bed_part);
        gc(verbose=TRUE, reset=TRUE);
        return(NULL);
      }

    lapply(c(1:(NROW(line.cutoff)-1)),FUN= cpu.fun)

    return( TRUE );
}

svm_predict_chrs<-function(impute_prefix, file.bw.plus, file.bw.minus, file.model.rdata, chrs, bigwig_compare=NULL, ncores=1, linear_scale=F  )
{
    model <- svm_load_model(file.model.rdata);       

    file.beds <- pred.chrs <- c();
    for( chr0 in chrs )
    {
cat("CHR=", chr0, "\n");

        file.imputed.bed <- paste(impute_prefix, "_", chr0, ".bed", sep="");
        bSucc <- svm_predict_chr( file.imputed.bed, file.bw.plus, file.bw.minus, model$model, model$gdm, chr0, ncores=ncores, linear_scale=linear_scale);
        if(bSucc) 
        {
            file.beds <- c(file.beds, file.imputed.bed );
            pred.chrs <- c(pred.chrs, chr0);
        }    
    }

    Rgtsvm::predict.unload( model$model );
    rm(model);
    
    r.cors <- c();
    for( file.bed in file.beds )
    {
          if( file.exists( file.bed ) && !is.null(bigwig_compare) )
        {
            r.cor <- compare_correlation(  file.bed, bigwig_compare )
            r.cors <- rbind(r.cors, r.cor);
        }
    }

    rownames(r.cors) <- pred.chrs;

    return(list(file.bed=file.beds, r.cors=r.cors));
}


svm_predict_all<-function(impute_prefix, file.bw.plus, file.bw.minus, file.model.rdata, bigwig_compare=NULL, linear_scale=F, ncores=1, gpu.idx=0)
{
    chrom.info.table <- get.chromosome.info( file.bw.plus, file.bw.minus );
    chrs.all <- c(); 
    for(chr in sort(as.character(chrom.info.table[,1] )) )
    {
        file.imputed.bed <- paste(impute_prefix, "_", chr, ".bed", sep="");
          if(!file.exists(file.imputed.bed))
        {
            chrs.all <- c(chrs.all, chr);
        }
        else
        {
            chr.size <- chrom.info.table[which(chrom.info.table[,1]==chr),2]
              chr.sample <- seq(1, chr.size-1, 10);
               tb.size <-  read.table(pipe(paste("wc -l ",file.imputed.bed)))$V1;
            if( tb.size < NROW(chr.sample)-1 )
                  chrs.all <- c(chrs.all, chr);
        }  
    }

    chrs.all <- setdiff(chrs.all, c("chrX", "chrM"))
    if(NROW(chrs.all)==0)
    {
         cat("***** All chromosomes are predicted.\n");
         return(0);
    }
    chrs.all <- sample(chrs.all);

print(chrs.all);
    
    cpu.fun<-function(gpu.cur.idx)
    {
        source("hist.svm.com.R")
        source("hist.svm.pred.R")
        source("hist.svm.main.R")
        
        library(dREG);
        library(bigWig);
        library(Rgtsvm);
        selectGPUdevice(gpu.cur.idx);

        task.idx <- c(0:ceiling(NROW(chrs.all)/NROW(gpu.idx)))*NROW(gpu.idx) + which(gpu.cur.idx==gpu.idx);
        task.idx <- task.idx[ task.idx <= NROW(chrs.all) ];

        model <- svm_load_model(file.model.rdata);       

        file.beds <- c();
        for(tix in task.idx )
        {
            file.imputed.bed <- paste(impute_prefix, "_", chrs.all[tix], ".bed", sep="");
cat("output file=", file.imputed.bed, "\n", file="debug.out",append=TRUE);

            bSucc <- svm_predict_chr( file.imputed.bed, file.bw.plus, file.bw.minus, model$model, model$gdm, chrs.all[tix], ncores=ncores, linear_scale=linear_scale);
            if(bSucc) file.beds <- c(file.beds, file.imputed.bed );
              
              gc(reset=TRUE);
        }  

        Rgtsvm::predict.unload( model$model );
        
        rm(model);
        gc();
        
        return(file.beds)
    }

    if(NROW(gpu.idx)>1)
    {
        library(snowfall);
        sfInit(parallel = TRUE, cpus = min( NROW(gpu.idx), NROW(chrs.all) ), type = "SOCK" )
        sfExport("file.model.rdata", "file.bw.plus", "file.bw.minus", "impute_prefix", "ncores", "chrs.all", "gpu.idx", "linear_scale" );

        fun <- as.function(cpu.fun);
        environment(fun)<-globalenv();

        file.beds<- unlist(sfLapply( gpu.idx, fun));
        sfStop();
    }
    else
       file.bed <- cpu.fun( gpu.idx );

    pred.cors <- NA;
    bigwig_compare <- NULL;
    if( !is.null(bigwig_compare) && file.exists(bigwig_compare) ) 
    {
          r.cors <- mclapply(sort(as.character(chrom.info.table[,1] )), function(chr0){
               file.imputed.bed <- paste(impute_prefix, "_", chr0, ".bed", sep="");
              if( file.exists( file.imputed.bed ) )
            {
                   r.cor<-compare_correlation( file.imputed.bed, bigwig_compare );
            }  
            else
                r.cor<-c(NA,NA);
        
             return(data.frame(chr=chr0, cor1=r.cor[1], cor2=r.cor[2]));
          }, mc.cores=ncores);

        pred.cors <- do.call("rbind", r.cors );
    }
    
    return(pred.cors);
}

svm_predict_chr_parallel<-function( impute_prefix, file.bw.plus, file.bw.minus, file.model.rdata, chr, bigwig_compare=NULL, ncores=1, linear_scale=F, gpu.idx=0 )
{
    chrom.info.table <- get.chromosome.info( file.bw.plus, file.bw.minus );
    seq.start <- 1;
    seq.end <- chrom.info.table[which(chr==chrom.info.table[,1]),2];

    cpu.fun<-function(gpu.cur.idx)
    {
        source("hist.svm.com.R")
        source("hist.svm.pred.R")
        source("hist.svm.main.R")
        
        library(dREG);
        library(bigWig);
        library(Rgtsvm);
        selectGPUdevice(gpu.cur.idx);

        model <- svm_load_model(file.model.rdata);       

        loci <- seq(seq.start, seq.end, 10);
        loci.block <- seq(1, NROW(loci), 1000*100);
        if( loci.block[NROW(loci.block)]<NROW(loci)) 
            loci.block <- c( loci.block, NROW(loci) );

        task.idx <- c(0:ceiling((NROW(loci.block)-1)/NROW(gpu.idx)))*NROW(gpu.idx) + which(gpu.cur.idx==gpu.idx);
        task.idx <- task.idx[ task.idx <= (NROW(loci.block)-1) ];

        ret.list <- list();
        i <- 1;
        for(tix in task.idx )
        { 
            idx.start <- loci.block [ tix] 
            idx.end   <- loci.block [ tix+1 ] - 1;
            cat(range(loci[idx.start:idx.end]), "\n");

            bed3 = data.frame(chr, loci[idx.start:idx.end], loci[idx.start:idx.end]+1 ); 
            ret.list[[i]] = run_svm_pred( model$gdm, bed3, file.bw.plus, file.bw.minus, model$model, NULL, ncores=ncores, linear_scale=linear_scale) 
            i = i+1;
              gc(reset=TRUE);
        }  

        Rgtsvm::predict.unload( model$model );
        
        rm(model);
        gc();
        
        return(do.call("rbind", ret.list) )
    }
    
    if(NROW(gpu.idx)>1)
    {
        library(snowfall);
        sfInit(parallel = TRUE, cpus = NROW(gpu.idx), type = "SOCK" )
        sfExport("file.model.rdata", "file.bw.plus", "file.bw.minus", "impute_prefix", "ncores", "chr", "seq.start", "seq.end", "gpu.idx", "linear_scale" );

        fun <- as.function(cpu.fun);
        environment(fun) <- globalenv();

        pred.bed <- rbindlist( sfLapply( gpu.idx, fun) );
        
        sfStop();
    }
    else
       pred.bed <- cpu.fun( gpu.idx );

    pred.bed <- pred.bed[order(pred.bed[,2]),]
    file.imputed.bed <- paste(impute_prefix, "_", chr, ".bed", sep="");
    write.table(pred.bed, file=file.imputed.bed, row.names=F, col.names=F, quote=F, sep="\t");
   
    r.cor <- NA;
    if( !is.null(bigwig_compare) && file.exists( file.imputed.bed ) )
    {
        r.cor <- compare_correlation( file.imputed.bed, bigwig_compare )
    }

    return(r.cor);

}

# different interface 
# output 5 columns including predict value and experiment value; 
svm_predict_bed2<-function( bed.pred, file.rdata.model, file.bw.histone, file.bw.plus, file.bw.minus, file.peak.pred.gz = NULL, ncores=10)
{
    if(is.null(file.peak.pred.gz)) 
       file.peak.pred.gz <- tempfile(fileext=".bed.gz");

cat("NROW(bed.pred)", NROW(bed.pred), "\n");
cat("Model:", file.rdata.model, "\n");
cat("Bed file:", file.peak.pred.gz, "\n");

    y_exp <- get_histone_read(bed.pred, file.bw.histone, block=500000, ncores=ncores );
    gc(reset=TRUE);

    model <- svm_load_model( file.rdata.model );   
    file.temp.bed <- tempfile(fileext=".bed");
    svm_predict_bed( file.temp.bed, file.bw.plus, file.bw.minus, model$model, model$gdm, bed.pred, ncores=10, linear_scale=F  );
    gc(reset=TRUE);

    tb.pred <- read.table(file.temp.bed);
    write.bed( data.frame(tb.pred, y_exp=y_exp), file.peak.pred.gz, compress=TRUE);
    rm(tb.pred);
    gc(reset=TRUE);
    
    Rgtsvm::predict.unload( model$model );
    rm(model);
    gc(reset=TRUE);
    
    return(file.peak.pred.gz);
}

