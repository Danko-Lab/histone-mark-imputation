library(dREG)
library(Rgtsvm)
library(bigWig)
library(data.table);
library(snowfall);
library(tools)
library(parallel)
library(MASS);

create_train_model <- function( path.proseq, file.bws.plus, file.bws.minus, 
     file.rdata.negative, file.rdata.positive, 
     path.histone, file.bw.histone, file.peak.histone, 
     bed.blacklist=NULL,
     file.all.histone=NULL,
     ratio = 0.1, samples=600000,strategy=0, 
     window.size = 0,
     include.bed=NULL,
     exclude=c("chrY", "chrM") )
{
    model <- list(
        src=list(path.proseq = path.proseq,
            file.bws.plus = if (is.null(path.proseq) || path.proseq=="") file.bws.plus else paste(path.proseq, file.bws.plus, sep="/"),
            file.bws.minus = if (is.null(path.proseq) || path.proseq=="") file.bws.minus else paste(path.proseq, file.bws.minus, sep="/"),
            file.rdata.negative = if (is.null(path.proseq) || path.proseq=="") file.rdata.negative else paste(path.proseq, file.rdata.negative, sep="/"),
            file.rdata.positive = if (is.null(path.proseq) || path.proseq=="") file.rdata.positive else paste(path.proseq, file.rdata.positive, sep="/"),
            path.histone = path.histone,
            file.bw.histone = if (is.null(path.histone) || path.histone=="") file.bw.histone else paste(path.histone, file.bw.histone, sep="/"),
            file.peak.histone = if (is.null(path.histone) || path.histone=="") file.peak.histone else paste(path.histone, file.peak.histone, sep="/"),
            file.all.histone = if (is.null(path.histone) || path.histone=="") file.all.histone else paste(path.histone, file.all.histone, sep="/")
            ),
        bed.blacklist = bed.blacklist,
        infp = list(),
        train = list(),
        test= list(),
        window.size =window.size,
        ratio = ratio,
        exclude = exclude,
        strategy = strategy,
        include.bed = include.bed,
        samples = samples );

    files1 <- c(file.bws.plus, file.bws.minus, file.rdata.negative, file.rdata.positive, file.bw.histone, file.peak.histone)
    for(f in files1)
        if( !file.exists(f) )
            stop("file doesn't exist", f, "\n")

    return(model);
}

build_train_model <- function( model )
{
    for(i in 1:NROW(model$src$file.bws.plus))
        if (length(model$infp)<i || is.null( model$infp[[i]]) )
        {
            file.bw.plus <- model$src$file.bws.plus[i];
            file.bw.minus <- model$src$file.bws.minus[i];
            infp <- get_informative_positions( file.bw.plus, file.bw.minus, depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE );
            infp <- infp[grep("_|chrM|chrY", infp[,1], invert=TRUE),]
            model$infp[[i]] <- infp;
        }

    for(i in 1:length(model$infp))
    {
        r <- select_train_sample( model, i );
        model$train[[i]] <- r$train;
        model$test[[i]]<- r$test;
    }

    return( model );
}

select_train_sample<-function( model, index )
{
    if(model$strategy==1)
        return(select_train_sample_S1(model,index))
    if(model$strategy==2)
        return(select_train_sample_S2(model,index))
    if(model$strategy==3)
        return(select_train_sample_S3(model,index, model$group.bw.histone, model$bed.black.region))
    if(model$strategy==4)
        return(select_train_sample_S4(model,index, model$group.bw.histone, model$bed.black.region, model$window.size))
    if(model$strategy==5)
        return(select_train_sample_S5(model,index))
    if(model$strategy==6)
        return(select_train_sample_S6(model,index))
}

get_histone_peak<-function( file.peak, file.histone)
{
   bw.hist <- load.bigWig(file.histone);
   bedA <- read.table(file.peak)[,c(1:3)];

   y <- try( unlist(bed.region.bpQuery.bigWig( bw=bw.hist, bed= bedA, op="max") ) );
   if(class(y)=="try-error")
   {
     stop("error");
     return(NULL);
   }
   unload.bigWig(bw.hist);
   return(data.frame(bedA, y));
}


get_proseq_reads<-function(bed.pos, file.bw.plus, file.bw.minus, half.width=50000)
{
    bed.pos[,2] <- bed.pos[,2] - half.width;
    bed.pos[,3] <- bed.pos[,3] + half.width;
    if( NROW(which(bed.pos[,2]<=0))>0 )
        bed.pos[which(bed.pos[,2]<=0), 2 ]<-0;
    
    bw.plus <- load.bigWig(file.bw.plus);
    bw.minus <-  load.bigWig(file.bw.minus);
    reads <- bed.region.bpQuery.bigWig(bw.plus, bed.pos, op = "sum", abs.value = TRUE) + bed.region.bpQuery.bigWig(bw.minus, bed.pos,  op = "sum", abs.value = TRUE);
    unload.bigWig(bw.plus);
    unload.bigWig(bw.minus);
    
    return(reads);
}

select_median_marks <- function(y_proseq, df_mark)
{
   for(j in 1:NCOL(df_mark))
      df_mark[,j] <- df_mark[,j]/max(df_mark[,j])*max(df_mark[,1]);

show(summary(y_proseq));
browser();

   yp_0.25 <- quantile(y_proseq, 0.25);
   yp_0.90 <- quantile(y_proseq, 0.90);
   
   ret <- c();
   for(i in 1:NROW(df_mark))
   {
      if(y_proseq[i]<yp_0.25)
         ret[i] <- min(df_mark[i,])
      else if(y_proseq[i]>yp_0.90)
         ret[i] <- max(df_mark[i,])
      else 
         ret[i] <- median(df_mark[i,])
   }   
   return(ret);
}

select_train_sample_S3<-function( model, index, files.bw.histone, bed.black.region=NULL )
{
    ## Pro-seq infp, positive/negative region, histone peak region
    ##
    ## (1*) Pro-seq infp & Positive (5%)
    ## (2 ) Pro-seq infp & Negative (85%)
    ## (3 ) Negative only (10%)

    train <- test <- list(pos=NULL, y=NULL);
    r1 <- 0.05; r2 <- 0.85; r3 <- 0.10;
    fmat0 <- fmat1 <- c();

    if( !file.exists( model$src$file.bw.histone) ) 
        stop(paste( "No histone file.", model$src$file.bw.histone, sep="") );

    if(file_ext(model$src$file.rdata.negative)=="rdata")
    {
        load(model$src$file.rdata.negative);
        load(model$src$file.rdata.positive);
    }
    else
    {
        positive_bed <- read.table(model$src$file.rdata.positive)[,c(1:3)]
        negative_bed <- read.table(model$src$file.rdata.negative)[,c(1:3)]
    }

    if(!is.null(bed.black.region))
    {
        positive_bed <- bedTools.subtract( positive_bed, bed.black.region);
        negative_bed <- bedTools.subtract( negative_bed, bed.black.region);
    }

    infp  <- model$infp[[index]];

    filter <- "";
    if( !is.null( model$exclude ) && model$exclude!="" )
        filter <- paste("_",  paste(model$exclude, collapse="|"), sep="|");

    if( filter != "" )
    {
        positive_bed <- positive_bed[grep(filter, positive_bed[,1], invert=TRUE),];
         negative_bed <- negative_bed[grep(filter, negative_bed[,1], invert=TRUE),];
         infp         <- infp[grep(filter, infp[,1], invert=TRUE),];
    }

    #part 1:
    tb1 <- bedTools.intersect( positive_bed, infp);
    tb1 <- unique(tb1[ sample( 1:NROW(tb1) ),c(1:3)]);
    colnames(tb1) <- c("chr", "start", "end");
    fmat0 <- rbind( fmat0, tb1[1:round(r1*model$samples),] );
    fmat1 <- rbind( fmat1, tb1[-c(1:round(r1*model$samples)),][1:round(r1*model$samples),]  );
cat("(1*)=", NROW(unique(fmat0)), "\n");

    #part 2:
    tb1 <- bedTools.intersect( negative_bed, infp);
    tb1 <- unique(tb1[ sample( 1:NROW(tb1) ),c(1:3)]);
    colnames(tb1) <- c("chr", "start", "end");
    fmat0 <- rbind( fmat0, tb1[1:round(r2*model$samples),] );
    fmat1 <- rbind( fmat1, tb1[-c(1:round(r2*model$samples)),][1:round(r2*model$samples),]  );
cat("(2 )=", NROW(unique(fmat0)), "\n");

    #part 3:
    tb1 <- bedTools.subtract( negative_bed, infp);
    tb1[,2] <- round((tb1[,2]+tb1[,3])/2);
    tb1[,3] <- tb1[,2] + 1;
    tb1 <- unique(tb1[ sample( 1:NROW(tb1) ),c(1:3)]);
    colnames(tb1) <- c("chr", "start", "end");
    fmat0 <- rbind( fmat0, tb1[1:(model$samples-NROW(unique(fmat0))),] );
    fmat1 <- rbind( fmat1, tb1[-c(1:(model$samples-NROW(unique(fmat0)))),][1:(model$samples-NROW(unique(fmat1))),]  );
cat("(3 )=", NROW(unique(fmat0)), "\n");

    fmat0 <- unique(fmat0);
    fmat1 <- unique(fmat1);
    na.idx.famt0 <- which( is.na(fmat0[,1]) |  is.na(fmat0[,2]) |  is.na(fmat0[,3]) ) ;
    na.idx.famt1 <- which( is.na(fmat0[,1]) |  is.na(fmat1[,2]) |  is.na(fmat1[,3]) ) ;
    if( NROW(na.idx.famt0)>0) fmat0 <- fmat0 [-na.idx.famt0,];
    if( NROW(na.idx.famt1)>0) fmat1 <- fmat1 [-na.idx.famt1,];


cat("(F )=", NROW(unique(fmat0)), "\n");

    train$pos <- fmat0;
    test$pos <- fmat1;

    ytrains <- ytests<-list();
    for(i in 1:NROW(model$src$file.all.histone))
    {
       ytrains[[i]] <- get_histone_read(train$pos, model$src$file.all.histone[i]);
       ytests[[i]]  <- get_histone_read(test$pos,  model$src$file.all.histone[i]);
    }

    yp_train <- get_proseq_reads(train$pos, model$src$file.bws.plus[index], model$src$file.bws.minus[index] )
    yp_test <- get_proseq_reads(test$pos,   model$src$file.bws.plus[index], model$src$file.bws.minus[index] )

    train$y <- select_median_marks(yp_train, do.call("cbind", ytrains) );
    test$y  <- select_median_marks(yp_test, do.call("cbind", ytests) );
   
    return(list(train=train, test=test));
}


select_train_sample_S4<-function( model, index, files.bw.histone, bed.black.region=NULL, window.size=10 )
{
    ## Pro-seq infp, positive/negative region, histone peak region
    ##
    ## (1*) Pro-seq infp & Positive (5%)
    ## (2 ) Pro-seq infp & Negative (93%)
    ## (3 ) Negative only (2%)

    train <- test <- list(pos=NULL, y=NULL);
    r1 <- 0.05; r2 <- 0.93; r3 <- 0.02;
    fmat0 <- fmat1 <- c();

    if( !file.exists(model$src$file.bw.histone ) )
        stop(paste( "No histone file.",  model$src$file.bw.histone) );

    if(file_ext(model$src$file.rdata.negative)=="rdata")
    {
        load(model$src$file.rdata.negative );
        load(model$src$file.rdata.positive );
    }
    else
    {
        positive_bed <- read.table( model$src$file.rdata.positive )[,c(1:3)]
        negative_bed <- read.table( model$src$file.rdata.negative )[,c(1:3)]
    }

    if(!is.null(bed.black.region))
    {
        positive_bed <- bedTools.subtract( positive_bed, bed.black.region);
        negative_bed <- bedTools.subtract( negative_bed, bed.black.region);
    }

    infp  <- model$infp[[index]];

    filter <- "";
    if( !is.null( model$exclude ) && model$exclude!="" )
        filter <- paste("_",  paste(model$exclude, collapse="|"), sep="|");

    if( filter != "" )
    {
        positive_bed <- positive_bed[grep(filter, positive_bed[,1], invert=TRUE),];
         negative_bed <- negative_bed[grep(filter, negative_bed[,1], invert=TRUE),];
         infp         <- infp[grep(filter, infp[,1], invert=TRUE),];
    }

    #part 1:
    tb1 <- bedTools.intersect( positive_bed, infp);
    tb1 <- unique(tb1[ sample( 1:NROW(tb1) ),c(1:3)]);
    colnames(tb1) <- c("chr", "start", "end");
    fmat0 <- rbind( fmat0, tb1[1:round(r1*model$samples),] );
    fmat1 <- rbind( fmat1, tb1[-c(1:round(r1*model$samples)),][1:round(r1*model$samples),]  );
cat("(1*)=", NROW(unique(fmat0)), "\n");

    #part 2:
    tb1 <- bedTools.intersect( negative_bed, infp);
    tb1 <- unique(tb1[ sample( 1:NROW(tb1) ),c(1:3)]);
    colnames(tb1) <- c("chr", "start", "end");
    fmat0 <- rbind( fmat0, tb1[1:round(r2*model$samples),] );
    fmat1 <- rbind( fmat1, tb1[-c(1:round(r2*model$samples)),][1:round(r2*model$samples),]  );
cat("(2 )=", NROW(unique(fmat0)), "\n");

    #part 3:
    tb1 <- bedTools.subtract( negative_bed, infp);
    tb1[,2] <- round((tb1[,2]+tb1[,3])/2);
    tb1[,3] <- tb1[,2] + 1;
    tb1 <- unique(tb1[ sample( 1:NROW(tb1) ),c(1:3)]);
    colnames(tb1) <- c("chr", "start", "end");
    fmat0 <- rbind( fmat0, tb1[1:(model$samples-NROW(unique(fmat0))),] );
    fmat1 <- rbind( fmat1, tb1[-c(1:(model$samples-NROW(unique(fmat0)))),][1:(model$samples-NROW(unique(fmat1))),]  );
cat("(3 )=", NROW(unique(fmat0)), "\n");

    fmat0 <- unique(fmat0);
    fmat1 <- unique(fmat1);
    na.idx.famt0 <- which( is.na(fmat0[,1]) |  is.na(fmat0[,2]) |  is.na(fmat0[,3]) ) ;
    na.idx.famt1 <- which( is.na(fmat0[,1]) |  is.na(fmat1[,2]) |  is.na(fmat1[,3]) ) ;
    if( NROW(na.idx.famt0)>0) fmat0 <- fmat0 [-na.idx.famt0,];
    if( NROW(na.idx.famt1)>0) fmat1 <- fmat1 [-na.idx.famt1,];

cat("(F )=", NROW(unique(fmat0)), "\n");

    train$pos <- fmat0;
    test$pos <- fmat1;

    bed3 <- train$pos;
    bed3[,2] <- bed3[,2] - window.size;
    bed3[,3] <- bed3[,3] + window.size;
    if (NROW(which(bed3[,2]<0))>0)
        bed3[which(bed3[,2]<0),2] <- 0;
    train$y <- (get_histone_read(bed3, model$src$file.bw.histone))/(2*window.size);
    
    bed3 <- test$pos;
    bed3[,2] <- bed3[,2] - window.size;
    bed3[,3] <- bed3[,3] + window.size;
    if (NROW(which(bed3[,2]<0))>0)
        bed3[which(bed3[,2]<0),2] <- 0;
    test$y  <- (get_histone_read(bed3, model$src$file.bw.histone))/(2*window.size);
   
    return(list(train=train, test=test));
}



select_train_sample_S1<-function( model, index )
{
    ## Pro-seq infp, positive/negative region, histone peak region
    ##
    ## (1*) Pro-seq infp & Positive (5%)
    ## (2 ) Pro-seq infp & Negative (93%)
    ## (3 ) Negative only (2%)

    train <- test <- list(pos=NULL, y=NULL);
    r1 <- 0.05; r2 <- 0.93; r3 <- 0.02;
    fmat0 <- fmat1 <- c();

    if( !file.exists( model$src$file.bw.histone ) )
        stop(paste( "No histone file.", model$src$file.bw.histone) );

    if(file_ext(model$src$file.rdata.negative)=="rdata")
    {
        load( model$src$file.rdata.negative );
        load( model$src$file.rdata.positive );
    }
    else
    {
        positive_bed <- read.table(model$src$file.rdata.positive )[,1:3]
        negative_bed <- read.table(model$src$file.rdata.negative )[,1:3]
    }

    infp  <- model$infp[[index]];

    filter <- "";
    if( !is.null( model$exclude ) && model$exclude!="" )
        filter <- paste("_",  paste(model$exclude, collapse="|"), sep="|");

    if( filter != "" )
    {
        positive_bed <- positive_bed[grep(filter, positive_bed[,1], invert=TRUE),];
         negative_bed <- negative_bed[grep(filter, negative_bed[,1], invert=TRUE),];
         infp         <- infp[grep(filter, infp[,1], invert=TRUE),];
    }

    #part 1:
    tb1 <- bedTools.intersect( positive_bed, infp);
    tb1 <- unique(tb1[ sample( 1:NROW(tb1) ),c(1:3)]);
    colnames(tb1) <- c("chr", "start", "end");
    fmat0 <- rbind( fmat0, tb1[1:round(r1*model$samples),] );
    fmat1 <- rbind( fmat1, tb1[-c(1:round(r1*model$samples)),][1:round(r1*model$samples),]  );
cat("(1*)=", NROW(unique(fmat0)), "\n");

    #part 2:
    tb1 <- bedTools.intersect( negative_bed, infp);
    tb1 <- unique(tb1[ sample( 1:NROW(tb1) ),c(1:3)]);
    colnames(tb1) <- c("chr", "start", "end");
    fmat0 <- rbind( fmat0, tb1[1:round(r2*model$samples),] );
    fmat1 <- rbind( fmat1, tb1[-c(1:round(r2*model$samples)),][1:round(r2*model$samples),]  );
cat("(2 )=", NROW(unique(fmat0)), "\n");

    #part 3:
    tb1 <- bedTools.subtract( negative_bed, infp);
    tb1[,2] <- round((tb1[,2]+tb1[,3])/2);
    tb1[,3] <- tb1[,2] + 1;
    tb1 <- unique(tb1[ sample( 1:NROW(tb1) ),c(1:3)]);
    colnames(tb1) <- c("chr", "start", "end");
    fmat0 <- rbind( fmat0, tb1[1:(model$samples-NROW(unique(fmat0))),] );
    fmat1 <- rbind( fmat1, tb1[-c(1:(model$samples-NROW(unique(fmat0)))),][1:(model$samples-NROW(unique(fmat1))),]  );
cat("(3 )=", NROW(unique(fmat0)), "\n");

    fmat0 <- unique(fmat0);
    fmat1 <- unique(fmat1);
    na.idx.famt0 <- which( is.na(fmat0[,1]) |  is.na(fmat0[,2]) |  is.na(fmat0[,3]) ) ;
    na.idx.famt1 <- which( is.na(fmat0[,1]) |  is.na(fmat1[,2]) |  is.na(fmat1[,3]) ) ;
    if( NROW(na.idx.famt0)>0) fmat0 <- fmat0 [-na.idx.famt0,];
    if( NROW(na.idx.famt1)>0) fmat1 <- fmat1 [-na.idx.famt1,];

cat("(F )=", NROW(unique(fmat0)), "\n");

    train$pos <- fmat0;
    test$pos <- fmat1;

    train$y <- get_histone_read(train$pos, model$src$file.bw.histone );
    test$y  <- get_histone_read(test$pos, model$src$file.bw.histone );

    return(list(train=train, test=test));
}

select_train_sample_S2<-function( model, index )
{
    ## Pro-seq infp, positive/negative region, histone peak region
    ##
    ## (1*) Pro-seq infp & Positive & histone peak (10%)
    ## (2*) Pro-seq infp & histone peak (5%)
    ## (3*) Pro-seq infp & Positive (5%)
    ## (4 ) Pro-seq infp & Negative (78%)
    ## (5 ) Negative only (2%)

    train <- test <- list(pos=NULL, y=NULL);
    r1 <- 0.1; r2 <- 0.05; r3 <- 0.05; r4 <- 0.78; r5<- 0.02;
    fmat0 <- fmat1 <- c();

    load( model$src$file.rdata.negative );
    load( model$src$file.rdata.positive );
    infp  <- model$infp[[index]];

    filter <- "";
    if( !is.null( model$exclude ) && model$exclude!="" )
        filter <- paste("_",  paste(model$exclude, collapse="|"), sep="|");

    if( filter != "" )
    {
        positive_bed <- positive_bed[grep(filter, positive_bed[,1], invert=TRUE),];
         negative_bed <- negative_bed[grep(filter, negative_bed[,1], invert=TRUE),];
         infp         <- infp[grep(filter, infp[,1], invert=TRUE),];
    }

    #part 1:
    tb1 <- bedTools.intersect( positive_bed, model$src$file.peak.histone );
    tb1 <- unique(tb1[ sample( 1:NROW(tb1) ), c(1:3) ]);
    tb1 <- tb1 [ tb1[,3] - tb1[,2] > 10, ];
    tb1 <- split.bed(tb1, 10);
    colnames(tb1) <- c("chr", "start", "end");
    fmat0 <- tb1[1:round(r1*model$samples),]
    fmat1 <- tb1[-c(1:round(r1*model$samples)),][1:round(r1*model$samples),];
cat("(1*)=", NROW(unique(fmat0)), "\n");

    #part 2:
    tb1 <- bedTools.subtract( model$src$file.peak.histone, positive_bed);
    tb1 <- bedTools.intersect( tb1[,c(1:3)], infp);
    tb1 <- unique(tb1[ sample( 1:NROW(tb1) ),c(1:3)]);
    colnames(tb1) <- c("chr", "start", "end");
    fmat0 <- rbind( fmat0, tb1[1:round(r2*model$samples),] );
    fmat1 <- rbind( fmat1, tb1[-c(1:round(r2*model$samples)),][1:round(r2*model$samples),]  );
cat("(2*)=", NROW(unique(fmat0)), "\n");

    #part 3:
    tb1 <- bedTools.subtract( positive_bed, model$src$file.peak.histone);
    tb1 <- bedTools.intersect( tb1, infp);
    tb1 <- unique(tb1[ sample( 1:NROW(tb1) ),c(1:3)]);
    colnames(tb1) <- c("chr", "start", "end");
    fmat0 <- rbind( fmat0, tb1[1:round(r3*model$samples),] );
    fmat1 <- rbind( fmat1, tb1[-c(1:round(r3*model$samples)),][1:round(r3*model$samples),]  );
cat("(3*)=", NROW(unique(fmat0)), "\n");

    #part 4:
    tb1 <- bedTools.intersect( negative_bed, infp);
    tb1 <- unique(tb1[ sample( 1:NROW(tb1) ),c(1:3)]);
    colnames(tb1) <- c("chr", "start", "end");
    fmat0 <- rbind( fmat0, tb1[1:round(r4*model$samples),] );
    fmat1 <- rbind( fmat1, tb1[-c(1:round(r4*model$samples)),][1:round(r4*model$samples),]  );
cat("(4 )=", NROW(unique(fmat0)), "\n");

    #part 5:
    tb1 <- bedTools.subtract( negative_bed, infp);
    tb1[,2] <- round((tb1[,2]+tb1[,3])/2);
    tb1[,3] <- tb1[,2] + 1;
    tb1 <- unique(tb1[ sample( 1:NROW(tb1) ),c(1:3)]);
    colnames(tb1) <- c("chr", "start", "end");
    fmat0 <- rbind( fmat0, tb1[1:(model$samples-NROW(unique(fmat0))),] );
    fmat1 <- rbind( fmat1, tb1[-c(1:(model$samples-NROW(unique(fmat0)))),][1:(model$samples-NROW(unique(fmat1))),]  );
cat("(5 )=", NROW(unique(fmat0)), "\n");

    fmat0 <- unique(fmat0);
    fmat1 <- unique(fmat1);
    na.idx.famt0 <- which( is.na(fmat0[,1]) |  is.na(fmat0[,2]) |  is.na(fmat0[,3]) ) ;
    na.idx.famt1 <- which( is.na(fmat0[,1]) |  is.na(fmat1[,2]) |  is.na(fmat1[,3]) ) ;
    if( NROW(na.idx.famt0)>0) fmat0 <- fmat0 [-na.idx.famt0,];
    if( NROW(na.idx.famt1)>0) fmat1 <- fmat1 [-na.idx.famt1,];


cat("(F )=", NROW(unique(fmat0)), "\n");

    train$pos <- fmat0;
    test$pos <- fmat1;

    train$y <- get_histone_read(train$pos, model$src$file.bw.histone);
    test$y  <- get_histone_read(test$pos, model$src$file.bw.histone);

    return(list(train=train, test=test));
}


summary_segnment_bed<-function( file_bw)
{
    read.bigwig.fast <- function(file.hist, tb.bed, op="sum" )
    {
      interval <- unique(c( seq( 1, NROW(tb.bed)+1, by = 1000*100 ), NROW(tb.bed)+1))

      ret <- do.call("c", mclapply(1:(length(interval)-1), function(x)
      {
        bw <- load.bigWig(file.hist);
        batch_indx<- c( interval[x]:(interval[x+1]-1) ) 
        dat <- bed.region.bpQuery.bigWig( bw, tb.bed[batch_indx, ], op=op);
        unload.bigWig(bw);
        return( c(dat) );
      }, mc.cores=24));

      return(ret);
    }


    find.peak.pos <- function(file.hist, tb.bed )
    {
      interval <- unique(c( seq( 1, NROW(tb.bed)+1, by = 1000*10 ), NROW(tb.bed)+1))

      ret <- do.call("rbind", mclapply(1:(length(interval)-1), function(x)
      {
        bw <- load.bigWig(file.hist);
        batch_indx<- c( interval[x]:(interval[x+1]-1) ) 
        dat <- bed.step.probeQuery.bigWig( bw, tb.bed[batch_indx, ], step=1, as.matrix=T );

        max.pos <- unlist(lapply(1:NROW(dat), function(i){which.max(dat[i,])}));
        min.pos <- unlist(lapply(1:NROW(dat), function(i){which.min(dat[i,])}));
        unload.bigWig(bw);
        return( data.frame(max=max.pos, min=min.pos ) );
      }, mc.cores=24));

      return(ret);
    }

    tb.chrom <- read.table("/fs/cbsudanko/storage/data/hg19/hg19.chromInfo");
    tb.chrom <- tb.chrom[grep("_|chrM|chrY|chrX|chr22", tb.chrom[,1], invert=TRUE),]

    tb.bed <-do.call("rbind", lapply(1:NROW(tb.chrom), function(i){
       return(data.frame(chr=tb.chrom[i,1], start=seq(1, tb.chrom[i,2]-100, 100), stop=seq(1, tb.chrom[i,2]-100, 100)+100 ));
    }));

    tb.bed$sum <- read.bigwig.fast( file_bw, tb.bed);
    tb.bed$max <- read.bigwig.fast( file_bw, tb.bed, op="max");
    tb.bed$signal = 0;
    tb.bed$max.pos = NA
    tb.bed$min.pos = NA

    q1 <- 0.0001
    q2 <- 0.05

    F = fitdistr(round(tb.bed$max), densfun="Poisson")
    Q2 = qpois(1-q1, F$estimate);
    tb.bed$signal[ round(tb.bed$max) >= Q2 & tb.bed$sum >= quantile(tb.bed$sum, 1-q2)] <- 1;
    df <- find.peak.pos(file_bw, tb.bed[tb.bed$signal==1,1:3] )
    tb.bed$max.pos[tb.bed$signal==1] <-df[,1];
    tb.bed$min.pos[tb.bed$signal==1] <-df[,2];
    
    return(tb.bed);
}


select_train_sample_S6<-function( model, index, nsample=600000 )
{
    if( !file.exists(model$src$file.bw.histone ) )
        stop(paste( "No histone file.", model$src$file.bw.histone) );

    df.bed <- summary_segnment_bed( model$src$file.bw.histone );

    train <- test <- list(pos=NULL, y=NULL);

    df1 <- df.bed[ df.bed$signal==1,]
    df0 <- df.bed[ df.bed$signal==0,]
     
    df1$offset <-  df1$max.pos - 1;
    df1$start <- df1$start + df1$offset;
    df1$stop <- df1$start + 1;
    
    df0$offset <-  round(runif(NROW(df0), 1, 100 ));
    df0$start  <- df0$start + df0$offset;
    df0$stop   <- df0$start + 1;
 
    df1.tr.idx <- sort(sample(1:NROW(df1))[1:round(NROW(df1)*0.5)]); 
    df1.tt.idx <- sort(c(1:NROW(df1))[-df1.tr.idx]); 
    df0.tr.idx <- sort(sample(1:NROW(df0))[1:round(NROW(df1)*0.5)]); 
    df0.tt.idx <- sort(sample(c(1:NROW(df0))[-df0.tr.idx])[1:NROW(df1.tt.idx)]); 

    train$pos <- rbind( df1[df1.tr.idx, 1:3],  df0[df0.tr.idx, 1:3] );
    test$pos  <- rbind( df1[df1.tt.idx, 1:3], df0[df0.tt.idx, 1:3] );

    train$pos <- train$pos[sort(sample(NROW(train$pos))[1:nsample]),];
    test$pos <- test$pos[sort(sample(NROW(test$pos))[1:nsample]),];

    train$y <- get_histone_read(train$pos, model$src$file.bw.histone);
    test$y  <- get_histone_read(test$pos, model$src$file.bw.histone);

print(head(train$y));
print(head(test$y));

    return(list(train=train, test=test));
}



select_train_sample_S5<-function( model, index, nsample=1000000 )
{
    if( !file.exists(model$src$file.bw.histone) )
        stop(paste( "No histone file.", model$src$file.bw.histone ) );

    train <- test <- list(pos=NULL, y=NULL);

    idx.train <- sample(1:NROW(model$include.bed))[1:round(NROW(model$include.bed)/2)]; 
    df.bed0 <- model$include.bed[idx.train,]
    df.bed1 <- model$include.bed[-idx.train,]
    
    train$pos<-do.call("rbind", lapply(1:NROW(df.bed0), function(i){ data.frame(chr=df.bed0[i,1], start=seq(df.bed0[i,2], df.bed0[i,3], 10), stop=seq(df.bed0[i,2], df.bed0[i,3], 10)+1)} ))
    test$pos <-do.call("rbind", lapply(1:NROW(df.bed1), function(i){ data.frame(chr=df.bed1[i,1], start=seq(df.bed1[i,2], df.bed1[i,3], 10), stop=seq(df.bed1[i,2], df.bed1[i,3], 10)+1)} ))

    train$pos <- train$pos[sample(NROW(train$pos))[1:nsample],];
    test$pos <- test$pos[sample(NROW(test$pos))[1:nsample],];
    
    train$y <- get_histone_read(train$pos, model$src$file.bw.histone);
    test$y  <- get_histone_read(test$pos, model$src$file.bw.histone);

    return(list(train=train, test=test));
}


svm_train_model <- function( model, gdm, file.rdata, ncores=1, skip.svm=FALSE )
{
    x_train <- y_train <- c();
    for(i in 1:length(model$train))
    {
        x <- extract_feature_matrix( model$train[[i]]$pos, model$src$file.bws.plus[i], model$src$file.bws.minus[i], gdm, linear_scale=F, ncores=ncores )
        x_train <- rbind(x_train, x$mat);
        y_train <- c(y_train, model$train[[i]]$y);
        rm(x);
        gc(reset=TRUE);
    }

if(!skip.svm)
{
    gamma<- 1/NCOL(x_train);
    cost <- 1;
    epsilon<- 0.1;
    scaled <- TRUE;

    bigmx <- attach.bigmatrix( data=as.matrix(x_train) );
    model$gdm <- gdm;
    rm(x_train);
    gc(reset=TRUE);

cat("start training...\n");
    model$svm <- Rgtsvm::svm(bigmx, y_train , type="eps-regression", cost=cost, gamma=gamma, epsilon=epsilon, scale=scaled, verbose=TRUE );
    save(model, file=file.rdata);

    rm(y_train);
    rm(bigmx);
    gc(reset=TRUE);

    x_test <- y_test <- c();
    for(i in 1:length(model$test))
    {
        x1 <- extract_feature_matrix( model$test[[i]]$pos, model$src$file.bws.plus[i], model$src$file.bws.minus[i], gdm, linear_scale=F, ncores=ncores )
        x_test <- rbind(x_test, x1$mat);
        y_test <- c(y_test, model$test[[i]]$y);
        rm(x1);
        gc();
    }

    model$y_pred <- Rgtsvm::predict.gtsvm(model$svm, x_test);
    model$y_test <- y_test;
    model$pred.cor <- cor( y_test, model$y_pred );
    save(model, file=file.rdata);

    cat("PRED COR=", model$pred.cor, "\n");
    rm(x_test);
    rm(y_test);
    gc(reset=TRUE);
}

    return(model);
}


svm_train_model2 <- function( model, gdm, file.rdata, uselog=FALSE, ncores=1 )
{
    x_train <- y_train <- c();

    for(i in 1:length(model$train))
    {
        x <- extract_feature_matrix( model$train[[i]]$pos, model$src$file.bws.plus[i], model$src$file.bws.minus[i], gdm, linear_scale=F, ncores=ncores )
        x_train <- rbind(x_train, x$mat);
        y_train <- c(y_train, model$train[[i]]$y);
        rm(x);
        gc();
    }

    gamma<- 1/NCOL(x_train);
    cost <- 1;
    epsilon<- 0.1;
    scaled <- TRUE;

    bigmx <- attach.bigmatrix( data=as.matrix(x_train) );
    model$gdm <- gdm;
    rm(x_train);
    gc();

    if(uselog) y_train <- log(y_train+1);

    model$svm <- Rgtsvm::svm(bigmx, y_train , type="eps-regression", cost=cost, gamma=gamma, epsilon=epsilon, scale=scaled, verbose=TRUE );
    save(model, file=file.rdata);

    rm(y_train);
    rm(bigmx);
    gc();

    x_test <- y_test <- df.test <- c();
    for(i in 1:length(model$test))
    {
        x1 <- extract_feature_matrix( model$test[[i]]$pos, model$src$file.bws.plus[i], model$src$file.bws.minus[i], gdm, linear_scale=F, ncores=ncores )
        x_test <- rbind(x_test, x1$mat);
        y_test <- c(y_test, model$test[[i]]$y);
        rm(x1);
        gc();
    }

    if(uselog) y_test <- log(y_test+1);

    model$uselog <- uselog;
    model$y_pred <- Rgtsvm::predict.gtsvm(model$svm, x_test);
    model$y_test <- y_test;
    model$pred.cor <- cor( y_test, model$y_pred );
    save(model, file=file.rdata);

    cat("PRED COR=", model$pred.cor, "\n");
    rm(x_test);
    rm(y_test);
    gc();

    return(model);
}

if(0)
{
load("H3k27ac.train.rdata")
model$svm <- NULL
model <- svm_train_model2 ( model, model$gdm, "H3k27ac.train.log.rdata", uselog=TRUE, ncores=15 )
}
