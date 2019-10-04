library(dREG)
library(Rgtsvm)
library(bigWig)
library(data.table);
library(snowfall);
library(tools)
library(parallel)
library(MASS);

create_h2h_model <- function( file.bw.train, file.bw.predict,
     file.peak.histone=NULL, 
     bed.blacklist=NULL,
     ratio = 0.1, samples=600000,strategy=1, 
     window.size = 0,
     include.chr="chr1",
     exclude.chr=NULL )
{
        model <- list(
        src=list(
            file.bw.train = file.bw.train,
            file.bw.predict = file.bw.predict,
            file.peak.histone = file.peak.histone ),
        bed.blacklist = bed.blacklist,
        infp = list(),
        train = list(),
        test= list(),
        window.size =window.size,
        ratio = ratio,
        exclude.chr = exclude.chr,
        strategy = strategy,
        include.chr = include.chr,
        samples = samples );

    return(model);
}

get_informative_region<-function(file.bw, chr,  step=50 )
{
    bw <- load.bigWig(file.bw[1]);
    chromsize <- bw$chromSizes[which(bw$chroms==chr)]

    infw <- step.bpQuery.bigWig(bw, chr, 1, floor(chromsize/step)*step, step)
    unload.bigWig(bw);

    if(NROW(file.bw)>1)
       for(j in 2:NROW(file.bw) )
       {
           bw0 <- load.bigWig(file.bw[j]);
           infw <- infw  + abs(step.bpQuery.bigWig(bw0, chr, 1, floor(chromsize/step)*step, step) )
           unload.bigWig(bw0);
       }   
    
    infw.bed<- data.frame(chr=chr, stop=seq(step+1, floor(chromsize/step)*step, step), reads=infw);
    infw.bed$start<- infw.bed$stop-step;
    infw.bed$stop<- infw.bed$start-1;
    infw.bed <- infw.bed[, c("chr", "start", "stop", "reads")];
    
    return(infw.bed);
}

build_h2h_model <- function( model )
{
    infp <- get_informative_region( model$src$file.bw.train, model$include.chr, step=50 );
    model$infp_train <- infp;

    infp <- get_informative_region( model$src$file.bw.predict, model$include.chr, step=50 );
    model$infp_predict <- infp;

    r <- select_train_sample_S1( model, i );
    model$train <- r$train;
    model$test  <- r$test;

    return( model );
}


select_train_sample_S1<-function( model, index )
{
    ## Pro-seq infp, positive/negative region, histone peak region
    ##
    ## (1*) Train Positive * Predict Positive (40%)
    ## (2 ) Train Positive * Predict Negative (20%)
    ## (3 ) Train Negative * Predict Positive (20%)
    ## (4 ) Train Negative * Predict Negative (20%)

    train <- test <- list(pos=NULL, y=NULL);
    r1 <- 1/3; r2 <- 1/3; r3 <- 1/3;
    
    m_train <- mean(model$infp_train[,4] )
    m_pred  <- mean(model$infp_predict[,4])
     
    idx.top1 <- which(model$infp_train[,4]>m_train | model$infp_predict[,4]>m_pred); 
    infp.pos1 <- model$infp_train[ sample(idx.top1, min(model$samples*r1*2, NROW(idx.top1) ) ), ];
    infp.pos1[,2] <- infp.pos1[,2] + round(runif(NROW(infp.pos1), 0, 49))
    infp.pos1[,3] <- infp.pos1[,2] + 1

    idx.top3 <- which(model$infp_train[,4]<m_train/5 & model$infp_predict[,4]<m_pred/5 ); 
    infp.pos3 <- model$infp_train[ sample(idx.top3, min(model$samples*r3*2, NROW(idx.top3) ) ), ];
    infp.pos3[,2] <- infp.pos3[,2] + round(runif(NROW(infp.pos3), 0, 49))
    infp.pos3[,3] <- infp.pos3[,2] + 1

    idx.top2 <- setdiff(1:NROW(model$infp_train), union(idx.top1, idx.top3));
    infp.pos2 <- model$infp_train[ sample(idx.top2, model$sample*2- NROW(infp.pos1)-NROW(infp.pos3) ), ];
    infp.pos2[,2] <- infp.pos2[,2] + round(runif(NROW(infp.pos2), 0, 49))
    infp.pos2[,3] <- infp.pos2[,2] + 1
    
    train$pos <- rbind(infp.pos1[1:(NROW(infp.pos1)/2),], infp.pos2[1:(NROW(infp.pos2)/2),], infp.pos3[1:(NROW(infp.pos3)/2),])
    test$pos <-  rbind(infp.pos1[-c(1:(NROW(infp.pos1)/2)),], infp.pos2[-c(1:(NROW(infp.pos2)/2)),], infp.pos3[-c(1:(NROW(infp.pos3)/2)),])

    train$pos <- train$pos[order(train$pos[,1], train$pos[,2]),]
    test$pos  <- test$pos[order(test$pos[,1], test$pos[,2]),]

    train$y <- get_histone_read(train$pos, model$src$file.bw.predict);
    test$y  <- get_histone_read(test$pos, model$src$file.bw.predict);

    return(list(train=train, test=test));
}

h2h_train_model <- function( model, gdm, file.rdata, ncores=1, skip.svm=FALSE )
{
    y_train <- model$train$y;

    if(NROW(model$src$file.bw.train)==1)
    {
       x_train <- extract_feature_matrix( model$train$pos, model$src$file.bw.train, model$src$file.bw.train, gdm, linear_scale=F, ncores=ncores)
       x_train <- x_train$mat[,c(1:(NCOL(x_train$mat)/2))]
    }
    else   
    {
       x_train <- data.frame(row.names=1:NROW(model$train$pos));
       for(i in 1:NROW(model$src$file.bw.train))
       {
          x_train2 <- extract_feature_matrix( model$train$pos, model$src$file.bw.train[i], model$src$file.bw.train[i], gdm, linear_scale=F, ncores=ncores)
          x_train <- cbind(x_train, x_train2$mat[,c(1:(NCOL(x_train2$mat)/2))]);
       }    
    }

    #save(model, gdm, file.rdata, x_train, y_train, x_test, y_test, file=paste(file.rdata, ".trainxy.rdata", sep="") );

if(!skip.svm)
{
    gamma<- 1/NCOL(x_train);
    cost <- 1;
    epsilon<- 0.1;
    scaled <- TRUE;

    bigmx <- attach.bigmatrix( data=as.matrix(x_train) );
    model$gdm <- gdm;
    rm(x_train);
    gc();

    model$svm <- Rgtsvm::svm(bigmx, y_train , type="eps-regression", cost=cost, gamma=gamma, epsilon=epsilon, scale=scaled, verbose=TRUE );
    save(model, file=file.rdata);

    rm(y_train);
    rm(bigmx);
    gc();

    y_test <- model$test$y;
    x_test <- data.frame(row.names=1:NROW(model$test$pos));
    if(NROW(model$src$file.bw.train)==1)
    {
       x_test <- extract_feature_matrix( model$test$pos, model$src$file.bw.train, model$src$file.bw.train, gdm, linear_scale=F, ncores=ncores)
       x_test <- x_test$mat[,c(1:(NCOL(x_test$mat)/2))]
    }
    else   
    {
       for(i in 1:NROW(model$src$file.bw.train))
       {
          x_test2 <- extract_feature_matrix( model$test$pos, model$src$file.bw.train[i], model$src$file.bw.train[i], gdm, linear_scale=F, ncores=ncores)
          x_test <- cbind(x_test, x_test2$mat[,c(1:(NCOL(x_test2$mat)/2))]);
       }    
    }

    model$y_pred <- Rgtsvm::predict.gtsvm(model$svm, x_test);
    model$y_test <- y_test;
    model$pred.cor <- cor( y_test, model$y_pred );
    save(model, file=file.rdata);

    cat("PRED COR=", model$pred.cor, "\n");
    rm(x_test);
    rm(y_test);
    gc();
}

    return(model);
}


h2h_predict_chr_parallel <- function( file.rdata.model, file.bw.train, file.bw.predict, chr="chr22", ncores=4, linear_scale=F)
{
    cat("Predict ", chr, "on", file.bw.train, "\n");

    bw <- load.bigWig(file.bw.train[1]);
    chromsize <- bw$chromSizes[which(bw$chroms==chr)]
    bed.infw<- data.frame(chr=chr, start=seq(1, floor(chromsize/10)*10, 10));
    bed.infw$stop <- bed.infw$start+1;
    bed.infw <- bed.infw[, c("chr", "start", "stop")];
    unload.bigWig(bw);
    rm(bw);
    
    y_test <- get_histone_read(bed.infw, file.bw.predict);

    load( file.rdata.model );
    model$svm <- Rgtsvm::predict.load( model$svm );

    sect <- unique(c(seq(1, NROW(bed.infw)+1, 500000), NROW(bed.infw)+1 ));
    y_pred <- c();
    for(i in 1:(NROW(sect)-1))
    {
       x_test <- data.frame(row.names=sect[i]:(sect[i+1]-1) );
       for(j in 1:NROW(model$src$file.bw.train))
       {
           x_test2 <- extract_feature_matrix( bed.infw[sect[i]:(sect[i+1]-1),], file.bw.train[j], file.bw.train[j], model$gdm, linear_scale=F, ncores=ncores )
           x_test <- cbind(x_test, x_test2$mat[,c(1:(NCOL(x_test2$mat)/2))]);
        }   

        y_pred <- c(y_pred, Rgtsvm::predict.gtsvm(model$svm, x_test));
        if(i%%3==0) gc();
     }
     
    Rgtsvm::predict.unload( model$svm );
    rm(model);
    rm(x_test);
    gc(reset=TRUE);

    pred.cor <- cor( y_pred, y_test );
cat("pred of chr.", chr, "=", pred.cor, "\n"); 

    return(data.frame(bed.infw, y_pred, y_test));
}
