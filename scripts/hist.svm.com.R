#utility functions

#get the chromosome length
get.chrom.length<-function(chr, chrom.info.table){
    return(chrom.info.table[(chrom.info.table[,1]==chr),2])
}

#extract to construct bedgraph and then to bigwig
tobigwig<-function(filename, temp.bg, chromInfo){
    #get bedgraph

    bedgraph.sorted=tempfile()

    options(scipen =99) # not to use scientific notation when writing out

    #write bedgraph formatted dataframes to tempfile
      #write.table(bedgraph,file= bedgraph.file,quote=F,sep="\t",col.names=F,row.names=F)
      command=paste("LC_ALL=C sort -k1,1 -k2,2n", temp.bg, "| uniq >", bedgraph.sorted,sep=" ")
      #browser()
       try(system(command))

      command=paste("bedGraphToBigWig", bedgraph.sorted, chromInfo ,filename, sep=" ")
      #cat(command,"\n")
    try(system(command))
    #unlink(bedgraph.file)
    unlink(bedgraph.sorted)
}

#betools merge
bedTools.merge<-function(bed)
{
  #create temp files
  a.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out

  #write bed formatted dataframes to tempfile
  write.table(bed,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)

  # create the command string and call the command using system()
  command=paste("LC_ALL=C sort -k1,1 -k2,2n",a.file,"| mergeBed -i stdin >",out,sep=" ")
  #cat(command,"\n")
  try(system(command))

  res=read.table(out,header=F)
  unlink(a.file);unlink(out)
  return(res)
}

#function that calls bedtools and operate on two bed dataframes
bedTools.2in<-function(functionstring="bedIntersect",bed1,bed2,opt.string="")
{
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out

  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)

  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  # cat(command,"\n")
  try(system(command))

  res=read.table(out,header=F)
  unlink(a.file);unlink(b.file);unlink(out)
  return(res)
}


#dREG HD peak calling functions
#scan for continuous peak that passes minimun peak length defined
#returns a matrix in the form of matrix(start, end)
scan_for_peak<-function(potential.peak.positions){
    length.tot<-length(potential.peak.positions)
    if (length.tot==0) return(cbind(numeric(0),numeric(0)))
    peak.starts<-c(1)
    peak.ends<-c()
    i<-1
    while (i< length.tot){
        if(potential.peak.positions[i+1]!=potential.peak.positions[i]+1){
        #break point found
        peak.starts<-c(peak.starts,i+1)
        peak.ends<-c(peak.ends,i)
        }
        i<-i+1
    }
    peak.ends<-c(peak.ends ,length.tot)
    return (cbind(potential.peak.positions [peak.starts[!is.na(peak.ends)]], (potential.peak.positions [peak.ends[!is.na(peak.ends)]]+1)))
}

split.bed.evenly<-function(bed, memory=0.1){
    memory <- 1;
    tot.num.examples <- as.integer(memory*1024*1024*1024/8/123)
    line.cutoff<-c(0)
    current.row<-1
    while(current.row<=nrow(bed)){
            current.num.examples <-bed[current.row,3]-bed[current.row,2]
        while(current.num.examples <= tot.num.examples && current.row<=nrow(bed)){
            current.row<-current.row+1
            current.num.examples<-current.num.examples+(bed[current.row,3]-bed[current.row,2])

        }
        line.cutoff<-c(line.cutoff,current.row-1)
    }

    return(line.cutoff)
}

get.chromosome.info <- function(file.plus, file.minus)
{
    bw.plus<-load.bigWig(file.plus);
    bw.minus<-load.bigWig(file.minus);

    chrom <- rbind( cbind( bw.plus$chroms, bw.plus$chromSizes), cbind( bw.minus$chroms, bw.minus$chromSizes) );
    chr.size <- unlist( lapply( unique(chrom[,1]), function(chr){max( as.numeric( chrom[which(chrom[,1]==chr),2])) } ) );

    unload.bigWig(bw.plus);
    unload.bigWig(bw.minus);

    return(data.frame( V1=unique(chrom[,1]), V2=chr.size ));
}

bedTools.intersect<-function(bedA, bedB, options="")
{
  options(scipen =99) # not to use scientific notation when writing out

  #create temp files
  if (is.data.frame(bedA))
  {
    file.a=tempfile(fileext=".bed")
    write.table(bedA, file = file.a, quote=F, sep="\t", col.names=F,row.names=F)
  }
  else
    file.a <- bedA;

  if (is.data.frame(bedB))
  {
    file.b=tempfile(fileext=".bed")
    write.table(bedB, file = file.b, quote=F, sep="\t", col.names=F,row.names=F)
  }
  else
    file.b <- bedB;

  res = read.table(pipe(paste("bedtools intersect -a", file.a, "-b", file.b, options )))

  if(is.data.frame(bedA)) unlink(file.a);
  if(is.data.frame(bedB)) unlink(file.b);

  return(res)
}

bedTools.subtract<-function(bedA, bedB)
{
  options(scipen =99) # not to use scientific notation when writing out

  #create temp files
  if (is.data.frame(bedA))
  {
    file.a=tempfile(fileext=".bed")
    write.table(bedA, file = file.a, quote=F, sep="\t", col.names=F,row.names=F)
  }
  else
    file.a <- bedA;

  if (is.data.frame(bedB))
  {
    file.b=tempfile(fileext=".bed")
    write.table(bedB, file = file.b, quote=F, sep="\t", col.names=F,row.names=F)
  }
  else
    file.b <- bedB;

  res = read.table(pipe(paste("bedtools subtract -a", file.a, "-b", file.b )))

  if(is.data.frame(bedA)) unlink(file.a);
  if(is.data.frame(bedB)) unlink(file.b);

  return(res)
}

split.bed <- function(bedA, interval=10)
{
    bedA <- data.frame(bedA, IDX=1:NROW(bedA));
    max.dist <- max(bedA[,3] - bedA[,2]);

    ret <- bedA[,c(1,2,4)];
    for(i in 1:ceiling(max.dist/interval))
    {
        tmp <- bedA[c(1,2,4)];
        tmp[,2] <- tmp[,2] + interval*i;
        tmp <- tmp[ tmp[,2] <= bedA[,3],];
        ret <- rbind(ret, tmp);
    }

    ret <- data.frame(ret[,c(1,2)], ret[,2]+1, ret[,3]);
    ret <- ret[order(ret[,4], ret[,3]),]

    return(ret[,c(1:3)]);
}


get_histone_read<-function(bedA, file.histone, block=500000, ncores=5)
{
   na.idx <- which( is.na(bedA[,1]) |  is.na(bedA[,2]) |  is.na(bedA[,3])  )
   if(length(na.idx)>0)
      stop("NA in bedA\n");

   cpu.fun <- function(idx) {

       library(bigWig);

       bed <- read.table(file.temp.beds[idx]);
       bw.hist <- load.bigWig(file.histone);
       y0 <- try( unlist(bed.region.bpQuery.bigWig( bw=bw.hist, bed= bed ) ) );
       if(class(y0)=="try-error")
       {
          browser();
          return(NULL);
       }
       unload.bigWig(bw.hist);
       unlink(file.temp.beds[idx]);
       
       return(y0);
   }
    
   file.temp.beds <- unlist(lapply(1:ceiling(NROW(bedA)/block), function(idx) {
       start <- 1+(idx-1)*block;
       end <- if(idx*block>NROW(bedA)) NROW(bedA) else idx*block;
       return( write.temp.bed(bedA[start:end, 1:3], compress=FALSE) );
   }));

   if(ncores>1)
   {
       library(snowfall);
       sfInit(parallel = TRUE, cpus = ncores, type = "SOCK" )
       sfExport("file.histone", "file.temp.beds", "block"  );

       fun <- as.function(cpu.fun);
       environment(fun)<-globalenv();

       y <- unlist (sfLapply(1:NROW(file.temp.beds), cpu.fun) )
       sfStop();
   }   
   else
       y <- unlist( lapply(1:ceiling(NROW(bedA)/block), cpu.fun ) );
   
   return(y);
}

bed.unintersect<-function( bedorfile1, bedorfile2 )
{
    bedfile1 <- bedorfile1;
    bedfile2 <- bedorfile2;

    if(class(bedorfile1)=="data.frame")
    {
        bedfile1 = tempfile();
        write.table(bedorfile1, file=bedfile1, quote=F, col.names=F, row.names=F, sep="\t");
    }

    if(class(bedorfile2)=="data.frame")
    {
        bedfile2 = tempfile();
        write.table(bedorfile2, file=bedfile2, quote=F, col.names=F, row.names=F, sep="\t");
    }

    #Only report those entries in A that have no overlap in B. Restricted by -f and -r.
    tbx <- NULL;
    try( tbx <- read.table(pipe(paste("bedtools intersect -a ", bedfile1, "-b ", bedfile2, " -v"))));

    if(class(bedorfile1)=="data.frame") unlink(bedfile1);
    if(class(bedorfile2)=="data.frame") unlink(bedfile2);

    return(tbx);
}

bed.intersect<-function( bedorfile1, bedorfile2 )
{
    bedfile1 <- bedorfile1;
    bedfile2 <- bedorfile2;

    if(class(bedorfile1)=="data.frame")
    {
        bedfile1 = tempfile();
        write.table(bedorfile1, file=bedfile1, quote=F, col.names=F, row.names=F, sep="\t");
    }

    if(class(bedorfile2)=="data.frame")
    {
        bedfile2 = tempfile();
        write.table(bedorfile2, file=bedfile2, quote=F, col.names=F, row.names=F, sep="\t");
    }


    tbx <- NULL;
    #Only report those entries in A that have no overlap in B. Restricted by -f and -r.
    try(tbx <- read.table(pipe(paste("bedtools intersect -a ", bedfile1, "-b ", bedfile2 ))));

    if(class(bedorfile1)=="data.frame") unlink(bedfile1);
    if(class(bedorfile2)=="data.frame") unlink(bedfile2);

    return(tbx);
}

write.bed<-function ( df.bed, file.bed, compress=FALSE )
{
    options("scipen"=100, "digits"=4);
    temp <- tempfile(fileext=".bed");
    write.table( df.bed, file=temp, quote=F, row.names=F, col.names=F, sep="\t");
    if(compress)
        system(paste0("sort-bed ", temp,  " | bgzip > ",  file.bed ))
    else
        system(paste0("sort-bed ", temp,  " > ",  file.bed ));
        
    invisible(unlink(temp));
}

write.temp.bed<-function ( df.bed, compress=FALSE )
{
    file.bed <- tempfile(fileext=".bed");

    options("scipen"=100, "digits"=4);
    temp <- tempfile(fileext=".bed");
    write.table( df.bed, file=temp, quote=F, row.names=F, col.names=F, sep="\t");
    if(compress)
        system(paste0("sort-bed ", temp,  " | bgzip > ",  file.bed ))
    else
        system(paste0("sort-bed ", temp,  " > ",  file.bed ));
        
    unlink(temp);

    return(file.bed);
}

extract_feature_matrix <- function( infp.bed, file.bw.plus, file.bw.minus, gdm, linear_scale=F, ncores=1, graininess=25 )
{
    zoom<- list(as.integer(gdm@window_sizes), as.integer(gdm@half_nWindows))

    bw.plus <- load.bigWig(file.bw.plus);
    bw.minus <-  load.bigWig(file.bw.minus);
    total <- abs(bw.plus$mean*bw.plus$basesCovered)+abs(bw.minus$mean*bw.minus$basesCovered);
    unload.bigWig(bw.plus);
    unload.bigWig(bw.minus);

    cpu.fun<-function(idx)
    {
        requireNamespace("dREG");
        require(bigWig)

        sect <- ceiling( NROW(infp.bed)/(ncores*graininess) );
        bed.idx <- c( (sect*(idx-1)+1) : (idx*sect) );
        if( min(bed.idx) > NROW(infp.bed) )
            return(NULL);

        if( max(bed.idx) > NROW(infp.bed) )
            bed.idx <- c(min(bed.idx): NROW(infp.bed) );

        #scaling using total read depth
        dat <- .Call("get_genomic_data_R",
                as.character(infp.bed[bed.idx,1]),
                as.integer(infp.bed[bed.idx,2]),
                as.character(file.bw.plus),
                as.character(file.bw.minus),
                zoom,
                as.logical( ifelse(linear_scale, FALSE, TRUE) ),
                PACKAGE= "dREG");

        if (linear_scale)
            dat <- unlist (dat)/(total/1E6)

        dreg_fet <- data.frame(chr=infp.bed[bed.idx,1], start=as.integer(infp.bed[bed.idx,2]), end=as.integer(infp.bed[bed.idx,2]+1), t(dat));
        return(dreg_fet);
    }

    sfInit(parallel = TRUE, cpus = ncores, type = "SOCK" )
    sfExport("infp.bed","zoom","file.bw.plus","file.bw.minus","total", "ncores", "linear_scale", "graininess");

    fun <- as.function(cpu.fun);
    environment(fun)<-globalenv();

    dat <- do.call( rbind.data.frame, sfClusterApplyLB( 1:(ncores*graininess), fun= cpu.fun))
    sfStop();

    #dat <- rbind.data.frame( lapply( 1:(ncores*5), cpu.fun));

    return(list(pos=dat[,c(1:3)], mat=dat[,-c(1:3)]));
}

get_chrom_info <- function( file.bw.plus )
{
    bw.plus <- load.bigWig(file.bw.plus);

    chrom <-  cbind( bw.plus$chroms, bw.plus$chromSizes)
    chr.size <- unlist( lapply( unique(chrom[,1]), function(chr){max( as.numeric( chrom[which(chrom[,1]==chr),2])) } ) );

    df.bed <- data.frame( V1=unique(chrom[,1]), V2=chr.size );
    return(df.bed);
}
