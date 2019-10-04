args <- commandArgs(trailingOnly=TRUE);

file_bed_gz  <- args[1]
file_output <- args[2]

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

tb <- read.table(file_bed_gz, stringsAsFactors=TRUE);

df <- data.frame(tb[-NROW(tb),1:3], tb[-1, 1:3]);

rem.idx <- which(df[,1]==df[,4] & df[,3]>df[,5])
if(NROW(rem.idx)==0)
 cat("No overlapped regions, quit.\n");

if(NROW(rem.idx)>0)
{
  tb <- tb[-rem.idx,]
  library(tools); 
  write.bed(tb, file_output, compress=if (file_ext(file_output)=="gz") TRUE else FALSE);  
}

  