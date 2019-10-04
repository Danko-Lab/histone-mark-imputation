source("../hist.svm.com.R")
source("../hist.svm.pred.R")
source("../hist.svm.main.R")
source("../hist.param.R");

bw.minus  <- list.files(pattern="*minus.bw");
df.bws <- data.frame(stringsAsFactors=F)
for(fn in bw.minus)
{
        r <- regexpr( "_minus.bw", fn)
        fid <- substr(fn, 1, r-1);
cat(fid, "\n");
        df.bws <- rbind(df.bws, data.frame(FID=fid, minus=fn, plus=paste(fid, "_plus.bw", sep=""), stringsAsFactors=F));
}


file.bed <- "all.merged.bdRM.centered1000.bed";
tb <- split.bed(read.table(file.bed, header=F), 10)
file.tmp.bed <- tempfile(fileext=".bed");
write.table(tb, file=file.tmp.bed, quote=F, row.names=F, col.names=F, sep="\t");

file.model <- "../H3K27ac-model.rdata";
file.histone <- "../histdata/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig"

#load model.hd, gdm.hd, DREG_THRES
load( file.model );
model <- svm;
rm(svm);

for(i in 1:NROW(df.bws))
{
        if(!file.exists(paste(df.bws[i,1], "_H3k27ac.bed", sep="")))
        {
                cat("Input File=", df.bws[i,1], df.bws[i,2], df.bws[i,3], "\n");

                pred <- svm_predict( paste(df.bws[i,1],"_tmp", sep=""),  df.bws[i,3], df.bws[i,2], model, gdm, file.bed=file.tmp.bed, bigwig_compare=file.histone,  ncores=15 );

                cat("Output File=", pred$file.bed, "\n");
                cat("COR=", pred$corr, "\n");

                file.rename( pred$file.bed, paste(df.bws[i,1], "_H3k27ac.bed", sep="") );
        }
}
