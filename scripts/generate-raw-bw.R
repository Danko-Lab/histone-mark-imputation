library(parallel);
library(data.table);
library(bigWig);

source("../scripts/hist.svm.com.R");
source("../scripts/hist.param.R");

convert_bed_to_bigwig<-function( input.path, prefix, output.path, q.threshold=NULL )
{
	file.beds <- list.files(path=input.path, pattern=glob2rx(paste0(prefix, "*.bed.gz")), full.names=T)

    if(!is.null(q.threshold))
    {
	    for(i in 1:NROW(file.beds))
	    {
	        tb <- read.table(file.beds[i]);
	        print(summary(tb));
	        tb[ tb[,4] <= quantile(tb[,4], q.threshold)*1.05, 4] <- 0;
	        f.temp <- write.temp.bed(tb, compress=TRUE);
	        file.beds[i] <- f.temp;
	        rm(tb);
	    }
    }

	file.tmp <- tempfile(fileext=".bed");
	system( paste("zcat", paste(file.beds, collapse="  " ), " > ", file.tmp ) );
	system( paste("bedGraphToBigWig", file.tmp, file.chromo.info, paste0(output.path, "/raw.", prefix, ".bw")));
	unlink(file.tmp);
	
}


if(0)
{
convert_bed_to_bigwig("../pred-gm/beds/", "H3k27ac.S1.V3.GM", "../pred-gm/");
convert_bed_to_bigwig("../pred-gm/beds/", "H3k27me3.S1.V3.GM", "../pred-gm/");
convert_bed_to_bigwig("../pred-gm/beds/", "H3k36me3.S1.V3.GM", "../pred-gm/");
convert_bed_to_bigwig("../pred-gm/beds/", "H3k4me1.S1.V2.GM",  "../pred-gm/");
convert_bed_to_bigwig("../pred-gm/beds/", "H3k4me2.S1.V2.GM",  "../pred-gm/");
convert_bed_to_bigwig("../pred-gm/beds/", "H3k4me3.S1.V3.GM",  "../pred-gm/");
convert_bed_to_bigwig("../pred-gm/beds/", "H3k9ac.S1.V2.GM",   "../pred-gm/");
convert_bed_to_bigwig("../pred-gm/beds/", "H3k9me3.S1.V2.GM",  "../pred-gm/");
convert_bed_to_bigwig("../pred-gm/beds/", "H4k20me1.S1.V3.GM", "../pred-gm/");
}

if(0)
{
convert_bed_to_bigwig("../pred-k562/beds/", "H3k27ac.S1.V3.G1", "../pred-k562/");
convert_bed_to_bigwig("../pred-k562/beds/", "H3k122ac.S1.V3.G1", "../pred-k562/");
convert_bed_to_bigwig("../pred-k562/beds/", "H3k27me3.S1.V3.G1", "../pred-k562/");
convert_bed_to_bigwig("../pred-k562/beds/", "H3k36me3.S1.V3.G1", "../pred-k562/");
convert_bed_to_bigwig("../pred-k562/beds/", "H3k4me1.S1.V2.G1",  "../pred-k562/");
convert_bed_to_bigwig("../pred-k562/beds/", "H3k4me2.S1.V2.G1",  "../pred-k562/");
convert_bed_to_bigwig("../pred-k562/beds/", "H3k4me3.S1.V3.G1",  "../pred-k562/");
convert_bed_to_bigwig("../pred-k562/beds/", "H3k9ac.S1.V2.G1",   "../pred-k562/");
convert_bed_to_bigwig("../pred-k562/beds/", "H3k9me3.S1.V2.G1",  "../pred-k562/");
convert_bed_to_bigwig("../pred-k562/beds/", "H4k20me1.S1.V3.G1", "../pred-k562/");
}

if(0)
{
convert_bed_to_bigwig("../pred-hct/beds/", "H3k122ac.S1.V3.HCT", "../pred-hct/");
convert_bed_to_bigwig("../pred-hct/beds/", "H3k27ac.S1.V3.HCT",  "../pred-hct/");
convert_bed_to_bigwig("../pred-hct/beds/", "H3k27me3.S1.V3.HCT", "../pred-hct/");
convert_bed_to_bigwig("../pred-hct/beds/", "H3k36me3.S1.V3.HCT", "../pred-hct/");
convert_bed_to_bigwig("../pred-hct/beds/", "H3k4me1.S1.V2.HCT",  "../pred-hct/");
convert_bed_to_bigwig("../pred-hct/beds/", "H3k4me2.S1.V2.HCT",  "../pred-hct/");
convert_bed_to_bigwig("../pred-hct/beds/", "H3k4me3.S1.V3.HCT",  "../pred-hct/");
convert_bed_to_bigwig("../pred-hct/beds/", "H3k9ac.S1.V2.HCT",   "../pred-hct/");
convert_bed_to_bigwig("../pred-hct/beds/", "H3k9me3.S1.V2.HCT",  "../pred-hct/");
convert_bed_to_bigwig("../pred-hct/beds/", "H4k20me1.S1.V3.HCT", "../pred-hct/");
}


if(0)
{ 
convert_bed_to_bigwig("../pred-cd4/beds/", "H3k122ac.S1.V3.CD4", "../pred-cd4/");
convert_bed_to_bigwig("../pred-cd4/beds/", "H3k27ac.S1.V3.CD4",  "../pred-cd4/");
convert_bed_to_bigwig("../pred-cd4/beds/", "H3k27me3.S1.V3.CD4", "../pred-cd4/");
convert_bed_to_bigwig("../pred-cd4/beds/", "H3k36me3.S1.V3.CD4", "../pred-cd4/");
convert_bed_to_bigwig("../pred-cd4/beds/", "H3k4me1.S1.V2.CD4",  "../pred-cd4/");
convert_bed_to_bigwig("../pred-cd4/beds/", "H3k4me2.S1.V2.CD4",  "../pred-cd4/");
convert_bed_to_bigwig("../pred-cd4/beds/", "H3k4me3.S1.V3.CD4",  "../pred-cd4/");
convert_bed_to_bigwig("../pred-cd4/beds/", "H3k9ac.S1.V2.CD4",   "../pred-cd4/");
convert_bed_to_bigwig("../pred-cd4/beds/", "H3k9me3.S1.V2.CD4",  "../pred-cd4/");
convert_bed_to_bigwig("../pred-cd4/beds/", "H4k20me1.S1.V3.CD4", "../pred-cd4/");
}


if(0)
{
convert_bed_to_bigwig("../pred-hela/beds/", "H3k27ac.S1.V3.HELA",  "../pred-hela/");
convert_bed_to_bigwig("../pred-hela/beds/", "H3k27me3.S1.V3.HELA", "../pred-hela/");
convert_bed_to_bigwig("../pred-hela/beds/", "H3k36me3.S1.V3.HELA", "../pred-hela/");
convert_bed_to_bigwig("../pred-hela/beds/", "H3k4me1.S1.V2.HELA",  "../pred-hela/");
convert_bed_to_bigwig("../pred-hela/beds/", "H3k4me2.S1.V2.HELA",  "../pred-hela/");
convert_bed_to_bigwig("../pred-hela/beds/", "H3k4me3.S1.V3.HELA",  "../pred-hela/");
convert_bed_to_bigwig("../pred-hela/beds/", "H3k9ac.S1.V2.HELA",   "../pred-hela/");
convert_bed_to_bigwig("../pred-hela/beds/", "H3k9me3.S1.V2.HELA",  "../pred-hela/");
convert_bed_to_bigwig("../pred-hela/beds/", "H4k20me1.S1.V3.HELA", "../pred-hela/");
}


if(0)
{
convert_bed_to_bigwig("/workdir/zw355/proj/prj15-histone/horse_don/pred_by_k562_model", "H3k27ac.S1.V3.G1", "/workdir/zw355/proj/prj15-histone/horse_don/");
convert_bed_to_bigwig("/workdir/zw355/proj/prj15-histone/horse_don/pred_by_k562_model", "H3k27me3.S1.V3.G1", "/workdir/zw355/proj/prj15-histone/horse_don/");
convert_bed_to_bigwig("/workdir/zw355/proj/prj15-histone/horse_don/pred_by_k562_model", "H3k4me1.S1.V2.G1", "/workdir/zw355/proj/prj15-histone/horse_don/");
convert_bed_to_bigwig("/workdir/zw355/proj/prj15-histone/horse_don/pred_by_k562_model", "H3k4me3.S1.V3.G1", "/workdir/zw355/proj/prj15-histone/horse_don/");
}


if(0)
{
convert_bed_to_bigwig("../pred-alex/", "Alex.H3K27ac.S1", "../pred-alex/");
convert_bed_to_bigwig("../pred-alex/", "Alex.H3k4me1.S1", "../pred-alex/");
convert_bed_to_bigwig("../pred-alex/", "Alex.H3k4me2.S1", "../pred-alex/");
convert_bed_to_bigwig("../pred-alex/", "Alex.H3k4me3.S1", "../pred-alex/");

convert_bed_to_bigwig("../pred-alex/", "Alex2.Mnase.H3k27ac.S1", "../pred-alex/");
convert_bed_to_bigwig("../pred-alex/", "Alex.Mnase.H3k4me1.S1", "../pred-alex/");
convert_bed_to_bigwig("../pred-alex/", "Alex.Mnase.H3k4me2.S1", "../pred-alex/");
convert_bed_to_bigwig("../pred-alex/", "Alex.Mnase.H3k4me3.S1", "../pred-alex/");
}

if(0)
{
#convert_bed_to_bigwig("../atac/", "ATAC.K562.S1", "../atac/");
convert_bed_to_bigwig("../atac/", "ATAC.exlchr22.K562.S1", "../atac/");
convert_bed_to_bigwig("../dnase/", "dNase.exlchr22.K562.S1", "../dnase/");

}


if(0)
{
convert_bed_to_bigwig("../pred-k562/", "H3K27ac.S3.V2", "../pred-k562/");
convert_bed_to_bigwig("../pred-k562/", "H3K27ac.S4.w10.V1", "../pred-k562/");
convert_bed_to_bigwig("../pred-k562/", "H3K27ac.S4.w50.V1", "../pred-k562/");
}

