library(parallel);
library(data.table);
library(bigWig);

source("../scripts/hist.svm.com.R");
source("../scripts/hist.param.R");

comvert_bed_to_bigwig<-function( input.path, prefix, output.path, rate=1 )
{
	file.beds <- list.files(path=input.path, pattern=glob2rx(paste0(prefix, "*.bed.gz")), full.names=T)

	file.tmp.beds <- unlist( mclapply(file.beds, function(file.bed)
	{
       tb <- read.table(file.bed);
       thres <- median(round(tb[,4],2));
cat(thres, file.bed, "\n")       
       tb[,4] <- tb[,4] - thres;
       tb <- tb[tb[,4]>0,];
	   file.tmp <- write.temp.bed(tb, compress=TRUE);
	   return(file.tmp);
	}, mc.cores=12, mc.preschedule = F));

	file.tmp <- tempfile(fileext=".bed");
	system( paste("zcat", paste(file.tmp.beds, collapse="  " ), " > ", file.tmp ) );
	system( paste("bedGraphToBigWig", file.tmp, file.chromo.info, paste0(output.path, "/raw.", prefix, ".bw")));
	unlink(file.tmp);
	unlink(file.tmp.beds);
}

if(1)
{
comvert_bed_to_bigwig("../pred-k562/beds/", "H3k122ac.S1.V3.G1", "../pred-k562/");
comvert_bed_to_bigwig("../pred-k562/beds/", "H3k27ac.S1.V3.G1",  "../pred-k562/");
comvert_bed_to_bigwig("../pred-k562/beds/", "H3k27me3.S1.V3.G1", "../pred-k562/");
comvert_bed_to_bigwig("../pred-k562/beds/", "H3k36me3.S1.V3.G1", "../pred-k562/");
comvert_bed_to_bigwig("../pred-k562/beds/", "H3k4me1.S1.V2.G1",  "../pred-k562/");
comvert_bed_to_bigwig("../pred-k562/beds/", "H3k4me2.S1.V2.G1",  "../pred-k562/");
comvert_bed_to_bigwig("../pred-k562/beds/", "H3k4me3.S1.V3.G1",  "../pred-k562/");
comvert_bed_to_bigwig("../pred-k562/beds/", "H3k9ac.S1.V2.G1",   "../pred-k562/");
comvert_bed_to_bigwig("../pred-k562/beds/", "H3k9me3.S1.V2.G1",  "../pred-k562/");
comvert_bed_to_bigwig("../pred-k562/beds/", "H4k20me1.S1.V3.G1", "../pred-k562/");
}

if(1)
{
comvert_bed_to_bigwig("../pred-gm/beds/", "H3k27ac.S1.V3.GM", "../pred-gm/");
comvert_bed_to_bigwig("../pred-gm/beds/", "H3k27me3.S1.V3.GM", "../pred-gm/");
comvert_bed_to_bigwig("../pred-gm/beds/", "H3k36me3.S1.V3.GM", "../pred-gm/");
comvert_bed_to_bigwig("../pred-gm/beds/", "H3k4me1.S1.V2.GM",  "../pred-gm/");
comvert_bed_to_bigwig("../pred-gm/beds/", "H3k4me2.S1.V2.GM",  "../pred-gm/");
comvert_bed_to_bigwig("../pred-gm/beds/", "H3k4me3.S1.V3.GM",  "../pred-gm/");
comvert_bed_to_bigwig("../pred-gm/beds/", "H3k9ac.S1.V2.GM",   "../pred-gm/");
comvert_bed_to_bigwig("../pred-gm/beds/", "H3k9me3.S1.V2.GM",  "../pred-gm/");
comvert_bed_to_bigwig("../pred-gm/beds/", "H4k20me1.S1.V3.GM", "../pred-gm/");
}


if(1)
{
comvert_bed_to_bigwig("../pred-hct/beds/", "H3k27ac.S1.V3.HCT",  "../pred-hct/");
comvert_bed_to_bigwig("../pred-hct/beds/", "H3k27me3.S1.V3.HCT", "../pred-hct/");
comvert_bed_to_bigwig("../pred-hct/beds/", "H3k36me3.S1.V3.HCT", "../pred-hct/");
comvert_bed_to_bigwig("../pred-hct/beds/", "H3k4me1.S1.V2.HCT",  "../pred-hct/");
comvert_bed_to_bigwig("../pred-hct/beds/", "H3k4me2.S1.V2.HCT",  "../pred-hct/");
comvert_bed_to_bigwig("../pred-hct/beds/", "H3k4me3.S1.V3.HCT",  "../pred-hct/");
comvert_bed_to_bigwig("../pred-hct/beds/", "H3k9ac.S1.V2.HCT",   "../pred-hct/");
comvert_bed_to_bigwig("../pred-hct/beds/", "H3k9me3.S1.V2.HCT",  "../pred-hct/");
comvert_bed_to_bigwig("../pred-hct/beds/", "H4k20me1.S1.V3.HCT", "../pred-hct/");
}


if(1)
{ 
comvert_bed_to_bigwig("../pred-cd4/beds/", "H3k27ac.S1.V3.CD4",  "../pred-cd4/");
comvert_bed_to_bigwig("../pred-cd4/beds/", "H3k27me3.S1.V3.CD4", "../pred-cd4/");
comvert_bed_to_bigwig("../pred-cd4/beds/", "H3k36me3.S1.V3.CD4", "../pred-cd4/");
comvert_bed_to_bigwig("../pred-cd4/beds/", "H3k4me1.S1.V2.CD4",  "../pred-cd4/");
comvert_bed_to_bigwig("../pred-cd4/beds/", "H3k4me2.S1.V2.CD4",  "../pred-cd4/");
comvert_bed_to_bigwig("../pred-cd4/beds/", "H3k4me3.S1.V3.CD4",  "../pred-cd4/");
comvert_bed_to_bigwig("../pred-cd4/beds/", "H3k9ac.S1.V2.CD4",   "../pred-cd4/");
comvert_bed_to_bigwig("../pred-cd4/beds/", "H3k9me3.S1.V2.CD4",  "../pred-cd4/");
comvert_bed_to_bigwig("../pred-cd4/beds/", "H4k20me1.S1.V3.CD4", "../pred-cd4/");
}


if(1)
{
comvert_bed_to_bigwig("../pred-hela/beds/", "H3k27ac.S1.V3.HELA",  "../pred-hela/");
comvert_bed_to_bigwig("../pred-hela/beds/", "H3k27me3.S1.V3.HELA", "../pred-hela/");
comvert_bed_to_bigwig("../pred-hela/beds/", "H3k36me3.S1.V3.HELA", "../pred-hela/");
comvert_bed_to_bigwig("../pred-hela/beds/", "H3k4me1.S1.V2.HELA",  "../pred-hela/");
comvert_bed_to_bigwig("../pred-hela/beds/", "H3k4me2.S1.V2.HELA",  "../pred-hela/");
comvert_bed_to_bigwig("../pred-hela/beds/", "H3k4me3.S1.V3.HELA",  "../pred-hela/");
comvert_bed_to_bigwig("../pred-hela/beds/", "H3k9ac.S1.V2.HELA",   "../pred-hela/");
comvert_bed_to_bigwig("../pred-hela/beds/", "H3k9me3.S1.V2.HELA",  "../pred-hela/");
comvert_bed_to_bigwig("../pred-hela/beds/", "H4k20me1.S1.V3.HELA", "../pred-hela/");
}


if(0)
{
comvert_bed_to_bigwig("../horse_don/pred_by_k562_model", "H3k27ac.S1.V3.G1",  "../horse_don/");
comvert_bed_to_bigwig("../horse_don/pred_by_k562_model", "H3k27me3.S1.V3.G1", "../horse_don/");
comvert_bed_to_bigwig("../horse_don/pred_by_k562_model", "H3k4me1.S1.V2.G1",  "../horse_don/");
comvert_bed_to_bigwig("../horse_don/pred_by_k562_model", "H3k4me3.S1.V3.G1",  "../horse_don/");
}


if(0)
{
comvert_bed_to_bigwig("../pred-alex/", "Alex.H3K27ac.Gx.S1", "../pred-alex/");
comvert_bed_to_bigwig("../pred-alex/", "Alex.H3k4me1.Gx.S1", "../pred-alex/");
comvert_bed_to_bigwig("../pred-alex/", "Alex.H3k4me2.Gx.S1", "../pred-alex/");
comvert_bed_to_bigwig("../pred-alex/", "Alex.H3k4me3.Gx.S1", "../pred-alex/");

comvert_bed_to_bigwig("../pred-alex/", "Alex.Mnase.H3k4me1.Gx.S1", "../pred-alex/");
comvert_bed_to_bigwig("../pred-alex/", "Alex.Mnase.H3k4me2.Gx.S1", "../pred-alex/");
comvert_bed_to_bigwig("../pred-alex/", "Alex.Mnase.H3k4me3.Gx.S1", "../pred-alex/");
comvert_bed_to_bigwig("../pred-alex/", "Alex2.Mnase.H3k27ac.Gx.S1", "../pred-alex/");
comvert_bed_to_bigwig("../pred-alex/", "Alex2.Mnase.H3k36me3.Gx.S1", "../pred-alex/");
comvert_bed_to_bigwig("../pred-alex/", "Alex2.Mnase.H3k79me3.Gx.S1", "../pred-alex/");
}

if(0)
{
comvert_bed_to_bigwig("../atac/", "ATAC.K562.S1", "../atac/");
comvert_bed_to_bigwig("../atac/", "ATAC.exlchr22.K562.S1", "../atac/");
comvert_bed_to_bigwig("../dnase/", "dNase.exlchr22.K562.S1", "../dnase/", rate=2000);
}

