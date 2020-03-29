source("hist.svm.com.R")
source("hist.svm.pred.R")
source("hist.svm.main.R")
source("hist.param.R");

#path.histone defined in "hist.param.R"
#path.proseq  defined in "hist.param.R"
#file.rdata.negative defined in "hist.param.R"
#file.rdata.positive defined in "hist.param.R"

file.mm9.mESC.proseq.minus <- c("/fs/cbsudanko/storage/data/mm9/esc/groseq/GSE48895_V6.5_untreated_Minus.bw",
                               "/fs/cbsudanko/storage/data/mm9/esc/groseq/GSE27037_MESC_GROseq_minus.bw",
                               "/workdir/zw355/proj/prj15-histone/train_mm9/GSM1422167_adleman_r1_QC_minus.mm9.bw",
                               "/workdir/zw355/proj/prj15-histone/train_mm9/GSM1422168_adleman_r2_QC_minus.mm9.bw");


file.mm9.mESC.proseq.plus  <- c("/fs/cbsudanko/storage/data/mm9/esc/groseq/GSE48895_V6.5_untreated_Plus.bw",
                                "/fs/cbsudanko/storage/data/mm9/esc/groseq/GSE27037_MESC_GROseq_plus.bw",
                                "/workdir/zw355/proj/prj15-histone/train_mm9/GSM1422167_adleman_r1_QC_plus.mm9.bw",
                                "/workdir/zw355/proj/prj15-histone/train_mm9/GSM1422168_adleman_r2_QC_plus.mm9.bw")

file.mm9.mESC.H3k4me3  <- "../mESC-mm9-2018/H3K4me3_WT_merge.bw"
file.mm9.mESC.H3k27me3 <- "../mESC-mm9-2018/H3K27me3_WT_merge.bw"

#file.mm9.mESC.dreg     <- "../pred-mm9-mESC/GSE48895-mm9.dREG.peak.score.bed.gz"
file.mm9.mESC.dreg     <- "../train_mm9/merge_dreg_peak.bed"

file.mm9.mESC.H3k4me3.peak  <- "../mESC-mm9-2018/H3K4me3_WT_merge.bed"
file.mm9.mESC.H3k27me3.peak <- "../mESC-mm9-2018/H3K27me3_WT_merge.bed"

file.H3k4me3.model  <- "../models/mESC.comb4.H3k4me3.S1.train.rdata";
file.H3k27me3.model <- "../models/mESC.comb4.H3k27me3.S1.train.rdata";

file.pos.bed.H3k4me3 <- "../mESC-mm9-2018/ATAC_WT_R.peak.ineresect.bed.gz"
file.neg.bed.H3k4me3 <- "../mESC-mm9-2018/ATAC_WT_R.peak.complement.bed.gz"

file.pos.bed.H3k27me3 <- "../mESC-mm9-2018/ATAC_WT_R.peak.ineresect.bed.gz"
file.neg.bed.H3k27me3 <- "../mESC-mm9-2018/ATAC_WT_R.peak.complement.bed.gz"

path.histone="";
path.proseq="";


#training H3k4me3 Model
if(1)
{
if(!file.exists(file.H3k4me3.model))
{
  tb <- rbind(read.table(file.mm9.mESC.dreg)[,c(1:3)], read.table(file.mm9.mESC.H3k4me3.peak)[,c(1:3)]);
  print(head(tb));
  file.tmp.pos.bed <- write.temp.bed(tb, compress=FALSE);
  print(file.tmp.pos.bed)

  tb <- read.table(pipe(paste("bedtools complement -i ",  file.tmp.pos.bed, " -g /fs/cbsudanko/storage/data/mm9/mm9.chromInfo.sorted | sort-bed -")));
  print(head(tb))
  file.tmp.neg.bed <- write.temp.bed(tb, compress=FALSE);
  print(head(tb))
  print(file.tmp.neg.bed)

  selectGPUdevice(0);

  model <- create_train_model( path.proseq, file.mm9.mESC.proseq.plus, file.mm9.mESC.proseq.minus, 
                file.pos.bed.H3k4me3, file.neg.bed.H3k4me3,
                ## file.tmp.neg.bed, file.mm9.mESC.dreg, 
                path.histone, file.mm9.mESC.H3k4me3, 
                file.mm9.mESC.H3k4me3.peak, ratio = 0.1, 
                samples=750000, strategy=1, 
                exclude=c("chr1", "chrX", "chrY", "chrM"))
  model <- build_train_model( model )
  model <- svm_train_model(model, gdm, file.H3k4me3.model, ncores=8);
  save(model, file = file.H3k4me3.model);
}
if(0){ 
pred1 <- svm_predict_chr_parallel( "mESC.H3k4me3.trainbyself.S1", 
                file.mm9.mESC.proseq.plus, 
                file.mm9.mESC.proseq.minus, 
                file.H3k4me3.model, 
                chr=c("chr1"), 
                bigwig_compare=file.mm9.mESC.H3k4me3, 
                ncores=5, linear_scale=F, gpu.idx=0  );
                
pred2 <- svm_predict_chr_parallel( "mESC.H3k4me3.trainbyself.S1", 
                file.mm9.mESC.proseq.plus, 
                file.mm9.mESC.proseq.minus, 
                file.H3k4me3.model, 
                chr=c("chr2"), 
                bigwig_compare=file.mm9.mESC.H3k4me3, 
                ncores=5, linear_scale=F, gpu.idx=0  );
                
}}

#training H3k27me3 Model
if(0) {
if(!file.exists(file.H3k27me3.model))
{
  tb <- rbind(read.table(file.mm9.mESC.dreg)[,c(1:3)], read.table(file.mm9.mESC.H3k27me3.peak)[,c(1:3)]);
  print(head(tb));
  file.tmp.pos.bed <- write.temp.bed(tb, compress=FALSE);

  tb <- read.table(pipe(paste("bedtools complement -i ",  file.tmp.pos.bed, " -g /fs/cbsudanko/storage/data/mm9/mm9.chromInfo.sorted | sort-bed -")));
  file.tmp.neg.bed <- write.temp.bed(tb, compress=FALSE);

  selectGPUdevice(1);

  model <- create_train_model( path.proseq, file.mm9.mESC.proseq.plus, file.mm9.mESC.proseq.minus, 
                #file.tmp.neg.bed, file.mm9.mESC.dreg, 
                file.pos.bed.H3k27me3, file.neg.bed.H3k27me3,
                path.histone, file.mm9.mESC.H3k27me3, 
                file.mm9.mESC.H3k27me3.peak, ratio = 0.1, 
                samples=400000, strategy=1, 
                exclude=c("chr1", "chrX", "chrY", "chrM"))
  model <- build_train_model( model )
  model <- svm_train_model(model, gdm, file.H3k27me3.model, ncores=8);
  save(model, file = file.H3k27me3.model);
}
 
pred1 <- svm_predict_chr_parallel( "mESC.H3k27me3.trainbyself.S1", 
                file.mm9.mESC.proseq.plus, 
                file.mm9.mESC.proseq.minus, 
                file.H3k27me3.model, 
                chr=c("chr1"), 
                bigwig_compare=file.mm9.mESC.H3k27me3, 
                ncores=5, linear_scale=F, gpu.idx=0  );
                
pred2 <- svm_predict_chr_parallel( "mESC.H3k27me3.trainbyself.S1", 
                file.mm9.mESC.proseq.plus, 
                file.mm9.mESC.proseq.minus, 
                file.H3k27me3.model, 
                chr=c("chr2"), 
                bigwig_compare=file.mm9.mESC.H3k27me3, 
                ncores=5, linear_scale=F, gpu.idx=0  );
}                
