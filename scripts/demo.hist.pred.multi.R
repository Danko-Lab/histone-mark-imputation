source("hist.svm.com.R")
source("hist.svm.pred.R")
source("hist.svm.main.R")
source("hist.svm.paralle.pred.R")
source("hist.param.R")

get_model_file <- function(str.hist)
{
	vec.hist <- c(
		"h3k122ac",  
		"h3k27ac",   
		"h3k27me3",  
		"h3k36me3", 
		"h3k4me1",  
		"h3k4me2",  
		"h3k4me3",  
		"h3k9ac",   
		"h3k9me3",  
		"h4k20me1");
	
	file.rdata <- c(
		file.H3k122ac.model,  
		file.H3k27ac.model,   
		file.H3k27me3.model,  
		file.H3k36me3.model, 
		file.H3k4me1.model,  
		file.H3k4me2.model,  
		file.H3k4me3.model,  
		file.H3k9ac.model,   
		file.H3k9me3.model,  
		file.H4k20me1.model);
		
	idx <- which(tolower(str.hist)==vec.hist);
	return(file.rdata[idx]);
}

args = commandArgs(trailingOnly=TRUE)
str.histone.type <- args[1];
file.bigwig.list <- args[2];
file.bed.region<- args[3];
str.prefix <- args[4];
gpu.count <- if(is.na(args[5])) 0 else as.numeric(args[5])

cat("Histone type = ", str.histone.type, "\n");
cat("bigWig list = ", file.bigwig.list, "\n");
cat("file of bed region = ", file.bed.region, "\n");
cat("Prefix = ", str.prefix, "\n");
cat("GPU count = ", gpu.count, "\n");

if(!file.exists(file.bigwig.list)) stop("No bigwig file!\n");
if(!file.exists(file.bed.region)) stop("No bed regions!\n");

file.model.rdata <- get_model_file( str.histone.type );
if(!file.exists(file.model.rdata)) stop("No model file!\n");

library(Rgtsvm)

tb.bigwig <- read.table(file.bigwig.list);
tb.bed.region <- read.table(file.bed.region);

svm_predict_multfile_parallel( str.prefix, file.model.rdata,  as.character(tb.bigwig[,1]), as.character(tb.bigwig[,2]), tb.bed.region, ncores=4, linear_scale=F, gpu.count=gpu.count, str.wd="../GBM20/" )

