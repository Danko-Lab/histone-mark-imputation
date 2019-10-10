# Computes enrichment around TSS and PAS of specific BED states.
setwd("/local/storage/projects/histone_imputation/chromHMM/")

require(bigWig)

tss <- read.table("/local/storage/projects/krab-enh-deletion/k562_annotations/tuSelecter/final_tus.txt", header=TRUE)
tss <- tss[tss$TXNAME != "Untranscribed",]
tss <- tss[tss$TXTYPE=="protein_coding",]
tss <- tss[tss$TXCHROM!="chrM",]
tss <- tss[tss$TXCHROM!="chrX",]
tss <- tss[(tss$TXEND - tss$TXSTART) > 20000,]
## 10198 protein-coding TSSs where TU > 10kb.

## Polyadenylation cleavage site.
pas <- tss
pas[pas$TXSTRAND=="+",2] <- pas[pas$TXSTRAND=="+",3]-1
pas[pas$TXSTRAND=="-",3] <- pas[pas$TXSTRAND=="-",2]+1

## Gene center.
gcen <- tss
gcen[,2] <- round((gcen$TXEND - gcen$TXSTART)/2) + gcen$TXSTART
gcen[,3] <- gcen[,2]+1

## Take 1 bp for the TSS
tss[tss$TXSTRAND=="+",3] <- tss[tss$TXSTRAND=="+",2]+1
tss[tss$TXSTRAND=="-",2] <- tss[tss$TXSTRAND=="-",3]-1

#tss <- tss[tss$TXSTRAND=="+",]## Limit to plus strand to check assymetry...

metaPlot <- function(bed, PREFIX, path="./", stp= 100, halfWindow = 10000, ...) {
        HP <- load.bigWig(paste(path, PREFIX, ".bw", sep=""))
		H_sum <- bed.step.bpQuery.bigWig(HP, center.bed(bed, halfWindow, halfWindow), step=stp, abs.value=TRUE)
        hmat <- matrix(unlist(H_sum), nrow= NROW(bed), byrow=TRUE)
		
		HR <- load.bigWig(paste(path, PREFIX, ".rand.bw", sep=""))
		H_sum <- bed.step.bpQuery.bigWig(HR, center.bed(bed, halfWindow, halfWindow), step=stp, abs.value=TRUE)
        hrand <- matrix(unlist(H_sum), nrow= NROW(bed), byrow=TRUE)

		## Compute seperately for the plus and minus strand. Have to reverse the minus strand.
		## Note that this logic will work ONLY when upstream and downstream distances are identical.
		enrich <- (colSums(hmat[bed$TXSTRAND == "+",])+rev(colSums(hmat[bed$TXSTRAND == "-",]))) / 
					(colSums(hrand[bed$TXSTRAND == "+",])+rev(colSums(hrand[bed$TXSTRAND == "-",])))

        N = length(enrich)
        x = 1:N*stp - stp/2 - halfWindow

        #plot(enrich~x, ...)
		return(enrich)
}

stp = 100; halfWindow = 10000; path="./"

org_tss<- lapply(1:18, function(x) { metaPlot(tss, paste("MS_K562_org_seg.state.", x, sep=""), path, stp, halfWindow) })
enew_tss<-lapply(1:18, function(x) { metaPlot(tss, paste("MS_K562_UWother_seg.state.", x, sep=""), path, stp, halfWindow) })
pred_tss<-lapply(1:18, function(x) { metaPlot(tss, paste("MS_K562_pred_seg.state.",x, sep=""), path, stp, halfWindow) })

# org_gcen<- lapply(1:18, function(x) { metaPlot(gcen, paste("MS_K562_org_seg.state.", x, sep=""), path, stp, halfWindow) })
# pred_gcen<-lapply(1:18, function(x) { metaPlot(gcen, paste("MS_K562_pred_seg.state.",x, sep=""), path, stp, halfWindow) })

# org_pas<- lapply(1:18, function(x) { metaPlot(pas, paste("MS_K562_org_seg.state.", x, sep=""), path, stp, halfWindow) })
# pred_pas<-lapply(1:18, function(x) { metaPlot(pas, paste("MS_K562_pred_seg.state.",x, sep=""), path, stp, halfWindow) })

cols=c("#FF0000", "#FF4500", "#FF4500", "#FF4500", "#008000", "#006400", "#006400", "#006400", "#FFC34D",  "#FFC34D", "#FFFF00", "#65CDAA", "#8A91D0", "#CD5C5C", "#BDB76B", "#808080", "#C0C0C0", "#DDDDDD")
labs=c("TssA", "TssFlnk", "TssFlnkU", "TssFlnkD", "Tx", "TxWk", "EnhG1", "EnhG2", "EnhA1",  "EnhA2", "EnhWk", "Znf/Rpts", "Het", "TssBiv", "EnhBiv", "ReprPC", "ReprPCWk", "Quies")

setwd("~/transfer/")

pdf("StateCompare.pdf", width=20/2, height=25/2)

par(mfrow=c(4, 5))

for(i in 1:18) {
 max_y <- max(c(org_tss[[i]], enew_tss[[i]], pred_tss[[i]])); if(max_y == Inf) {max_y = 15}
 plot(-100000,-1, xlim=c(min(x), max(x)), ylim=c(0, max_y), main= labs[i], xlab= "Distance from gene center [bp]", ylab= "Fold enrichment")
 lines(org_tss[[i]]~x, col=cols[i], lwd=0.5)
 lines(enew_tss[[i]]~x, col=cols[i], lwd=0.5, lty="dotted")
 lines(pred_tss[[i]]~x, col=cols[i], lwd=3)
}

dev.off()

pdf("MetaPlots.pdf", width=5, height=5)

par(mfrow=c(1, 1))

x = 1:length(org_tss[[1]])*stp - stp/2 - halfWindow
plot(-100000,-1, xlim=c(min(x), max(x)), ylim=c(0, 150), main="Experimental", xlab= "Distance from TSS [bp]", ylab= "Fold enrichment")
for(i in 1:18) {
 lines(org_tss[[i]]~x, col=cols[i])
}

plot(-100000,-1, xlim=c(min(x), max(x)), ylim=c(0, 150), main="Experimental (alternate datasets)", xlab= "Distance from TSS [bp]", ylab= "Fold enrichment")
for(i in 1:18) {
 lines(enew_tss[[i]]~x, col=cols[i])
}

plot(-100000,-1, xlim=c(min(x), max(x)), ylim=c(0, 150), main="Predicted", xlab= "Distance from TSS [bp]", ylab= "Fold enrichment")
for(i in 1:18) {
 lines(pred_tss[[i]]~x, col=cols[i])
}

# plot(-100000,-1, xlim=c(min(x), max(x)), ylim=c(0, 10), main="Experimental", xlab= "Distance from gene center [bp]", ylab= "Fold enrichment")
# for(i in 1:18) {
 # lines(org_gcen[[i]]~x, col=cols[i])
# }

# plot(-100000,-1, xlim=c(min(x), max(x)), ylim=c(0, 10), main="Predicted", xlab= "Distance from gene center [bp]", ylab= "Fold enrichment")
# for(i in 1:18) {
 # lines(pred_gcen[[i]]~x, col=cols[i])
# }

# plot(-100000,-1, xlim=c(min(x), max(x)), ylim=c(0, 10), main="Experimental", xlab= "Distance from polyadenylation cleavage site [bp]", ylab= "Fold enrichment")
# for(i in 1:18) {
 # lines(org_pas[[i]]~x, col=cols[i])
# }

# plot(-100000,-1, xlim=c(min(x), max(x)), ylim=c(0, 10), main="Predicted", xlab= "Distance from polyadenylation cleavage site [bp]", ylab= "Fold enrichment")
# for(i in 1:18) {
 # lines(pred_pas[[i]]~x, col=cols[i])
# }

dev.off()
