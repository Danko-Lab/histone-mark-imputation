##
## Draws a beeswarm plot comparing molecular assays.

load("/workdir/zw355/proj/prj15-histone/h2h-model/summary.h2h.train.L1.rdata")
rdOne  <- rDist

load("/workdir/zw355/proj/prj15-histone/h2h-model/summary.hx2h.train.L1.rdata")
rdMult <- rDist

require(reshape2)
require(RColorBrewer)
require(beeswarm)
rd <- melt(rdOne[order(rowMeans(rdOne)),], value.name="L1.Norm", varnames=c("train", "predict"))
rd <- rbind(rd, melt(rdMult, value.name="L1.Norm", varnames=c("train", "predict")))
rd <- rd[rd$L1.Norm > 0,]

cols <- brewer.pal(10, "Paired")

pdf("BeeswarmPlot.pdf", width = 10, height = 5)

beeswarm(L1.Norm~train, data = rd, pch=16, cex=1.5, pwcol = cols[as.factor(predict)], las=2)
bxplot(L1.Norm~train, data=rd, probs = 0.5, add=TRUE)
bxplot(L1.Norm~train, data=rd[rd$predict!="H3k27me3" & rd$predict!="H3k9me3",], probs = 0.5, add=TRUE, col="dark blue")
legend("topright", legend = levels(rd$predict), title = "Predicted mark", pch=16, col = cols)

dev.off()

