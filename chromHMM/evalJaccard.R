# First build a matrix where rows indicate state in V5, the fraction of times observed in V4. 

## On CBSUDanko.
cd /local/ftp/pub/hub/histone/ChromHMM/
bedtools unionbedg -i MS_K562_pred_adjust2_seg.bed.gz MS_K562_org_seg.bed.gz MS_K562_UWother_seg.bed.gz MS_K562_pred_seg.bed.gz | gzip > ~/tmp.chomHMMcomparison.bed.gz


library(ggplot2)
require(reshape2) 
require(RColorBrewer)
require(gridExtra) 

tm <- read.table("~/tmp.chomHMMcomparison.bed.gz")

## Compute Jaccard using definition here: https://bedtools.readthedocs.io/en/latest/content/tools/jaccard.html. 

## Column is experimental, row is predicted
k562 <- NULL
for(i in 1:18) {
 k562 <- rbind(k562, sapply(c(1:18), 
    function(x) {
        sum((tm$V3-tm$V2)[tm$V5 == i & tm$V4 == x]) / (sum((tm$V3-tm$V2)[tm$V5 == i | tm$V4 == x]) - sum((tm$V3-tm$V2)[tm$V5 == i & tm$V4 == x]))
 }))
}

k562dr <- NULL
for(i in 1:18) {
 k562dr <- rbind(k562dr, sapply(c(1:18), 
    function(x) {
        sum((tm$V3-tm$V2)[tm$V5 == i & tm$V6 == x]) / (sum((tm$V3-tm$V2)[tm$V5 == i | tm$V6 == x]) - sum((tm$V3-tm$V2)[tm$V5 == i & tm$V6 == x]))
 }))
}

k562pr1 <- NULL
for(i in 1:18) {
 k562pr1 <- rbind(k562pr1, sapply(c(1:18), 
    function(x) {
        sum((tm$V3-tm$V2)[tm$V5 == i & tm$V7 == x]) / (sum((tm$V3-tm$V2)[tm$V5 == i | tm$V7 == x]) - sum((tm$V3-tm$V2)[tm$V5 == i & tm$V7 == x]))
 }))
}
 
k562prvdr <- NULL
for(i in 1:18) {
 k562prvdr <- rbind(k562prvdr, sapply(c(1:18), 
    function(x) {
        sum((tm$V3-tm$V2)[tm$V6 == i & tm$V4 == x]) / (sum((tm$V3-tm$V2)[tm$V6 == i | tm$V4 == x]) - sum((tm$V3-tm$V2)[tm$V6 == i & tm$V4 == x]))
 }))
}

## melt -> 1st variable is row, 2nd is column
p1 <- ggplot(melt(k562, value.name="Jaccard", varnames=c("Experimental", "Predicted.adjust2_seg")), aes(x = Predicted.adjust2_seg, y = Experimental)) + 
 # scale_y_discrete(name="", labels= c(1:18)) +
  geom_tile(color = "white")+
  geom_raster(aes(fill=Jaccard)) + 
  scale_fill_gradient2(low = "#FFFFFF", high = "#012345", limit = c(0,1), space = "Lab", name="Jaccard\nDistance") + #mid = "orange", midpoint = 0.125, 
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))

p2 <- ggplot(melt(k562dr, value.name="Jaccard", varnames=c("Experimental.ENCODE", "Experimental.Alt")), aes(x = Experimental.Alt, y = Experimental.ENCODE)) + 
 # scale_y_discrete(name="", labels= c(1:18)) +
  geom_tile(color = "white")+
  geom_raster(aes(fill=Jaccard)) + 
  scale_fill_gradient2(low = "#FFFFFF", high = "#012345", limit = c(0,1), space = "Lab", name="Jaccard\nDistance") + #mid = "orange", midpoint = 0.125, 
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))

p3 <- ggplot(melt(k562pr1, value.name="Jaccard", varnames=c("Experimental", "Predicted.pred_seg")), aes(x = Predicted.pred_seg, y = Experimental)) + 
  geom_tile(color = "white")+
  geom_raster(aes(fill=Jaccard)) + 
  scale_fill_gradient2(low = "#FFFFFF", high = "#012345", limit = c(0,1), space = "Lab", name="Jaccard\nDistance") + #mid = "orange", midpoint = 0.125, 
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))

pdf("Jaccard.pdf", width=9, height=3)
grid.arrange(p1, p3, p2, nrow= 1)
dev.off()
  
sapply(1:NROW(k562), function(i) {k562prvdr[i,i]})


cols=c("#FF0000", "#FF4500", "#FF4500", "#FF4500", "#008000", "#006400", "#006400", "#006400", "#FFC34D",  "#FFC34D", "#FFFF00", "#65CDAA", "#8A91D0", "#CD5C5C", "#BDB76B", "#808080", "#C0C0C0", "#000000")

pdf("JaccardScatterplot.pdf", width=5, height=5)

plot( sapply(1:NROW(k562), function(i) {k562dr[i,i]}), sapply(1:NROW(k562), function(i) {k562pr1[i,i]}), xlim=c(0,2.2), ylim=c(0,2.2), pch=19, col=cols, xlab="Jaccard(Experimental, Experimental)", ylab="Jaccard(Predicted,Experimental)", cex=3)
#abline(0,1, lty="dotted", col="light gray")
#points(k562dr[18,18], k562pr1[18,18], cex=3)

plot( c(k562dr), c(k562pr1), xlim=c(0,2.2), ylim=c(0,2.2), xlab="Jaccard(Experimental, Experimental)", col="light gray", ylab="Jaccard(Predicted,Experimental)")
#abline(0,1, lty="dotted", col="light gray")
points(sapply(1:NROW(k562), function(i) {k562dr[i,i]}), sapply(1:NROW(k562), function(i) {k562pr1[i,i]}), pch=19, col=cols)
points(k562dr[18,18], k562pr1[18,18])

dev.off()

cor( sapply(1:NROW(k562), function(i) {k562dr[i,i]}), sapply(1:NROW(k562), function(i) {k562pr1[i,i]}))
# 0.9133127 # Spearman
# 0.9158804 # Pearson

mean(sapply(1:NROW(k562), function(i) {k562pr1[i,i]})-sapply(1:NROW(k562), function(i) {k562dr[i,i]}))
#[1] -0.07964478

mean((sapply(1:NROW(k562), function(i) {k562[i,i]})-sapply(1:NROW(k562), function(i) {k562pr1[i,i]}))/ (sapply(1:NROW(k562), function(i) {k562pr1[i,i]})))

median(sapply(1:NROW(k562), function(i) {k562pr1[i,i]})-sapply(1:NROW(k562), function(i) {k562dr[i,i]}))
#[1] -0.11899

mean(sapply(11:18, function(i) {k562pr1[i,i]})-sapply(11:18, function(i) {k562dr[i,i]}))
  
###################################################################################
bedtools unionbedg -i MS_GM12878_pred_seg.bed.gz MS_GM12878_org_seg.bed.gz MS_K562_org_seg.bed.gz | gzip > ~/tmp.chomHMMcomparison.bed.gz

library(ggplot2)
require(reshape2) 
require(RColorBrewer)
require(gridExtra) 

tm <- read.table("~/tmp.chomHMMcomparison.bed.gz")

## Column is experimental, row is predicted
gm <- NULL
for(i in 1:18) {
 gm <- rbind(gm, sapply(c(1:18), 
    function(x) {
        sum((tm$V3-tm$V2)[tm$V5 == i & tm$V4 == x]) / (sum((tm$V3-tm$V2)[tm$V5 == i | tm$V4 == x]) - sum((tm$V3-tm$V2)[tm$V5 == i & tm$V4 == x]))
 }))
}

gm.k562 <- NULL
for(i in 1:18) {
 gm.k562 <- rbind(gm.k562, sapply(c(1:18), 
    function(x) {
        sum((tm$V3-tm$V2)[tm$V5 == i & tm$V6 == x]) / (sum((tm$V3-tm$V2)[tm$V5 == i | tm$V6 == x]) - sum((tm$V3-tm$V2)[tm$V5 == i & tm$V6 == x]))
 }))
}

p1 <- ggplot(melt(gm, value.name="Jaccard", varnames=c("Experimental", "Predicted")), aes(x = Predicted, y = Experimental)) + 
 # scale_y_discrete(name="", labels= c(1:18)) +
  geom_tile(color = "white")+
  geom_raster(aes(fill=Jaccard)) + 
  scale_fill_gradient2(low = "#FFFFFF", high = "#012345", limit = c(0,1), space = "Lab", name="Jaccard\nDistance") + #mid = "orange", midpoint = 0.125, 
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))

p2 <- ggplot(melt(gm.k562, value.name="Jaccard", varnames=c("Experimental", "Predicted")), aes(x = Predicted, y = Experimental)) + 
 # scale_y_discrete(name="", labels= c(1:18)) +
  geom_tile(color = "white")+
  geom_raster(aes(fill=Jaccard)) + 
  scale_fill_gradient2(low = "#FFFFFF", high = "#012345", limit = c(0,1), space = "Lab", name="Jaccard\nDistance") + #mid = "orange", midpoint = 0.125, 
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))

grid.arrange(p1, p2, nrow= 1)
