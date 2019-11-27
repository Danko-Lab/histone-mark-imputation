## Merges plus and minus strand annotations following the format 
## by Tome and Tippens.

pl <- read.table("maxPl.bed.gz", fill=TRUE)
mn <- read.table("maxMn.bed.gz", fill=TRUE)

## Exclude TIDs without TSS on plus and minus strand.
exclude <- is.na(pl$V9) | is.na(mn$V9) | pl$V1 == "chrX" | pl$V1 == "chrY" | pl$V1 == "chrM"
pl <- pl[!exclude,]
mn <- mn[!exclude,]

## Include a unique ID. Sanity checked and the names are identical.
pl[,4] <- paste(pl$V7,pl$V8,pl$V9, sep="-")
mn[,4] <- paste(mn$V7,mn$V8,mn$V9, sep="-")

## Merge so that every other line is interspersed.
allData <- data.frame()
for(x in 1:NROW(pl)) {
 allData <- rbind(allData, mn[x,], pl[x,]) 
}

## Write it out!
options(scipen=999)
#write.table(pl, "maxPl.filt.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE) 
#write.table(mn, "maxMn.filt.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(allData, "hg19.gm12878.TIDs.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
