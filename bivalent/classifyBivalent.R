## Classifies bivalent promoters using imputed data.
## 

biv <- read.table("poised-plosgen.mm8.bed", fill=TRUE)
coord <- read.table("bivalent.mm9.bed")

biv <- biv[biv$V6 %in% coord$V4,]

biv[,2] <- coord[match(biv$V6, coord$V4),2]
biv[,3] <- coord[match(biv$V6, coord$V4),3]

biv <- biv[biv[,1] == "chr1",]

summary(biv$V8)

## Next, count reads in windows near each gene.
## Build a RF model to classify based on these values.

# Get bigWigs.
require(bigWig)
K4 <- load.bigWig("raw.MM9.mEsc.H3k4me3.S1.pred.bw")
K27<- load.bigWig("raw.MM9.mEsc.H3k27me3.S1.pred.bw")

K4e<- load.bigWig("mm9.GSM307618_ES.H3K4me3.bw")
K27e<- load.bigWig("mm9.GSM307619_ES.H3K27me3.bw")

bed <- biv[,c(1:3,6,8,4)]
K4_CountMatrix <- matrix(unlist(bed.step.bpQuery.bigWig(K4, center.bed(bed, 1000, 1000), step=250, abs.value=TRUE, op="sum")), nrow= NROW(bed), byrow=TRUE)
K27_CountMatrix<- matrix(unlist(bed.step.bpQuery.bigWig(K27,center.bed(bed, 60000, 60000), step=15000, abs.value=TRUE, op="sum")), nrow= NROW(bed), byrow=TRUE)
K4e_CountMatrix <- matrix(unlist(bed.step.bpQuery.bigWig(K4e, center.bed(bed, 1000, 1000), step=250, abs.value=TRUE, op="sum")), nrow= NROW(bed), byrow=TRUE)
K27e_CountMatrix<- matrix(unlist(bed.step.bpQuery.bigWig(K27e,center.bed(bed, 60000, 60000), step=15000, abs.value=TRUE, op="sum")), nrow= NROW(bed), byrow=TRUE)

## Have to reverse minus strand windows?
for(x in which(bed$V4 == "-")) { K4_CountMatrix[x,]  <- rev(K4_CountMatrix[x,]); K4e_CountMatrix[x,]  <- rev(K4e_CountMatrix[x,]) }
for(x in which(bed$V4 == "-")) { K27_CountMatrix[x,] <- rev(K27_CountMatrix[x,]); K27e_CountMatrix[x,] <- rev(K27e_CountMatrix[x,]) }

## Now build the random forest.
require(randomForest)
df <- data.frame(class= bed$V8, K4= K4_CountMatrix, K27= K27_CountMatrix)
dfe<- data.frame(class= bed$V8, K4= K4e_CountMatrix, K27= K27e_CountMatrix)

rf <- randomForest(class ~ ., data=df)
rf

rfe <- randomForest(class ~ ., data=dfe)
rfe

## Narrow to two classes ... bivalent, not bivalent.
tc <- rep("nbv", NROW(df)); tc[bed$V8 == "K4+K27"] <- "bv" ## Set up as predicting bivalent from everything else.
df[,1] <- as.factor(tc)
dfe[,1]<- as.factor(tc) 

indxTrain <- c(sample(which(bed$V8 != "K4+K27"), 100), sample(which(bed$V8 == "K4+K27"), 100)) ## Matched training set.

rf <- randomForest(class ~ ., data=df, subset=indxTrain); rf ## Train the RF model. Print confusion matrix.
rfe<- randomForest(class ~ ., data=dfe, subset=indxTrain); rfe

lr <- glm(class ~ ., data=df, subset=indxTrain, family=binomial)

## Test random forest.
pd <- data.frame(predict= predict(rf, df[-indxTrain,]), exp= df[-indxTrain,"class"]) ## Predict everything not in training set ...
pde<- data.frame(predict= predict(rfe,dfe[-indxTrain,]),exp= dfe[-indxTrain,"class"])
indxTest <- c(sample(which(pd$exp == "nbv"), 100), sample(which(pd$exp == "bv"), 100)) ## Form a matched set.
pd <- pd[indxTest,]
pde<- pde[indxTest,] 
sum(pd$predict == pd$exp)/NROW(pd) ## Compute total accuracy.
sum(pde$predict == pde$exp)/NROW(pde) ## Compute total accuracy.

data.frame(bv= c(sum(pd$exp == "bv" & pd$predict == "bv"), 
                        sum(pd$exp == "bv" & pd$predict == "nbv")), 
                nbv= c(sum(pd$exp == "nbv" & pd$predict == "bv"), 
                        sum(pd$exp == "nbv" & pd$predict == "nbv"))) ## Confusion matrix.
importance(rf)

## Plot ROC curve.
library(ROCR)
pred <- predict(rf, df[-indxTrain,], type = "prob")
pred <- prediction(pred[,1][indxTest], df[-indxTrain,"class"][indxTest] == "bv")

prede <- predict(rfe, dfe[-indxTrain,], type = "prob")
prede <- prediction(prede[,1][indxTest], dfe[-indxTrain,"class"][indxTest] == "bv")

roc <- performance(pred, "tpr", "fpr")
prc <- performance(pred, "prec", "rec")

roce <- performance(prede, "tpr", "fpr")
prce <- performance(prede, "prec", "rec")

 plot(prc, ylim=c(0,1))
 plot(prce, add=TRUE, col="green")
 abline(h=0.5, col="gray")

pdf("ROC.PRC.pdf")
 plot(roc); abline(0,1)
 plot(roce, add=TRUE, col="green")
 
 plot(prc, ylim=c(0,1))
 plot(prce, add=TRUE, col="green")
 abline(h=0.5, col="gray")
dev.off()

## Test logistic regression.
pd_v <- predict(lr, df[-indxTrain,])
pd <- data.frame(predict= pd_v, exp= df[-indxTrain,"class"]) ## Predict everything not in training set ...
indxTest <- c(sample(which(pd$exp == "nbv"), 100), sample(which(pd$exp == "bv"), 100)) ## Form a matched set.
pd[pd_v >=0,"predict"] <- "nbv"
pd[pd_v < 0,"predict"] <- "bv"
pd <- pd[indxTest,] 
sum(pd$predict == pd$exp)/NROW(pd) ## Compute total accuracy.
data.frame(bv= c(sum(pd$exp == "bv" & pd$predict == "bv"), 
			sum(pd$exp == "bv" & pd$predict == "nbv")), 
		nbv= c(sum(pd$exp == "nbv" & pd$predict == "bv"), 
			sum(pd$exp == "nbv" & pd$predict == "nbv"))) ## Confusion matrix.

