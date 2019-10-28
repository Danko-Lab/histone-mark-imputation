
if(0)
{
  load("calc-comb-GM-K562.rdata");
  df.allcomb <- rbind( df.cor.10k$K562,  df.cor.10k$GM);
  
  load("draw-cor-chr22.rdata");
  df.K562.chr22<- df.G1.chr22[df.G1.chr22$win.size==10000,] 
  df.GM.chr22 <- df.GM.chr22[df.GM.chr22$win.size==10000,] 
  df.Alex.chr22 <- df.Alex.Mnase.chr22[df.Alex.Mnase.chr22$win.size==10000,] 
  
  df.K562.chr22$Mark<- unlist(lapply( as.character(df.K562.chr22$plot), function(str.plot) return(strsplit(str.plot, "\\.")[[1]][2]) ))
  df.GM.chr22$Mark<- unlist(lapply( as.character(df.GM.chr22$plot), function(str.plot) return(strsplit(str.plot, "\\.")[[1]][2]) ))
  df.Alex.chr22$Mark<- unlist(lapply( as.character(df.Alex.chr22$plot), function(str.plot) return(strsplit(str.plot, "\\.")[[1]][3]) ))

  df.points <- rbind( data.frame(df.K562.chr22, cols="red"), 
                      data.frame(df.GM.chr22, cols="blue"), 
                      data.frame(df.Alex.chr22, cols="darkgreen") );
  df.points <- df.points[df.points$Mark!="H3k122ac",]
#  df.points <- df.points[df.points$Mark!="H3k9me3",]
#  df.points <- df.points[df.points$Mark!="H4k20me1",]
#  df.points <- df.points[df.points$Mark!="H3k27me3",]

  pdf("boxplot-histone-allcomb-10kb.pdf");
  
  boxplot( pearson~Mark,data=df.allcomb, main="K562+GM12878",  xlab="Histone", ylab="Pearson Correlation", cex.axis=0.5,ylim=c(0,1), range=0) 
  points(factor(df.points$Mark), df.points$pearson, pch=15, cex=1, col=as.character(df.points$cols))
  legend(7, 0.25, c("K562","GM12878","MNase-seq"), fill=c("red","blue","darkgreen"), col=c("red","blue","darkgreen")); 

  boxplot( spearman~Mark,data=df.allcomb, main="K562+GM12878",  xlab="Histone", ylab="Spearman Correlation", cex.axis=0.5,ylim=c(0,1), range=0)
  points(factor(df.points$Mark), df.points$spearman, pch=15, cex=1, col=as.character(df.points$cols))
  legend(7, 0.25, c("K562","GM12878","MNase-seq"), fill=c("red","blue","darkgreen"), col=c("red","blue","darkgreen"));

  boxplot( JSD~Mark,data=df.allcomb, main="K562+GM12878",  xlab="Histone", ylab="Jensen-Shannon Divergence", cex.axis=0.5, range=0, ylim=c(0,0.5))
  points(factor(df.points$Mark), df.points$JSD, pch=15, cex=1, col=as.character(df.points$cols))
  legend(7, 0.25, c("K562","GM12878","MNase-seq"), fill=c("red","blue","darkgreen"), col=c("red","blue","darkgreen"));


  boxplot( mad~Mark,data=df.allcomb, main="K562+GM12878",  xlab="Histone", ylab="Mean Absolute Deviation", cex.axis=0.5, range=0)
  points(factor(df.points$Mark), df.points$mad, pch=15, cex=1, col=as.character(df.points$cols))
  legend(7.25, 15000, c("K562","GM12878","MNase-seq"), fill=c("red","blue","darkgreen"), col=c("red","blue","darkgreen"));

  dev.off();
  
}


if(1)
{
    # loading correlations based on raw signals for chr 22. 
    load(file="draw-cor-chr22.rdata")
    load("calc-comb-GM-K562.rdata");
    df.cor <- rbind(df.cor.10k$K562[,-11], df.cor.10k$GM[,-11], 
              df.cor.1k$K562[,-11],  df.cor.1k$GM[,-11], 
              df.cor.100$K562[,-11], df.cor.100$GM[,-11], 
              df.cor.10$K562[,-11],  df.cor.10$GM[,-11]);

    df.G1.chr22$Cell <- unlist(lapply( as.character(df.G1.chr22$plot), function(str.plot) return(strsplit(str.plot, "\\.")[[1]][1]) ))
    df.G1.chr22$Mark <- unlist(lapply( as.character(df.G1.chr22$plot), function(str.plot) return(strsplit(str.plot, "\\.")[[1]][2]) ))
    df.GM.chr22$Cell <- unlist(lapply( as.character(df.GM.chr22$plot), function(str.plot) return(strsplit(str.plot, "\\.")[[1]][1]) ))
    df.GM.chr22$Mark <- unlist(lapply( as.character(df.GM.chr22$plot), function(str.plot) return(strsplit(str.plot, "\\.")[[1]][2]) ))
#    df.Alex.Mnase.chr22$Mark <- unlist(lapply( as.character(df.Alex.Mnase.chr22$plot), function(str.plot) return(strsplit(str.plot, "\\.")[[1]][3]) ))
#    df.Alex.Mnase.chr22$Cell <- "G1";

#    df.chr22 <- rbind(df.G1.chr22, df.GM.chr22, df.Alex.Mnase.chr22)[,-c(9:10)]
    df.chr22 <- rbind(df.G1.chr22, df.GM.chr22)[,-c(9:10)]
    df.cor <- rbind(data.frame(df.cor, Pred=0),data.frame(df.chr22,Pred=1))
    df.cor <- df.cor[df.cor$Mark!="H3k122ac",]
    df.cor <- df.cor[df.cor$Mark!="H3k9me3",]
    df.cor <- df.cor[df.cor$Mark!="H4k20me1",]
    df.cor <- df.cor[df.cor$Mark!="H3k27me3",]
    
    df.cor$win.size<- as.numeric(df.cor$win.size);
    df.cor$pearson <- as.numeric(df.cor$pearson);
    df.cor$spearman<- as.numeric(df.cor$spearman);
    df.cor$mad<- as.numeric(df.cor$mad);
    df.cor$JSD<- as.numeric(df.cor$JSD);

    #There is one NA in the result
    df.cor <- df.cor[!is.na(df.cor$pearson),]

    df.cor$mark_idx  <- as.numeric(df.cor$Mark)
    df.cor$x  <-  log10(df.cor$win.size) - 0.4 - 0.8/18 + 0.8/9*df.cor$mark_idx  

    df.plot<- df.cor[df.cor$Pred==0, c("pearson", "Mark", "win.size")];
    colnames(df.plot) <- c("Cor", "Mark", "Size" );                    

    df.points <- df.cor [df.cor$Pred==1 ,c("x", "pearson",  "Mark" )];
    colnames(df.points) <- c("x", "Cor", "Mark"  );                    

    # check the order of plot factor in two data frames. 
    all(levels(df.cor$Mark) == levels(df.points$Mark))

    df.plot$Size <- factor(df.plot$Size)
    df.plot$Mark<- factor(df.plot$Mark)

    cols <-  c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#1111d6")    
    df.points$cols <- cols[as.factor(df.points$Mark)];
    
    
    pdf("boxplot-cor-comb-xkb.pdf", width=8, height=5);
    
library(ggplot2)
   ggplot(data = df.plot, aes(x = Size, y = Cor)) + 
     geom_boxplot(aes(fill = Mark), width = 0.8,outlier.shape = NA) +
     geom_point(data = df.points, size = 1.5, shape=21, colour="black",  show.legend=FALSE, aes(x=x, y=Cor, fill=Mark) ) +
     scale_fill_manual(values=cols)+
     xlab("Window size") + 
     ylab("Correlation (Pearson)") + 
     #scale_x_discrete(breaks=c(1,2,3, 4), labels=c("10", "100", "1000", "10000")) + 
     theme_bw()

dev.off()

}

