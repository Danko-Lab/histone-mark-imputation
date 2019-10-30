# Col3: pearson
# Col4: spearman
# Col5: MAD

get_marker_cell_matrix<-function(win.size=1000, idx.col=4)
{
   load("draw-cor-chr22.rdata");

   df.K562.chr22  <- df.G1.chr22[df.G1.chr22$win.size==win.size,];
   df.GM.chr22    <- df.GM.chr22[df.GM.chr22$win.size==win.size,];
   df.CD4.chr22   <- df.CD4.chr22[df.CD4.chr22$win.size==win.size,];
   df.HCT.chr22   <- df.HCT.chr22[df.HCT.chr22$win.size==win.size,];
   df.HELA.chr22  <- df.HELA.chr22[df.HELA.chr22$win.size==win.size,];
   df.Horse.chr22 <- df.Horse.chr22[df.Horse.chr22$win.size==win.size,];
   df.Alex.chr22  <- df.Alex.Mnase.chr22[df.Alex.Mnase.chr22$win.size==win.size,];
   df.mESC.chr1   <- df.mESC.chr1[df.mESC.chr1$win.size==win.size,];

   df.K562.chr22$plot <- unlist(lapply( as.character(df.K562.chr22$plot), function(str.plot) return(strsplit(str.plot, "\\.")[[1]][2]) ))
   df.GM.chr22$plot   <- unlist(lapply( as.character(df.GM.chr22$plot), function(str.plot) return(strsplit(str.plot, "\\.")[[1]][2]) ))
   df.CD4.chr22$plot  <- unlist(lapply( as.character(df.CD4.chr22$plot), function(str.plot) return(strsplit(str.plot, "\\.")[[1]][2]) ))
   df.HCT.chr22$plot  <- unlist(lapply( as.character(df.HCT.chr22$plot), function(str.plot) return(strsplit(str.plot, "\\.")[[1]][2]) ))
   df.HELA.chr22$plot <- unlist(lapply( as.character(df.HELA.chr22$plot), function(str.plot) return(strsplit(str.plot, "\\.")[[1]][2]) ))
   df.Horse.chr22$plot<- unlist(lapply( as.character(df.Horse.chr22$plot), function(str.plot) return(strsplit(str.plot, "\\.")[[1]][3]) ))
   df.Alex.chr22$plot <- unlist(lapply( as.character(df.Alex.chr22$plot), function(str.plot) return(strsplit(str.plot, "\\.")[[1]][3]) ))
   df.mESC.chr1$plot  <- unlist(lapply( as.character(df.mESC.chr1$plot), function(str.plot) return(strsplit(str.plot, "\\.")[[1]][2]) ))
   
   L <- do.call("rbind", lapply(as.character(df.K562.chr22$plot), function(str.mark) {
       C1=df.K562.chr22 [which(df.K562.chr22$plot == str.mark), idx.col];
       C2=df.GM.chr22   [which(df.GM.chr22$plot   == str.mark), idx.col];
       C3=df.HCT.chr22  [which(df.HCT.chr22$plot  == str.mark), idx.col];
       C4=df.HELA.chr22 [which(df.HELA.chr22$plot == str.mark), idx.col];
       C5=df.CD4.chr22  [which(df.CD4.chr22$plot  == str.mark), idx.col];
       C6=df.Horse.chr22[which(df.Horse.chr22$plot== str.mark), idx.col];
       C7=df.Alex.chr22 [which(df.Alex.chr22$plot == str.mark), idx.col];
       C8=df.mESC.chr1  [which(df.mESC.chr1$plot == str.mark),  idx.col];

       return(c(if(NROW(C1)==0)NA else C1,
                if(NROW(C2)==0)NA else C2,
                if(NROW(C3)==0)NA else C3,
                if(NROW(C4)==0)NA else C4,
                if(NROW(C5)==0)NA else C5,
                if(NROW(C6)==0)NA else C6,
                if(NROW(C7)==0)NA else C7,
                if(NROW(C8)==0)NA else C8))
       }));
       
   colnames(L) <- c("K562", "GM", "HCT", "HELA", "CD4", "Horse", "Alex", "mESC");
   rownames(L) <- as.character(df.K562.chr22$plot)
   
   return(L);
}

colmap <- c("#e5f5f9", "#99d8c9", "#2ca25f");
colmap2<- c("#efedf5", "#bcbddc", "#756bb1")

if(1)
{
  library(gplots)
  library(RColorBrewer)

  mat_data <- round(get_marker_cell_matrix(10000, 3),2)

  # creates a own color palette from red to green
  my_palette <- colorRampPalette(colmap)(n = 299)

  # (optional) defines the color breaks manually for a "skewed" color transition
  col_breaks = c(seq(0,0.33,length=100),   # for red
                 seq(0.34,0.66,length=100),            # for yellow
                 seq(0.67,1,length=100))              # for green

# creates a 5 x 5 inch image
  pdf("draw-cormat-pearson-10k.pdf")

 heatmap.2(mat_data,
  cellnote = mat_data,  # same data set for cell labels
  main = "Pearson Correlation", # heat map title
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="row",     # only draw a row dendrogram
  Colv="NA",            # turn off column clustering
           sepwidth=c(0.005,0.005),
           sepcolor="black",
           colsep=0:ncol(mat_data),
           rowsep=0:nrow(mat_data))

  dev.off()
}

if(1)
{
  library(gplots)
  library(RColorBrewer)

  # Col4: pearson
  # Col5: spearman
  # Col6: MAD
  mat_data <- round(get_marker_cell_matrix(10000, 4),2)

  # creates a own color palette from red to green
  my_palette <- colorRampPalette(colmap2)(n = 299)

  # (optional) defines the color breaks manually for a "skewed" color transition
  col_breaks = c(seq(0,0.33,length=100),   # for red
                 seq(0.34,0.66,length=100),            # for yellow
                 seq(0.67,1,length=100))              # for green

# creates a 5 x 5 inch image
  pdf("draw-cormat-spearman-10k.pdf")

 heatmap.2( mat_data,
  cellnote = mat_data,  # same data set for cell labels
  main = "Spearman Correlation", # heat map title
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="row",     # only draw a row dendrogram
  Colv="NA",            # turn off column clustering
           sepwidth=c(0.005,0.005),
           sepcolor="black",
           colsep=0:ncol(mat_data),
           rowsep=0:nrow(mat_data))

  dev.off()
}

if(1)
{
  library(gplots)
  library(RColorBrewer)

  mat_data <- round(get_marker_cell_matrix(10000, 5),2)

  # creates a own color palette from red to green
  my_palette <- rev(colorRampPalette(colmap2)(n = 299))

  # (optional) defines the color breaks manually for a "skewed" color transition
  col_breaks = c(seq(0,700,length=100),   # for red
                 seq(701,1400,length=100),# for yellow
                 seq(1401,2100,length=100))   # for green

# creates a 5 x 5 inch image
  pdf("draw-cormat-mad-10k.pdf")

 heatmap.2(mat_data,
  cellnote = mat_data,  # same data set for cell labels
  main = "MAD",         # heat map title
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="row",     # only draw a row dendrogram
  Colv="NA",            # turn off column clustering
           sepwidth=c(0.005,0.005),
           sepcolor="black",
           colsep=0:ncol(mat_data),
           rowsep=0:nrow(mat_data))

  dev.off()
}


if(1)
{
  library(gplots)
  library(RColorBrewer)

  mat_data <- round(get_marker_cell_matrix(10000, 6),2)

  # creates a own color palette from red to green
  my_palette <- rev(colorRampPalette(colmap2)(n = 299))

  # (optional) defines the color breaks manually for a "skewed" color transition
  col_breaks = c(seq(0,700,length=100),   # for red
                 seq(701,1400,length=100),# for yellow
                 seq(1401,2100,length=100))   # for green

# creates a 5 x 5 inch image
  pdf("draw-cormat-JSD-10k.pdf")

 heatmap.2(mat_data,
  cellnote = mat_data,  # same data set for cell labels
  main = "MAD",         # heat map title
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="row",     # only draw a row dendrogram
  Colv="NA",            # turn off column clustering
           sepwidth=c(0.005,0.005),
           sepcolor="black",
           colsep=0:ncol(mat_data),
           rowsep=0:nrow(mat_data))

  dev.off()
}


colmap <- c("#e5f5f9", "#99d8c9", "#2ca25f");
colmap2<- c("#efedf5", "#bcbddc", "#756bb1")

if(1)
{
  library(gplots)
  library(RColorBrewer)

  # Col4: pearson
  # Col5: spearman
  # Col6: MAD
  mat_data <- round(get_marker_cell_matrix(1000, 3),2)

  # creates a own color palette from red to green
  my_palette <- colorRampPalette(colmap)(n = 299)

  # (optional) defines the color breaks manually for a "skewed" color transition
  col_breaks = c(seq(0,0.33,length=100),   # for red
                 seq(0.34,0.66,length=100),            # for yellow
                 seq(0.67,1,length=100))              # for green

# creates a 5 x 5 inch image
  pdf("draw-cormat-pearson-1k.pdf")

 heatmap.2(mat_data,
  cellnote = mat_data,  # same data set for cell labels
  main = "Pearson Correlation", # heat map title
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="row",     # only draw a row dendrogram
  Colv="NA",            # turn off column clustering
           sepwidth=c(0.005,0.005),
           sepcolor="black",
           colsep=0:ncol(mat_data),
           rowsep=0:nrow(mat_data))

  dev.off()
}

if(1)
{
  library(gplots)
  library(RColorBrewer)

  # Col4: pearson
  # Col5: spearman
  # Col6: MAD
  mat_data <- round(get_marker_cell_matrix(1000, 4),2)

  # creates a own color palette from red to green
  my_palette <- colorRampPalette(colmap2)(n = 299)

  # (optional) defines the color breaks manually for a "skewed" color transition
  col_breaks = c(seq(0,0.33,length=100),   # for red
                 seq(0.34,0.66,length=100),            # for yellow
                 seq(0.67,1,length=100))              # for green

# creates a 5 x 5 inch image
  pdf("draw-cormat-spearman-1k.pdf")

 heatmap.2( mat_data,
  cellnote = mat_data,  # same data set for cell labels
  main = "Spearman Correlation", # heat map title
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="row",     # only draw a row dendrogram
  Colv="NA",            # turn off column clustering
           sepwidth=c(0.005,0.005),
           sepcolor="black",
           colsep=0:ncol(mat_data),
           rowsep=0:nrow(mat_data))

  dev.off()
}
