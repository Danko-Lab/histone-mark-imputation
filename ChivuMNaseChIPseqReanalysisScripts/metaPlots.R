library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(IRanges)

## Get TSS data.
tss_data <- read.table("/local/storage/data/hg19/k562/copro/promoter.maxTSS.K562.bed", header=TRUE)
gr_tss <- GRanges(seqnames = tss_data$seqnames, ranges = IRanges(start = tss_data$start, end = tss_data$end))

#tss_data <- read.table("/local/storage/data/hg19/k562/dreg/G1.dREG.peak.full.bed.gz")
#midpoints <- round((tss_data$V2 + tss_data$V3) / 2)
#gr_tss <- GRanges(seqnames = tss_data$V1, ranges = IRanges(start = midpoints, end = midpoints))

# Extend each TSS to 2kb upstream and downstream
extended_tss <- resize(gr_tss, width = 4000, fix = "center")

# Create non-overlapping 100 bp windows
windows <- slidingWindows(extended_tss, width = 100, step = 100)
flattened_windows <- unlist(windows)

countBam <- function(bam_file) {
	# Create a read coverage object
	cov <- summarizeOverlaps(flattened_windows, bam_file, mode="Union", singleEnd=TRUE, ignore.strand=TRUE)


	## Get the count vector.
	count_matrix <- assays(cov)$counts

	# Reshape the matrix
	# Number of TSSs
	num_tss <- length(windows)

	# Number of windows per TSS (assuming each TSS has the same number of windows)
	num_windows_per_tss <- length(flattened_windows) / num_tss

	# Reshape the matrix
	reshaped_matrix <- matrix(count_matrix, nrow = num_tss, byrow = TRUE)

	return(colSums(reshaped_matrix))
}


plotNormalizedCounts <- function(combined_data, output_pdf_path, title) {

    # Trim the first and last two windows from each dataset. Avoids displaying edge effect in counting.
    combined_data <- combined_data[-c(1:2, (nrow(combined_data)-1):nrow(combined_data)),]
    print(combined_data)

    # Determine the overall y-axis range across all vectors
    y_range <- range(combined_data, na.rm = TRUE)
    
    # Colors for the different time points (assuming 6 columns: 0h_r1, 0h_r2, 1h_r1, 1h_r2, 4h_r1, 4h_r2)
    colors <- c("red2", "goldenrod1", "purple4")

    # Calculate base pair positions for the x-axis (assuming 100 bp per window)
    window_size <- 100  # window size in base pairs
    locus_radius <- 2000 - 2 * window_size # * window_size # size of the expanded window
    x_axis_labels <- (1:nrow(combined_data)) * window_size - (window_size / 2) - locus_radius
    print(x_axis_labels)

    # Open PDF device
    pdf(output_pdf_path)

    # Create the plot
    matplot(combined_data, 
            type = "l", 
            lty = 1,
	    lwd = 3, 
            col = colors, 
            ylim = y_range, 
            xlab = "Distance to coPRO TSS", 
            ylab = "Normalized Counts", 
            main = title,
            cex.main = 1.5,
            cex.lab = 1.5,
            cex.axis = 1.5,
            xaxt = "n")

    # Add custom X axis
    axis(1, at = 1:nrow(combined_data), labels = x_axis_labels)

    # Add a legend
    legend("topright", 
           legend = c("0h", "1h", "4h"), 
           col = colors,
           lty = 1,
	   lwd = 3,
           cex = 1.5)

    # Close the PDF device
    dev.off()
}



## Count reads, divide by spike-in count.
k27ac_0h_r1 <- countBam(BamFile("K27ac_0h_r1_merged.sorted.bam")) / 4470
k27ac_1h_r1 <- countBam(BamFile("K27ac_1h_r1_merged.sorted.bam")) / 12306
k27ac_4h_r1 <- countBam(BamFile("K27ac_4h_r1_merged.sorted.bam")) / 11241

k27ac_0h_r2 <- countBam(BamFile("K27ac_0h_r2_merged.sorted.bam")) / 1594
k27ac_1h_r2 <- countBam(BamFile("K27ac_1h_r2_merged.sorted.bam")) / 10269
k27ac_4h_r2 <- countBam(BamFile("K27ac_4h_r2_merged.sorted.bam")) / 10939

## Take means of replicates.
k27ac_0h <- colSums(rbind(k27ac_0h_r1, k27ac_0h_r2))
k27ac_1h <- colSums(rbind(k27ac_1h_r1, k27ac_1h_r2))
k27ac_4h <- colSums(rbind(k27ac_4h_r1, k27ac_4h_r2))

## Print metaplots.
plotNormalizedCounts(cbind(k27ac_0h, k27ac_1h, k27ac_4h), "~/chivu_2022/k27ac.pdf", "H3K27ac")

## Now do K4me3.
k4me3_0h_r1 <- countBam(BamFile("K4me3_0h_r1_merged.sorted.bam")) / 4889
k4me3_1h_r1 <- countBam(BamFile("K4me3_1h_r1_merged.sorted.bam")) / 4054
k4me3_4h_r1 <- countBam(BamFile("K4me3_4h_r1_merged.sorted.bam")) / 5810

k4me3_0h_r2 <- countBam(BamFile("K4me3_0h_r2_merged.sorted.bam")) / 1294
k4me3_1h_r2 <- countBam(BamFile("K4me3_1h_r2_merged.sorted.bam")) / 2309
k4me3_4h_r2 <- countBam(BamFile("K4me3_4h_r2_merged.sorted.bam")) / 2898

## Take means of replicates.
k4me3_0h = colSums(rbind(k4me3_0h_r1, k4me3_0h_r2))
k4me3_1h = colSums(rbind(k4me3_1h_r1, k4me3_1h_r2))
k4me3_4h = colSums(rbind(k4me3_4h_r1, k4me3_4h_r2))

## Print metaplot.
plotNormalizedCounts(cbind(k4me3_0h, k4me3_1h, k4me3_4h), "~/chivu_2022/k4me3.pdf", "H3K4me3")




