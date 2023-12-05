## Chivu reanalysis

These script were written to reanalyzes the MNase ChIP-seq data from the SRA files. 

getData.py -- Downloads fastq files from the short read archive, using sraids.txt.  
trim_fastq.bsh -- Trims Illumina adapters for all fastq files.  
align_all(_notrim).bsh -- Aligns reads using BWA.  
merge_bams(_notrim).bsh -- Merges BAM files for two separate sequencing runs for the same sample.   
getSpikeCounts.bsh -- Counts reads mapping to the human or D Iulia genome.  
