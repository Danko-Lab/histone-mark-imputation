mcmc_polii - simulates a gene under a fixed initiation and pause release rate.

mcmc_polii.c should compile fine w/: gcc -o mcmc_polii mcmc_polii.c

makeBed.bsh runs the imputation over a range of initiation and release parameters in INIT and RELE variables.
The resulting data is output to bigWig files in hg19 coordinate system. The gene is on the chromosome chr6_ssto_hap7. It starts at position 1000000 and ends at 25000.
