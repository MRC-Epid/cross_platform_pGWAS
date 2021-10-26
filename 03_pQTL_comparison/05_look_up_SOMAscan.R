#!/usr/bin/env Rscript

## script to search for regional lead variant
## in Olink files - Maik Pietzner  11/11/2019
rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)

## assign names of the proteins to be used
snp   <- args[1]
sl_n  <- args[2]
## chromosome
chr.s <- as.numeric(args[3])
## start position of the gene
pos.s <- as.numeric(args[4])

## define start and end of the signal 
pos.e <- pos.s + 5e5
pos.s <- pos.s - 5e5

## read the relevant data
res.soma        <- paste0("zcat ~/rds/rds-rjh234-mrc-epid/Studies/People/Ellie/SomaLogic/GWAS_discovery_subsets/Meta-Analysis/Fenland_auto_chrX_filtered/output/",
                          sl_n,"_Fenland_MA_auto_chrX_filtered.txt.gz",
                          " | awk -v chr=", chr.s, " -v low=", pos.s, " -v upp=", pos.e, 
                          " '$18 == chr && $19 >= low && $19 <= upp {print $0}' -")
res.soma        <- data.table::fread(cmd = res.soma, sep = "\t", header = F, data.table = F)
## add names
names(res.soma) <- c("rsid", "MarkerName", "Allele1", "Allele2", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Effect", "StdErr", "Pvalue", "Direction", "HetSq", "HetChiSq",
                     "HetDf", "HetPVal", "TotalSampleSize", "chr", "pos")

print(dim(res.soma))
## restrict to MAF > 1%
res.soma$MAF    <- ifelse(res.soma$Freq1 > .5, 1-res.soma$Freq1, res.soma$Freq1)
res.soma        <- subset(res.soma, MAF >= .01)

## get the strongest association
ii       <- which.max(abs(res.soma$Effect/res.soma$StdErr))   
res.soma <- res.soma[ii,] 

## add some names to ease mapping afterwards
res.soma$pheno            <- sl_n
res.soma$Markername.Olink <- snp

## store the results
write.table(res.soma, paste0("/home/mdp50/rds/rds-rjh234-mrc-epid/Studies/People/Maik/compare_SomaLogic_Olink/01_pQTLs_Fenland/output/", sl_n, ".", snp, ".regional"), sep="\t", row.names=F, quote = F)

