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
ol_n  <- args[2]
## chromosome
chr.s <- as.numeric(args[3])
## start position of the gene
pos.s <- as.numeric(args[4])

## define start and end of the signal 
pos.e <- pos.s + 5e5
pos.s <- pos.s - 5e5

## import results files
require(data.table)
ol <- paste0("zcat /home/mdp50/rds/rds-rjh234-mrc-epid/Studies/People/Ellie/Olink/GWAS/output/formatted/combined/", gsub("res_invn", "invn_res", ol_n), "_withMarkerName.out.gz")
ol <- data.frame(fread(cmd=ol, sep=" ", colClasses=c("numeric", "character", "numeric", "character", "character", rep("numeric", 6), "character")))

## subset to region of interest (500kb region around the signal)
ol <- subset(ol, chr == chr.s & (pos >= pos.s & pos <= pos.e))
gc()

## rename to ease mapping afterwards
names(ol)[8:11] <- c("beta", "se", "tval", "log10p")

## subset to MAF of at least 1%
ol$MAF <- ifelse(ol$af > .5, 1-ol$af, ol$af)
ol     <- subset(ol, MAF >= .01)

## sort by strongest p-value and store
ii     <- which.max(abs(ol$log10p))
ol     <- ol[ii,]

## add some information to ease mapping afterwards
ol$pGWAS_olink     <- ol_n
ol$markername.soma <- snp 

## store the results
write.table(ol, paste0("/home/mdp50/rds/rds-rjh234-mrc-epid/Studies/People/Maik/compare_SomaLogic_Olink/01_pQTLs_Fenland/output/", ol_n, ".", snp, ".regional"), sep="\t", row.names=F, quote = F)
