#!/usr/bin/env Rscript

## script to run confirmative phenotype coloc for SomaScan and Olink
## Maik Pietzner 06/05/2021
rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)

setwd("/rds/project/rjh234/rds-rjh234-mrc-epid/Studies/People/Maik/compare_SomaLogic_Olink/14_protein_net/")

## --> import parameters <-- ##

## aptamer
soma  <- args[1]
## Olink
olink <- args[2]
## trait
trait <- args[3]
## chromosome
chr.s <- as.numeric(args[4])
## start position of the region
pos.s <- as.numeric(args[5])
## end position of the region
pos.e <- as.numeric(args[6])


#-----------------------------------------#
##-- 	       import protein data       --##
#-----------------------------------------#

## read the relevant data
res.soma         <- paste0("zcat ~/rds/rds-rjh234-mrc-epid/Studies/People/Ellie/SomaLogic/GWAS_discovery_subsets/Meta-Analysis/Fenland_auto_chrX_filtered/output/",
                           soma,"_Fenland_MA_auto_chrX_filtered.txt.gz",
                           " | awk -v chr=", chr.s, " -v low=", pos.s, " -v upp=", pos.e, 
                           " '$18 == chr && $19 >= low && $19 <= upp {print $0}' -")
res.soma         <- data.table::fread(cmd = res.soma, sep = "\t", header = F, data.table = F)
## add names
names(res.soma)  <- c("rsid", "MarkerName", "Allele1", "Allele2", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Effect", "StdErr", "Pvalue", "Direction", "HetSq", "HetChiSq",
                      "HetDf", "HetPVal", "TotalSampleSize", "chr", "pos")
## keep only needed
res.soma         <- res.soma[, c("rsid", "MarkerName", "Freq1", "Allele1", "Allele2", "Effect", "StdErr", "Pvalue", "TotalSampleSize", "chr", "pos")]

## apply MAF filter
res.soma$MAF     <- ifelse(res.soma$Freq1 > .5, 1 - res.soma$Freq1, res.soma$Freq1)

## make Alleles to upper
res.soma$Allele1 <- toupper(res.soma$Allele1)
res.soma$Allele2 <- toupper(res.soma$Allele2)

#-----------------------------------------#
##-- 	       import protein data       --##
#-----------------------------------------#

## read the relevant data
res.olink                 <- paste0("zcat /home/mdp50/rds/rds-rjh234-mrc-epid/Studies/People/Ellie/Olink/GWAS/output/formatted/combined/",
                                    olink,"_withMarkerName.out.gz",
                                    " | awk -v chr=", chr.s, " -v low=", pos.s, " -v upp=", pos.e, 
                                    " '$1 == chr && $3 >= low && $3 <= upp {print $0}' -")
res.olink                 <- data.table::fread(cmd = res.olink, sep = " ", header = F, data.table = F)
## add names
names(res.olink)          <- c("chr", "rsid", "pos", "Allele2", "Allele1", "Freq1", "info", "Effect", "StdErr", "tval", "log10p", "MarkerName")
## add sample size and Pvalue
res.olink$TotalSampleSize <- 485
res.olink$Pvalue          <- 10^-res.olink$log10p
## keep only needed
res.olink                 <- res.olink[, c("rsid", "MarkerName", "Freq1", "Allele1", "Allele2", "Effect", "StdErr", "Pvalue", "TotalSampleSize", "chr", "pos")]

#-----------------------------------------#
##-- 	         import trait data       --##
#-----------------------------------------#

## is included in the Open GWAS data base already
require(ieugwasr)

## obtain summary stats and trait info
reg              <- paste0(chr.s, ":", pos.s,"-", pos.e)
res.trait        <- associations(reg, trait)
## just to make sure everything will work fine
res.trait        <- subset(res.trait, !is.na(beta))
## make unique, since some data processing error
res.trait        <- unique(res.trait)

## get meta information (careful: sometimes missed)
tr.info         <- tryCatch({
  gwasinfo(trait)
}, error=function(e){
  return(data.frame(trait=NA))
})
print(tr.info)

## recode INDELS
res.trait[, c("ea", "nea")] <- t(apply(res.trait[, c("ea", "nea")], 1, function(x){
  if(nchar(x[1])>1){
    return(c("I", "D"))
  }else if(nchar(x[2])>1){
    return(c("D", "I"))
  }else{
    return(x)
  }
}))

## create MarkerName to have a common ID to merge on
res.trait$MarkerName <- apply(res.trait[, c("chr", "position", "ea", "nea")], 1, function(x){
  return(paste0("chr", as.numeric(x[1]), ":", as.numeric(x[2]), "_", paste(sort(x[3:4]), collapse = "_")))
})

## transform to simple data frame
res.trait      <- data.frame(res.trait)

## omit some columns
res.trait$id   <- NULL
res.trait$rsid <- NULL

#-----------------------------------------#
##--       import genotype data        --##
#-----------------------------------------#

## read the dosage file
require(data.table)
tmp           <- fread(paste("tmp_input//tmp", soma, olink, trait, chr.s, pos.s, pos.e, "dosage", sep="."), sep=" ", header=T, data.table=F)
## transpose
rownames(tmp) <- tmp$rsid
## store allele information to realign effect estimates afterwards
tmp.info      <- tmp[, 1:6]
tmp           <- t(tmp[,-c(1:6)])
## retransform to data set (keep X in mind for variable names)
tmp           <- data.frame(ID_1=rownames(tmp), tmp)

## create another column to info to map names from the SNP data set
tmp.info$id         <- sapply(tmp.info$rsid, function(x) ifelse(substr(x, 1, 2) == "rs", x, paste0("X", gsub(":", ".", x))))
## edit some IDs (X-chromosome)
tmp.info$id         <- gsub("XX", "X", tmp.info$id)
tmp.info$id         <- gsub("XAffx-", "Affx.", tmp.info$id)
## create MarkerName column as well
tmp.info$MarkerName <- apply(tmp.info, 1, function(x){
  paste0("chr", as.numeric(x[1]), ":", as.numeric(x[4]), "_", paste(sort(x[5:6]), collapse = "_"))
})

## rename
pheno <- tmp
rm(tmp); gc()

################################################
####    combine everything and run coloc    ####
################################################

## create combined data set
res.protein               <- merge(res.soma, res.olink, by=c("rsid", "MarkerName", "chr", "pos"), suffixes = c(".soma", ".olink"))
## simplify
res.protein$Effect.olink  <- ifelse(res.protein$Allele1.soma == res.protein$Allele1.olink, res.protein$Effect.olink, -res.protein$Effect.olink)
## drop what is no longer needed
res.protein$Allele1.olink <- res.protein$Allele2.olink <- res.protein$Freq1.olink <- NULL
## edit some names
names(res.protein)        <- gsub("Allele1.soma", "Allele1", names(res.protein))
names(res.protein)        <- gsub("Allele2.soma", "Allele2", names(res.protein))
names(res.protein)        <- gsub("Freq1.soma", "Freq1", names(res.protein))

## add the trait data
res.all                   <- merge(res.protein, res.trait, 
                                   by.x=c("MarkerName", "chr", "pos"), 
                                   by.y=c("MarkerName", "chr", "position"))

## remove non-biallelelic variants
ii                        <- table(res.all$rsid)
res.all                   <- subset(res.all, rsid %in% names(ii[ii==1]))

## align effect estimates
res.all$beta.aligned      <- ifelse(res.all$Allele1 == res.all$ea, res.all$beta, -res.all$beta)

## add to ease mapping of LD matrix
res.all                   <- merge(res.all, tmp.info[, c("rsid", "id")], by="rsid", suffix=c(".trait", ".snp"))

## order by position
res.all                   <- res.all[order(res.all$pos),]

#-----------------------------------------#
##--       import genotype data        --##
#-----------------------------------------#

require(coloc)

res.coloc <- lapply(c("soma", "olink"), function(x){
 
  #-----------------------------------------#
  ##-- 	         sanity check            --##
  #-----------------------------------------#
  
  ## to ease some typing
  b     <- paste("Effect", x, sep=".")
  s     <- paste("StdErr", x, sep=".")
  
  print(c(b,s))
  
  ## top signal for the protein in all
  it    <- res.all$MarkerName[which.max(abs(res.all[, b]/res.all[, s]))]
  ## top signal SOMAscan
  if(x == "soma"){
    is    <- res.soma$MarkerName[which.max(abs(res.soma$Effect/res.soma$StdErr))] 
  }else{
    ## top signal Olink
    is    <- res.olink$MarkerName[which.max(abs(res.olink$Effect/res.olink$StdErr))] 
  }
            
  ## get the top SNP for the outcome
  io    <- res.all$MarkerName[which.max(abs(res.all$beta/res.all$se))]
  
  ## keep names
  isnps <- sapply(c(it, is, io), function(x){
    if(x %in% tmp.info$MarkerName){
      tmp.info$id[which(tmp.info$MarkerName == x)]
    }else{
      return(NA)
    }
  })
  
  print(isnps)
  
  ## get the LDs for each of the hits
  if(!is.na(isnps[2])){
    is <- cor(pheno[, isnps[1]], pheno[, isnps[2]])^2
  }else{
    is <- NA
  }
  
  io <- cor(pheno[, isnps[1]], pheno[, isnps[3]])^2
  
  #-----------------------------------------#
  ##-- 	            run coloc            --##
  #-----------------------------------------#
  
  ## prepare input
  D1          <- list(beta=res.all[, b], varbeta=res.all[, s]^2, 
                      type="quant", 
                      N=max(res.all[, paste("TotalSampleSize", x, sep=".")], na.rm=T), 
                      sdY=1,
                      MAF=res.all$MAF,
                      snp=res.all$id,
                      position=1:nrow(res.all))
  
  ## try out
  if(!("ncase" %in% names(tr.info))){
    print("use quantitative")
    
    ## define N
    if(!is.na(tr.info$trait)){
      n <- tr.info$sample_size
    }else{
      n <- max(res.all$n, na.rm=T)
    }
    
    D2          <- list(beta=res.all$beta.aligned, varbeta=res.all$se^2, 
                        type="quant", 
                        N=n,
                        sdY=1,
                        MAF=res.all$MAF,
                        snp=res.all$id,
                        position=1:nrow(res.all))
  }else{
    
    ## binary outcome
    D2          <- list(beta=res.all$beta.aligned, varbeta=res.all$se^2, 
                        type="cc",
                        s=tr.info$ncase/(tr.info$ncontrol+tr.info$ncase), 
                        N=tr.info$sample_size,
                        MAF=res.all$MAF,
                        snp=res.all$id,
                        position=1:nrow(res.all))
  }
  
  ## do naive coloc as well
  naive.coloc <- coloc.signals(D1, D2, method="single", p12=1e-6)
  
  ## get the strongest pQTL
  ii          <- which.max(abs(res.all[, b]/res.all[, s]))
  
  ## create output
  res         <- data.frame(platform = x, 
                            naive.coloc$summary,
                            ld.check.top  = io,
                            ld.check.sens = is,
                            res.all[ii, c("MarkerName", "rsid", "chr", "pos", "Freq1", "Allele1", "Allele2", grep(x, names(res.all), value=T), "beta.aligned", "se", "n", "p")])
  ## edit names to ease merging afterwards
  names(res)  <- gsub(paste0("\\.", x), "", names(res))
  ## return
  return(res)
  
})
res.coloc <- do.call(rbind, res.coloc)

## reshape and keep what is needed
res.coloc <- reshape(res.coloc, idvar = "chr", timevar = "platform", direction="wide")

## add some last coloums to ease mapping
res.coloc$pheno.soma  <- soma
res.coloc$pheno.olink <- olink
res.coloc$id.ieu      <- trait
# res.coloc$trait       <- tr.info$trait
res.coloc$trait       <- res.all$trait[1]

## write to file
write.table(res.coloc, paste("output/results", soma, olink, trait, chr.s, pos.s, pos.e, sep="."), sep="\t", row.names=F)

#-----------------------------------------#
##--          stacked locus plot       --##
#-----------------------------------------#

## create simple label data set
lab <- data.frame(name=c("soma", "olink", "trait"),
                  label=c(paste0("SomaScan (", soma, ")"), 
                          paste0("Olink (", olink, ")"), 
                          tr.info$trait))

source("scripts/plot_hyprcoloc.R")

png(paste0("graphics/", soma, ".", olink, ".", trait, ".", chr.s, ".", pos.s, ".", pos.e, ".png"), width=8, height=8, units="cm", res=200)
par(mar=c(.1,1.5,.75,.5), cex.axis=.5, cex.lab=.6, bty="l", tck=-.01, mgp=c(.6,0,0), xaxs="i", lwd=.5)
layout(matrix(1:4,4,1), heights = c(.3,.3,.3,.15))
plot.hyprcoloc(res.all, res.coloc, lab, pheno, tmp.info)
dev.off()


