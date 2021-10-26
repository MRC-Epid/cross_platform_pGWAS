################################################
#### naive cis-pQTL coloc PheWAS            ####
#### Maik Pietzner               11/11/2020 ####
################################################

rm(list=ls())
setwd("/rds/project/rjh234/rds-rjh234-mrc-epid/Studies/People/Maik/compare_SomaLogic_Olink/09_coloc_cis_phewas/input/")
options(stringsAsFactors = F)
load(".RData")

##############################################
####          load Olink targets          ####
##############################################

## import prep from COVID work
olink.targets       <- data.frame(readxl::read_excel("../../01_pQTLs_Fenland/input/Olink.mapping.genomic.position.xlsx", 1))

## add GWAS id
olink.targets$pheno <- paste0("invn_res_X", tolower(olink.targets$varname))

## check chromosome definition
table(olink.targets$chromosome_name)

## define regions
olink.targets$region_start <- as.numeric(olink.targets$start_position) - 5e5
olink.targets$region_end   <- as.numeric(olink.targets$end_position) + 5e5
## edit some negative boundaries
olink.targets$region_start <- ifelse(olink.targets$region_start < 0, 0, olink.targets$region_start)

## make chromosome numeric
olink.targets$chr          <- as.numeric(olink.targets$chromosome_name)

## check region width
olink.targets$region.width <- olink.targets$region_end - olink.targets$region_start

##############################################
####        create input for coloc        ####
##############################################

## subset to promising regions
coloc.regions        <- subset(olink.targets, !is.na(chr))

## create input file
write.table(coloc.regions[, c("pheno", "chr", "region_start", "region_end")], "coloc.regions.txt", sep="\t", row.names=F, col.names=F, quote=F)

##############################################
####            collect output            ####
##############################################

ii        <- dir("../output/") ## 1 region missing
res.coloc <- lapply(ii, function(x){
  ## read results
  tmp              <- read.table(paste0("../output/", x), sep="\t", header=T)
  ## add some identifier
  tmp$region_start <- as.numeric(strsplit(x, "\\.")[[1]][5]) 
  tmp$region_end   <- as.numeric(strsplit(x, "\\.")[[1]][6])
  tmp$chr          <- as.numeric(strsplit(x, "\\.")[[1]][4])
  ## indicate no coloc
  if(is.na(tmp$id.ieu[1])){
    tmp$id.ieu <- "no evidence"
  }
  return(tmp)
})
require(plyr)
res.coloc <- do.call(rbind.fill, res.coloc)

## combine to see results
res.coloc <- merge(coloc.regions, res.coloc, by=c("pheno", "chr", "region_start", "region_end"), all = T, suffixes = c(".discovery", ".olink"))

##############################################
####            filter results            ####
##############################################

## drop missing runs
res.coloc           <- subset(res.coloc, !is.na(id.ieu) & id.ieu != "no evidence")
## drop bad stats (missing top protein hit)
res.coloc           <- subset(res.coloc, ld.check.sens > .9 & ld.lead.soma > .8)

## flag interesting results
res.coloc$PP.cand   <- ifelse(res.coloc$PP.H4.abf > .8, 1, 0)
res.coloc$LD.cand   <- ifelse(res.coloc$ld.check.top > .8, 1, 0)

##############################################
####           All for SCALLOP            ####
##############################################

## import SCALLOP results (do not care about cis/trans here)
res.cvd1 <- read.table("../../05_variant_annotation/output/Regional.sentinels.pGWAS.SCALLOP.CVD1.annotated.20200907.txt", sep="\t", header=T)

## create input for coloc script
coloc.cvd1              <- unique(res.cvd1[, c("rsid", "MarkerName", "Olink_varname", "Protein", "chr", "pos")])
## add 500kb regions to signals
coloc.cvd1$region_start <- coloc.cvd1$pos - 5e5
coloc.cvd1$region_end   <- coloc.cvd1$pos + 5e5

## write data to file
write.table(coloc.cvd1[, c("Protein", "chr", "region_start", "region_end")], "coloc.regions.scallop.txt", sep="\t", row.names=F, col.names=F, quote=F)

##############################################
####          collect results             ####
##############################################

ii             <- dir("../output_scallop//") ## 1 region missing
res.coloc.cvd1 <- lapply(ii, function(x){
  ## read results
  tmp              <- read.table(paste0("../output_scallop//", x), sep="\t", header=T)
  ## add some identifier
  tmp$region_start <- as.numeric(strsplit(x, "\\.")[[1]][5]) 
  tmp$region_end   <- as.numeric(strsplit(x, "\\.")[[1]][6])
  tmp$chr          <- as.numeric(strsplit(x, "\\.")[[1]][4])
  ## indicate no coloc
  if(is.na(tmp$id.ieu[1])){
    tmp$id.ieu <- "no evidence"
  }
  return(tmp)
})
require(plyr)
res.coloc.cvd1 <- do.call(rbind.fill, res.coloc.cvd1)

## combine to see results
res.coloc.cvd1 <- merge(coloc.cvd1, res.coloc.cvd1, 
                        by.x=c("Protein", "chr", "region_start", "region_end"), 
                        by.y=c("pheno", "chr", "region_start", "region_end"),
                        all = T, suffixes = c(".discovery", ".scallop"))

#-------------------------------------#
##--          filter results       --##
#-------------------------------------#

## drop missing runs
res.coloc.cvd1           <- subset(res.coloc.cvd1, !is.na(id.ieu) & id.ieu != "no evidence")
## drop bad stats (missing top protein hit)
res.coloc.cvd1           <- subset(res.coloc.cvd1, ld.check.sens > .9 & ld.lead.soma > .8)

## flag interesting results
res.coloc.cvd1$PP.cand   <- ifelse(res.coloc.cvd1$PP.H4.abf > .8, 1, 0)
res.coloc.cvd1$LD.cand   <- ifelse(res.coloc.cvd1$ld.check.top > .8, 1, 0)
