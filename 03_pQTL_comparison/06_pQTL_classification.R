####################################################
#### try to build a classifier for pQTLs        ####
#### Maik Pietzner                   17/11/2020 ####
####################################################

rm(list=ls())
setwd("v:/Programme1_DiabetesAetiology/People/Maik/Cross-platform proteomics/04_classifier_consistency/data/")
options(stringsAsFactors = F)
load(".RData")

##############################################
####          load results files          ####
##############################################

## import annotation from genetic results
dat.fenland            <- read.table("../../02_collate_results_genetics/data/pQTL.classification.set.Olink.20201117.txt", sep="\t", header=T)
dat.scallop            <- read.table("../../02_collate_results_genetics/data/pQTL.classification.set.SCALLOP.20201118.txt", sep="\t", header=T)

## slight edition to dat.scallop
dat.scallop$cis_trans  <- ifelse(dat.scallop$cis_trans == "cis", "cis", "trans")

## add a mapping for few slightly different rsIDs
tmp                    <- read.table("../../02_collate_results_genetics/data/mapping.Olink.conditional.SomaScan.Fenland.txt", sep="\t", header=T)
## add to the data
dat.fenland            <- merge(dat.fenland, tmp, by.x="rsid", by.y="rsid.cond", all.x = T)
dat.fenland$rsid.olink <- ifelse(!is.na(dat.fenland$rsid.olink), dat.fenland$rsid.olink, dat.fenland$rsid)
## rename
names(dat.fenland)     <- gsub("rsid\\.olink", "rsid.fused", names(dat.fenland))

##############################################
####     import phenotypic follow-up      ####
##############################################

## --> GWAS catalog data <-- ##
res.gwas.soma  <- read.table("Lookup.GWAS.catalogue.all.pQTLs.20200814.txt", sep="\t", header=T)
res.gwas.olink <- read.table("GWAS.catalog.linkage.pQTLs.Olink.SCALLOP.txt", sep="\t", header=T)

## add to the data 
dat.fenland$pheno.GWAS <- ifelse(dat.fenland$rsid.fused %in% res.gwas.soma$rsid_query | dat.fenland$rsid.fused %in% res.gwas.olink$rsid_query, "yes", "no")
## how many phenotypes
dat.fenland$num.GWAS   <- sapply(dat.fenland$rsid.fused, function(x){
  ii <- subset(res.gwas.soma, rsid_query == x)
  io <- subset(res.gwas.olink, rsid_query == x)
  return(max(c(nrow(ii),nrow(io)), na.rm=T))
})

## add to the data 
dat.scallop$pheno.GWAS <- ifelse(dat.scallop$rsid %in% res.gwas.soma$rsid_query | dat.scallop$rsid %in% res.gwas.olink$rsid_query, "yes", "no")
## how many phenotypes
dat.scallop$num.GWAS   <- sapply(dat.scallop$rsid, function(x){
  ii <- subset(res.gwas.soma, rsid_query == x)
  io <- subset(res.gwas.olink, rsid_query == x)
  return(max(c(nrow(ii),nrow(io)), na.rm=T))
})

## --> IEU data <-- ##

## import the results from coloc for Olink in Fenland
res.ieu.cis.olink <- read.table("Results.naive.coloc.cis.regions.Olink.txt", sep="\t", header=T)
res.ieu.cis.olink <- subset(res.ieu.cis.olink, PP.cand == 1)
## discard eQTL data
res.ieu.cis.olink <- subset(res.ieu.cis.olink, substr(id.ieu, 1, 4) != "eqtl")

## import the results from SCALLOP (cis and trans)
res.ieu.scallop   <- read.table("Results.naive.coloc.all.regions.SCALLOP.txt", sep="\t", header=T)
res.ieu.scallop   <- subset(res.ieu.scallop, PP.cand == 1)
## discard eQTL data
res.ieu.scallop   <- subset(res.ieu.scallop, substr(id.ieu, 1, 4) != "eqtl")

## cis SomaScan
res.ieu.cis.soma  <- read.table("../../../pGWAS SomaLogic/07_PheWAS_v2//data/Results.naive.coloc.cis.regions.SomaLogic.txt", sep="\t", header=T)
res.ieu.cis.soma  <- subset(res.ieu.cis.soma, PP.cand == 1)
## discard eQTL data
res.ieu.cis.soma  <- subset(res.ieu.cis.soma, substr(id.ieu, 1, 4) != "eqtl")

## trans SomaScan overlapping and excluding MHC region
res.ieu.trans.soma <- read.table("Results.naive.coloc.trans.regions.SOMAscan.txt", sep="\t", header=T)
res.ieu.trans.soma <- subset(res.ieu.trans.soma, PP.cand == 1)
## discard eQTL data
res.ieu.trans.soma <- subset(res.ieu.trans.soma, substr(id.ieu, 1, 4) != "eqtl")

## add to the data (careful can't be done on variant basis, since conditional SomaScan signals do not map 1:1 to lead Olink signals)
dat.fenland$pheno.IEU <- ifelse(dat.fenland$rsid.fused %in% c(res.ieu.cis.olink$rsid.lead, res.ieu.cis.soma$rsid.lead, res.ieu.trans.soma$rsid.lead), "yes", "no")
## how many phenotypes
dat.fenland$num.IEU   <- apply(dat.fenland[, c("rsid.fused", "pheno.soma", "pheno.olink")], 1, function(x){
  ii <- subset(res.ieu.cis.olink, rsid.lead == x[1] & pheno == x[3])
  io <- subset(res.ieu.cis.soma, rsid.lead == x[1] & pheno == x[2])
  it <- subset(res.ieu.trans.soma, rsid.lead == x[1] & pheno.soma == x[2])
  return(max(c(nrow(ii), nrow(io), nrow(it)), na.rm=T))
})

## add to the data 
dat.scallop$pheno.IEU <- ifelse(dat.scallop$rsid %in% c(res.ieu.scallop$rsid.discovery, res.ieu.cis.soma$rsid.lead, res.ieu.trans.soma$rsid.lead), "yes", "no")
## how many phenotypes
dat.scallop$num.IEU   <- apply(dat.scallop[, c("rsid", "pheno.soma", "pheno.scallop")], 1, function(x){
  ii <- subset(res.ieu.scallop, rsid.lead == x[1] & Protein == x[3])
  io <- subset(res.ieu.cis.soma, rsid.lead == x[1] & pheno == x[2])
  it <- subset(res.ieu.trans.soma, rsid.lead == x[1] & pheno.soma == x[2])
  return(max(c(nrow(ii), nrow(io), nrow(it)), na.rm=T))
})


##############################################
####           import GTEx-results        ####
##############################################

## --> import eQTL data <-- ##
res.eqtl.soma  <- read.table("../../../pGWAS SomaLogic/08_GTEx/data/GTEx.eQTL.mapping.SomaScan.txt", sep="\t", header=T)
res.eqtl.olink <- read.table("GTEx.eQTL.mapping.Olink.SCALLOP.txt", sep="\t", header=T)

## add cis-eQTL mapping
dat.fenland$cis.eqtl <- apply(dat.fenland[, c("rsid.fused", "pheno.soma", "Olink_varname")], 1, function(x){
  ## test for mapping eQTLs
  ii <- subset(res.eqtl.soma, rsid == x[1] & pheno == x[2] & cis_eQTL == "yes")
  io <- subset(res.eqtl.olink, rsid == x[1] & Olink_varname == x[3] & cis_eQTL == "yes")
  ## report back
  if(nrow(ii) > 0 | nrow(io) > 0){
    return("yes")
  }else{
    return("no")
  }
})

## add any eQTL
dat.fenland$any.eqtl <- apply(dat.fenland[, c("rsid.fused", "pheno.soma", "Olink_varname")], 1, function(x){
  ## test for mapping eQTLs
  ii <- subset(res.eqtl.soma, rsid == x[1] & pheno == x[2] & any_eQTL == "yes")
  io <- subset(res.eqtl.olink, rsid == x[1] & Olink_varname == x[3] & any_eQTL == "yes")
  ## report back
  if(nrow(ii) > 0 | nrow(io) > 0){
    return("yes")
  }else{
    return("no")
  }
})

## add cis-eQTL mapping
dat.scallop$cis.eqtl <- apply(dat.scallop[, c("rsid", "pheno.soma", "Olink_varname")], 1, function(x){
  ## test for mapping eQTLs
  ii <- subset(res.eqtl.soma, rsid == x[1] & pheno == x[2] & cis_eQTL == "yes")
  io <- subset(res.eqtl.olink, rsid == x[1] & Olink_varname == x[3] & cis_eQTL == "yes")
  ## report back
  if(nrow(ii) > 0 | nrow(io) > 0){
    return("yes")
  }else{
    return("no")
  }
})

## add any eQTL
dat.scallop$any.eqtl <- apply(dat.scallop[, c("rsid", "pheno.soma", "Olink_varname")], 1, function(x){
  ## test for mapping eQTLs
  ii <- subset(res.eqtl.soma, rsid == x[1] & pheno == x[2] & any_eQTL == "yes")
  io <- subset(res.eqtl.olink, rsid == x[1] & Olink_varname == x[3] & any_eQTL == "yes")
  ## report back
  if(nrow(ii) > 0 | nrow(io) > 0){
    return("yes")
  }else{
    return("no")
  }
})

##############################################
####          additional features         ####
##############################################

## incorporate the same features as in the observational correlation
tmp.feat             <- read.table("../../01_observational_correlations/data/Meta.data.protein.targets.SomaScan.Olink.txt", sep="\t", header=T)

## edit phenotype names for Olink
tmp.feat$pheno.olink <- gsub("invn_res_Xtdgf1_nex", "invn_res_Xtdgf1_nex_ipc", tmp.feat$pheno.olink)
tmp.feat$pheno.olink <- gsub("invn_res_Xx_5__nt_onc2", "invn_res_X_5__nt_onc2", tmp.feat$pheno.olink)

## add to the classification data
dat.fenland          <- merge(dat.fenland, tmp.feat[, -13])
dat.scallop          <- merge(dat.scallop, tmp.feat[,-13])

##############################################
####           univariate analysis        ####
##############################################

## import feature mapping
feat          <- read.table("features.classification.txt", sep="\t", header = T)
## susbet to those available so far
feat          <- subset(feat, variable_id %in% names(dat.fenland))
## count numbers
feat$fen.miss <- sapply(feat$variable_id, function(x) sum(is.na(dat.fenland[,x]))/nrow(dat.fenland))
feat$sca.miss <- sapply(feat$variable_id, function(x) sum(is.na(dat.scallop[,x]))/nrow(dat.scallop))

#----------------------#
##-- create factors --##
#----------------------#

## recode factors
for(j in feat$variable_id[which(feat$type == "factor")]){
  dat.fenland[,j] <- as.factor(dat.fenland[,j])
  ## set reference
  dat.fenland[,j] <- relevel(dat.fenland[,j], ref=feat$ref[which(feat$variable_id == j)])
}
## same for outcome
dat.fenland$consistency <- as.factor(dat.fenland$consistency)
dat.fenland$consistency <- relevel(dat.fenland$consistency, ref="inconsistent")

## recode factors
for(j in feat$variable_id[which(feat$type == "factor")]){
  dat.scallop[,j] <- as.factor(dat.scallop[,j])
  ## set reference
  dat.scallop[,j] <- relevel(dat.scallop[,j], ref=feat$ref[which(feat$variable_id == j)])
}
## same for outcome
dat.scallop$consistency <- as.factor(dat.scallop$consistency)
dat.scallop$consistency <- relevel(dat.scallop$consistency, ref="inconsistent")

#------------------------------#
##-- transform some outcome --##
#------------------------------#

## transform selected outcomes
for(j in feat$variable_id[which(feat$transform == "sqrt")]){
  ## define transformed variable
  dat.fenland[,j] <- sqrt(dat.fenland[,j])
  ## define transformed variable
  dat.scallop[,j] <- sqrt(dat.scallop[,j])
}

#----------------------#
##--  edit outcome  --##
#----------------------#


## import newly stratified reason
tmp.lab            <- read.table("fused.reaons.txt", sep="\t", header=T)

## edit SCALLOP to map better to fused reasons
dat.scallop$reason <- gsub("unique to SCALLOP", "unique to Olink", dat.scallop$reason)

## add to the data
dat.fenland        <- merge(dat.fenland, tmp.lab)
dat.scallop        <- merge(dat.scallop, tmp.lab)

## prune some as reason for condiotional SNPs differs
dat.fenland$fused.reason <- ifelse(dat.fenland$consistency == "consistent", "consistent", dat.fenland$fused.reason) 
dat.scallop$fused.reason <- ifelse(dat.scallop$consistency == "consistent", "consistent", dat.scallop$fused.reason) 
table(dat.fenland$fused.reason)
table(dat.scallop$fused.reason)

## recode as factors
dat.fenland$fused.reason <- as.factor(dat.fenland$fused.reason) 
dat.fenland$fused.reason <- relevel(dat.fenland$fused.reason, ref = "consistent") 
dat.scallop$fused.reason <- as.factor(dat.scallop$fused.reason) 
dat.scallop$fused.reason <- relevel(dat.scallop$fused.reason, ref = "consistent") 

##############################################
####          stratify by reason          ####
##############################################

## load function to do the testing
source("../scripts/simple_logistic_regression.R")

#--------------------------------------#
##--            Fenland             --##
#--------------------------------------#

## run for each reason separately
res.fenland <- lapply(c("distinct", "cond_snp", "direct_incon", "unique to Olink", "unique to SOMAscan"), function(x){
  ## define new reference level to align to consistency
  dat.fenland$fused.reason <- relevel(dat.fenland$fused.reason, ref=x)
  res                      <- run.logisitc(subset(dat.fenland, fused.reason %in% c(x,"consistent")), feat$variable_id, "fused.reason")
  res$reason               <- x
  return(res)
})
res.fenland <- do.call(rbind, res.fenland)

#--------------------------------------#
##--            SCALLOP             --##
#--------------------------------------#

## run for each reason separately
res.scallop <- lapply(c("distinct", "direct_incon", "unique to Olink", "unique to SOMAscan"), function(x){
  ## define new reference level to align to consistency
  dat.scallop$fused.reason <- relevel(dat.scallop$fused.reason, ref=x)
  res                      <- run.logisitc(subset(dat.scallop, fused.reason %in% c(x,"consistent")), feat$variable_id, "fused.reason")
  res$reason               <- x
  return(res)
})
res.scallop <- do.call(rbind, res.scallop)

##############################################
####             create table             ####
##############################################

## prepare label for plotting
feat$srt <- 1:nrow(feat)

## add label
res.fenland               <- merge(res.fenland, feat, by.x = "exposure", by.y= "variable_id")
## sort again
res.fenland               <- res.fenland[order(res.fenland$srt, res.fenland$level),]
## create new plotting label
res.fenland$plot_label    <- apply(res.fenland[, c("label", "level")], 1, function(x){
  if(x[2] == ""){
    return(x[1])
  }else{
    ## add reference level compared to reference
    paste(x[1], "-", x[2])
  }
})

## add label
res.scallop               <- merge(res.scallop, feat, by.x = "exposure", by.y= "variable_id")
## sort again
res.scallop               <- res.scallop[order(res.scallop$srt, res.scallop$level),]
## create new plotting label
res.scallop$plot_label    <- apply(res.scallop[, c("label", "level")], 1, function(x){
  if(x[2] == ""){
    return(x[1])
  }else{
    ## add reference level compared to reference
    paste(x[1], "-", x[2])
  }
})


## create column for Odds ratio
res.fenland$OR.col <- paste0(sprintf("%.2f", exp(res.fenland$beta)), " (", sprintf("%.2f", exp(res.fenland$beta - 1.96*res.fenland$se)),
                             "; ", sprintf("%.2f", exp(res.fenland$beta + 1.96*res.fenland$se)), ")") 

## now SCALLOP
res.scallop$OR.col <- paste0(sprintf("%.2f", exp(res.scallop$beta)), " (", sprintf("%.2f", exp(res.scallop$beta - 1.96*res.scallop$se)),
                             "; ", sprintf("%.2f", exp(res.scallop$beta + 1.96*res.scallop$se)), ")") 

## create common table
res.table          <- merge(subset(res.fenland, priority == 1), subset(res.scallop, priority == 1),
                            by=c("exposure", "reason", "level", "label", "type", "ref", "priority", "comment", "srt", "plot_label", "sca.miss", "fen.miss"),
                            suffixes = c(".fenland", ".scallop"), all.x = T) 
## reduce to the reasons of interest
res.table          <- subset(res.table, reason %in% c("distinct", "unique to Olink", "unique to SOMAscan"))

## order
res.table          <- res.table[order(res.table$srt),]

write.table(res.table, "Results.pQTL.consistency.factors.20210208.txt", sep="\t", row.names = F)