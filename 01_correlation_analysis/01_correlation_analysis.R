####################################################
#### Comparison SomaLogic and Olink             ####
#### Maik Pietzner                   06/08/2020 ####
####################################################

rm(list=ls())
setwd("v:/Programme1_DiabetesAetiology/People/Maik/Cross-platform proteomics/01_observational_correlations/data/")
options(stringsAsFactors = F)
load(".RData")

##################################################
#### load the data to compare measured values ####
##################################################

#--------------------#
##--     label    --##
#--------------------#

require(readxl)
label <- data.frame(read_excel("V:/Programme1_DiabetesAetiology/Data/SomaLogic/Annotations_Sep2019/merged_files/proteomics_druggablegenome_secreted_5Dec2019.xlsx", 1))

## edit two protein names
label$Olink_varname[which(label$Olink_varname == "_4E_BP1_INF")] <- "X_4E_BP1_INF"
label$Olink_varname[which(label$Olink_varname == "_5__NT_ONC2")] <- "X_5__NT_ONC2"


#--------------------#
##--   SomaLogic  --##
#--------------------#

## load the SomaLogic data
require(readstata13)
sl <- read.dta13("V:/Programme1_DiabetesAetiology/Data/Fenland/Somalogic/Fenland_main_03Jul2018/Cleaned_AMN_normalised_Oct2018/All_AMN_proteinAsAnalyzed-2018-10-04_R205.dta")

#--------------------#
##--    Olink     --##
#--------------------#

## load the data 
olink             <- read.csv("V:/Programme1_DiabetesAetiology/Data/Fenland/O-Link/including_below_LOD_cleaned_26Feb2019/Olink_including_below_LOD_cleaned_26Feb2019.csv")
## undo log2-transformation
ol.prots          <- unique(label$Olink_varname)
## delete one
ol.prots          <- sort(ol.prots)[-1]
olink[, ol.prots] <- 2^olink[, ol.prots]

## add identifier to map both
label$slid        <- paste("SeqId", gsub("-", "_", label$SeqId), sep="_")
sum(label$slid %in% names(sl)) ## N = 4,979

## what is in the overlap
tmp               <- subset(label, Olink_varname != "" & SeqId != "" & !is.na(version_v4))
## N = 937 (caveat, includes multiple somamers targeting the same olink protein)

#--------------------#
##--    combine   --##
#--------------------#

## merge both data sets
dat <- merge(sl[, c("SampleId", tmp$slid)], olink[, c("SampleID", unique(tmp$Olink_varname))], by.x="SampleId", by.y="SampleID")
## N = 485 samples

rm(sl); gc()

#####################################################
#### calculate different measures of correlation ####
#####################################################

## simple spearman, pearsons on raw, log10 and inverse normal data
res        <- tmp[, c("slid", "Olink_varname")]

## simple correlation
res[, 3:4] <- t(apply(res[, 1:2], 1, function(x){
  s <- cor(dat[, x[1]], dat[, x[2]], m="s", u="p")
  p <- cor(log2(dat[, x[1]]), log2(dat[, x[2]]), m="p", u="p")
  return(c(s,p))
}))

## assign names
names(res) <- c("slid", "Olink_varname", "s.raw", "p.log2")

## order by spearman correlation
res <- res[order(res$s.raw, decreasing = T), ]

## add annoation of targets
res <- merge(res, label[, c("slid", "Target", "TargetFullName", "EntrezGeneSymbol")], by="slid")

#####################################################
####  try to find factors explaining correlation ####
#####################################################

##-- missingness, panel measured, type of protein (characteristics -- enrichment using DAVID)
##-- genetics (use cis and trans-pQTLs from Ellies list)

## calculate vlaues below LOD for SOMAscan


## load meta data from Olink
olinkMetaInfo <- read.csv("V:/Programme1_DiabetesAetiology/Data/Fenland/O-Link/interim_cleaned_31Jan2019/Olink_protein_metadata.csv")
## rename two proteins
olinkMetaInfo$varname[which(olinkMetaInfo$varname == "_4E_BP1_INF")] <- "X_4E_BP1_INF"
olinkMetaInfo$varname[which(olinkMetaInfo$varname == "_5__NT_ONC2")] <- "X_5__NT_ONC2"

## add to the results from comparison
res           <- merge(res, olinkMetaInfo, by.x="Olink_varname", by.y="varname", all.x=T)
## two proteins missing

## import outlier measurement for Olink
require(readstata13)
tmp            <- read.dta13("V:/Programme1_DiabetesAetiology/Data/Fenland/O-Link/including_below_LOD_cleaned_26Feb2019/outliers_by_protein_19Jan2021.dta") 
## add postfix
names(tmp)[-1] <- paste0(names(tmp)[-1], ".olink") 
## edit some names to ease mapping
tmp$varname    <- gsub("_5__NT_ONC2", "X_5__NT_ONC2", tmp$varname) 
tmp$varname    <- gsub("_4E_BP1_INF", "X_4E_BP1_INF", tmp$varname) 
tmp$varname    <- gsub("IgG_Fc_receptor_II_b", "IgG_Fc_receptor_II_b_CAR2", tmp$varname) 

## add to the data
res            <- merge(res, tmp, by.x="Olink_varname", by.y="varname", all.x=T)

#------------------------------------#
##-- import LOD data for SOMAscan --##
#------------------------------------#

require(readstata13)
soma.lod <- read.dta13("V:/Programme1_DiabetesAetiology/Data/Fenland/Somalogic/Fenland_main_03Jul2018/Cleaned_ANML_normalised_Dec2019/V4-18-039.ANML.FINAL_somamer_QC_corrected16Jan.dta")

## add information to the results
res      <- merge(res, soma.lod[, c("MRC_seqid", "eLOD", "percent_below_lod", "percent_outliers", "ColCheck")], by.x="slid", by.y="MRC_seqid")

#------------------------------------#
##--      protein annotations     --##
#------------------------------------#

## import annotation from Julia for the proteins
require(readxl)
prots.annotation <- read.csv("V:/Programme1_DiabetesAetiology/People/Julia/PhD/Crossplatform_comparisons/data_output/protein_annotation_master_database.csv") 
## one entry looks odd
table(sapply(prots.annotation$Entry, nchar))
## some manual tweaking is needed
res$Entry        <- res$Uniprot_ID
res$Entry        <- gsub("O43521-2", "O43521", res$Entry)
res$Entry        <- gsub("Q8NEV9,Q14213", "Q8NEV9", res$Entry, fixed = T)
res$Entry        <- gsub("P29460,P29459", "P29459", res$Entry, fixed = T)
## two proteins are missing, one with no Uniprot_ID (NTpro_BNP and one w/o annotation in the data base FGF_5_INF)
grep("P29460|P29459", prots.annotation$Entry, value=T)

## create new data set for prediction
res.pred         <- merge(res, prots.annotation, by = "Entry", all.x = T)
res.pred$Mass    <- as.numeric(gsub(",", "", res.pred$Mass))

## combine with the data for prediction
res.pred         <- merge(res.pred, label[, c("slid", "ApparentKdM", "CharacterizationInfo", "MassSpecConfirmationinMatrix", "TotalCVPlasma", "Concentrationpgl")], 
                          by="slid", all.x = T)
## transform one measure
res.pred$ApparentKdM      <- -log10(res.pred$ApparentKdM)
res.pred$Concentrationpgl <- -log10(res.pred$Concentrationpgl)

## create list of features to be investigated
feat.pred <- c(names(prots.annotation)[-c(1,3,25,27,28)], "dilution", "panelshortname", "ApparentKdM", "TotalCVPlasma", "Concentrationpgl", "LOD", "Missing_data_percent",
               "percent_below_lod", "eLOD", "percent_outliers", "percent_outliers.olink")

#-----------------------------------------#
##--   test each feature separately    --##
#-----------------------------------------#

## use linear regression models
source("../scripts/run_regression_exposure.R")

## restict to numeric variables for now [sapply(feat.pred, function(x) is.numeric(res.pred[,x]))]
res.fac.cor        <- run.regression(res.pred, feat.pred, "s.raw", "")

## create dataset without missing values
tmp                <- na.omit(res.pred[, c("s.raw", feat.pred[-which(feat.pred %in% c("Concentrationpgl", "Human.secretome.annotation"))])])
tmp$dilution       <- as.factor(tmp$dilution)
tmp$panelshortname <- as.factor(tmp$panelshortname)

## use Boruta feature selection as an alternative
require(Boruta)

## apply to the set of features
res.boruta         <- Boruta(tmp[, -1], tmp[,"s.raw"], maxRuns = 500)
plot(res.boruta)
res.boruta         <- TentativeRoughFix(res.boruta)

## combine variable importance with p-Value graph and reporte achieved R2 in the paper
tmp                <- attStats(res.boruta)
tmp$exposure       <- rownames(tmp)
res.fac.cor        <- merge(res.fac.cor, tmp)

##################################################
####          correlation raw data            ####
##################################################

## obtain data prepared by Julia
soma.raw <- data.table::fread("Fenland_raw_SL_V4_proteome_measures.txt", sep="\t", header=T)
## keep only those overlapping with Olink
soma.raw <- subset(soma.raw, SampleId %in% dat.ivn$SampleId)

## add the data needed
soma.raw <- merge(dat.ivn, soma.raw, by="SampleId", suffixes = c(".norm", ".raw"))

## create new results data frame
res.raw  <- res

## add correlation coefficients
names(res.raw)             <- gsub("s\\.raw", "spearman.norm.olink", names(res.raw))
res.raw$spearman.raw.olink <- apply(res.raw[, 1:2], 1, function(x){
  cor(soma.raw[,paste0(x[2], ".raw")], soma.raw[,x[1]], m="s", u="p")
})
## raw vs. correlated
res.raw$spearman.raw.norm <- sapply(res.raw$slid, function(x){
  cor(soma.raw[, paste0(x, ".raw")], soma.raw[, paste0(x, ".norm")], m="s", u="p")
})

#---------------------------------------------------#
##--      Factors influencing the difference     --##
#---------------------------------------------------#

## add prediciton factors
ii                 <- names(res.pred)[!(names(res.pred) %in% names(res.raw))]
res.raw            <- merge(res.raw, res.pred[, c("Olink_varname", "slid", ii)])

## restict to numeric variables for now [sapply(feat.pred, function(x) is.numeric(res.pred[,x]))]
res.fac.raw        <- run.regression(res.raw, c(feat.pred, "beta.prot_PC", "absz.pc", "beta.prot_PC.RAW", "absz.pc.raw", grep("p_locus", names(tmp.trans), value=T)), "diff.cor", "")

## create a temporary data set with no missing values
tmp                <- na.omit(res.raw[, c("diff.cor", "beta.prot_PC", "absz.pc", "beta.prot_PC.RAW", "absz.pc.raw",
                                          feat.pred[-which(feat.pred %in% c("Concentrationpgl", "Human.secretome.annotation"))],
                                          grep("p_locus", names(tmp.trans), value=T))])
tmp$dilution       <- as.factor(tmp$dilution)
tmp$panelshortname <- as.factor(tmp$panelshortname)

## use Boruta feature selection as an alternative
require(Boruta)

## apply to the set of features
res.boruta.raw     <- Boruta(tmp[, -1], tmp[,"diff.cor"], maxRuns = 500)
plot(res.boruta.raw)
res.boruta.raw     <- TentativeRoughFix(res.boruta.raw)

## combine variable importance with p-Value graph and reporte achieved R2 in the paper
tmp                <- attStats(res.boruta.raw)
tmp$exposure       <- rownames(tmp)
res.fac.raw        <- merge(res.fac.raw, tmp)