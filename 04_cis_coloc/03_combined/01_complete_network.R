##############################################
#### fused protein - phenotype network    ####
#### Maik Pietzner             06/05/2021 ####
##############################################

rm(list=ls())
setwd("/rds/project/rjh234/rds-rjh234-mrc-epid/Studies/People/Maik/compare_SomaLogic_Olink/14_protein_net/data/")
options(stringsAsFactors = F)
load(".RData")

#############################################
#### import fused SomaScan/Olink network #### 
#############################################

## import fused coloc network
fused.network              <- read.table("Fused.cis.protein.phenotype.results.SomaScan.Olink.20210505.txt", sep = "\t", header = T)

## define consistent regional defintions for coloc
fused.network$chr          <- ifelse(!is.na(fused.network$chr.soma), fused.network$chr.soma, fused.network$chr.olink)
fused.network$region_start <- ifelse(!is.na(fused.network$region_start.soma), fused.network$region_start.soma, fused.network$region_start.olink)
fused.network$region_end   <- ifelse(!is.na(fused.network$region_end.soma), fused.network$region_end.soma, fused.network$region_end.olink)


## write to file to run jobs
write.table(fused.network[, c("pheno.soma", "pheno.olink", "id.ieu", "chr", "region_start", "region_end")], "cis.coloc.phenotypes.SomaScan.Olink.txt",
            sep="\t", col.names=F, row.names=F, quote=F)

## some are redundant
fused.network$id <- apply(fused.network[, c("pheno.soma", "pheno.olink", "id.ieu", "chr", "region_start", "region_end")], 1, function(x){
  ## generate file name
  return(paste(x[1], x[2], x[3], as.numeric(x[4]), as.numeric(x[5]), as.numeric(x[6]), sep="."))
})

#############################################
####            import results           #### 
#############################################

ii <- dir("../output/")

## import the results
res.net <- lapply(ii, function(x){
  ## import table
  res              <- read.table(paste0("../output/", x), sep="\t", header=T)
  ## add regional boundaries
  x                <- strsplit(x, "\\.")[[1]]
  res$region_start <- as.numeric(x[6])
  res$region_end   <- as.numeric(x[7])
  return(res)
})
res.net <- do.call(rbind, res.net)

#############################################
####       prepare clean output file     #### 
#############################################

## recode some alleles
res.net$Allele1.soma  <- gsub("TRUE", "T", res.net$Allele1.soma)
res.net$Allele2.soma  <- gsub("TRUE", "T", res.net$Allele2.soma)
res.net$Allele1.olink <- gsub("TRUE", "T", res.net$Allele1.olink)
res.net$Allele2.olink <- gsub("TRUE", "T", res.net$Allele2.olink)

## edit some names
names(res.net)        <- gsub("beta", "Effect.trait", names(res.net))
names(res.net)        <- gsub("^se", "StdErr.trait", names(res.net))
names(res.net)        <- gsub("^p\\.", "Pvalue.trait.", names(res.net))

## recode everything to the protein increasing allele
res.net$Effect.trait.soma  <- sign(res.net$Effect.soma)*res.net$Effect.trait.aligned.soma
res.net$EA.soma            <- ifelse(res.net$Effect.soma > 0, res.net$Allele1.soma, res.net$Allele2.soma)
res.net$NEA.soma           <- ifelse(res.net$Effect.soma > 0, res.net$Allele2.soma, res.net$Allele1.soma)
res.net$Effect.soma        <- abs(res.net$Effect.soma)
## same for Olink
res.net$Effect.trait.olink <- sign(res.net$Effect.olink)*res.net$Effect.trait.aligned.olink
res.net$EA.olink           <- ifelse(res.net$Effect.olink > 0, res.net$Allele1.olink, res.net$Allele2.olink)
res.net$NEA.olink          <- ifelse(res.net$Effect.olink > 0, res.net$Allele2.olink, res.net$Allele1.olink)
res.net$Effect.olink       <- abs(res.net$Effect.olink)

## keep relevant information
res.net                    <- res.net[, c("pheno.soma", "pheno.olink", "chr", "region_start", "region_end", "trait", "id.ieu", "PP.H4.abf.soma", "ld.check.top.soma", "ld.check.sens.soma", "MarkerName.soma",
                                          "EA.soma", "NEA.soma", "Effect.soma", "StdErr.soma", "Pvalue.soma", "Effect.trait.soma", "StdErr.trait.soma",
                                          "Pvalue.trait.soma",  "PP.H4.abf.olink", "ld.check.top.olink", "ld.check.sens.olink", "MarkerName.olink",
                                          "EA.olink", "NEA.olink", "Effect.olink", "StdErr.olink", "Pvalue.olink", "Effect.trait.olink", "StdErr.trait.olink",
                                          "Pvalue.trait.olink")]
## write to file
write.table(res.net, "Fused.cis.protein.phenotype.results.SomaScan.Olink.20210507.filled.txt", sep="\t", row.names=F)
