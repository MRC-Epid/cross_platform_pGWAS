#################################################################
####  compare genetic estimates between SomaLogic and Olink  ####
####  Maik Pietzner                              26/08/2020  ####
#################################################################

rm(list=ls())
setwd("/rds/project/rjh234/rds-rjh234-mrc-epid/Studies/People/Maik/compare_SomaLogic_Olink/01_pQTLs_Fenland/input/")
options(stringsAsFactors = F)
load(".RData")

##############################################
####       load results and label         ####
##############################################

## results from SOMAscan discovery
res.soma  <- read.table("../../../pGWAS_SomaLogic/00_collect_results/data/Regional.sentinels.pGWAS.SomaLogic.20200807.txt", sep="\t", header=T)
cond.soma <- read.table("../../../pGWAS_SomaLogic/00_collect_results/data/Conditional.results.pGWAS.SomaLogic.20200807.txt", sep="\t", header=T)
## discard regional sentinels
cond.soma <- subset(cond.soma, hit != 1)
  
## results from Olink
res.olink     <- read.table("combined_olink_sentinel.txt", sep="\t", header=T)
## discard rare variants
res.olink$MAF <- ifelse(res.olink$af > .5, 1-res.olink$af, res.olink$af)
res.olink     <- subset(res.olink, MAF >= .01)
## N = 301 signals

## merge regions!

## sort by protein and position
res.olink <- res.olink[order(res.olink$pheno, res.olink$chr, res.olink$pos),]

## treat MHC region different
res.olink$signal.start[which(res.olink$chr == 6 & res.olink$pos >= 28477797 & res.olink$pos <= 33448354)] <- 28477797 
res.olink$signal.end[which(res.olink$chr == 6 & res.olink$pos >= 28477797 & res.olink$pos <= 33448354)]   <- 33448354 

## data frame to store merged results
source("../scripts/merge_regions.R")
res.olink <- merge.region(res.olink)
## 287

## --> add cis/trans annotation for OLINK <-- ##
require(readxl)
## has been edited no include only unique assignments (if possible absed on cis-pQTLs)
olink.genes         <- data.frame(read_excel("Olink.mapping.genomic.position.xlsx", 1))
## add identifier
olink.genes$pheno   <- paste0("invn_res_X", tolower(olink.genes$varname))
## add to the data
res.olink           <- merge(res.olink, olink.genes, by="pheno")
## add cis/trans information
res.olink$cis_trans <- ifelse(res.olink$chr != res.olink$chromosome_name, "trans",
                              ifelse(abs(res.olink$pos - as.numeric(res.olink$start_position)) <= 5e5, "cis", "trans"))
table(res.olink$cis_trans)
# cis trans 
# 246    41

## load label for mapping
require(readxl)
label       <- data.frame(read_excel("../../../COVID_19_targets/proteomics_druggablegenome_secreted_5Dec2019.xlsx"), sep="\t")
label       <- subset(label, !is.na(SomaId_v4))
label$pheno <- paste0("res_invn_X", gsub("-", "_", label$proteomicsid))

## import results from observational correlation to ease mapping
res.obs             <- read.table("SL.OL.observational.correlation.20200806.txt", sep="\t", header=T)
## add GWAS ID to ease mapping
res.obs$pheno.soma  <- gsub("SeqId_", "res_invn_X", res.obs$slid)
res.obs$pheno.olink <- paste0("invn_res_X", tolower(res.obs$Olink_varname))
## edit one entry
res.obs$pheno.olink <- gsub("invn_res_Xx_5__nt_onc2", "invn_res_X_5__nt_onc2", res.obs$pheno.olink)
res.obs$pheno.olink <- gsub("invn_res_Xtdgf1_nex", "invn_res_Xtdgf1_nex_ipc", res.obs$pheno.olink)

######################################
####     prepare look-up list     ####
######################################

## edit one name
res.soma$Olink_varname <- gsub("_5__NT_ONC2", "X_5__NT_ONC2", res.soma$Olink_varname)

## break down to overlapping targets
res.soma            <- merge(res.soma, res.obs, 
                             by.x=c("pheno", "Olink_varname", "Target", "TargetFullName", "EntrezGeneSymbol"), 
                             by.y=c("pheno.soma", "Olink_varname", "Target", "TargetFullName", "EntrezGeneSymbol")) 
## 1,796

## restrict to autosomes
res.soma            <- subset(res.soma, chr != 23)
## 1,780

## now Olink
res.olink           <- merge(res.olink, res.obs, by.x="pheno", by.y="pheno.olink") ## 264
## do some renaming to ease SNP-list
names(res.soma)[1]  <- "pheno.soma"
names(res.olink)[1] <- "pheno.olink"

## create common list to query results
query               <- unique(rbind(res.soma[, c("pheno.soma", "pheno.olink", "MarkerName")],
                                    res.olink[, c("pheno.soma", "pheno.olink", "MarkerName")])) 


## write to file to look up statistics
write.table(query, "Reciprocal.Look.up.SOMAscan.Olink.txt", sep="\t", row.names = F, col.names = F, quote = F)

## list of rsIDs, chr, pos for LD-files
tmp <- subset(cond.soma, pheno %in% res.soma$pheno.soma) ## 744
tmp <- unique(rbind(res.soma[, c("MarkerName", "chr", "pos", "rsid")],
                    tmp[, c("MarkerName", "chr", "pos", "rsid")],
                    res.olink[, c("MarkerName", "chr", "pos", "rsid")])) ## 2,038
write.table(tmp, "../../02_ld_files/input/Fenland.pQTLs.SomaLogic.Olink.txt", sep = "\t", row.names=F)

## --> create simple results table <-- ##
write.table(res.soma, "Results.SOMAscan.sentinel.pQTLs.Olink.overlap.20200826.txt", sep="\t", row.names=F)
write.table(subset(cond.soma, pheno %in% res.soma$pheno.soma), "Results.SOMAscan.secondary.pQTLs.Olink.overlap.20200826.txt", sep="\t", row.names=F)
write.table(res.olink, "Results.Olink.sentinel.pQTLs.SOMAscan.overlap.20200826.txt", sep="\t", row.names=F)

######################################
####        import look-up        ####
######################################

ii   <- dir("../output/")
## separate into SOMAscan and Olink
ii.s <- grep("res_invn_X[0-9]", ii, value=T) 
ii.o <- grep("res_invn_X[0-9]", ii, value=T, invert = T) 
## some are missing --> naming issue

## head for SOMAscan and Olink
s.names <- c("rsid", "MarkerName", "Allele1", "Allele2", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Effect", "StdErr", "P-value",
             "Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal", "TotalSampleSize", "chr", "pos")
## rename allele coding for Olink to match with coding from METAL output
o.names <- c("chr", "rsid", "pos", "Allele2", "Allele1", "Freq1", "info", "Effect", "StdErr", "tval", "log10p", "MarkerName")

#----------------------------#
##-- create combined file --##
#----------------------------#

## SOMAscan
ii.s <- lapply(ii.s, function(x){
  ## load the file
  tmp            <- read.table(paste0("../output/", x), sep="\t")
  ## assign names
  names(tmp)     <- s.names
  ## add protein
  tmp$pheno.soma <- strsplit(x, "\\.")[[1]][1]
  return(tmp)
})
ii.s <- do.call(rbind, ii.s)
## edit some names
names(ii.s)[-c(1:2,18:20)] <- paste0(names(ii.s)[-c(1:2,18:20)],".soma")

## Olink
ii.o <- lapply(ii.o, function(x){
  ## load the file (some files do not contain any lines)
  tryCatch({tmp <- read.table(paste0("../output/",x), sep=" ")
  ## assign names
  names(tmp)      <- o.names
  ## add protein
  tmp$pheno.olink <- strsplit(x, "\\.")[[1]][1]
  return(tmp)
  }, error=function(e){
    cat("no lines in", x, "\n")
    return(NA)
  })
})
ii.o <- do.call(rbind, ii.o)
ii.o <- na.omit(ii.o)

## edit some names
names(ii.o)[4:11] <- paste0(names(ii.o)[4:11], ".olink")

## create final look-up file
res.lookup <- merge(query, ii.s) 
res.lookup <- merge(res.lookup, ii.o, by=c("pheno.olink", "MarkerName"), suffixes = c(".soma", ".olink"))

## recode alleles
res.lookup$Allele1.olink[which(res.lookup$Allele1.olink == T)] <- "T"
res.lookup$Allele2.olink[which(res.lookup$Allele2.olink == T)] <- "T"

## align effect estimates from Olink
res.lookup$Effect.olink.aligned <- ifelse(toupper(res.lookup$Allele1.soma) == res.lookup$Allele1.olink, res.lookup$Effect.olink, -res.lookup$Effect.olink)

## indicate discovery cohort
res.lookup$discovery            <- sapply(res.lookup$MarkerName, function(x){
  if(x %in% res.soma$MarkerName & x %in% res.olink$MarkerName){
    return("both")
  }else if(x %in% res.soma$MarkerName & !(x %in% res.olink$MarkerName)){
    return("SOMAscan")
  }else{
    return("Olink")
  }
})
table(res.lookup$discovery)

## edit one name
names(res.lookup) <- gsub("rsid.soma", "rsid", names(res.lookup), fixed = T)
names(res.lookup) <- gsub("chr.soma", "chr", names(res.lookup), fixed = T)
names(res.lookup) <- gsub("pos.soma", "pos", names(res.lookup), fixed = T)

#--------------------------------------------------#
##--              apply LD-pruning              --##
#--------------------------------------------------#

## --> create LD matrix <-- ##

## files for LD
ii <- dir("../../02_ld_files//output/")
ii <- grep("\\.ld", ii, value=T)

## create LD-matrix for all SNPs
require(data.table)
ld.matrix <- lapply(ii, function(x){
  tmp        <- data.frame(fread(paste0("../../02_ld_files/output/", x), sep=" ", header=T))
  ## do two subsets to get only SNPs of interest
  tmp        <- subset(tmp, SNP_A %in% res.lookup$rsid)
  print(summary(tmp))
  return(tmp)
})
ld.matrix <- do.call(rbind, ld.matrix)

## generate list of SNPs in high LD
ex.snps <- subset(ld.matrix, R2 >= .8 & SNP_B %in% res.lookup$rsid)
ex.snps <- subset(ex.snps, SNP_A != SNP_B)
## 1724 high LD combinations

## ease mapping
ex.snps$SNP.1 <- apply(ex.snps[, c("SNP_A", "SNP_B")], 1, function(x){
  return(sort(x)[1])
})
ex.snps$SNP.2 <- apply(ex.snps[, c("SNP_A", "SNP_B")], 1, function(x){
  return(sort(x)[2])
})

## subset again
ex.snps <- unique(ex.snps[, c("SNP.1", "SNP.2", "R2")])

#----------------------------------------------#
##-- drop individual SNPs if they are in LD --##
##-- for a given pair of interest           --##
#----------------------------------------------#

## create a mapping identifier
res.lookup$tmp1 <- paste(res.lookup$pheno.soma, res.lookup$pheno.olink, sep="$")
## order the data
res.lookup      <- res.lookup[order(res.lookup$tmp1, res.lookup$chr, res.lookup$pos),]

## collect the possible duplicates
tmp <- lapply(unique(res.lookup$tmp1), function(x){
  
  ## get the sub data set
  tmp <- subset(res.lookup, tmp1 == x)
  
  ## do only with at least two SNPs
  if(nrow(tmp) > 1){
    ## get all pairwise combinations of SNPs (sort already)
    ii <- apply(combn(tmp$rsid, 2), 2, sort)
    ## test if those pairs are in the pruned LD list
    jj <- apply(ii, 2, function(x){ifelse(any(which(ex.snps$SNP.1 == x[1] & ex.snps$SNP.2 == x[2])), T, F)})
    ## select SNPs if present
    if(sum(jj) > 0){
      ## subset to pairs of interest
      ii <- ii[,jj,drop=F]
      for(k in 1:ncol(ii)){
        if (k == 1){
          kk <- data.frame(rsid.x = ii[1,k], rsid.y=ii[2,k], tmp1=x)
        }else{
          kk <- rbind(kk, data.frame(rsid.x = ii[1,k], rsid.y=ii[2,k], tmp1=x))
        }
      }
      return(kk)
    }else{
      return(NA)
    }
  }else{
    return(NA)
  }
})
tmp <- na.omit(do.call(rbind, tmp))
##-- 78 pairs in total --##

## do on dummy data frame
foo <- res.lookup

## go through the list
for(k in 1:nrow(tmp)){
  ## identify positions in the results data frame
  ii  <- which(foo$tmp1 == tmp[k,"tmp1"] & foo$rsid %in% unlist(tmp[k,c("rsid.x", "rsid.y")]))
  ## delete only, if both originate from different GWAS (ie. Olink and SomaLogic)
  if(nrow(foo[ii,])==2){
    # print(foo$cis_trans[ii])
    ## decide which has the lower z-score
    jj  <- which.min(abs(foo$Effect.soma[ii]/foo$StdErr.soma[ii]))
    ## define novel entry for platform column
    foo$discovery[ii] <- "both_by_LD"
    ## delete from the data set
    cat("\ndelete from the data\n")
    print(foo[ii[jj],])
    cat("\nout of \n")
    print(foo[ii,])
    foo <- foo[-ii[jj],]
  }
}

##-- apply exclusion --##
res.lookup <- foo

## add cis/trans information
res.lookup$cis_trans <- apply(res.lookup[, c("MarkerName", "pheno.soma", "pheno.olink")], 1, function(x){
  if(x[1] %in% subset(res.soma, pheno.soma == x[2] & cis_trans == "cis")$MarkerName | x[1] %in% subset(res.olink, pheno.olink == x[3] & cis_trans == "cis")$MarkerName){
    return("cis")
  }else{
    return("trans")
  }
})

#######################################################
####              collapse into regions            ####
#######################################################

## run through results and collapse signals if those are less than 1MB apart from each other
## for a given SOMAmer - Olink - genetic locus

## store names to ease regional mapping
n.sl <- c("pheno.soma", "Target", "rsid", "MarkerName", "chr", "pos", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P.value")
n.ol <- c("pheno.olink", "rsid", "MarkerName", "chr", "pos", "a_0", "a_1", "af", "beta", "se", "p")

## use observational pairings as start
res.regions <- lapply(1:nrow(res.obs), function(x){
  
  ## get the names for SOMAmer and Olink
  sl     <- res.obs$pheno.soma[x]
  ol     <- res.obs$pheno.olink[x]
  
  # print(x)
  
  ## get the subsets of relevant hits
  tmp.sl <- subset(res.soma, pheno.soma == sl)
  tmp.ol <- subset(res.olink, pheno.olink == ol)
  ## make unique
  tmp.ol <- unique(tmp.ol[, n.ol])
  
  ## do only if at least one hit in either
  if(nrow(tmp.sl) > 0 | nrow(tmp.ol) > 0){
    
    ## data frame to store the merged results
    res        <- array(data=NA, dim=c(nrow(tmp.sl)+nrow(tmp.ol),length(c(n.sl, n.ol))))
    res        <- as.data.frame(res)
    names(res) <- c(paste0(n.sl, ".soma"), paste0(n.ol, ".olink"))
    
    ## row counter
    kk         <- 0
    
    ## map SomaLogic data
    if(nrow(tmp.sl) > 0){
      for(j in 1:nrow(tmp.sl)){
        ## get the marker position
        ii.chr <- tmp.sl$chr[j]
        ii.pos <- tmp.sl$pos[j]
        ## search for Olink hits in close proximity (<= 1MB)
        ii.ol  <- which(tmp.ol$chr == ii.chr & tmp.ol$pos >= ii.pos - 1e6 & tmp.ol$pos <= ii.pos + 1e6)
        if(length(ii.ol) == 1){
          res[j,]               <- c(tmp.sl[j, n.sl], tmp.ol[ii.ol, n.ol])
          ## delete from Olink data set
          tmp.ol                <- tmp.ol[-ii.ol, n.ol]
          ## increase counter
          kk                    <- kk+1
        }else if(length(ii.ol) == 0){
          ## add the Olink name just to ensure unique mapping
          res[j,1:(length(n.sl)+1)] <- c(tmp.sl[j, n.sl], ol)  
          ## increase counter
          kk                         <- kk+1
        }else{
          print(tmp.ol[ii.ol, n.ol])
          cat("mapped multiple regions at", tmp.sl$MarkerName[j], "\n")
        }
      }
    }
    ## do Olink in case something is left
    if(nrow(tmp.ol) > 0){
      for(j in 1:nrow(tmp.ol)){
        res[kk+j, c(1,(1:length(n.ol)) + length(n.sl))] <- c(sl,tmp.ol[j, n.ol])
        ## increase counter
        kk <- kk+1
      }
    }
    return(res[1:kk,])
  }else{
    return(NA)
  }
})
res.regions <- do.call(rbind, res.regions)
res.regions <- unique(res.regions)
## delete complete NA entries
res.regions <- res.regions[apply(res.regions, 1, function(x) sum(is.na(x)) != ncol(res.regions)),]

## edit some names
names(res.regions) <- gsub(".soma.soma", ".soma", names(res.regions), fixed = T)
names(res.regions) <- gsub(".olink.olink", ".olink", names(res.regions), fixed = T)

## create column indicating whether two signals, SOMAlogic only, Olink only
res.regions$platform <- apply(res.regions[, c("MarkerName.soma", "MarkerName.olink")], 1, function(x){
  ifelse(!is.na(x[1]) & !is.na(x[2]), "Both", ifelse(!is.na(x[1])& is.na(x[2]), "SOMAscan", "Olink"))
})
table(res.regions$platform)     

#######################################################
####          fill gaps with look up data          ####
#######################################################

## do regional look-up for each of the pairs: e.g search for an signal in <= 500kb distance
tmp <- subset(res.regions, is.na(MarkerName.soma) | is.na(MarkerName.olink))

## loop to get the needed information
foo <- lapply(1:nrow(tmp), function(c){
  ## check what is present
  if(!is.na(tmp$MarkerName.soma[c])){
    return(data.frame(MarkerName=tmp$MarkerName.soma[c], pheno.soma = tmp$pheno.soma[c], pheno.olink = tmp$pheno.olink[c],
                      chr=tmp$chr.soma[c], pos=tmp$pos.soma[c], platform.lookup="Olink"))
  }else{
    return(data.frame(MarkerName=tmp$MarkerName.olink[c], pheno.soma = tmp$pheno.soma[c], pheno.olink = tmp$pheno.olink[c],
                      chr=tmp$chr.olink[c], pos=tmp$pos.olink[c], platform.lookup="SOMAscan"))
  }
})
foo <- do.call(rbind, foo)

## write to file to do look-up on HPC
write.table(foo, "Regional.look.up.SOMAscan.Olink.20200827.txt", sep="\t", row.names=F, quote=F, col.names = F)
## test for completeness
foo$id     <- apply(foo, 1, function(c){
  if(c[6] == "Olink"){
    return(paste(c[3], c[1], "regional", sep = "."))
  }else{
    return(paste(c[2], c[1], "regional", sep = "."))
  }
})


## look-up in SOMAscan
ii.s <- unique(subset(foo, platform.lookup == "SOMAscan")$id)
ii.o <- unique(subset(foo, platform.lookup == "Olink")$id)

## loop over SOMAscan files
ii.s <- lapply(ii.s, function(x){
  tmp <- read.table(paste0("../output/", x), sep="\t", header=T)
})
ii.s <- do.call(rbind, ii.s)

## loop over Olink files
ii.o <- lapply(ii.o, function(x){
  tmp <- read.table(paste0("../output/", x), sep="\t", header=T)
})
ii.o <- do.call(rbind, ii.o)
## add p-value
ii.o$p <- 10^-(ii.o$log10p)

#----------------------------------#
##--       fill the gaps        --##
#----------------------------------#

back.up <- res.regions

## fill the missing data on SOMAscan
r.names.sl   <- grep("soma", names(res.regions), value=T)
## omit one
r.names.sl   <- r.names.sl[-2]
## which names to replace with
g.names.sl   <- gsub("\\.soma", "", r.names.sl)
names(ii.s)  <- gsub("Pvalue", "P.value", names(ii.s))

## now replace matching entries
for(j in 1:nrow(ii.s)){
  ## find the gap
  ii                          <- which(res.regions$MarkerName.olink == ii.s$Markername.Olink[j] & res.regions$pheno.soma == ii.s$pheno[j])
  print(ii)
  ## now replace with missing entries
  res.regions[ii, r.names.sl] <- ii.s[j, g.names.sl]
}

## now replace the Olink gaps
r.names.ol <- grep("olink", names(res.regions), value=T)
g.names.ol <- c("pGWAS_olink", "rsid", "MarkerName", "chr", "pos", "a_0", "a_1", "af", "beta", "se", "p")

## now replace matching entries
for(j in 1:nrow(ii.o)){
  ## find the gap
  ii <- which(res.regions$MarkerName.soma == ii.o$markername.soma[j] & res.regions$pheno.olink == ii.o$pGWAS_olink[j])
  print(ii)
  ## now replace with missing entries
  res.regions[ii, r.names.ol] <- ii.o[j, g.names.ol]
}

## check if anything is missing
summary(res.regions)

## recode alleles
res.regions$a_0.olink[which(res.regions$a_0.olink == T)] <- "T"
res.regions$a_1.olink[which(res.regions$a_1.olink == T)] <- "T"

#--------------------------------------#
##-- check for LD between sentinels --##
#--------------------------------------#

res.regions$R2_sentinels <- apply(res.regions[, c("rsid.soma", "rsid.olink")], 1, function(x){
  ## test whether the same signal
  if(x[1] == x[2]){
    return(1)
  }else{
    ## identify proxies from both
    tmp <- subset(ld.matrix, (SNP_A == x[1] & SNP_B == x[2]) | (SNP_A == x[2] & SNP_B == x[1]))
    if(nrow(tmp) > 0){
      return(tmp$R2[1])
    }else{
      ## code lower values as NA
      return(NA)
    }
  }
})

#---------------------------------------#
##-- check for LD with cond. signals --##
#---------------------------------------#

## check for LD with conditional signals from SomaLogic
tmp <- lapply(1:nrow(res.regions), function(x){
  ## test whether included
  if(res.regions$rsid.olink[x] %in% subset(cond.soma, MarkerName_sentinel == res.regions$MarkerName.soma[x] & pheno == res.regions$pheno.soma[x])$rsid){
    ii <- which(cond.soma$rsid == res.regions$rsid.olink[x] & cond.soma$pheno == res.regions$pheno.soma[x])
    # print(ii)
    return(data.frame(rsid.cond=cond.soma$rsid[ii], pheno.olink=res.regions$pheno.olink[x], pheno.soma=res.regions$pheno.soma[x],
                      rsid.olink=res.regions$rsid.olink[x],
                      LD.cond=1))
  }
  ## test whether any proxy might be in LD with a conditional hit
  pr <- subset(ld.matrix, R2 >= .8 & (SNP_A == res.regions$rsid.olink[x] | SNP_B == res.regions$rsid.olink[x]))
  ll <- unique(c(pr$SNP_A, pr$SNP_B))
  # print(ll)
  if(length(ll > 0)){
    ## test whether any of the proxies is in the conditional list
    if(sum(ll %in% subset(cond.soma, MarkerName_sentinel == res.regions$MarkerName.soma[x] & pheno == res.regions$pheno.soma[x])$rsid) > 0){
      ii <- which(cond.soma$rsid %in% ll & cond.soma$pheno == res.regions$pheno.soma[x])
      print(ii)
      ## get the LD information between both variants
      ld <- subset(pr, SNP_B == cond.soma$rsid[ii])
      print(ld)
      return(data.frame(rsid.cond=cond.soma$rsid[ii], pheno.olink=res.regions$pheno.olink[x], pheno.soma=res.regions$pheno.soma[x],
                        rsid.olink=res.regions$rsid.olink[x],
                        LD.cond=ld$R2))
    }
  }
})
tmp <- do.call(rbind, tmp)

## add to the data
res.regions <- merge(res.regions, tmp, all.x=T)

#---------------------------------------#
##--  check directional concordance  --##
#---------------------------------------#

## create a list of all SNPs of interest and corresponding proteins for look-up (one comprehensive data set)
tmp1        <- res.regions[, c("pheno.soma", "pheno.olink", "MarkerName.soma")]
tmp2        <- res.regions[, c("pheno.soma", "pheno.olink", "MarkerName.olink")]
names(tmp1) <- names(tmp2) <- c("pheno.soma", "pheno.olink", "MarkerName")
tmp         <- unique(rbind(tmp1, tmp2))

## write to file to do look up
write.table(tmp, "SOMAscan.Olink.complete.lookup.20200828.txt", sep="\t", row.names=F, col.names=F, quote=F)

## --> import look-up and create column for directional concordance <-- ##

ii   <- dir("../output/")
ii   <- grep("lookup", ii, value=T)
## separate into SOMAscan and Olink
ii.s <- grep("res_invn_X[0-9]", ii, value=T) 
ii.o <- grep("res_invn_X[0-9]", ii, value=T, invert = T) 

## head for SOMAscan and Olink
s.names <- c("rsid", "MarkerName", "Allele1", "Allele2", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Effect", "StdErr", "P-value",
             "Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal", "TotalSampleSize", "chr", "pos")
## rename allele coding for Olink to match with coding from METAL output
o.names <- c("chr", "rsid", "pos", "Allele2", "Allele1", "Freq1", "info", "Effect", "StdErr", "tval", "log10p", "MarkerName")

#----------------------------#
##-- create combined file --##
#----------------------------#

## SOMAscan (some are not found)
ii.s <- lapply(ii.s, function(x){
  ## test for empty files
  if(file.info(paste0("../output/", x))$size != 0){
    ## load the file
    tmp            <- read.table(paste0("../output/", x), sep="\t")
    ## assign names
    names(tmp)     <- s.names
    ## add protein
    tmp$pheno.soma <- strsplit(x, "\\.")[[1]][1]
    return(tmp)
  }else{
    cat("no lines for", x, "\n")
  }
})
# no file for res_invn_X3284_75.chr16:17564311_A_C.lookup 
ii.s <- do.call(rbind, ii.s)
## edit some names
names(ii.s)[-c(1:2,18:20)] <- paste0(names(ii.s)[-c(1:2,18:20)],".soma")

## Olink
ii.o <- lapply(ii.o, function(x){
  ## load the file (some files do not contain any lines)
  tryCatch({tmp <- read.table(paste0("../output/",x), sep=" ")
  ## assign names
  names(tmp)      <- o.names
  ## add protein
  tmp$pheno.olink <- strsplit(x, "\\.")[[1]][1]
  return(tmp)
  }, error=function(e){
    cat("no lines in", x, "\n")
    return(NA)
  })
})
## no complains
ii.o <- do.call(rbind, ii.o)
ii.o <- na.omit(ii.o)

## edit some names
names(ii.o)[4:11] <- paste0(names(ii.o)[4:11], ".olink")

## create final look-up file
tmp <- merge(tmp, ii.s, all.x=T) 
tmp <- merge(tmp, ii.o, by=c("pheno.olink", "MarkerName"), suffixes = c(".soma", ".olink"), all.x=T)

## find what is missing an rerun
subset(tmp, is.na(log10p.olink) | is.na(Effect.soma))
## one not captured (misses SOMAscan equivalent)

## recode alleles
tmp$Allele1.olink[which(tmp$Allele1.olink == T)] <- "T"
tmp$Allele2.olink[which(tmp$Allele2.olink == T)] <- "T"

## align effect estimates from Olink
tmp$Effect.olink.aligned <- ifelse(toupper(tmp$Allele1.soma) == tmp$Allele1.olink, tmp$Effect.olink, -tmp$Effect.olink)

#---------------------------------------#
##-- assign directional concordance  --##
#---------------------------------------#

res.regions$directional.concordant <- apply(res.regions[, c("MarkerName.soma", "MarkerName.olink", "pheno.soma", "pheno.olink", "R2_sentinels", "rsid.cond")], 1, function(x){
  ## test whether SNPs in LD
  if(!is.na(as.numeric(x[5]))){
    ## get the SOMAmer entry
    ii <- subset(tmp, MarkerName == x[1] & pheno.soma == x[3] & pheno.olink == x[4])
    if(sign(ii$Effect.soma) == sign(ii$Effect.olink.aligned)){
      return("consistent")
    }else{
      return("inconsistent")
    }
  }else if(!is.na(x[6])){
    ## test for conditional hit
    ## get the Olink entry
    ii <- subset(tmp, MarkerName == x[2] & pheno.soma == x[3] & pheno.olink == x[4])
    if(sign(ii$Effect.soma) == sign(ii$Effect.olink.aligned)){
      return("consistent")
    }else{
      return("inconsistent")
    }
  }else{
    ## distinct SNPs
    return("N/A")
  }
})


#---------------------------------------#
##--          add cis/trans          --##
#---------------------------------------#

res.regions$cis_trans <- ifelse(res.regions$MarkerName.soma %in% subset(res.soma, cis_trans == "cis")$MarkerName | res.regions$MarkerName.olink %in% subset(res.olink, cis_trans == "cis")$MarkerName,
                                "cis", "trans")

#---------------------------------------#
##--        assign replication       --##
#---------------------------------------#

res.regions$reason <- apply(res.regions[, c("power", "platform", "R2_sentinels", "P.value.soma", "p.olink", "rsid.cond", "directional.concordant", "pev.soma")], 1, function(x){
  if(x[1] == "yes"){
    ## test whether on both platforms
    if(x[2] == "Both"){
      ## regional sentinels from both studies
      if(!is.na(x[3]) & x[3] > .6){
        if(x[7] == "consistent"){
          return("genome-wide sig. & LD sentinel SNP")
        }else{
          return("genome-wide sig. but directionally inconsistent")
        }
      }else if(!is.na(x[6])){
        if(x[7] == "consistent"){
          return("genome-wide sig. & LD conditional SNP")
        }else{
          return("genome-wide sig. & LD conditional SNP but directionally inconsistent")
        }
      }else{
        return("distinct genome-wide SNPs")
      }
    }else{
      if(!is.na(x[3]) & x[3] > .6){
        if(x[7] == "consistent"){
          return("moderate sig. level & LD sentinel SNP")
        }else{
          return("moderate sig. level & LD sentinel SNP but directionally inconsistent")
        }
      }else if(!is.na(x[6])){
        if(x[7] == "consistent"){
          return("moderate sig. level & LD conditional SNP")
        }else{
          return("moderate sig. level & LD conditional SNP but directional inconsistent")
        }
      }else{
        if(as.numeric(x[4]) > 1e-5 & as.numeric(x[5]) < 4.5e-11){
          return("unique to Olink")
        }else if(as.numeric(x[5]) > 1e-5 & as.numeric(x[8]) > .05){
          return("unique to SOMAscan")
        }else{
          return("uncertain")
        }
      }
    }
    
  }else{
    return("uncertain")
  }
})
table(res.regions$reason, useNA = "always")
head(subset(res.regions, is.na(reason)))

########################################################################################################################################
####                                                          END OF SCRIPT                                                         ####
########################################################################################################################################
