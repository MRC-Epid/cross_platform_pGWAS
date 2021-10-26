########################################
## function to plot hyprcoloc results

plot.hyprcoloc <- function(sum.stat, res, lab, snp.dat, snp.lab){
  
  ## 'sum.stat' -- summary statistics for plotting
  ## 'res'      -- results from coloc
  ## 'lab'      -- label for phenotypes
  ## 'snp.dat'  -- individual levels SNP data
  ## 'snp.lab'  -- infomation on the SNP to ease merging 
  
  #-----------------------#
  ##--  perpare input  --##
  #-----------------------#
  
  ## get the variant to be annotated
  tr <- c("soma", "olink", "trait")
  
  ## create log10p
  for(j in tr){
    if(j == "trait"){
      sum.stat[, paste0("log10p.",j)] <- -pchisq((sum.stat$beta/sum.stat$se)^2, df=1, lower.tail=F, log.p=T)/log(10)
    }else{
      sum.stat[, paste0("log10p.",j)] <- -pchisq((sum.stat[, paste0("Effect.",j)]/sum.stat[, paste0("StdErr.",j)])^2, df=1, lower.tail=F, log.p=T)/log(10)
    }
  }
  
  ## add LD information
  if(res$PP.H4.abf.soma > .8){
    rr <- res$rsid.soma
  }else{
    rr <- res$rsid.olink
  }
  
  print(rr)
  
  ## compute correlation
  ld.mat    <- cor(pheno[, tmp.info$id[which(tmp.info$rsid== rr)]], pheno[, tmp.info$id])^2
  ld.mat    <- data.frame(id=colnames(ld.mat), R2=t(ld.mat))
  ## ease merging
  ld.mat    <- merge(ld.mat, tmp.info)
  ## add to the data
  sum.stat  <- merge(sum.stat, ld.mat, all.x=T)
  
  ## define a colour gradient
  tmp.col   <- colorRampPalette(c("grey80", "firebrick"))(100) 
  
  #-----------------------#
  ##--    plot stats   --##
  #-----------------------#
  
  ## now draw the plot in a loop
  for(j in 1:length(tr)){
    
    ## get maximum
    mp <- max(sum.stat[, paste("log10p", tr[j], sep=".")])
    ## define boundaries for plotting
    if(mp < 7){
      ylim <- c(0,7)
    }else{
      ylim <- c(0, mp*1.1)
    }
    
    ## empty plot
    plot(sum.stat[, "pos"], sum.stat[, paste("log10p", tr[j], sep=".")], type = "n", xlab="", ylab=expression(-log[10]("p-value")),
         xaxt="n", yaxt="n", ylim=ylim)
    axis(2, lwd=.5)
    
    ## add varaints not in LD
    tmp <- subset(sum.stat, R2 < .2)
    points(tmp[, "pos"], tmp[, paste("log10p", tr[j], sep=".")], pch=21, lwd=.1, cex=.3, 
           bg=tmp.col[ceiling(tmp$R2*100)],
           col="grey70")
    
    ## now add the important stuff
    tmp <- subset(sum.stat, R2 >= .2)
    points(tmp[, "pos"], tmp[, paste("log10p", tr[j], sep=".")], pch=21, lwd=.1, cex=.4, col= "grey20", 
           bg=tmp.col[ceiling(tmp$R2*100)], 
           type="p")
    
    ## add trait
    if(tr[j] != "trait"){
      legend("topleft", pch=NA, lty=0, cex=.5, 
             legend = paste0(lab$label[which(lab$name == tr[j])], ": H3=", sprintf("%.1f", res[, paste0("PP.H3.abf.", tr[j])]*100),"% | ",
                             "H4=", sprintf("%.1f", res[, paste0("PP.H4.abf.", tr[j])]*100), "%"), 
             bty="n")
    }else{
      legend("topleft", pch=NA, lty=0, cex=.5, legend = lab$label[which(lab$name == tr[j])], bty="n")
    }
    
    
    ## add legend for LD-colouring
    if(j == 1){
      
      ## get parameter to ease plotting
      pm <- par("usr")
      ## legend for LD
      l <- seq(pm[1]+(pm[2]-pm[1])*.7, pm[1]+(pm[2]-pm[1])*.9, length.out = 100)
      rect(l-(l[2]-l[1])/2, pm[3]+(pm[4]-pm[3])*.7, l+(l[2]-l[1])/2, pm[3]+(pm[4]-pm[3])*.8, border=NA, col=tmp.col)
      ## box
      rect(l[1]-(l[2]-l[1])/2, pm[3]+(pm[4]-pm[3])*.7, l[100]+(l[2]-l[1])/2, pm[3]+(pm[4]-pm[3])*.8, border="black", col=NA, lwd=.3)
      ## simple axis
      text(l[c(1,50,100)], pm[3]+(pm[4]-pm[3])*.7, labels=c(0,.5,1), pos=1, cex=.4, offset = .1)
      ## add header
      text(pm[1]+(pm[2]-pm[1])*.8, pm[3]+(pm[4]-pm[3])*.8, cex=.5, labels = paste("r2 with", rr), pos=3,
           offset = .15)
    }
  }
  
  #-----------------------#
  ##--    plot genes   --##
  #-----------------------#
  
  #------------------------------------#
  ##--        gene assignment       --##
  #------------------------------------#
  
  ## incorporate gene annotation
  library(biomaRt)
  
  ## get data on build 37
  gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37) 
  
  ## get all genes in the region
  tmp.genes <- getBM(attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype'),
                     filters = c('chromosome_name','start','end'),
                     values = list(sum.stat$chr[1], min(sum.stat$pos, na.rm=T), max(sum.stat$pos, na.rm=T)),
                     mart = gene.ensembl)
  
  ## restrict to protein encoding genes for now
  tmp.genes <- subset(tmp.genes, gene_biotype == "protein_coding")
  ## sort by start
  tmp.genes <- tmp.genes[order(tmp.genes$start_position),]
  
  ## dummay for the line in the plot
  tmp.genes$line <- NA
  
  ## start sorting
  l <- 0
  
  ## loop over everything, collect genes row-wise
  while(sum(is.na(tmp.genes$line)) > 0){
    ## increase line 
    l <- l+1
    e <- 0
    for(j in 1:nrow(tmp.genes)){
      ## test whether new start is larger than the current end
      if(tmp.genes$start_position[j] > e+4e4 & is.na(tmp.genes$line[j])){
        ## assign line to be drawn
        tmp.genes$line[j] <- l
        ## assign new end
        e                 <- tmp.genes$end_position[j]
      }
    }
  }
  
  print(head(tmp.genes))
  
  ## now plot it just below
  par(mar=c(1.5,1.5,.1,.5), tck=-.02, bty="o", lwd=.5)
  plot(range(sum.stat$pos), c(0,max(tmp.genes$line)+.5), yaxt="n", xaxt="n", ylab="", 
       xlab=paste("Genomic position on chromosome", sum.stat$chr[1]), type="n",
       ylim=rev(c(0,max(tmp.genes$line)+.5)))
  axis(1, lwd=.5)
  ## position of the lead variant
  ii <- which(sum.stat$rsid == rr)
  abline(v=sum.stat$position[ii], lwd=.3, col="red1")
  ## add genes
  arrows(tmp.genes$start_position, tmp.genes$line, tmp.genes$end_position, tmp.genes$line, lwd=.5, length = 0)
  ## add Gene Names on top
  text(tmp.genes$start_position+(tmp.genes$end_position - tmp.genes$start_position)/2, tmp.genes$line, cex=.3, font=3,
       xpd=NA, labels=tmp.genes$external_gene_name, pos=3, offset=.1)
}