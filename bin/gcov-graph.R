#!/usr/bin/env Rscript
#
args <- commandArgs(TRUE)

if (length(args) < 1 | "--help" %in% args) {
  cat("Usage: graph-gcov.R <gcov.gz file> [ gcov.gz file ... ]\n", file=stderr())
  q()
}

window <- 1000
dither <- 1000
for(i in 1:length(args)) {
  if (file.exists(args[i])) {
    if (grepl("gcov.gz$", args[i])) {
      cov <- read.table(gzfile(args[i]))
      output <- sub("gz$", "png", args[i])
      run <- sub(".gcov.gz$", "", args[i])
    } else {
      cov <- read.table(args[i])
      output <- paste(args[i], ".png", sep="")
    }
    png(output, width=1024, height=512)
    ref <- unique(cov$V1)
    refsize <- rep(0, length(ref))
    ref <- ref[order(refsize, decreasing=T)]
    refsize <- refsize[order(refsize, decreasing=T)]
    cov$filter = 0;
    for(j in 1:length(ref)) {
      refsize[j] <- max(cov$V2[cov$V1==ref[j]])
      window <- min(1000, refsize[j]-1)
      cov$filter[cov$V1==ref[j]] <- filter(cov$V3[cov$V1==ref[j]], rep(1/(window+1), window+1), sides=2)
    }
    size <- round(log10(refsize))
    size[1] <- size[1] + 1.5
    ifelse(length(ref) > 1, linepos <- rep(1, length(ref)), linepos <- 0.5)
    layout(matrix(1:length(ref), 1, length(ref)), size)
    par(oma=c(3,4,6,3))
    for(j in 1:length(ref)) {
      ifelse(j==1, par(mar=c(0,2,0,0)), ifelse(j==length(ref), par(mar=c(0,0,0,2)), par(mar=c(0,0,0,0))))
      tx <- cov$V2[cov$V1==ref[j]]
      ty <- cov$filter[cov$V1==ref[j]]
      toplot <- seq(1, length(tx), dither)
      plot(tx[toplot], log10(ty[toplot]+1), type='l', lwd=1, ylim=log10(range(cov$V3)+1), xaxt="n", yaxt="n", axes=F)
       if (j %% 2 == 1) {
         mtext(paste(sep="", "Reference: \n", ref[j], "(", max(tx), " nt)"), side=3, cex=1.25, line=0.5)
         axis(1, xlim=range(tx), cex.axis=1)
       } else {
         mtext(paste(sep="", "Reference: \n", ref[j], "(", max(tx), " nt)"), side=1, cex=1.25, line=1)
         axis(3, xlim=range(tx), cex.axis=1)
       }
       ifelse(j==1, axis(2, ylim=log10(range(cov$V3)+1), cex.axis=1), ifelse(j==length(ref), axis(4, ylim=log10(range(cov$V3)+1), cex.axis=1), NA))
    }
    ydetail <- paste(sep="", "\n(sliding window of size ", window, " bp, step size ", dither, " bp)")
    mtext(paste(sep="", "Log10 coverage", ydetail), side=2, line=1, cex=1, outer=TRUE)
    if(length(ref) > 1) mtext("Log10 coverage", side=4, line=1.5, cex=1, outer=TRUE)
    mtext(paste(sep="", "Coverage plot (", run, ")"), side=3, line=3.5, outer=T, cex=1.75)
    dev.off()
  } else {
    cat("Can't find args[i], skipping...", file=stderr())
  }
}
