suppressWarnings(suppressMessages(library(genlasso)))
options("scipen"=100, "digits"=4)

adjustMedian <- function(remotedir,seriesName,arrayName,workingdir,chipType){
  segfile <- read.table(file.path(remotedir,seriesName,arrayName,'provenance/segments,cn.tsv'),header = F,stringsAsFactors = F)
  colnames(segfile) <- c('sample_id','chromosome','start','end','value','probes')
  probefile <- read.table(file.path(remotedir,seriesName,arrayName,'probes,cn.tsv'),header = T,stringsAsFactors = F)
  med <- median(segfile[,5])
  segfile[,5] <-segfile[,5]-med
  probefile[,4] <- probefile[,4] -med
  write.table(segfile,file.path(remotedir,seriesName,arrayName,'segments,cn.tsv'), sep="\t", quote=FALSE,row.names=FALSE)
  write.table(probefile,file.path(remotedir,seriesName,arrayName,'probes,cn.tsv'), sep="\t", quote=FALSE,row.names=FALSE)
  logfile <- file.path(remotedir,seriesName,arrayName,'cnseg,log.txt')
  cat(as.character(Sys.time()), sprintf("adjusted probe and segment values with median: %s \n",med), file=logfile,append = F)
}

getLmd <- function(remotedir,seriesName,arrayName,workingdir,chipType){
  segfile <- read.table(file.path(remotedir,seriesName,arrayName,'provenance/segments,cn.tsv'),header = T,stringsAsFactors = F)
  noseg <- nrow(segfile)
  if (noseg < 200) {
    l <- 0
  } else if(noseg <500){
    l <- 0.3
  } else if(noseg <1000){
    l <- 0.5
  } else {
    l <- 1
  }
  return (l)
}

getGP <- function(remotedir,seriesName,arrayName,workingdir,chipType) {"1e4"}

stepFilter <- function(remotedir,seriesName,arrayName,workingdir,chipType,gp,lmd){
  log <- vector()
  segStep <- stepwiseDifference(as.numeric(gp),file.path(remotedir,seriesName,arrayName,'provenance/segments,cn.tsv'))
  # segStep <- stepwiseDifference(as.numeric(gp),file.path(remotedir,seriesName,arrayName,'segments,cn,5_sdundo_1.tsv'))
  segStep <- segStep[,c(2:ncol(segStep))]

  x = segStep[,5]
  idx <- calculateLasso(segStep,lmd)[[2]]
#  print(segStep[idx,])
  lmd <- calculateLasso(segStep,lmd)[[3]]
  wm <- calculateWeightedMean(idx,segStep)
  ##chipType CHRO chr_start chr_end probes
  chrPos <- read.table(sprintf("%s/PlatformInfo/%s_chrPos.tab",workingdir,chipType),header = T)
  # chrPos <- read.table(sprintf("%s/PlatformInfo/%s_chrPos.tab",'/Users/pgweb/arraydata/aroma/hg19/',chipType),header = T)
  newseg <- data.frame()
  idx <- c(idx,nrow(segStep)) ##each idx is the end of segment and add the last row of all segments
  for (i in 1:length(idx)){
    if (i==1){
      countProbes <- 0 ## count probes that are in the hidden segments
      if(segStep[idx[i],2]!=1){
        log <- c(log,'stepFilter involves chromosomes merging on chr1!')
        chrPassed <- 1 : segStep[idx[i],2]
        for (j in 1:(length(chrPassed)-1)){
            tmpseg <- segStep[idx[1],]
            tmpseg[,2] <- chrPassed[j]
            tmpseg[,c(3,4)] <- chrPos[chrPos[,2]==chrPassed[j],c(3,4)]
            tmpseg[,5] <- wm[i]
            tmpseg[,6] <- chrPos[chrPos[,2]==chrPassed[j],5]
            newseg <- rbind(newseg, tmpseg)
            countProbes <- countProbes + tmpseg[,6]
          }
        }

      tmpseg <- vector()
      tmpseg <- segStep[idx[i],]
      tmpseg[,3] <- segStep[1,3]
      tmpseg[,5] <- wm[i]
      tmpseg[,6] <- sum(segStep[c(1:idx[i]),6])-countProbes
      newseg<- rbind(newseg, tmpseg)

    }
    ## new segment and previous segment on the same chromosome
    else if (segStep[idx[i],2] == segStep[idx[i-1],2]){
      tmpseg <- vector()
      tmpseg <- segStep[idx[i],]
      tmpseg[,3] <- segStep[idx[i-1]+1,3] #start of last row +1's start
      tmpseg[,5] <- wm[i] # weighted mean
      tmpseg[,6] <- sum(segStep[c((idx[i-1]+1):idx[i]),6]) # sum of probe numbers since last row +1
      newseg<- rbind(newseg, tmpseg)
    }
    ## new segment starts with a new chromosome
    else if (segStep[idx[i],2] != segStep[idx[i-1],2]){
      ## new segment starts with another chromosome instead of the same one
        countProbes <- 0 ## count probes that are in the hidden segments
        chrPassed <- segStep[idx[i-1],2] : segStep[idx[i],2]
        usedProbes <- sum(newseg[newseg[,2]==chrPassed[1],6]) ## probes used in the segments in that chromosome
        tmpseg <- segStep[idx[i-1],]
        tmpseg[,3] <- tmpseg[,4]
        tmpseg[,4] <- chrPos[chrPos[,2]==chrPassed[1],4]
        tmpseg[,5] <- wm[i]
        tmpseg[,6] <- chrPos[chrPos[,2]==chrPassed[1],5] - usedProbes
        newseg <- rbind(newseg, tmpseg)
        countProbes <- countProbes + tmpseg[,6]
        ## the chromosome is not the immediate next one
        if (length(chrPassed) > 2){
          log <- c(log,sprintf('stepFilter involves chromosomes merging between chr%s and chr%s!',segStep[idx[i-1],2],segStep[idx[i],2]))
          for (j in 2:(length(chrPassed)-1)){
            tmpseg <- segStep[idx[i-1],]
            tmpseg[,2] <- chrPassed[j]
            tmpseg[,c(3,4)] <- chrPos[chrPos[,2]==chrPassed[j],c(3,4)]
            tmpseg[,5] <- wm[i]
            tmpseg[,6] <- chrPos[chrPos[,2]==chrPassed[j],5]
            newseg <- rbind(newseg, tmpseg)
            countProbes <- countProbes + tmpseg[,6]
          }
        }
        tmpseg <- segStep[idx[i],]
        tmpseg[,3] <- chrPos[chrPos[,2]==segStep[idx[i],2],3]
        tmpseg[,5] <- wm[i]
        tmpseg[,6] <- sum(segStep[(idx[i-1]+1):idx[i],6]) -countProbes
        newseg <- rbind(newseg, tmpseg)

      }

  }
  write.table(newseg, file.path(remotedir,seriesName,arrayName,'segments,cn.tsv'), sep="\t", quote=FALSE,row.names=FALSE)
  logfile <- file.path(remotedir,seriesName,arrayName,'cnseg,log.txt')
  cat(as.character(Sys.time()), 'stepFilter warning:', paste(log, collapse='; '), '\n',sep='\t', file=logfile,append = T)
  cat(as.character(Sys.time()), sprintf("filtered segments with gap size = %s, lambda = %s, previous %s segments, now %s segments, %s indices \n",gp,lmd,length(x),nrow(newseg),length(idx)),sep='\t',file=logfile,append = T)
}

calculateLasso <- function(segmentData, lmd) {
  set.seed(1)
  x <- segmentData$value
  out = fusedlasso1d(x)
  minl <- min(out$lambda)
  lmd <- max(minl,lmd)
  beta1 = coef(out, lambda=lmd)$beta

  ## calculate weighted mean in each segment
  beta1 <- round(beta1,4)
  idx = which(abs(diff(beta1)) > 1e-4)
  return(list(beta1,idx,lmd))
}

calculateWeightedMean <- function(idx, segmentData) {
  lasti=1
  wm=rep(0,length(idx)+1)
  for(i in 1:(length(idx)+1)){
    wm[i] <- ifelse(i!=length(idx)+1,
                    weighted.mean(segmentData$value[lasti:idx[i]],w=segmentData$end[lasti:idx[i]]-segmentData$start[lasti:idx[i]]),
                    weighted.mean(segmentData$value[lasti:nrow(segmentData)],w=segmentData$end[lasti:nrow(segmentData)]-segmentData$start[lasti:nrow(segmentData)]))

    lasti <- idx[i]+1}
  wm <- round(wm, 4)
  return (wm)
}

stepwiseDifference <- function(gapSize,segmentFile){
  seg <- read.table(segmentFile,stringsAsFactors = F,header = T)
  print(nrow(seg))
  accumulatePos <- 0
  newseg <-data.frame()
  rmProbes <- 0
  for (row in 1:nrow(seg)){ ###work on this to add the probe number from the removed small segments
    segLen <- seg[row,4]-seg[row,3]
    if ( segLen> gapSize) {
      accumulatePos <- accumulatePos + segLen
      tmpseg<-cbind(accumulatePos,seg[row,c(1:6)])
      tmpseg[,7] <- tmpseg[,7]+rmProbes
      newseg<-rbind(newseg,tmpseg)
      rmProbes <- 0
    } else{
      rmProbes <- rmProbes + seg[row,6]
    }
  }
  return(newseg)
}

statsSmallSeg <- function(gapSizes,segmentFile){
  seg <- read.table(segmentFile,stringsAsFactors = F,header = T)
  smallSeg <- lapply(gapSizes, function(gapSize){
    newseg <-data.frame()
    for (row in 1:nrow(seg)){
      segLen <- seg[row,4]-seg[row,3]
      if ( segLen< gapSize) {
        newseg<-rbind(newseg,data.frame(length=seg[row,4]-seg[row,3],probes=seg[row,5]))
      }
    }
    return (newseg)
  })
  names(smallSeg) <- gapSizes
  return(smallSeg)
}

rmGaps <- function(remotedir,seriesName,arrayName,workingdir,chipType){
  fn <- file.path(remotedir,seriesName,arrayName,'segments,cn.tsv')
  file <- read.table(fn,header = T, stringsAsFactors = F)
  gapfile <- read.table(sprintf("%s/PlatformInfo/%s_GapPos.tab",workingdir,chipType),header = T)
  newfile <- data.frame()
  for (chr in 1:23) {
    subfile <- subset(file,file$chromosome == chr)
    gapstart <- gapfile$Gap_start[chr]
    gapend <- gapfile$Gap_end[chr]
    for (row in 1:nrow(subfile)) {
      if ((subfile[row,]$start) < gapstart & (subfile[row,]$end) > gapend) {
        newrow <- subfile[row,]
        newrow$start <- gapend
        subfile[row,]$end  <- gapstart
        n <- subfile[row,]$probes
        RatioBefAft <- (gapstart - subfile[row,]$start) / (newrow$end - gapend)
        subfile[row,]$probes <- round(RatioBefAft/(RatioBefAft+1) * n)
        newrow$probes <- round(1/(RatioBefAft+1) * n)
        if (row < nrow(subfile)) {
            subfile <- rbind(subfile[1:row,], newrow, subfile[(row+1):nrow(subfile),])
            next
          }
        else {
          subfile <- rbind(subfile[1:row,], newrow)
          next
        }
      }
      else if((subfile[row,]$start) < gapstart & subfile[row,]$end <= gapend &subfile[row,]$end > gapstart) {
        subfile[row,]$end <-gapstart
      }
      else if ((subfile[row,]$start) >=gapstart & (subfile[row,]$start) < gapend & subfile[row,]$end > gapend) {
        subfile[row,]$start <-gapend
      } else next
    }
    newfile <- rbind(newfile,subfile)
  }
  newfile <- newfile[newfile$probes!=0,]
  write.table(newfile, file=fn, sep="\t", quote=FALSE, row.names=FALSE, col.names=T)
  logfile <- file.path(remotedir,seriesName,arrayName,'cnseg,log.txt')
  cat(as.character(Sys.time()), sprintf("removed centromere gaps, now %s segments", nrow(newfile)), file=logfile,append = T)
}
