suppressWarnings(suppressMessages(library(genlasso)))
options("scipen"=100, "digits"=4)

getLmd <- function(remotedir,seriesName,arrayName,workingdir,chipType,filename){
  name <- chooseFile(filename,1) [[1]]
  segfile <- read.table(file.path(remotedir,seriesName,arrayName,name),header = T,stringsAsFactors = F)
  noseg <- nrow(segfile)
  lmdtable <- data.frame('cn'=c(0,0.3,0.5,1),'fracb'=c(0,0.05,0.075,0.2))
  if (noseg < 200) {
    l <- lmdtable[filename][1,]
  } else if(noseg <500){
    l <- lmdtable[filename][2,]
  } else if(noseg <1000){
    l <- lmdtable[filename][3,]
  } else {
    l <- lmdtable[filename][4,]
  }
  return (l)
}

getGP <- function(remotedir,seriesName,arrayName,workingdir,chipType,filename) {"1e4"}

stepFilter <- function(remotedir,seriesName,arrayName,workingdir,chipType,gp,lmd,filename){
  log <- vector()
  name <- chooseFile(filename,1) [[1]]
  valuecol <- chooseFile(filename,1) [[2]]
  probecol <- chooseFile(filename,1) [[3]]
  segStep <- stepwiseDifference(as.numeric(gp),file.path(remotedir,seriesName,arrayName,name),filename)
  segStep <- segStep[,c(2:ncol(segStep))]

  x = segStep[,valuecol]
  idx <- calculateLasso(segStep,lmd,filename)[[2]]
  lmd <- calculateLasso(segStep,lmd,filename)[[3]]
  wm <- calculateWeightedMean(idx,segStep,filename)
  ##chipType CHRO chr_start chr_end probes
  if (filename == 'cn') {
    chrPos <- read.table(sprintf("%s/PlatformInfo/%s_chrPos.tab",workingdir,chipType),header = T)
  } else if(filename =='fracb'){
    chrPos <- read.table(sprintf("%s/PlatformInfo/%s_chrPos_SNP.tab",workingdir,chipType),header = T)
  }

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
            tmpseg[,valuecol] <- wm[i]
            tmpseg[,probecol] <- chrPos[chrPos[,2]==chrPassed[j],5]
            newseg <- rbind(newseg, tmpseg)
            countProbes <- countProbes + tmpseg[,probecol]
          }
        }

      tmpseg <- vector()
      tmpseg <- segStep[idx[i],]
      tmpseg[,3] <- segStep[1,3]
      tmpseg[,valuecol] <- wm[i]
      tmpseg[,probecol] <- sum(segStep[c(1:idx[i]),probecol])-countProbes
      newseg<- rbind(newseg, tmpseg)

    }
    ## new segment and previous segment on the same chromosome
    else if (segStep[idx[i],2] == segStep[idx[i-1],2]){
      tmpseg <- vector()
      tmpseg <- segStep[idx[i],]
      tmpseg[,3] <- segStep[idx[i-1]+1,3] #start of last row +1's start
      tmpseg[,valuecol] <- wm[i] # weighted mean
      tmpseg[,probecol] <- sum(segStep[c((idx[i-1]+1):idx[i]),probecol]) # sum of probe numbers since last row +1
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
        tmpseg[,valuecol] <- wm[i]
        tmpseg[,probecol] <- chrPos[chrPos[,2]==chrPassed[1],5] - usedProbes
        newseg <- rbind(newseg, tmpseg)
        countProbes <- countProbes + tmpseg[,probecol]
        ## the chromosome is not the immediate next one
        if (length(chrPassed) > 2){
          log <- c(log,sprintf('stepFilter involves chromosomes merging between chr%s and chr%s!',segStep[idx[i-1],2],segStep[idx[i],2]))
          for (j in 2:(length(chrPassed)-1)){
            tmpseg <- segStep[idx[i-1],]
            tmpseg[,2] <- chrPassed[j]
            tmpseg[,c(3,4)] <- chrPos[chrPos[,2]==chrPassed[j],c(3,4)]
            tmpseg[,valuecol] <- wm[i]
            tmpseg[,probecol] <- chrPos[chrPos[,2]==chrPassed[j],5]
            newseg <- rbind(newseg, tmpseg)
            countProbes <- countProbes + tmpseg[,6]
          }
        }
        tmpseg <- segStep[idx[i],]
        tmpseg[,3] <- chrPos[chrPos[,2]==segStep[idx[i],2],3]
        tmpseg[,valuecol] <- wm[i]
        tmpseg[,probecol] <- sum(segStep[(idx[i-1]+1):idx[i],probecol]) -countProbes
        newseg <- rbind(newseg, tmpseg)

      }

  }
  newname <-chooseFile(filename,0) [[1]]
  write.table(newseg, file.path(remotedir,seriesName,arrayName,newname), sep="\t", quote=FALSE,row.names=FALSE)
  logfile <- file.path(remotedir,seriesName,arrayName,sprintf('%sseg,log.txt',filename))
  cat(as.character(Sys.time()), 'stepFilter warning:', paste(log, collapse='; '), '\n',sep='\t', file=logfile,append = T)
  cat(as.character(Sys.time()), sprintf("filtered segments with gap size = %s, lambda = %s, previous %s segments, now %s segments, %s indices \n",gp,lmd,length(x),nrow(newseg),length(idx)),sep='\t',file=logfile,append = T)
}

calculateLasso <- function(segmentData, lmd, filename) {
  valuecol <- chooseFile(filename,1)[[2]]
  set.seed(1)
  x <- segmentData[,valuecol]
  out = fusedlasso1d(x)
  minl <- min(out$lambda)
  lmd <- max(minl,lmd)
  beta1 = coef(out, lambda=lmd)$beta

  ## calculate weighted mean in each segment
  beta1 <- round(beta1,4)
  idx = which(abs(diff(beta1)) > 1e-4)
  return(list(beta1,idx,lmd))
}

calculateWeightedMean <- function(idx, segmentData, filename) {
  valuecol <- chooseFile(filename,1)[[2]]
  lasti=1
  wm=rep(0,length(idx)+1)
  for(i in 1:(length(idx)+1)){
    wm[i] <- ifelse(i!=length(idx)+1,
                    weighted.mean(segmentData[,valuecol][lasti:idx[i]],w=segmentData[,4][lasti:idx[i]]-segmentData[,3][lasti:idx[i]]),
                    weighted.mean(segmentData[,valuecol][lasti:nrow(segmentData)],w=segmentData[,4][lasti:nrow(segmentData)]-segmentData[,3][lasti:nrow(segmentData)]))

    lasti <- idx[i]+1}
  wm <- round(wm, 4)
  return (wm)
}

stepwiseDifference <- function(gapSize,segmentFile,filename){
  probecol <- chooseFile(filename,1)[[3]]
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
      tmpseg[,probecol+1] <- tmpseg[,probecol+1]+rmProbes
      newseg<-rbind(newseg,tmpseg)
      rmProbes <- 0
    } else{
      rmProbes <- rmProbes + seg[row,probecol]
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
