##########################
### baf analysis final ###
##########################

suppressWarnings(suppressMessages(library(plyr)))

BafAnalysis <- function(seriesName,chipType,workingdir,sourcedir,memory) {
  options("scipen"=100, "digits"=4)

  # cids <- gsub(".CEL","",list.files(paste(getwd(),'rawData',seriesName,chipType,sep="/")))
  remoteworkingdir = "/Volumes/arraymapMirror/arraymap/hg19/"
  remoterawdir = "/Volumes/arraymapIncoming/aroma/aromaRaw/"
  cids <- gsub(".CEL","",list.files(paste(remoterawdir,seriesName,chipType,sep="/")))
  # cids <- cids[1]
  for (cid in cids){
    cat("Processing sample:",cid,'\n')
    filelist <- list.files(paste(getwd(),"processed",seriesName,cid,sep="/"))
    cnfile <- filelist[grep("fracB,chr",filelist)]
    allfracb <- data.frame()
    for (j in cnfile) {
      allfracb <- rbind(allfracb,read.table(paste(workingdir,"processed",seriesName,cid,j,sep="/"),header=TRUE))
    }

    ## read the segments,fracB file; pre-filtering the segments.
    allseg <- read.table(sprintf("%s/%s/%s/fracBseg.tab",remoteworkingdir,seriesName,cid),header=T)
    #dynamic filtering, compare each time the new segment mean with the next one
    ss4 <- data.frame()
    for (chr in 1:23) {
      rows <- nrow(ss4)
      ss2 <- subset(allseg,allseg$chrom==chr)
      #ss2 <- subset(ss2,ss2$num.mark>20)
      count=0
      ss4 <- rbind(ss4,ss2[1,])
      if (nrow(ss2) >1) {
        for (i in 1:(nrow(ss2)-1)){
          lines <- abs(ss2$seg.mean[i+1] - ss4$seg.mean[rows+i-count]) > 0.035
          if (lines == 1) ss4 <- rbind(ss4,ss2[i+1,])
          if (lines == 0) {
            ss4$loc.end [rows+i-count] <- ss2$loc.end[i+1]
            ss4$seg.mean[rows+i-count] <- ss2$seg.mean[i+1]
            ss4$num.mark[rows+i-count] <- ss4$num.mark[rows+i-count]+ss2$num.mark[i+1]
            #ss4[rows+i-1-count,4:ncol(ss2)] <- ss2[i,4:ncol(ss2)]
            count= count+1
          }
        }
      }

    }
    allseg <- ss4


    Out <- sprintf("%s/processed/%s/%s/segments,fracb.tsv",workingdir,seriesName,cid)
    cat("ID","chr","loc.start","loc.end","fracB\n",sep="\t",file=Out,append=FALSE)
    for (chr in 1:23){
      #cat(sprintf("Chr%d \n",chr),file=Out,append=TRUE)
      seg <- subset(allseg, allseg$chrom ==chr & allseg$loc.end - allseg$loc.start > 1000000 & allseg$num.mark >0)
      ##remove extreme values
      fracb <- subset(allfracb, allfracb$CHRO ==chr)

      for (j in 1:nrow(seg)){
        range <- c(seg[j,"loc.start"], seg[j,"loc.end"])
        id <- which(fracb$BASEPOS>=range[1] & fracb$BASEPOS< range[2])
        subfracb <- fracb[id,]
        if (nrow(subfracb) <= 5) next
        nbin <- 50
        bin <- seq(0,1,length.out = nbin)
        k <- vector()
        for (i in 1:(nbin-1)) {
          k[i] <- sum(subfracb$VALUE < bin[i+1] & subfracb$VALUE > bin[i],na.rm = T)
        }

        period <- 10
        kn <- vector()
        for(i in 1:(nbin-period)){
          kn[i] <- mean(k[i:(i+period-1)])
        }
        kn <- log2(kn+15)

        ext <- nbin/10
        l <- length(kn)
        kn2 <- kn[-c(1:ceiling(ext/2),(l-ceiling(ext/2)+1):l)]

        m<- mean(kn2)
        Vec <- which(kn2>m)
        if (length(Vec) == 0 & mean(kn) > m) peak <- c(0.025,0.975)
        else{
            Breaks <- sort(c(1, which(diff(Vec) != 1),which(diff(Vec) != 1)+1, length(Vec)))
            Ranges <- matrix(Vec[Breaks],ncol=2,byrow=T)

            # find midpoint
            peak <- apply(Ranges,1,mean)
            peak <- peak/(nbin/10*9-period)
            peak <- sort(c(peak,1-peak))
            peak <- rm.near.dy(peak)
        }
        if (length(peak) == 0) next
        for (line in 1:length(peak)){
          cat(cid,chr,range,peak[line],'\n',sep="\t",file=Out,append=TRUE)
        }
      }
    }
    # dev.off()
  }
}


  ##dynamic remove
rm.near.dy <- function(x,distance=0.15) {
    count = 0
    i = 1
    #i <= (length(x)+count-1)
    while (i <= (length(x)+count-1)){
      c <- 0
      while (x[i+1-count]-x[i-c-count] <= distance) {
        #x[i+1-count]-x[i-c-count] <= distance
        i <- i+1
        c <- c+1
        if (i+1-count > length(x)) break
        #i+1-count > length(x)
      }
      x[i-c-count]<- mean(x[(i-c-count):(i-count)])
      x
      if (c>0) x<- x[-c((i-c-count+1):(i-count))]
      count <- count+c
      i <- i+1
    }
    x
}