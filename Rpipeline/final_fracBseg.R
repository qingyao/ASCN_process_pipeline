########################
### BAF segmentation ###
########################

fracBseg <- function(seriesName,chipType,workingdir,sourcedir,memory){
  setOption(aromaSettings, "memory/ram", memory)
  suppressWarnings(suppressMessages(library(DNAcopy)))
  cids <- gsub(".CEL","",list.files(paste(getwd(),'rawData',seriesName,chipType,sep="/")))
  #cids <- cids[c(2,3)]
  for (cid in cids){
    fn <- paste0(workingdir,"/processed/",seriesName,"/",cid, '/fracBseg.tab')
    cat("ID","chrom","loc.start", "loc.end", "num.mark", "seg.mean","seg.sd","seg.median", "seg.mad\n",sep="\t",file=fn,append = F)
    for (chrname in 1:23){
      data<- read.table(sprintf("%s/processed/%s/%s/fracB,chr%d.tab",workingdir,seriesName,cid,chrname),header=T)
      #cat("Working on", cid, "\n", file=stderr())
      #for (chrname in c('X', 'Y', 22:1)) {
      #cat("\tProcessing", chrname, "\n", file=stderr())
      posn <-data$BASEPOS
      baf1 <-data$VALUE

      ### mBAF mirrored
      mbaf <- unlist(sapply(baf1, function(x) {
        if (is.na(x) == FALSE) {
          if (x < 0.5) x <- 1-x
          if (x >= 0.5) x<-x
        }
        else x<-x
      }))
      mbaf[abs(mbaf-0.5) >0.35] <- NA
      yy <- mbaf
      # segment to estimate both LOH and BAF
      cna1 <- CNA(as.matrix(yy),
                  chrom=rep(chrname, length(yy)),
                  maploc=posn,
                  data.type="logratio",
                  sampleid=cid)

      smoo1 <- tryCatch({
        smooth.CNA(cna1, smooth.region=10,smooth.SD.scale=4)
      }, error = function(err) {
        # error handler picks up where error was generated
        print(paste("MY_ERROR:  ",err))
        mbaf <- unlist(sapply(baf1, function(x) {
          if (is.na(x) == FALSE) {
            if (x < 0.5) x <- 1-x
            if (x >= 0.5) x<-x
          }
          else x<-x
        }))
        yy <- mbaf
        cna1 <- CNA(as.matrix(yy),
                    chrom=rep(chrname, length(yy)),
                    maploc=posn,
                    data.type="logratio",
                    sampleid=cid)
         smooth.CNA(cna1, smooth.region=10,smooth.SD.scale=4)
      },finally = {
        message(paste("Processed fracB segmentation for sample:",cid,"Chr:",chrname))
        })

      seg1 <- segment(smoo1, min.width=5, verbose=0)
      ss1 <- segments.summary(seg1)
      write.table(ss1, file=fn, sep="\t", quote=FALSE,
                  append = T, row.names=FALSE, col.names=F)

    }
  }

  gc()
  ### remove the gaps of chromosomes ####
  for (cid in cids){
    fn <- paste0(workingdir,"processed/",seriesName,"/",cid, '/fracBseg.tab')

    file <- read.table(fn,header = T)
    gapfile <- read.table(sprintf("%s/PlatformInfo/%s_GapPos.tab",workingdir,chipType),header = T)
    newfile <- data.frame()
    for (chr in 1:23) {
      subfile <- subset(file,file$chrom == chr)
      gapstart <- gapfile$Gap_start[chr]
      gapend <- gapfile$Gap_end[chr]
      for (row in 1:nrow(subfile)) {
        if ((subfile[row,]$loc.start) < gapstart & (subfile[row,]$loc.end) > gapend) {
          emptyrow <- subfile[row,]
          emptyrow$loc.start <- gapstart
          emptyrow$loc.end <- gapend
          emptyrow[,5:ncol(emptyrow)] <- rep(0,ncol(emptyrow)-4)
          newrow <- subfile[row,]
          newrow$loc.start <- gapend
          subfile[row,]$loc.end  <- gapstart
          n <- subfile[row,]$num.mark
          subfile[row,]$num.mark <- n/2
          newrow$num.mark <- n/2
          if (row < nrow(subfile)) {

            subfile <- rbind(subfile[1:row,],emptyrow, newrow,subfile[(row+1):nrow(subfile),])
            next
          }
          else {
            subfile <- rbind(subfile[1:row,],emptyrow, newrow)
            next
          }
        }
        else if((subfile[row,]$loc.start) < gapstart & subfile[row,]$loc.end <= gapend &subfile[row,]$loc.end > gapstart) {
          emptyrow <- subfile[row,]
          emptyrow$loc.start <- gapstart
          emptyrow$loc.end <- gapend
          emptyrow[,5:ncol(emptyrow)] <- rep(0,ncol(emptyrow)-4)
          subfile[row,]$loc.end <-gapstart
          subfile <- rbind(subfile[1:row,],emptyrow, subfile[(row+1):nrow(subfile),])
        }
        else if ((subfile[row,]$loc.start) >=gapstart & (subfile[row,]$loc.start) < gapend & subfile[row,]$loc.end > gapend) {
          emptyrow <- subfile[row,]
          emptyrow$loc.start <- gapstart
          emptyrow$loc.end <- gapend
          emptyrow[,5:ncol(emptyrow)] <- rep(0,ncol(emptyrow)-4)
          subfile[row,]$loc.start <-gapend
          subfile <- rbind(subfile[1:(row-1),],emptyrow, subfile[(row):nrow(subfile),])
        } else next
      }
      newfile <- rbind(newfile,subfile)
    }
    write.table(newfile, file=fn, sep="\t", quote=FALSE, row.names=FALSE, col.names=T)
    ####

    rownames(file) <- 1:nrow(file)
    rownames(newfile) <- 1:nrow(newfile)

  }


}
