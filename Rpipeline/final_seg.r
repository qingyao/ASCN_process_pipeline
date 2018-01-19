suppressWarnings(suppressMessages(library(DNAcopy)))

########################
### CN segmentation ####
########################


cnsegPerArray <- function(workingdir,seriesName, arrayName, undosd, chipType){
  dir.create(file.path(workingdir,'processed',seriesName,arrayName), showWarnings = FALSE)
  dir.create(file.path(workingdir,'processed',seriesName,arrayName,'provenance'), showWarnings = FALSE)
  fn <- file.path(workingdir,'processed',seriesName,arrayName,'segments,cn.tsv')
  cat("sample_id","chromosome","start", "end", "value", "probes\n",sep="\t",file=fn,append = F)
  fp <- file.path(workingdir,'processed',seriesName,arrayName,'provenance','segments,cn.tsv')
  cat("sample_id","chromosome","start", "end", "value", "probes\n",sep="\t",file=fp,append = F)
  for (chrname in 1:23){
    data<- read.table(sprintf('/Volumes/arraymapMirror/arraymap/hg19/%s/%s/probes,cn,chr%d.tsv',seriesName,arrayName,chrname),header=T)
    cna1 <- CNA(genomdat=data$VALUE,
                chrom=rep(chrname, length(data$VALUE)),
                maploc=data$BASEPOS,
                data.type="logratio",
                sampleid=arrayName)
    smoo1 <- smooth.CNA(cna1, smooth.region=5,smooth.SD.scale=2)
    message(paste("Processed CN segmentation for sample:",arrayName,"Chr:",chrname))

    seg1 <- segment(smoo1, min.width=5, verbose=0, undo.splits='sdundo',undo.SD=undosd)
    ss1 <- segments.summary(seg1)[c(1:4,6,5)]
    write.table(ss1, file=fn, sep="\t", quote=FALSE,
                append = T, row.names=FALSE, col.names=F)
    write.table(ss1, file=fp, sep="\t", quote=FALSE,
                append = T, row.names=FALSE, col.names=F)

    }
  logfile <- file.path(workingdir,'processed',seriesName,arrayName,'cnseg,log.txt')
  cat(Sys.time(), sprintf("undosd = %s",undosd)) ## add more information later ??not written to log file, and no new line
  rmGaps(workingdir,seriesName,arrayName,chipType)
}

############################
#### fracb segmentation ####
############################

fracBsegperArray <- function(seriesName,arrayName, chipType,workingdir){

    fn <- file.path(workingdir,"processed",seriesName,arrayName,'fracBseg.tab')
    cat("ID","chrom","loc.start", "loc.end", "num.mark", "seg.mean","seg.sd","seg.median", "seg.mad\n",sep="\t",file=fn,append = F)
    for (chrname in 1:23){
      data<- read.table(sprintf("%s/processed/%s/%s/fracB,chr%d.tab",workingdir,seriesName,arrayName,chrname),header=T)
      #cat("Working on", arrayName, "\n", file=stderr())
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
                  sampleid=arrayName)


      smoo1 <- smooth.CNA(cna1, smooth.region=10,smooth.SD.scale=4)

      message(paste("Processed fracB segmentation for sample:",arrayName,"Chr:",chrname))

      seg1 <- segment(smoo1, min.width=5, verbose=0)
      ss1 <- segments.summary(seg1)
      write.table(ss1, file=fn, sep="\t", quote=FALSE,
                  append = T, row.names=FALSE, col.names=F)

    }
  }

rmGaps <- function(workingdir,seriesName,arrayName,chipType){
  fn <- file.path(workingdir,'processed',seriesName,arrayName,'segments,cn.tsv')
  file <- read.table(fn,header = T)
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
        RatioBefAft <- (gapstart - subfile[row,]$start) / (subfile[row,]$end - gapend)
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
  write.table(newfile, file=fn, sep="\t", quote=FALSE, row.names=FALSE, col.names=T)

}
