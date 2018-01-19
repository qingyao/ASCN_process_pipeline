########################
### CN segmentation ####
########################


cnsegPerArray <- function(workingdir,seriesName, cid, undosd, chipType){
  dir.create(file.path(workingdir,'processed',seriesName,cid), showWarnings = FALSE)
  fn <- file.path(workingdir,'processed',seriesName,cid, 'segments,cn.tsv')
  cat("sample_id","chromosome","start", "end", "value", "probes\n",sep="\t",file=fn,append = F)

  for (chrname in 1:23){
    data<- read.table(sprintf('/Volumes/arraymapMirror/arraymap/hg19/%s/%s/probes,cn,chr%d.tsv',seriesName,cid,chrname),header=T)
    cna1 <- CNA(genomdat=data$VALUE,
                chrom=rep(chrname, length(data$VALUE)),
                maploc=data$BASEPOS,
                data.type="logratio",
                sampleid=cid)
    smoo1 <- smooth.CNA(cna1, smooth.region=5,smooth.SD.scale=2)
    message(paste("Processed CN segmentation for sample:",cid,"Chr:",chrname))

    seg1 <- segment(smoo1, min.width=5, verbose=0, undo.splits='sd.undo',undo.SD=1)
    ss1 <- segments.summary(seg1)[c(1:4,6,5)]
    write.table(ss1, file=fn, sep="\t", quote=FALSE,
                append = T, row.names=FALSE, col.names=F)
    }
  rmGaps(workingdir,seriesName,cid,chipType)
}

cnseg <- function(seriesName,arrayName=NULL,chipType,workingdir,sourcedir,memory,sdforce){
  if (sdforce == 0) {
    if (chipType %in% c('Mapping10K_Xba142')){
      undosd = 2
    } else if (chipType %in% c('Mapping50K_Hind240','Mapping50K_Xba240','Mapping250K_Nsp','Mapping250K_Sty','CytoScan750K_Array','GenomeWideSNP_5')){
      undosd = 3
    } else if (chipType %in% c('GenomeWideSNP_6','CytoScanHD_Array')){
      undosd = 5
    }
  } else {
    undosd = sdforce
  }
  print(paste(seriesName,arrayName,chipType,workingdir,sourcedir,memory,sdforce,undosd))
  setOption(aromaSettings, "memory/ram", memory)
  suppressWarnings(suppressMessages(library(DNAcopy)))
  cids <- gsub(".CEL","",list.files(paste('/Volumes/arraymapIncoming/aroma/aromaRaw',seriesName,chipType,sep="/")))
  dir.create(file.path(workingdir,'processed',seriesName), showWarnings = FALSE)
  if (is.null(arrayName)) {
    for (cid in cids){
      cnsegPerArray(workingdir,seriesName, cid, undosd)
    }
  } else{
    cnsegPerArray(workingdir,seriesName, arrayName, undosd)
  }

  gc()


  ### remove the gaps of chromosomes ####
  for (cid in cids){
    rmGaps(workingdir,seriesName,cid,chipType)
    ####
  }

}

fracBseg <- function(seriesName,chipType,workingdir,sourcedir,memory){
  setOption(aromaSettings, "memory/ram", memory)
  suppressWarnings(suppressMessages(library(DNAcopy)))

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
    rmGaps(workingdir,seriesName,cid,chipType)
  }

rmGaps <- function(workingdir,seriesName,cid,chipType){
  fn <- file.path(workingdir,'processed',seriesName,cid, 'segments,cn.tsv')

  file <- read.table(fn,header = T)
  gapfile <- read.table(sprintf("%s/PlatformInfo/%s_GapPos.tab",workingdir,chipType),header = T)
  newfile <- data.frame()
  for (chr in 1:23) {
    subfile <- subset(file,file$chromosome == chr)
    gapstart <- gapfile$Gap_start[chr]
    gapend <- gapfile$Gap_end[chr]
    for (row in 1:nrow(subfile)) {
      if ((subfile[row,]$start) < gapstart & (subfile[row,]$end) > gapend) {
        # emptyrow <- subfile[row,]
        # emptyrow$start <- gapstart
        # emptyrow$end <- gapend
        # emptyrow[,5:ncol(emptyrow)] <- rep(0,ncol(emptyrow)-4)
        newrow <- subfile[row,]
        newrow$start <- gapend
        subfile[row,]$end  <- gapstart
        n <- subfile[row,]$probes
        subfile[row,]$probes <- n/2
        newrow$probes <- n/2
        if (row < nrow(subfile)) {

          subfile <- rbind(subfile[1:row,],newrow,subfile[(row+1):nrow(subfile),])#emptyrow,
          next
        }
        else {
          subfile <- rbind(subfile[1:row,], newrow)#emptyrow,
          next
        }
      }
      else if((subfile[row,]$start) < gapstart & subfile[row,]$end <= gapend &subfile[row,]$end > gapstart) {
        # emptyrow <- subfile[row,]
        # emptyrow$start <- gapstart
        # emptyrow$end <- gapend
        # emptyrow[,5:ncol(emptyrow)] <- rep(0,ncol(emptyrow)-4)
        subfile[row,]$end <-gapstart
        subfile <- rbind(subfile[1:row,], subfile[(row+1):nrow(subfile),])#emptyrow,
      }
      else if ((subfile[row,]$start) >=gapstart & (subfile[row,]$start) < gapend & subfile[row,]$end > gapend) {
        # emptyrow <- subfile[row,]
        # emptyrow$start <- gapstart
        # emptyrow$end <- gapend
        # emptyrow[,5:ncol(emptyrow)] <- rep(0,ncol(emptyrow)-4)
        subfile[row,]$start <-gapend
        subfile <- rbind(subfile[1:(row-1),], subfile[(row):nrow(subfile),])#emptyrow,
      } else next
    }
    newfile <- rbind(newfile,subfile)
  }
  write.table(newfile, file=fn, sep="\t", quote=FALSE, row.names=FALSE, col.names=T)

}
