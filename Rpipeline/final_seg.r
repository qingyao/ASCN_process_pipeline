suppressWarnings(suppressMessages(library(DNAcopy)))

########################
### CN segmentation ####
########################


cnsegPerArray <- function(workingdir,seriesName, arrayName, undosd, chipType){
  dir.create(file.path(workingdir,'processed',seriesName,arrayName), showWarnings = FALSE)
  fn <- file.path(workingdir,'processed',seriesName,arrayName,'segments,cn.tsv')
  cat("sample_id","chromosome","start", "end", "value", "probes\n",sep="\t",file=fn,append = F)
  fp <- file.path(workingdir,'processed',seriesName,arrayName,'segments,cn,provenance.tsv')
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
  cat(as.character(Sys.time()), sprintf("undosd = %s",undosd), file=logfile, append=T)
}

############################
#### fracb segmentation ####
############################

fracbsegPerArray <- function(seriesName,arrayName, chipType,workingdir){

    fn <- file.path(workingdir,"processed",seriesName,arrayName,'fracbseg.tsv')
    cat("ID","chrom","loc.start", "loc.end", "num.mark", "seg.mean","seg.sd","seg.median", "seg.mad\n",sep="\t",file=fn,append = F)
    fp <- file.path(workingdir,"processed",seriesName,arrayName,'fracbseg,provenance.tsv')
    cat("ID","chrom","loc.start", "loc.end", "num.mark", "seg.mean","seg.sd","seg.median", "seg.mad\n",sep="\t",file=fp,append = F)
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
