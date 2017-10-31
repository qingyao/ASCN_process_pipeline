########################
### CN segmentation ####
########################

cnseg <- function(seriesName,chipType,workingdir,memory,sourcedir,sdforce){
  if (sdforce == 0) {
    if (chipType %in% c('Mapping10k_Xba142')){
      undosd = 2
    } else if (chipType %in% c('Mapping50K_Hind240','Mapping50K_Xba240','Mapping250K_Nsp','Mapping250K_Sty','CytoScan750K','GenomeWideSNP_5')){
      undosd = 3
    } else if (chipType %in% c('GenomeWideSNP_6','CytoScanHD')){
      undosd = 5
    }
  } else undosd = sdforce
  setOption(aromaSettings, "memory/ram", memory)
  suppressWarnings(suppressMessages(library(DNAcopy)))
  cids <- gsub(".CEL","",list.files(paste('/Volumes/arraymapIncoming/aroma/aromaRaw',seriesName,chipType,sep="/")))
  dir.create(file.path(workingdir,'processed',seriesName), showWarnings = FALSE)
  for (cid in cids){
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
      smoo1 <- smooth.CNA(cna1, smooth.region=10,smooth.SD.scale=undosd)
      message(paste("Processed CN segmentation for sample:",cid,"Chr:",chrname))

      seg1 <- segment(smoo1, min.width=5, verbose=0)
      ss1 <- segments.summary(seg1)[c(1:4,6,5)]
      write.table(ss1, file=fn, sep="\t", quote=FALSE,
                  append = T, row.names=FALSE, col.names=F)
      }
  }


  gc()


  ### remove the gaps of chromosomes ####
  for (cid in cids){
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
          emptyrow <- subfile[row,]
          emptyrow$start <- gapstart
          emptyrow$end <- gapend
          emptyrow[,5:ncol(emptyrow)] <- rep(0,ncol(emptyrow)-4)
          newrow <- subfile[row,]
          newrow$start <- gapend
          subfile[row,]$end  <- gapstart
          n <- subfile[row,]$probes
          subfile[row,]$probes <- n/2
          newrow$probes <- n/2
          if (row < nrow(subfile)) {

            subfile <- rbind(subfile[1:row,],emptyrow, newrow,subfile[(row+1):nrow(subfile),])
            next
          }
          else {
            subfile <- rbind(subfile[1:row,],emptyrow, newrow)
            next
          }
        }
        else if((subfile[row,]$start) < gapstart & subfile[row,]$end <= gapend &subfile[row,]$end > gapstart) {
          emptyrow <- subfile[row,]
          emptyrow$start <- gapstart
          emptyrow$end <- gapend
          emptyrow[,5:ncol(emptyrow)] <- rep(0,ncol(emptyrow)-4)
          subfile[row,]$end <-gapstart
          subfile <- rbind(subfile[1:row,],emptyrow, subfile[(row+1):nrow(subfile),])
        }
        else if ((subfile[row,]$start) >=gapstart & (subfile[row,]$start) < gapend & subfile[row,]$end > gapend) {
          emptyrow <- subfile[row,]
          emptyrow$start <- gapstart
          emptyrow$end <- gapend
          emptyrow[,5:ncol(emptyrow)] <- rep(0,ncol(emptyrow)-4)
          subfile[row,]$start <-gapend
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
