########################
### CN segmentation ####
########################


cnsegPerArray <- function(workingdir,seriesName, cid, chrs, minw, smrg, smsd, undosplit, undopr, undosd){
  dir.create(file.path(workingdir,'processed',seriesName,cid), showWarnings = FALSE)
  fn <- file.path(workingdir,'processed',seriesName,cid, sprintf('segments,cn,%s_%s_%s_%s_%s.tsv',minw,smrg, smsd, undosplit, undopr,undosd))
  cat("sample_id","chromosome","start", "end", "value", "probes\n",sep="\t",file=fn,append = F)

  for (chrname in chrs){
    data<- read.table(sprintf('/Volumes/arraymapMirror/arraymap/hg19/%s/%s/probes,cn,chr%d.tsv',seriesName,cid,chrname),header=T)
    cna1 <- CNA(genomdat=data$VALUE,
                chrom=rep(chrname, length(data$VALUE)),
                maploc=data$BASEPOS,
                data.type="logratio",
                sampleid=cid)
    time <- system.time(smoo1 <- smooth.CNA(cna1, smooth.region=smrg,smooth.SD.scale=smsd)
    message(paste("Processed CN segmentation for sample:",cid,"Chr:",chrname))

    seg1 <- segment(smoo1, min.width= minw, verbose=0, undo.splits = undosplit, undo.prune = undoprune,undo.SD = undosd))
    ss1 <- segments.summary(seg1)[c(1:4,6,5)]
    write.table(ss1, file=fn, sep="\t", quote=FALSE,
                append = T, row.names=FALSE, col.names=F)
    logf <- file.path(workingdir,'processed',seriesName,cid, 'log.txt')
    cat(nrow(ss1),sum(segments.p(seg1)$pval < 0.1), time, minw,smrg, smsd, undosplit, undopr,undosd, file = logf,append = TRUE)
    }
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


  aft <- c(90,125,90,125,125,100,90,85,70,80,75,85,50,75,75,40,50,55,25,30,30,30,75)
  bef <- c(110,80,80,40,40,50,55,37,37,37,40,30,15,15,30,20,12,20,20,5,30,)

  }

}
