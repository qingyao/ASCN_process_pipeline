### aroma affymetrix processing pipeline ###

future::plan("multiprocess")
args = commandArgs(trailingOnly = TRUE)
workingdir <- args[1]
setwd(workingdir)
seriesName <- args[2]
sourcedir <- args[3]
force <- as.numeric(args[4])
filetype <- args[5]
Pipelinedir <- file.path(sourcedir,'Rpipeline')
source(file.path(Pipelinedir,'common_function.r'))
source(file.path(Pipelinedir,'final_seg.r'))

log <- vector()
localdir <- paste0(workingdir,"/rawData/",seriesName)
processedlocalpath <- paste0(workingdir,"/processed/",seriesName)
if (!dir.exists(processedlocalpath)) dir.create(processedlocalpath)
remotedir <- paste0("/Volumes/arraymapIncoming/aroma/aromaRaw/",seriesName)
chipTypes <- list.files(remotedir)
time <- format(Sys.time(), "%m%d-%H%M%S")
for (chipType in chipTypes){
  if (!chipType %in% list.files(file.path(getwd(),"annotationData","chipTypes"))) break

  if (dir.exists(localdir) == FALSE) dir.create(localdir)
  localpath <- paste0(localdir,"/",chipType)
  if (dir.exists(localpath) == FALSE) dir.create(localpath)
  remotepath <- paste0(remotedir,"/",chipType)

  files<- list.files(remotepath)
  remoteProcessPath <- file.path("/Volumes/arraymapMirror/arraymap/hg19",seriesName)
  cids <- gsub(".CEL","",files)

  settings = list (
    seriesName = seriesName,
    arrayName = NULL,
    chipType = chipType,
    workingdir = workingdir,
    undosd = 1
  )

  #check if all segments,cn.tsv are there for all arrays.
  i <- 1
  incomplete <- force
  notStarted <- 0
  if (force == 0) {
    newcids = vector()
    for (i in 1:length(cids)){
      if(!file.exists(file.path(remoteProcessPath,cids[i],chooseFile(filetype,1)[[1]] )) {
        newcids <- c(newcids,cids[i])
      }
    }
    cids <- newcids
    if (length(cids) > 0) incomplete <- 1
  }


  if (incomplete == 1) {
      for (cid in cids) {
        localsettings <- settings
        localsettings[['arrayName']] <- cid
        log <- c(log,tryCatch({do.call(get(sprintf('%ssegperArray',filetype)),localsettings)},error=function(e){
            message("Here's the original error message:")
            message(e,"\n")
            return(paste0("Error\t",format(Sys.time(), "%y-%m-%d %H:%M:%S"),"\t",sprintf("%s segmentation\t",toupper(filetype)),seriesName,"\t",e))}))
      }
    } else {
      next
    }

  }
write.table(paste0(log),paste0(workingdir,"/processed/aroma_",format(Sys.time(), "%y-%m-%d"),".log"),quote=F,row.names = F,col.names = F,append=T)
