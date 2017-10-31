### aroma affymetrix processing pipeline ###

future::plan("multiprocess")
args = commandArgs(trailingOnly = TRUE)
workingdir <- args[1]
setwd(workingdir)
seriesName <- args[2]
cleanup <- args[3]
sourcedir <- args[4]
memory <- args[5]
force <- as.numeric(args[6])
sdforce <- as.numeric(args[7])
Pipelinedir <- file.path(sourcedir,'Rpipeline')
sourcefiles <- list.files(Pipelinedir)
sapply(file.path(Pipelinedir,sourcefiles),source,.GlobalEnv)

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
    chipType = chipType,
    workingdir = workingdir,
    sourcedir = sourcedir,
    memory = memory,
    sdforce = sdforce
  )


    log <- c(log,tryCatch({do.call(cnseg,settings)},error=function(e){
        message("Here's the original error message:")
        message(e,"\n")
        return(paste0("Error\t",format(Sys.time(), "%y-%m-%d %H:%M:%S"),"\t","CNsegmentation\t",seriesName,"\t",e))}))




}
write.table(paste0(log),paste0(workingdir,"/processed/aroma_",format(Sys.time(), "%y-%m-%d"),".log"),quote=F,row.names = F,col.names = F,append=T)
