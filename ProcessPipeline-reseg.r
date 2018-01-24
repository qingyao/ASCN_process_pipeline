### aroma affymetrix processing pipeline ###

future::plan("multiprocess")
args = commandArgs(trailingOnly = TRUE)
workingdir <- args[1]
setwd(workingdir)
seriesName <- args[2]
sourcedir <- args[3]
filename <- args[4]
Pipelinedir <- file.path(sourcedir,'Rpipeline')
source(file.path(Pipelinedir,'common_function.r'))
source(file.path(Pipelinedir,'seg_eval.r'))
source(file.path(Pipelinedir,'final_BafAnalysis.r'))

log <- vector()
remotedir <- paste0("/Volumes/arraymapIncoming/aroma/aromaRaw/",seriesName)
chipTypes <- list.files(remotedir)
print(chipTypes)
time <- format(Sys.time(), "%m%d-%H%M%S")
for (chipType in chipTypes){
  if (!chipType %in% list.files(file.path(getwd(),"annotationData","chipTypes"))) break

  remotepath <- paste0(remotedir,"/",chipType)

  files<- list.files(remotepath)
  cids <- gsub(".CEL","",files)

  settings = list (
    seriesName = seriesName,
    chipType = chipType,
    workingdir = workingdir,
    remotedir = "/Volumes/arraymapMirror/arraymap/hg19",
    filename = filename
  )


  for (cid in cids) {
    localsettings <- settings
    localsettings[['arrayName']] <- cid
    if (filename == 'cn') {
        log <- c(log,tryCatch({
#          do.call(adjustMedian,localsettings)
          lmd <- do.call(getLmd,localsettings)
          gp <- do.call(getGP,localsettings)
          filtersettings <- c(localsettings,lmd=lmd,gp=gp)
          do.call(stepFilter,filtersettings)
          do.call(rmGaps,localsettings)
        }, error=function(e){
          message('Error! ',e,'\n')
          return(paste0("Error\t",format(Sys.time(), "%y-%m-%d %H:%M:%S"),"\t","\t",seriesName,"\t",e))}
        )
      )}
    else if (filename == 'fracb') {
      log <- c(log,tryCatch({
        do.call(BafAnalysis,localsettings)
      }, error=function(e){
        message('Error! ',e,'\n')
        return(paste0("Error\t",format(Sys.time(), "%y-%m-%d %H:%M:%S"),"\t","\t",seriesName,"\t",e))}
      )
    )
    }
  }
}
write.table(paste0(log),paste0(workingdir,"/processed/aroma_",format(Sys.time(), "%y-%m-%d"),".log"),quote=F,row.names = F,col.names = F,append=T)
