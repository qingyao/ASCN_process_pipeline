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


remotedir <- paste0("/Volumes/arraymapIncoming/aroma/aromaRaw/",seriesName)
chipTypes <- list.files(remotedir)
time <- format(Sys.time(), "%m%d-%H%M%S")
for (chipType in chipTypes){
  settings = list (
    seriesName = seriesName,
    arrayName = NULL,
    chipType = chipType,
    workingdir = workingdir,
    sourcedir = sourcedir,
    memory = memory,
    sdforce = sdforce
  )
  print(settings)
  cnseg(seriesName,NULL,chipType,workingdir,sourcedir,memory,sdforce)
}
