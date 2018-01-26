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
Pipelinedir <- file.path(sourcedir,'Rpipeline')
source(file.path(Pipelinedir,'final_ACNE_hg19.R'))
source(file.path(Pipelinedir,'final_doCRMAv2.R'))

log <- vector()
localdir <- paste0(workingdir,"/rawData/",seriesName)
processedlocalpath <- paste0(workingdir,"/processed/",seriesName)
remotedir <- paste0("/Volumes/arraymapIncoming/aroma/aromaRaw/",seriesName)
chipTypes <- list.files(remotedir)
time <- format(Sys.time(), "%m%d-%H%M%S")
for (chipType in chipTypes){
  if (!chipType %in% list.files(file.path(getwd(),"annotationData","chipTypes"))) next
  if (!chipType %in% c('Mapping50K_Xba240','Mapping50K_Hind240')) next
  if (dir.exists(localdir) == FALSE) dir.create(localdir)
  localpath <- paste0(localdir,"/",chipType)
  if (dir.exists(localpath) == FALSE) dir.create(localpath)
  remotepath <- paste0(remotedir,"/",chipType)

  files<- list.files(remotepath)
  remoteProcessPath <- file.path("/Volumes/arraymapMirror/arraymap/hg19",seriesName)
  cids <- gsub(".CEL","",files)

### check if processed probe files are complete
  checkIncomplete <- function(namestructure) {
    incomplete <- 0
    for(i in 1:length(cids)){
      if (incomplete == 0){
        for (chr in 1:23){
          if(!file.exists(file.path(remoteProcessPath,cids[i],sprintf(namestructure,chr)))) {
            incomplete <- 1
            break
          }
        }

      }
    }
    return(incomplete)
  }
  fracbIncomplete <- force | checkIncomplete('fracB,chr%s.tab')
  cnIncomplete <- force | checkIncomplete('probes,cn,chr%s.tsv')
  message("fracb",fracbIncomplete)
  message("cn",cnIncomplete)
  if (fracbIncomplete | cnIncomplete){
    for (file in files){
      print(file)
      file.copy(from=paste0(remotepath,"/",file), to=localpath, overwrite = FALSE, recursive = FALSE, copy.mode = TRUE)
    }
  }

  settings = list (
    seriesName = seriesName,
    chipType = chipType,
    workingdir = workingdir,
    sourcedir = sourcedir,
    memory = memory
  )

  if (fracbIncomplete) {

    log <- c(log,tryCatch({do.call(ACNE,settings)},error=function(e){
    message("Here's the original error message:")
    message(e,"\n")
    return(paste0("Error\t",format(Sys.time(), "%y-%m-%d %H:%M:%S"),"\t","ACNE\t",seriesName,"\t",e))}))

  }

  if (cnIncomplete) {

    log <- c(log,tryCatch({do.call(CRMAv2,settings)},error=function(e){
    message("Here's the original error message:")
    message(e,"\n")
    return(paste0("Error\t",format(Sys.time(), "%y-%m-%d %H:%M:%S"),"\t","CRMAv2\t",seriesName,"\t",e))}))

  }



  for (file in files){
    file.remove(paste0(localpath,"/",file))
  }

  if (cleanup==1) {
    tmpfiledir <- paste0(workingdir,c("/plmData/","/probeData/"))
    for (dir in tmpfiledir) {
      tmpfiles <- list.files(dir)
      tmpfiles <- tmpfiles[grep(seriesName,tmpfiles)]
      for (tmpfile in tmpfiles){
        unlink(file.path(tmpfiledir,tmpfile),recursive = TRUE)
      }
    }
  }


}
if (dir.exists(paste0(workingdir,"/processed/",seriesName)) == F) dir.create(paste0(workingdir,"/processed/",seriesName))
write.table(paste0(log),paste0(workingdir,"/processed/aroma_",format(Sys.time(), "%y-%m-%d"),".log"),quote=F,row.names = F,col.names = F,append=T)
