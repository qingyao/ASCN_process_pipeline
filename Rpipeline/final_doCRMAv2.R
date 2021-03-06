####################################
### Total copy number extraction ###
####################################
suppressWarnings(suppressMessages(library(aroma.affymetrix)))

CRMAv2  <- function(seriesName,chipType,workingdir,sourcedir,memory) {
  setOption(aromaSettings, "memory/ram", memory)
  future::plan("multiprocess")
  log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE)
  # Don't display too many decimals.
  options(digits=4)

  setwd(workingdir)
  ## get sampletype from metadata with python
  command <- "python3"
  path2script = sprintf("%s/getmeta.py",sourcedir)
  output <- system2(command,args = c(path2script,seriesName),stdout = TRUE)

  m <- gregexpr("'[A-Za-z]+'",output)
  output <- unlist(regmatches(output,m))
  output <- gsub("'","",output)

  if (length(which(output == "Normal")) >= 10) {
    useExtRef <- FALSE
    refIdx <- which(output=="Normal")
  } else useExtRef <- TRUE

  if (chipType %in% c("GenomeWideSNP_6","GenomeWideSNP_5")){cdf <- AffymetrixCdfFile$byChipType(chipType, tags="Full")
  }else {cdf <- AffymetrixCdfFile$byChipType(chipType)}
  #print(cdf)

  gi <- getGenomeInformation(cdf)
  #print(gi)

  si <- getSnpInformation(cdf)
  #print(si)

  acs <- AromaCellSequenceFile$byChipType(getChipType(cdf, fullname=FALSE))
  #print(acs)


  csR <- AffymetrixCelSet$byName(seriesName, cdf=cdf)
  #print(csR)

  acc <- AllelicCrosstalkCalibration(csR, model="CRMAv2")

  csC <- process(acc, verbose=verbose)

  bpn <- BasePositionNormalization(csC, target="zero")
  # print(bpn)

  csN <- process(bpn, verbose=verbose)
  # print(csN)

  plm <- AvgCnPlm(csN, mergeStrands=TRUE, combineAlleles=TRUE)
  #print(plm)

  if (length(findUnitsTodo(plm)) > 0) {
    # Fit CN probes quickly (~5-10s/array + some overhead)
    units <- fitCnProbes(plm, verbose=verbose)
#    str(units)
    # int [1:945826] 935590 935591 935592 935593 935594 935595 ...

    # Fit remaining units, i.e. SNPs (~5-10min/array)
    units <- fit(plm, verbose=verbose)
#    str(units)
  }

  ces <- getChipEffectSet(plm)
  #print(ces)

  fln <- FragmentLengthNormalization(ces, target="zero")
  #print(fln)

  cesN <- process(fln, verbose=verbose) #26h for GSE19949
  #print(cesN)

  if (useExtRef == FALSE) {ceR <- getAverageFile(cesN[refIdx], verbose=verbose)
  #save(ceR,file=file.path(getwd(),"ReferenceFile",chipType,"ceR.RData"))
  }else {load(file.path(getwd(),"ReferenceFile",chipType,"ceR.RData")) }
  #print(ceR)

  # ce <- cesN[[3]]  # Array #3
  # theta <- extractTheta(ce, unit=987)
  # thetaR <- extractTheta(ceR, unit=987)
  # C <- 2*theta/thetaR
  # print(C)
  #
  # print(c(theta, thetaR))


  path <- file.path(getwd(),"processed",seriesName)
  if (dir.exists(path) == 0) dir.create(path)

  # Write out logR in a file
  for (chr in 1:23) {
    # Define the chromosome
    cdf <- getCdf(cesN)
    gi <- getGenomeInformation(cdf)
    units <- getUnitsOnChromosome(gi, chromosome=chr)
    pos <- getPositions(gi, units=units)
    pos <- as.numeric(pos)
    unitnames <- getUnitNames(cdf,units=units)

    # Retrieving CN estimates of the reference in the above region
    if (useExtRef == TRUE) {
      thetaR <- extractTheta(saveLoadReference,units=units)
    } else {
      thetaR <- extractTheta(ceR, units=units)
    }

    # Retrieving the corresponding estimates of samples
    sampleIDs <- gsub(".CEL","",list.files(paste(getwd(),'rawData',seriesName,chipType,sep="/")))
    for (i in 1:length(sampleIDs)) {
      ce <- cesN[[indexOf(cesN, sampleIDs[i])]]
      theta <- extractTheta(ce, units=units)
      # Calculate the relative CNs
      R <- theta/thetaR
      logR <- round(log2(R),4)
      out <- data.frame(cbind(unitnames,rep(chr,length(pos)),pos,logR),stringsAsFactors = F)
      out <- out[order(as.numeric(out$pos)),]
      colnames(out) <- c("PROBEID","CHRO","BASEPOS","VALUE")
      samplepath <- file.path(path,sampleIDs[i])
      if (dir.exists(samplepath) == 0) dir.create(samplepath)
      write.table(out,file.path(samplepath,sprintf("probes,cn,chr%s.tsv",chr)),sep = "\t", quote = F,row.names = F)
      cat(sprintf('%s/%s\n',seriesName,sampleIDs[i]),file=file.path(getwd(),"processed",'logs','log_CRMAv2.txt'),append=T)
    }

  }

}
