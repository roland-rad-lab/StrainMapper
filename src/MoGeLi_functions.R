suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(splitstackshape))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(tidytext))

# STRAIN - FUNCTIONS ----------

# this ended up being some premium spaghetti code, sorry about that :)

# load data from file, this is specific to our local files
LoadData <- function(samplePath, inputType, strainDT, selectFiles="TumorNormal", minCoverage=1, minAF=0.2, minMQ=58, minBQ=20, combineTN="all", removeClusteredEvents=T, maxNeighbors=1){
  if(inputType=="MoCaSeq"){
    # classic MoCaSeq results (this will expects the defined folder structure of MoCaSeq)
    
    samplename <- basename(samplePath)
    
    if(selectFiles == "TumorNormal"){
      tumor.file <- paste(samplePath,"/results/Mutect2/",samplename,".Tumor.Mutect2.Positions.txt", sep="")
      normal.file <- paste(samplePath,"/results/Mutect2/",samplename,".Normal.Mutect2.Positions.txt", sep="")
      if(!file.exists(tumor.file)){stop(paste0("File not found: ", tumor.file))}
      if(!file.exists(normal.file)){stop(paste0("File not found: ", tumor.file))}
      
      tumor <- fread(tumor.file, header=T, sep="\t")
      normal <- fread(normal.file, header=T, sep="\t")
      
      # Tumor and Normal can be combined in 2 ways: 
      # all mutations from both (more noise but also more data/coverage) or 
      # only if found in both (less noise + cleam picture, but maybe information loss)
      if(combineTN=="all"){
        sampleSNPs <- merge(tumor, normal, all=T, by=c("CHROM", "POS", "REF", "ALT"))
      } else if(combineTN=="shared"){
        sampleSNPs <- merge(tumor, normal, all=F, by=c("CHROM", "POS", "REF", "ALT"))
      }
      
      sampleSNPs.unfiltered <- sampleSNPs[, .(Chr=CHROM, Pos=POS, Ref=REF, Alt=ALT)]
      
      # set filters on read coverage, allele frequency and mapping quality
      # in this approach, if ONE OF BOTH passes the filter criteria, the SNP is used (this is less stringend)
      sampleSNPs <- sampleSNPs[`GEN[Tumor].AD[1]` >= minCoverage | `GEN[Normal].AD[1]` >= minCoverage]
      sampleSNPs <- sampleSNPs[`GEN[Tumor].AF` >= minAF | `GEN[Normal].AF` >= minAF]
      sampleSNPs <- sampleSNPs[`MMQ[1].x` >= minMQ | `MMQ[1].y` >= minMQ]
      
    } else if(selectFiles == "tumor" | selectFiles == "Tumor"){
      
      tumor.file <- paste(samplePath,"/results/Mutect2/",samplename,".Tumor.Mutect2.Positions.txt", sep="")
      if(!file.exists(tumor.file)){stop(paste0("File not found: ", tumor.file))}
      
      sampleSNPs <- fread(tumor.file, header=T, sep="\t")
      sampleSNPs.unfiltered <- sampleSNPs[, .(Chr=CHROM, Pos=POS, Ref=REF, Alt=ALT)]
      
      # set filters on read coverage, allele frequency and mapping quality
      sampleSNPs <- sampleSNPs[`GEN[Tumor].AD[1]` >= minCoverage]
      sampleSNPs <- sampleSNPs[`MMQ[1]` >= minMQ]
      sampleSNPs <- sampleSNPs[`GEN[Tumor].AF` >= minAF]
    } else if(selectFiles == "normal" | selectFiles == "Normal"){
      normal.file <- paste(samplePath,"/results/Mutect2/",samplename,".Normal.Mutect2.Positions.txt", sep="")
      if(!file.exists(normal.file)){stop(paste0("File not found: ", normal.file))}
      
      sampleSNPs <- fread(normal.file, header=T, sep="\t")
      sampleSNPs.unfiltered <- sampleSNPs[, .(Chr=CHROM, Pos=POS, Ref=REF, Alt=ALT)]
      
      # set filters on read coverage, allele frequency and mapping quality
      sampleSNPs <- sampleSNPs[`GEN[Normal].AD[1]` >= minCoverage]
      sampleSNPs <- sampleSNPs[`MMQ[1]` >= minMQ]
      sampleSNPs <- sampleSNPs[`GEN[Normal].AF` >= minAF]
    }
    
    sampleSNPs <- sampleSNPs[, .(Chr=CHROM, Pos=POS, Ref=REF, Alt=ALT)]
    
  } else if(inputType == "vcf"){
    # custom VCF approach (e.g. lcWGS) --> better use TXT if available (much faster)
    
    # this can be either a single file or the partial path for Tumor+Normal (PATH/SAMPLE.Tumor.vcf.gz and PATH/SAMPLE.Normal.vcf.gz should be given as PATH/SAMPLE)
    if(selectFiles == "TumorNormal" & !grepl("\\.vcf", samplePath)){
      samplename <- basename(samplePath)
      
      tumor.file <- paste(samplePath,".Tumor.vcf.gz", sep="")
      if(!file.exists(tumor.file)){stop(paste0("File not found: ", tumor.file))}
      normal.file <- paste(samplePath,".Normal.vcf.gz", sep="")
      if(!file.exists(normal.file)){stop(paste0("File not found: ", tumor.file))}
      
      tumor.vcfDat <- as.data.table(read.table(tumor.file, header = F, comment.char = "#"))
      tumor.vcfDat <- tumor.vcfDat[, .(Chr=V1, Pos=V2, Ref=V4, Alt=V5, Info.Tumor=V8)]
      
      normal.vcfDat <- as.data.table(read.table(normal.file, header = F, comment.char = "#"))
      normal.vcfDat <- normal.vcfDat[, .(Chr=V1, Pos=V2, Ref=V4, Alt=V5, Info.Normal=V8)]
      
      # Tumor and Normal can be combined in 2 ways: 
      # all mutations from both (more noise but also more data/coverage) or 
      # only if found in both (less noise + cleam picture, but maybe information loss)
      if(combineTN=="all"){
        print("Processing with all SNPs from Tumor and Normal, this is very slow! Using mpileup on Tumor.bam and Normal.bam to generate a single combined VCF is recommended.")
        vcfDat <- merge(tumor.vcfDat, normal.vcfDat, all=T, by=c("Chr", "Pos", "Ref", "Alt"))
      } else if(combineTN=="shared"){
        vcfDat <- merge(tumor.vcfDat, normal.vcfDat, all=F, by=c("Chr", "Pos", "Ref", "Alt"))
      }
      
      sampleSNPs.unfiltered <- vcfDat[, .(Chr, Pos, Ref, Alt)]
      
      # FILTERING
      vcfDat[, MQ.Tumor := gsub(".*MQ=(.*?);.*", "\\1",Info.Tumor)]
      vcfDat[, DP4.Tumor := gsub(".*DP4=(.*?);.*", "\\1",Info.Tumor)]
      vcfDat[, MQ.Normal := gsub(".*MQ=(.*?);.*", "\\1",Info.Normal)]
      vcfDat[, DP4.Normal := gsub(".*DP4=(.*?);.*", "\\1",Info.Normal)]
      
      vcfDat <- suppressWarnings(cSplit(vcfDat, c("DP4.Tumor", "DP4.Normal"), ",", type.convert = T))
      
      # DP4_1 and DP4_2 are the reference counts (on forward and reverse reads), therefore DP4_3 and DP4_4 are the alternative counts
      vcfDat[, DP.Tumor := DP4.Tumor_1+DP4.Tumor_2+DP4.Tumor_3+DP4.Tumor_4]
      vcfDat[, AF.Tumor := (DP4.Tumor_3+DP4.Tumor_4)/DP.Tumor]
      vcfDat[, DP.Normal := DP4.Normal_1+DP4.Normal_2+DP4.Normal_3+DP4.Normal_4]
      vcfDat[, AF.Normal := (DP4.Normal_3+DP4.Normal_4)/DP.Normal]
      
      vcfDat <- vcfDat[DP.Tumor >= minCoverage | DP.Normal >= minCoverage]
      vcfDat <- vcfDat[AF.Tumor >= minAF | AF.Normal >= minAF]
      vcfDat <- vcfDat[MQ.Tumor >= minMQ | MQ.Normal >= minMQ]
      
      sampleSNPs <- vcfDat[, .(Chr, Pos, Ref, Alt)]
      
    } else {
      if(!file.exists(samplePath)){stop(paste0("File not found: ", samplePath))}
      samplename <- gsub(".vcf.gz", "", basename(samplePath))

      vcfDat <- fread(samplePath, skip = "#CHROM", fill = T)
      sampleSNPs.unfiltered <- vcfDat[, .(Chr=V1, Pos=V2, Ref=V4, Alt=V5)]
      
      vcfDat <- vcfDat[(grepl("DP4=", V8) & grepl("MQ=", V8))]
      
      if(nrow(vcfDat) == 0){stop("No SNPs found, probably the required DP4 and MQ values in the VCF file were not found")}
      
      vcfDat[, MQ := gsub(".*MQ=(.*?);.*", "\\1",V8)]
      vcfDat[, DP4 := gsub(".*DP4=(.*?);.*", "\\1",V8)]
      suppressWarnings(vcfDat <- cSplit(vcfDat, "DP4", ",", type.convert = T))
      
      # DP4_1 and DP4_2 are the reference counts (on forward and reverse reads), therefore DP4_3 and DP4_4 are the alternative counts
      vcfDat[, DP := DP4_1+DP4_2+DP4_3+DP4_4]
      vcfDat[, AF := (DP4_3+DP4_4)/DP]
      
      vcfDat <- vcfDat[MQ >= minMQ]
      vcfDat <- vcfDat[DP >= minCoverage]
      vcfDat <- vcfDat[AF >= minAF]
      sampleSNPs <- vcfDat[, .(Chr=V1, Pos=V2, Ref=V4, Alt=V5)]
    }
  } else if(inputType == "txt"){
    
    # this can be either a single file or the partial path for Tumor+Normal (PATH/SAMPLE.Tumor.vcf.gz and PATH/SAMPLE.Normal.vcf.gz should be given as PATH/SAMPLE)
    if(selectFiles == "TumorNormal" & !grepl("\\.txt", samplePath)){
      
      samplename <- basename(samplePath)
      
      tumor.file <- paste(samplePath,".Tumor.txt.gz", sep="")
      if(!file.exists(tumor.file)){stop(paste0("File not found: ", tumor.file))}
      normal.file <- paste(samplePath,".Normal.txt.gz", sep="")
      if(!file.exists(normal.file)){stop(paste0("File not found: ", normal.file))}
      
      tumor <- fread(tumor.file, header=T, sep="\t")
      tumor <- tumor[POS != 0] # remove # header lines
      setnames(tumor, c("DP","DP4", "MQ"), c("DP.Tumor","DP4.Tumor", "MQ.Tumor"))
      normal <- fread(normal.file, header=T, sep="\t")
      normal <- normal[POS != 0] # remove # header lines
      setnames(normal, c("DP", "DP4", "MQ"), c("DP.Normal", "DP4.Normal", "MQ.Normal"))
      
      # Tumor and Normal can be combined in 2 ways: 
      # all mutations from both (more noise but also more data/coverage) or 
      # only if found in both (less noise + cleam picture, but maybe information loss)
      if(combineTN=="all"){
        print("Processing with all SNPs from Tumor and Normal, this is very slow! Using mpileup on Tumor.bam and Normal.bam to generate a single combined TXT is recommended.")
        txtDat <- merge(tumor, normal, all=T, by=c("CHROM", "POS", "REF", "ALT"))
      } else if(combineTN=="shared"){
        txtDat <- merge(tumor, normal, all=F, by=c("CHROM", "POS", "REF", "ALT"))
      }
      
      sampleSNPs.unfiltered <- txtDat[, .(Chr=CHROM, Pos=POS, Ref=REF, Alt=ALT)]
      
      suppressWarnings(txtDat <- cSplit(txtDat, c("DP4.Tumor", "DP4.Normal"), ",", type.convert = T))
      txtDat[, AF.Tumor := (DP4.Tumor_3+DP4.Tumor_4)/DP.Tumor]
      txtDat[, AF.Normal := (DP4.Normal_3+DP4.Normal_4)/DP.Normal]
      
      txtDat <- txtDat[DP.Tumor >= minCoverage | DP.Normal >= minCoverage]
      txtDat <- txtDat[AF.Tumor >= minAF | AF.Normal >= minAF]
      txtDat <- txtDat[MQ.Tumor >= minMQ | MQ.Normal >= minMQ]
      
      sampleSNPs <- txtDat[, .(Chr=CHROM, Pos=POS, Ref=REF, Alt=ALT)]
    } else {
      samplename <- gsub(".txt.gz", "", basename(samplePath))
      txtDat <- fread(samplePath)
      txtDat <- txtDat[POS != 0] # remove # header lines
      sampleSNPs.unfiltered <- txtDat[, .(Chr=CHROM, Pos=POS, Ref=REF, Alt=ALT)]
      
      suppressWarnings(txtDat <- cSplit(txtDat, "DP4", ",", type.convert = T))
      txtDat[, AF := (DP4_3+DP4_4)/DP]
      
      txtDat <- txtDat[MQ >= minMQ]
      txtDat <- txtDat[DP >= minCoverage]
      txtDat <- txtDat[AF >= minAF]
      sampleSNPs <- txtDat[, .(Chr=CHROM, Pos=POS, Ref=REF, Alt=ALT)]
    }
    
    
  } else if(inputType == "data"){
    samplename <- "custom"
    sampleSNPs <- copy(sampleSNPs) # this requires a global variable!
    sampleSNPs.unfiltered <- copy(sampleSNPs)
    
    if(!identical(colnames(sampleSNPs), c("Chr", "Pos", "Ref", "Alt"))){
      stop("Wrong input data format! If you use inputType==data the global variable sampleSNPs should have the 4 tab-separated columns (this order and naming, chromosomes without chr): 
           Chr Pos Ref Alt
           1 105336683 T C
           X 87691002 G A")
    }
    
  } else if(inputType == "BRC"){
    samplename <- gsub(".strains.brc.tsv","",basename(samplePath))
    
    if(!file.exists(samplePath)){
      stop("For the inputType=BRC the given samplePath should be the path as well as filename to the X.strain.brc.tsv from PrintCommandsBRC()")
    }
    BRC.unfiltered <- ProcessBRC(samplePath, minMQ=0, minReads=1, minBQ=0, minAF=0) # this can be used for the plot data itself
    
    if(is.null(BRC.unfiltered)){
      return(NULL)
    }
    
    sampleSNPs.unfiltered <- BRC.unfiltered[Ref != Alt, .(Chr, Pos, Ref, Alt)]
    
    BRC <- ProcessBRC(samplePath, minMQ=minMQ, minReads=minCoverage, minBQ=minBQ, minAF=minAF) # this can be used for the plot data itself
    sampleSNPs <- BRC[Ref != Alt, .(Chr, Pos, Ref, Alt)]
  } else {
    stop("Invalid inputType argument!")
  }
  
  # sub-select the MNVs (2x or Nx substitution), split those into separate lines and remerge into SNP dt
  MNVs <- sampleSNPs[nchar(Ref) == nchar(Alt) & nchar(Alt) != 1]
  if(nrow(MNVs) > 0){
    MNVs[, Ref := gsub("(?<=.)(?=.)", ",", Ref, perl = TRUE)]
    MNVs[, Alt := gsub("(?<=.)(?=.)", ",", Alt, perl = TRUE)]
    MNVs[, Pos := paste0(Pos,",",Pos+1)] # also increase the second base by 1
    MNVs <- cSplit(MNVs, c("Ref", "Alt", "Pos"), sep=",", direction = "long", type.convert = F)
    MNVs <- MNVs[!is.na(Ref)]
    MNVs <- MNVs[!is.na(Pos)]
    MNVs[, Pos := as.numeric(Pos)]
    sampleSNPs <- rbind(sampleSNPs, MNVs)
  }
  
  # now delete all the indels (and MNVs which were split above)
  sampleSNPs <- sampleSNPs[nchar(Ref) + nchar(Alt) == 2]
  
  # redo for unfiltered, to find all neighbors for "ClusteredEvents"
  MNVs.unfiltered <- sampleSNPs.unfiltered[nchar(Ref) == nchar(Alt) & nchar(Alt) != 1]
  sampleSNPs.unfiltered <- sampleSNPs.unfiltered[!(nchar(Ref) == nchar(Alt) & nchar(Alt) != 1)] # already remove from input list (but keep indels which should be counted as neighbors)
  if(nrow(MNVs.unfiltered) > 0){
    MNVs.unfiltered[, Ref := gsub("(?<=.)(?=.)", ",", Ref, perl = TRUE)]
    MNVs.unfiltered[, Alt := gsub("(?<=.)(?=.)", ",", Alt, perl = TRUE)]
    MNVs.unfiltered[, Pos := paste0(Pos,",",Pos+1)] # also increase the second base by 1
    MNVs.unfiltered <- cSplit(MNVs.unfiltered, c("Ref", "Alt", "Pos"), sep=",", direction = "long", type.convert = F)
    MNVs.unfiltered <- MNVs.unfiltered[!is.na(Ref)]
    MNVs.unfiltered <- MNVs.unfiltered[!is.na(Pos)]
    MNVs.unfiltered[, Pos := as.numeric(Pos)]
    sampleSNPs.unfiltered <- rbind(sampleSNPs.unfiltered, MNVs.unfiltered)
  }
  
  # removeClusteredEvents can not be done for empty data
  if(nrow(sampleSNPs.unfiltered) == 0){removeClusteredEvents <- F}
  
  if(removeClusteredEvents){
    blacklistedIDs <- BlacklistClusteredEvents(sampleSNPs.unfiltered, strainDT, maxNeighbors = maxNeighbors) 
  } else {
    blacklistedIDs <- c()
  }
  
  return(list(sampleSNPs, samplename, blacklistedIDs))
}

# read files and merge with strain signatures
MergeDataWithStrains <- function(samplePath, strainDT, GenomeBins, genomeSNPs, returnData=F, inputType="vcf", selectFiles="TumorNormal", minCoverage=1, minAF=0.2, minMQ=58, minBQ=20, removeClusteredEvents=T, maxNeighbors=1, minBinSNPS=2, combineTN="all", addData=NULL){

  GenomeBinsGR <- makeGRangesFromDataFrame(GenomeBins)
  
  snpDT <- LoadData(samplePath, inputType, strainDT, selectFiles, minCoverage, minAF, minMQ, minBQ, combineTN, removeClusteredEvents, maxNeighbors)

  # highly customized runs
  if(!is.null(addData)){
    snpDT[[1]] <- rbind(snpDT[[1]], addData[, .(Chr, Pos, Ref, Alt)]) # debug test
  }
  
  # will be returned instead of NULL
  emptyReturn <- list(strainBinsAll=NULL, strainBinsRaw=NULL, sampleName=snpDT[[2]])
  
  if(is.null(snpDT)){
    print("No SNPs were found, please double check your data.")
    return(emptyReturn)
  } else if(nrow(snpDT[[1]]) == 0){
    print("No SNPs were found, please double check your data.")
    return(emptyReturn)
  }
  
  blacklistedIDs <- snpDT[[3]]
  samplename <- snpDT[[2]]
  snpDT <- snpDT[[1]]
  
  strainDT.filtered <- copy(strainDT) # this is called filtered because a million years ago there was some filtering, now this is just legacy naming convention
  
  if(removeClusteredEvents){
    strainDT.filtered <- strainDT.filtered[!ID %in% blacklistedIDs]
  }
  
  snpDT[, Chr := as.character(Chr)]
  snpDT <- merge(snpDT, strainDT.filtered, by=c("Chr", "Pos", "Ref", "Alt")) # identify snps specific for a strain
  
  if(is.null(snpDT)){
    print("No SNPs were found after filtering, please double check your data.")
    return(emptyReturn)
  } else if(nrow(snpDT) == 0){
    print("No SNPs were found after filtering, please double check your data.")
    return(emptyReturn)
  }
  
  if(nrow(snpDT) == 0){
    stop(paste0("No strain SNPs found for : ", dir))
  }
  
  snpGR <- makeGRangesFromDataFrame(snpDT[, .(chr=Chr, start=Pos, end=Pos)])
  hits <- findOverlaps(GenomeBinsGR, snpGR)
  
  # prepare some QC data (same as below but just with all the invidivual SNPs)
  hitStats <- NULL
  if(returnData){
    hitStats <- cbind(GenomeBins[queryHits(hits), -c("plotstart", "plotend")], snpDT[subjectHits(hits), .(Pos, Ref, Alt, ID, strainGroup, strain)])
    hitStats[, nFound := .N, by=.(UnfiedPosition, strain)]
    hitStats <- merge(hitStats, genomeSNPs[, .(UnfiedPosition, nExpected=nTotal, strain)], by=c("UnfiedPosition", "strain"), all.x=T, sort=F)
    hitStats[, binValue := nFound / nExpected]
    hitStats[, sample := samplename]
  }
  
  # calculate the value per bin
  strainBins <- cbind(GenomeBins[queryHits(hits)], snpDT[subjectHits(hits), .(strain)])
  
  strainBins <- strainBins[, .N, by=.(UnfiedPosition, strain)]
  
  # filtering out bins with only a small number of SNPs
  strainBins <- strainBins[N >= minBinSNPS]
  
  # normalize by number of total snps of that strain in that bin 
  strainBins <- merge(strainBins, genomeSNPs[, .(UnfiedPosition, nTotal, strain)], by=c("UnfiedPosition", "strain"), all.x=T, sort=F)
  strainBins[, value := N / nTotal]
  strainBins[, sample := samplename]
  
  return(list(strainBinsAll=strainBins, strainBinsRaw=hitStats, sampleName=samplename))
}

# will "count" the number of strain SNPs per bin, based on all, best (only best strain per bin is kept), second, both
BinScoring <- function(strainBinsAll, keepHits="all"){
  # select and score bins
  if(keepHits == "best"){
    # get only the best hit for each bin
    sampleResults <- strainBinsAll[strainBinsAll[, .I[value == max(value)], by=.(UnfiedPosition, sample)]$V1]
  } else if(keepHits == "second"){
    sampleResults <- strainBinsAll[strainBinsAll[, .I[value == sort(value, decreasing = T)[2]], by=.(UnfiedPosition, sample)]$V1]
  } else if(keepHits == "both"){
    # second and first best hit
    strainBinsTop1 <- strainBinsAll[strainBinsAll[, .I[value == sort(value, decreasing = T)[1]], by=.(UnfiedPosition, sample)]$V1]
    strainBinsTop2 <- strainBinsAll[strainBinsAll[, .I[value == sort(value, decreasing = T)[2]], by=.(UnfiedPosition, sample)]$V1]
    sampleResults <- rbind(strainBinsTop1, strainBinsTop2)
  } else if(keepHits == "all"){
    # all hits
    sampleResults <- copy(strainBinsAll)
  }
  sampleResults[, BinScoring := keepHits]
  return(sampleResults)
}

# do some postprocessing and calculations
# minGES = the minimal percentage of genome (i.e. number of bins) to be shown (e.g. 1 bin would be 0.004%, which is removed)
PostProcessStrainBins <- function(sampleResults, samplename, strainGroupAssignment, nBins, minBES, minGES=1, keepEmptyStrains=F, roundDigits=0){
  keepHits <- sampleResults[, unique(BinScoring)]
  
  # rename strain groups into meaningful names
  sampleResults <- merge(sampleResults, strainGroupAssignment, by="strain", all.x=T)
  sampleResults[, strain := groupname]
  
  # filter by a minimal value score
  sampleResults <- sampleResults[value > minBES]
  
  # calculate "strain-ness" score for each sample
  # naive: percentage of bins belonging to a strain
  spDT <- sampleResults[, .N, by=.(sample, strain)]
  spDT <- spDT[, .(SP=(N/nBins*100), strain), by=sample] # strain percent
  
  nEmptyBins <- sampleResults[, length(unique(UnfiedPosition)), by=sample]
  nEmptyBins[, V1 := ((nBins-V1) / nBins) * 100]
  nEmptyBins[, strain := NA]
  setnames(nEmptyBins, "V1", "SP")
  spDT <- rbind(spDT, nEmptyBins)
  
  # also add the "NA" strain if there are 0 hits found for a sample
  if(nrow(sampleResults) == 0){
    missingSamples <- data.table(sample = samplename, SP=100, strain=NA)
    spDT <- rbind(spDT, missingSamples)
    keepHits <- "skip"
  }
  
  spDT[, SP := round(SP, digits = roundDigits)] # round now, to prevent rounding errors
  spDT[, SPtext := paste0(SP, "%")]
  spSummary <- copy(spDT)
  spSummary[SP < 5, strain := "others"]
  spSummary <- spSummary[, .(SP=sum(SP)), by=.(sample, strain)]
  setorder(spSummary, sample, -SP)
  
  # remove bins with a single bin in the entire genome (0% if rounded) and set bin to NA
  strain0 <- spDT[SP <= minGES, paste0(strain, sample)]
  sampleResults[, tmpID := paste0(strain, sample)]
  sampleResults <- sampleResults[!tmpID %in% strain0]
  sampleResults[, tmpID := NULL]
  spDT[, tmpID := paste0(strain, sample)]
  
  # calculate how much is removed and add it to NA
  if(keepHits %in% c("best", "second")){
    addToNA <- spDT[tmpID %in% strain0, sum(SP)]
    spDT <- spDT[!tmpID %in% strain0]
    spDT[, tmpID := NULL]
    spDT[is.na(strain), SP := SP + addToNA]
    
  } else {
    # get the updated number of empty bins
    nEmptyBins <- sampleResults[, length(unique(UnfiedPosition)), by=sample]
    nEmptyBins[, V1 := ((nBins-V1) / nBins) * 100]
    nEmptyBins[, strain := NA]
    nEmptyBins[, SP := round(V1, digits = 0)]
    
    addNA <- merge(spDT[is.na(strain), -c("SP")], nEmptyBins[, .(sample, SP)], by="sample", all.x=T)
    addNA[is.na(SP), SP := 100]
    
    spDT <- spDT[!tmpID %in% strain0] # remove strain if only supported by x (default 1) bins
    spDT <- rbind(spDT[!is.na(strain)], addNA)
    spDT[, tmpID := NULL]
  }
  spDT[, SPtext := paste0(SP, "%")]
  
  # this will add empty facets to the plot
  if(keepEmptyStrains){
    spAdd <- data.table()
    for(mysample in spSummary[, unique(sample)]){
      
      allStrains <- strainDT[, unique(strain)] 
      missingStrains <- allStrains[!allStrains %in% spDT[sample == mysample, unique(strain)]]
      
      tmp <- data.table(sample=mysample, SP=0, strain=missingStrains, SPtext="0%")
      spAdd <- rbind(spAdd, tmp)
    }
    spDT <- rbind(spDT, spAdd)
  }
  
  spDT <- merge(spDT, strainGroupAssignment[, .(strain, strainGroup)], all.x=T, by="strain")
  
  return(list(sampleResults, spDT, spSummary))
}

# Plotting procedures for bubbles
PlotStrainBins <- function(sampleResults, spDT, samplename, GenomeBins, strainGroupColors, annotationLoci=NULL, returnData=F, excludeNA=F){
  
  # some data for plotting
  chromEnds <- GenomeBins[, max(plotend), by=chr]$V1 # equal to the original "Len" variable
  chromEnds <- c(0, chromEnds)
  chromnames <- GenomeBins[, paste0("chr", unique(chr))]
  chromLabelPos <- diff(chromEnds)/2
  names(chromLabelPos) <- chromnames
  chromLabelPos <- chromEnds+c(chromLabelPos, 0) # with pseuso 0 after last chromosome
  GenomeBinsGR <- makeGRangesFromDataFrame(GenomeBins)
  
  plotDT <- merge(GenomeBins, sampleResults[, .(UnfiedPosition, strain, N, value, sample, strainGroup)], all.x=T, by="UnfiedPosition")
  
  if(excludeNA){
    plotDT <- plotDT[!is.na(strain)]
    spDT <- spDT[!is.na(strain)]
  }
  
  spDT.plot <- copy(spDT) # work copy
  # remove NA if the entire genome is covered by SNPs (this would break the % text)
  if(all(spDT.plot[is.na(strain), SP] %in% 0)){
    spDT.plot <- spDT.plot[!is.na(strain)]
  }
  
  # assign the max possible bubble size to the NA bubbles (if to prevent warnings for full NA samples)
  if(nrow(plotDT[!is.na(value)]) == 0){
    maxValue <- 1 # if all bins are NA
  } else {
    maxValue <- plotDT[, max(value, na.rm = TRUE)]
  }
  
  # legacy code
  naClasses <- NA
  plotDT[is.na(strain), value := maxValue] 
  
  sub.spDT <- spDT.plot[!strain %in% naClasses]
  
  sub.spDT <- sub.spDT[!grepl("\\(", strain)]
  setorder(sub.spDT, SP)
  strainOrder <- as.character(sub.spDT$strain)
  
  todos.strains <- unique(spDT.plot[grepl("\\(", strain), gsub(" \\(.*\\)", "", strain)])
  for(strainname in todos.strains){
    addIndex <- which(strainOrder==strainname)-1
    strainOrder <- append(strainOrder, paste0(strainname, " (filtered)"), after=addIndex)
    strainOrder <- append(strainOrder, paste0(strainname, " (no coverage)"), after=addIndex)
  }
  
  plotDT[, strain := factor(strain, levels = c(strainOrder, naClasses))]
  spDT.plot[, strain := factor(strain, levels = c(strainOrder, naClasses))]
  
  setorder(sub.spDT, -SP)
  straingroupOrder <- unique(as.character(sub.spDT$strainGroup))
  plotDT[, strainGroup := factor(strainGroup, levels = c(straingroupOrder, naClasses))]
  spDT.plot[, strainGroup := factor(strainGroup, levels = c(straingroupOrder, naClasses))]
  
  # breaks if there are no empty bins  
  if(nrow(plotDT[is.na(strain)]) == 0){
    spDT.plot <- spDT.plot[!is.na(strain)]
  }
  
  pointshapeGroups <- c(dot=16, cross=4)
  plotDT[, pshape := "dot"]
  plotDT[grepl("\\(", strain), pshape := "cross"] # ye ... I know this is dirty ...
  
  outplot <- ggplot(plotDT, aes(plotstart, strain, color=strainGroup)) +
    geom_point(aes(size=value, shape=pshape), alpha=0.8)  + 
    geom_vline(xintercept = chromEnds, linetype=3, color="grey") +
    facet_grid(strainGroup ~ ., scales="free", space="free") +
    theme_grey() +
    scale_x_continuous(breaks = chromLabelPos, labels = names(chromLabelPos), expand = c(0.05,0.05)) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid")) +
    scale_shape_manual(values = pointshapeGroups, guide="none") + 
    scale_size_continuous("BES") + 
    scale_color_manual("Strain Group", values = strainGroupColors, na.value = "darkgrey", guide="none") +
    geom_text(data = spDT.plot, mapping = aes(x = Inf, y = strain, label = SPtext), hjust = 1.2)
  
  # annotate specific loci of interest (e.g. Pdx-cre on chr8 from FVB mice)
  if(!is.null(annotationLoci)){
    annotationLoci.gr <- makeGRangesFromDataFrame(annotationLoci)
    hits <- findOverlaps(annotationLoci.gr, GenomeBinsGR)
    
    markedNotFound <- annotationLoci[!queryHits(hits)]
    if(nrow(markedNotFound) > 0){
      warning(paste0("One or more annotation loci are not valid mm10 genomic locations: ", markedNotFound[, paste0(name, collapse = ",")]))
    }
    
    annotationLoci.plot <- cbind(annotationLoci[queryHits(hits)], GenomeBins[subjectHits(hits), .(plotstart, plotend)])
    annotationLoci.plot <- annotationLoci.plot[, .(plotstart=min(plotstart), plotend=max(plotend)), by=.(strainGroup, name, strain)]
    annotationLoci.plot[, center := plotstart + ((plotend-plotstart)/2)]
    
    # 1. convert "all" to all possible strain groups
    plotGroups <- plotDT[!is.na(strainGroup), paste0(unique(strainGroup), collapse = ",")]
    annotationLoci.plot[strainGroup == "all", strainGroup := plotGroups]
    annotationLoci.plot <- cSplit(annotationLoci.plot, "strainGroup", ",", direction = "long", type.convert = F)
    
    # 2. determine the max "height" the rect can go for each straingroup based on the number of identified substrains, i.e. the index
    groupN <- plotDT[!is.na(strainGroup), uniqueN(strain), by=strainGroup]
    setnames(groupN, "V1", "maxN")
    annotationLoci.plot <- merge(annotationLoci.plot, groupN)
    
    # 3. if strain is not "all", the rect should only be around this specific strain, so we need to find its index (=height boundaries)
    annotationLoci.plot[, minN := 1]
    setorder(plotDT, -strain)
    strainIndex <- unique(plotDT[!is.na(strainGroup), .(strain, strainGroup)])
    strainIndex[ , I := .N:1, by = strainGroup] # in decreasing order
    
    #strainIndex[, lineN := .N-1, by=gsub(" \\(.*\\)", "", strain)] # weird error
    strainIndex[, tmpstrain := gsub(" \\(.*\\)", "", strain)]
    strainIndex[, lineN := .N-1, by=tmpstrain]
    strainIndex[, tmpstrain := NULL]
    
    annotationLoci.plot <- merge(annotationLoci.plot, strainIndex[, .(strain, I, lineN)], by="strain", all.x=T)
    annotationLoci.plot[!is.na(I), minN := I-lineN]
    annotationLoci.plot[!is.na(I), maxN := I]
    annotationLoci.plot[, strainGroup := factor(strainGroup)] # without this, the ordering will be broken, no idea why
    
    # text offset (not sure how this will behave for very large plots with many many strains)
    text.offset <- 0.5
    outplot <- outplot + 
      geom_rect(data=annotationLoci.plot,aes(xmin=plotstart,xmax=plotend,ymin=minN-text.offset,ymax=maxN+(text.offset-0.15)),inherit.aes=FALSE, fill=NA, color="black") +
      geom_text(data=annotationLoci.plot, aes(x=center, y = maxN, label = name), color="black", size=3, vjust = -2)
  }
  
  if(returnData){
    return(list(plotBins=plotDT, plotStats=spDT.plot, plotObj=outplot))
  } else {
    return(list(plotObj=outplot))
  }
}

# this is a (very very quick and dirty) version of the plotting function, just for individual chromosomes returned as a list
PlotStrainBinsByChrom <- function(sampleResults, spDT, samplename, GenomeBins, strainGroupColors, annotationLoci=NULL){
  GenomeBinsGR <- makeGRangesFromDataFrame(GenomeBins)
  
  plotDT <- merge(GenomeBins, sampleResults[, .(UnfiedPosition, strain, N, value, sample, strainGroup)], all.x=T, by="UnfiedPosition")
  
  spDT.plot <- copy(spDT) # work copy
  # remove NA if the entire genome is covered by SNPs (this would break the % text)
  if(all(spDT.plot[is.na(strain), SP] %in% 0)){
    spDT.plot <- spDT.plot[!is.na(strain)]
  }
  
  # assign the max possible bubble size to the NA bubbles (if to prevent warnings for full NA samples)
  if(nrow(plotDT[!is.na(value)]) == 0){
    maxValue <- 1 # if all bins are NA
  } else {
    maxValue <- plotDT[, max(value, na.rm = TRUE)]
  }
  
  # legacy code snippets
  naClasses <- NA
  plotDT[is.na(strain), value := maxValue] 
  
  sub.spDT <- spDT.plot[!strain %in% naClasses]
  
  sub.spDT <- sub.spDT[!grepl("\\(", strain)]
  setorder(sub.spDT, SP)
  strainOrder <- as.character(sub.spDT$strain)
  
  todos.strains <- unique(spDT.plot[grepl("\\(", strain), gsub(" \\(.*\\)", "", strain)])
  for(strainname in todos.strains){
    addIndex <- which(strainOrder==strainname)-1
    strainOrder <- append(strainOrder, paste0(strainname, " (filtered)"), after=addIndex)
    strainOrder <- append(strainOrder, paste0(strainname, " (no coverage)"), after=addIndex)
  }
  
  plotDT[, strain := factor(strain, levels = c(strainOrder, naClasses))]
  spDT.plot[, strain := factor(strain, levels = c(strainOrder, naClasses))]
  
  setorder(sub.spDT, -SP)
  straingroupOrder <- unique(as.character(sub.spDT$strainGroup))
  plotDT[, strainGroup := factor(strainGroup, levels = c(straingroupOrder, naClasses))]
  spDT.plot[, strainGroup := factor(strainGroup, levels = c(straingroupOrder, naClasses))]
  
  # breaks if there are no empty bins  
  if(nrow(plotDT[is.na(strain)]) == 0){
    spDT.plot <- spDT.plot[!is.na(strain)]
  }
  
  pointshapeGroups <- c(dot=16, cross=4)
  plotDT[, pshape := "dot"]
  plotDT[grepl("\\(", strain), pshape := "cross"] # ye ... I know this is dirty ...
  
  outplot.chrlist <- list()
  for(mychrom in levels(plotDT[, unique(chr)])){
    chrom.length <- plotDT[chr == mychrom, max(end)]
    chrom.cuts <- seq(0,chrom.length, 10000000)
    names(chrom.cuts) <- chrom.cuts / 1000000
    
    outplot <- ggplot(plotDT[chr == mychrom], aes(start, strain, color=strainGroup)) +
      geom_point(aes(size=value, shape=pshape), alpha=0.8)  + 
      facet_grid(strainGroup ~ ., scales="free", space="free") +
      theme_grey() +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid")) +
      scale_shape_manual(values = pointshapeGroups, guide="none") + 
      scale_size_continuous("BES") + 
      scale_color_manual("Strain Group", values = strainGroupColors, na.value = "darkgrey", guide="none") +
      scale_x_continuous(breaks = chrom.cuts, labels = names(chrom.cuts), expand = c(0.05,0.05)) +
      geom_vline(xintercept = chrom.cuts, linetype=3, color="grey") +
      xlab(paste0("chr", mychrom))
    
    # annotate specific loci of interest (e.g. Pdx-cre on chr8 from FVB mice)
    if(!is.null(annotationLoci)){

      # for the corresponding bin
      annotationLoci.gr <- makeGRangesFromDataFrame(annotationLoci)
      hits <- findOverlaps(annotationLoci.gr, GenomeBinsGR)
      
      markedNotFound <- annotationLoci[!queryHits(hits)]
      if(nrow(markedNotFound) > 0){
        warning(paste0("One or more annotation loci are not valid mm10 genomic locations: ", markedNotFound[, paste0(name, collapse = ",")]))
      }
      
      annotationLoci.plot <- cbind(annotationLoci[queryHits(hits)], GenomeBins[subjectHits(hits), .(plotstart=start, plotend=end)])
      annotationLoci.plot <- annotationLoci.plot[, .(chr=chr, plotstart=min(plotstart), plotend=max(plotend)), by=.(strainGroup, name, strain)]
      annotationLoci.plot[, center := plotstart + ((plotend-plotstart)/2)]
      
      # reduce to relevant chromosome
      annotationLoci.plot <- annotationLoci.plot[chr == mychrom]
      
      if(nrow(annotationLoci.plot) > 0){
        # addition: only for specific chrom
        keepstrain <- plotDT[chr==mychrom, unique(strain)]
        
        # 1. convert "all" to all possible strain groups
        plotGroups <- plotDT[!is.na(strainGroup), paste0(unique(strainGroup), collapse = ",")]
        annotationLoci.plot[strainGroup == "all", strainGroup := plotGroups]
        annotationLoci.plot <- cSplit(annotationLoci.plot, "strainGroup", ",", direction = "long", type.convert = F)
        
        # 2. determine the max "height" the rect can go for each straingroup based on the number of identified substrains, i.e. the index
        groupN <- plotDT[strain %in% keepstrain, uniqueN(strain), by=strainGroup]
        setnames(groupN, "V1", "maxN")
        annotationLoci.plot <- merge(annotationLoci.plot, groupN)
        
        # 3. if strain is not "all", the rect should only be around this specific strain, so we need to find its index (=height boundaries)
        annotationLoci.plot[, minN := 1]
        setorder(plotDT, -strain)
        strainIndex <- unique(plotDT[!is.na(strainGroup), .(strain, strainGroup)])
        strainIndex <- strainIndex[strain %in% keepstrain]
        strainIndex[ , I := .N:1, by = strainGroup] # in decreasing order
        strainIndex[, lineN := .N-1, by=gsub(" \\(.*\\)", "", strain)]
        annotationLoci.plot <- merge(annotationLoci.plot, strainIndex[, .(strain, I, lineN)], by="strain", all.x=T)
        annotationLoci.plot[!is.na(I), minN := I-lineN]
        annotationLoci.plot[!is.na(I), maxN := I]
        annotationLoci.plot[, strainGroup := factor(strainGroup)] # without this, the ordering will be broken, no idea why
        
        # text offset (not sure how this will behave for very large plots with many many strains)
        text.offset <- 0.5
        
        outplot <- outplot + 
          geom_rect(data=annotationLoci.plot,aes(xmin=plotstart,xmax=plotend,ymin=minN-text.offset,ymax=maxN+(text.offset-0.15)),inherit.aes=FALSE, fill=NA, color="black") +
          geom_text(data=annotationLoci.plot, aes(x=center, y = maxN, label = name), color="black", size=3, vjust = -2)
      }
    }
    outplot.chrlist[[mychrom]] <- outplot
  }
  return(outplot.chrlist)
}

# for each mutation, count the number of neighboring mutations (+-5bp) to find clustered events
# this has to be tested, but for some examples this clearly reduced the noise (e.g. JAX-FVB 10% 129 was removed)
BlacklistClusteredEvents <- function(sampleSNPs, strainDT, window=5, maxNeighbors=1, preserveStrains=F){
  blacklistedSNPs <- copy(sampleSNPs)
  blacklistedSNPs[, Chr := as.character(Chr)]
  blacklistedSNPs <- merge(blacklistedSNPs, strainDT[, .(Chr, Pos, Ref, Alt, strainGroup, strain)], by=c("Chr", "Pos", "Ref", "Alt"), all.x=T)
  
  bl.tmp <- unique(blacklistedSNPs[, .(Chr, start=Pos, end=Pos)])
  gr <- makeGRangesFromDataFrame(bl.tmp)
  gr.extended <- makeGRangesFromDataFrame(bl.tmp[, .(Chr, start=start-window, end=end+window)])
  
  counts <- countOverlaps(gr, gr.extended)-1 # -1 to remove "self-hit"
  bl.tmp[, nNeighbors := counts]
  
  # optional, remove those neighbors from the counting if they are form the same strand group. this can be better or not and has to be tested further
  if(preserveStrains){
    bl.tmp2 <- unique(blacklistedSNPs[, .(Chr=paste0(Chr, "-", strainGroup), start=Pos, end=Pos)])
    gr2 <- makeGRangesFromDataFrame(bl.tmp2)
    gr2.extended <- makeGRangesFromDataFrame(bl.tmp2[, .(Chr, start=start-window, end=end+window)])
    counts2 <- countOverlaps(gr2, gr2.extended)-1 # -1 to remove "self-hit"
    bl.tmp2[, counts2 := counts2]
    bl.tmp2 <- bl.tmp2[!grepl("-NA", Chr)]
    bl.tmp2[, Chr := gsub("-.*", "", Chr)]
    bl.tmp <- merge(bl.tmp, bl.tmp2, by=c("Chr", "start", "end"), all.x=T)
    bl.tmp[, nNeighbors := nNeighbors - counts2]
    bl.tmp[, counts2 := NULL]
  }
  
  blacklistedSNPs <- merge(blacklistedSNPs, bl.tmp, by.x=c("Chr", "Pos"),by.y=c("Chr", "start"))
  blacklistedSNPs <- unique(blacklistedSNPs[, sum(nNeighbors), by=.(Chr, Pos)])
  setnames(blacklistedSNPs, "V1", "nNeighbors")
  blacklistedSNPs <- merge(strainDT, blacklistedSNPs, by=c("Chr", "Pos"))
  blacklistedIDs <- blacklistedSNPs[nNeighbors >= maxNeighbors, ID]
  
  return(blacklistedIDs)
}

# now we can classify each position into:
# - B6 SNP = no strain specific SNP found for this position
# - not B6 SNP = another strain SNP was found for this position (the strain SNP is unspecific, since all strains have it except B6)
# - no coverage = position is not covered or the SNPs did not pass the MAPQ and min-reads cutoffs (or the bin does not have any overlapping B6 SNPs even with coverage)
DetermineBlack6 <- function(B6DT, b6Counts, GenomeBinsB6, minReadsPerSNP=1){
  
  GenomeBinsB6.gr <- makeGRangesFromDataFrame(GenomeBinsB6)
  
  # add additional information to the pool of B6 candidate SNPs
  B6DT.sample.snps <- merge(B6DT, b6Counts[TotalCount > minReadsPerSNP, .(Chr=as.character(Chr), Pos, Cov=TotalCount, isB6)], by=c("Chr", "Pos"), all.x=T) 
  
  B6DT.sample.snps[is.na(Cov), Cov := 0]
  B6DT.sample.snps[is.na(isB6), isB6 := "no"] # SNPs without read coverage
  
  B6DT.sample.snps.gr <- makeGRangesFromDataFrame(B6DT.sample.snps[, .(chr=Chr, start=Pos, end=Pos)])
  hits <- findOverlaps(GenomeBinsB6.gr, B6DT.sample.snps.gr)
  B6DT.sample.snps <- cbind(GenomeBinsB6[queryHits(hits), -c("plotstart", "plotend")], B6DT.sample.snps[subjectHits(hits), .(Pos, Ref, Alt, ID, strainGroup, strain, isB6, Cov)])
  
  # add bins without any Black6 SNPs (those bins will be added to "not covered")
  missingBins <- GenomeBinsB6[!UnfiedPosition %in% B6DT.sample.snps$UnfiedPosition, ]
  missingBins[, c("Pos", "Ref", "Alt", "ID", "strainGroup", "strain", "isB6", "Cov") := NA]
  missingBins <- missingBins[, -c("plotstart", "plotend")]
  missingBins[, Cov := 0]
  B6DT.sample.snps <- rbind(B6DT.sample.snps, missingBins)
  
  # for now a naive approach: more yes or no will define the bin status
  B6DT.sample.snps[Cov > 0, isB6 := ifelse(sum(isB6 == "no") > sum(isB6 == "yes"), "no", "yes"), by=UnfiedPosition]
  
  # can be either Cov=0, Cov!=0 but not B6, or real B6
  notB6 <- B6DT.sample.snps[, all(isB6 == "no"), by=UnfiedPosition]
  notB6 <- notB6[V1 == "TRUE", UnfiedPosition]
  B6DT.sample.snps[UnfiedPosition %in% notB6, strain := "not C57BL_6J"]
  
  notCov <- B6DT.sample.snps[, sum(Cov) == 0, by=UnfiedPosition]
  notCov <- notCov[V1 == "TRUE", UnfiedPosition]
  B6DT.sample.snps[UnfiedPosition %in% notCov, strain := "no coverage"]
  
  B6DT.sample.bins <- B6DT.sample.snps[, .N, by=.(UnfiedPosition, strain)]
  
  B6DT.sample.bins <- merge(B6DT.sample.bins, GenomeBinsB6[, .(UnfiedPosition, plotstart, plotend)], by="UnfiedPosition", all.x=T)
  B6DT.sample.bins[, value := NA]
  B6DT.sample.bins[, strainGroup := "custom"]
  
  spDT.plot.add <- data.table(B6DT.sample.bins[, table(strain) / nrow(GenomeBinsB6) * 100])
  spDT.plot.add <- spDT.plot.add[, .(strain, sample=as.character(NA), SP=round(N), SPtext=paste0(round(N), "%"), strainGroup="custom")]
  
  # add the % also regarding covered regions (i.e. black6 vs NOT black6)
  spDT.plot.add[strain != "no coverage", SPcovered := round(SP / sum(SP) * 100)]
  spDT.plot.add[strain != "no coverage", SPtext := paste0(SPtext, "\n  (",SPcovered,"%)")]
  spDT.plot.add[, SPcovered := NULL]
  
  B6DT.sample.bins[, strain := factor(strain, levels=c("no coverage", "not C57BL_6J", "C57BL_6J"))]
  
  return(list(B6DT.sample.bins=B6DT.sample.bins, spDT.plot.add=spDT.plot.add))
}

PlotBlack6 <- function(B6.res, GenomeBins, strainGroupColors, samplename=NULL, sampleResults=NULL, spDT=NULL, mode="withmain", annotationLoci=NULL, returnData=F){

  # some data for plotting
  chromEnds <- GenomeBins[, max(plotend), by=chr]$V1 # equal to the original "Len" variable
  chromEnds <- c(0, chromEnds)
  chromnames <- GenomeBins[, paste0("chr", unique(chr))]
  chromLabelPos <- diff(chromEnds)/2
  names(chromLabelPos) <- chromnames
  chromLabelPos <- chromEnds+c(chromLabelPos, 0) # with pseuso 0 after last chromosome
  GenomeBinsGR <- makeGRangesFromDataFrame(GenomeBins)
  
  if(mode == "withmain" & (is.null(sampleResults) | is.null(spDT) | is.null(samplename))){
    stop("To use withmain, the main function has to be executed first and the 3 objects needs to be given to the function")
  }
  
  # plot as standalone
  if(mode == "standalone"){
    plotDT2 <- B6.res$B6DT.sample.bins
    plotDT2[is.na(value), value := 1]
    plotDT2[, value := as.numeric(value)]
    spDT.plot2 <- B6.res$spDT.plot.add
    plotDT2[, pshape := "cross"]
  } else if(mode == "withmain"){
    # or combined with main plot (add below)
    plotList <- PlotStrainBins(sampleResults, spDT, samplename, GenomeBins, strainGroupColors, annotationLoci = annotationLoci, returnData = T, excludeNA=T)
    plotDT2 <- rbind(plotList$plotBins, B6.res$B6DT.sample.bins, fill=T)
    
    if(all(is.na(plotDT2$value))){
      maxVal <- 1
    } else {
      maxVal <- plotDT2[, max(value, na.rm = T)]
    }
    
    plotDT2[is.na(value), value := maxVal]
    plotDT2[is.na(pshape), pshape := "cross"]
    spDT.plot2 <- rbind(plotList$plotStats, B6.res$spDT.plot.add)
    plotDT2[is.na(strainGroup), strainGroup := "NA"]
    spDT.plot2[is.na(strainGroup), strainGroup := "NA"]
  }
  
  pointshapeGroups <- c(dot=16, cross=4)
  
  outplot <- ggplot(plotDT2, aes(plotstart, strain, color=strainGroup)) +
    geom_point(aes(size=value, shape=pshape), alpha=0.8)  + 
    geom_vline(xintercept = chromEnds, linetype=3, color="grey") +
    facet_grid(strainGroup ~ ., scales="free", space="free") +
    theme_grey() +
    scale_x_continuous(breaks = chromLabelPos, labels = names(chromLabelPos), expand = c(0.05,0.05)) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid")) +
    scale_shape_manual(values = pointshapeGroups) + 
    scale_color_manual("Strain Group", values = strainGroupColors, na.value = "darkgrey") +
    geom_text(data = spDT.plot2, mapping = aes(x = Inf, y = strain, label = SPtext), hjust = 1.2)

  # this is copy-paste from the main function PlotStrainBins
  # annotate specific loci of interest (e.g. Pdx-cre on chr8 from FVB mice)
  if(!is.null(annotationLoci)){
    annotationLoci.gr <- makeGRangesFromDataFrame(annotationLoci)
    hits <- findOverlaps(annotationLoci.gr, GenomeBinsGR)
    
    markedNotFound <- annotationLoci[!queryHits(hits)]
    if(nrow(markedNotFound) > 0){
      warning(paste0("One or more annotation loci are not valid mm10 genomic locations: ", markedNotFound[, paste0(name, collapse = ",")]))
    }
    
    annotationLoci.plot <- cbind(annotationLoci[queryHits(hits)], GenomeBins[subjectHits(hits), .(plotstart, plotend)])
    annotationLoci.plot <- annotationLoci.plot[, .(plotstart=min(plotstart), plotend=max(plotend)), by=.(strainGroup, name, strain)]
    annotationLoci.plot[, center := plotstart + ((plotend-plotstart)/2)]
    
    # 1. convert "all" to all possible strain groups
    plotGroups <- plotDT2[!is.na(strainGroup) & strainGroup != "custom", paste0(unique(strainGroup), collapse = ",")]
    annotationLoci.plot[strainGroup == "all", strainGroup := plotGroups]
    annotationLoci.plot <- cSplit(annotationLoci.plot, "strainGroup", ",", direction = "long", type.convert = F)
    
    # 2. determine the max "height" the rect can go for each straingroup based on the number of identified substrains, i.e. the index
    groupN <- plotDT2[!is.na(strainGroup) & strainGroup != "custom", uniqueN(strain), by=strainGroup]
    setnames(groupN, "V1", "maxN")
    annotationLoci.plot <- merge(annotationLoci.plot, groupN)
    
    # 3. if strain is not "all", the rect should only be around this specific strain, so we need to find the its index (=height boundaries)
    annotationLoci.plot[, minN := 1]
    setorder(plotDT2, -strain)
    strainIndex <- unique(plotDT2[!is.na(strainGroup) & strainGroup != "custom", .(strain, strainGroup)])
    strainIndex[ , I := .N:1, by = strainGroup] # in decreasing order
    strainIndex[, nLines := .N, by=gsub(" \\(.*\\)", "", strain)] # to account for added X lines
    annotationLoci.plot <- merge(annotationLoci.plot, strainIndex[, .(strain, I, nLines)], by="strain", all.x=T)
    annotationLoci.plot[!is.na(I), minN := I]
    annotationLoci.plot[!is.na(I), maxN := I]
    annotationLoci.plot[, strainGroup := factor(strainGroup)] # without this, the ordering will be broken, no idea why
    
    # if X for no-cov or filtered were added, we need to adjust this value
    annotationLoci.plot[, minN := minN-(nLines-1)]
    
    # text offset (not sure how this will behave for very large plots with many many strains)
    text.offset <- 0.5
    
    outplot <- outplot + 
      geom_rect(data=annotationLoci.plot,aes(xmin=plotstart,xmax=plotend,ymin=minN-text.offset,ymax=maxN+(text.offset-0.15)),inherit.aes=FALSE, fill=NA, color="black") +
      geom_text(data=annotationLoci.plot, aes(x=center, y = maxN, label = name), color="black", size=3, vjust = -2)
  }

  if(returnData){
    return(list(plotBins=plotDT2, plotStats=spDT.plot2, plotObj=outplot))
  } else {
    return(list(plotObj=outplot))
  }
}

# all BRC commands together
PrintCommandsBRC <- function(strainResources, mhcResources, samplename, bamfile, workfolder, reffile, append=T){
  PrintCommandsBRC.STRAIN(strainResources, samplename, bamfile, workfolder, reffile, append=append)
  PrintCommandsBRC.MHC(mhcResources, samplename, bamfile, workfolder, reffile, append=T)
  PrintCommandsBRC.MHC.B6(mhcResources, samplename, bamfile, workfolder, reffile, append=T)
}

# outputs BASH script: get the coverage and the allele information for chosen SNPs (once for B6 and once for all strain-specific ones)
# this requires Docker!
# this will append by default, so one can do this once for all samples
# # EXAMPLE (on WS4)
# workfolder <- "/home/rad/Downloads/" # /home/rad/media/rad/HDD2/temp_niklas/mousestrains_bam_read_count
# bamfile=paste0("/mnt/WS6/media/rad/HDD2/MoCaSeq_runs/JAX_mouseStrains/",samplename,"/results/bam/",samplename,".Normal.bam")
# bamfile=paste0("/mnt/WS6/media/rad/HDD2/MoCaSeq_runs/JAX_mouseStrains/",samplename,"/results/Mutect2/",samplename,".Normal.m2.bam")
# reffile="/media/rad/SSD1/MoCaSeq_ref/GRCm38.p6/GRCm38.p6.fna"
# PrintCommandsBRC(samplename, bamfile, workfolder, reffile, append=F)
PrintCommandsBRC.STRAIN <- function(strainResources, samplename, bamfile, workfolder, reffile, append=T){
  # VERSION: bam-readcount version: 0.8.0-unstable-50-f5f5ed0 (commit f5f5ed0)
  # this will write 2 commands into the file, once line for black6 and one for the other strains (it is not slower and cleaner this way)
  
  # load files from resource data
  B6DT <- strainResources$B6STRAIN
  strainDT <- strainResources$variantsStrain
  
  bamPath <- dirname(bamfile)
  bamFile <- basename(bamfile)
  
  refPath <- dirname(reffile)
  refFile <- basename(reffile)
  
  # FOR BLACK6
  outfile <- paste0(workfolder, "/", samplename, ".b6.brc.tsv.gz")
  
  # write the lookup bed to the workdir
  b6BED <- unique(B6DT[, .(Chr, Pos, Pos)])
  bedfile <- paste0(workfolder, "/b6_positions.bed")
  if(!file.exists(bedfile)){  fwrite(b6BED, bedfile, sep="\t", col.names = F)}
  bedfile <- basename(bedfile)
  
  # print commands to file
  command <- paste0("time docker run --user $(id -u):$(id -g) -v ",
                    refPath,":/reference/ -v ",bamPath,":/bam/ -v ",workfolder,":/workdir/ mgibio/bam-readcount -f /reference/",refFile,
                    " /bam/",bamFile," -l /workdir/",bedfile," -w 0 | cut -f1,2,3,4,6,7,8,9 | gzip > ", outfile)
  
  commandFile <- paste0(workfolder, "/get_read_counts_from_bam.sh")
  
  echostring <- paste0("echo \"Running bam-readcount (black6) for: \"", samplename)
  command <- paste0(echostring, " & ", command)
  cat("\n",file=commandFile,append=append)
  cat(command,file=commandFile,append=TRUE)
  
  
  # FOR ALL OTHER STRAINS
  outfile <- paste0(workfolder, "/", samplename, ".strains.brc.tsv.gz")
  
  # write the lookup bed to the workdir
  strainBED <- unique(strainDT[, .(Chr, Pos, Pos)])
  bedfile <- paste0(workfolder, "/strain_positions.bed")
  if(!file.exists(bedfile)){fwrite(strainBED, bedfile, sep="\t", col.names = F)}
  bedfile <- basename(bedfile)
  
  # print commands to file
  command <- paste0("time docker run --user $(id -u):$(id -g) -v ",
                    refPath,":/reference/ -v ",bamPath,":/bam/ -v ",workfolder,":/workdir/ mgibio/bam-readcount -f /reference/",refFile,
                    " /bam/",bamFile," -l /workdir/",bedfile," -w 0 | cut -f1,2,3,4,6,7,8,9 | gzip > ", outfile)
  
  commandFile <- paste0(workfolder, "/get_read_counts_from_bam.sh")

  echostring <- paste0("echo \"Running bam-readcount (strain specific) for: \"", samplename)
  command <- paste0(echostring, " & ", command)
  cat("\n",file=commandFile,append=TRUE)
  cat(command,file=commandFile,append=TRUE)
}

# BRC (bam read count from mgibio/bam-readcount)
ProcessBRC <- function(brc.file, minMQ=60, minReads=1, minBQ=20, minAF=0.2, isB6=F, isHLA=F, B6data=B6DT){
  
  # B6data: B6DT for strain, B6HLA for MHC
  
  #check if file is empty
  if(file.size(brc.file) <= 20 | is.na(file.size(brc.file))){
    print(paste0("BRC file is empty: ", brc.file))
    return(NULL)
  }
  
  brcCounts <- fread(brc.file, sep="\t")
  
  brcCounts[, V3 := toupper(V3)]
  
  # doing columns to rows and then working on it is actually slower, so we do it like this
  
  # capture the nucleotide information + mapping and base quality
  keepIndices <- c(2:4)
  brcCounts[, c("A", "A.MAPQ", "A.BASEQ") := tstrsplit(V5, ":", fixed=TRUE, keep=keepIndices)]
  brcCounts[, c("C", "C.MAPQ", "C.BASEQ") := tstrsplit(V6, ":", fixed=TRUE, keep=keepIndices)]
  brcCounts[, c("G", "G.MAPQ", "G.BASEQ") := tstrsplit(V7, ":", fixed=TRUE, keep=keepIndices)]
  brcCounts[, c("T", "T.MAPQ", "T.BASEQ") := tstrsplit(V8, ":", fixed=TRUE, keep=keepIndices)]
  
  brcCounts <- brcCounts[, -c(paste0("V", 5:8))] # remove unwanted columns
  brcCounts <- unique(brcCounts)
  
  setnames(brcCounts, c(paste0("V", 1:4)), c("Chr", "Pos", "Ref", "TotalCount"))
  setcolorder(brcCounts, c("Chr", "Pos", "Ref", "TotalCount", "A", "C", "T", "G"))
  
  brcCounts <- data.table(gather(brcCounts, Alt, AltCount, A:T))
  brcCounts[, AltCount := as.numeric(AltCount)]
  brcCounts[, TotalCount := as.numeric(TotalCount)]
  
  brcCounts <- brcCounts[AltCount >= minReads] # additional depth filter
  
  brcCounts[, AF := AltCount / TotalCount]
  brcCounts <- brcCounts[AF >= minAF] # additional frequency filter (strain SNP should be hom)
  
  brcCounts[, ID := paste(Chr, Pos, Ref, Alt, sep="-")]
  
  # summarize the MAPQ by nucleotide
  brcCounts[Alt == "A", MAPQ := as.numeric(A.MAPQ)]
  brcCounts[Alt == "C", MAPQ := as.numeric(C.MAPQ)]
  brcCounts[Alt == "G", MAPQ := as.numeric(G.MAPQ)]
  brcCounts[Alt == "T", MAPQ := as.numeric(T.MAPQ)]
  
  # summarize the BASEQ by nucleotide
  brcCounts[Alt == "A", BASEQ := as.numeric(A.BASEQ)]
  brcCounts[Alt == "C", BASEQ := as.numeric(C.BASEQ)]
  brcCounts[Alt == "G", BASEQ := as.numeric(G.BASEQ)]
  brcCounts[Alt == "T", BASEQ := as.numeric(T.BASEQ)]
  
  # filter by mapping quality
  brcCounts <- brcCounts[MAPQ >= minMQ]
  brcCounts <- brcCounts[BASEQ >= minBQ]
  
  # remove unwanted columns (we do not select the wanted because we probably need to filter by additional columns later at some point)
  brcCounts <- brcCounts[, -grepl("\\.", colnames(brcCounts)), with=F]
  
  setcolorder(brcCounts, c("ID", "Chr", "Pos", "Ref", "Alt", "MAPQ", "TotalCount", "AltCount"))
  
  brcCounts <- unique(brcCounts)
  
  # this is needed to process the black6 "anti-SNPs"
  if(isB6){
    brcCounts[, nSNPs := .N, by=.(Chr, Pos)] # these are homozygous positions
    brcCounts[nSNPs == 1 & Ref==Alt, isB6 := "yes"]
    brcCounts[is.na(isB6), isB6 := "no"]
    
    # now find positions which seem to be "not-B6" because there is some mutation/SNP --> check if this SNP is indeed a strain specific SNP else keep
    # if nSNPs != 1, check if the SNP is a strain signature SNP (keep if it is not, which means it is a mutation or technical artifact)
    rescueHet <- brcCounts[nSNPs != 1]
    rescueHet[, ID2 := paste(Chr, Pos, Ref, Alt, sep="-")]
    rescueHet <- rescueHet[!ID2 %in% B6data$ID2]
    rescueHet[, N:=.N, by=.(Chr, Pos)] # all with N = 1 are still invalid, since there is a strain specific mutation (the missing line will be isB6=no in the other set)
    rescueHet <- rescueHet[N != 1]
    brcCounts[ID %in% rescueHet$ID, isB6 := "yes"]
  }
  
  # this is needed to process the HLA black6 "anti-SNPs"
  if(isHLA){
    brcCounts[, nSNPs := .N, by=.(Chr, Pos)] # these are homozygous positions
    brcCounts[nSNPs == 1 & Ref==Alt, isB6 := "yes"]
    brcCounts[is.na(isB6), isB6 := "no"]
    
    # now find positions which seem to be "not-B6" because there is some mutation/SNP --> check if this SNP is indeed a strain specific SNP else keep
    # if nSNPs != 1, check if the SNP is a strain signature SNP (keep if it is not, which means it is a mutation or technical artifact)
    rescueHet <- brcCounts[nSNPs != 1]
    rescueHet[, ID2 := paste(Chr, Pos, Ref, Alt, sep="-")]
    rescueHet <- rescueHet[!ID2 %in% B6data$ID2]
    rescueHet[, N:=.N, by=.(Chr, Pos)] # all with N = 1 are still invalid, since there is a strain specific mutation (the missing line will be isB6=no in the other set)
    rescueHet <- rescueHet[N != 1]
    brcCounts[ID %in% rescueHet$ID, isB6 := "yes"]
  }
  
  return(brcCounts)
}

# this will add additional annotations to the bubble plot: 2 new lines with X for the filtered and not-covered bins (strains are select by percent with addCovMinGES)
# this requires additional pre-calculation by the user (mgibio/bam-readcount) for the coverage at SNP positions
# addCovMinGES: for all strains with a percentage coverage (number of the right side of the plot) >= "value", the additional coverage/filter with X will be shown (if addCovX=T)
BubblesAddCovX <- function(brcFile, sampleResults, spDT, GenomeBins, strainDT, CovFile, addCovMinGES, strainGroupAssignment){
  GenomeBinsGR <- makeGRangesFromDataFrame(GenomeBins)
  
  # Process the BRC (bam read count) file
  BRC.unfiltered <- ProcessBRC(brcFile, minMQ=0, minReads=0, minBQ = 0, minAF=0) # this will be used for the X in the plot
  BRC.unfiltered <- merge(strainDT, BRC.unfiltered[, .(ID, AltCount, TotalCount)], all.x=T)
  BRC.unfiltered[is.na(AltCount), AltCount := 0]
  BRC.unfiltered[is.na(TotalCount), TotalCount := 0]
  
  # generate the bins for "is covered or not"
  BRC.unfiltered.GR <- makeGRangesFromDataFrame(BRC.unfiltered[, .(chr=Chr, start=Pos, end=Pos)])
  hits <- findOverlaps(GenomeBinsGR, BRC.unfiltered.GR)
  countBins <- cbind(GenomeBins[queryHits(hits)], BRC.unfiltered[subjectHits(hits), .(strain, AltCount, TotalCount)])
  countBins <- countBins[, .(.N, TotalCounts=sum(TotalCount), AltCounts=sum(AltCount)), by=.(UnfiedPosition, strain)]
  
  # there can be missings bins (e.g. small end of chromosomes), we re-add them with 0 counts
  initBins <- strainGroupAssignment[, .(UnfiedPosition=rep(GenomeBins$UnfiedPosition)), by=strain]
  countBins <- merge(countBins, initBins, by=c("UnfiedPosition", "strain"), all.y=T)
  countBins[is.na(countBins)] <- 0 # replace all NAs with 0
  countBins <- merge(countBins, strainGroupAssignment, by="strain")
  
  # for each strain with enough strain percent (GES) add the missing bins with annotations 
  foundStrains <- spDT[SP >= addCovMinGES, strain]
  foundStrains <- foundStrains[!is.na(foundStrains)]
  for(astrain in foundStrains){
    missingBins <- GenomeBins[!UnfiedPosition %in% sampleResults[strain == astrain, UnfiedPosition], UnfiedPosition]
    
    addLines <- countBins[strain == astrain & UnfiedPosition %in% missingBins, .(strain=astrain, UnfiedPosition, N, nTotal=as.numeric(NA), value=as.numeric(NA), sample="FVB_NJ-1", 
                                                                                 BinScoring="all", strainGroup, groupname, plotcolor, groupplotcolor, TotalCounts, AltCounts)]
    # what is defined as "not covered"?
    addLinesNoCov <- addLines[TotalCounts == 0 & AltCounts == 0]
    addLinesNoCov[, strain := paste0(strain, " (no coverage)")]
    addLinesNoCov.sp <- round(nrow(addLinesNoCov) / nrow(GenomeBins) * 100) # this will be added to the summary at the side
    addLinesNoCov.sp.lines <- spDT[strain == astrain]
    addLinesNoCov.sp.lines[, strain := paste0(strain, " (no coverage)")]
    addLinesNoCov.sp.lines[, SP := addLinesNoCov.sp]
    
    # the others are covered but not used, i.e. filtered due to things like MAPQ, BASQ, AF, min coverage, other strain ...
    addLinesFiltered <- addLines[!UnfiedPosition %in% addLinesNoCov$UnfiedPosition]
    addLinesFiltered[, strain := paste0(strain, " (filtered)")]
    addLinesFiltered.sp <- round(nrow(addLinesFiltered) / nrow(GenomeBins) * 100) # this will be added to the summary at the side
    addLinesFiltered.sp.lines <- spDT[strain == astrain]
    addLinesFiltered.sp.lines[, strain := paste0(strain, " (filtered)")]
    addLinesFiltered.sp.lines[, SP := addLinesFiltered.sp]
    
    # combine and add
    addLines <- rbind(addLinesNoCov, addLinesFiltered)
    sampleResults <- rbind(sampleResults, addLines, fill=T)
    maxVal <- sampleResults[, mean(value, na.rm = T)]
    sampleResults[is.na(value), value := maxVal]
    
    # also add to sp annotation
    addLines2 <- rbind(addLinesNoCov.sp.lines, addLinesFiltered.sp.lines)
    
    # due to rounding errors (3x), this will will occasionally not be 100% --> we add those numbers to filtered, because they should be there anyways
    totalGES <- addLines2[, sum(SP)] + spDT[strain == astrain, SP]
    
    if(totalGES != 100){
      
      missingP <- 100-totalGES
      addLines2[grepl("filtered", strain), SP := SP + missingP]
      
      # in case filtered not found (this should never happen)
      if(nrow(addLines2[grepl("filtered", strain)]) == 0){addLines2[grepl("no coverage", strain), SP := SP + missingP]}
    }
    
    addLines2[, SPtext := paste0(SP, "%")]
    spDT <- rbind(spDT, addLines2)
  }
  return(list(sampleResults=sampleResults, spDT=spDT))
}



# STRAIN - MAIN FUNCTION ----------

# PARAMETERS EXPLAINED
# minCoverage = minimum number of reads for each SNP to be included (>=X, from raw mutation calling data)
# minAF = minimum allele frequency for each SNP to be included (>=X, from raw mutation calling data)
# minMQ = minimum mean mapping quality of reads for each SNP position to be included (>=X, from raw mutation calling data)
# minBQ = minimum base quality of reads for each SNP position to be included (>=X, from raw mutation calling data) (! only used for BRC files)
# removeClusteredEvents = remove SNPs with any neighboring (+-5bp) other SNPs
# maxNeighbors = number of allowed neighboring SNPs (>=X)
# minGES = the minimum percentage of genome explained by a strain to be shown (>X, e.g. 1 bin would be 0.004% explained for 10Mb bins, which is removed)
# --> decreasing this to -1 (so no bin filtering) can help to identify small transgenic loci (i.e. neighboring bins of the same strand) 
#     but will also increase the number of noise (i.e. individual strain bins scattered across the genome)
# minBinSNPS = the minimum number of SNPs per bin to assign a strain to that bin (>=X)
# minBES = the minimum enrichment score for each bin to be included (>X)
# keepHits = for each bin, choose whether to keep all, only the best or best+second best strain [DEFAULT = all]
# annotationLoci = data.table with information about which regions to highlight and annotate (can be specific for a strain but also all strains, see example)
# returnData = returns a list with different objects, including information about the bins and SNPs
# inputType = input data can be BRC, MoCaSeq, data
# selectFiles = TumorNormal, Tumor or Normal
# combineTN = combine Tumor+Normal by all hits or just the shared ones found in both

# NA-categories
# lowCov: no valid reads for all SNPs in a bin (relates to minCoverage, minAF and minMQ; not valid means either 0 reads or reads below those 3 cutoffs)
# lowSNPs: valid SNPs but the number of SNPs per bin is too low (relates to the minBinSNPs cutoff)
# lowES: SNPs are found, but the bin or the genomic enrichment score are below the chosen cutoffs (this relates to minBES and minGES, reduce these cutoffs to include lowES-bins into the analysis)

CheckParameters <- function(samplePath, inputType, seqType, selectFiles, combineTN, removeClusteredEvents, minBinSNPS, keepHits, minBES, 
                            minGES, minCoverage, minAF, minMQ, annotationLoci, returnData, keepEmptyStrains){
  
  if(!keepHits %in% c("best", "second", "both", "all")){
    stop(paste0("Invalid parameter used for keepHits (\"",keepHits,"\"). Please use one of the following: ", paste0(c("best", "second", "both", "all"), collapse = ", ")))
  }
  
  if(!inputType %in% c("MoCaSeq", "BRC", "data")){
    stop(paste0("Invalid parameter used for inputType (\"",inputType,"\"). Please use one of the following: ", paste0(c("MoCaSeq", "BRC", "data"), collapse = ", ")))
  }
  
}

RunAnalysis <- function(strainResources, samplePath, inputType="MoCaSeq", seqType="WES", selectFiles="TumorNormal", combineTN="all", removeClusteredEvents=T, maxNeighbors=1,
                        brcFile=NULL, b6brcFile=NULL, addCovMinGES=100, minBinSNPS=1, keepHits="all", minBES=NA, minGES=1, minCoverage=1, minAF=0.2, minMQ=58, minBQ=20,
                        annotationLoci=NULL, returnData=F, keepEmptyStrains=F, plotByChrom=F){

  
  # load the resource data from lookup object
  GenomeBins <- strainResources$GenomeBins
  genomeSNPs <- strainResources$genomeSNPs
  strainDT <- strainResources$variantsStrain
  strainGroupAssignment <- strainResources$strainGroupAssignment
  strainGroupColors <- strainResources$strainGroupColors
  nBins <- GenomeBins[, length(unique(UnfiedPosition))]
  B6DT <- strainResources$B6STRAIN
  GenomeBinsB6 <- strainResources$GenomeBinsB6
  
  # check if parameters are valid
  CheckParameters(samplePath, inputType, seqType, selectFiles, combineTN, removeClusteredEvents, minBinSNPS, keepHits, minBES, 
                  minGES, minCoverage, minAF, minMQ, annotationLoci, returnData, keepEmptyStrains)
  
  # some preset cutoffs for each sequencing type
  if(seqType=="WES" & is.na(minBES)){
    minBES <- 0.01 # maybe this cutoff --> see JAX-FVB2
  } else if(seqType=="WGS" & is.na(minBES)){
    minBES <- 0.1
  } else if(seqType=="lcWGS" & is.na(minBES)){
    minBES <- 0.01
  }
  
  print("Reading input data")
  strainBinsAllList <- MergeDataWithStrains(samplePath, strainDT, GenomeBins, genomeSNPs, returnData, inputType, selectFiles, 
                                            minCoverage, minAF, minMQ, minBQ, removeClusteredEvents, maxNeighbors, minBinSNPS, combineTN)
  
  # stop if no SNPs were found
  if(is.null(strainBinsAllList$strainBinsAll)){
    # but at least try to do the B6 plot is possible 
    if(!is.null(b6brcFile)){
      samplename <- strainBinsAllList$sampleName
      print("Checking black6 (C57BL_6J plot)")
      B6BRC <- ProcessBRC(b6brcFile, minMQ=58, minReads=2, minBQ = 20, minAF=0.3, isB6 = T, B6data = B6DT)
      
      if(is.null(B6BRC)){
        return(NULL)
      }
      
      # stop if no SNPs were found
      if(nrow(B6BRC) > 0){
        B6.res <- DetermineBlack6(B6DT, B6BRC, GenomeBinsB6)
        # by adding this, we can remove the original NA track
        output <- PlotBlack6(B6.res, GenomeBins, strainGroupColors, samplename, sampleResults=data.table(), spDT=data.table(), mode = "standalone", annotationLoci, returnData=returnData)
        return(output)
      }
    }
    return(NULL)
  }
  
  strainBinsAll <- strainBinsAllList$strainBinsAll
  strainBinsRaw <- strainBinsAllList$strainBinsRaw
  samplename <- strainBinsAllList$sampleName
  
  print("Calculating mouse genetic background")
  sampleResults <- BinScoring(strainBinsAll, keepHits=keepHits)
  
  resList <- PostProcessStrainBins(sampleResults, samplename, strainGroupAssignment, nBins, minBES, minGES = minGES, keepEmptyStrains=keepEmptyStrains)
  sampleResults <- resList[[1]]
  spDT <- resList[[2]]
  spSummary <- resList[[3]]
  
  # optional: add X to the plot based on "no-coverage" or "filtered" for each strain (depending on addCovMinGES)
  if(!is.null(brcFile)){
    print("Adding coverage/filtering X annotations to plot")
    AddCovRes <- BubblesAddCovX(brcFile, sampleResults, spDT, GenomeBins, strainDT, "", addCovMinGES, strainGroupAssignment)
    
    # overwrite the original results
    sampleResults <- AddCovRes$sampleResults
    spDT <- AddCovRes$spDT
  }
  
  # optional: add an additional annotation track for the "anti-strain; black6" hits
  if(!is.null(b6brcFile)){
    print("Adding black6 (C57BL_6J track) annotation to plot")
    B6BRC <- ProcessBRC(b6brcFile, minMQ=58, minReads=2, minBQ = 20, minAF=0.3, isB6 = T, B6data = B6DT)

    # stop if no SNPs were found
    if(nrow(B6BRC) > 0){
      B6.res <- DetermineBlack6(B6DT, B6BRC, GenomeBinsB6)
      # by adding this, we can remove the original NA track
      output <- PlotBlack6(B6.res, GenomeBins, strainGroupColors, samplename, sampleResults, spDT, mode = "withmain", annotationLoci, returnData=returnData)
    } else {
      print("No black6-SNPs were found, this will be NA in the plot.")
      
      # do it like below with the classic approach
      output <- PlotStrainBins(sampleResults, spDT, samplename, GenomeBins, strainGroupColors, annotationLoci, returnData=returnData)
    }
    
  } else {
    # classic plot
    output <- PlotStrainBins(sampleResults, spDT, samplename, GenomeBins, strainGroupColors, annotationLoci, returnData=returnData)
  }
  
  # classic but split by chromosome (optional and not default)
  if(plotByChrom){
    output.allchroms <- PlotStrainBinsByChrom(sampleResults, spDT, samplename, GenomeBins, strainGroupColors, annotationLoci)
    output[["plotObjByChrom"]] <- strainBinsRaw
  }
  
  if(returnData){
    output[["strainSNPs"]] <- strainBinsRaw
  }
  
  return(output)
}


# MHC - FUNCTIONS ------
RunAnalysisHLA <- function(mhcResources=mhcResources, samplePath, minCoverage=1, minAF=0.2, minMQ=58, minBQ=20, addTypeB=T){

  # default arguments (do not change these, they are only required in the strain analysis but the functions expect something)
  # TODO: CHECK WHICH ONES ARE NEEDED
  keepEmptyStrains=F
  returnData <- F
  minBES <- 0
  minGES <- 0
  minBinSNPS <- 1
  brcFile <- NULL
  b6brcFile <- NULL
  removeClusteredEvents <- F
  maxNeighbors <- 1
  combineTN <- "all" 
  selectFiles <- "TumorNormal"
  seqType <- "WES"
  annotationLoci <- NULL 
  
  # load the resource data from lookup object
  HLA.assignment <- mhcResources$HLA.assignment
  HLA.assignment2 <- mhcResources$HLA.assignment2
  HLA.cluster <- mhcResources$HLA.cluster
  HLA.colors <- mhcResources$HLA.colors
  HLA.K <- mhcResources$HLA.K
  HLA.A <- mhcResources$HLA.A
  HLA.E <- mhcResources$HLA.E
  HLA.D <- mhcResources$HLA.D
  HLA.Q <- mhcResources$HLA.Q
  HLA.T <- mhcResources$HLA.T
  GenomeBins <- mhcResources$GenomeBins
  genomeSNPs <- mhcResources$genomeSNPs
  mhcDT <- mhcResources$variantsMHC # strainDT was renamed to mhcDT in this function
  mhcDT.addedB6 <- mhcResources$variantsMHCB6
  B6HLA <- mhcResources$B6HLA
  strainGroupAssignment <- mhcResources$strainGroupAssignment # legacy, not really needed
  nBins <- GenomeBins[, length(unique(UnfiedPosition))]
  
  # load BRC
  sample.HLA.BRC.counts <- NULL
  if(addTypeB){
    samplePath2 <- gsub("MHC.brc","B6MHC.brc",samplePath)
    sample.HLA.BRC.counts <- ProcessBRC(brc.file=samplePath2, minReads=minCoverage, minMQ=minMQ, minBQ=minBQ, minAF=minAF, isHLA=T, B6data=B6HLA)
  }
  
  # want to add b6 negative snps?
  if(addTypeB){
    mhcDT <- rbind(mhcDT, mhcDT.addedB6)
  }
  
  print("Reading input data")
  
  if(addTypeB & !is.null(sample.HLA.BRC.counts)){
    # with hlatype b added
    strainBinsAllList <- MergeDataWithStrains(samplePath, mhcDT, GenomeBins, genomeSNPs, returnData, inputType="BRC", selectFiles, minCoverage, minAF, minMQ, minBQ, removeClusteredEvents, maxNeighbors, minBinSNPS, combineTN, addData = sample.HLA.BRC.counts)
  } else {
    # classic
    strainBinsAllList <- MergeDataWithStrains(samplePath, mhcDT, GenomeBins, genomeSNPs, returnData, inputType="BRC", selectFiles, minCoverage, minAF, minMQ, minBQ, removeClusteredEvents, maxNeighbors, minBinSNPS, combineTN)
  }
  
  if(is.null(strainBinsAllList)){return(NULL)}   # stop if no SNPs were found
  if(is.null(strainBinsAllList$strainBinsAll)){return(NULL)}   # stop if no SNPs were found
  strainBinsAll <- strainBinsAllList$strainBinsAll
  strainBinsRaw <- strainBinsAllList$strainBinsRaw
  samplename <- strainBinsAllList$sampleName
  
  print("Calculating mouse HLA background")
  sampleResults <- BinScoring(strainBinsAll, keepHits="all")
  
  resList <- PostProcessStrainBins(sampleResults, samplename, strainGroupAssignment, nBins, minBES, minGES, keepEmptyStrains=keepEmptyStrains, roundDigits = 5) # for HLA
  sampleResults <- resList[[1]]
  spDT <- resList[[2]]
  spSummary <- resList[[3]]
  
  res1 <- PlotV1HLA(sampleResults, spDT, samplename, HLA.cluster, HLA.assignment2, GenomeBins, HLA.colors, includeEmptyStrains = F)
  #res1
  
  res1b <- ScoreHLA(res1$data, GenomeBins, HLA.cluster, HLA.assignment2)
  #res1b
  
  return(list(res1=res1, res1b=res1b))
}

PlotV1HLA <- function(sampleResults, spDT, samplename, HLA.cluster, HLA.assignment2, GenomeBins, HLA.colors, annotationLoci=NULL, includeEmptyStrains=T){
  HLA.cluster.gr <- makeGRangesFromDataFrame(HLA.cluster[, .(chr, start, end)])
  
  plotDT <- merge(GenomeBins, sampleResults[, .(UnfiedPosition, strain, N, value, sample, strainGroup)], all.x=T, by="UnfiedPosition")
  plotDT <- plotDT[!is.na(strain)]
  
  plotGR <- makeGRangesFromDataFrame(plotDT)
  hits <- findOverlaps(plotGR, HLA.cluster.gr)
  plotDT[queryHits(hits), HLAlocus := HLA.cluster[subjectHits(hits), Cluster]]
  
  # assign loci
  plotDT <- plotDT[!is.na(HLAlocus)]
  
  plotDT <- merge(plotDT, HLA.assignment2, by=c("HLAlocus", "strain"))
  
  plotDT[, HLAlocus := factor(HLAlocus, levels = levels(HLA.assignment2$HLAlocus))]
  setorder(plotDT, HLAlocus)
  HLAtype.sorting <- plotDT[, sum(value), by=.(HLAlocus, HLAtype)]
  setorder(HLAtype.sorting, -V1)
  plotDT[, HLAtype := factor(HLAtype, levels = unique(HLAtype.sorting$HLAtype))]
  
  # I want every locus to be at least plotted with 3 bins
  # I also want to add 1-2 bins at the start or end of every region
  GenomeBins.tmp <- copy(GenomeBins)
  GenomeBins.gr <- makeGRangesFromDataFrame(GenomeBins.tmp)
  hits <- findOverlaps(GenomeBins.gr, HLA.cluster.gr)
  GenomeBins.tmp[queryHits(hits), HLAlocus := HLA.cluster[subjectHits(hits), Cluster]]
  GenomeBins.tmp1 <- GenomeBins.tmp[, .(plotstart=min(plotstart)-20000), by=HLAlocus]
  GenomeBins.tmp2 <- GenomeBins.tmp[, .(plotstart=max(plotend)+20000), by=HLAlocus]
  atype <- plotDT[, unique(HLAtype)][1] # just select any found type, else NA is fine
  GenomeBins.tmp1[, HLAtype := atype]
  GenomeBins.tmp2[, HLAtype := atype]
  plotDT <- rbind(plotDT, GenomeBins.tmp1, fill=T)
  plotDT <- rbind(plotDT, GenomeBins.tmp2, fill=T)
  
  if(includeEmptyStrains){
    missingStrains <- HLA.assignment[!Strains %in% plotDT[, unique(strain)], Strains]
    dummyDT <- HLA.assignment2[strain %in% missingStrains]
    dummyDT <- cbind(dummyDT, plotDT[1, -c("strain", "HLAtype", "HLAlocus")])
    setcolorder(dummyDT, colnames(plotDT))
    dummyDT[, value := NA]
    dummyDT[, plotstart := as.numeric(NA)]
    plotDT <- rbind(plotDT, dummyDT)
  }
  
  plotDT[, tmp := min(plotstart, na.rm = T), by=.(HLAlocus)]
  plotDT[is.na(plotstart), plotstart := tmp]
  plotDT <- plotDT[!is.infinite(plotstart)]
  plotDT <- plotDT[!is.na(HLAtype)]
  
  plotDT[, HLAtype := factor(HLAtype, levels = names(HLA.colors))]
  
  borderAdds <- plotDT[is.na(strain)] # needed to add empty boxes for missing loci (especially Q)
  plotDT <- plotDT[!is.na(strain)]
  
  outplot <- ggplot(plotDT, aes(plotstart, strain, color=HLAtype)) +
    geom_point(aes(size=value), alpha=0.8)  + 
    theme_grey() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid")) +
    scale_size_continuous("BES") + 
    facet_grid(HLAtype ~ HLAlocus, scales="free", space="free") +
    scale_color_manual("HLA type", values = HLA.colors, na.value = "darkgrey", guide="none") +
    geom_blank(data=borderAdds[is.na(strain)], aes(plotstart, strain)) #+ theme(legend.position="none", text = element_text(size=20)) 
  
  return(list(plot=outplot, data=plotDT))
}
PlotV2HLA <- function(sampleResults, spDT, samplename, HLA.cluster, HLA.assignment2, GenomeBins, HLA.colors, annotationLoci=NULL, includeEmptyStrains=T){
  HLA.cluster.gr <- makeGRangesFromDataFrame(HLA.cluster[, .(chr, start, end)])
  
  plotDT <- merge(GenomeBins, sampleResults[, .(UnfiedPosition, strain, N, value, sample, strainGroup)], all.x=T, by="UnfiedPosition")
  plotDT <- plotDT[!is.na(strain)]
  
  plotGR <- makeGRangesFromDataFrame(plotDT)
  hits <- findOverlaps(plotGR, HLA.cluster.gr)
  plotDT[queryHits(hits), HLAlocus := HLA.cluster[subjectHits(hits), Cluster]]
  
  plotDT <- merge(plotDT, HLA.assignment2, by=c("HLAlocus", "strain"), all.x=T)
  plotDT[, HLAtmp := paste(HLAlocus, HLAtype, sep="-")]
  
  strain.order <- plotDT[, sum(value), by=strain]
  setorder(strain.order, V1)
  
  plotDT[, strain := factor(strain, levels = strain.order$strain)]
  
  plotDT <- plotDT[!is.na(HLAtype)]
  plotDT[, HLAtype := factor(HLAtype, levels = names(HLA.colors))]
  
  outplot <- ggplot(plotDT, aes(plotstart, strain, color=HLAtype)) +
    geom_point(aes(size=value), alpha=0.8)  + 
    theme_grey() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid")) +
    coord_cartesian(xlim = c(2247000000, 2253000000)) +
    scale_size_continuous("BES") +
    scale_color_manual("HLA type", values = HLA.colors, na.value = "darkgrey", guide="none")
  
  if(!is.null(annotationLoci)){
    annotationLoci.gr <- makeGRangesFromDataFrame(annotationLoci)
    hits <- findOverlaps(annotationLoci.gr, GenomeBinsGR)
    
    markedNotFound <- annotationLoci[!queryHits(hits)]
    if(nrow(markedNotFound) > 0){
      warning(paste0("One or more annotation loci are not valid mm10 genomic locations: ", markedNotFound[, paste0(name, collapse = ",")]))
    }
    
    annotationLoci.plot <- cbind(annotationLoci[queryHits(hits)], GenomeBins[subjectHits(hits), .(plotstart, plotend)])
    annotationLoci.plot <- annotationLoci.plot[, .(plotstart=min(plotstart), plotend=max(plotend)), by=.(strainGroup, name, strain)]
    annotationLoci.plot[, center := plotstart + ((plotend-plotstart)/2)]
    
    outplot <- outplot + 
      geom_rect(data=annotationLoci.plot,aes(xmin=plotstart,xmax=plotend,ymin=-Inf,ymax=Inf),inherit.aes=FALSE, fill="black", alpha=0.05) +
      geom_text(data=annotationLoci.plot, aes(x=center, y = 0, label = name), color="black", size=3, vjust = -2)
    
    # reorder layers (i.e. put rect in background)
    outplot$layers <- outplot$layers[c(2,1,3)]
  }
  
  return(list(plot=outplot, data=plotDT))
}
ScoreHLA <- function(plotDT, GenomeBinsTmp, HLA.cluster, HLA.assignment2){
  plotDT <- plotDT[!is.na(strain)]
  
  HLA.cluster.gr <- makeGRangesFromDataFrame(HLA.cluster[, .(chr, start, end)])
  
  obsVals <- plotDT[, .(nBinsObs=uniqueN(UnfiedPosition), nStrainsObs=uniqueN(strain)), by=.(HLAlocus, HLAtype)]
  
  expVals1 <- HLA.assignment2[, .(nStrainsExp=uniqueN(strain)), by=.(HLAlocus, HLAtype)]
  
  GenomeBinsTmp.gr <- makeGRangesFromDataFrame(GenomeBinsTmp)
  hits <- findOverlaps(GenomeBinsTmp.gr, HLA.cluster.gr)
  GenomeBinsTmp[queryHits(hits), HLAlocus := HLA.cluster[subjectHits(hits), Cluster]]
  
  GenomeBinsTmp <- GenomeBinsTmp[!is.na(HLAlocus)]
  expVals2 <- GenomeBinsTmp[, .(nBinsExp=uniqueN(UnfiedPosition)), by=HLAlocus]
  
  scoreVals <- merge(obsVals, expVals1, by=c("HLAlocus", "HLAtype"))
  scoreVals <- merge(scoreVals, expVals2, by="HLAlocus")
  
  scoreVals[, scoreP := (nBinsObs*nStrainsObs) / (nBinsExp*nStrainsExp) * 100]
  scoreVals[, scoreP := round(scoreP, digits = 0)]
  
  scoreVals[, HLAlocus := factor(HLAlocus, levels = levels(HLA.assignment2$HLAlocus))]
  setorder(scoreVals, HLAlocus, -scoreP)
  
  bestVals <- scoreVals[scoreVals[, .I[scoreP == max(scoreP)], by=HLAlocus]$V1]
  
  allLoci <- HLA.assignment2[, unique(HLAlocus)]
  missingLoci <- allLoci[!allLoci %in% bestVals$HLAlocus]
  
  bestVals <- rbind(bestVals, data.table(HLAlocus=missingLoci), fill=T)
  
  setorder(bestVals, HLAlocus)
  
  # combine if multiple ones have the same score
  bestVals[, HLAtype := paste0(HLAtype, collapse = "/"), by=HLAlocus]
  bestVals <- unique(bestVals[, .(HLAlocus, HLAtype)])
  hlatype <- bestVals[, paste0(HLAtype, collapse = "-")]
  
  return(list(data=scoreVals, bestHaplotype=hlatype))
}
PrintCommandsBRC.MHC <- function(mhcResources, samplename, bamfile, workfolder, reffile, append=T){
  # VERSION: bam-readcount version: 0.8.0-unstable-50-f5f5ed0 (commit f5f5ed0)
  # this will write 2 commands into the file, once line for black6 and one for the other strains (it is not slower and cleaner this way)
  
  # load from resources
  HLA.K <- mhcResources$HLA.K
  HLA.A <- mhcResources$HLA.A
  HLA.E <- mhcResources$HLA.E
  HLA.D <- mhcResources$HLA.D
  HLA.Q <- mhcResources$HLA.Q
  HLA.T <- mhcResources$HLA.T
  
  bamPath <- dirname(bamfile)
  bamFile <- basename(bamfile)
  
  refPath <- dirname(reffile)
  refFile <- basename(reffile)
  
  outfile <- paste0(workfolder, "/", samplename, ".MHC.brc.tsv.gz")
  
  # write the lookup bed to the workdir
  
  HLA.locilist <- list(HLA.K, HLA.A, HLA.E, HLA.D, HLA.Q, HLA.T)
  HLA.BED <- do.call(rbind, HLA.locilist)
  HLA.BED <- unique(HLA.BED[, .(Chr, Pos, Pos)])
  bedfile <- paste0(workfolder, "/mhc_positions.bed")
  if(!file.exists(bedfile)){  fwrite(HLA.BED, bedfile, sep="\t", col.names = F)}
  bedfile <- basename(bedfile)
  
  # print commands to file
  command <- paste0("time docker run --user $(id -u):$(id -g) -v ",
                    refPath,":/reference/ -v ",bamPath,":/bam/ -v ",workfolder,":/workdir/ mgibio/bam-readcount -f /reference/",refFile,
                    " /bam/",bamFile," -l /workdir/",bedfile," -w 0 | cut -f1,2,3,4,6,7,8,9 | gzip > ", outfile)
  
  #commandFile <- paste0(workfolder, "/get_hla_read_counts_from_bam.sh")
  commandFile <- paste0(workfolder, "/get_read_counts_from_bam.sh")
  
  echostring <- paste0("echo \"Running bam-readcount (MHC, strain specific) for: \"", samplename)
  command <- paste0(echostring, " & ", command)
  cat("\n",file=commandFile,append=append)
  cat(command,file=commandFile,append=TRUE)
}
PrintCommandsBRC.MHC.B6 <- function(mhcResources, samplename, bamfile, workfolder, reffile, append=T){
  # VERSION: bam-readcount version: 0.8.0-unstable-50-f5f5ed0 (commit f5f5ed0)
  # this will write 2 commands into the file, once line for black6 and one for the other strains (it is not slower and cleaner this way)
  
  # load from resource
  strainselect.KAEDQT <- mhcResources$HLA.KAEDQT
  
  bamPath <- dirname(bamfile)
  bamFile <- basename(bamfile)
  
  refPath <- dirname(reffile)
  refFile <- basename(reffile)
  
  # FOR BLACK6
  outfile <- paste0(workfolder, "/", samplename, ".B6MHC.brc.tsv.gz")
  
  # write the lookup bed to the workdir
  hlab6BED <- unique(strainselect.KAEDQT[, .(Chr, Pos, Pos)])
  bedfile <- paste0(workfolder, "/b6mhc_positions.bed")
  if(!file.exists(bedfile)){  fwrite(hlab6BED, bedfile, sep="\t", col.names = F)}
  bedfile <- basename(bedfile)
  
  # print commands to file
  command <- paste0("time docker run --user $(id -u):$(id -g) -v ",
                    refPath,":/reference/ -v ",bamPath,":/bam/ -v ",workfolder,":/workdir/ mgibio/bam-readcount -f /reference/",refFile,
                    " /bam/",bamFile," -l /workdir/",bedfile," -w 0 | cut -f1,2,3,4,6,7,8,9 | gzip > ", outfile)
  
  #commandFile <- paste0(workfolder, "/get_hla_read_counts_from_bam.sh")
  commandFile <- paste0(workfolder, "/get_read_counts_from_bam.sh")
  
  echostring <- paste0("echo \"Running bam-readcount (MHC, black6) for: \"", samplename)
  command <- paste0(echostring, " & ", command)
  cat("\n",file=commandFile,append=append)
  cat(command,file=commandFile,append=TRUE)
}