#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#-#-#-#-#-#-#-#-#-#-#-#-#-# MoGeLi (mouse strains and MHC analysis) #-#-#-#-#-#-#-#-#-#-#-#-#-#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


# #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# 0. LOAD FUNCTIONS AND DATA #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
setwd("GITCLONEFOLDER/") # navigate the the git clone folder (of MoGeLi repository)

source("src/MoGeLi_functions.R") # load all functions 
strainResources <- readRDS("data/strainResources.RDS") # load strain SNP data
mhcResources <- readRDS("data/mhcResources.RDS") # load MHC SNP data
reffile <- "GRCm38.p6.fna" # should be downloaded, only used for bam-readcounts


# #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# 1. PREPROCESSING #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# extract all variant information directly from the BAM file

# the example is a sample from the JAX MMRDB WES sequencing (MMR-48497) 
# the sample is from a mixed strain background (A/J + B6C3Fe/C3H), additional can be found here information: https://www.jax.org/strain/1432
# the BAM file can not be provided, but it can be generated based on the available FASTQ files (SRR1784009)
samplename <- "MMR48497"
bamfile <- "MMR48497.bam" 
workfolder <- "/media/rad/HDD1/MoGeLi/BRC/" # for the bam-readcount output


# first we need to extract specific SNPs from the BAM file using a bash script calling bamreadcounts (docker)
# (skip this for the MMR-48497 example and continue with step 2 using the provided .brc files in the folder example/)
# alternatively, you can use your own custom set of called variants (e.g. from Mutect2, see below)

PrintCommandsBRC(strainResources, mhcResources, samplename, bamfile, workfolder, reffile, append=F)
# INPUT: BAM file
# OUTPUT (of this function): 1x command file (bash script with docker) and 4x bed files in the workfolder
# RUN (outside of R): parallel -j10 -a get_read_counts_from_bam.sh
# DESCRIPTION: 
  # this will print a bash file (get_read_counts_from_bam.sh) and 4 bed files into the workfolder
  # for each sample there will be 4 lines (strain, strain B6, mhc, mhc B6) 
  # this can be done for a large number of samples, by putting the above command in a loop and setting append=T
  # the command uses docker, so it can run on most systems without major dependencies
# OUTPUT (of bash file): 4x .brc.tsv.gz files (these are required for the strain and MHC analysis below)
# NOTE: if the sequencing depth is insufficient, this can be done for the Tumor and Normal and the resulting BRC files can be merged (only of both are from the same mouse)

# #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# 2. STRAIN ANALYSIS #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# INPUT: BRC files from PREPROCESSING
# OUTPUT: strain analysis bubble plot and data

# important parameters

sequencingType <- "WES" 
# this internally defines a set of default cutoffs for: WES, lcWGS or WGS 
# (the individual cutoffs can also be manually defined, see below)

# no need to change these (see step 1. PREPROCESSING)
brcfile.strain <- paste0(workfolder, "/", samplename, ".strains.brc.tsv.gz")
brcfile.strainB6 <- paste0(workfolder, "/", samplename, ".b6.brc.tsv.gz")

# run the main function
results.strain <- RunAnalysis(strainResources, brcfile.strain, inputType="BRC", seqType=sequencingType, 
                    brcFile = brcfile.strain, b6brcFile = brcfile.strainB6) 

# results plot
results.strain$plotObj
#ggsave(paste0(workfolder, "/", samplename, "_STRAIN.png"), results.strain$plotObj, width = 16, height = 9)


# WITH CUSTOM FILTERING (every sequencing run is different, here we can change some parameters to improve the results)
returnData <- F # output only the plot without additional data (SNP or bin assignment and percentage statistics)
# minBES <- 0.01 # minimal bin enrichment score (i.e. "how large should a bubble be to be included"). 
# Only change this, if there are too many small noisy bubbles (internally defined by sequencingType)
minGES <- 5 # minimal global enrichment score (i.e. "how much strain % (GES) is needed to be shown as track")
minMQ <- 60 # minimal mapping quality of reads to be used
minBQ <- 20 # minimal base quality of reads to be used
minAF <- 0.2 # minimal variant allele frequency of variants to be used
minCoverage <- 1 # minimal coverage of variants in the tumor to be used
addCovMinGES <- 100 # will add an additional track line with XXXX for filtered/not-covered reads (only for advanced QC, 100=deactivated)

# results.strain <- RunAnalysis(strainResources, brcfile.strain, inputType="BRC", seqType=sequencingType, returnData = returnData, 
#                     brcFile = brcfile.strain, b6brcFile = brcfile.strainB6,
#                     minGES = minGES, minMQ = minMQ, minBQ = minBQ, minAF = minAF, minCoverage = minCoverage, addCovMinGES=addCovMinGES) 

# if returnData=T, additional text outputs will be available
# results.strain$plotStats # percent values shown in the plot
# results.strain$plotBins # assignment and score (value) of each genomic bin to one/more specific strain(s)
# results.strain$strainSNPs # assignment and score (value) of each SNP to a bin and one/more specific strain(s)


# ALTERNATIVE DATA INPUT
if(F){
  # custom input data can be loaded and used
  # this is especially useful if mutation calling data (e.g. Mutect2) is available, which will reduce the overall noise in this analysis
  # the vcf output should be converted into a tab text file with the columns chr, pos, ref, alt
  # filtering should be performed following GATK guidelines and at your own discretion
  txtPath <- "example/MMR48497_Mutect2.txt.gz"
  txtRawData <- fread(txtPath)
  txtRawData <- txtRawData[, .(Chr=CHROM, Pos=POS, Ref=REF, Alt=ALT)] # required format and column names
  sampleSNPs <- copy(txtRawData)
  # sampleSNPs should be a global variable
  # FORMAT: data.table with 4 columns (Chr, Pos, Ref, Alt)
  results.strain <- RunAnalysis(strainResources, samplePath=NA, inputType="data", seqType="WES", returnData = T,
                      minGES = 3, minMQ = 60, minBQ = 30, minAF = 0.2, minCoverage = 5)
  
  results.strain$plotObj
}



# #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# 3. MHC ANALYSIS #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# INPUT: BRC files from PREPROCESSING
# OUTPUT: MHC analysis bubble plot and haplotype estimation

# no need to change these (see step 1. PREPROCESSING)
brcfile.mhc <- paste0(workfolder, "/", samplename, ".MHC.brc.tsv.gz")
brcfile.mhcB6 <- paste0(workfolder, "/", samplename, ".B6MHC.brc.tsv.gz")

results.mhc <- RunAnalysisHLA(mhcResources, samplePath=brcfile.mhc, addTypeB = T)
results.mhc$res1$plot # plot 
results.mhc$res1b$bestHaplotype # final haplotype prediction
# results.mhc$res1$data # assignment and score (value) of each SNP to a bin and one/more specific MHC types
# results.mhc$res1b$data # summary statistics of all MHC types 

ggsave(paste0(workfolder, "/", samplename, "_MHC.png"), results.mhc$res1$plot, width = 16, height = 9)

# the following parameters can be modified and added to the function 
minCoverage <- 1 # minimal coverage of variants in the tumor to be used
minMQ <- 60 # minimal mapping quality of reads to be used
minBQ <- 20 # minimal base quality of reads to be used
minAF <- 0.2 # minimal variant allele frequency of variants to be used
addTypeB <- T # the haplotype "b" is a special group and are added/calculated based on a very different approach as the other strains (works great, just letting you know)

# NOTE: custom data can not be used for this step, as most mutation calling tools do not output the highly variable MHC loci variants