## WARNING: DO NOT PERFORM OPERATIONS WITH SHORTREAD AND RSUBREAD (EXCEPT FOR FEATURECOUNTS)
## DURING THE WORKSHOP, AS IT IS MEMORY-INTENSIVE AND WILL TAKE A VERY LONG TIME

## INSTEAD, USE THE FILES PROVIDED OUTSIDE OF THE FOLDERS. FOR THE WORKSHOP, START WITH 
## THE LINE 'BAMS'

library(ShortRead)

## Define locations for raw files
paths <- list.dirs()
paths <- list.dirs()[2:length(paths)][c(F,T)]

## Access the raw file folders with readFastq()
for (x in paths){
  # Set the output name: We want the name of the folder where the fastq files are 
  # contained. Here it is the 3rd element of the 1st element ([[1]][3]) of the string x 
  # after we've split it (strsplit()) by the '/'. We want to add the extension .fastq, 
  # so we paste '.fastq' to the end of the name
  outname <- paste(strsplit(x,'/')[[1]][3],'.fastq', sep = '')
  # Get all the names of the fastqs in each folder 
  flist <- list.files(x, full.names = T)
  # Create another loop that gets the files in flist to add them into one file
  for (y in flist){
    # Read file with readFastq()
    yfile <- readFastq(y)
    writeFastq(yfile, outname, mode = 'a', compress = T)
  }
}

library(Rsubread)

## Build index, make sure it's gapped! Total time: 45.7 minutes using 10G Memory
buildindex('h_index','GRCh38.p13.genome.fa.gz', gappedIndex = T, memory = 10000)

## Make .bam files with align()
afiles <- list.files(pattern = '.fastq')
for (x in afiles){
  align('h_index',x)
}

## Map all the reads to genomic "features"(genes) with featureCounts
# All files with .BAM as the *end*
bams <- list.files(pattern = '.BAM$')
# Save the counts as a matrix
fc <- featureCounts(files=bams, annot.ext = 'gencode.v32.annotation.gtf.gz',
                    isGTFAnnotationFile = T)

## Rename observations according to metadata
meta <- read.csv('reshu_meta.csv')
colnames(fc$counts) <- meta$samples
## Check out what your matrix looks like for the entire experiment.
# Your column names should be the short form of the data files (I.e. 'Control-1')
# Rows should be the Ensemble gene ID's (I.e. 'ENSGxxxxx.5')
fc$counts

## Speed things up! Parallelize everything
library("Biocparallel")
register(SnowParam(8))




