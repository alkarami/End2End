## Everything to install
install.packages('BiocManager')
BiocManager::install(c('ShortRead','Rsubread','DESeq2','EnhancedVolcano','genefilter','pheatmap','Rsamtools'), update = T, ask = F)

## Speed things up! Parallelize everything. I have 8 cores, so I'll use them all.
library("BiocParallel")
register(SnowParam(8))

library(ShortRead)

## Define locations for raw files
paths <- list.dirs()
paths <- list.dirs()[2:length(paths)]

## Access the raw file folders with readFastq()
for (x in paths){
  # Set the output name: We want the name of the folder where the fastq files are 
  # contained. Here it is the 2nd element of the 1st element ([[1]][2]) of the string x 
  # after we've split it (strsplit()) by the '/'. We want to add the extension .fastq, 
  # so we paste '.fastq' to the end of the name
  outname_1 <- paste(strsplit(x,'/')[[1]][2],'_R1.fastq', sep = '')
  outname_2 <- paste(strsplit(x,'/')[[1]][2],'_R2.fastq', sep = '')
  # Get all the names of the fastqs in each folder 
  flist_1 <- list.files(x, full.names = T, pattern = '_R1_')
  flist_2 <- list.files(x, full.names = T, pattern = '_R2_')
  # Create another loop that gets the files in flist to add them into one file
  for (y in flist_1){
    # Read file with readFastq()
    yfile <- readFastq(y)
    writeFastq(yfile, outname_1, mode = 'a', compress = T)
  }
  for (y in flist_2){
    # Read file with readFastq()
    yfile <- readFastq(y)
    writeFastq(yfile, outname_2, mode = 'a', compress = T)
  }
}

library(Rsubread)

## SKIP THIS PART IF YOU ALREADY HAVE THE INDEX! 
## Build index, make sure it's gapped! Total time: 45.7 minutes using 10G Memory
buildindex('m2','GRCm38.p6.genome.fa.gz', gappedIndex = T)

################ DAY 2 ##############################################################

## Make .bam files with align()
afiles <- list.files(pattern = '_R1.fastq')
bfiles <- list.files(pattern = '_R2.fastq')
align('m2',afiles,readfile2 = bfiles, nthreads = 4)

# Run featurecounts to count reads mapped to each feature (gene)
bams <- list.files(pattern = '.BAM$')
fc <- featureCounts(files=bams, annot.ext = 'gencode.vM24.annotation.gtf.gz', isPairedEnd= T, 
                    GTF.attrType="gene_name", isGTFAnnotationFile = T, nthreads = 4)

# Create metadata, instructions on how to group the samples
samples <- colnames(fc$counts)
ed.samples <- c()
for (x in samples){
  ed.samples <- append(ed.samples,strsplit(x,split = 'R1')[[1]][1])
}
colnames(fc$counts) <- ed.samples
CytoID <- c('Neg','Pos','Neg','Pos','Neg','Pos','Neg','Pos')
Sex <- c('F','F','F','F','M','M','M','M')
meta <- as.data.frame(cbind(ed.samples,CytoID,Sex))

library(DESeq2)

# Create DESeqDataSet 
CytoIDds <- DESeqDataSetFromMatrix(countData = fc$counts, colData = meta, design = ~Sex+CytoID)

# Run DESeq2, which is the core of the code. It is *very* fast.
CytoIDdds <- DESeq(CytoIDds)

# PCA
rld <- rlog(CytoIDdds)
p1 <- plotPCA(rld,intgroup = 'CytoID')
p2 <- plotPCA(rld,intgroup = 'Sex')
plot_grid(p1,p2, ncol = 1)

# Get results for DEGs
CytoID_HvL <- results(CytoIDdds, contrast = c('CytoID','Pos','Neg'), alpha = 0.01)
Sex_FvM <- results(CytoIDdds, contrast = c('Sex','F','M'), alpha = 0.05)

# Make a heatmap
library(genefilter)
library(pheatmap)

topVarGenes <- head(order(-rowVars(assay(rld))),20)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c('CytoID','Sex')])
pheatmap(mat, annotation_col=df, cluster_cols = T, cluster_rows = T)

