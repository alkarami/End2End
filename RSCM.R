## Everything to install
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

## Make .bam files with align()
afiles <- list.files(pattern = '_R1.fastq')
bfiles <- list.files(pattern = '_R2.fastq')
for (x in 1:length(afiles)){
  f1 = afiles[x]
  f2 = bfiles[x]
  align('m2',f1,readfile2 = f2, nthreads = 4)
}

################ DAY 2 ##############################################################

## Map all the reads to genomic "features"(genes) with featureCounts
# All files with .BAM as the *end*
bams <- list.files(pattern = '.BAM$')
# Save the counts as a matrix
fc <- featureCounts(files=bams, annot.ext = 'gencode.v32.annotation.gtf.gz', isPairedEnd= T, 
                    GTF.attrType="gene_name", isGTFAnnotationFile = T, nthreads = 4)

## Rename observations according to metadata
meta <- read.csv('reshu_meta.csv')
colnames(fc$counts) <- meta$samples
## Check out what your matrix looks like for the entire experiment.
# Your column names should be the short form of the data files (I.e. 'Control-1')
# Rows should be the Ensemble gene ID's (I.e. 'ENSGxxxxx.5')
fc$counts

## Now we start with DESeq2
library(DESeq2)

# Create DESeqDataSet 
IL13ds <- DESeqDataSetFromMatrix(countData = fc$counts, colData = meta, design = ~group)

# Run DESeq2, which is the core of the code. It is *very* fast.
IL13dds <- DESeq(IL13ds)

# Run PCA to see how similar the samples are
rld <- rlog(IL13dds)
plotPCA(rld,intgroup = 'group')

# Let's check out results
IL13_3v1 <- results(IL13dds, contrast = c('group','3d','cont'), alpha = 0.05)
IL13_7v1 <- results(IL13dds, contrast = c('group','7d','cont'), alpha = 0.05)

# Add gene symbols to the results
library(stringr)
gene_ids <- str_replace(row.names(IL13dds),
                        pattern = ".[0-9]+$",
                        replacement = "")

BiocManager::install(c('AnnotationDbi','org.Hs.eg.db'), ask = F, update = T)

library("AnnotationDbi")
library("org.Hs.eg.db")
symbols <- mapIds(org.Hs.eg.db,
                  keys=gene_ids,
                  column="SYMBOL",
                  keytype="ENSEMBL",
                  multiVals="first")

IL13_3v1$Symbols <- symbols
IL13_7v1$Symbols <- symbols

# Check it out with a volcano plot
library(EnhancedVolcano)
EnhancedVolcano(IL13_3v1, IL13_7v1$Symbols,'log2FoldChange','padj')
EnhancedVolcano(IL13_7v1, IL13_7v1$Symbols,'log2FoldChange','padj')

# Make a heatmap
library(genefilter)
library(pheatmap)

topVarGenes <- head(order(-rowVars(assay(rld))),20)
mat <- assay(rld)[ topVarGenes, ]
heat_ids <- str_replace(row.names(mat),
                        pattern = ".[0-9]+$",
                        replacement = "")
heatsyms <- mapIds(org.Hs.eg.db,
                              keys=heat_ids,
                              column="SYMBOL",
                              keytype="ENSEMBL",
                              multiVals="first")
rownames(mat) <- heatsyms
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c('samples','group')])
pheatmap(mat, annotation_col=df, cluster_cols = F, cluster_rows = T)

# Cutoff the FC
sig3v1ds <- IL13_3v1[which(IL13_3v1$log2FoldChange>=abs(1)),]

# Cutoff the pvalue (padj)
sig7v1ds <- IL13_7v1[which(IL13_7v1$padj<=0.01),]

# Write results
write.csv(as.data.frame(IL13_3v1[,c(7,2,6)],row.names = NULL), file = 'RSCM_IL13_3v1.csv')
write.csv(as.data.frame(sig7v1ds[,c(7,2)],row.names = NULL), file = 'RSCM_IL13_7v1.csv')
