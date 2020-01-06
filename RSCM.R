if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c('ShortRead','Rsubread','DESeq2','EnhancedVolcano','genefilter','pheatmap','Rsamtools'), update = T, ask = F)

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

################ DAY 2 ##############################################################

## Speed things up! Parallelize everything. I have 8 cores, so I'll use them all.
library("BiocParallel")
register(SnowParam(8))

## Map all the reads to genomic "features"(genes) with featureCounts
# All files with .BAM as the *end*
bams <- list.files(pattern = '.BAM$')
# Save the counts as a matrix
fc <- featureCounts(files=bams, annot.ext = 'gencode.v32.annotation.gtf.gz',
                    isGTFAnnotationFile = T, nthreads = 8)

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

rld <- rlog(IL13dds)
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
pheatmap(mat, annotation_col=df, cluster_cols = F)

# Cutoff the FC 
sig7v1ds <- IL13_7v1[which(IL13_7v1$log2FoldChange>=abs(1)),]

# Cutoff the pvalue (padj)
sig7v1ds <- sig7v1ds[which(sig7v1ds$padj<=0.1),]

# Write results
write.csv(as.data.frame(sig7v1ds[,c(7,2)],row.names = NULL), file = 'RSCM_IL13_7v1.csv')

# Test no FC cutoff
sig7v1ds <- IL13_7v1[which(IL13_7v1$padj<=0.1),]
write.csv(as.data.frame(sig7v1ds[,c(7,2)],row.names = NULL), file = 'RSCM_IL13_7v1.csv')
