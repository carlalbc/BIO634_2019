# BIO634 - Next generation sequencing (NGS) II. Transcriptomes, Variant Calling and Biological Interpretation
## June 3-4th 2019 
## University of Zürich (UZH) 

![alt text](https://github.com/carlalbc/URPP_tutorials/blob/master/img/Logo_URPP_kl2.png)

## URPP Evolution in action

# Day 2.- RNA sequencing: Transcriptomes and differential gene expression analyses

# I. Aligning transcriptomes with STAR

## Step 1: Read mapping with STAR aligner:

It is recommended using the STAR aligner for all genomes where there are no alternative alleles. For genomes such as hg38 that have alt alleles, hisat2 should be used as it handles the alts correctly and STAR does not yet. Use Tophat2 only if you do not have enough RAM available to run STAR (about 30 GB). The documentation for STAR is available [here](https://github.com/alexdobin/STAR/raw/master/doc/STARmanual.pdf)

Today we will work with data from the Zebrafish at different stages of differentiation. The data files are contained in the subdirectory called data and are the following:

- 2cells_1.fastq and 2cells_2.fastq: these files are based on RNA-seq data of a 2-cell zebrafish embryo, and
- 6h_1.fastq and 6h_2.fastq: these files are based on RNA-seq data of zebrafish embryos 6h post fertilization.

## Step 2: Install STAR aligner

```
sudo apt install rna-star
```

- The documentation for STAR is available [here](https://github.com/alexdobin/STAR/raw/master/doc/STARmanual.pdf)

## Step 3: Prepare all the directories and download the data

- Remember you can go from one directory to the next using ***cd***. Now let's create a new directory called STARGenome.

```
mkdir danRer 
mkdir danRer/reference danRer/data
cd danRer
wget ftp://ftp.ebi.ac.uk/pub/training/Train_online/RNA-seq_exercise/2cells_1.fastq ftp://ftp.ebi.ac.uk/pub/training/Train_online/RNA-seq_exercise/2cells_2.fastq ftp://ftp.ebi.ac.uk/pub/training/Train_online/RNA-seq_exercise/6h_1.fastq ftp://ftp.ebi.ac.uk/pub/training/Train_online/RNA-seq_exercise/6h_2.fastq -P data
wget ftp.ensembl.org/pub/release-93/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.toplevel.fa.gz -P reference
wget ftp://ftp.ensembl.org/pub/release-93/gtf/danio_rerio/Danio_rerio.GRCz11.93.gtf.gz -P reference
```

- If you feel adventurous you can run FastQC on all the reads in the **data** folder we just created using **cd** and find out the read length and quality:

```
fastqc *.gz
firefox *.html
```

## Step 4: Generate the genome index with STAR

- Run STAR in "genomeGenerate" mode

```
 STAR  [options]... --genomeDir REFERENCE   --readFilesIn R1.fq R2.fq
```

```
mkdir STARindex
gunzip reference/*.gz

STAR --runMode genomeGenerate --readFilesCommand zcat --genomeDir STARindex --genomeFastaFiles danRer11.fa --sjdbGTFfile danRer11.gtf --sjdbOverhang 75 --runThreadN 4

```

## Step 5: Perform the alignment with STAR:

- Make a folder to store the STAR output in it

```
mkdir alignment_STAR
```
- You can align the fastq files to the genome. The commands for this are explained in section 3.1 of the manual.12

Now align the pair of files from the 2cells sample to the genome, using the following parameters:


```
STAR --runThreadN 4 --genomeDir danRer11.fa --readFilesIn ../data/2cells_1.fastq  ../data/2cells_2.fastq --outSAMtype BAM SortedByCoordinate
```

Do the same for the other pair of reads (6h_1.fastq and 6h_2.fastq)

# II. Exploration of airway library: 

- To start let's install some R packages. In the terminal write the following:

## Step 1: Install R packages 

```
sudo R
source("https://bioconductor.org/biocLite.R")
biocLite(c("VennDiagram", "DESeq","edgeR", "Matrix", "airway", "Rsamtools", "pasilla", "GenomicFeatures", "GenomicAlignments","BiocParallel", "Rsubread"))
```
## Step 2: Open **rstudio** by typing ***rstudio*** in the command-line

Now go to the following [link](https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html
) to where it says *"Locating BAM files and the sample table"* and start from there.


# III. Differential analysis: Comparison between DESEq and edgeR


## Step 1 and 2: Open **rstudio** by typing ***rstudio*** in the command-line and install the packages like in the previous part.

- You can always copy&paste, but it is important that you understand what you are doing... go through chunks of the modules instead of everything at once. Each module is separated by a *#* which denotes a "comment" on the script. 

- Start by loading the libraries:

```
library(DESeq) 
library(edgeR) 
library(VennDiagram)


# Read in data ------------------------------------------------------------

## Use pasilla data (from Drosophila)
datafile = system.file( "extdata/pasilla_gene_counts.tsv", package="pasilla" )
datafile


## Read in the data making the row names the first column
counttable <- read.table(datafile, header=T, row.names=1)
head(counttable)
```

## Make metadata data.frame

```
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("untreated", "untreated", "untreated", "untreated", "treated", "treated", "treated"),
  libType=c("single", "single", "paired", "paired", "single", "paired", "paired"))
meta$condition <- relevel(meta$condition, ref="untreated")
meta

## Independent filtering?
# keep_cpm <- rowSums(cpm(counttable)>2) >=2
# keep_quantile <- rowSums(counttable)>quantile(rowSums(counttable), probs=.5)
# addmargins(table(keep_cpm, keep_quantile))
# counttable <- counttable[keep_cpm, ]

# DESeq -------------------------------------------------------------------

## Make a new countDataSet
d <- newCountDataSet(counttable, meta)

## Estimate library size and dispersion
d <- estimateSizeFactors(d)
d <- estimateDispersions(d)
plotDispEsts(d, main="DESeq: Per-gene dispersion estimates")

## Principal components biplot on variance stabilized data, color-coded by condition-librarytype
print(plotPCA(varianceStabilizingTransformation(d), intgroup=c("condition", "libType")))

## Fit full and reduced models, get p-values
dfit1 <- fitNbinomGLMs(d, count~libType+condition)
dfit0 <- fitNbinomGLMs(d, count~libType)
dpval <- nbinomGLMTest(dfit1, dfit0)
dpadj <- p.adjust(dpval, method="BH")

## Make results table with pvalues and adjusted p-values
dtable <- transform(dfit1, pval=dpval, padj=dpadj)
dtable <- dtable[order(dtable$padj), ]
head(dtable)

# Now with edgeR

# edgeR -------------------------------------------------------------------

## Make design matrix
condition <- relevel(factor(meta$condition), ref="untreated")
libType <- factor(meta$libType)
edesign <- model.matrix(~libType+condition)

## Make new DGEList, normalize by library size, and estimate dispersion allowing possible trend with average count size
e <- DGEList(counts=counttable)
e <- calcNormFactors(e)
e <- estimateGLMCommonDisp(e, edesign)
e <- estimateGLMTrendedDisp(e, edesign) 
e <- estimateGLMTagwiseDisp(e, edesign)

## MDS Plot
plotMDS(e, main="edgeR MDS Plot")

## Biological coefficient of variation plot
plotBCV(e, cex=0.4, main="edgeR: Biological coefficient of variation (BCV) vs abundance")

## Fit the model, testing the coefficient for the treated vs untreated comparison
efit <- glmFit(e, edesign)
efit <- glmLRT(efit, coef="conditiontreated")

## Make a table of results
etable <- topTags(efit, n=nrow(e))$table
etable <- etable[order(etable$FDR), ]
head(etable)

## ~MA Plot
with(etable, plot(logCPM, logFC, pch=20, main="edgeR: Fold change vs abundance"))
with(subset(etable, FDR<0.05), points(logCPM, logFC, pch=20, col="red"))
abline(h=c(-1,1), col="blue")

# Comparison --------------------------------------------------------------

head(etable)
head(dtable)

addmargins(table(sig.edgeR=etable$FDR<0.05, sig.DESeq=dtable$padj<0.05))

merged <- merge(etable, dtable, by='row.names')
with(                     merged, plot(logFC, conditiontreated, xlab="logFC edgeR", ylab="logFC DESeq", pch=20, col="black", main="Fold change for DESeq vs edgeR"))
with(subset(merged, FDR<0.05),  points(logFC, conditiontreated, xlab="logFC edgeR", ylab="logFC DESeq", pch=20, col="red"))
with(subset(merged, padj<0.05), points(logFC, conditiontreated, xlab="logFC edgeR", ylab="logFC DESeq", pch=20, col="green"))
legend("topleft", xjust=1, yjust=1, legend=c("FDR<0.05 edgeR only", "FDR<0.05 DESeq & edgeR", "FDR>0.05"), pch=20, col=c("red", "green", "black"), bty="n")

```

- Taken from "A survey of best practices for RNA-seq data analysis" https://doi.org/10.1186/s13059-016-0881-8
