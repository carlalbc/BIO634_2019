# BIO634 - Next generation sequencing (NGS) II. Transcriptomes, Variant Calling and Biological Interpretation
## June 3-4th 2019 
## University of Zürich (UZH) 

![alt text](https://github.com/carlalbc/URPP_tutorials/blob/master/img/Logo_URPP_kl2.png)

## URPP Evolution in action

# Day 2.- RNA sequencing: Transcriptomes and differential gene expression analyses

# I. Aligning transcriptomes with Salmon

Today we will use **[Salmon](https://combine-lab.github.io/salmon/)** to align a transcriptome :fish:

[Salmon](https://combine-lab.github.io/salmon/) is a tool for quantifying the expression of transcripts using RNA-seq data. Salmon uses new algorithms (specifically, coupling the concept of *quasi-mapping* with a two-phase inference procedure) to provide accurate expression estimates very quickly (i.e. *wicked-fast*) and while using little memory. In Salmon it's all about quantification! 

:information_source: You could also use the STAR aligner, it is particularly good for all genomes where there are no alternatives alleles. For genomes such as hg38 that have alt alleles, hisat2 should be used as it handles the alts correctly and STAR does not yet. Use Tophat2 only if you do not have enough RAM available to run STAR (about 30 GB). The documentation for STAR is available [here](https://github.com/alexdobin/STAR/raw/master/doc/STARmanual.pdf).

### Salmon installation

Installation from source: 

```sh
# Make software directory in your home folder
mkdir ~/software
# Get the source file
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.14.0/salmon-0.14.0_linux_x86_64.tar.gz
# Uncompress the file
tar xzvf salmon-0.14.0_linux_x86_64.tar.gz
# Export salmon into your PATH variable
export PATH=$PATH:~/software/salmon-latest_linux_x86_64/bin/
```
:information_source: If you had conda installed in your system, you could also install it with conda which would be much faster.

- Run salmon in the terminal:
```
$ salmon -h
salmon v0.14.0

Usage:  salmon -h|--help or 
        salmon -v|--version or 
        salmon -c|--cite or 
        salmon [--no-version-check] <COMMAND> [-h | options]

Commands:
     index Create a salmon index
     quant Quantify a sample
     alevin single cell analysis
     swim  Perform super-secret operation
     quantmerge Merge multiple quantifications into a single file

```


## Step 1: Analyzing RNA-seq data with Salmon

In order to quantify transcript-level abundances, Salmon requires a target transcriptome. This transcriptome is given to Salmon in the form of a (possibly compressed) multi-FASTA file, with each entry providing the sequence of a transcript. For this example, we’ll be analyzing some Arabidopsis thaliana data, so we’ll download and index the A. thaliana transcriptome. First, create a directory where we’ll do our analysis, let’s call it `salmon`: 



```
# Make a working directory and go to it
mkdir salmon
cd salmon
```

Here, we’ve used a reference transcriptome for Arabadopsis. However, one of the benefits of performing quantification directly on the transcriptome (rather than via the host genome), is that one can easily quantify assembled transcripts as well (obtained via software such as StringTie for organisms with a reference or Trinity for de novo RNA-seq experiments).

Next, we’re going to build an index on our transcriptome. The index is a structure that salmon uses to quasi-map RNA-seq reads during quantification. The index need only be constructed once per transcriptome, and it can then be reused to quantify many experiments. We use the index command of salmon to build our index:


```
salmon index -t athal.fa.gz -i athal_index
```

More info on parameters for indexing [here](https://salmon.readthedocs.io/en/latest/)


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




:information_source: Today we will work with the Jupyter notebook

## Jupyter notebook and R installation

- For installing under Debian (command under Ubuntu are analogous), as root:

```sh
#Become root
sudo su        #Give password "URPP$2019" for the student account and press enter

#In case python 3 it's not installed
apt-get install build-essential python3-dev python3-pip
```

- Then, after the installation (again, as root) run:

`pip3 install jupyter`

:information_source: If the last command was launched without root privileges, Jupyter is not available in the command line, but it can be launched:   `~/.local/bin/jupyter-notebook`

:information_source: If Jupyter was installed for all the users (with root priviledges), then simply type in the terminal:
`jupyter-notebook`

- The R installation is, with root privileges, as easy as:

`apt-get install r-base r-base-dev libssl-dev libcurl3-dev curl`

- The installation of the R kernel for Jupyter is performed under R command line (as from the GitHub page of the project, https://github.com/IRkernel/IRkernel):

```sh
#Enter R with root priviledges in the terminal, give the same password than before, URPP$2019 and type:
sudo R
```
- Run the following in R:
```r
install.packages(c('pbdZMQ', 'repr', 'devtools')) 
devtools::install_github('IRkernel/IRkernel') 
IRkernel::installspec()
```
Now the kernel is installed. Launch Jupyter (by typing `jupyter notebook`) and open the file given [here]() from the Jupyter notebook






