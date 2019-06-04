# BIO634 - Next generation sequencing (NGS) II. Transcriptomes, Variant Calling and Biological Interpretation
## June 3-4th 2019 
![alt text](https://github.com/carlalbc/URPP_tutorials/blob/master/img/Logo_URPP_kl2.png)
### University of Zürich (UZH) & URPP Evolution in action
----------
## Day 2.- RNA-seq and gene expression analyses
### I. Analyzing RNA-seq data with Salmon

*(Tutorial modified from the [Salmon](https://combine-lab.github.io/salmon/) official documentation)*

Today we will use **[Salmon](https://combine-lab.github.io/salmon/)** to align a transcriptome :fish:

[Salmon](https://combine-lab.github.io/salmon/) is a tool for quantifying the expression of transcripts using RNA-seq data. Salmon uses new algorithms (specifically, coupling the concept of *[quasi-mapping](https://academic.oup.com/bioinformatics/article/32/12/i192/2288985)* with a two-phase inference procedure) to provide accurate expression estimates very quickly (i.e. *wicked-fast*) and while using little memory. In Salmon it's all about quantification! 

:information_source: You could also use the STAR aligner, it is particularly good for all genomes where there are no alternatives alleles. For genomes such as hg38 that have alt alleles, hisat2 should be used as it handles the alts correctly and STAR does not yet. Use Tophat2 only if you do not have enough RAM available to run STAR (about 30 GB). The documentation for STAR is available [here](https://github.com/alexdobin/STAR/raw/master/doc/STARmanual.pdf).

### Salmon installation

Installation from source: 

```sh
# Make software directory in your home folder
mkdir ~/software
#Go to directory
cd ~/software
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

### a) Obtaining a transcriptome and building and index 
In order to quantify transcript-level abundances, Salmon requires a target transcriptome. 

This transcriptome is given to Salmon in the form of a (possibly compressed) multi-FASTA file, with each entry providing the sequence of a transcript. 

For this example, we’ll be analyzing some *Arabidopsis thaliana* data, so we’ll download and index the *A. thaliana* transcriptome. 
- First, create a directory where we’ll do our analysis, let’s call it `salmon`: 


```sh
# Make a working directory and go to it
mkdir salmon
cd salmon
```

- Download the transcriptome:

```sh
$ curl ftp://ftp.ensemblgenomes.org/pub/plants/release-28/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.28.cdna.all.fa.gz -o athal.fa.gz
```

Here, we’ve used a reference transcriptome for Arabadopsis. However, one of the benefits of performing quantification directly on the transcriptome (rather than via the host genome), is that one can easily quantify assembled transcripts as well (obtained via software such as StringTie for organisms with a reference or Trinity for de novo RNA-seq experiments).

Next, we’re going to build an index on our transcriptome. The index is a structure that salmon uses to quasi-map RNA-seq reads during quantification. The index need only be constructed once per transcriptome, and it can then be reused to quantify many experiments. 

- Now we use the index command of salmon to build our index:


```
salmon index -t athal.fa.gz -i athal_index
```
More info on parameters for indexing [here](https://salmon.readthedocs.io/en/latest/)

### b) Obtaining sequencing data

In addition to the index, salmon obviously requires the RNA-seq reads from the experiment to perform quantification. In this tutorial, we’ll be analyzing data from this [4-condition experiment](https://www.ebi.ac.uk/ena/data/view/PRJDB2508)[accession PRJDB2508]. You can use the following shell script to obtain the raw data and place the corresponding read files in the proper locations. Here, we’re simply placing all of the data in a directory called `data`, and the left and right reads for each sample in a sub-directory labeled with that sample’s ID (i.e. `DRR016125_1.fastq.gz` and `DRR016125_2.fastq.gz` go in a folder called `data/DRR016125`).

```sh
#!/bin/bash
mkdir data
cd data
for i in `seq 25 28`; 
do 
  mkdir DRR0161${i}; 
  cd DRR0161${i}; 
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR016/DRR0161${i}/DRR0161${i}_1.fastq.gz; 
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR016/DRR0161${i}/DRR0161${i}_2.fastq.gz; 
  cd ..; 
done
cd .. 
```
- Place this commands in a script called `download_reads.sh`. You can use Gedit:
`gedit download_reads.sh`
- Copy and paste the contents above, save the file with `Ctrl+S` and Exit with `Alt+F4`
- To download the data, run the scrip and wait for it to complete:
```
bash download_reads.sh
```

### c) Quantifying the samples

Now that we have our index built and all of our data downloaded, we’re ready to quantify our samples. Since we’ll be running the same command on each sample, the simplest way to automate this process is, again, a simple shell script (`quant_samples.sh`):

```sh
#!/bin/bash
for fn in data/DRR0161{25..28};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i athal_index -l A \
         -1 ${fn}/${samp}_1.fastq.gz \
         -2 ${fn}/${samp}_2.fastq.gz \
         -p 8 --validateMappings -o quants/${samp}_quant
done 
```

- Place this commands in a script called `quant_samples.sh`. You can use Gedit:
`gedit quant_samples.sh`
- Copy and paste the contents above, save the file with `Ctrl+S` and Exit with `Alt+F4`
- To download the data, run the scrip and wait for it to complete:
```sh
bash quant_samples.sh
```

This script simply loops through each sample and invokes `salmon` using fairly barebone options. The `-i` argument tells salmon where to find the index `-l A` tells salmon that it should automatically determine the library type of the sequencing reads (e.g. stranded vs. unstranded etc.). The `-1` and `-2` arguments tell salmon where to find the left and right reads for this sample (notice, salmon will accept gzipped FASTQ files directly). Finally, the `-p 8` argument tells salmon to make use of 8 threads and the `-o` argument specifies the directory where salmon’s quantification results sould be written. Salmon exposes *many* different options to the user that enable extra features or modify default behavior. However, the purpose and behavior of all of those options is beyond the scope of this introductory tutorial. You can read about salmon’s many options in the [documentation](http://salmon.readthedocs.io/en/latest/).

After the salmon commands finish running, you should have a directory named `quants`, which will have a sub-directory for each sample. These sub-directories contain the quantification results of salmon, as well as a lot of other information salmon records about the sample and the run. The main output file (called `quant.sf`) is rather self-explanatory. For example, take a peek at the quantification file for sample `DRR016125` in `quants/DRR016125/quant.sf` and you’ll see a simple TSV format file listing the name (`Name`) of each transcript, its length (`Length`), effective length (`EffectiveLength`) (more details on this in the documentation), and its abundance in terms of Transcripts Per Million (`TPM`) and estimated number of reads (`NumReads`) originating from this transcript.

### d) After quantification

That’s it! Quantifying your RNA-seq data with salmon is that simple (and fast). Once you have your quantification results you can use them for downstream analysis with differential expression tools like [DESeq2](https://bioconductor.org/packages/DESeq2), [edgeR](https://bioconductor.org/packages/edgeR), [limma](https://bioconductor.org/packages/limma), or [sleuth](http://pachterlab.github.io/sleuth/). Using the [tximport](http://bioconductor.org/packages/tximport) package, you can import salmon’s transcript-level quantifications and optionally aggregate them to the gene level for gene-level differential expression analysis. You can read more about how to import salmon’s results into DESeq2 by reading the tximport section of the excellent [DESeq2 vignette](https://bioconductor.org/packages/DESeq2). For instructions on importing for use with edgeR or limma, see the [tximport vignette](http://bioconductor.org/packages/tximport). For preparing salmon output for use with sleuth, see the [wasabi](https://github.com/COMBINE-lab/wasabi) package.

### II. Exploration of airway library: 

- To start let's install some R packages. In the terminal write the following:

#### Step 1: Install R packages 

```
sudo R
source("https://bioconductor.org/biocLite.R")
biocLite(c("VennDiagram", "DESeq","edgeR", "Matrix", "airway", "Rsamtools", "pasilla", "GenomicFeatures", "GenomicAlignments","BiocParallel", "Rsubread"))
```
## Step 2: Open **rstudio** by typing ***rstudio*** in the command-line

Now go to the following [link](https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html
) to where it says *"Locating BAM files and the sample table"* and start from there.


### III. Differential analysis: Comparison between DESEq and edgeR


#### Step 1 and 2: Open **rstudio** by typing ***rstudio*** in the command-line and install the packages like in the previous part.

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

### Make metadata data.frame

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


# Useful workflows for RNA-seq data analyses

- [RNA-seq workflow: gene-level exploratory analysis and differential expression, 2018](https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html)
- Importing transcript abundance datasets with [tximport](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html), 2019.

-  [Analysis of Genomics Data with R/Bioconductor](https://ivanek.github.io/analysisOfGenomicsDataWithR/06_RNAseqFromFASTQtoCountData_html.html)
- [ARMOR (Automated Reproducible MOdular RNA-seq) a 2019 Snakemake RNA-seq workflow](https://github.com/csoneson/ARMOR).
- [RNA-seq workflow - gene-level exploratory analysis and differential expression, 2016](https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html).







