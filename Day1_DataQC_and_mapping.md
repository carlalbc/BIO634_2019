# BIO634 - Next generation sequencing (NGS) II. Transcriptomes, Variant Calling and Biological Interpretation

## June 3 - 4th, 2019
### University of Zürich (UZH) & URPP "Evolution in action"

![alt text](https://github.com/carlalbc/URPP_tutorials/blob/master/img/Logo_URPP_kl2.png)

## Part I.- Data Quality Control (QC), pre-processing and mapping genomes 

Before we get started --remember, Linux is your friend :penguin:

:information_source: A brief reminder of useful commands that will come in handy for Today's workflow:

| Command | Function |
| ------ | ------ |
| ls   | listing files
less, more | view files
head | view first 10 lines of the file
cd | change directory (enter)
cd ..   |  change one directory down (exit)
mv      |             moving file/ dir
cp      |            copying file / dir
mkdir   |           make directory
rmdir   | remove directory
wget    |          download files from the web  


## I) Quality Assesment: Pre-processing of the reads

We will be working with data from the [Lenski lab](http://myxo.css.msu.edu/ecoli/genomicsdat.html) and Desai. They are known for making long-term evolution experiments in *E. coli* since the early 00’s. The strain we will work with today is an *E. coli* from a long-term evolution experiment (LTEE). The twelve LTEE populations have been serially propagated
in the same medium for more than 60,000 generations, with samples preserved every 500 generations.  

- The sequencing data is from: [Good, B. H., McDonald, M. J., Barrick, J. E., Lenski, R. E., & Desai, M. M. (2017). The dynamics of molecular evolution over 60,000 generations. Nature, 551(7678), 45–50. doi:10.1038/nature24287](https://www.nature.com/articles/nature24287
).

Let's get started! 

## Step 1: Downloading raw sequencing FASTQ files from a database.

- There are two main databases, the **Sequence Read Archive** (SRA, US based) and the **European nucleotide archive** (ENA, EU based). 

- Today we will download the raw reads from the ENA. The project accession name is [PRJNA380528](https://www.ebi.ac.uk/ena/data/view/PRJNA380528). Please, follow the next steps before downloading the fastq files:
 
####  1. Make the following directories and go to the main directory (remember it's good to keep things tidy!):


```sh
# 1) Make directory called mapping:
mkdir mapping

# 2) Make subdirectory called fastq:
mkdir mapping/fastq           

# 3) Make another subdirectory for the fastq files called SRR6170103
mkdir mapping/fastq/SRR6170103       

# 4) Go to the SRR6170103 directory inside the mapping/fastq directories
cd mapping/fastq/SRR6170103 

# 5)  Make yet another directory at your current location to store the QC results later
mkdir FastQC                        
``` 

#### 2. Download the FASTQ files from ENA - we will work with paired-end (PE) reads from *E. coli*:

- Verify that you are the directory `mapping/fastq/SRR6170103` by using `pwd`. If you are continue with the below commands to download the files that we will work on today.

```sh 
# Get both FASTQ PE read files from ENA and store it in the subdirectory we just created:

wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR617/003/SRR6170103/SRR6170103_1.fastq.gz && wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR617/003/SRR6170103/SRR6170103_2.fastq.gz 
```

- Let’s check the FASTQ files:

```sh 
less SRR6170103_1.fastq.gz        # Exit with Ctrl+Z
# shows the first 10 lines
head SRR6170103_1.fastq.gz        # In this case it's a compressed file so it will be non-human readable!
``` 

- You can also use `tail` to see the end of the file.

Now that you have seen the files, continue with the rest of the workflow.

## Step 2: Quality check of the FASTQC files by using FastQC

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a quality control tool for high throughput sequence data.

When running new software it is always useful to understand it first. A quick glimpse to different options can be obtained by looking at the in-built help:

```sh
# Find help in the help manual of FastQC

fastqc --help	

Usage: fastqc seqfile1 seqfile2 .. seqfileN
```

#### 1) Assess the quality of the FASTQC files containing PE sequencing reads by running FastQC:

- We can do this either by using the graphical user interface of the program (GUI) or through the command-line (recommended). 
 
##### Graphical User Interface (GUI) - option 1:

```sh
#Open the GUI of FastQC by typing fastqc in the command-line (don't forget the ampersand)

fastqc &
```

That will open FastQC and you will be able to open the fastq files directly with the program. If you prefer to use the command-line (recommended) do the following:

##### Command-line -  option 2:

```sh  
#Run FastQC on both files and wait till it's done running
fastqc SRR6170103_1.fastq.gz SRR6170103_2.fastq.gz      

# Keep it tidy by moving the resulting files to the FastQC folder we created at the beginning
mv *.zip *.html FastQC         

# Go to the FastQC folder
cd FastQC
``` 

:information_source: **Reminder:** You can always check where you are in the terminal using `pwd`


#### 2) Open the FastQC results with your favorite html visualizer (i.e firefox, chrome, etc.) or if you prefer it, you can open the file through your GUI by directly clicking on it.

```sh 
firefox SRR6170103_1_fastqc.html
```

- You will see the following:
  
![alt text](https://github.com/carlalbc/BIO694_2018/blob/master/img/fastqc_report1.png)

Pretty good quality reads! :heavy_check_mark: :octocat:



- :information_source: Sometimes you can get very **bad** quality reads. See the example below:

![alt text](https://github.com/carlalbc/BIO694_2018/blob/master/img/fastqc_bad1.png)

Those are quite bad quality reads! :heavy_multiplication_x: :warning: :bangbang: :octocat:

:diamond_shape_with_a_dot_inside: Go through this [MANUAL](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf) to understand each of the results you have gotten in your report. And answer the following:

### :beginner: Questions:


1. What’s the total number of sequences in each of the paired-end fastq files?  (Hint: look at the Basic statistics module)

- SRR6170103_1 Total number of sequences  _________

- SRR6170103_2 Total number of sequences  _________

2. What’s the type of encoding used?
3. What’s the length of the reads? 
4. What are the warnings we get for each of the fastq files? 
5. Which sequence seems to be overrepresented? 

There is usually an expected drop in quality at the 3’ end of the sequences as well as an expected overrepresentation of adaptor sequences. We will learn how to trim low quality ends and remove adaptors.

**Do not close the FastQC window as we will compare the original files to the ones we will produce after adapter removal and quality filtering.**

## Step 3: Use Trimmomatic for trimming, removing adapter sequences  and low quality reads

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is a java tool for performing a range of trimming tasks on Illumina paired end and single end read data. Here's the [full manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

The parameters used for Trimmomatic are defined as follows:

| Parameter | Definition |
| ------ | ------ |
 **PE** | Data is paired end
**phred33** | Quality scores are 33 offset
**threads 1** | Number of threads to use
**trimlog logfile** | Name of logfile for summary information
**paired_end1.fastq** | Name of input fastq file for left reads
**paired_end2.fastq** | Name of input fastq file for right reads
**Left_paired.fastq** | Paired trimmed output fastq file for left reads
**Left_unpaired.fastq** | Unpaired trimmed output fastq file for left reads
**Right_paired.fastq** | Paired trimmed output fastq file for right reads
**Right_unpaired.fastq** | Unpaired trimmed output fastq file for right reads
**ILLUMINACLIP** | Parameters for the adapter clipping
**TruSeq3-PE-2.fa** | Text file of adapter sequences to search for
**:2:40:15** | Adapter-read alignment settings 
**MINLEN:36** | Delete reads trimmed below length MINLEN

#### a) Trimming PE reads and removing adapters:

```sh
# Trimming and removing Illumina adapters
java -jar /usr/share/java/trimmomatic.jar PE -phred33  SRR6170103_1.fastq.gz SRR6170103_2.fastq.gz SRR6170103_1_paired.fastq.gz SRR6170103_1_unpaired.fastq.gz SRR6170103_2_paired.fastq.gz SRR6170103_2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### :beginner: Questions: 

1. What does :2:40:15 means?      (Hint: It has something to do with the adapters, check Trimmomatic’s manual)
2. How many reads survived after Trimmomatic? (Hint: Check the messages left in the terminal after Trimmomatic)


:information_source: **NOTE**: Trimmomatic only deletes reads if the length after trimming of adapter sequences is less than MINLEN (which we set to 36bp)

#### b) Use Trimmomatic to filter low quality reads

Next, we can remove low quality reads of the sequences by trimming the bases at the 3' end of the reads with the following command:

```sh
#Filtering low quality reads
java -jar /usr/share/java/trimmomatic.jar PE -phred33 -threads 1 -trimlog logfile2 SRR6170103_1_paired.fastq.gz SRR6170103_2_paired.fastq.gz SRR6170103_1_trim_paired.fastq SRR6170103_1_unpaired.fastq SRR6170103_2_trim_paired.fastq SRR6170103_2_trim_unpaired.fastq SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36
```

### :beginner: Questions:

- What was the screenlog of Trimmomatic this time? How many reads were removed? Whas there any?
- What does the SLIDINGWINDOW:4:15 means? Check the manual.

## Step 4: Correcting sequencing errors by using SOAPec 

[SOAPec](https://omictools.com/soapec-tool) is a tool is a read correction tool specifically designed for illumina short reads, the manual can be found [here](http://soap.genomics.org.cn/about.html). The command below will take the trimmed fastq files generated in step 4 and correct sequencing errors in a two step process using tools called KmerFreq_AR and Corrector_AR.

:diamond_shape_with_a_dot_inside: If not installed already, download SOAPec

```sh
#Make directory
mkdir software
# Get SOAPec
wget http://sourceforge.net/projects/soapdenovo2/files/ErrorCorrection/SOAPec_v2.01.tar.gz -P software
# Go to directory
cd software
#Untar
tar -zxf SOAPec_v2.01.tar.gz
#Get out of folder
cd ..
```

Now we need to create a file with the fastq files location for SOAPec to be able to work. Use your favorite text editor (ie. gedit) and write the path of our fastq files. If you can't remember, use ***pwd**.

```
gedit files.txt &
```

It will open up gedit where you have to write the path of the files that you created in Step 4. 
In **my** case it will look like this:

```
/home/vega/URPP_2018/BIO634-2019/PartI/fastq/SRR6170103/SRR6170103_1_trim_paired.fastq
/home/vega/URPP_2018/BIO634-2019/PartI/fastq/SRR6170103/SRR6170103_2_trim_paired.fastq
```
- Save the file (Ctrl+S) and exit (Alt+F4). 
- Run the `KmerFreq_AR` command below and when it finishes run the `Corrector_AR` command.

:warning:Depending on your location and where you downloaded SOAPEc *the path will change slightly*. Remember to use **pwd** to know where you are :warning:

```
#Substitute the path below for your own
/home/vega/URPP_2018/BIO634-2019/software/SOAPec_v2.01/bin/KmerFreq_AR -k 16 -t 1 -q 33 -p Error_Corr files.txt > kmerfreq16.log 2> kmerfreq16.err
/home/vega/URPP_2018/BIO634-2019/software/SOAPec_v2.01/bin/Corrector_AR -k 16 -Q 33 -t 1 -o 3 Error_Corr.freq.cz Error_Corr.freq.cz.len files.txt > Corr16.log 2>Corr16.err
```
### :beginner: Questions

1. Inspect the output of Corrector_AR output file SRR6170103_1_trim_paired.fastq.cor.stat (Hint, use ***more*** / ***less***)

For SRR6170103_1:

- How many bases were corrected by the "Fast method"?  ______________

- How many bases were corrected by the "BBtree method"?  ______________

What about for SRR6170103_2?

- How many bases were corrected by the "Fast method"?  ______________

- How many bases were corrected by the "BBtree method"?  ______________

- **Are there any differences between the two? What conclusions can you make?**

_____________________________________________________________________________________


**If you closed FastQC before, open it and run it for the trimmed paired-end reads, do they look very different?**


Congratulations you now can continue with the next step of our workflow! :octocat:

NOTE: Exit the fastq folder using ***cd ..*** until you get to your main directory

## II) Mapping sequencing reads to a reference genome using the Burrows-Wheeler Aligner (BWA) tool


BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. You can find the manual [here](http://bio-bwa.sourceforge.net/bwa.shtml)
 
- In this part of the hands-on session we will map the cleaned reads from the previous steps to the reference genome of *E. coli*.

## Step 1: Download the reference genome and its annotation file

Go to the main folder (we called it mapping) and download the files.

```
cd mapping
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz
gunzip *
```

## Step 2: Create indices and dictionaries for bwa, samtools and picard.

Indices are necessary for quick access to specific information in very large files.

As you go creating the indices check for the files with ***ls***

- bwa index
```
bwa index -a is GCF_000005845.2_ASM584v2_genomic.fna
```
***-a is*** Sets the algorithm to be used to construct a suffix array. This is suitable for
databases smaller than 2GB.

- samtools index
```
samtools faidx GCF_000005845.2_ASM584v2_genomic.fna
```
- picard index
```
java -jar /home/student/APPL/PICARD/picard.2.18.0.jar CreateSequenceDictionary R=GCF_000005845.2_ASM584v2_genomic.fna O=GCF_000005845.2_ASM584v2_genomic.dict
```

# Step 3: Align reads to the Reference Genome using BWA

```
bwa mem GCF_000005845.2_ASM584v2_genomic.fna fastq/SRR6170103/SRR6170103_1_trim_paired.fastq fastq/SRR6170103/SRR6170103_2_trim_paired.fastq > SRR6170103.sam
```

See bwa manual [here](http://bio-bwa.sourceforge.net/bwa.shtml) for more options.

- Convert the new sam file to bam format (bam is a binary version of the sam format):
```
samtools view -b SRR6170103.sam -o SRR6170103.bam
```
- Sort the bam file
```
samtools sort SRR6170103.bam -o SRR6170103_sorted.bam
```
- Index the sorted bam for fast access
```
samtools index SRR6170103_sorted.bam 
```
- View the sorted bam file 

```
samtools view -H SRR6170103_sorted.bam 
```

## Questions

1. What is the index of the new sorted bam file?
2. How many chromosomes are present and which version of the SAM is it?


## Step 4: BAM Refinement: duplicate removal with Picard

Picard is a set of command line tools for manipulating high-throughput sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF. More information can be found [here](https://broadinstitute.github.io/picard/command-line-overview.html#Overview)

```
java -jar /home/student/APPL/PICARD/picard.2.18.0.jar MarkDuplicates INPUT=SRR6170103_sorted.bam OUTPUT=SRR6170103_final.bam METRICS_FILE=dupl_metrics.txt
```

Let's do a quick BAMQC by running samtools:
```
samtools flagstat SRR6170103_sorted.bam > SRR6170103_sorted.flagstat
samtools flagstat SRR6170103_final.bam > SRR6170103_final.flagstat
```

## Questions

1. What's the percentage of duplicated reads? (Hint: look at the file dupl_metrics.txt)
2. Do you observe any diferences between the file with and without duplicates? (Hint: look at the flagstat files)


## BAM visualization

We will use the java web start version of IGV.
Follow this link: ​ https://www.broadinstitute.org/software/igv/download
Register, and you’ll find the IGV Java Web start.​ ​ Launch IGV with 750 MB.

- Load the bam file using *E. coli*'s reference genome.

Or run IGV from the /home/student/APLL/


