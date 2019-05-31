
## Variant Calling 2

### Evolution in Action
### Stefan Wyder
  
![URPP logo](../Logo_URPP_kl2.png)  
  
#### Topics
- Variant calling using FreeBayes
- Filtering variants
- Visualizing aligned reads
- Variant calling using GATK
- Variant Annotation using SnpEff & SnpSift
  

This training gives an introduction to the use of some popular NGS analysis at the command line. It presents several tools and some code, which can be used as a basis for your own workflows. However, be aware that for your own analysis different parameters than the ones used here might be better suited. Read carefully the documentation and take care that you understand what every parameter does.  

- **Sanity check your results frequently!**  
- **Use more than one SNP caller and compare their result**  


The approaches presented here work for Illumina data.


### Download data

```
wget https://www.dropbox.com/s/hgucirf8a5iys3k/DATA_NGS2_SNP.zip?dl=0
mv DATA_NGS2_SNP.zip?dl=0 DATA_NGS2_SNP.zip
unzip DATA_NGS2_SNP.zip
```


## SNP Calling using FreeBayes

Here we use FreeBayes, which is a fast easy-to-run variant caller for short variants (Garrison & Marth, http://arxiv.org/abs/1207.3907). It identifies small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), small indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment. 

### Simple Example with *E.coli*

We will identify short variants relative to the genome of the ***E.coli* DB10 strain** (Ncbi NC_010473). The experiment was done with the same strain the reference genome was made so we expect only a very small number of genetic differences. It is a relatively easy task – *E. coli* is haploid, the genome does not contain long repeats and we have a homogenous coverage across almost the whole genome. In many cases it will be much harder to identify SNPs reliably e.g. for diploid eukaryotic genomes containing lots of repetitive elements, pseudogenes and close paralogs.

FreeBayes expects a cleaned BAM file as input with marked duplicates. It calls variants for any number of individuals of a population, we can combine data from multiple individuals into a single BAM by attaching Read Group (RG) to the alignments.

How we can get information about freebayes options:
```
freebayes --help
```

The reads were already mapped to the *E.coli* DH10B genome, you will find a BAM file `MiSeq_Ecoli_DH10B_110721_PF_subsample.bam`. Let us rename the file to `Ecoli_DH10B.bam` to save typing. 


Get information about the BAM file doing:
```
samtools flagstat Ecoli_DH10B.bam
```
In our sample we observe very few anomalous read pairs, so we will skip for now the filtering of anomalous read pairs (You find a script with a more extended pipeline at the bottom of this page). 

We remove duplicate reads using samtools:
```
samtools rmdup Ecoli_DH10B.bam Ecoli_DH10B-rmdup.bam
[bam_rmdup_core] 1200 / 650578 = 0.0018 in library
```
Our file only contains 0.18% of duplicate reads.
(We could also just mark duplicate reads using Picard tools, reads marked as duplicates in the BAM file will then be ignored by freebayes).


That's it, now we have an analysis-ready BAM.


As always we have to index the BAM file:
```
samtools index Ecoli_DH10B-rmdup.bam
```
The index command creates an additional file allowing faster access. BAM index files are required for the next command.


Now we can run an analysis without any parameter optimization just setting the ploidy to 1 (this will take a few minutes):
```
freebayes -p 1 -f EcoliDH10B.fa Ecoli_DH10B-rmdup.bam > Ecoli_DH10B-rmdup.vcf
```

### Indexing vcf files

Many downstream applications only work with compressed & indexed vcf files using `bgzip` and `tabix` allowing fast access.

```
bgzip Ecoli_DH10B-rmdup.vcf   # makes Ecoli_DH10B-rmdup.vcf.gz
tabix -p vcf Ecoli_DH10B-rmdup.vcf.gz
```
As we have to do this for each vcf file, we write a script that does the job. Save the script below as `indexVCF.sh` and make it executable. It sorts, compresses and indexes a vcf file by doing: `./indexVCF.sh Ecoli_DH10B-rmdup.vcf` 
```
#!/bin/bash
# processes a vcf file for downstream analysis: sorting by chromosome and position, compressing, and indexing a vcf file
# from http://wiki.bits.vib.be/index.php/NGS_variant_analysis_custom_functions_and_code

# accepts both piped content or a vcf file as input
# writes header through and sorts remaining lines
# compresses with bgzip, indexes with tabix 
# and cleans up temp file (piped version)
    if [ $# -eq 1 ]; then
        title=$(basename "$1" ".vcf");
        ( grep ^"#" $1 | perl -pi -e "s/$title.*$/$title/";
        grep -v ^"#" $1 | sort -k 1V,1 -k 2n,2 ) | \
        	bgzip -c > $1".gz" && tabix -p vcf $1".gz";
    else
        cat > .tmpf;
        outfile="sorted-data.vcf.gz";
        ( grep ^"#" .tmpf;
        grep -v ^"#" .tmpf | sort -k 1V,1 -k 2n,2 ) | \
        	bgzip -c > ${outfile} && ( rm .tmpf;
        tabix -p vcf ${outfile} );
    fi
```

### Take a peek using `vt`

[vt](http://genome.sph.umich.edu/wiki/vt) is a toolkit for variant annotation and manipulation. In addition to other methods, it provides a nice method, vt peek, to determine basic statistics about the variants in a VCF file.

Add the directory that contain the vt binary to the PATH like this:
```
export PATH=$PATH:~/APPL/vt
```

Now we can get a summary like so:
```
vt peek Ecoli_DH10B-rmdup.vcf.gz
```

Inspect also the vcf output file and try to understand what the different columns mean. You will find a very nice poster summary under http://vcftools.sourceforge.net/VCF-poster.pdf . The detailed vcf format specifications you find [here](http://samtools.github.io/hts-specs/)

Different metrics can be used to filter the variants, the FreeBayes manual states: 
```
Of primary interest to most users is the QUAL field, which estimates the probability that there is a polymorphism at the loci described by the record. In freebayes, this value can be understood as 1 - P(locus is homozygous given the data). It is recommended that users use this value to filter their results, rather than accepting anything output by freebayes as ground truth. By default, records are output even if they have very low probability of variation, in expectation that the VCF will be filtered using tools such as vcftools.
```

- Sort the records according to quality and check some proposed variants using the IGV genome browser. Hint: you can sort the vcf file using `zgrep -v "#" File.vcf.gz | sort -k6,6gr | head`  
- Check what are the options of FreeBayes. Use the FreeBayes manual (https://github.com/ekg/freebayes) to explore and change parameters. Important parameters are `--min-mapping-quality` and `--min-base-quality`. How do different parameters influence your results?


### Filtering variants

With any variant caller you use, there will be an inherent trade-off between sensitivity and specificity. Typically, you would carry forward as much data as practical at each processing step, deferring final judgement until later so you have more information. For example, you might not want to exclude low coverage regions in one sample from your raw data because you may be able to infer a genotype based on family information later.

This typically leads to NGS pipelines that maximize sensitivity, leading to a **substantial number of false positives.** Filtering after variant calling can be useful to eliminate false positives based on all the data available after numerous analyses.

Like this we remove all low-quality variants:
```
sudo apt install bcftools # install bcftools
bcftools filter -O z -o Ecoli_DH10B-rmdup-filtered.vcf.gz -i'%QUAL>10' Ecoli_DH10B-rmdup.vcf.gz
vt peek Ecoli_DH10B-rmdup.vcf.gz
```
All but 1 variant have been filtered out as low-quality.

<details>
<summary>(Optional) An alternative graphical summary with more information using bcftools:</summary>
```
bcftools stats -F EcoliDH10B.fa -s - Ecoli_DH10B-rmdup.vcf.gz > Ecoli_DH10B-rmdup.vcf.gz.stats
sudo apt-get install python-matplotlib
sudo apt-get install texlive-latex-base exlive-latex-recommended exlive-latex-extra
plot-vcfstats -p bcftools_plots/ Ecoli_DH10B-rmdup.vcf.gz.stats
```
</details>


### Visualize the aligned reads

Seeing is believing! One should always have at look at the data to get a feeling about the error rate, coverage heterogeneity, ...  

Go to the Integrative Genome Viewer (IGV) website http://www.broadinstitute.org/igv/ 

1. Open a new tab in the terminal (File | Open Tab)
1. Go to IGV directory doing `cd ~/APPL/IGV/IGV_2.4.10/`
2. Launch IGV doing `./igv.sh`
3. Load the fasta file of the genome: File | Load Genome from File...  
  then choose the file `EcoliDH10B.fa`
4. Load the BAM file: File | Load from File...   
  then choose the BAM file
5. Load the genome annotation (gff or bed): File | Load from File...  
  then choose the file `EcoliDH10B.gff


- Do you think some records are SNPs? Or rather sequencing errors? Look up the base call qualities at the varying positions by moving the cursor to a nucleotide.


## BAM file preprocessing

BAM File preprocessing often improves SNP calling accuracy. Possible steps include (For details check the GATK documentation):

1. Filtering anomalous read pairs (e.g. with unexpected insert size, read pairs where only 1 read could be mapped, ...)
2. Mark Duplicate Removal
3. Base Recalibration

Picard, samtools and bamtools are the most widely used tools for BAM file preprocessing.

Depending on the SNP Caller we have to perform all or just a subset of preprocessing steps. Above we have just done steps 1 + 2 before having used freebayes. Luckily, freebayes does the indel realignment step internally and does not require base recalibration. We just have to do steps 1 + 2.


In contrast, for GATK it is recommended to use all 3 BAM preprocessing steps. Let's look at the GATK workflow.


### GATK


- GATK is powerful, but it requires relatively complex pipelines 
- GATK has been developed for human data and they provide files with known SNPs, indels etc. only for human. Some of the preprocessing steps can not be run for other organisms.
- Even if you do not work with human data, you may find valuable ideas by looking over mature pipelines from large genome projects.
- Notably base quality recalibration can be done for any organism with a workaround. Check this [publication] (http://www.g3journal.org/content/5/4/655.full)
- GATK is slow and resource-hungry. If you want to run it for pooled or population data, prepare to run in on a cluster.
- GATK is under constant development. Check the website from time to time.


### Installing Oracle Java

GATK doesn't work with JDK Java - it requires Oracle Java
```
sudo add-apt-repository ppa:webupd8team/java
sudo apt-get install oracle-java8-installer
sudo update-alternatives --config java
Select Oracle Java from the list
```


The following script `Run_GATK_Ecoli.sh` runs a GATK pipeline for the *E.coli* BAM file. Change the variable BAMFILE if you have renamed the BAM file earlier.
```
#!/bin/bash
set -o nounset
set -o errexit

 # BAM processing using PICARD, SNP calling using GATK

#SAMTOOLS=~/software/SAMTOOLS/samtools-1.3/samtools
SAMTOOLS=samtools
GATK=~/APPL/GATK/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
PICARD=~/APPL/PICARD/picard.2.18.0.jar

 # IDEAS to improve the script:
 # add variable with memory parameter to java -Xmx6g
 # add 2> to capture errors
 # add input and output folders


REF_FILE=EcoliDH10B.fa
BAM_FILE=MiSeq_Ecoli_DH10B_110721_PF_subsample.bam
 # outfolder=gatk_variants

 # build indexes for Genome and BAM files
$SAMTOOLS faidx $REF_FILE
$SAMTOOLS index $BAM_FILE
java -jar $PICARD CreateSequenceDictionary R=$REF_FILE O=${REF_FILE%%.*}.dict

 # sort mappings by name to keep paired reads together
java -Xmx1g -jar $PICARD SortSam \
	I=$BAM_FILE \
	O=${BAM_FILE%%.*}_sorted.bam \
	SO=queryname \
	VALIDATION_STRINGENCY=LENIENT \
	2>SortSam_queryname.err

 # fix mate information and sort
java -Xmx1g -jar $PICARD FixMateInformation \
	I=${BAM_FILE%%.*}_sorted.bam \
	O=${BAM_FILE%%.*}_fixmate.bam \
	SO=coordinate \
	VALIDATION_STRINGENCY=LENIENT \
	2>FixMateInformation.err


 # mark duplicates
java -Xmx1g -jar $PICARD MarkDuplicates \
	I=${BAM_FILE%%.*}_fixmate.bam \
	O=${BAM_FILE%%.*}_rdup.bam \
	M=duplicate_metrics.txt \
	VALIDATION_STRINGENCY=LENIENT \
	#2>MarkDuplicates.err


 # Add Read Groups @RG
java -Xmx1g -jar $PICARD AddOrReplaceReadGroups \
	I=${BAM_FILE%%.*}_rdup.bam \
	O=${BAM_FILE%%.*}_rdup-rg.bam \
	RGID="NA18507" \
	RGLB="lib-NA18507" \
	RGPL="ILLUMINA" \
	RGPU="unkn-0.0" \
	RGSM="MiSeq_Ecoli_DH10B" \
	VALIDATION_STRINGENCY=LENIENT \
	# 2>AddOrReplaceReadGroups.err

 # Build BAM index for fast access
java -Xmx1g -jar $PICARD BuildBamIndex \
	I=${BAM_FILE%%.*}_rdup-rg.bam \
	2>BuildBamIndex.err


 # identify regions for indel local realignment of the selected chromosome
java -Xmx1g -jar $GATK -T RealignerTargetCreator \
	-R $REF_FILE \
	-I ${BAM_FILE%%.*}_rdup-rg.bam \
	-o target_intervals.list

 # perform indel local realignment of the selected chromosome
java -Xmx1g -jar $GATK -T IndelRealigner \
	-R $REF_FILE \
	-I ${BAM_FILE%%.*}_rdup-rg.bam \
	-targetIntervals target_intervals.list \
	-o ${BAM_FILE%%.bam}_realigned.bam

 # analyze patterns of covariation
 # only possible with known SNPs 
 # for non-model organisms it is possiple to interpret the high quality fraction of called SNPs as knownSites and perform multiple rounds of BaseRecalibration and SNPcalling
 #java -Xmx1g -jar $GATK -T BaseRecalibrator -R $REF_FILE -I ${BAM_FILE%%.bam}_realigned.bam
 #	-knownSites ${dbsnp} -knownSites ${gold_indels} -o recal_data.table 2>BaseRecalibrator.err

 #1. Run UnifiedGenotyper
 # this takes approx. 1.8 min
java -Xmx1g -jar $GATK -T UnifiedGenotyper \
	-R $REF_FILE \
	-I ${BAM_FILE%%.bam}_realigned.bam \
	-ploidy 1 \
	-glm BOTH \
	-stand_call_conf 30 \
	-stand_emit_conf 10 \
	-mbq 10 \
	-o raw_variants_UG.vcf
 
 #2. Run HaplotypeCaller
java -Xmx1g -jar $GATK -T HaplotypeCaller \
	-R $REF_FILE \
	-I ${BAM_FILE%%.bam}_realigned.bam \
	--genotyping_mode DISCOVERY \
	-ploidy 1 \
	-stand_call_conf 30 \
	-stand_emit_conf 10 \
	-o raw_variants_HC.vcf
	# -L 20
```

- Try to understand all the steps.
- We run both UnifiedGenotyper and HaplotypeCaller. What is the difference?
- A good page how to interpret vcf is [here](https://www.broadinstitute.org/gatk/guide/article?id=1268)  
 
#### GATK with a human sample

Next we will work on the low coverage data set of the [1000 genomes project](http://www.1000genomes.org/). The human sample (HG00154) was sequenced with ILLUMINA and the reads were aligned using bwa ([see details](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/README.alignment_data)) against the human reference genome hg19.

We will download only a subset of the original data (to safe disk space and execution time), 1 MB of chromosome 20.

- First create a new folder and download this human sample (The download may take some minutes - in the meanwhile, open a new tab of your terminal and proceed with the next steps). Download the sequence of chromosome 20 and the bam file from Dropbox:
```
wget https://www.dropbox.com/s/rc7rdbm764mm89s/human_example.zip
```
 
The BAM file was prepared for you like this (you don't have to do it): 
  ```
  samtools view -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00154/alignment/HG00154.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam 20:1000000-2000000 > HG00154.low_coverage.chr20.1000000-2000000.sam
samtools view -bS HG00154.low_coverage.chr20.1000000-2000000.sam > HG00154.low_coverage.chr20.1000000-2000000.bam
```
    
- Index the BAM file
- Copy the script, adapt it (BAM_FILE, REF_FILE, ploidy) and run GATK on it
- Compress & index the resulting vcf files
- Filter out low-quality variants and have a look at the filtered file:
```
bcftools filter -O z -o VCF.vcf-filtered.gz -i'%QUAL>10' VCF.vcf.gz
vt peek VCF.vcf-filtered.gz
``` 
- Check some variants in IGV



## Variant Annotation using SnpEff & SnpSift

There are many software tools available for annotating variants. Many tools only work with Human and mouse variants where  lots of functional annotation are available. [SnpEff & SnpSift](http://snpeff.sourceforge.net/SnpEff.html) is widely used and is also useful for work with other organisms. SnpEff provides functions to annotate variants and predict the effects of variants on genes (such as amino acid changes), while SnpSift allows you to filter and manipulate annotated files.


First download snfEff via wget (or alternatively, from the [download](http://snpeff.sourceforge.net/download.html) page):
```
 # Download using wget
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
 # Install
$ unzip snpEff_latest_core.zip 
```
Install it in the `~/APPL` directory.
 

### Download databases

As of May 2015 it comprises > 35'000 databases for many organisms (including plants, bacteria) and genome assembly versions.

Here we check how many databases are available for the Arabidopsis genus. 

```
java -jar snpEff.jar databases | grep Arabidopsis
TAIR10.26                                                   	Arabidopsis_thaliana                                        	          	                              	http://downloads.sourceforge.net/project/snpeff/databases/v4_1/snpEff_v4_1_TAIR10.26.zip
alyrata1                                                    …
...                                              	Arabidopsis_lyrata                                          	          	                              	http://downloads.sourceforge.net/project/snpeff/databases/v4_1/snpEff_v4_1_v.1.0.26.zip
```

We see that there are databases available for A. thaliana (both TAIR9 and TAIR10 annotation) and A. lyrata. 
SnpEff can also help you to build a database from reference genome files (Fasta, GFF, etc.).

We can download and install a selected database like:
```
java -jar snpEff.jar download -v athalianaTair10
```

### Annotate using SnpEff

```
java -Xmx4g -jar snpEff.jar -stats summary.html athalianaTair10 test.Ath.vcf > test.Ath.ann.vcf
```

```
##SnpEffVersion="4.1e (build 2015-05-02), by Pablo Cingolani"
##SnpEffCmd="SnpEff  athalianaTair10 test.Ath.vcf "
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
##INFO=<ID=LOF,Number=.,Type=String,Description="Predicted loss of function effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected' ">
##INFO=<ID=NMD,Number=.,Type=String,Description="Predicted nonsense mediated decay effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected' ">
Chr1	117511	.	C	T	.	.	28,27,14,82,;ANN=T|synonymous_variant|LOW|AT1G01300|AT1G01300|transcript|AT1G01300.1|Coding|1/1|c.447C>T|p.Gly149Gly|569/1822|447/1458|149/485||,T|synonymous_variant|LOW|Exon_1_117065_118522|Gene_AT1G01300.1-Protein|transcript|AT1G01300.1-Protein|Coding|1/1|c.447C>T|p.Gly149Gly|447/1458|447/1458|149/485||,T|upstream_gene_variant|MODIFIER|AT1G01305|AT1G01305|transcript|AT1G01305.1|Coding||c.-33C>T|||||1886|,T|upstream_gene_variant|MODIFIER|Exon_1_119429_119854|Gene_AT1G01305.1-Protein|transcript|AT1G01305.1-Protein|Coding||c.-1C>T|||||1918|,T|upstream_gene_variant|MODIFIER|AT1G01310|AT1G01310|transcript|AT1G01310.1|Coding||c.-68C>T|||||2643|,T|upstream_gene_variant|MODIFIER|Exon_1_120221_120946|Gene_AT1G01310.1-Protein|transcript|AT1G01310.1-Protein|Coding||c.-1C>T|||||2710|,T|downstream_gene_variant|MODIFIER|AT1G01280|AT1G01280|transcript|AT1G01280.1|Coding||c.*1560C>T|||||3564|,T|downstream_gene_variant|MODIFIER|Exon_1_112290_113195|Gene_AT1G01280.1-Protein|transcript|AT1G01280.1-Protein|Coding||c.*1533C>T|||||3606|,T|downstream_gene_variant|MODIFIER|AT1G01290|AT1G01290|transcript|AT1G01290.1|Coding||c.*826C>T|||||1962|,T|downstream_gene_variant|MODIFIER|Exon_1_114299_114433|Gene_AT1G01290.1-Protein|transcript|AT1G01290.1-Protein|Coding||c.*813C>T|||||2215|,T|downstream_gene_variant|MODIFIER|Exon_1_114299_114433|Gene_AT1G01290.2-Protein|transcript|AT1G01290.2-Protein|Coding||c.*813C>T|||||2215|,T|downstream_gene_variant|MODIFIER|AT1G01320|AT1G01320|transcript|AT1G01320.1|Coding||c.*5643G>A|||||3613|,T|downstream_gene_variant|MODIFIER|Exon_1_130039_130099|Gene_AT1G01320.1-Protein|transcript|AT1G01320.1-Protein|Coding||c.*5394G>A|||||4071|,T|downstream_gene_variant|MODIFIER|Exon_1_130039_130099|Gene_AT1G01320.2-Protein|transcript|AT1G01320.2-Protein|Coding||c.*5364G>A|||||4071|,T|downstream_gene_variant|MODIFIER|AT1G01290|AT1G01290|transcript|AT1G01290.2|Coding||c.*826C>T|||||1403|,T|downstream_gene_variant|MODIFIER|AT1G01320|AT1G01320|transcript|AT1G01320.2|Coding||c.*5620G>A|||||3556|
```
See the documentation for more options and functions.


### Exercise: VCFtools

VCFtools are specialized tools for working with VCF files: validating, filtering merging, comparing and calculate some basic population genetic statistics (see the documentation http://vcftools.sourceforge.net/docs.html). There are many other tools with similar functionality available for the same purpose.
￼￼￼￼￼￼
### Exercise: How many low-coverage regions?

Often some regions of the genome are low coverage only (or even without any aligned reads) consequently we cannot tell whether polymorphisms exist in these regions. We want to identify such regions and find out whether they overlap with genes.
The bedtools utilities are convenient for working with genomic coordinates for example bedtools allows one to intersect, merge, count, complement, and shuffle genomic intervals from multiple files in widely-used genomic file formats such as BAM, BED, GFF/GTF, VCF. You will find the documentation under http://bedtools.readthedocs.org/en/latest/index.html


Use `MiSeq_Ecoli_DH10B_110721_PF_subsample.bam` from the first exercise.  

- Try to find out how many nucleotides in the *E. coli* genome do not reach a minimal read coverage of 5, 10 or 20 reads (Hint: Use samtools depth, direct the output into a file and then use awk or R to process the file)
- Do some low coverage nucleotides overlap with annotated genes? (Hint: use the command intersect from bedtools)
- We now counted individual nucleotides but actually we would like to know whether there are some regions with low coverage(Hint: use mergeBed from bedtools)
- Are there any regions with unexpected high coverage? As mapping problems are likely (due to repetitive regions) they are often excluded from the analysis.

#### Solution
```
samtools depth MiSeq_Ecoli_DH10B_110721_PF_subsample.bam | awk 'BEGIN{OFS="\t"} $3<10 {print $1,$2-1,$2}' | bedtools merge -i - | head
EcoliDH10B.fa	0	17
EcoliDH10B.fa	15618	15749
EcoliDH10B.fa	16379	16496
EcoliDH10B.fa	19963	19964
EcoliDH10B.fa	19981	20114
```

## Mixed *E.coli* population

This exercise was taken from a NGS course at the University of Texas at Austin by Jeff Barrick https://wikis.utexas.edu/pages/viewpage.action?pageId=33949671

**Example 1**
A mixed population of E. coli from an evolution experiment was sequenced at several different time points ([PMID: 19838166](http://www.ncbi.nlm.nih.gov/pubmed/19838166) , [PMID:19776167](http://www.ncbi.nlm.nih.gov/pubmed/19776167)). At generation 0 it consisted of a clone (cells grown from a colony with essentially no genetic variation), then additional samples were taken after 20K and 40K generations (roughly 10 and 20 years) of culture in the laboratory. By these later times new mutations had arisen. Some had taken over (swept to fixation in) the population and were present in 100% of the individuals. Some mutations were only present in certain subpopulations and thus they were are intermediate frequencies (like 25%) in the sample.

### Data

File Name | Description
SRR030252.fastq.gz | Illumina reads, 0K generation individual clone from population
SRR032374.fastq.gz | Illumina reads, 20K generation mixed population
SRR032376.fastq.gz | Illumina reads, 40K generation mixed population
NC_012967.1.fasta.gz | E. coli B str. REL606 genome



The read files were downloaded from the [ENA SRA study](http://www.ebi.ac.uk/ena/data/view/SRP001569).
So that you can treat all the data as single-ended for simplicity, we concatenated two separate FASTQ (paired-end) files for sample SRR030252 using this command
```
cat SRR030252_1.fastq SRR030252_2.fastq > SRR030252.fastq
```
Alternatively, you could map that data set as paired-end.

The reference genome file was downloaded from the [NCBI Genomes page](http://www.ncbi.nlm.nih.gov/genome/167?project_id=58803).
We renamed the FASTA sequence from "gi|254160123|ref|NC_012967.1|" to "NC_012967" by changing the first line of the NC_012967.1.fasta sequence using a text editor. It's just easier to deal with the shorter name.

### Map Reads

Choose an appropriate program and map the reads. As for other variant callers, convert the mapped reads to BAM format, then sort and index the BAM file.

### Additional exercises

- Determine the approximate depth of mapped read coverage for each sequencing data set.
- Try using different mappers or changing the default alignment settings to find more variants.


### Run FreeBayes

FreeBayes can be used to treat the sample as a mixture of pooled samples. (In our case it is actually a mixture of >1 million bacteria, but we have nowhere near that level of coverage, so we give an arbitrary mixed ploidy of 100, which means we use a statistical model that predicts variants only with frequencies of 1%, 2%, 3%, ... 98%, 99%, 100%). This command runs pretty fast, so you can do it in interactive mode.
Example command for running freebayes
```
$ freebayes --min-alternate-total 3 --ploidy 100 --pooled --vcf SRR032374.vcf \
        --fasta-reference NC_012967.1.fasta SRR032374.sorted.bam
```
Nice! We now have a familiar VCF file.

#### Additional exercises

- Write a script or use a linux command to filter the output files to only contain variants that are predicted to have frequencies > 0.10 and < 0.90 or scores > 100.
- What does the `--min-alternate-total` option mean? Experiment with using other options and examine how they affect speed versus what is predicted. Beware that some choices of options can cause crashes.
- Compare the variants predicted in samples SRR032374 and SRR032376. You could view them in IGV.

### Sources & Links

- [Tutorial from the creator of freebayes](https://github.com/ekg/alignment-and-variant-calling-tutorial)
- [Extensive course on variant analysis](http://wiki.bits.vib.be/index.php/Hands-on_introduction_to_NGS_variant_analysis)

### More information

- [Detailed GATK Presentations](https://www.broadinstitute.org/gatk/guide/presentations?id=4768) (2-day workshop) and [Presentation PDFs](https://www.dropbox.com/sh/uw2cor3eix3w47e/AADCAkp7l51Wz9uCIjUboyU2a?dl=0)
- Beware format differences: 1-based VCF format and 0-based BED!
  The UCSC Genome Browser has nice short [descriptions](http://genome.ucsc.edu/FAQ/FAQformat#format10.1) of widely used formats.
 
 
### Example Scripts
 
#### SNP calling using freebayes with BAM preprocessing (fix mate pairs, mark duplicates)

```
#SAMTOOLS=~/software/SAMTOOLS/samtools-1.3
SAMTOOLS=samtools
PICARD=~/APPL/PICARD/picard.2.18.0.jar
FREEBAYES=~/APPL/FREEBAYES/freebayes/bin/freebayes

memory="1G"
threads=1

cp -p MiSeq_Ecoli_DH10B_110721_PF_subsample.bam Ecoli_DH10B.bam

 # sort reads by name required by the next command
$SAMTOOLS sort -n -m $memory -@ $threads -O bam -T /tmp/sort1 -o Ecoli_DH10B-sorted.bam Ecoli_DH10B.bam

 # fix mate information and output to bam
$SAMTOOLS fixmate -O bam Ecoli_DH10B-sorted.bam Ecoli_DH10B-fixmate.bam

 # sort the bam data by position
$SAMTOOLS sort -m $memory -@ $threads -O bam -T /tmp/sort2 -o Ecoli_DH10B-pos-sorted.bam Ecoli_DH10B-fixmate.bam

 # mark duplicates
 # Note that we use samtools version 0.1.18 here as this function is not yet implemented in samtools v1.2
samtools rmdup Ecoli_DH10B-pos-sorted.bam Ecoli_DH10B-fixmate-mdup.bam
 #$SAMTOOLS rmdup Ecoli_DH10B-pos-sorted.bam Ecoli_DH10B-fixmate-mdup.bam

 # index BAM file$
$SAMTOOLS index Ecoli_DH10B-fixmate-mdup.bam

 # call variants
$FREEBAYES -f EcoliDH10B.fa -p 1 Ecoli_DH10B-fixmate-mdup.bam > Ecoli_DH10B-fixmate-mdup.vcf
```


### Install software



#### Install vt

Go to the `~/APPL` folder and build & install vt:

```
git clone https://github.com/atks/vt
cd vt
make
```

Export the path to vt like this:
```
export PATH=~/APPL/vt:$PATH
```

#### Add directory to the PATH permanently

An export command works only temporarily
```
export PATH=~/APPL/vt:$PATH
```

To make it permanent we have to add the line to the file `~/.bashrc`

In order to immediately reflect changes execute `source ~/.bashrc` 

#### An alternative graphical summary with more information using bcftools:
```
bcftools stats -F EcoliDH10B.fa -s - Ecoli_DH10B-rmdup.vcf.gz > Ecoli_DH10B-rmdup.vcf.gz.stats
sudo apt-get install python-matplotlib
sudo apt-get install texlive-latex-base
plot-vcfstats -p bcftools_plots/ Ecoli_DH10B-rmdup.vcf.gz.stats
```
