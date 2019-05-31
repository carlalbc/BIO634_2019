## BIO634 Next-Generation Sequencing 2 – Transcriptomes, Variant Calling and Biological Interpretation

## September 17/18 2018


### University of Zurich
### URPP Evolution in Action
![URPP logo](Logo_URPP_kl2.png)

Stefan Wyder & Carla Bello

stefan.wyder@uzh.ch  
carla.bello@ieu.uzh.ch


## Table of Content

### Day 1


&nbsp;   | &nbsp; | &nbsp;
-------- | --- | --- 
9.30 - 9.40 | **Welcome & Introduction** | SW&CB
9.40 - 12.30 | **QC and Mapping** <br /> [Presentation](QC_MAPPING/QC_and_mapping.pdf) \| [Hands-on](https://github.com/carlalbc/BIO634_2018/blob/master/Day1_DataQC_and_mapping.md) | CB
&nbsp; | &nbsp; | &nbsp;
13.30 - 15.45 | **Variant Calling 2** <br /> [Presentation](VARIANT_CALLING/Variant_Calling2.pdf)  \| [Hands-on](VARIANT_CALLING/Exercises_Variant_Calling.md) | SW
16.00 - 17.30 | *Talk:* Dr. Martin Fischer (ETH): Detection of the genomic signature of selection [pdf](TALKS/NGS_Selection_Fischer.pdf) |



### Day 2

&nbsp;   | &nbsp; | &nbsp;
-------- | --- | --- 
9.30 - 12.00 | **RNA-seq** <br /> [Presentation](RNAseq/RNAseq.pdf) \| [Hands-on](https://github.com/carlalbc/BIO634_2018/blob/master/Day2_RNAseq.md) | CB
12.00 - 13.00 | *Talk:* Dr. Charlotte Soneson (UZH): RNA-seq analysis [pdf](TALKS/RNAseq_Soneson.pdf) |
&nbsp; | &nbsp; | &nbsp; 
14.00 - 15.00 | *Talk:* Dr. Jean-Claude Walser (ETH): RNA-seq in ecology and evolutionary biology [pdf](TALKS/RNAseq_Walser.pdf) | 
15.00 - 17.30 | **Making sense of gene lists** <br /> [Presentation](GENE_LISTS/MakingSenseOfGeneLists.pdf)  \| [Hands-on](GENE_LISTS/Exercises_MakingSenseOfGeneLists.md) | SW


## Prerequisites for the course

- Basic command line 
- Basic knowledge about NGS data structure: reads, alignments (BAM), quality scores

or attendance of 'BIO609 Introduction to Linux and Bash Scripting' and  
`BIO610 Next-Generation Sequencing 1 – Introductory Course: Assembly, Mapping, and Variant Calling


## Directory structure

```
(Unix System hierarchical tree)

|--HOME--|
      ├─APPL
           ├── IGV/...
           ├── FREEBAYES/...
           ├── PICARD/...
           ├── GATK/...
      ├─NGS2
	   ├─RNA-Seq
	   ├─VariantCalling2
        	├── Final_Results (with pre-computed results)
            ├── ref (reference files like genomes)
            └── scripts

( more folders will be added in this level during the workshop )
```


## Installation Instructions for the Virtual Machine

We will reuse the Virtual Machine (VM) of the Linux course BIO609 and NGS course BIO610. However, we need to download the data for this course.


### Download the data for this course

- Start the VM and login
- Open a terminal
- Download the file by typing `wget url` (Note to myself: Use URL shortener)
- Unzip the file: `unzip data_NGS2.zip`


### If you have *not* yet installed the VM 
- Download the virtual machine manager VirtualBox, from [virtualbox.org](https://www.virtualbox.org/). Make sure you pick the right operation system for your laptop. 
- Install VirtualBox on your machine
- Download the VM image (~4 GB) from dropfiles.uzh. Ask for the link.
- Run VirtualBox and do `File | Import Appliance` from the menu. Choose the VM image you just downloaded (file with extension .OVA). This will trigger a menu where you can change the Appliance settings. We recommend giving the VM as much memory as you can given your local machine (about 2/3 of the total memory, but between 2-4 GB). Start the import process.

Now you can start the VM by selecting it in the list and clicking on the Start button. Login and proceed with the instructions 

### Copy/paste doesn't work in the command window

To copy, select text and click the title bar and go to Edit->Copy. You can use Edit->Paste to paste.  
  
  
The usual shortcuts for copy/paste don't work in the virtual machine. For Macs use ctrl+shift+C for copying and ctrl+shift+v for pasting.  
  
You can also setup shared folders between the VM and the host system. See below or ask us how to do it.

### Adding Shared Folders

It is possible to share data (read/write/copy files) between Virtual Machine and your host operating system (e.g. Mac OS or Windows). By default no access is granted, we have to configure it: 

1. Share a folder in your operating system  
   With a running VM go to the menu of the VirtualBox software `Devices | Shared Folder Settings…` and add 1 or multiple folders
2. In Ubuntu VM, go to the terminal and type:  
`usermod -aG vboxsf student`           
 (student is our username)
3. Reboot the VM

On Ubuntu the shared folder is located in the path `/media`


## Recommended books (Practical Computing Skills)

- [Haddock & Dunn. Practical Computing for Biologists. Sinauer Associates 2011.](http://practicalcomputing.org)  
  A good book that covers the shell/command line, programming in python & bash, databases, regular expressions. 
  Suitable for self-study and as a reference book.

- [Vince Buffalo. Bioinformatics Data Skills. O'reilly 2015](http://shop.oreilly.com/product/0636920030157.do)  
  This practical book teaches the skills that scientists need for turning large sequencing datasets into reproducible and robust biological findings.
  Also covers methods on Sequence and Alignment Data. 
  More advanced than Haddock & Dunn and progresses with faster pace.


## Recommended websites

**General**  
- <http://software-carpentry.org/>  
  Scientific Computing Resources for learning bash shell, programming in python, R, …]  
- [SEQanswers](http://seqanswers.com/) the NGS community (Questions&Answers, protocols, software lists, news)   
- [BioStars](https://www.biostars.org/) for questions about biocomputing and scripting for biologists  
- [stackoverflow](http://stackoverflow.com/) for questions related to coding

**Linux/Shell**  
- [Cheatsheet intermediate](http://www.cheatography.com/davechild/cheat-sheets/linux-command-line/pdf/)  
- [SIB e-learning: UNIX fundamentals](http://edu.isb-sib.ch/pluginfile.php/2878/mod_resource/content/3/couselab-html/content.html)
