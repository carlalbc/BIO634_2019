## BIO634 Next-Generation Sequencing 2 – Transcriptomes, Variant Calling and Biological Interpretation

## June 3-4th 2019


### University of Zurich
### URPP Evolution in Action
![URPP logo](Logo_URPP_kl2.png)

Carla Bello & Gregor Rot
 
carla.bello@ieu.uzh.ch
gregor.rot@uzh.ch 


## Table of Content

### Day 1
&nbsp; | &nbsp; | &nbsp;
-------- | --- | --- 
9.30 - 9.40 | **Welcome & Introduction** | CB&GR
9.40 - 11:00 | **QC and Mapping** <br /> [Presentation](https://github.com/carlalbc/BIO634_2019/blob/master/Day1_QC_and_mapping/Day1_DataQC_and_BWAmapping.pdf) \| [Hands-on](https://github.com/carlalbc/BIO634_2019/blob/master/Day1_QC_and_mapping/Day1_DataQC_and_mapping.md) | CB
11:00 - 11.20 | *Coffee break*
11.20 - 12.30 | **QC and Mapping: Continuation** | CB
12.30 - 13.30 | *Lunch at the cafeteria*
13.30 - 15.45 | **Variant Calling 2** <br /> [Presentation](variant_calling/variant_calling_presentation.pdf)  \| [Hands-on](variant_calling/variant_calling_exercises.md) | GR
15.45 - 16.00 | *Coffee break*
16.00 - 17.30 | *Talk:* Dr. Jean-Claude Walser (ETH): RNA-seq in ecology and evolutionary biology [pdf](https://github.com/carlalbc/BIO634_2019/blob/master/UniZH_Bio634_JCW_190603.pdf)



### Day 2
&nbsp; | &nbsp; | &nbsp;
-------- | --- | --- 
9.30 - 10.45 | **RNA-seq** <br /> [Presentation](https://github.com/carlalbc/BIO634_2019/blob/master/Day2_RNAseq/Day2_RNAseq%20.pdf) \| [Hands-on](https://github.com/carlalbc/BIO634_2019/blob/master/Day2_RNAseq/Day2_RNAseq.md) | CB
10.45 - 11:00 | *Coffee break*
11.00 - 12.00 |  **Continuation: RNA-seq**  | CB
12.00 - 13.00 | *Lunch at the cafeteria* 
13:00 - 14:30 |  **Making sense of gene lists** <br /> [Presentation](gene_lists/gene_lists_presentation.pdf)  \| [Hands-on](gene_lists/gene_lists_exercises.md) | GR
14.30 - 14.50| *Coffee break* |
14.50 - 15.45 | **Making sense of gene lists: Continuation**  | GR
16.00 - 17.30 | **Talk** TBD


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
