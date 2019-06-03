# Making sense of gene lists

### URPP Evolution in Action
### Stefan Wyder


We will practise the 3 different ways of making sense of a gene list we have discussed in the presentation.



## Exercise 1: Over-Representation / Enrichment Analysis

Over-representation analysis is simply a 2x2 contigency table separately for each gene. 

There are many online tools available for enrichment analysis, but many tools are not recommendable. Before you use one, make sure the tool  
- lets you define the background (i.e. all measured genes)  
- uses an up-to-date annotation  
- does multiple testing correction (e.g. FDR)  

A good choice is e.g. [GOrilla](http://cbl-gorilla.cs.technion.ac.il/) for 6 different species.  
  
In R/Bioconductor, I like the [topGO](http://www.bioconductor.org/packages/release/bioc/html/topGO.html) package. 
 
**Some tips**  
Generally ignore categories  
- with very few genes  
- with large overlap with other categories 

Often you get back a huge list of significantly enriched terms which is time-consuming to digest. [REViGO](http://revigo.irb.hr/) 
([PMID:21789182](http://www.ncbi.nlm.nih.gov/pubmed/21789182)) takes long lists of Gene Ontology terms and summarizes them by
removing redundant GO terms. It also makes summary plots.


### Exercise

Try out GOrilla with the list of genes overexpressed in Adult mice relative to Newborns.  
  
First we have to prepare the input files (1 gene name per line). We take the 734 genes overexpressed in Adult relative to Newborn mice
(from the RNA-seq tutorial) with a FDR < 5%. 
```
awk '$5<1 && $8<=0.05 {print $1}' diffExp_N-A.txt  | sed 's/"//g' > diffExp_N_smallerThan_A.FDR5perc.txt
awk 'NR>1 {print $1}' diffExp_N-A.txt  | sed 's/"//g' > diffExp_N-A.allGenes.txt
```

Run an over-representation analyis using [GOrilla](http://cbl-gorilla.cs.technion.ac.il/). Use the basic functionality <Two unranked lists of genes (target and background lists)> and upload the 2 files we just prepared. In step 4 choose <Process> and run the analysis.

Look at the table output below the figure. Click at some significantly enriched processes for the full description. It is well visible from the graph the parents of highly significantly GO categories often also become significant (Parents contain all the genes of their children).  
  
Try out GOrilla with the list of differentially expressed genes from the RNA-seq session.


## Gene Set Enrichment Analysis (GSEA)

Gene Set Enrichment Analysis (GSEA) is a more sensitive and robust alternative to Over-Representation Analysis. It does not require a (somewhat
arbitrary) significance cut-off but it uses the ranked list of all genes.
We want to find out using GSEA which pathways are enriched in the Newborn - Adult comparison from the previous exercise. 
Here we will be using the original GSEA software that has a nice Graphical User Interface providing its output as simple webpages. 
However, many derivatives of the GSEA approach are now available in Bioconductor with more precise statistics (check http://www.bioconductor.org/packages/release/BiocViews.html#___Pathways).

### Exercise 

First we need to provide a **list of ranked genes**. Here we take the probabilities for differential expression between the 2 groups, 
Newborn vs Adults, and we separate the 2 possible directions by multiplying FDR with -1 for genes which are underexpressed in 
Newborns relative to Adults (we can take either p-value or FDR which give the same result).
```
awk '$5<1 {factor=-1} NR>1 {print $1"\t"$8*factor;factor=1}' diffExp_N-A.txt | sed 's/"//g' > diffExp_N-A.4GSEA.rnk 
```
We also have to provide *gene sets*. The GSEA gene set files are available for human only therefore I have prepared a file 
c2.cp.v2.5.symbols_MmProjection_ensembl.gmt which projects the human pathway annotation to mouse genes via orthology.

Now we are ready to run GSEA. Launch GSEA in the virtual machine:
```
java -Xmx1000m -jar /software/GSEA/gsea2-2.0.14.jar
```

Next we load the data. 
- Click on the on the left <Load data> and <Browse for files...> then select the 2 required files.
- From the menu select <Tools | GseaPreranked>.
- Fill out the required fields: For Gene Sets Database click on <...> and press the right arrow until you reach <Gene matrix local gmx/gmt>
then click on the file and press <OK> 
- Set <Collapse dataset to gene symbols> to <false> 
- Under <Basic Fields> and <Advanced Fields> one can change the analysis parameters. We don't change anything.
- At the bottom of the window change to <Normal> CPU Usage.
- Run the analysis by pressing the <Run> button at the lower bottom of the window
- Once finished, click on the green <Success> in the lower left corner of GSEA. A webpage will then open in your browser.  


Go to <Detailed enrichment results in html format> and on the following page go to the <Details> for some pathways. 
Try to understand the different columns.  
Documentation of the GSEA method (User guide, tutorials, file formats etc) is available [here](http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Main_Page)  
  
Alternatively, you can access the GSEA output from the web browser (the folder containing the date might be called differently - use the <Tab> autocompletion!)
```
firefox /home/student/gsea_home/output/may18/my_analysis.GseaPreranked.1401102143051/index.html
```

**Interpretation**  
- You can look up the pathway description on [MSigDB](http://www.broadinstitute.org/gsea/msigdb/index.jsp) (Registration required).  
- Explore the most significant pathways.  
- Are they all plausible?  
- Pathway analysis is easier to interpret and often gives better biological insight than Gene Ontology (but a larger fraction of genes is annotated with Gene Ontology terms than with patways.


**Non-model organisms**  
The GSEA websites provides only annotation files for human, but you can prepare gene set (.gmt) files for any species from Gene Ontology,or Pathway annotation. Check the [GSEA documentation (http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Main_Page) how to format gene set (.gmt) files.

Check the scope of available pathways in KEGG on http://www.genome.jp/kegg/pathway.html



### Network analysis with STRING

As an example for network analysis we try out STRING.  
  
STRING stands for *S*earch *T*ool for the *R*etrieval of *I*nteracting *G*enes/Proteins. It is a webtool and database for retrieving networks of interacting genes. The interactions include direct (physical) and indirect (functional) associations. STRING integrates previous knowledge from known experimental interactions and pathway databases, it does automated text-mining of the literature but it also covers the less well understood portion of gene interactions as predicted *de novo* by a number of algorithms using genomic information. Importantly, each association has a probabilistic confidence score.  
 
The current version v10.5 comprises more than 2'000 species while interaction information is transferred between species.


### Exercise

The gene list we use here is a list of significantly up-regulated genes in the *met1* background in *Arabidopsis thaliana*
defective in methylation maintenance (PMID:18423832).  
(obtained from Table S2 from the publication [PMID:23146178](http://onlinelibrary.wiley.com/doi/10.1111/tpj.12070/abstract)
The gene list is not significantly enriched in any Gene Ontology terms or in a KEGG pathway.

------------

1. Open http://string-db.org in your browser
2. STRING accepts different types of input, often we want to lookup a gene list, but sometimes only a single gene.
  Here we have a gene list so click on `multiple proteins` and copy&paste the gene list into the search field.
3. Type Arabidopsis in the `organism` and choose `Arabidopsis thaliana`
  Press the `GO!` button.
4. Next you see a page with suggested matches. Press `Continue ->` to accept all suggestions.
5. Now we get to a page with the result network. Each dot is a gene/protein of the input list. The lines connecting genes denote
known or predicted interactions with different line colors representing the types of evidence for the association. 
6. Click on a gene to get more information about it (and links out to other databases). 
  Hover and click on some lines and find out where the evidence comes from. Right-Click on <Show> to see the origin.
  Click on some green lines and right-click on <Show> to show the relevant publications co-mentioning the gene pair.  
7. If you scroll down further you see <Your input> - the input gene(s) with descriptions.
  Below the network there is a <Settings> tab where you can switch on/off individual prediction methods and change
the required minimum confidence score.
  Change the confidence, press <Update> and check the effect on the resulting networks.
8. Below the network is a taskbar with icons `+` and `-` Here you can grow the network by 10 more partners. By pressing the `+ (more)` button 10 genes not present in your list are added to the network so that
your input genes are maximally connected. These newly added genes (in white color) are candidates for working in the same process/pathway
as your input genes (Find more information about them in the <Predicted Functional Partners:> window at the bottom of the page) 
  Press the `+ (more)` button multiple times and see how the network grows (the same functionality you can get from the <Settings> tab: `max number of interactors to show`).  
9. In the tab <Analysis>, you can see enrichment analysis with Gene Ontology and KEGG pathways.   
10. In the tab <Clusters>, you can apply different cluster algorithms on your network.    

  Make an Enrichment analysis with <GO Biological Processes> and <KEGG pathways>
  Remove genes not connected to any other genes by going to <Settings> | <hide disconnected nodes in the network>



## Tips

- Often you get back a huge list of significantly enriched terms which is difficult to digest. [REViGO](http://revigo.irb.hr/) ([PMID:21789182](http://www.ncbi.nlm.nih.gov/pubmed/21789182)) takes long lists of Gene Ontology terms and summarizes them by removing redundant GO terms. It also makes summary plots.
 
- Interesting tool seen during the preparation of the course: https://github.com/tanghaibao/goatools
