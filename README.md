## Description:
This repository contains scripts and sample code for creating and analyzing gene-level networks of phage diversity, as in our recent preprint, "Gene networks provide a high-resolution view of bacteriophage ecology," available via [biorxiv](http://www.biorxiv.org/content/early/2017/06/11/148668).  
A detailed guide (Walkthrough.txt) is available for download with instructions for going from genome sequences for a collection of phage genomes to building networks and generating host predictions for novel genomes. It can be used to recreate the analysis in our manuscript and can also be applied to new data. Please refer to the folders "R Scripts" for necessary code and to "Output Data" for examples of output. Necessary input files are contained in "Input Data."    

## Required Software 
This code is intended to be run in a UNIX environment and requires:  
  -[R](https://cran.r-project.org/) (and the packages 'Matrix' and 'igraph'): Required for most analysis.  
  -[usearch] (http://drive5.com/usearch/): Used to cluster homologous genes (Note: only 32-bit usearch is free. Larger projects might not be possible without buying a 64-bit license.)
  -[MCL] (https://micans.org/mcl/index.html): used to identify graphical clusters in networks.  
  -[Cytoscape] (http://www.cytoscape.org): used to visualize networks.    

### *Contact information*
Please contact Jason Shapiro at jshapiro2@luc.edu with any questions.  
 
