# HIVdynamite (BETA)
dynamic identification of transmission clusters in phylogenies

**This project is still in beta and as such may contain major bugs and we do not guarantee proper results/functionality**

## Requirements
R 3.6.0 or higher

**Libraries**

* phytools
* data.tree
* tidytree
* rlist
* familyR
* tidyverse
* ggplot2
* gridExtra
* ggtree
* drc
* quantmod
* remotes
* growthrates
* phangorn
* parallel
* lubridate

## Execution 

```R
source("HIVdynamite.R")
processTree("newick_tree_file")
```

Important functions (with defaults):

```R
lsd2(inputTree=sub_tree,inputDate=sts,estimateRoot="as", constraint=TRUE,variance=1, ZscoreOutlier = 3,outFile = "lsd2_results",seqLen = seqLen,nullblen=-1)
  
 skygrowth.map(tree, res = round((present_date-start_date)*12*4), tau0 = .1) 
 
 pickClust()
                 
```

## Output

HIVdynamite will result in the following output:


### Data Table
result.rds file containing classification information (e.g., growth, decay), relevant values (e.g., R0 and TMRCA), and trait distribution data over time (e.g., age, gender) for each identified cluster. This file can be opened with \url{need info from NANA analytics}.

 
### Tree
result.nexus will contain a nexus tree file with annotated information regarding cluster association, time, and traits inferred using ancestral state reconstruction.


## Notes
Reference article and describe DYNAMITE here.


![Alt text](./Images/numbering_example.png?raw=true "Numbered Tree")

This example image may be generated using:

```R
library(ape)

tree <- read.tree(text = "(((A,B),(C,D)),(E,F));")
plot(tree, edge.width = 2, label.offset = 0.1)
nodelabels()
tiplabels()
```


## References (not finished)
Paradis E, Claude J & Strimmer K (2004). APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289-290.

Yu G, Smith D, Zhu H, Guan Y, Lam TT (2017). ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods in Ecology and Evolution, 8, 28-36.

Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research, 43(7), e47. 

Revell, LJ (2012). phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol. 3 217-223.

Kembel SW, Cowan PD, Helmus MR, Cornwell WK, Morlon H, Ackerly DD, Blomberg SP, and Webb CO (2010). Picante: R tools for integrating phylogenies and ecology. Bioinformatics 26:1463-1464. 

