# DYNAMITE (BETA)
DYNAMic Identification of Transmission Epicenters within phylogenies

**This project is still in beta and as such may contain major bugs and we do not guarantee proper results/functionality**

## Requirements
R 3.6.0 or higher

**Libraries**

* optparse
* remotes
* phytools
* data.tree
* tidytree
* rlist
* familyR
* tidyverse
* ggtree
* parallel
* geiger
* tibble
* skygrowth
* treedater

**File format requirements**

DYNAMITE takes a tree (scaled in substitutions/site) in newick or nexus format with support values at the nodes. Support values can be provided by any (single) method - i.e., combined result of 2 or more methods will not be accepted. These values can be scaled from 0-1 or 0-100.

Additionally, a metadata file with sampling dates and traits of interest is required. Column headers are needed. The first column must be named "id" and correspond to the taxa names, the second column must be named "date" and correspond to date information (must be consistently either numeric or Date format). Remaining columns can consist of any number of traits of interest. Specific header names are not required, as only the order is called.

## Execution 

```R
--help (-h) = provides helpful information for running DYNAMITE
--tree (-t) = path to treefile [default= .nwk file]
--metadata (-m) = path to metadata file [default= .tab file]
--seqLen (-s)  = sequence length [default=10000]
--cluster (-c) = choice of cluster algorithm from 'c' (Phylopart's cladewise) or 'b' (DYNAMITE's branchwise) [default= b]
--range (-r) = range of branch length threshold quantiles used to determine the optimal cluster branch length threshold [default=30]
--asr (-asr) = option of ancestral state reconstruction for each cluster [default= N]
```

**Important functions (with defaults):**

```R
treedater::dater(tree, dates, s=seqLen, ncpu=numCores, omega0=8E-04)

branchwise()

phylopart()

ancStateRecon(ace(model="ER", type="discrete"))
                 
```

## Output

DYNAMITE will result in the following output:


### Data Table

"result_data.rds" file containing the following information for each cluster:

*	metadata distributions over time
*	tree statistics (e.g., Pybus's gamma, yule)
*	relevant values (e.g., R0, TMRCA, timespan)

This file can be opened with \url{need info from NANA analytics} for interactive viewing.

 
### Tree
"result_tree.tree" containing a nexus timed tree file with annotated information regarding cluster association. This file is meant to be viewed alongside the result_data.rds in \url{need info from NANA analytics} but can also be opened in various alternative applications including (but not limited to) BEAST, Figtree, and NextStrain. Taxa in the tree have been modofied to supply sampling dates from metadata file following "|". 


## References
Prosperi MCF (2011). A novel methodology for large-scale phylogeny partition. Nat Commun. 2: 321.

Paradis E, Claude J & Strimmer K (2004). APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289-290.

Revell, LJ (2012). phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol. 3:217-223.

Yu G, Smith DK, Zhu H, Guan Y, Lam TT (2016). GGTREE: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods Ecol. Evol. 8:28-36.

Pennell MW, Eastman JM, Slater GJ, Brown JW, Uyeda JC, Fitzjohn RG, Alfaro ME, Harmon LJ (2014). geiger v2.0: an expanded suite of methods for fitting macroevolutionary models to phylogenetic trees. Bioinformatics. 30:2216-2218.

Erik Volz (2020). treedater: Fast Molecular Clock Dating of Phylogenetic Trees with Rate Variation. R package version 0.5.0.  https://CRAN.R-project.org/package=treedater

