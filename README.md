# DYNAMITE (BETA)
DYNAMic Identification of Transmission Epicenters within phylogenies

**This project is still in beta and as such may contain major bugs and we do not guarantee proper results/functionality**

## Requirements
R 3.6.0 or higher

**Libraries**

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
* rlsd2

**File format requirements**

DYNAMITE requires a tree (scaled in substitutions/site) in newick or nexus format with support values at the nodes. Support values can be provided by any (single) method - i.e., combined result of 2 or more methods will not be accepted. These values can be scaled from 0-1 or 0-100.

Additionally, a metadata file with sampling dates and traits of interest is required. Column headers are needed. The first column must correspond to the taxa names, the second column dates (must be consistently either numeric or Date format), and remainig columns any number of traits of interest. Specific header names are not required, as only the order is called.


## Execution 

```R
arg1 <- path_to_tree_file
arg2 <- path_to_metadata_file
arg3 <- sequence_length
system(paste("dynamite.R", arg1, arg2, arg3))
```

**Important functions (with defaults):**

```R
lsd2(estimateRoot="as", constraint=TRUE, variance=1, ZscoreOutlier = 3, seqLen = arg3, nullblen=-1)
  
phylopart(phylopart.threshold=0.10)

ancStateRecon(ace(model="ER", type="discrete"))

calculateNe(skygrowth.map(res = 10, tau0 = 0.1))

calculateRe(conf.level=0.95, s=100, psi=rnorm(s, 14, 5))

yule()

ltt()                 
```

## Output

DYNAMITE will result in the following output:


### Data Table
"result_data.rds" file containing metadata distributions over time, tree statistics (e.g., Pybus's gamma), and relevant values (e.g., R0 and TMRCA) for each identified cluster. This file can be opened with \url{need info from NANA analytics} for interactive viewing.

 
### Tree
"result_tree.tree" containing a nexus timed tree file with annotated information regarding cluster association. This file is meant to be viewed alongside the result_data.rds in \url{need info from NANA analytics} but can also be opened in various alternative applications including (but not limited to) BEAST, Figtree, and NextStrain. Taxa in the tree have been modofied to supply sampling dates from metadata file following "|". 


## References
Paradis E, Claude J & Strimmer K (2004). APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289-290.

Yu G, Smith D, Zhu H, Guan Y, Lam TT (2017). ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods Ecol. Evol. 8:28-36.

Revell, LJ (2012). phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol. 3:217-223.

Yu G, Smith DK, Zhu H, Guan Y, Lam TT (2016). GGTREE: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods Ecol. Evol. 8:28-36.

Pennell MW, Eastman JM, Slater GJ, Brown JW, Uyeda JC, Fitzjohn RG, Alfaro ME, Harmon LJ (2014). geiger v2.0: an expanded suite of methods for fitting macroevolutionary models to phylogenetic trees. Bioinformatics. 30:2216-2218.

Volz EM & Didelot X (2018). Modeling the growth and decline of pathogen effective population size provides insight into epidemic dynamics and drivers of antimicrobial resistance. Syst. Biol. 67:719-728.

To T-H, Jung M, Lycett S, Gascuel O (2016) Fast dating using least-squares criteria and algorithms. Syst. Biol. 65:82-97.

Pybus OG, Harvey PH (2000). Testing macro-evolutionary models using incomplete molecular phylogenies. Proc. Roy. Soc. Lond. B: Biol. Sci. 267:2267–2272.

Holmes EC, Nee S, Rambaut A, Garnett GP, Harvey PH (1995). Revealing the history of infectious disease epidemics through phylogenetic trees. Philos. Trans. Roy. Soc. B: Biol. Sci. 349:33.


