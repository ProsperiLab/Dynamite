# DYNAMITE (BETA)
DYNAMic Identification of Transmission Epicenters within phylogenies

## Description

DYNAMITE is a depth-first search algorithm that is not restricted to monophyletic clades, allowing for improved sensitivity to the detection of dynamic, minimal evolution chains (MECs) that comprise transmission patterns of interest for public health relevance. A more detailed description of the DYNAMITE algorithm and its applications can be found at https://www.biorxiv.org/content/10.1101/2021.01.21.427647v1. Since the release of this article, we have added a time component to the algorithm, allowing the user to restrain the inclusion of individuals in clusters based on a specified serial interval. This serial interval is defined as the period of time between a primary case-patient (infector) with symptom onset and a secondary case-patient, which has been observed to be 5-6 days for SARS-CoV-2 Rai et al., 2021. Direct transmission chains are intuitively represented by individuals separated by both minimal genetic evolution and time difference within this serial interval.

**This project is still in beta and as such may contain major bugs and we do not guarantee proper results/functionality**

## Requirements
R 3.6.0 or higher

**Libraries**

* optparse
* remotes
* phytools
* data.tree
* dplyr
* tidytree
* lubridate
* rlist
* familyR
* tidyverse
* ggtree
* parallel
* foreach
* geiger
* tibble
* treedater
* Rlsd2
* treeio

**File format requirements**

DYNAMITE takes a tree (scaled in substitutions/site) in newick or nexus format with support values at the nodes. Support values can be provided by any (single) method - i.e., combined result of 2 or more methods will not be accepted. These values can be scaled from 0-1 or 0-100.

Additionally, a metadata file with sequence IDs, sampling dates, and traits of interest is required. Column headers are needed and must contain the strings "ID" and "Date" (specific case format not required). Date information must be consistently either numeric or Date format. Remaining columns can consist of any number of traits of interest.

It is important to note that large trees (>30,000 sequences) cannot be scaled in time using the current implementation of treedater. While we have supplied a the least-squares dating function from Rlsd2 as an alternative, this function is also not entirely scaleable for large datasets. DYNAMITE has been successful with datasets of up to 10,000 sequences, though we are continuing to search for improved options for larger datasets.

## Execution 

```R
--help (-h) = provides helpful information for running DYNAMITE
--tree (-t) = path to treefile [default = .nwk file]
--metadata (-m) = path to metadata file [default= .csv file]
--cluster (-c) = choice of cluster algorithm from 'c' (Phylopart's cladewise) or 'b' (DYNAMITE's branchwise) [default= b]
--timetree (-q) = option (Y/N) for molecular clock calibration and time tree output/statistics [default=Y]
--seqLen (-s) = sequence length used in molecular clock calibration [default=30000]
--threshold (-l) = threshold for cluster determination, which can be numeric or "median" [default= 0.05]
--serial (-i) = serial interval (in days) for cluster filtering [default= 6]
--asr (-asr) = option (Y/N) of ancestral state reconstruction for each cluster [default= N]
```

**Important functions (with defaults):**

```R
treedater::dater(tree, dates, s=1000, ncpu=numCores, omega0=8E-04)

rlsd2::lsd2(inputTree=sub_tree,
                     inputDate=decimal_date(sts),
                     estimateRoot="as",
                     constraint=TRUE,
                     variance=1,
                     ZscoreOutlier = 3,
                     outFile = "lsd2_results",
                     seqLen = seqLen,
                     nullblen=-1

branchwise()

phylopart()

ancStateRecon(ace(model="ER", type="discrete"))
                 
```

## Output

DYNAMITE will result in the following output:


### Data Table

"trait_distributions_<original tree name>.csv" file containing the following information for each cluster:

*	metadata distributions over time
*	tree statistics (Infection rate [Oster], phylogenetic diversity [PD])
*	temporal information (TMRCA, timespan) if timetree option chosen

This file can be opened with \url{need info from NANA analytics} for interactive viewing.

 
### Optimal branch length
"branch_length_limit.txt" file containing the optimal threshold (% branch length distribution) and corresponding branch length

### Trees
"dynamite_subtree_<original tree name>.tree" containing a nexus tree file scaled in substitutions/site with annotated information regarding cluster association. 

"dynamite_timetree_<original tree name>.tree" containing a nexus tree file scaled in substitutions/site with annotated information regarding cluster association if the timetree option is chosen. 


These tree files can also be viewed alongside the output data table in \url{need info from NANA analytics} but can also be opened in various alternative applications including (but not limited to) BEAST, Figtree, and NextStrain. Taxa in the tree have been modofied to supply sampling dates from metadata file following "|". 

### Tree statistics summary table
A summary file detailing tree statistics such as the estimated infection rate, described by Oster et al. (2018) and overall phylogenetic diversity (sum of branch lengths), described by Faith (1992). These statistics are calculated for each cluster as well as for the background population.

### Run time
For a single tree with with 362 tips (see /simulations/ folder), the DYNAMITE branchwise algorithm with no serial time restraints completed in 5.571 seconds in R v4.0 on a MacBook Air (2 GHz Intel Core i7, 8GB).

### Visualization
An R script found in COVEME_2021/DEMO referred to as dynaviz.R can also be used to visualize the distributions (both interactively using the shiny R package and in PDF form), as well the timed tree annotated according to cluster and metadata information in the form of a heatmap.


## References
Rai B, Shukla A, Kant Dwivedi L (2021). Estimates of serial interval for COVID-19: A systematic review andmeta-analysis. Clin Epidemiol Glob Health. 9: 157-61.

Prosperi MCF (2011). A novel methodology for large-scale phylogeny partition. Nat Commun. 2: 321.

Paradis E, Claude J & Strimmer K (2004). APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289-290.

Revell, LJ (2012). phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol. 3:217-223.

Yu G, Smith DK, Zhu H, Guan Y, Lam TT (2016). GGTREE: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods Ecol. Evol. 8:28-36.

Pennell MW, Eastman JM, Slater GJ, Brown JW, Uyeda JC, Fitzjohn RG, Alfaro ME, Harmon LJ (2014). geiger v2.0: an expanded suite of methods for fitting macroevolutionary models to phylogenetic trees. Bioinformatics. 30:2216-2218.

Erik Volz (2020). treedater: Fast Molecular Clock Dating of Phylogenetic Trees with Rate Variation. R package version 0.5.0.  https://CRAN.R-project.org/package=treedater

