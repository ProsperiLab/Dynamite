# HIVdynamite (BETA)
dynamic identification of transmission clusters in phylogenies

<b>This project is still in beta and as such may contain major bugs and we do not guarantee proper results/functionality</b>

## Requirements
R 3.3 or higher

<b>Libraries</b>

* phytools
* picante
* ggtree
* colorspace

## Execution 

```R
source("HIVdynamite.R")
processTree(tree_file.nwk)
```

Important functions (with defaults):

```R
processTree(input_tree, slice_count=10, bootstrap=0.70, min_leaves=15, perc_dist=0.05)
plotClusters(percentiles, node_tree_file="./treeFigures/node_tree.nwk",
             processed_tree_table="./treeTables/processedTree.csv",
             width = 16, height = 16, units = "in", dpi = 600)
```

## Output

For example, a 10 slice tree will result in:

### treeTables
* <b>HIVdynamite.csv</b>                 Final table with cluster paths and repeats identified

### treeFigures
* <b>losalamos_small.nwk.LTT.pdf</b>: Lineage-through-time plot with extant and extinct lineages -- conducts γ-test of Pybus & Harvey (2000).
* <b>node_tree.nwk</b>: Tree with node identifiers
* <b>node_tree.pdf</b>: Plot of tree with node identifiers
* <b>tree.nwk</b>: Original tree with polytomies resolved to dichotomies
* <b>tree.pdf</b>: Plot of dichotomous tree
* <b>cluster_tree.pdf</b>: Plot of tree with colored clusters

### treeSlices (intermediate files)
* phylopart0.1.csv                Results from PhyloPart at 10th percentile
* phylopart0.2.csv                20th percentile
* phylopart0.3.csv                ...
* phylopart0.4.csv                ...
* phylopart0.5.csv                ...
* phylopart0.6.csv                ...
* phylopart0.7.csv                ...
* phylopart0.8.csv                ...
* phylopart0.9.csv                ...
* treeSlice0.1.nwk                Tree sliced at 10th percentile
* treeSlice0.1.pdf                Plot of tree sliced at 10th percentile
* treeSlice0.2.nwk                20th percentile
* treeSlice0.2.pdf                ...
* treeSlice0.3.nwk                ...
* treeSlice0.3.pdf                ...
* treeSlice0.4.nwk                ...
* treeSlice0.4.pdf                ...
* treeSlice0.5.nwk                ...
* treeSlice0.5.pdf                ...
* treeSlice0.6.nwk                ...
* treeSlice0.6.pdf                ...
* treeSlice0.7.nwk                ...
* treeSlice0.7.pdf                ...
* treeSlice0.8.nwk                ...
* treeSlice0.8.pdf                ...
* treeSlice0.9.nwk                ...
* treeSlice0.9.pdf                ...
* renamed_leaves0.1.csv           Old_leaf_name to New_leaf_name table for 10th percentile
* renamed_leaves0.2.csv           ...
* renamed_leaves0.3.csv           ...
* renamed_leaves0.4.csv           ...
* renamed_leaves0.5.csv           ...
* renamed_leaves0.6.csv           ...
* renamed_leaves0.7.csv           ...
* renamed_leaves0.8.csv           ...
* renamed_leaves0.9.csv           ...
* renamed_leaves1.csv             ...

## Notes
### Node enumeration
By convention, the tips of the tree are numbered 1 through n for n tips; and the nodes are numbered n + 1 through n + m for m nodes. m = n - 1 for a fully bifurcating tree. So, if there are 6 leaves, then the root number is 7. From the root (7), nodes are enumerated in a depth-first manner (This follows the convention of the R package ape numbering scheme).

![Alt text](./Images/numbering_example.png?raw=true "Numbered Tree")

This example image may be generated using:

```R
library(ape)

tree <- read.tree(text = "(((A,B),(C,D)),(E,F));")
plot(tree, edge.width = 2, label.offset = 0.1)
nodelabels()
tiplabels()
```

