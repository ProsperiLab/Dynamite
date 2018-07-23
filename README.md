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
processTree("losalamos_small.nwk")
```

Important functions (with defaults):

```R
processTree(input_tree, slice_count=10, bootstrap=0.70, min_leaves=15, perc_dist=0.05)
plotClusters(percentiles, node_tree_file="./treeFigures/node_tree.nwk",
             processed_tree_table="./treeTables/processedTree.csv",
             width = 16, height = 16, units = "in", dpi = 600)
```

## Output

For example, a 10 slice tree will result in 3 directories containing:

### treeTables
* <b>HIVdynamite.csv</b>: Primary output consisting of a final table with cluster paths and repeats identified

 | TreeLevel | MedianOfDistances | ClusterName | SequencesPerCluster | Bootstrap | Path | PybusGamma | CutHeight | MRCANode | MRCAHeight | SpeciationRate | PhyloDiversity | RepeatCluster | Leaves |
 | --------- | ----------------- | ----------- | ------------------- | --------- | ---- | ---------- | --------- | -------- | ---------- | -------------- | -------------- | ------------- | ------ |
 | 6 | 0.00349160000000001 | c6.1 | 5 | 0.9999 | c6.1 | 4.06535504585726 | 0.0770163 | 232 | 0.0752705 | 93.4310414402262 | 0.0810695 | - | 234 | 235 | 238 | 244 | 251 |
 | 7 | 0.0152437 | c7.1 | 2 | 0.99 | c7.1 | -3.50603898980592 | 0.08985235 | 133 | 0.0822305 | 57.8392248849795 | 0.0974742 | - | 84@57_0$0.00762185-133 | 122@95_0$0.00762185-133 |
 | 7 | 0.0222637 | c7.2 | 8 | 0.998 | c7.2 | -3.50603898980592 | 0.08985235 | 149 | 0.0781005 | 57.8392248849795 | 0.15892715 | - | 158 | 4@4_2$0.00702185-162 | 15@12_2$0.00264185-163 | 150 | 160 | 13@11_2$0.00943185-161 | 8@8_2$0.00028185-159 | 107@80_0$0.01647185-217 |
 | 7 | 0.00118370000000001 | c7.3 | 2 | 1 | c7.3 | -3.50603898980592 | 0.08985235 | 166 | 0.0892605 | 57.8392248849795 | 0.0904442 | - | 167 | 200 |
 | 7 | 0.0217837 | c7.4 | 2 | 0.988 | c7.4 | -3.50603898980592 | 0.08985235 | 209 | 0.0789605 | 57.8392248849795 | 0.1007442 | - | 1@1_1$0.01089185-209 | 2@2_1$0.01089185-209 |
 | 7 | 0.0216437 | c7.5 | 2 | 0.946 | c7.5 | -3.50603898980592 | 0.08985235 | 211 | 0.0790305 | 57.8392248849795 | 0.1006742 | - | 79@52_0$0.01082185-211 | 121@94_0$0.01082185-211 |
 | 7 | 0.0196037 | c7.6 | 2 | 0.902 | c7.6 | -3.50603898980592 | 0.08985235 | 223 | 0.0800505 | 57.8392248849795 | 0.0996542 | - | 87@60_0$0.00980185-223 | 102@75_0$0.00980185-223 |
 | 7 | 0.0215837 | c7.7 | 3 | 0.905 | c7.7 | -3.50603898980592 | 0.08985235 | 230 | 0.0790605 | 57.8392248849795 | 0.11014605 | - | 68@41_0$0.01079185-230 | 81@54_0$0.00950185-231 | 101@74_0$0.00950185-231 |
 | 7 | 0.0274237 | c7.8 | 15 | 0.9999 | c6.1 | c7.8 | -3.50603898980592 | 0.08985235 | 232 | 0.0752705 | 57.8392248849795 | 0.20808825 | - | 92@65_0$0.00534185-238 | 243 | 119@92_0$0.01066185-235 | 89@62_0$0.01066185-235 | 19@15_3$0.00025185-246 | 67@40_0$0.01230185-234 | 120@93_0$0.01166185-251 | 112@85_0$0.00927185-244 | 124@97_0$0.00221185-248 | 83@56_0$0.01230185-234 | 249 | 247 | 240 | 64@37_0$0.01166185-251 | 242 |

### treeFigures
* <b>losalamos_small.nwk.LTT.pdf</b>: Lineage-through-time plot with extant and extinct lineages -- conducts γ-test of Pybus & Harvey (2000).
* <b>node_tree.nwk</b>: Tree with node identifiers
* <b>node_tree.pdf</b>: Plot of tree with node identifiers
* <b>tree.nwk</b>: Original tree with polytomies resolved to dichotomies
* <b>tree.pdf</b>: Plot of dichotomous tree
* <b>cluster_tree.pdf</b>: Plot of tree with colored clusters

<b>Figure cluster_tree.pdf</b>
![Alt text](./Images/cluster_tree.png?raw=true "Clustered Tree")

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

