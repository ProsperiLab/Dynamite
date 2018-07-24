# HIVdynamite (BETA)
dynamic identification of transmission clusters in phylogenies

**This project is still in beta and as such may contain major bugs and we do not guarantee proper results/functionality**

## Requirements
R 3.3 or higher

**Libraries**

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
* **HIVdynamite.csv**: Primary output consisting of a final table with cluster paths and repeats identified

This table (depicted below) contains 14 columns:

- **TreeLevel** - the percentile cut height at which the current cluster occurs
- **MedianOfDistance** - median of patristic distances in the cluster
- **ClusterName** - the name of the current cluster
- **SequencesPerCluster** - number of leaves in the current cluster
- **Bootstrap** - the bootstrap score of the current cluster
- **Path** - the clusters at lower cuts (closer to root) containing the current cluster
- **PybusGamma** - result of the γ-test of Pybus & Harvey (2000).
- **CutHeight** - the height of the cut along the full-length tree
- **MRCANode** - the node representing the most recent common ancestor (MRCA) of the current cluster
- **MRCAHeight** - the height from the MRCANode
- **SpeciationRate** - measure of speciation rate generated using **yule** from the ape package
- **PhyloDiversity** - measure of phylogenetic diversity  generated using **pd** from the picante package
- **RepeatCluster** - indicates identical clusters occurring at lower cuts
- **Leaves** - names of the leaves within the current cluster

 | TreeLevel | MedianOfDistances | ClusterName | SequencesPerCluster | Bootstrap | Path | PybusGamma | CutHeight | MRCANode | MRCAHeight | SpeciationRate | PhyloDiversity | RepeatCluster | Leaves |
 | --------- | ----------------- | ----------- | ------------------- | --------- | ---- | ---------- | --------- | -------- | ---------- | -------------- | -------------- | ------------- | ------ |
 | 6 | 0.00349160000000001 | c6.1 | 5 | 0.9999 | c6.1 | 4.06535504585726 | 0.0770163 | 232 | 0.0752705 | 93.4310414402262 | 0.0810695 | - | 235;244;238;234;251 |
 | 7 | 0.0152437 | c7.1 | 2 | 0.99 | c7.1 | -3.50603898980592 | 0.08985235 | 133 | 0.0822305 | 57.8392248849795 | 0.0974742 | - | 84@57_0$0.00762185-133;122@95_0$0.00762185-133 |
 | 7 | 0.0222637 | c7.2 | 8 | 0.998 | c7.2 | -3.50603898980592 | 0.08985235 | 149 | 0.0781005 | 57.8392248849795 | 0.15892715 | - | 13@11_2$0.00943185-161;150;107@80_0$0.01647185-217;160;4@4_2$0.00702185-162;158;8@8_2$0.00028185-159;15@12_2$0.00264185-163 |
 | 7 | 0.00118370000000001 | c7.3 | 2 | 1 | c7.3 | -3.50603898980592 | 0.08985235 | 166 | 0.0892605 | 57.8392248849795 | 0.0904442 | - | 167;200 |
 | 7 | 0.0217837 | c7.4 | 2 | 0.988 | c7.4 | -3.50603898980592 | 0.08985235 | 209 | 0.0789605 | 57.8392248849795 | 0.1007442 | - | 1@1_1$0.01089185-209;2@2_1$0.01089185-209 |
 | 7 | 0.0216437 | c7.5 | 2 | 0.946 | c7.5 | -3.50603898980592 | 0.08985235 | 211 | 0.0790305 | 57.8392248849795 | 0.1006742 | - | 79@52_0$0.01082185-211;121@94_0$0.01082185-211 |
 | 7 | 0.0196037 | c7.6 | 2 | 0.902 | c7.6 | -3.50603898980592 | 0.08985235 | 223 | 0.0800505 | 57.8392248849795 | 0.0996542 | - | 87@60_0$0.00980185-223;102@75_0$0.00980185-223 |
 | 7 | 0.0215837 | c7.7 | 3 | 0.905 | c7.7 | -3.50603898980592 | 0.08985235 | 230 | 0.0790605 | 57.8392248849795 | 0.11014605 | - | 68@41_0$0.01079185-230;81@54_0$0.00950185-231;101@74_0$0.00950185-231 |
 | 7 | 0.0274237 | c7.8 | 15 | 0.9999 | c6.1;c7.8 | -3.50603898980592 | 0.08985235 | 232 | 0.0752705 | 57.8392248849795 | 0.20808825 | - | 240;120@93_0$0.01166185-251;247;64@37_0$0.01166185-251;119@92_0$0.01066185-235;249;92@65_0$0.00534185-238;242;112@85_0$0.00927185-244;124@97_0$0.00221185-248;243;67@40_0$0.01230185-234;83@56_0$0.01230185-234;89@62_0$0.01066185-235;19@15_3$0.00025185-246 |
 | 8 | 0.0409158 | c8.1 | 2 | 0.99 | c7.1;c8.1 | -3.29374770216132 | 0.1026884 | 133 | 0.0822305 | 45.9755454590621 | 0.1231463 | - | 84@57_0$0.0204579-133;122@95_0$0.0204579-133 |
 | 8 | 0.0371679 | c8.2 | 13 | 0.998 | c7.2;c8.2 | -3.29374770216132 | 0.1026884 | 149 | 0.0781005 | 45.9755454590621 | 0.2222842 | - | 3@3_2;97@70_0;98@71_0;7@7_2;5@5_2;6@6_2;8@8_2$0.0131179-159;4@4_2;152;12@11_2;15@12_2$0.0154779-163;10@10_2;13@11_2 |
 | 8 | 0.0258379 | c8.3 | 18 | 0.912 | c7.3;c8.3 | -3.29374770216132 | 0.1026884 | 165 | 0.0724805 | 45.9755454590621 | 0.1706211 | - | 69@42_0;55@34_5$0.0005079-202;35@24_5$0.0003779-203;50@33_5$0.0002679-195;60@35_5;40@25_5$0.0005679-187;33@24_5$0.0096779-200;169;28@20_5$0.0013179-197;192;62@35_5$0.0002679-195;31@22_5$0.0014779-199;191;48@32_5$0.0013179-197;41@26_5$0.0049879-201;34@24_5$0.0003779-203;43@28_5;45@30_5 |
 | 8 | 0.032 | c8.4 | 2 | 0.988 | c7.4;c8.4 | -3.29374770216132 | 0.1026884 | 209 | 0.0789605 | 45.9755454590621 | 0.1109605 | - | 1@1_1;2@2_1 |
 | 8 | 0.0368179 | c8.5 | 2 | 0.946 | c7.5;c8.5 | -3.29374770216132 | 0.1026884 | 211 | 0.0790305 | 45.9755454590621 | 0.1158484 | - | 79@52_0$0.0236579-211;121@94_0 |
 | 8 | 0.0179079 | c8.6 | 2 | 1 | c8.6 | -3.29374770216132 | 0.1026884 | 224 | 0.0923905 | 45.9755454590621 | 0.1102984 | - | 93@66_0$0.0102979-224;94@67_0 |
 | 8 | 0.0363558 | c8.7 | 4 | 0.998 | c6.1;c7.8;c8.7 | -3.29374770216132 | 0.1026884 | 238 | 0.0845105 | 45.9755454590621 | 0.1523221 | - | 242;240;243;92@65_0$0.0181779-238 |
 | 8 | 0.0348579 | c8.8 | 7 | 0.966 | c6.1;c7.8;c8.8 | -3.29374770216132 | 0.1026884 | 244 | 0.0805805 | 45.9755454590621 | 0.1654721 | - | 19@15_3;250;112@85_0$0.0221079-244;113@86_0;20@15_3;124@97_0$0.0150479-248;18@14_3$0.0027079-249 |
 | 9 | 0.01584395 | c9.1 | 3 | 1 | c7.2;c8.2;c9.1 | -3.29374770216132 | 0.11552445 | 152 | 0.1067005 | 47.7620532715653 | 0.1282884 | - | 16@13_2;154;14@11_2$0.00574395-153 |
 | 9 | 0.00328000000000001 | c9.2 | 2 | 0.976 | c7.2;c8.2;c9.2 | -3.29374770216132 | 0.11552445 | 158 | 0.0912005 | 47.7620532715653 | 0.0944805 | - | 3@3_2;10@10_2 |
 | 9 | 0 | c9.3 | 2 | 0.9999 | c7.2;c8.2;c9.3 | -3.29374770216132 | 0.11552445 | 160 | 0.1011005 | 47.7620532715653 | 0.1011005 | - | 6@6_2;97@70_0 |
 | 9 | 0 | c9.4 | 2 | 0.9999 | c7.2;c8.2;c9.4 | -3.29374770216132 | 0.11552445 | 164 | 0.0968405 | 47.7620532715653 | 0.0968405 | - | 7@7_2;98@71_0 |
 | 9 | 0.01577395 | c9.5 | 34 | 1 | c7.3;c8.3;c9.5 | -3.29374770216132 | 0.11552445 | 167 | 0.0989305 | 47.7620532715653 | 0.25527395 | - | 59@35_5;52@33_5$0.00288395-185;29@21_5$0.00556395-179;39@25_5;40@25_5;58@35_5$0.00447395-178;30@21_5;32@23_5$0.00556395-179;31@22_5;56@35_5$0.00702395-180;36@25_5$0.00702395-180;37@25_5$0.00288395-185;54@33_5;28@20_5;53@33_5;57@35_5$0.00388395-175;24@19_5$0.00538395-184;45@30_5;61@35_5;62@35_5;26@19_5;25@19_5;49@33_5$0.00447395-178;51@33_5;47@31_5;42@27_5;27@19_5;43@28_5;48@32_5;38@25_5;46@31_5$0.00678395-183;44@29_5;50@33_5;60@35_5 |
 | 9 | 0.01332 | c9.6 | 4 | 0.977 | c7.3;c8.3;c9.6 | -3.29374770216132 | 0.11552445 | 201 | 0.0977005 | 47.7620532715653 | 0.1189505 | - | 41@26_5;55@34_5;35@24_5;34@24_5 |
 | 9 | 0 | c9.7 | 2 | 0.9999 | c6.1;c7.8;c8.7;c9.7 | -3.29374770216132 | 0.11552445 | 242 | 0.1063205 | 47.7620532715653 | 0.1063205 | - | 22@17_4;95@68_0 |
 | 9 | 0 | c9.8 | 2 | 0.9999 | c6.1;c7.8;c8.8;c9.8 | -3.29374770216132 | 0.11552445 | 247 | 0.1011705 | 47.7620532715653 | 0.1011705 | - | 20@15_3;113@86_0 |
 | 9 | 0 | c9.9 | 2 | 0.9999 | c6.1;c7.8;c8.8;c9.9 | -3.29374770216132 | 0.11552445 | 250 | 0.1064205 | 47.7620532715653 | 0.1064205 | - | 17@14_3;100@73_0 |
 | 10 | 0.01612 | c10.1 | 5 | 1 | c7.2;c8.2;c9.1;c10.1 | -3.29374770216132 | 0.1283605 | 152 | 0.1067005 | 47.9902161881836 | 0.1335105 | - | 9@9_2;11@11_2;99@72_0;14@11_2;16@13_2 |
 | 10 | 0.00328000000000001 | c10.2 | 2 | 0.976 | c7.2;c8.2;c9.2;c10.2 | -3.29374770216132 | 0.1283605 | 158 | 0.0912005 | 47.9902161881836 | 0.0944805 | c9.2 | 3@3_2;10@10_2 |
 | 10 | 0 | c10.3 | 2 | 0.9999 | c7.2;c8.2;c9.3;c10.3 | -3.29374770216132 | 0.1283605 | 160 | 0.1011005 | 47.9902161881836 | 0.1011005 | c9.3 | 6@6_2;97@70_0 |
 | 10 | 0 | c10.4 | 2 | 0.9999 | c7.2;c8.2;c9.4;c10.4 | -3.29374770216132 | 0.1283605 | 164 | 0.0968405 | 47.9902161881836 | 0.0968405 | c9.4 | 7@7_2;98@71_0 |
 | 10 | 0.01703 | c10.5 | 34 | 1 | c7.3;c8.3;c9.5;c10.5 | -3.29374770216132 | 0.1283605 | 167 | 0.0989305 | 47.9902161881836 | 0.2746605 | - | 29@21_5;37@25_5;51@33_5;39@25_5;38@25_5;58@35_5;48@32_5;43@28_5;61@35_5;45@30_5;59@35_5;50@33_5;31@22_5;26@19_5;53@33_5;47@31_5;62@35_5;57@35_5;27@19_5;54@33_5;36@25_5;24@19_5;40@25_5;42@27_5;28@20_5;30@21_5;52@33_5;60@35_5;32@23_5;49@33_5;44@29_5;46@31_5;25@19_5;56@35_5 |
 | 10 | 0.01332 | c10.6 | 4 | 0.977 | c7.3;c8.3;c9.6;c10.6 | -3.29374770216132 | 0.1283605 | 201 | 0.0977005 | 47.9902161881836 | 0.1189505 | c9.6 | 41@26_5;55@34_5;35@24_5;34@24_5 |
 | 10 | 0.01882 | c10.7 | 2 | 1 | c8.6;c10.7 | -3.29374770216132 | 0.1283605 | 224 | 0.0923905 | 47.9902161881836 | 0.1112105 | - | 93@66_0;94@67_0 |
 | 10 | 0 | c10.8 | 2 | 0.9999 | c6.1;c7.8;c8.7;c10.8 | -3.29374770216132 | 0.1283605 | 240 | 0.1208705 | 47.9902161881836 | 0.1208705 | - | 21@16_4;96@69_0 |
 | 10 | 0 | c10.9 | 2 | 0.9999 | c6.1;c7.8;c8.7;c9.7;c10.9 | -3.29374770216132 | 0.1283605 | 242 | 0.1063205 | 47.9902161881836 | 0.1063205 | c9.7 | 22@17_4;95@68_0 |
 | 10 | 0 | c10.10 | 2 | 0.9999 | c6.1;c7.8;c8.7;c10.10 | -3.29374770216132 | 0.1283605 | 243 | 0.1194505 | 47.9902161881836 | 0.1194505 | - | 23@18_4;105@78_0 |
 | 10 | 0.01782 | c10.11 | 3 | 0.98 | c6.1;c7.8;c8.8;c9.8;c10.11 | -3.29374770216132 | 0.1283605 | 246 | 0.0896005 | 47.9902161881836 | 0.1074205 | - | 19@15_3;20@15_3;113@86_0 |
 | 10 | 0.0176 | c10.12 | 3 | 0.999 | c6.1;c7.8;c8.8;c9.9;c10.12 | -3.29374770216132 | 0.1283605 | 249 | 0.0999805 | 47.9902161881836 | 0.1175805 | - | 17@14_3;100@73_0;18@14_3 |

### treeFigures
* **losalamos_small.nwk.LTT.pdf**: Lineage-through-time plot with extant and extinct lineages -- conducts γ-test of Pybus & Harvey (2000).
* **node_tree.nwk**: Tree with node identifiers
* **node_tree.pdf**: Plot of tree with node identifiers
* **tree.nwk**: Original tree with polytomies resolved to dichotomies
* **tree.pdf**: Plot of dichotomous tree
* **cluster_tree.pdf**: Plot of tree with colored clusters

**Figure cluster_tree.pdf**. Grey vertical lines represent the percentile cut heights along the tree. Clusters sharing a path are represented by different shades of the same color. Clade name represents the cluster path for the corresponding cluster. Repeated clusters are not represented. Coloring begins at the MRCA of the cluster not at the slice height.
![Alt text](./Images/cluster_tree.png?raw=true "Clustered Tree")

### treeSlices (intermediate files)
* phylopart0.1.csv: Results from PhyloPart at 10th percentile
* phylopart0.2.csv: 20th percentile
* phylopart0.3.csv: ...
* phylopart0.4.csv: ...
* phylopart0.5.csv: ...
* phylopart0.6.csv: ...
* phylopart0.7.csv: ...
* phylopart0.8.csv: ...
* phylopart0.9.csv: ...
* treeSlice0.1.nwk: Tree sliced at 10th percentile
* treeSlice0.1.pdf: Plot of tree sliced at 10th percentile
* treeSlice0.2.nwk: 20th percentile
* treeSlice0.2.pdf: ...
* treeSlice0.3.nwk: ...
* treeSlice0.3.pdf: ...
* treeSlice0.4.nwk: ...
* treeSlice0.4.pdf: ...
* treeSlice0.5.nwk: ...
* treeSlice0.5.pdf: ...
* treeSlice0.6.nwk: ...
* treeSlice0.6.pdf: ...
* treeSlice0.7.nwk: ...
* treeSlice0.7.pdf: ...
* treeSlice0.8.nwk: ...
* treeSlice0.8.pdf: ...
* treeSlice0.9.nwk: ...
* treeSlice0.9.pdf: ...
* renamed_leaves0.1.csv: Old_leaf_name to New_leaf_name table for 10th percentile
* renamed_leaves0.2.csv: ...
* renamed_leaves0.3.csv: ...
* renamed_leaves0.4.csv: ...
* renamed_leaves0.5.csv: ...
* renamed_leaves0.6.csv: ...
* renamed_leaves0.7.csv: ...
* renamed_leaves0.8.csv: ...
* renamed_leaves0.9.csv: ...
* renamed_leaves1.csv: ...

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

### Leaf names
Since clusters are generated from slices taken across the original tree, leaf names must be edited accordingly. For example, looking at some leaves from cluster c7.2:

**13@11_2$0.00943185-161;150;107@80_0$0.01647185-217;160;**

Leaf names may be changed to the corresponding node (as seen by leaf 160) or if the cut is along the branch leading to the leaf, the name becomes *LeafName$HeightFromNode-MRCANode* to indicate that the leaf in that cluster represents an ancestral state.

### Java error
Phylogenetic tree formatting may cause issues when running HIVdynamite. An example error:
```Java
Exception in thread "main" java.lang.Exception: Error while parsing tree: parsed tree is different form the one in the file.
```

One possible solution involves checking the parentheses around the tree. 

**Fails**: (122:0.2245,148:0.2245):1.000;
**Works**: ((122:0.2245,148:0.2245):1.000);

Otherwise, be sure the tree format matches the description in PhyloPart_README.txt

## References
Paradis E, Claude J & Strimmer K (2004). APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289-290.

Ihaka R, Murrell P, Hornik K, Fisher JC, Zeileis A (2016). colorspace: Color Space Manipulation. R package version 1.3-2. URL https://CRAN.R-project.org/package=colorspace

Yu G, Smith D, Zhu H, Guan Y, Lam TT (2017). ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods in Ecology and Evolution, 8, 28-36.

Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research, 43(7), e47. 

Revell, LJ (2012). phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol. 3 217-223.

Kembel SW, Cowan PD, Helmus MR, Cornwell WK, Morlon H, Ackerly DD, Blomberg SP, and Webb CO (2010). Picante: R tools for integrating phylogenies and ecology. Bioinformatics 26:1463-1464. 

