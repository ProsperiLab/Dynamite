VERSION 2.1 (2016)

Contact: Dr. Iuri Fanti <iuri.fanti@gmail.com>

Usage: double click on the "PhyloPart_v2.0.jar" file
or (for MS Windows users) double click on the "phyloPart_v2.0_GUI.bat" file
or type "java [-Xmx{1,2,3...}G] -jar PhyloPart_v2.0.jar {GUI | treeFilePath {percentileThreshold | -DIST | -LEVELS} -oOutputFile}" from a terminal/command line prompt

Examples of input trees (standard NEWICK format):
	a.	"((pippo_blabla:0.09,pluto_blabla:0.007)0.89:0.076,topolino_blabla:0.81);"
	b.	"((leaf1@01_:0.09,leaf2@02_:0.007)0.89:0.076,leaf3@01_:0.81);"

	Rules:
	- The tree is surrounded by round brackets and ends with a semicolon.
	- At the right of the leaves there is a colon and the branch_length.
	- At the right of the intermediate nodes there is NO COLON, the bootstrap_value, a colon and the branch_length.
	- Intermediate nodes are surrounded by round brackets.
	- All nodes are separated by comma.
	- Each leaf name should be formatted in this way: "idx_anythingelse" (example a), if you have multiple observations per patient, you can set up the identifier as "idx@patientid_anythingelse" (example b).


The program has THREE working modalities (both from the command line and GUI interface):


1. CALCULATE A PARTITION OF A PHYLOGENETIC TREE WITH A GIVEN PERCENTILE DISTANCE THRESHOLD (5th perc. for instance)
	"java -jar phyloPart.jar yourTree.nwk 0.05 -oOUT.txt"
	Note: for large trees, when using the command-line, the distance calculation is limited to a maximum of 1000000 pairwise calculations (use always at least -Xmx2G for trees with >10000 leaves, better -Xmx4G). The user can augment this limit using the GUI interface (and presumably setting up at least -Xmx8G for trees with >10000 leaves).
	1.1 CALCULATE A SET OF PARTITIONS OF A PHYLOGENETIC TREE WITH A GIVEN PERCENTILE DISTANCE THRESHOLD INTERVAL
		Same as above, but from the GUI a user can run the program using a set of thresholds (specifying min, max and steps)

2. CALCULATE PAIRWISE PATRISTIC DISTANCES BETWEEN LEAVES OF THE TREE
	"java -jar phyloPart.jar yourTree.nwk -DIST -oOUT.txt"

3. CALCULATE TREE TOPOLOGY STATISTICS
	"java -jar phyloPart.jar yourTree.nwk -LEVELS -oOUT.txt"


All output files are plain texts in .csv format.


-----------
WALKTHROUGH 
-----------

The program is easy to run using the graphical interface. Just double click on the "phyloPart_GUI.bat" file, otherwise open a terminal and run the phyloPart.jar with the GUI option (see syntax above).
After the program starts, you can select a rooted phylogenetic tree file (it must be in newick format, with bootstrap/posterior probability/reliability values). 
If you want to perform clustering, select the Depth-first option and choose a percentile threshold value between 0 and 1 (this is related to the tree patrisitc distance and can make clusters bigger or smaller, I suggest you to try 0.05 first, which corresponds to the 5th percentile). 
If your tree is big (>1000 or >10000 leaves), you might decide to set up a sample distance limit, say 50000/100000, or increase the RAM usage (modifying the -Xmx option in the .bat file, which now is set to 2Gb by default, for instance using -Xmx4G or -Xmx8G). You need to specify an output file (write any name, like result.txt). Then click "run".

The output file is written in csv format with the following columns:

clustername;bootstrap; leafname;branchPath; medianOfDistances; sequencesperCluster

clustername = name of the cluster (can be a number, but the number has no meaning, sequences in the same cluster have the same number. The special number "0" means that the sequence was not placed in any cluster)
bootstrap = the bootstrap value (or any reliability/posterior probability value) of the cluster, it will be always >90%
leafname = name of the sequence
branchPath = root-to-tip distance of the sequence
medianOfDistances = median of patristic distances in the cluster
sequencesperCluster = number of sequences per cluster

In the main program window you will see that many lines are written. If you scroll down, you will find a special keyword "threshold:" (repeated many times) that tells you which is the corresponding nucleotide substitution per site value of the threshold used in the input (so if you change the percentile, this value will change).
