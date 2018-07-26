#!/usr/bin/env Rscript

tryInstall <- function(missing_lib){

    write(paste(missing_lib, "not installed"), stderr())
    write("Installing required packages", stderr())
    install = ""
    while ((install != "y") || (install != "n")){
        install <- readline(prompt=paste("Install ", missing_lib, "? [y/n]\n", sep=""))
        if (install == "y"){
            install.packages(missing_lib)
        } else if (install == "n"){
            quit(status=0)
        }# end if-else
    }# end while
}# end tryInstall

tryCatch({
            write("Loading required package: phytools",stderr())
            library(phytools)
            write("Loading required package: picante", stderr())
            library(picante)
            write("Loading required package: ggtree", stderr())
            library(ggtree)
            write("Loading required package: colorspace", stderr())
            library(colorspace)
         }, error = function(err){
            if (grepl("picante", err)){
                tryInstall("picante")
            }
            if (grepl("phytools", err)){
                tryInstall("phytools")
            }
            if (grepl("ggtree", err)){
                tryInstall("ggtree")
            }
            if (grepl("colorspace", err)){
                tryInstall("colorspace")
            }
            quit(status=0)
         })# end tryCatch

# Function to read Nexus files or Newick files
checkAndConvertTree <- function(filename, percentiles){
    # Check tree format -- either Newick or Nexus -- using suffices
    if (endsWith(filename, "nwk") || endsWith(filename,"newick")){
        write("Newick file detected", stderr())
        return(filename)
    } else if (endsWith(filename, "nex") || endsWith(filename,"nexus") || endsWith(filename, "nxs")){
        write("Nexus file detected",stderr())
        # If file is Nexus format, the Java program can't read it.
        # So, we use the library 'ape' to read Nexus and write Newick
        nexus_tree <- read.nexus(filename)
        nexus_name <- paste(removeExt(c(filename)),"converted.nwk",sep=".")
        write.tree(nexus_tree, file=nexus_name)
        return(nexus_name)
    } else {
        # Neither Newick nor Nexus identified -- quit
        write("Cannot identify tree format. Use nexus (nex,nxs,nexus) or newick (nwk,newick) formats.",stderr())
        quit(status=1)
    }# end if-else
}# end check_and_convert_tree

# Remove extensions -- originally from limma by Gordon Smyth
removeExt <- function (x, sep = ".") 
{
    x <- as.character(x)
    n <- length(x)
    if (length(grep(sep, x)) < n) 
        return(x)
    sep <- protectMetachar(sep)
    RegExpr <- paste0("(.*)", sep, "(.*)$")
    ext <- sub(RegExpr, "\\2", x)
    if (all(ext[1] == ext)) 
        return(sub(RegExpr, "\\1", x))
    else return(x)
}# end removeExt

# Function necessary for removing extensions -- originally from limma by Gordon Smyth
protectMetachar <- function(x) 
{
    x <- gsub("\\.", "\\\\.", x)
    x <- gsub("\\|", "\\\\|", x)
    x <- gsub("\\(", "\\\\(", x)
    x <- gsub("\\)", "\\\\)", x)
    x <- gsub("\\[", "\\\\[", x)
    x <- gsub("\\{", "\\\\{", x)
    x <- gsub("\\^", "\\\\^", x)
    x <- gsub("\\$", "\\\\$", x)
    x <- gsub("\\*", "\\\\*", x)
    x <- gsub("\\+", "\\\\+", x)
    x <- gsub("\\?", "\\\\?", x)
    x
}# end protectMetachar

# Function to plot trees
treePlot <- function(tree, percentiles){
    # Generate PDF image
    write.tree(tree, "./treeFigures/tree.nwk", tree.names=TRUE)
    gg <- ggtree(tree, ladderize=T)
    # Add percentile slices across the plot
    for (i in 1:length(percentiles)){
        gg <- gg + geom_vline(xintercept = max(gg$data$x) * percentiles[i], color = "grey", size = 0.3)
    }
    tree_height <- max(nodeHeights(tree))
    gg <- gg + geom_nodelab(cex=2)
    gg <- gg + geom_tiplab(cex=2) + xlim(0, max(gg$data$x) + 0.15*max(gg$data$x))
    ggsave(gg, file="./treeFigures/tree.pdf")
    
    node_tree <- treeSlice(tree, tree_height+1, orientation="rootwards")
    write.tree(node_tree, "./treeFigures/node_tree.nwk", tree.names=TRUE)
    # Plot node_tree with slices
    gg <- ggtree(node_tree, ladderize=T)
    for (i in 1:length(percentiles)){
        gg <- gg + geom_vline(xintercept = max(gg$data$x) * percentiles[i], color = "grey", size = 0.3)
    }
    gg <- gg + geom_nodelab(cex=2) 
    gg <- gg + geom_tiplab(cex=2) + xlim(0, max(gg$data$x) + 0.1*max(gg$data$x))
    ggsave(gg, file="./treeFigures/node_tree.pdf")
}# end treePlot

# Function to generate percentiles based upon the number of desired slices
percentilesFromSliceCount <- function(num_slices){
    return(seq((1/num_slices),1,(1/num_slices)))
}# end percentiles_from_slice_count

# Function to run previously build JAR file that performs actual clustering
runJar <- function(jar_path="PhyloPart_v2.1.jar",filename, percentile, bootstrap, perc_dist){
    # Run JAR in system call
    command <- paste("java -jar ", jar_path, " ", filename," ", perc_dist,
                     " -b", bootstrap, " -ophylopart", percentile, ".csv", sep='')
    write(command, stderr())
    sys_out <- try(system(command, intern=TRUE))
}# end run_jar

# Function to extend LETTERS for clusters > 26
extend <- function(alphabet) function(i) {
   base10toA <- function(n, A) {
      stopifnot(n >= 0L)
      N <- length(A)
      j <- n %/% N 
      if (j == 0L) A[n + 1L] else paste0(Recall(j - 1L, A), A[n %% N + 1L])
   }   
   vapply(i-1L, base10toA, character(1L), alphabet)
}# end LETTERS extension

# Function to generate alphabetical cluster IDs from unique IDs
getAlphaIds <- function(uniq_ids){
    num_uniq_ids <- length(uniq_ids)
    # Generate list for dictionary style access
    alpha_cluster_dict <- vector(mode="list", length=num_uniq_ids)
    MORELETTERS <- extend(LETTERS)

    # Generate alphabetical cluster IDs 
    # Set keys
    names(alpha_cluster_dict) <- as.character(uniq_ids)
    letter_range <- MORELETTERS(1:num_uniq_ids)
    for (x in 1:num_uniq_ids){
        # Set values
        alpha_cluster_dict[[x]] <- letter_range[x]
    }# end for
    return(alpha_cluster_dict)
}# end get_alpha_ids

getAlphaNumIds <- function(uniq_ids, tree_level){
    num_uniq_ids <- length(uniq_ids)
    alpha_cluster_dict <- vector(mode="list", length=num_uniq_ids)

    names(alpha_cluster_dict) <- as.character(uniq_ids)
    for (x in 1:num_uniq_ids){
        alpha_cluster_dict[[x]] <- as.character(x)
    }
    return(alpha_cluster_dict)
}# end getAlphaNumIds

expandNodes <- function(cluster_leaves, all_leaves, tree){
    clust_leaves <- c()
    cluster <- c()
    # Extend node labels into leaves
    for (leaf in 1:length(cluster_leaves$LeafName)){
        if (! as.character(cluster_leaves$LeafName[leaf]) %in% all_leaves){
            tips <- getDescendants(tree, as.character(cluster_leaves$LeafName[leaf]))
            keep <- c(length(tips))
            # Remove descendants that are nodes and keep leaves
            for (tip in 1:length(tips)){
                if (tips[tip] < length(all_leaves)){
                    keep[tip] <- TRUE
                } else {
                    keep[tip] <- FALSE
                }
            } # end inner for
            tips <- tips[as.logical(keep)]
            # Add all leaves (descendants) of nodes
            clust_leaves <- c(clust_leaves, tree$tip.label[tips])
            # For each newly add leaf -- add a cluster ID
            for (i in 1:length(tips)){
                cluster <- c(cluster, as.character(cluster_leaves$ClusterName[leaf]))
            }
        } else {
            clust_leaves <- c(clust_leaves, as.character(cluster_leaves$LeafName[leaf]))
            cluster <- c(cluster, as.character(cluster_leaves$ClusterName[leaf]))
        }# end if
    }# end outer for
    if (length(clust_leaves) != length(cluster)){
        write("Length mismatch in leaves/cluster IDs",stderr())
    }
    return(list("leaves"=clust_leaves, "clusters"=cluster))
}

getPaths <- function(final_long_df, final_short_df, tree){
    # For each row of the short Dataframe, determine the path at >= TreeLevel for each leaf 
    leaves <- tree$tip.label
    paths <- c(length(nrow(final_short_df)))
    clust_mrca_node <- list()
    clust_mrca_height <- list()
    for (x in 1:nrow(final_short_df)){
        row_num <- nrow(final_short_df)-(x-1) # start from last row (last cluster)
        # Get all leaves matching cluster ID
        cluster_leaves <- final_long_df[final_long_df$OldCluster == final_short_df[row_num,10],]
        # and where TreeLevel matches -- as cluster ID may be repeated in other cuts
        cluster_leaves <- cluster_leaves[cluster_leaves$TreeLevel == final_short_df[row_num,1],]
        # For nodes in cluster, get all downstream leaves
        cluster_leaves <- expandNodes(cluster_leaves, leaves, tree)
        clust_mrca_node[[row_num]] <- findMRCA(tree, as.character(cluster_leaves$leaves), type="node")
        clust_mrca_height[[row_num]] <- findMRCA(tree, as.character(cluster_leaves$leaves), type="height")
        # Get all leaves lower than current TreeLevel
        all_lower_leaves <- final_long_df[as.integer(final_long_df$TreeLevel) <= as.integer(final_short_df[row_num,1])-1,]
        if (nrow(all_lower_leaves) == 0){
            next
        }
        all_lower_leaves <- expandNodes(all_lower_leaves, leaves, tree)
        # Start indexing for local_paths (list) and initialize
        index = 1
        local_paths <- list()
        # For all rows in all_lower_leaves -- if lower leaves match current cluster leaves
        # add them to the local_paths list
        for (i in 1:length(all_lower_leaves$leaves)){
            if (as.character(all_lower_leaves$leaves[i]) %in% cluster_leaves$leaves){
                local_paths[index] <- all_lower_leaves$clusters[i]
                index = index + 1
            } 
        }
        local_paths[index] <- final_short_df[row_num, 2]
        local_paths <- as.character(local_paths)
        paths[row_num] <- paste(unique(local_paths), collapse=';')
    }# end outer for
    # Set Path column
    lowest_cut <- final_short_df[as.integer(final_short_df$TreeLevel) == min(as.integer(final_short_df$TreeLevel)),]
    paths[1:length(lowest_cut$ClusterName)] <- as.character(lowest_cut$ClusterName)
    final_short_df$Path <- paths
    final_short_df$MRCANode <- as.character(clust_mrca_node)
    final_short_df$MRCAHeight <- as.character(clust_mrca_height)
    # Reorder and cut unnecessary columns
    final_short_df <- final_short_df[c("TreeLevel","MedianOfDistances","ClusterName",
                           "SequencesPerCluster","Bootstrap","Path","PybusGamma","CutHeight",
                           "MRCANode","MRCAHeight","SpeciationRate")]
    return(final_short_df)
}# end getPath

getClusterPaths <- function(clusters){
    levels <- unique(clusters$TreeLevel)
    for (x in 1:length(levels)){
        rev_index <- length(levels)-(x-1)
        clusters_at_level <- clusters[clusters$TreeLevel == levels[rev_index],]$Path
        if (x == 1){
            # If first level analyzed, set all paths to current levels
            paths <- as.character(clusters_at_level)
        }
        if (x==length(levels)-1){
            return(paths)
        }
        for (i in 1:length(paths)){
            # Determine all paths in the next level
            next_level <- clusters[clusters$TreeLevel == levels[rev_index-1],]$Path
            # Get the split path
            split_string <- strsplit(as.character(paths[i]),';')
            # for existing paths, find preceeding cluster and append to path
            if (length(split_string[[1]]) >= 2){
                for (j in 1:length(next_level)){
                    if (identical(split_string[[1]][1], next_level[j])){
                        next_cluster <- strsplit(as.character(next_level[j]),';')[[1]][1]
                        if (!startsWith(paths[i], next_cluster)){
                            paths[i] <- paste(next_cluster, paths[i], sep=";")
                        }
                    }
                }# end j-for
                # for new/opening paths, add path to paths
                for (j in 1:length(next_level)){
                    first_occur <- strsplit(as.character(next_level[j]), ';')[[1]]
                    if (! is.na(first_occur[2])){
                        if (! any( grepl( first_occur[2], paths) ) ){
                            paths <- c(paths, as.character(next_level[j]))
                        }
                    } else if(is.na(first_occur[2])){
                        if (! any( grepl(first_occur[1], paths))){
                            paths <- c(paths, first_occur[1])
                        }
                    }
                }
            }# end if-else
        }# end i-for
    }# end x-for
}# end getClusterPaths

reduceClusterPaths <- function(paths, clusters, percentiles, tree){
    # Read all sliced trees
    trees <- list()
    for (i in 1:length(percentiles)){
        trees[[i]] <- read.newick(paste("./treeSlices/treeSlice",percentiles[i],".nwk",sep=""))
    }
    new_paths <- c(length(paths))
    cluster_mrca_dict <- vector(mode="list", length=length(clusters$MRCANode))
    names(cluster_mrca_dict) <- as.character(clusters$MRCANode)
    # Populate dictionary
    for (i in 1:length(clusters$ClusterName)){
        cluster_mrca_dict[[i]] <- as.character(clusters$ClusterName)[i]
    }
    # For each node in a path, determine the MRCAs and if they repeat
    for (i in 1:length(paths)){
        # For path in paths, get clusters
        clust_vector <- strsplit(paths[i],';')[[1]]
        if (length(clust_vector) > 1){
            temp_mrcas <- c(length(clust_vector))
            # For cluster in path, get MRCAs
            for (j in 1:length(clust_vector)){
                index <- match(clust_vector[j], clusters$ClusterName) 
                temp_mrcas[j] <- clusters$MRCANode[index]
            }
            # Get list of duplicated mrcas
            dup_mrcas <- unique(temp_mrcas[duplicated(temp_mrcas)]) # uniq set of duplicated mrcas
            # if no duplicates, pass path and go to next
            if (identical(integer(0), dup_mrcas)){
                new_paths[i] <- paths[i]
                next
            }
            bad_clusters <- c()
            for (j in 1:length(dup_mrcas)){
                # Determine indices of duplicated mrcas
                dup_boolean <- grepl(paste("^",dup_mrcas[j],"$",sep=""), temp_mrcas)
                # Get clusters from dup_mrca indices
                dup_clusters <- clust_vector[dup_boolean]
                # Get slice number for duplicated cluster
                slices <- c(length(dup_clusters))
                for (k in 1:length(dup_clusters)){
                    slices[k] <- as.integer(strsplit(strsplit(dup_clusters[k], '[.]')[[1]][1], 'c')[[1]][2])
                }
                # For mrca descendants in sliced tree, are descendants the same
                actual_dup_boolean <- c(length(slices))
                actual_dup_boolean[1] <- FALSE
                leaves <- getDescendants(tree, dup_mrcas[j])
                # Get leaves from mrca
                keep <- c(length(leaves))
                for (tip in 1:length(leaves)){
                    if (leaves[tip] < length(tree$tip.label)){
                        keep[tip] <- TRUE
                    } else {
                        keep[tip] <- FALSE
                    }
                } # end inner for
                leaves <- leaves[as.logical(keep)]
                leaves <- tree$tip.label[leaves]
                for (k in 2:length(slices)){
                    actual_dup_boolean[k] <- as.logical(all(leaves %in% trees[[slices[k-1]]]$tip.label) & all(leaves %in% trees[[slices[k]]]$tip.label))
                }
                bad_clusters <- c(bad_clusters, dup_clusters[as.logical(actual_dup_boolean)])
            }
            clean_path <- clust_vector[! clust_vector %in% bad_clusters]
        } else {
            clean_path <- paths[i]
        }
        new_paths[i] <- paste(clean_path, collapse = ";")
    }
    return(unique(new_paths))
}# end reduceClusterPaths

darken <- function(color, factor=1.4){
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue = 255)
    return(col)
}# end darken

plotClusters <- function(percentiles, node_tree_file="./treeFigures/node_tree.nwk", 
                         processed_tree_table="./treeTables/processedTree.csv",
                         width = 16, height = 16, units = "in", dpi = 600){
    tree <- read.newick(node_tree_file)
    clusters <- read.csv(processed_tree_table)
    paths <- getClusterPaths(clusters)
    paths <- reduceClusterPaths(paths, clusters, percentiles, tree)
    n <- length(paths)
    clust_count <- 0
    for (i in 1:n){
        clust_count <- clust_count + length(strsplit(paths[i],';')[[1]])
    }
    cols <- palette(rainbow(n))
    new_cols <- c(clust_count)
    all_node_paths <- c(clust_count)
    index <- 1
    for (i in 1:n){
        each_node <- strsplit(paths[i],';')[[1]]
        for (j in 1:length(each_node)){
            all_node_paths[index] <- each_node[j]
            if (j == 1){
                new_cols[index] <- darken(cols[i])
                index <- index + 1
            } else {
                new_cols[index] <- darken(new_cols[index-1])
                index <- index + 1
            }
        }
    }
    gg <- ggtree(tree, ladderize=T)
    clust_vector <- as.character(clusters$MRCANode)
    clust_height <- as.character(clusters$MRCAHeight)
    highlighted_clusters <- c()
    index <- 1
    # For path in paths
    for (i in 1:n){
        # For cluster in path
        for (j in 1:length(strsplit(paths[i],';')[[1]])){
            clust_vector <- clusters[as.character(clusters$ClusterName) == strsplit(paths[i],';')[[1]][j],]
            if (clust_vector$MRCANode %in% highlighted_clusters){
                index <- index + 1
                next
            }
            gg <- gg + geom_hilight(clust_vector$MRCANode, fill = new_cols[index], alpha = 0.2)
            highlighted_clusters <- c(highlighted_clusters, clust_vector$MRCANode)
            index <- index + 1
        }
    }
    repeat_node <- c(n)
    index <- 1
    for (i in 1:n){
        nodes <- strsplit(paths[i],';')[[1]]
        last_node <- nodes[length(nodes)]
        clust_vector <- clusters[as.character(clusters$ClusterName) == last_node,]
        if (any(grepl(pattern = clust_vector$MRCANode, x = repeat_node))){
            gg <- gg + geom_cladelabel(node = clust_vector$MRCANode, label = paths[i],
                                       align = T, fontsize = 2, offset.text = max(gg$data$x)*0.1)
        } else{
            gg <- gg + geom_cladelabel(node = clust_vector$MRCANode, label = paths[i],
                                       align = T, fontsize = 2)
        }
        repeat_node[i] <- clust_vector$MRCANode
    }
    
    for (i in 1:length(percentiles)){
        gg <- gg + geom_vline(xintercept = max(gg$data$x) * percentiles[i], color = "grey", size = 0.3)
    }
    gg <- gg +  geom_nodelab(cex=2) + geom_tiplab(cex=2) + ggplot2::xlim(0, max(gg$data$x) + 0.15*max(gg$data$x))
    ggsave(gg, file="./treeFigures/cluster_tree.pdf", dpi = dpi, width = width, height = height, units = units)
    # Generate paths reflecting reduction
    remaining_nodes <- c(sum(lengths(regmatches(paths, gregexpr("c", paths)))))
    for (i in 1:length(paths)){
        remaining_nodes <- c(remaining_nodes, strsplit(paths[i], ";")[[1]]) 
    }
    repeat_column <- c(length(clusters$ClusterName))
    old_paths <- as.character(clusters$Path)
    for (i in 1:length(clusters$Path)){
        if (clusters$ClusterName[i] %in% remaining_nodes){
            repeat_column[i] <- "-"
        } else {
            split_path <- strsplit(old_paths[i], ";")[[1]]
            if (length(split_path) > 1){
                repeat_column[i] <- split_path[length(split_path)-1]
            } else {
                repeat_column[i] <- old_paths[i]
            }
        }
    }
    clusters$RepeatCluster <- repeat_column
    return(clusters)
}

generateLeafColumn <- function(percentiles, final_clusters){
    tree_levels <- as.integer(final_clusters$TreeLevel)
    start_level <- min(tree_levels)
    stop_level <- max(tree_levels)
    all_leaves <- c(length(final_clusters$TreeLevel))
    index = 1
    for (i in start_level:stop_level){
        clusters <- as.character(final_clusters[final_clusters$TreeLevel == i,]$ClusterName)
        filename <- paste("./treeSlices/phylopart",percentiles[i],".csv", sep="")
        temp_df <- read.csv(filename, header=TRUE)
        for (j in 1:length(clusters)){
            all_leaves[index] <- paste(sort(as.character(temp_df[temp_df$clustername == clusters[j],]$leafname)), sep="", collapse=";")
            index <- index + 1
        }
    }
    final_clusters$Leaves <- all_leaves
    return(final_clusters)
}

processTree <- function(input_tree, slice_count=10, bootstrap=0.70, min_leaves=15, perc_dist=0.05){
    # Get percentiles from slice count
    percentiles <- percentilesFromSliceCount(slice_count)
    # Check tree type and covert if necessary
    filename <- checkAndConvertTree(input_tree, percentiles)
    # Check tree is rooted
    tree <- read.newick(filename)
    if (!is.rooted(tree)){
        write("Tree not rooted. Must provide rooted tree.", stderr())
        quit(status=1)
    }
    # Check tree has node labels
    if (is.null(tree$node.label)){
        tree$node.label <- rep_len(1, tree$Nnode)
    }    
    # Resolve Multichotomies
    tree <- multi2di(tree, random = FALSE)
    # Make directory for treeSlices
    dir.create(file.path(".","treeSlices"))
    dir.create(file.path(".","treeFigures"))
    dir.create(file.path(".","treeTables"))
    # Plot main tree and node tree
    treePlot(tree, percentiles)
    # Check tree height
    height <- max(nodeHeights(tree))
    # Initialize list of quality DataFrames
    quality_percentiles <- list()
    # Plot Lineage-through-time plot
    pdf(paste("./treeFigures/",filename, ".LTT.pdf", sep=""))
    ltt(tree, gamma=TRUE)
    dev.off()
    
    phylogenetic_diversity = c()
    # Check JAR system call using percentiles -- may need to alter JAR path
    for (x in 1:length(percentiles)){
        # Cut the tree at the designated percentile, keep rootward, write new file
        cut_tree <- treeSlice(tree, height*percentiles[x], orientation="rootwards")
        # Get the index for each nodel label from tree edge -- use index to get tree node label
        new_labels = list()
        unique_edges = unique( tree$edge[,1])
        for (i in 1:length(cut_tree$node.label)){
            index <- match(cut_tree$node.label[i],unique_edges)
            new_labels[[i]] <- tree$node.label[index]
        }
        cut_tree$node.label <- as.character(new_labels)
        cut_tree_file <- paste("./treeSlices/treeSlice",percentiles[x],".nwk", sep="")

        write.tree(cut_tree, file=cut_tree_file, tree.names=TRUE)
        # Calculate speciation rate
        Speciation_rate <- yule(cut_tree)
        # Calculate Pybus gamma without plot
        Pybus_gamma <- ltt(cut_tree, plot=FALSE, gamma=TRUE)
        # Run jar
        runJar("PhyloPart_v2.1.jar", filename=cut_tree_file, percentile=percentiles[x], 
               bootstrap=bootstrap, perc_dist=perc_dist)
        # Read output into DataFrame
        temp_file <- paste("phylopart",percentiles[x],".csv", sep='')
        df <- read.csv(temp_file)
        old_tips <- cut_tree$tip.label
        # Need to rename leaves that aren't full length
        for (i in 1:length(cut_tree$tip.label)){
            if (cut_tree$tip.label[i] %in% tree$tip.label){                
                index <- match(i, cut_tree$edge[,2])
                old_tip_index <- match(cut_tree$tip.label[i], tree$tip.label)
                old_edge_index <- match(old_tip_index, tree$edge[,2])
                if ( abs(cut_tree$edge.length[index] - tree$edge.length[old_edge_index]) >= 1e-5 ){
                    new_id <- paste(cut_tree$tip.label[i],"$",format(cut_tree$edge.length[index],digits=8), 
                                    "-", getParent(tree, match(cut_tree$tip.label[i], tree$tip.label)), sep="")
                    cut_tree$tip.label[i] <- new_id
                }
            }
        }
        # Plot the cut trees
        gg <- ggtree(cut_tree, ladderize=T)
        gg <- gg +  geom_nodelab(cex=2) + geom_tiplab(cex=2) + xlim(0, max(gg$data$x) + 0.15*max(gg$data$x))
        ggsave(gg, file=paste("./treeSlices/treeSlice",percentiles[x],".pdf", sep=""))
        # Move phylopart files to treeSlices
        file.remove(temp_file)
        write.tree(cut_tree, file=cut_tree_file, tree.names=TRUE)
        if (nrow(df) < min_leaves){
            write(paste(percentiles[x]*100,
                        "th percentile has too few intersections along the cut",
                        sep=''),
                  stderr())
            next
        }
        # Filter unclustered leaves
        df <- df[ df$clustername != 0,]
        # Make sure DataFrame isn't empty. If it is, go to next percentile
        if (nrow(df) == 0){
            write(paste(percentiles[x]*100,
                        "th percentile has zero leaves after filtering.", 
                        sep=''),
                  stderr())
            next
        }# end if
        tip_table <- data.frame(old_tips, cut_tree$tip.label)
        write.csv(tip_table, paste("./treeSlices/renamed_leaves",percentiles[x],".csv",sep=""), row.names = FALSE)
        new_leaves <- c(length(df$leafname))
        for (i in 1:length(df$leafname)){
            new_leaf_index <- grepl(paste("^",df$leafname[i], "$", sep=""), tip_table$cut_tree.tip.label)
            new_leaves[i] <- as.character(tip_table$cut_tree.tip.label)[as.logical(new_leaf_index)][1]
        }
        new_df <- df
        new_df$leafname <- new_leaves
        # Identify unique cluster IDs
        uniq_ids <- unique(df$clustername)
        # Convert random ID numbers to alpha characters
        alpha_ids <- getAlphaNumIds(uniq_ids)
        # Assign alpha_ids over numeric
        new_ids <- list()
        for (i in 1:length(df$clustername)){
            clust_name <- as.character(df$clustername[i])
            new_ids[[i]] <- paste("c",x,".",alpha_ids[[clust_name]], sep='')
        }# end for
        new_df$clustername <- as.character(unlist(new_ids))
        write.csv(new_df, paste("./treeSlices/",temp_file,sep=""), row.names=FALSE)
        df$oldcluster <- df$clustername
        df$clustername <- as.character(unlist(new_ids))
        df$TreeLevel <- rep(as.character(x),times=length(df$clustername))
        df$PybusGamma <- rep(as.character(Pybus_gamma$gamma), times=length(df$clustername))
        df$CutHeight <- rep(as.character(height*percentiles[x]), times=length(df$clustername))
        df$SpeciationRate <- rep(as.character(Speciation_rate$lambda), times=length(df$clustername))
        # Rename Dataframe headers
        # clustername,bootstrap,leafname,branchPath,medianOfDistances,sequencesperCluster,TreeLevel,Pybus_gamma
        names(df) <- c("ClusterName","Bootstrap","LeafName","RootToTipDist",
                       "MedianOfDistances","SequencesPerCluster","OldCluster",
                       "TreeLevel","PybusGamma","CutHeight","SpeciationRate")
        # Reorder Dataframe
        df <- df[c("TreeLevel","ClusterName","SequencesPerCluster",
                   "Bootstrap","LeafName","RootToTipDist","MedianOfDistances",
                   "PybusGamma","CutHeight","OldCluster","SpeciationRate")]
        # For cluster in DataFrame, calculate phylogenetic diversity
        uniq_ids <- unique(df$ClusterName)
        community_matrix <- matrix(data=rep(0,length(uniq_ids)*length(cut_tree$tip.label)), 
                                   nrow=length(uniq_ids), ncol=length(cut_tree$tip.label))
        rownames(community_matrix) <- uniq_ids
        colnames(community_matrix) <- cut_tree$tip.label
        for (i in 1:length(unique(uniq_ids))){
            cluster_leaves <- df[df$ClusterName == uniq_ids[i],]$LeafName
            for (j in 1:length(cut_tree$tip.label)){
                for (k in 1:length(cluster_leaves)){
                    if (identical(cluster_leaves[k], cut_tree$tip.label[j])){
                        community_matrix[i,j] <- 1
                    }
                }
            }
        }
        phylogenetic_diversity <- c(phylogenetic_diversity, pd(community_matrix, cut_tree, include.root=TRUE)$PD)
        quality_percentiles[[as.character(percentiles[x])]] = df
    }# end for
    # Merge Dataframes in list into single large Dataframe
    final_long_df <- Reduce(function(df1, df2) rbind(df1, df2),quality_percentiles)
    write.csv(final_long_df, "./treeTables/long_format.csv", row.names=FALSE)
    # Reindex final_long_df
    rownames(final_long_df) <- 1:nrow(final_long_df)
    # Filter Dataframe using rows with unique information for specified columns
    final_short_df <- final_long_df[!duplicated(final_long_df[,c("TreeLevel","ClusterName",
                                        "SequencesPerCluster","Bootstrap","MedianOfDistances","PybusGamma","CutHeight")]),]
    # Get paths for each row and filter columns
    final_short_df <- getPaths(final_long_df, final_short_df, tree)
    final_short_df$PhyloDiversity <- phylogenetic_diversity
    rownames(final_short_df) <- 1:nrow(final_short_df)
    write.csv(final_short_df, file="./treeTables/processedTree.csv", row.names=FALSE)
    #file.rename(from="final_paths.csv", "./treeTables/final_paths.csv")  <-- consider recursively removing
    # Generate plot of clusters
    final_clusters <- plotClusters(percentiles)
    final_clusters <- generateLeafColumn(percentiles, final_clusters)
    write.csv(final_clusters, file="./treeTables/HIVdynamite.csv", row.names=FALSE)
    file.remove("./treeTables/processedTree.csv")
    file.remove("./treeTables/long_format.csv")
    return(final_clusters)
}# end function processTree

testProcessTree <- function(){
    input_tree <- "losalamos_small.nwk"
    slice_count <- 10
    bootstrap <- 0.90
    min_leaves <- 15
    perc_dist <- 0.05

    output <- processTree(input_tree, slice_count, bootstrap, min_leaves, perc_dist)
    #percentiles <- percentilesFromSliceCount(slice_count)
    #print(plotClusters(percentiles))
}
