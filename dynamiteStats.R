birthCountClusteredLeafCount <- function(df) {
  # Given data.frame, returns birth count and other info
  births <- c()
  clust_leaves <- c()
  median_seq_per_clust <- c()
  median_pd <- c()
  for (x in 0:9) {
    birth_count <- 0
    for (i in df[df$TreeLevel == x + 1, ]$Path) {
      if (length(strsplit(i, ';')[[1]]) == 1) {
        birth_count <- birth_count + 1
      }
    }
    births <- c(births, birth_count)
    clust_leaves <- c(clust_leaves,
                      sum(df[df$TreeLevel == x + 1, ]$SequencesPerCluster))
    med <- median(df[df$TreeLevel == x + 1, ]$SequencesPerCluster)
    if (length(med) == 0) {
      median_seq_per_clust <- c(median_seq_per_clust, 0)
    }  else {
      median_seq_per_clust <- c(median_seq_per_clust, med)
    }
    med_pd <- median(df[df$TreeLevel == x + 1, ]$PhyloDiversity)
    if (length(med_pd) == 0) {
      median_pd <- c(median_pd, 0)
    } else {
      median_pd <- c(median_pd, med_pd)
    }
  }
  results_list <- list(
    "Births" = births,
    "ClustLeaves" = clust_leaves,
    "MedianSeqPerClust" = median_seq_per_clust,
    "MedianPD" = median_pd
  )
  return(results_list)
}


deathCount <- function(df) {
  # Returns death count
  clusters <- df$ClusterName
  lvl10 <- df[df$TreeLevel == 10, ]
  dead_clusters = c()
  for (cluster in clusters) {
    cluster_in <- c()
    for (i in lvl10$Path) {
      if (cluster %in% strsplit(i, ';')[[1]]) {
        cluster_in <- c(cluster_in, TRUE)
      } else {
        cluster_in <- c(cluster_in, FALSE)
      }
    }
    if (!any(cluster_in)) {
      dead_clusters <- c(dead_clusters, cluster)
    }
  }
  # Determine where clusters die -- for TreeLevel count
  death_counts = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  dead_path = c()
  for (cluster in dead_clusters) {
    for (x in 0:8) {
      lvl = df[df$TreeLevel == x + 2, ]
      
      cluster_level = as.integer(strsplit(gsub("c", "", cluster), "\\.")[[1]][1])
      path_in <- c()
      for (i in lvl$Path) {
        if (cluster %in% strsplit(i, ';')[[1]]) {
          path_in <- c(path_in, TRUE)
        } else {
          path_in <- c(path_in, FALSE)
        }
      }
      if (((x + 2) > cluster_level) & !any(path_in)) {
        # if cluster is also not in a path already identified as dead
        cluster_in_path <- c()
        for (i in dead_path) {
          if (cluster %in% strsplit(i, ';')[[1]]) {
            cluster_in_path <- c(cluster_in_path, TRUE)
          } else {
            cluster_in_path <- c(cluster_in_path, FALSE)
          }
        }
        if (!any(cluster_in_path)) {
          # for any paths from the previous level containing that cluster
          paths <- c()
          for (i in df[df$TreeLevel == x + 1, ]$Path) {
            if (cluster %in% strsplit(i, ';')[[1]]) {
              paths <- c(paths, i)
            }
          }
          # append them
          dead_path <- c(dead_path, paths)
          # count the dead
          death_counts[x + 2] <- death_counts[x + 2] + length(paths)
        }
      }
    }
  }
  return(death_counts)
}

singletonCount <- function(df) {
  # Returns singleton count
  singletons <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  # only clusters not in the last level may be singletons
  lvlNot10 <- df[df$TreeLevel != 10, ]
  clusters <- lvlNot10$ClusterName
  # for all clusters, if the cluster is in only 1 path and that path has 1 cluster -- it's a singleton
  for (cluster in clusters) {
    paths_with_cluster <- c()
    for (i in df$Path) {
      if (cluster %in% strsplit(i, ';')[[1]]) {
        paths_with_cluster <- c(paths_with_cluster, TRUE)
      } else {
        paths_with_cluster <- c(paths_with_cluster, FALSE)
      }
    }
    if (sum(paths_with_cluster) == 1) {
      hit_single = as.character(df$Path[paths_with_cluster])
      #hit_single.index = c(0)
      level = df$TreeLevel[paths_with_cluster][1]
      if (length(strsplit(hit_single[1], ';')) == 2) {
        singletons[level] = singletons[level] + 1
      }
    }
  }
  return(singletons)
}


getPybusSpeciationClustcount <- function(df) {
  # Process Pybus gamma, speciationRate and cluster count from dataframe
  df2 = df[, c("TreeLevel", "SpeciationRate", "PybusGamma")]
  df2 = unique(df2)
  df2$ClusterCount <- sort(table(df$TreeLevel))
  #df2$TreeNum <- rep(tree_count, nrow(df2))
  #df2$LeafCount <- rep(leaf_count, nrow(df2))
  return(df2)
}


getBirthDeathSingleClust <- function(df) {
  # births, clust_leaves, median_seq_per_clust, median_pd
  bcclc <- birthCountClusteredLeafCount(df)
  bcclc$Deaths <- deathCount(df)
  bcclc$Singletons <- singletonCount(df)
  return(data.frame(bcclc))
}


fillEmptyTreeLevels <- function(df) {
  # If there are no clusters on a tree level, fill empty values
  for (x in 0:9) {
    lvl = df[df$TreeLevel == x + 1, ]
    if (nrow(lvl) == 0) {
      # df = rbind(df, c(x+1,NA,NA,0,tree_count, leaf_count))
      df <- rbind(df, c(x + 1, NA, NA, 0))
    }
  }
  df <- df[order(df$TreeLevel), ]
  rownames(df) <- 1:nrow(df)
  return(df)
}


totalLeafCount <- function(slices) {
  counts = c()
  for (i in slices) {
    subtree_leaf_count <-
      ape::read.tree(paste("./treeSlices/treeSlice", i, ".nwk", sep = ""))$tip.label
    counts <- c(counts, subtree_leaf_count)
  }
  return(list("TotalLeaves" = counts))
}


clusterSizeByTreeSize <- function(df) {
  # Get data form ClusterSize/TreeSize boxplot with tertile
  df <- df[, c("TreeLevel",
               "SequencesPerCluster",
               "PhyloDiversity",
               "MedianOfDistances")]
  #df$TreeCount <- rep(tree_count, nrow(df))
  #df$LeafCount <- rep(leaf_count, nrow(df))
  return(df)
}

# totalBranchLengthFromRoot <- function(){
#   # For each tree slice, read the sum(branchLengthFromRoot)
#   # and add to appropriate slice in df
#   slices <- c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1")
#   df <- data.frame()
#   for (level in slices){
#     df <- rbind(df, read.csv(paste("/treeSlices/treeSliceTotalBranchLength",
#                                   level, ".csv", sep="")))
#   }
#   colnames(df) <- c("TotalBranchLength")
#   return(df$TotalBranchLength)
# }


rootTipDist <- function() {
  scores <- c()
  for (j in c(seq(0.1, 0.9, 0.1), "1")) {
    filename <- paste("./treeSlices/treeSlice",
                      j, ".nwk", sep = "")
    tree <- ape::read.tree(filename)
    tree$node.label <- seq(1, length(tree$node.label))
    scores <- c(scores, sum(adephylo::distRoot(tree)))
  }
  return(scores)
}

dynaStats <- function() {
  df <- read.csv("./treeTables/HIVdynamite.csv")
  
  cluster_level_df <- clusterSizeByTreeSize(df)
  
  df_psc <- getPybusSpeciationClustcount(df)
  df_psc <- fillEmptyTreeLevels(df_psc)
  mid <- getBirthDeathSingleClust(df)
  df_psc <- merge(df_psc, mid, by = "row.names", all = TRUE)
  df_psc$Row.names <- NULL
  df_psc$TotalBranchLength <- rootTipDist()
  df_psc <- df_psc[order(df_psc$TreeLevel), ]
  
  stability_set <- c(0)
  shrink_set <- c(0)
  growth_set <- c(0)
  founder_set <- c(0)
  
  for (i in 1:9) {
    shrink <- 0
    growth <- 0
    stability <- 0
    founder <- 0
    previous_clust_name = as.character(df[df$TreeLevel == i + 1, ]$ClusterName)
    previous_seqs_per_clust = as.character(df[df$TreeLevel == i + 1, ]$SequencesPerCluster)
    current_path = as.character(df[df$TreeLevel == i + 2, ]$Path)
    current_seqs_per_clust = data.frame(df[df$TreeLevel == i + 2, ]$SequencesPerCluster)
    colnames(current_seqs_per_clust) <- c("SequencesPerCluster")
    
    if (length(previous_clust_name) == 0) {
      stability_set <- c(stability_set, stability)
      shrink_set <- c(shrink_set, shrink)
      growth_set <- c(growth_set, growth)
      founder_set <- c(founder_set, founder)
      next
    }
    current_path_with_prev_clust = c()
    
    used <- c()
    for (m in current_path) {
      used <- c(used, F)
    }
    for (x in 1:length(previous_clust_name)) {
      contained_indices <- c()
      for (k in current_path) {
        if (previous_clust_name[x] %in% strsplit(k, ';')[[1]]) {
          contained_indices <- c(contained_indices, TRUE)
        } else {
          FALSE
        }
      }
      
      not_used <- c()
      for (m in used) {
        not_used <- c(not_used,!m)
      }
      contained_indices <- contained_indices & not_used
      contained_counts <- sum(contained_indices)
      if (any(contained_indices)) {
        #rownames(previous_seqs_per_clust) <- 1:length(previous_seqs_per_clust)
        
        if (contained_counts > 1) {
          founder = founder + contained_counts - 1
          for (clust_count in current_seqs_per_clust[contained_indices, ]) {
            if (previous_seqs_per_clust[x] == clust_count) {
              stability <- stability + 1
            } else if (previous_seqs_per_clust[x] < clust_count) {
              growth <- growth + 1
            } else if (previous_seqs_per_clust[x] > clust_count) {
              shrink <- shrink + 1
            }
          }
        } else if (contained_counts == 1) {
          if (previous_seqs_per_clust[x] < sum(current_seqs_per_clust[contained_indices, ])) {
            growth <- growth + 1
          } else if (previous_seqs_per_clust[x] > sum(current_seqs_per_clust[contained_indices, ])) {
            shrink <- shrink + 1
          } else if (previous_seqs_per_clust[x] == sum(current_seqs_per_clust[contained_indices, ])) {
            stability <- stability + 1
          }
        }
      }
      #useable = pd.Series([not m for m in contained_indices])
      used <- contained_indices | used
    }
    stability_set <- c(stability_set, stability)
    shrink_set <- c(shrink_set, shrink)
    growth_set <- c(growth_set, growth)
    founder_set <- c(founder_set, founder)
  }
  stability_set <- list("Stability" = stability_set)
  shrink_set <- list("Shrink" = shrink_set)
  growth_set <- list("Growth" = growth_set)
  founder_set <- list("Founder" = founder_set)
  more <-
    data.frame(c(stability_set, shrink_set, growth_set, founder_set))
  
  df_psc <- merge(df_psc, more, by = "row.names", all = TRUE)
  df_psc$Row.names <- NULL
  return(df_psc[order(df_psc$TreeLevel),])
}
