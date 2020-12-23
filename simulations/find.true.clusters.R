find.true.clusters <- function(family_tree, clades) {
  true_clusters <- list()
  for (clade in seq_along(clades)) { # subtrees (and whole tree) dataframes correspond to output dataframes from get_children function in familyR (see below)
    anc <- tail(familyR::get_parents(family_tree, clades[[clade]]$from[1], return_nodes=T)[[2]]$label, n=1)
    for (state in cluster_dynamics$state[cluster_dynamics$state %notin% c('A')]) { # states correspond to any letter that is not A, which is background
      if ((length(grep(state, clades[[clade]]$label)) >= 0.95*nrow(clades[[clade]])) & # If 90% of subtree is comprised of one particular state
          (length(grep(state, clades[[clade]]$label[1])) == 1) & # And if the first child in the subtree belongs to that state
          (length(grep(state, tail(familyR::get_parents(family_tree, clades[[clade]]$from[1],
                                                        return_nodes=T)[[2]]$label, n=2))) == 1)) { # AND if parent (indexed from whole tree) of that parent is not that state
        true_clusters[[clade]] <- data.frame(cluster_id = paste0(state,clade), # This will be some arbitrary name you assign to the cluster, such as C1, C2, C3, etc.
                                             taxa = toString(c(anc, clades[[clade]]$label)), # list of taxa
                                             mean_R0 = cluster_dynamics$mean_R0[cluster_dynamics$state == state],
                                             upper_R0 = cluster_dynamics$upper_R0[cluster_dynamics$state == state],
                                             lower_R0 = cluster_dynamics$lower_R0[cluster_dynamics$state == state])
                                             #dynamic = cluster_dynamics$dynamic[cluster_dynamics$state == state],
                                             #birth = cluster_dynamics$birth[cluster_dynamics$state == state]) # dynamic corresponds to the dynamic state (growth, etc.) assigned to the cluster state
      } else{NULL}
    } # End for states DNE A (the code below is for when births are being tested (e.g., D -> E))
    ## Now need to account for D (beause not going to be 95%)
    # state <- cluster_dynamics$state[cluster_dynamics$state %in% c('D','E')]
    #   if ((length(grep(state[1], clades[[clade]]$label)) +
    #        length(grep(state[2], clades[[clade]]$label)) >= 0.95*nrow(clades[[clade]])) &
    #       (length(grep(state[1], clades[[clade]]$label[1])) == 1) &
    #       (length(grep(state[1], tail(familyR::get_parents(family_tree, clades[[clade]]$from[1],
    #                                                        return_nodes=T)[[2]]$label, n=2))) == 1)) {
    #     true_clusters[[clade]] <- data.frame(cluster_id = paste0(state[1],clade), # This will be some arbitrary name you assign to the cluster, such as C1, C2, C3, etc.
    #                                          #taxa = toString(clades[[clade]]$label[grepl(state[1], clades[[clade]]$label)]), # list of taxa
    #                                          taxa = toString(c(anc, clades[[clade]]$label)), # list of taxa
    #                                          mean_R0 = cluster_dynamics$mean_R0[cluster_dynamics$state == state[1]],
    #                                          upper_R0 = cluster_dynamics$upper_R0[cluster_dynamics$state == state[1]],
    #                                          lower_R0 = cluster_dynamics$lower_R0[cluster_dynamics$state == state[1]],
    #                                          #dynamic = cluster_dynamics$dynamic[cluster_dynamics$state == state[1]],
    #                                          birth = cluster_dynamics$birth[cluster_dynamics$state == state[1]]) # dynamic corresponds to the dynamic state (growth, etc.) assigned to the cluster state
    #   } else{NULL} # End for state 'D'
  # state <- "E"
  # if ((length(grep(state, clades[[clade]]$label)) >= 0.95*nrow(clades[[clade]])) &
  #     (length(grep(state, clades[[clade]]$label[1])) == 1) &
  #     (length(grep(state, tail(familyR::get_parents(family_tree, clades[[clade]]$from[1],
  #                                                   return_nodes=T)[[2]]$label, n=2))) == 1)) {
  #   true_clusters[[clade]] <- data.frame(cluster_id = paste0(state,clade), # This will be some arbitrary name you assign to the cluster, such as C1, C2, C3, etc.
  #                                        taxa = toString(c(anc, clades[[clade]]$label)), # list of taxa
  #                                        mean_R0 = cluster_dynamics$mean_R0[cluster_dynamics$state == state],
  #                                        upper_R0 = cluster_dynamics$upper_R0[cluster_dynamics$state == state],
  #                                        lower_R0 = cluster_dynamics$lower_R0[cluster_dynamics$state == state],
  #                                        #dynamic = cluster_dynamics$dynamic[cluster_dynamics$state == state],
  #                                        birth = cluster_dynamics$birth[cluster_dynamics$state == state]) # dynamic corresponds to the dynamic state (growth, etc.) assigned to the cluster state
  # } else{NULL} # End for state 'E'
} # End loop along clades
  return(Filter(Negate(function(y) is.null(unlist(y))), true_clusters))
} # End function
merge.overlap.clust.true <- function(true_clusters) {
  copy <- true_clusters
  unwanted <- list()
  unwanted2 <- list()
  for (ct in seq_along(true_clusters)) { # Do not want to merge cluster 'F' with origin cluster 'E'
    true_clusters[[ct]]$taxa <- gsub(" ", "", true_clusters[[ct]]$taxa) # Need to get rid of spaces
    for (cc in seq_along(copy)) {
      copy[[cc]]$taxa <- gsub(" ", "", copy[[cc]]$taxa) # Need to get rid of spaces
      # if (is.na(true_clusters[[ct]]$birth) &
      #     is.na(copy[[cc]]$birth) &
      #     isTRUE(all(unlist(str_split(true_clusters[[ct]]$taxa, ',')) %in% unlist(str_split(copy[[cc]]$taxa, ','))) &
      #            true_clusters[[ct]]$cluster_id != copy[[cc]]$cluster_id)) {
      #   unwanted[[ct]] <- true_clusters[[ct]]
      # } else {NULL}
      # if (!is.na(true_clusters[[ct]]$birth) &
      #     !is.na(copy[[cc]]$birth) &
      #     isTRUE(all(unlist(str_split(true_clusters[[ct]]$taxa, ',')) %in% unlist(str_split(copy[[cc]]$taxa, ','))) &
      #            true_clusters[[ct]]$cluster_id != copy[[cc]]$cluster_id)) {
      #   unwanted2[[ct]] <- true_clusters[[ct]]
      # } else {NULL}
    } # End loop along copy
  } # End loop along true
  result <- setdiff(true_clusters, unwanted)
  result <- setdiff(result, unwanted2)
  return(result)
}