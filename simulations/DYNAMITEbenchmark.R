
## Initialize stable working environment and store time of intiation
rm(list=ls())



# List of packages for session
.packages <-  c("optparse", "remotes", "phytools", "treeio", "dplyr", "plyr", "tidyr", "tidytree", "data.table", 
                "parallel", "stringr", "rlist", "paleotree", "phylodyn",
                "ggtree", "drc", "growthrates", "TreeTools", "geiger") # May need to incorporate code for familyR (https://rdrr.io/github/emillykkejensen/familyR/src/R/get_children.R) i fno longer supported.
.github_packages <- c("emillykkejensen/familyR", "mrc-ide/skygrowth") 

# # Install CRAN packages (if not already installed)
# .inst <- .packages %in% installed.packages()
# if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# .inst_github <- .packages %in% installed.packages()
# ## Install GitHub packages(if not already installed)
# if(length(github_packages[!.inst_github]) > 0) try(remotes::install_github(github_packages[!.inst_github]))
# if(length(github_packages[!.inst_github]) > 0) try(devtools::install_github(github_packages[!.inst_github]))

# Load packages into session 
lapply(.packages, require, character.only=TRUE)
lapply(gsub(".+\\/(.+)", "\\1", .github_packages), require, character.only=TRUE)
numCores <- detectCores()

option_list = list(
  make_option(c("-s", "--sim_index"), type="numeric"),
  make_option(c("-a", "--cluster"), type="character", default="b", 
              help="choice of cluster algorithm from c (Phylopart's cladewise) or b (DYNAMITE's branchwise) [default= phylopart]", metavar="character"),
  make_option(c("-l", "--leaves"), type="character", default="", 
              help="choice of transformation to tree from bifurcating or addLeaves [default=empty]", metavar="character"),
  make_option(c("-t", "--threshold"), type="numeric", default=0.05, 
              help="branch length threshold [default= 0.05]", metavar="numeric")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$s)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

### Functions ##############################################################################################################

`%notin%` <- Negate(`%in%`) # Just plain useful

true.cluster.dyn <- function(conf.level=0.95){
  cluster_dynamics <- data.frame(state=rbind('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'),
                                 dynamic = rbind('background', 'static', 'static', 'static', 'static', 'growth', 'decay', 'static', 'static'),
                                 birth = rbind(NA,NA,NA,NA,NA,NA,NA,NA, 'birth')
                                 #death = rbind(NA,NA,NA,NA,NA), stringsAsFactors = F
                                 )
  
  ## Estimate R0 for each cluster based on paramaters used in simulation
  s <- 100 # sample size of 100
  n_contacts_A <- rnorm(s, 16, 1)
  n_contacts_BC <- rnorm(s, 4, 1)
  n_contacts_D <- rnorm(s, 6, 1)
  n_contacts_EHI <- rnorm(s, 4, 1)
  n_contacts_F <- as.numeric(2)
  n_contacts_G <- as.numeric(6)
 
  t_incub <- rnorm(s, 5, 2)
  t_exit <- rnorm(s, 14, 2)
  p_trans_A <- 0.015
  p_trans_BD <- 0.1
  p_trans_CFGHI <- 0.15
  p_trans_E <- 0.2

  
  Z=qnorm(0.5*(1 + conf.level))
  
  R0_A = sample(n_contacts_A*p_trans_A*(t_exit-t_incub),
                size=s, replace=F)
  R0_B = sample(n_contacts_BC*p_trans_BD*(t_exit-t_incub),
                size=s, replace=F)
  R0_C = sample(n_contacts_BC*p_trans_CFGHI*(t_exit-t_incub),
                  size=s, replace=F)
  R0_D = sample(n_contacts_D*p_trans_BD*(t_exit-t_incub),
                size=s, replace=F)
  R0_E = sample(n_contacts_EHI*p_trans_E*(t_exit-t_incub),
                 size=s, replace=F)
  R0_F = sample(n_contacts_F*p_trans_CFGHI*(t_exit-t_incub),
                size=s, replace=F)
  R0_G = sample(n_contacts_G*p_trans_CFGHI*(t_exit-t_incub),
                size=s, replace=F)
  R0_H = sample(n_contacts_EHI*p_trans_CFGHI*(t_exit-t_incub),
                size=s, replace=F)
  R0_I = sample(n_contacts_EHI*p_trans_CFGHI*(t_exit-t_incub),
                size=s, replace=F)
  
  logR0_A = log(R0_A)
  logR0_B = log(R0_B)
  logR0_C = log(R0_C)
  logR0_D = log(R0_D)
  logR0_E = log(R0_E)
  logR0_F = log(R0_F)
  logR0_G = log(R0_G)
  logR0_H = log(R0_H)
  logR0_I = log(R0_I)
  
  SElogR0_A = sd(logR0_A, na.rm=T)/sqrt(s) # standard deviation of 2, sample size of 10
  SElogR0_B = sd(logR0_B, na.rm=T)/sqrt(s) # standard deviation of 2, sample size of 10
  SElogR0_C = sd(logR0_C, na.rm=T)/sqrt(s) # standard deviation of 2, sample size of 10
  SElogR0_D = sd(logR0_D, na.rm=T)/sqrt(s) # standard deviation of 2, sample size of 10
  SElogR0_E = sd(logR0_E, na.rm=T)/sqrt(s) # standard deviation of 2, sample size of 10
  SElogR0_F = sd(logR0_F, na.rm=T)/sqrt(s) # standard deviation of 2, sample size of 10
  SElogR0_G = sd(logR0_G, na.rm=T)/sqrt(s) # standard deviation of 2, sample size of 10
  SElogR0_H = sd(logR0_H, na.rm=T)/sqrt(s) # standard deviation of 2, sample size of 10
  SElogR0_I = sd(logR0_I, na.rm=T)/sqrt(s) # standard deviation of 2, sample size of 10
  
  lower_A = exp(mean(logR0_A, na.rm=T) - Z*SElogR0_A) 
  upper_A = exp(mean(logR0_A, na.rm=T) + Z*SElogR0_A) 
  lower_B = exp(mean(logR0_B, na.rm=T) - Z*SElogR0_B) 
  upper_B = exp(mean(logR0_B, na.rm=T) + Z*SElogR0_B) 
  lower_C = exp(mean(logR0_C, na.rm=T) - Z*SElogR0_C) 
  upper_C = exp(mean(logR0_C, na.rm=T) + Z*SElogR0_C) 
  lower_D = exp(mean(logR0_D, na.rm=T) - Z*SElogR0_D) 
  upper_D = exp(mean(logR0_D, na.rm=T) + Z*SElogR0_D) 
  lower_E = exp(mean(logR0_E, na.rm=T) - Z*SElogR0_E) 
  upper_E = exp(mean(logR0_E, na.rm=T) + Z*SElogR0_E) 
  lower_F = exp(mean(logR0_F, na.rm=T) - Z*SElogR0_F) 
  upper_F = exp(mean(logR0_F, na.rm=T) + Z*SElogR0_F) 
  lower_G = exp(mean(logR0_G, na.rm=T) - Z*SElogR0_G) 
  upper_G = exp(mean(logR0_G, na.rm=T) + Z*SElogR0_G) 
  lower_H = exp(mean(logR0_H, na.rm=T) - Z*SElogR0_H) 
  upper_H = exp(mean(logR0_H, na.rm=T) + Z*SElogR0_H) 
  lower_I = exp(mean(logR0_I, na.rm=T) - Z*SElogR0_I) 
  upper_I = exp(mean(logR0_I, na.rm=T) + Z*SElogR0_I) 
  
  #R0_A <- read.table(paste0("R0_", sim_index, ".tab"), header=T)
  #cluster_dynamics$mean_R0[1] = as.numeric(R0_A[1]) # Can't use exact R0 because standard deviation includes the high R0 of E and 0.
  cluster_dynamics$mean_R0[1] = mean(R0_A)
  cluster_dynamics$upper_R0[1] = upper_A
  cluster_dynamics$lower_R0[1] = lower_A
  
  cluster_dynamics$mean_R0[2] = mean(R0_B)
  cluster_dynamics$upper_R0[2] = upper_B
  cluster_dynamics$lower_R0[2] = lower_B
  
  cluster_dynamics$mean_R0[3] = mean(R0_C)
  cluster_dynamics$upper_R0[3] = upper_C
  cluster_dynamics$lower_R0[3] = lower_C
  
  cluster_dynamics$mean_R0[4] = mean(R0_D)
  cluster_dynamics$upper_R0[4] = upper_D
  cluster_dynamics$lower_R0[4] = lower_D

  cluster_dynamics$mean_R0[5] = mean(R0_E)
  cluster_dynamics$upper_R0[5] = upper_E
  cluster_dynamics$lower_R0[5] = lower_E
  
  cluster_dynamics$mean_R0[6] = mean(R0_F)
  cluster_dynamics$upper_R0[6] = upper_F
  cluster_dynamics$lower_R0[6] = lower_F
  
  cluster_dynamics$mean_R0[7] = mean(R0_G)
  cluster_dynamics$upper_R0[7] = upper_G
  cluster_dynamics$lower_R0[7] = lower_G
  
  cluster_dynamics$mean_R0[8] = mean(R0_H)
  cluster_dynamics$upper_R0[8] = upper_H
  cluster_dynamics$lower_R0[8] = lower_H
  
  cluster_dynamics$mean_R0[9] = mean(R0_I)
  cluster_dynamics$upper_R0[9] = upper_I
  cluster_dynamics$lower_R0[9] = lower_I
  
  
  return(cluster_dynamics)
  
}# End function
define.clades <- function(sub_tree) {
  ## Grab tree scaled in substitutions/site from time_tree ("intree") so that we know that node numbers line up later on when obtaining temporal information
  sub_tree <- sub_tree
  ## Grab only subtrees that are supported by bootstrap values >90
  ## We may need to change this in case people have other support values
  ## Note that subtrees() function in ape renumbers nodes, so useless here, since at the end we wish to recombine the information
  family_tree <- tidytree::as_tibble(sub_tree)
  ## Need to relabel columns so that "parent" and "node" are "from" and "to" for familyR::get_children function
  colnames(family_tree)[1:2] <- c("from", "to")
  ## The dataframe needs to be transformed back into a data.table for future analyses
  family_tree <- data.table::as.data.table(family_tree)
  assign("family_tree", family_tree, envir = globalenv())
  
  supported_nodes <- family_tree# All of family_tree for simulations
  
  ## Creating a list for all sub_trees using the familyR package (get_children function)
  clades <- list()
  for (node in unique(supported_nodes$to)) {
    clades[[node]] <- familyR::get_children(family_tree, node)
  } # End loop along family_tree_supported
  ## This function introduces several null items in the list, which can be removed by the following:
  clades<-plyr::compact(clades)
  ## Merge information across 'nodes' and 'edges' dataframes within each nested list
  for (i in seq_along(clades)) {
    clades[[i]] <- merge(clades[[i]]$edges,clades[[i]]$nodes, by.x = "to", by.y = "id")
  } # End loop along clades
  ## Restructure list of clades for easy visualization and manipulation and remove
  ## Zeroth level (contains root branch length) from first clade (full tree) if exists
  clades <- lapply(clades, function(x) {
    dplyr::select(x, from, everything()) %>%
      arrange(level) %>%
      filter(!level==0)
  }) %>%
    list.filter(from[1] != rootnode(sub_tree)) # Whole tree is included because support is alwasy 100% for root, so discard# End loop along clades
  
  # Save to global environment or merge next function.
  assign("clades", clades, envir=globalenv())
  
  
} # End defineClades function
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
add.test.clust <- function(true_clusters) {
  dirtrans_clusters <- readRDS(file=paste0("dirtrans_clusters_", opt$sim_index, ".rds"))
  dirtrans_clusters <- lapply(dirtrans_clusters, function(x) {
    x <- x %>%
      mutate(label = paste(label, state, sep="_"))
  })
  n <- length(true_clusters)
  state <- list()
  for (i in seq_along(dirtrans_clusters)){
  state[[i]] <- dirtrans_clusters[[i]]$state[1]
  true_clusters[[n+i]] <- data.frame(cluster_id = paste0(state[[i]], "_direct"), # This will be some arbitrary name you assign to the cluster, such as C1, C2, C3, etc.
                                   taxa = str_replace_all(toString(unique(dirtrans_clusters[[i]]$label)), fixed(" "), ""), # list of taxa 
                                   mean_R0 = cluster_dynamics$mean_R0[cluster_dynamics$state == state[[i]]],
                                   upper_R0 = cluster_dynamics$upper_R0[cluster_dynamics$state == state[[i]]],
                                   lower_R0 = cluster_dynamics$lower_R0[cluster_dynamics$state == state[[i]]])
                                   #dynamic = cluster_dynamics$dynamic[cluster_dynamics$state == state[[i]]],
                                   #birth = cluster_dynamics$birth[cluster_dynamics$state == state[[i]]], stringsAsFactors = F) # dynamic corresponds to the dynamic state (growth, etc.) assigned to the cluster state
  }
  return(true_clusters)
}
merge.nested.clust <- function(clusters) {
  copy <- clusters
  result <- list()
  unwanted <- list()
  for (ct in seq_along(clusters)) {
    for (cc in seq_along(copy)) {
      if (isTRUE(all(clusters[[ct]]$label %in% copy[[cc]]$label) &
                 length(clusters[[ct]]$label) != length(copy[[cc]]$label))) {
        unwanted[[ct]] <- clusters[[ct]]
      } else{NULL}
    } # End loop along copy
  } # End loop along true
  result <- setdiff(clusters, unwanted)
  for (j in seq_along(result)) {
    names(result)[[j]] <- paste0("c", j)
  }
  return(result)
}  
merge.overlap.clust <- function(clusters) {
  copy <- clusters
  unwanted <- list()
  result <- list()
  for (ct in seq_along(clusters)) {
    for (cc in seq_along(copy)) {
      if (isTRUE(sum(copy[[cc]]$label %in% clusters[[ct]]$label) > 0.01*length(copy[[cc]]$label)) &
          isTRUE(names(copy)[[cc]] != names(clusters)[[ct]])) {
        unwanted[[cc]] <- copy[[cc]]
        clusters[[ct]] <- full_join(copy[[cc]], clusters[[ct]], by=c("from", "to", "branch.length", "label"))
      } else{clusters[[ct]] <- clusters[[ct]]}
    } # End loop along copy
  } # End loop along true
  result <- setdiff(clusters, Filter(Negate(function(x) is.null(unlist(x))), unwanted)) %>%
    lapply(., function(x){
      dplyr::select(x, from, to, branch.length, label) %>%
        dplyr::arrange(from,to) 
    }) %>%
    unique()
  
  return(result)
}  # In case you want to remove this and consider only fully nested clusters
branchLengthLimit <- function(tree) {
  p.dist.mat.leaves <- cophenetic(tree)
  get.node.leaf.MPPD <- function(node,tree,distmat){
    nlist <- tips(tree,node)
    foo <- distmat[nlist,nlist]
    return(median(foo[upper.tri(foo,diag=FALSE)]))
  } ## Given a node, tree, and distance matrix, return median   pairwise patristic distance (MPPD) of its leaves
  get.node.full.MPPD <- function(node,tree,distmat){
    nlist <- tips(tree, node)
    elist <- tree$edge[which.edge(tree,nlist),2]
    foo <- distmat[elist,elist]
    return(median(foo[upper.tri(foo,diag=FALSE)]))
  } ## Given a node, tree, and distance matrix, return median pairwise patristic distance (MPPD) of all of its decendants
  pdist.clusttree <- function(tree,distmat, mode=c('leaf', 'all')){
    mode <- match.arg(mode)
    ntips<- Ntip(tree)
    nint <- tree$Nnode # Number of internal nodes
    node_num <- (ntips+2):(ntips+nint)
 
    if(mode=='leaf'){
      distmat <-  p.dist.mat.leaves
      MPPD <- sapply(node_num,get.node.leaf.MPPD,tree,distmat)
      return(data.frame(node_num=node_num, MPPD=MPPD))
    }
    else{
      distmat <-  dist.nodes(tree)
      MPPD <- sapply(node_num,get.node.full.MPPD,tree,distmat)
      return(data.frame(node_num=node_num, MPPD=MPPD))
    }
  } ## Given a tree and (optionally) a distance matrix, return a vector giving the median pairwise patristic distance of the subtree under each internal node
  pdist.clades <- function(clades, tree, distmat, mode=c('leaf', 'all')){
    mode <- match.arg(mode)
   
    if(mode=='leaf'){
      distmat <-  p.dist.mat.leaves
      mclapply(clades, function(x) {
        get.node.leaf.MPPD(x$from[1], tree, distmat)
      }, mc.cores=numCores)
    } else{
      distmat <-  dist.nodes(tree)
      mclapply(clades, function(x) {
        get.node.full.MPPD(x$from[1], tree, distmat)
      }, mc.cores=numCores)
    }
  } ## Determine MPPD for all well-supported clades
  ### Create matrix of each pairwise patristic distance for external leaves using the following
  # leaves <- expand.grid(tree$tip.label,tree$tip.label)
  # p.dist.leaves <- sapply(seq_len(nrow(leaves)), ## Create list of all pairwise combinations of IDs using expand.grid()
  #                         function(k) { #future_sapply actually slower here!
  #                           i <- leaves[k,1]
  #                           j <- leaves[k,2]
  #                           fastDist(tree, i,j)
  #                         })
  # p.dist.mat.leaves <- matrix(p.dist.leaves,
  #                             nrow=Ntip(tree), ncol=Ntip(tree),
  #                             dimnames=list(leaves,leaves))
  
  
  ## Create a vector of MPPDs for plotting and determining branch length limit
  distvec <- pdist.clusttree(tree, mode='leaf')
  hist(distvec$MPPD)
  
  ## Determine MPPDs for all well-supported clades
  clade_MPPD <- pdist.clades(clades, tree, mode='leaf')
  assign("clade_MPPD", clade_MPPD, envir=globalenv())
  
  branch_length_limit <- quantile(distvec$MPPD, opt$threshold)
  return(branch_length_limit)
}
branchWise <- function(tree, branch_length_limit, make_tree) {
#   findThreshold <- function(tree) {
#     sub_tree <- as.phylo(tree)
#     
#     sub_tree <- multi2di(sub_tree)
# #    fit <- tryCatch(skygrowth.map(sub_tree, res=round(180/7)), error=function(e) NULL)
# #    p <- data.frame(time=fit$time, nemed=fit$ne)
# #    fit <- fit_easylinear(p$time, p$nemed, quota=0.99) #Foound this to be the best, except when prob_exit and p_trans are both high
# #    plot(fit)
# #    time_epi_peak <- max(nodeHeights(sub_tree, root.edge=TRUE)) - min(fit@fit$model$x)
#     
#     
#     
#      b0 <- BNPR(multi2di(tree))
#      plot_BNPR( b0 )
#      p <- data.frame(time=rev(b0$summary$time), ne=b0$summary$mean)
#      mL <- try(drm(data=p, ne~time, fct = LL.3(), type = "continuous"), silent=T)
#      p2 <- plot(mL)
#      
#     # dY <- diff(p2$`1`)/diff(p2$time)  # the derivative of your function
#     # dX <- rowMeans(embed(p2$time,2)) # centers the X values for plotting
#     # plot(dX,dY,type="l",main="Derivative") #check
#     # d1 <- data.frame(x=dX, y=dY)
#     # 
#     # time_epi_peak <- max(nodeHeights(tree, root.edge=TRUE)) - d1$x[which.max(d1$y)]
#     
#     fit <- fit_easylinear(p2$`time`, p2$`1`, quota=0.90) #Foound this to be the best, except when prob_exit and p_trans are both high
#     plot(fit)
# #    time_epi_peak <- max(nodeHeights(sub_tree, root.edge=TRUE)) - min(fit@fit$model$x)
#     time_epi_peak <- min(fit@fit$model$x)
#     
#     time.tree.exp <- paleotree::timeSliceTree(sub_tree, time_epi_peak)
#     #  bl_rescaled <- time.tree.exp$edge.length*(time_tree_data$mean.rate)
#     bl_rescaled <- time.tree.exp$edge.length
#     plot(density(bl_rescaled))
#     
#     branch_length_limit <- median(bl_rescaled) ## Ideally will need to multiply branch lengths by estimated rate.
#     #branch_length_limit <- mean(bl_rescaled) ## Ideally will need to multiply branch lengths by estimated rate. # Sometimes median works best
#     ## Maybe we should take both and whichever is larger?
#     #branch_length_limit <- mean(time.tree.exp$edge.length) + qnorm(.95)*(sd(time.tree.exp$edge.length)/sqrt(Ntip(time.tree.exp)))
#     
#     return(branch_length_limit)
#   } # End findThreshold()
  pickClust <- function(clade){
    
    ## Need to add mean_bl column  to original clade list
    clade$mean_bl <- rep(Inf, nrow(clade))
    
    
    ## Set initial values
    current_level <- as.numeric(2)
    current_mean_bl <- as.numeric(-Inf)
    
    
    ## Initiate subclade using first two edges connected to root of clade (level=1)
    sub_clade <- filter(clade, level==1,
                        branch.length <= branch_length_limit)
    
    
    
    while(isTRUE(current_mean_bl <= branch_length_limit)){
      
      # Create a vector of nodes sampled from the subsequent level
      # Each iteration chooses amongst nodes that are connected to the current sub-clade:
      
      next_level_nodes <- filter(clade, level == current_level,
                                 from %in% sub_clade$to)
      
      # A shortlist of possible enlargements of the sub-clade is kept to be able
      # to compare each potential enlargement of the sub-clade and always keep the enlargement
      # if the mean branch length is under the limit
      #
      # The shortlist is enlarged by vertices that are:
      #  1) adjacent to the most recent added node(s)
      #  2) not already IN the sub_clade
      new_node <- dplyr::setdiff(next_level_nodes, sub_clade)
      sub_clade <- rbind(sub_clade,new_node, fill=T)
      
      
      # The branch length is NOT calculated by the branch length of an individual
      # edges leading to nodes in the shortlist BUT on the mean of the nodes in the previous level
      # and added node.
      for (x in 1:nrow(sub_clade)) {
        if (isTRUE(sub_clade$level[x] < current_level &
                   sub_clade$mean_bl[x] != Inf)) {
          sub_clade$mean_bl[x] <- sub_clade$mean_bl[x]
        } else{
          sub_clade$mean_bl[x] <- sum(c(sub_clade$branch.length[x],
                                        subset(sub_clade$branch.length, sub_clade$level < sub_clade$level[x])))/
            length(c(sub_clade$branch.length[x],
                     subset(sub_clade$branch.length, sub_clade$level < sub_clade$level[x])))
        }
      }
      
      # Identify nodes with mean branch length > current mean branch length limit
      # and remove from shortlist
      unwanted_edges_at_current_level <- filter(sub_clade, level==current_level, mean_bl >= branch_length_limit)
      sub_clade <- setdiff(sub_clade, unwanted_edges_at_current_level)
      
      ## Redefine current level and mean branch length, taking into account whether
      ## or not all nodes belonging to the current level have been filtered out
      if (max(sub_clade$level)==current_level) {
        #current_mean_bl <- min(sub_clade$mean_bl[level=current_level])
        current_mean_bl <- mean(sub_clade$branch.length)
        current_level <- current_level+1
      } else {
        current_mean_bl = Inf
        current_level <- current_level+1
      } # End ifelse statment
    } # End tree traversal (while loop)
    if (nrow(sub_clade) == 0) {
      return(NULL)
    }else {
      return(sub_clade)
    }
  } # End pickClust function
  bifurcate <- function(tree, clusters) {
    results <- list()# Copy clusters list for ease
    for (c in seq_along(clusters)) {
      results[[c]] <- extract.clade(tree, clusters[[c]]$from[1])
      results[[c]] <- drop.tip(results[[c]], subset(results[[c]]$tip.label, results[[c]]$tip.label %notin% clusters[[c]]$label),
                               trim.internal = FALSE,
                               collapse.singles = TRUE)
    } #End loop along clusters
    # for (j in seq_along(results)) {
    #   names(results)[[j]] <- paste0("c", j)
    # }
    return(results)
  } # End function
  pullClade <- function(tree, clusters) {
    results <- list()# Copy clusters list for ease
    for (c in seq_along(clusters)) {
      results[[c]] <- extract.clade(sub_tree, clusters[[c]]$from[1])
    } #End loop along clusters
    # for (j in seq_along(results)) {
    #   names(results)[[j]] <- paste0("c", j)
    # }
    return(results)
  } # End function
  addLeaves <- function(tree, clusters) {
    sub_tree <- as_tibble(tree)
    x <- clusters  # Copy clusters list for ease
    for (c in seq_along(x)) {
      x[[c]] <- x[[c]] %>% dplyr::rename(parent=from, node=to) # comparison of clusters and subtree is based on same column names
      for (node in 2:nrow(x[[c]])) {
        if (length(x[[c]]$parent[x[[c]]$parent==x[[c]]$parent[node]]) == 1) { # For each row in a cluster, find parent nodes
          ## that only have one edge (need two for downstream analysis)
          p <- x[[c]]$parent[node]
          n <- x[[c]]$node[node]
          l <- as.data.frame(filter(sub_tree, parent == p & node != n)) ## Find the other edge in sub_tree
          l$branch.length =
            x[[c]]$branch.length[x[[c]]$parent == l$parent &
                                   x[[c]]$node != l$node] ## Replace the current branch.length with the branch.length of the
          ## other edge in the current list
          x[[c]] <- rbind(x[[c]], l) ## Add the new edge to the current cluster
        } # End if statement ## If no lonely edges, let the cluster list be itself
      } # End loop along row
    } # End loop along clusters
    return(x)
  } # End function; this added articial leaves to create bifurcating tree, but changed to the bifurcate f(x) above, which just involes dropping tips.
  
#  branch_length_limit <- findThreshold(sub_tree)
  
  clusters <- mclapply(clades, pickClust, mc.cores=numCores) %>%
    compact() %>%
    merge.nested.clust() %>%
    merge.overlap.clust() %>%
    merge.nested.clust()
  ## Remove singleton nodes (non-bifurcating branches) by using bifurcate() function or extracting entire clade.
  
  if (make_tree=="bifurcate") {
  clusters <- list.filter(clusters, length(label) >= 4)
  clusters <- bifurcate(tree, clusters)
  clusters <- mclapply(clusters, as_tibble, mc.cores=numCores)
  } else {
    if (make_tree=="pullClade") {
      clusters <- pullClade(tree, clusters)
    } else {
      if (make_tree=="addLeaves") {
        clusters <- addLeaves(tree, clusters)
        clusters <- mclapply(clusters, as_tibble, mc.cores=numCores)
        }
      }
    }

  clusters <-  list.filter(clusters, sum(label %in% tree$tip.label) >= 5)
  
  return(clusters)
}
phylopart <- function(tree, branch_length_limit) {
 
  
  clusters <- list()
  for (clade in seq_along(clades)) {
    if (isTRUE(clade_MPPD[[clade]] <= branch_length_limit)) {
      clusters[[clade]] <- clades[[clade]]
    } else{NULL}
  }
  clusters <- compact(clusters) %>%
    merge.nested.clust() %>%
    #  merge.overlap.clust() %>% # Will not work well with support in labels
    #  merge.nested.clust() %>%
    list.filter(length(label) >= 5)
  return(clusters)
}
benchmark <- function(clusters, true_clusters) {
# First benchmark true clusters
compareClusters <- function(x,y) {
   if (isTRUE(
     names(which.max(table(gsub(".+\\_([A-Z])", "\\1", unlist(str_split(x$taxa, ',')))))) == 
     names(which.max(table(gsub(".+\\_([A-Z])", "\\1", y$label)))) &
     sum(y$label %in% unlist(str_split(x$taxa, ','))) >= 0.40*length(unlist(str_split(x$taxa, ','))))) {
     #length(y$label))) {
     state <- names(which.max(table(gsub(".+\\_([A-Z])", "\\1", unlist(str_split(x$taxa, ','))))))
     proportion <- sum(y$label %in% unlist(str_split(x$taxa, ',')))/
       length(unlist(str_split(x$taxa, ','))) 
     }
   else{ 
     state <- names(which.max(table(gsub(".+\\_([A-Z])", "\\1", unlist(str_split(x$taxa, ','))))))
     proportion=0

   }
  taxa <- x$taxa
   return(data.frame(state=state, proportion=proportion, taxa=taxa, stringsAsFactors = F))
 }
 
 performance <- mclapply(true_clusters, function(x) lapply(clusters, function(y) compareClusters(x,y)), mc.cores = numCores)
 performance <- mclapply(performance, rbindlist, mc.cores = numCores)
 ## If multiple clusters as part of the same true cluster are identified, add proportions together (because distinct)
 performance <- lapply(performance, function(x) {
   state=x$state[1]
   proportion=sum(x$proportion)
   taxa=x$taxa
   return(data.frame(state=state, proportion, taxa=taxa, stringsAsFactors = F))
 })
 performance <- rbindlist(performance) %>% dplyr::distinct()
 return(performance)
}
calculateNe <- function(cluster, tree) {
  if ("phylo" %in% class(cluster)) {  
    fit <- ltt(tree, plot=F)
    } else {
      if ("tbl" %in% class(cluster))
      x <- extract.clade(tree, cluster$label)
     #  res=round(max(nodeHeights(x))/days) # Total Ne set to every two weeks, but setting clusters up for every week
      #  fit <- tryCatch(skygrowth.map(x, res=res, tau0=0.1), error=function(e) NULL)
      #  p <- data.frame(time=fit$time, nemed=fit$ne, ne_ci=fit$ne_ci)
      fit <- tryCatch(ltt(x, plot=F), error=function(e) NULL)
    }
  return(fit)
}
modelNe <- function() {
  require(adephylo)
  
  clusterMetadata <- function(tree, clusters) {
    sts <- distRoot(tree@phylo, tips = "all", method="patristic")
    sts <- data.frame(time = sts, ID=names(sts))
    
    
    phylo_list <- lapply(clusters, function(x) {
      extract.clade(tree@phylo, phytools::findMRCA(tree@phylo, x$label[x$label %in% tree@phylo$tip.label]))
    })#, mc.cores=numCores) # End creation of phylo_list
    names(phylo_list) <- names(clusters)
    #  saveRDS(phylo_list, "clusters_as_trees.rds")
    assign("phylo_list", phylo_list, envir = globalenv())
    
    tbl_list <- mclapply(phylo_list, tidytree::as_tibble, mc.cores=numCores)
    tbl_list <- mclapply(tbl_list, function(x) {
      mrsd <- max(sts$time[sts$ID %in% x$label], na.rm=T)
      x <- cbind(x, mrsd=mrsd,
                 cluster_size = length(as.phylo(x)$tip.label),
                 root_age = mrsd-max(node.depth.edgelength(as.phylo(x))) 
      )}, mc.cores=numCores) # End creation of tbl_list
    
    
    
    for (i in seq_along(tbl_list)) {
      tbl_list[[i]]$cluster_id <- names(tbl_list)[[i]]
    }
    
    
    clust_metadata <- rbindlist(tbl_list, fill=T) %>%
      dplyr::select(-parent, -node, -branch.length) # No longer need parenta and node information, since not unique
    names(clust_metadata)[1] <- "ID"
    
    
    
    metadata <- merge(sts, clust_metadata, by="ID", all.x=T)
    
    return(metadata)
  }
  metadata <- clusterMetadata(tree, clusters)
  assign("metadata", metadata, envir = globalenv())
  
  
  total_Ne <- calculateNe(tree@phylo, tree@phylo) # Set to every two weeks
  total_gamma <- total_Ne$gamma
  assign("total_gamma", total_gamma, envir = globalenv())
  
  total_Ne_df <- data.frame(ne=total_Ne$ltt, time=total_Ne$times)
  
  clusters_Ne <- lapply(clusters, function(x) calculateNe(x, tree@phylo))
  clusters_gamma <- mclapply(clusters_Ne, function(x) return(x$gamma))
  assign("clusters_gamma", clusters_gamma, envir = globalenv())
  
  clusters_Ne_df <- mclapply(clusters_Ne, function(x) data.frame(ne=x$ltt, time=x$times), mc.cores=numCores)
  
  c_root_age <- NA
  for (i in seq_along(clusters_Ne_df)) {
    c_root_age[i] <- unique(metadata$root_age[metadata$cluster_id == names(clusters_Ne_df)[i] & !is.na(metadata$cluster_id)])
    #  x1[i] <- max(metadata$time[metadata$cluster_id == names(clusters_Ne_df)[i]], na.rm=T)
    clusters_Ne_df[[i]]$time <- clusters_Ne_df[[i]]$time + c_root_age[i]
  }
  
  reslm_total <- function(Ne_df) {
    result <- lm(Ne_df$ne ~ sin(2*pi/365*Ne_df$time)+cos(2*pi/365*Ne_df$time) +
                   sin(2*pi/365*2*Ne_df$time)+cos(2*pi/365*2*Ne_df$time) +
                   sin(2*pi/365*3*Ne_df$time)+cos(2*pi/365*3*Ne_df$time))
    return(result)
  }
  
  reslm_clusters <- function(Ne_df) {
    result <- lm(Ne_df$ne ~ sin(2*pi/365*Ne_df$time)+cos(2*pi/365*Ne_df$time) +
                   sin(2*pi/365*2*Ne_df$time)+cos(2*pi/365*2*Ne_df$time))
    return(result)
  }
  
  reslm_total <- reslm_total(total_Ne_df)
  reslm_cluster <- mclapply(clusters_Ne_df, reslm_clusters, mc.cores=numCores)
  
  # plot(total_Ne_df$ne~total_Ne_df$time)
  # lines(total_Ne_df$time,reslm_total$fitted,col=2)
  
  # plot(clusters_Ne_df[[2]]$ne~clusters_Ne_df[[2]]$time)
  # lines(clusters_Ne_df[[2]]$time,reslm_cluster[[2]]$fitted,col=2)
  
  
  model.coef <- NULL
  model.coef.total <- rbind.data.frame(model.coef,
                                       data.frame(Beta=reslm_total$coefficients[1], 
                                                  k1=reslm_total$coefficients[2], 
                                                  k2=reslm_total$coefficients[3],
                                                  k3=reslm_total$coefficients[4],
                                                  k4=reslm_total$coefficients[5],
                                                  k5=reslm_total$coefficients[6],
                                                  k6=reslm_total$coefficients[7]))
  model.coef.cluster <- mclapply(reslm_cluster, function(x) {
    rbind.data.frame(model.coef,
                     data.frame(Beta=x$coefficients[1], 
                                k1=x$coefficients[2], 
                                k2=x$coefficients[3],
                                k3=x$coefficients[4],
                                k4=x$coefficients[5]))
  }, mc.cores=numCores)
  
  fun_fs1 = function(x, k1, k2, k3, k4, k5, k6, Beta){
    y = k1*sin(2*pi/365*x)+k2*cos(2*pi/365*x) +
      k3*sin(2*pi/365*2*x)+k4*cos(2*pi/365*2*x) +
      k5*sin(2*pi/365*3*x)+k6*cos(2*pi/365*3*x)  + Beta
    for (i in seq_along(y)) {
      if (isTRUE(y[i]<0)) {
        y[i] = 0
      } else {y[i] = y[i]}}
    return(y)}
  
  fun_fs2 = function(x, k1, k2, k3, k4,  Beta){
    y = k1*sin(2*pi/365*x)+k2*cos(2*pi/365*x) +
      k3*sin(2*pi/365*2*x)+k4*cos(2*pi/365*2*x)  + Beta
    for (i in seq_along(y)) {
      if (isTRUE(y[i]<0)) {
        y[i] = 0
      } else {y[i] = y[i]}}
    return(y)}
  
  
  fun_fs_total <- function(x) fun_fs1(x, model.coef.total$k1, model.coef.total$k2,
                                      model.coef.total$k3, model.coef.total$k4,
                                      model.coef.total$k5, model.coef.total$k6, model.coef.total$Beta)
  
  total_Ne_model <- fun_fs_total(seq(1,max(nodeHeights(tree@phylo)), 1))
  total_Ne_model <- data.frame(time=seq(1,max(nodeHeights(tree@phylo)), 1), Ne=total_Ne_model)
  assign("total_Ne_model", total_Ne_model, envir = globalenv())
  
  dates <- list()
  clusters_Ne_model <- list()
  fun_fs_cluster <- list()
  for (i in seq_along(clusters_Ne_df)) {
    fun_fs_cluster[[i]] <- function(x) fun_fs2(x, model.coef.cluster[[i]]$k1, model.coef.cluster[[i]]$k2,
                                               model.coef.cluster[[i]]$k3, model.coef.cluster[[i]]$k4,model.coef.cluster[[i]]$Beta)
    dates[[i]] <- seq(min(clusters_Ne_df[[i]]$time), max(clusters_Ne_df[[i]]$time), 1)
    clusters_Ne_model[[i]] <- fun_fs_cluster[[i]](dates[[i]])
    clusters_Ne_model[[i]] <- data.frame(time=dates[[i]], Ne=clusters_Ne_model[[i]])
  }
  assign("clusters_Ne_model", clusters_Ne_model, envir = globalenv())
  
  
  clusters_Re <- mclapply(clusters_Ne_model, function(x) {
    calculateRe(x, conf.level = 0.95)}, mc.cores = numCores)
  assign("clusters_Re", clusters_Re, envir = globalenv())
} # End modelNe function
calculateRe <- function(Ne, conf.level=0.95) {
  s <- 100 # sample size of 100
#  psi <- rnorm(s, 0.038, 0.014) #Duration of infection around 14 days 
  psi <- rnorm(s, 14, 5) #Duration of infection around 14 days -incubation period of 5 days
#  psi <- rnorm(s, 9, 5) #Duration of infection around 14 days -incubation period of 5 days
  Z=qnorm(0.5*(1 + conf.level))
  Re <- list()
  if (isTRUE(nrow(Ne) >2)) {
    for (i in 2:nrow(Ne)) {
      time = Ne$time[i-1]
      Re.dist = sample(1+psi*(Ne$Ne[i]-Ne$Ne[i-1])/((Ne$time[i]-Ne$time[i-1])*Ne$Ne[i-1]),
                       size=s, replace=F)
      logRe = log(Re.dist)
      SElogRe = sd(logRe)/sqrt(s) # standard deviation of 2, sample size of 10
      LCL = exp(mean(logRe) - Z*SElogRe) 
      UCL = exp(mean(logRe) + Z*SElogRe) 
      Re[[i-1]] <- data.frame(time = time, mean_Re=mean(Re.dist), conf.int=paste0("(", LCL, "," ,UCL, ")"))
    }
  } else {Re[[i-1]] <- NULL}
  Re <- do.call("rbind",Re)
}
growthMetrics <- function(clusters, clusters_Ne_model) {
  
  Ne_peak_t <- total_Ne_model$time[total_Ne_model$Ne==max(total_Ne_model$Ne)][1]
  max_t <- max(total_Ne_model$time)
  p_time_growth <- Ne_peak_t/max_t
  Ne_growth <- total_Ne_model[total_Ne_model$time <= Ne_peak_t,]
  Ne_lm <- lm(Ne~time, data=Ne_growth)
  total_growth_rate <- coef(Ne_lm)[2]
  
  
  
  growth_criteria <- NULL
  ## Calculate slopes for Ne and ltt, and record R0
  for (i in seq_along(clusters)) {
    cluster_state <- names(which.max(table(gsub(".+([A-Z])$", "\\1", clusters[[i]]$label))))
    cluster_taxa = paste0(clusters[[i]]$label[clusters[[i]]$label!="NA"], collapse=",")
    cluster_size = length(clusters[[i]]$label)
    ## Find amount of time spent in growth phase (based on maximum)
    Ne_peak_t <- clusters_Ne_model[[i]]$time[clusters_Ne_model[[i]]$Ne==max(clusters_Ne_model[[i]]$Ne)][1]
    max_t <- max(clusters_Ne_model[[i]]$time)
    p_time_growth <- Ne_peak_t/max_t
    Ne_growth <- clusters_Ne_model[[i]][clusters_Ne_model[[i]]$time <= Ne_peak_t,]
    ## Find growth rate
    Ne_lm <- lm(Ne~time, data=Ne_growth)
    growth_rate <- coef(Ne_lm)[2]
    rel_growth_rate <- growth_rate/total_growth_rate
    rife_metric <- p_time_growth*rel_growth_rate
    ## Report whether growth occurred 
    mean_R0 <- clusters_Re[[i]]$mean_Re[1]
    lower_R0 <- as.numeric(gsub("\\((.+)\\,.+", "\\1", as.character(clusters_Re[[i]]$conf.int[1])[1]))
    upper_R0 <- as.numeric(gsub(".+\\,(.+)\\)", "\\1", as.character(clusters_Re[[i]]$conf.int[1])[1]))
    gamma <- clusters_gamma[[i]]
    growth_criteria <- rbind(growth_criteria, 
                             data.frame(state = cluster_state,
                                        taxa = cluster_taxa,
                                        size = cluster_size,
                                        p_time_growth = p_time_growth,
                                        growth_rate = growth_rate,
                                        rel_growth_rate = rel_growth_rate,
                                        pybus_gamma = gamma,
                                        mean_R0 = mean_R0,
                                        lower_R0 = lower_R0,
                                        upper_R0 = upper_R0, 
                                        rife_metric = rife_metric,
                                        dynamic=NA,
                                        priority=NA,
                                        stringsAsFactors = F))
    } # End loop along clusters

  # Label as growing, not growing, or dead according to ltt slope
  for (i in 1:nrow(growth_criteria)) {
    if (isTRUE(growth_criteria$p_time_growth[i] == 1)){
      growth_criteria$dynamic[i] <- "growth"
    } else {
    if (isTRUE(growth_criteria$p_time_growth[i] >= 0.5 & growth_criteria$p_time_growth[i] < 1)){
      growth_criteria$dynamic[i] <- "static"
    } else {
      if (isTRUE(growth_criteria$p_time_growth[i] < 0.5)){
        growth_criteria$dynamic[i] <- "decay"
      } 
    }
    }
    if (isTRUE(growth_criteria$rife_metric[i] >= 0.5)){
      growth_criteria$priority[i] <- "high priority"
    } else {
      if (isTRUE(growth_criteria$rife_metric[i] < 0.5)){
        growth_criteria$priority[i] <- "low priority"
    }
    }
  }
  return(growth_criteria)
}
getContempTaxa <- function(sub_tree) {
  tip_heights <- NULL
  sub_tree <- as.phylo(sub_tree)
  for (i in sub_tree$tip.label) {
    tip_heights[i] <- fastHeight(sub_tree, i, i)}
  tip_heights <- plyr::round_any(tip_heights, accuracy = 0.01)
  contemp_taxa <- names(tip_heights[tip_heights <= max(tip_heights) & tip_heights >= (max(tip_heights)-0.02)])[-1] #Sampled individuals within a week
  #of the max time are still considered
  return(contemp_taxa)
}
connectClust <- function(sub_tree, clusters) {
  
  cluster_data <- clusters
  dup_cluster_data <- cluster_data
  full_tree <- as_tibble(sub_tree)
  
  for (i in seq_along(cluster_data)) {
    cluster_data[[i]]$birth_origin <- NA
    for (j in seq_along(dup_cluster_data[-i])) {
      ## If the parent of origin node of one cluster (found by referencing full tree) is a child node in one of the other clusters... 
      if (isTRUE(tidytree::parent(full_tree, dup_cluster_data[[j]]$parent[1])$parent %in% cluster_data[[i]]$node |
          ## Or if the parent of the parent of the origin node of one cluster (found by referencing full tree) is a child node in one of the other clusters... (meaning clusters are allowed to be separated by up to two branches)
          tidytree::parent(full_tree, 
                           tidytree::parent(full_tree, dup_cluster_data[[j]]$parent[1])$parent)$parent %in% cluster_data[[i]]$node)) {
        cluster_data[[j]]$birth_origin == paste0("cluster_", cluster_data[[i]]$parent[1] )
      } #End if statement
    } # End for loop along duplicate list
  } # End for loop along original list
  return(cluster_data)
} #End connectClust function

###########################################################################################################################


tree_list = list.files(pattern=paste0("sim_", opt$sim_index, "_.+\\.tree$"))
print("reading in tree...")
tree = lapply(tree_list, read.beast)[[1]]
sub_tree <- tree@phylo



## Set known dynamics ##################################################################################
cluster_dynamics <- true.cluster.dyn(conf.level=0.95)

####### Define clades and find true clusters ###########################################################
define.clades(sub_tree)

print("Finding true clusters...")
true_clusters <- find.true.clusters(family_tree, clades) %>%
  merge.overlap.clust.true() %>%
  list.filter(length(unlist(str_split(taxa, ','))) >= 6)
## Add test "A" clusters, which is a clade of sequences taken directly from the full transmission tree
true_clusters <- add.test.clust(true_clusters)
states_present <- lapply(true_clusters, function(x) {
  (gsub(".+_([A-Z]).+", "\\1", x$taxa))
})
states_present <- data.frame(sim=opt$sim_index, state = do.call("rbind", states_present), stringsAsFactors = F) %>%
  group_by(sim, state) %>%
  dplyr::summarise(state_count = n())

#######Pick clusters ##################################################################################
write("Picking clusters based on user choice of algorithm....")
branch_length_limit <- branchLengthLimit(sub_tree)
if (opt$cluster == "b") {
  clusters <- branchWise(sub_tree, branch_length_limit, make_tree=opt$leaves)
} else {
  if (opt$cluster == "c") {
    clusters <- phylopart(sub_tree, branch_length_limit)
    clusters <- mclapply(clusters, function(x) dplyr::rename(x, parent=from, node=to), mc.cores=numCores)
  } else {
    write("Incorrect cluster_picking algorithm choice. Please choose between 'b' (branch-wise) or 'c' (clade-wise) and run script again.")
  }
}


### Benchmarking ###################################################################################

print("Benchmarking performance")
performance <- benchmark(clusters, true_clusters)


## Calculate growth metrics ##########################################################################
print("Calculating Ne(t) and Re(t)...")
modelNe()

### Relative Ne ########################################################################################

relative_Ne <- mclapply(clusters_Ne_model, function(x) {
  Ne_merged <- merge(total_Ne_model, x, by="time", all.x=T)
  Ne_merged$rel_Ne <- Ne_merged$Ne.y/Ne_merged$Ne.x
  return(dplyr::select(Ne_merged, time, rel_Ne))
}, mc.cores=numCores)

## Populate growth metric table ########################################################################

growth_criteria <- growthMetrics(clusters, clusters_Ne_model)


## Classify births ######################################################################################

clusters <- connectClust(sub_tree, clusters) # This needs to be the original clusters, rather than the
## filled clusters (same for below)

for (i in seq_along(clusters)) {
  growth_criteria$birth <- NA
  if(!is.na(clusters[[i]]$birth_origin[1])) {
    growth_criteria$birth[i] <- "birth"
  } else {NULL}
}

# ## Classify deaths ###################################################################################

## Grab contemporaneous taxa to classify deaths
contemp_taxa <- getContempTaxa(sub_tree)

for (i in seq_along(clusters)) {
  if (isTRUE(any(clusters[[i]]$label %in% contemp_taxa))) {
    growth_criteria$death[i] <- NA
  } else{growth_criteria$death[i] <- "death"} # End if-else statement
} # End for loop

## Now populate performance table based on classification ###############################################

for (i in 1:nrow(growth_criteria)) {
  for (j in 1:nrow(cluster_dynamics)) {
    
    # R0_diff
    if (isTRUE(growth_criteria$state[i] == cluster_dynamics$state[j])) {
    growth_criteria$R0_diff[i] <- abs(growth_criteria$mean_R0[i] - cluster_dynamics$mean_R0[j])
    } else {NULL}
    
    # R0 intervals           
    if (isTRUE(growth_criteria$state[i] == cluster_dynamics$state[j] &
               (growth_criteria$lower_R0[i] <= cluster_dynamics$lower_R0[j] &
               growth_criteria$upper_R0[i] >= cluster_dynamics$upper_R0[j]) |
               (growth_criteria$lower_R0[i] >= cluster_dynamics$lower_R0[j] &
                growth_criteria$upper_R0[i] <= cluster_dynamics$upper_R0[j]))){
    growth_criteria$R0_overlap[i] <- "included"
    } else {
      if (isTRUE(growth_criteria$state[i] == cluster_dynamics$state[j] &
               (growth_criteria$lower_R0[i] <= cluster_dynamics$lower_R0[j] & 
               growth_criteria$upper_R0[i] >= cluster_dynamics$lower_R0[j]) |
               (growth_criteria$upper_R0[i] >= cluster_dynamics$upper_R0[j] &
               growth_criteria$lower_R0[i] <= cluster_dynamics$upper_R0[j]))){# If true cluster dynamic classification is contained within the reported cluster dynamic,
      growth_criteria$R0_overlap[i] <- "overlapping"
      } else {growth_criteria$R0_overlap[i] <- "non-overlapping"}
    }
  }
}

## Merge with performance metrics ###############################################
`%!=na%` <- function(e1, e2) (is.na(e1) & is.na(e2))

final_df <- NULL
for (i in 1:nrow(growth_criteria)) {
  for (j in 1:nrow(performance)) {
    for (k in 1:nrow(cluster_dynamics)) {
    if (isTRUE(growth_criteria$state[i] == performance$state[j] &
        any(unlist(str_split(performance$taxa[j], ',')) %in% 
            unlist(str_split(growth_criteria$taxa[i], ','))) &
        performance$state[j] == cluster_dynamics$state[k])) {
      state <- performance$state[j]
      proportion <- performance$proportion[j]
      taxa <- performance$taxa[j]
      dynamic <- if_else(growth_criteria$dynamic[i] == cluster_dynamics$dynamic[k] |
                           (growth_criteria$dynamic[i] == "static" &
                           cluster_dynamics$dynamic[k] == "background"), "TRUE", "FALSE")
      birth <- if_else(growth_criteria$birth[i] == cluster_dynamics$birth[k] |
                         growth_criteria$birth[i] %!=na% cluster_dynamics$birth[k], "TRUE", "FALSE")
      priority <- if_else((cluster_dynamics$dynamic[k] == "background" &
                            growth_criteria$priority[i] == "low priority") |
                            (cluster_dynamics$dynamic[k] != "background" &
                               growth_criteria$priority[i] == "high priority"), "TRUE", "FALSE")
      R0_diff <- growth_criteria$R0_diff[i]
      R0_overlap <- growth_criteria$R0_overlap[i]

    } else {
      if (isTRUE(performance$proportion[j] == 0)) {
      state <- performance$state[j]
      proportion <- performance$proportion[j]
      taxa <- performance$taxa[j]
      dynamic <- NA
      birth <- NA
      priority <- NA
      R0_diff <- NA
      R0_overlap <- NA
      }
    }
       final_df <- rbind(final_df,
                        data.frame(state=state, proportion=proportion, taxa=taxa, dynamic=dynamic,
                             birth=birth, priority=priority, R0_diff=R0_diff, R0_overlap=R0_overlap)) %>%
        distinct()
    }
  }
}


## Save performance metrics for each run ################################################################
if (opt$leaves == "") {
  setwd(paste0("/blue/salemi/brittany.rife/dynamite/simulations/", 
               opt$cluster, "/perc_",
               opt$threshold))
} else {
  setwd(paste0("/blue/salemi/brittany.rife/dynamite/simulations/", 
             opt$cluster, "/", 
             opt$leaves, "/perc_",
             opt$threshold))
}

saveRDS(clusters, file=paste0("clusters_", opt$sim_index, ".rds"))
saveRDS(true_clusters, file=paste0("true_clusters_", opt$sim_index, ".rds"))

growth_criteria$sim <- opt$sim_index
final_df$sim <- opt$sim_index

write.table(states_present, file=paste0('num_true_clusters_', opt$sim_index, ".tab"), sep='\t', quote=F, row.names = F)
write.table(growth_criteria, file=paste0('sim_growth_stats_', opt$sim_index, ".tab"), sep='\t', quote=F, row.names = F) 
write.table(final_df, file=paste0('sim_performance_', opt$sim_index, ".tab"), sep='\t', quote=F, row.names=F) 

#########################################################################################


