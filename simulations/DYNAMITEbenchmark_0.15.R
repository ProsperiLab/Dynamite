
## Initialize stable working environment and store time of intiation
rm(list=ls())

## Apply benchmarking to each simulation
args = commandArgs(trailingOnly=TRUE)
sim_index = as.numeric(args[1]) # arg 1 fraction new
#sim_index <- "test" # test

`%notin%` <- Negate(`%in%`) # Just plain useful


# List of packages for session
.packages <-  c("treeio", "phytools", "remotes", "dplyr", "plyr", "tidyr", "tidytree", "data.table", 
                "parallel", "stringr", "rlist", "paleotree", "phylodyn",
                "ggtree", "drc", "growthrates", "TreeTools", "geiger") # May need to incorporate code for familyR (https://rdrr.io/github/emillykkejensen/familyR/src/R/get_children.R) i fno longer supported.
github_packages <- c("emillykkejensen/familyR", "mrc-ide/skygrowth") 

# # Install CRAN packages (if not already installed)
# .inst <- .packages %in% installed.packages()
# if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# .inst_github <- .packages %in% installed.packages()
# ## Install GitHub packages(if not already installed)
# if(length(github_packages[!.inst_github]) > 0) try(remotes::install_github(github_packages[!.inst_github]))
# if(length(github_packages[!.inst_github]) > 0) try(devtools::install_github(github_packages[!.inst_github]))

# Load packages into session 
lapply(.packages, require, character.only=TRUE)
lapply(gsub(".+\\/(.+)", "\\1", github_packages), require, character.only=TRUE)
numCores <- detectCores()



### Functions ##############################################################################################################


true.cluster.dyn <- function(conf.level=0.95){
  cluster_dynamics <- data.frame(state=rbind('A','B','C','D','E')
                                 #dynamic = rbind('background', 'static', 'growth', 'static', 'static'),
                                 #birth = rbind(NA,NA,NA,NA, 'birth'),
                                 #death = rbind(NA,NA,NA,NA,NA), stringsAsFactors = F
                                 )
  
  ## Estimate R0 for each cluster based on paramaters used in simulation
  s <- 100 # sample size of 100
  n_contacts_A <- rnorm(s, 20, 1)
  n_contacts_BC <- rnorm(s, 4, 1)
  n_contacts_D <- rnorm(s, 8, 1)
  n_contacts_E <- rnorm(s, 4, 1)
  t_incub <- rnorm(s, 5, 2)
  t_exit <- rnorm(s, 14, 3)
  p_trans_A <- 0.015
  p_trans_B <- 0.125
  p_trans_C <- 0.175
  p_trans_D <- 0.0875
  p_trans_E <- 0.275
  
  
  Z=qnorm(0.5*(1 + conf.level))
  
  R0_A = sample(n_contacts_A*p_trans_A*(t_exit-t_incub),
                size=s, replace=F)
  R0_B = sample(n_contacts_BC*p_trans_B*(t_exit-t_incub),
                size=s, replace=F)
  R0_C = sample(n_contacts_BC*p_trans_C*(t_exit-t_incub),
                  size=s, replace=F)
  R0_D = sample(n_contacts_D*p_trans_D*(t_exit-t_incub),
                size=s, replace=F)
  R0_E = sample(n_contacts_E*p_trans_E*(t_exit-t_incub),
                 size=s, replace=F)
  
  logR0_A = log(R0_A)
  logR0_B = log(R0_B)
  logR0_C = log(R0_C)
  logR0_D = log(R0_D)
  logR0_E = log(R0_E)
  
  SElogR0_A = sd(logR0_A, na.rm=T)/sqrt(s) # standard deviation of 2, sample size of 10
  SElogR0_B = sd(logR0_B, na.rm=T)/sqrt(s) # standard deviation of 2, sample size of 10
  SElogR0_C = sd(logR0_C, na.rm=T)/sqrt(s) # standard deviation of 2, sample size of 10
  SElogR0_D = sd(logR0_D, na.rm=T)/sqrt(s) # standard deviation of 2, sample size of 10
  SElogR0_E = sd(logR0_E, na.rm=T)/sqrt(s) # standard deviation of 2, sample size of 10

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
  dirtrans_clusters <- readRDS(file=paste0("dirtrans_clusters_", sim_index, ".rds"))
  n <- length(true_clusters)
  state <- list()
  for (i in seq_along(dirtrans_clusters)){
  state[[i]] <- dirtrans_clusters[[i]]$state[1]
  true_clusters[[n+i]] <- data.frame(cluster_id = paste0(state[[i]], "_direct"), # This will be some arbitrary name you assign to the cluster, such as C1, C2, C3, etc.
                                   taxa = str_replace_all(toString(dirtrans_clusters[[i]]$host), fixed(" "), ""), # list of taxa 
                                   mean_R0 = cluster_dynamics$mean_R0[cluster_dynamics$state == state[[i]]],
                                   upper_R0 = cluster_dynamics$upper_R0[cluster_dynamics$state == state[[i]]],
                                   lower_R0 = cluster_dynamics$lower_R0[cluster_dynamics$state == state[[i]]])
                                   #dynamic = cluster_dynamics$dynamic[cluster_dynamics$state == state[[i]]],
                                   #birth = cluster_dynamics$birth[cluster_dynamics$state == state[[i]]], stringsAsFactors = F) # dynamic corresponds to the dynamic state (growth, etc.) assigned to the cluster state
  }
  return(true_clusters)
}
fastDist<-function(tree,sp1,sp2){
  fastHeight(tree,sp1,sp1)+fastHeight(tree,sp2,sp2)-
    2*fastHeight(tree,sp1,sp2)}
get.node.leaf.MPPD <- function(node,tree,distmat){
  nlist <- tips(tree,node)
  foo <- distmat[nlist,nlist]
  return(median(foo[upper.tri(foo,diag=FALSE)]))
} ## Given a node, tree, and distance matrix, return median pairwise patristic distance (MPPD) of its leaves
get.node.full.MPPD <- function(node,tree,distmat){
  nlist <- tips(tree, node)
  elist <- tree$edge[which.edge(tree,nlist),2]
  foo <- distmat[elist,elist]
  return(median(foo[upper.tri(foo,diag=FALSE)]))
} ## Given a node, tree, and distance matrix, return median pairwise patristic distance (MPPD) of all of its decendants
pdist.clusttree <- function(tree,distmat=NULL,mode=c('leaf','all')){
  mode <- match.arg(mode)
  if(is.null(distmat)){
    if(mode=='leaf'){ distmat <-  p.dist.mat.leaves}
    else{ distmat <-  dist.nodes(tree) }
  }
  ntips<- Ntip(tree)
  nint <- tree$Nnode # Number of internal nodes
  node_num <- (ntips+2):(ntips+nint)
  #  node_sbsmpl <- sample((ntips+1):(ntips+nint), 0.50*tree$Nnode)
  if(mode=='leaf'){
    MPPD <- sapply(node_num,get.node.leaf.MPPD,tree,distmat)
    return(data.frame(node_num=node_num, MPPD=MPPD))
  }
  else{
    #    return(sapply(node_sbsmpl,get.node.full.MPPD,tree,distmat))
    MPPD <- sapply(node_num,get.node.full.MPPD,tree,distmat)
    return(data.frame(node_num=node_num, MPPD=MPPD))
  }
} ## Given a tree and (optionally) a distance matrix, return a vector giving the median pairwise patristic distance of the subtree under each internal node
pdist.clades <- function(clades, tree, distmat=NULL, mode=c('leaf', 'all')){
  mode <- match.arg(mode)
  if(is.null(distmat)){
    if(mode=='leaf'){ distmat <-  p.dist.mat.leaves}
    else{ distmat <-  dist.nodes(tree) }
  }
  if(mode=='leaf'){
    mclapply(clades, function(x) {
      get.node.leaf.MPPD(x$from[1], tree, distmat)
    }, mc.cores=numCores)
  } else{
    mclapply(clades, function(x) {
      get.node.full.MPPD(x$from[1], tree, distmat)
    }, mc.cores=numCores)
  }
} ## Determine MPPD for all well-supported clades
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
      if (isTRUE(sum(copy[[cc]]$label %in% clusters[[ct]]$label) > 0.05*length(copy[[cc]]$label)) &
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
benchmark <- function(clusters, true_clusters) {
  clusters <- lapply(clusters, as_tibble)
  performance <- data.frame(state=rep(NA, length(clusters)), proportion=rep(NA, length(clusters)))
  for (ct in seq_along(true_clusters)) {
    for (cc in seq_along(clusters)) {
      if (sum(unlist(str_split(true_clusters[[ct]]$taxa, ',')) %in% clusters[[cc]]$label)>0.70* #>80% of true cluster sequences in identified cluster
          length(unlist(str_split(true_clusters[[ct]]$taxa, ','))) |
          sum(unlist(str_split(true_clusters[[ct]]$taxa, ',')) %in% clusters[[cc]]$label)>0.70* # Or... #>80% of identified cluster sequences in true cluster
          length(clusters[[cc]]$label) ) {
        ## Change state so that state is the majority state present, rather than root
        performance$state[cc] = names(which.max(table(gsub(".+\\_([A-Z])", "\\1", clusters[[cc]]$label))))
        performance$proportion[cc] <- sum(clusters[[cc]]$label %in% unlist(str_split(true_clusters[[ct]]$taxa, ',')))/
          length(unlist(str_split(true_clusters[[ct]]$taxa, ',')))
      }
      if (is.na(performance$proportion[cc])) {
        performance$state[cc] <-  names(which.max(table(gsub(".+\\_([A-Z])", "\\1", clusters[[cc]]$label))))
        performance$proportion[cc] <- 0
      } # End if-else statement
    } # End loop along copy
  } # End loop along true
  return(performance)
} 

addLeaves <- function(sub_tree, clusters) {
  sub_tree <- as_tibble(sub_tree)
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
calculateNe <- function(cluster) {
#  heights <- sapply(cluster$node, function(x) nodeheight(as.phylo(sub_tree), x))
#  times <- sapply(cluster$node, function(x) nodeheight(as.phylo(time_tree), x))
#  rtt <- summary(lm(heights~times))$adj.r.squared
#  if (rtt >= 0.10) {
    x <- as_tibble(cluster)   
    x <- arrange(x, parent, node)
    class(x) = c("tbl_tree", class(x))
    x <- as.phylo(x)
    x <- rtree(n = Ntip(x), tip.label = x$tip.label, br = x$edge.length) # Have to create new tree because nodes need to be numbered sequentially
 # } else{NULL}
  #b0 <- BNPR(x)
  #plot_BNPR( b0 )
  #p <- data.frame(time=rev(b0$x), Ne=b0$effpopmean)
  fit <- tryCatch(skygrowth.map(x, res=10), error=function(e) NULL)
  gr <- growth.plot(fit)
  p <- data.frame(time=fit$time, nemed=fit$ne)
  return(p)
}
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
      Re.dist = sample(1+psi*(Ne$nemed[i]-Ne$nemed[i-1])/((Ne$time[i]-Ne$time[i-1])*Ne$nemed[i-1]),
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
growthMetrics <- function(clusters) {
  growth_criteria <- NULL
  ## Calculate slopes for Ne and ltt, and record R0
  for (i in seq_along(clusters)) {
    clusters[[i]] <- as_tibble(clusters[[i]])   
    clusters[[i]] <- arrange(clusters[[i]], parent, node)
    class(clusters[[i]]) = c("tbl_tree", class(clusters[[i]]))
    clusters[[i]] <- as.phylo(clusters[[i]])
    cluster_state <- names(which.max(table(gsub(".+([A-Z])$", "\\1", clusters[[i]]$tip.label))))
    cluster_taxa = paste0(clusters[[i]]$tip.label[clusters[[i]]$tip.label!="NA"], collapse=",")
    cluster_size = length(clusters[[i]]$tip.label)
    ## Find slope of Ne after peak growth
    #  Ne <- data.frame(time = clusters_Ne[[i]]$time, nemed = clusters_Ne[[i]]$nemed)
    #  Ne_growth <- fit_easylinear(Ne$time, Ne$nemed, quota=0.90)
    #  Ne_peak_t <- max(Ne_growth@fit$model)
    # dY <- diff(Ne$nemed)/diff(Ne$time)  # the derivative of your function
    # dX <- rowMeans(embed(Ne$time,2)) # centers the X values for plotting
    # d1 <- data.frame(x=dX, y=dY)
    # Ne_peak_t <- d1$x[which.max(d1$y)]
    # Ne_post_peak <- subset(Ne, Ne$time >= Ne_peak_t)
#    Ne_50 <- subset(Ne, Ne$time >= 0.50*Ne$time)
#    Ne_50_slope <- tryCatch(as.numeric(coef(lm(Ne_50$nemed~Ne_50$time))[2]), error=function(e) NA)
#    Ne_pr <- tryCatch(summary(lm(Ne_50$nemed~Ne_50$time))$coefficients[2,4], error=function(e) NA)
    ## Do the same for ltt  
#    ltt_data <- data.frame(time = clusters_ltt[[i]]$times, ltt = clusters_ltt[[i]]$ltt) %>%
#      dplyr::group_by(time) %>%
#      dplyr::summarize(lineages=mean(ltt))
    #     ltt_growth <- fit_easylinear(ltt_data$time, ltt_data$lineages, quota=0.90)
    #     ltt_peak_t <- max(ltt_growth@fit$model)
    # dY <- diff(ltt_data$lineages)/diff(ltt_data$time)  # the derivative of your function
    # dX <- rowMeans(embed(ltt_data$time,2)) # centers the X values for plotting
    # d1 <- data.frame(x=dX, y=dY)
    #    ltt_peak_t <- d1$x[which.max(d1$y)]
    #    ltt_peak_t <- d1$x[which.max(d1$y)]
    #    ltt_post_peak <- subset(ltt_data, ltt_data$time >= ltt_peak_t)
#    ltt_50 <- subset(ltt_data, ltt_data$time >= 0.50*ltt_data$time)
#    ltt_50_slope <- as.numeric(coef(lm(ltt_50$lineages~ltt_50$time))[2])
    # if (isTRUE(is.na(ltt_post_peak_slope) & 
    #     ltt_data$time[ltt_data$lineages==ltt_peak_t] == max(ltt_data$time))) { # If slope is NA because last time point is the peak ltt,
    #   ltt_post_peak_slope <- Inf
    #   }
#    ltt_pr <- summary(lm(ltt_50$lineages~ltt_50$time))$coefficients[2,4]
    mean_R0 <- clusters_Re[[i]]$mean_Re[1]
    lower_R0 <- as.numeric(gsub("\\((.+)\\,.+", "\\1", as.character(clusters_Re[[i]]$conf.int[1])[1]))
    upper_R0 <- as.numeric(gsub(".+\\,(.+)\\)", "\\1", as.character(clusters_Re[[i]]$conf.int[1])[1]))
    growth_criteria <- rbind(growth_criteria, 
                             data.frame(state = cluster_state,
                                        taxa = cluster_taxa,
                                        size = cluster_size,
                                        #Ne_50_slope = Ne_50_slope,
                                        #Ne_pr = Ne_pr,
                                        #ltt_50_slope = ltt_50_slope,
                                        #ltt_pr = ltt_pr,
                                        mean_R0 = mean_R0,
                                        lower_R0 = lower_R0,
                                        upper_R0 = upper_R0, stringsAsFactors = F))
    } # End loop along clusters

  ## Label as growing, not growing, or dead according to ltt slope
  # for (i in 1:nrow(growth_criteria)) {
  #   if (isTRUE(growth_criteria$ltt_pr[i] > 0.10)){
  #     growth_criteria$dynamic[i] <- "static or shrinking"
  #   } else {
  #   if (isTRUE(growth_criteria$ltt_50_slope[i] <= 0)){
  #     growth_criteria$dynamic[i] <- "static or shrinking"
  #   } else {
  #     if (isTRUE(growth_criteria$ltt_50_slope[i] > 0)){
  #       growth_criteria$dynamic[i] <- "growth"
  #       } 
  #     }
  #   } # End if-else statement
  # } # End for loop
  ## Don't really need all of the info in the saved output
  #growth_criteria <- dplyr::select(growth_criteria, cluster, R0, R0.conf.int.hi, R0.conf.int.lo, dynamic)
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
      if (isTRUE(tidytree::parent(full_tree, dup_cluster_data[[j]]$from[1])$parent %in% cluster_data[[i]]$to |
          ## Or if the parent of the parent of the origin node of one cluster (found by referencing full tree) is a child node in one of the other clusters... (meaning clusters are allowed to be separated by up to two branches)
          tidytree::parent(full_tree, 
                           tidytree::parent(full_tree, dup_cluster_data[[j]]$from[1])$parent)$parent %in% cluster_data[[i]]$to)) {
        cluster_data[[j]]$birth_origin == paste0("cluster_", cluster_data[[i]]$from[1] )
      } #End if statement
    } # End for loop along duplicate list
  } # End for loop along original list
  return(cluster_data)
} #End connectClust function

###########################################################################################################################

#setwd("/Users/macbook/Dropbox (UFL)/DYNAMITE/HIVdynamite/nosoi_simulations/different_R0")
setwd("/blue/salemi/brittany.rife/dynamite/simulations/")
tree_list = list.files(pattern=paste0("sim_", sim_index, "_.+\\.tree$"))
print("reading in tree...")
tree = lapply(tree_list, read.beast)[[1]]



## Set known dynamics ##################################################################################
cluster_dynamics <- true.cluster.dyn(conf.level=0.95)

### Pick clusters #######################################################################################
sub_tree <- tree@phylo

####### Compare clusters between DYNAMITE and simulation ##################################
define.clades(sub_tree)

print("Finding true clusters...")
true_clusters <- find.true.clusters(family_tree, clades) %>%
  merge.overlap.clust.true() %>%
  list.filter(length(unlist(str_split(taxa, ','))) >= 9)
## Add test "A" clusters, which is a clade of sequences taken directly from the full transmission tree
true_clusters <- add.test.clust(true_clusters)
states_present <- lapply(true_clusters, function(x) {
  (gsub(".+_([A-Z]).+", "\\1", x$taxa))
})
states_present <- data.frame(sim=sim_index, state = do.call("rbind", states_present), stringsAsFactors = F) %>%
  group_by(sim, state) %>%
  dplyr::summarise(state_count = n())


### DYNAMITE way ######################################################################################
# branch_length_limit <- findThreshold(sub_tree)
# 
# clusters <- mclapply(clades, pickClust, mc.cores=numCores) %>%
#   compact() %>%
#   merge.nested.clust() %>%
#   merge.overlap.clust() %>%
#   merge.nested.clust() %>% ## Why am I having to run this again????
#   list.filter(length(label) >= 5) 
#   #list.filter(sum(label %in% sub_tree$tip.label) >= 5) 


### Phylopart way ######################################################################################

## Leaves method ####################################################################################
## Function to retrieve patristic distances from phytools

### Create matrix of each pairwise patristic distance for external leaves using the following
leaves <- sample(sub_tree$tip.label, 0.50*length(sub_tree$tip.label))
leaves <- expand.grid(leaves,leaves)
p.dist.leaves <- sapply(seq_len(nrow(leaves)), ## Create list of all pairwise combinations of IDs using expand.grid()
                        function(k) { #future_sapply actually slower here!
                          i <- leaves[k,1]
                          j <- leaves[k,2]
                          fastDist(sub_tree, i,j)
                        })
p.dist.mat.leaves <- matrix(p.dist.leaves,
                            nrow=Ntip(sub_tree), ncol=Ntip(sub_tree),
                            dimnames=list(sub_tree$tip.label,sub_tree$tip.label))


## Create a vector of MPPDs for plotting and determining branch length limit
distvec <- pdist.clusttree(sub_tree, mode='all')
hist(distvec$MPPD)

## Determine MPPDs for all well-supported clades
clade_MPPD <- pdist.clades(clades, sub_tree, mode='all')

phylopart.threshold <- 0.15
branch_length_limit <- quantile(distvec$MPPD, phylopart.threshold)

clusters <- list()
for (clade in seq_along(clades)) {
  if (isTRUE(clade_MPPD[[clade]] <= branch_length_limit)) {
    clusters[[clade]] <- clades[[clade]]
  } else{NULL}
}
clusters <- compact(clusters) %>%
  merge.nested.clust() %>%
  merge.overlap.clust() %>%
  merge.nested.clust() %>%
  list.filter(length(label) >= 5)

### Benchmarking ###################################################################################

## Remove singleton nodes (non-bifurcating branches) by using bifurcate() function or extracting entire clade.
#clusters <- bifurcate(sub_tree, clusters)
#clusters <- pullClade(sub_tree, clusters)
clusters <- addLeaves(sub_tree, clusters)
# Pulling the entire clade can result in a change in the majority state present, which needs to be reflected in 
## the performance table

performance <- benchmark(clusters, true_clusters)


## Calculate growth metrics ##########################################################################
print("Calculating Ne(t) and Re(t)...")
clusters_Ne <- mclapply(clusters, function(x) {
  calculateNe(x)
}, mc.cores = numCores)
## When using the bifurcate() function, clusters could be reduced to as few as 3 nodes, for which
## Ne estimation will not work. This results in an empty dataframe, which needs to be filled
## with 'NA' in order to populate growth_criteria dataframe downstream.
print(" Note that when forcing clusters to be bifurcating trees through reduction, clusters could be reduced to as few as 3 nodes, for which Ne estimation will not produce results.")
for (i in seq_along(clusters_Ne)) {
  if (nrow(clusters_Ne[[i]])==0) {
    clusters_Ne[[i]] <- data.frame(time=NA,nemed=NA)
  } # end if-else statement
} # end for loop


clusters_Re <- mclapply(clusters_Ne, function(x) {
  calculateRe(x, conf.level = 0.95)
}, mc.cores = numCores)
## When using the bifurcate() function, clusters could be reduced to as few as 3 nodes, for which
## Ne estimation will not work. This results in an empty dataframe, which needs to be filled
## with 'NA' in order to populate growth_criteria dataframe downstream.
for (i in seq_along(clusters_Re)) {
  if (is.null(clusters_Re[[i]])) {
    clusters_Re[[i]] <- data.frame(time=NA,mean_Re=NA)
  } # end if-else statement
} # end for loop

# clusters_ltt <- mclapply(clusters, function(x) {
#   x <- as_tibble(x)   
#   x <- arrange(x, parent, node)
#   class(x) = c("tbl_tree", class(x))
#   x <- as.phylo(x)
#   x <- rtree(n = Ntip(x), tip.label = x$tip.label, br = x$edge.length) # Have to create new tree because nodes need to be numbered sequentially
#   ltt(x, plot=F)
#   #  png(file=paste0("ltt.cluster_c", i, "_", sim_index, ".png"))
#   #  ltt(multi2di(clusters[[i]]))
#   #  dev.off()
# }, mc.cores=numCores)

## Populate growth metric table ########################################################################

growth_criteria <- growthMetrics(clusters)

## Classify births ######################################################################################

# clusters <- connectClust(sub_tree, clusters) # This needs to be the original clusters, rather than the 
# ## filled clusters (same for below)
# 
# for (i in seq_along(clusters)) {
#   growth_criteria$birth <- NA
#   if(!is.na(clusters[[i]]$birth_origin[1])) {
#     growth_criteria$birth[i] <- "birth"
#   } else {NULL}
# }
# 
# ## Classify deaths ###################################################################################
# 
# ## Grab contemporaneous taxa to classify deaths 
# contemp_taxa <- getContempTaxa(sub_tree)
# 
# for (i in seq_along(clusters)) {
#   if (isTRUE(any(clusters[[i]]$tip.label %in% contemp_taxa))) {
#     growth_criteria$death[i] <- NA
#   } else{growth_criteria$death[i] <- "death"} # End if-else statement
# } # End for loop

## Now populate performance table based on classification ###############################################
performance$R0 <- NA
performance$R0_diff <- NA

for (i in 1:nrow(growth_criteria)) {
  for (j in 1:nrow(cluster_dynamics)) {
    
    # R0_diff
    if (isTRUE(growth_criteria$state[i] == cluster_dynamics$state[j])) {
    performance$R0_diff[i] <- abs(growth_criteria$mean_R0[i] - cluster_dynamics$mean_R0[j])
    } else {NULL}
    
    # R0 intervals           
    if (isTRUE(growth_criteria$state[i] == cluster_dynamics$state[j] &
               cluster_dynamics$mean_R0[j] >= growth_criteria$lower_R0[i] &
               cluster_dynamics$mean_R0[j] <= growth_criteria$upper_R0[i])) {
    performance$R0[i] <- "included"
    } 
    if (isTRUE((is.na(performance$R0[i]) | performance$R0[i] != "included") &
               growth_criteria$state[i] == cluster_dynamics$state[j] &
               ((growth_criteria$lower_R0[i] >= cluster_dynamics$lower_R0[j] & # For any cluster state...
               growth_criteria$lower_R0[i] <= cluster_dynamics$upper_R0[j]) |
               (growth_criteria$upper_R0[i] >= cluster_dynamics$lower_R0[j]&
               growth_criteria$upper_R0[i] <= cluster_dynamics$upper_R0[j])))){# If true cluster dynamic classification is contained within the reported cluster dynamic,
      performance$R0[i] <- "overlapping"
    }
    if (isTRUE(performance$R0[i] %notin% c("included","overlapping"))) {
      performance$R0[i] <- "non-overlapping"}
# Dynamics
  # if (isTRUE(growth_criteria$state[i] == cluster_dynamics$state[j] & # For any cluster state...
  #     grepl(cluster_dynamics$dynamic[j], growth_criteria$dynamic[i]))) {# If true cluster dynamic classification is contained within the reported cluster dynamic,
  # performance$dynamic[i] <- "correct"
  # } else {
  #   if (isTRUE(growth_criteria$state[i] == cluster_dynamics$state[j]) & # For any cluster state...
  #            !isTRUE(grepl(cluster_dynamics$dynamic[j], growth_criteria$dynamic[i]))) {
  #     performance$dynamic[i] <- "incorrect"
  #   }
  # }
# Births
#     if (isTRUE(growth_criteria$state[i] == cluster_dynamics$state[j] & # For any cluster state...
#                cluster_dynamics$birth[j] == "birth" &
#                growth_criteria$birth[i] == "birth")) {# If true cluster dynamic classification is contained within the reported cluster dynamic,
#       performance$birth[i] <- "correct"
#     } else {
#       if (isTRUE(growth_criteria$state[i] == cluster_dynamics$state[j] & # For any cluster state...
#                  cluster_dynamics$birth[j] == "birth" &
#                  growth_criteria$birth[i] == "NA")) {
#         performance$birth[i] <- "incorrect"
#       } else{
#         if (isTRUE(growth_criteria$state[i] == cluster_dynamics$state[j] & # For any cluster state...
#                    cluster_dynamics$birth[j] == "NA" &
#                    growth_criteria$birth[i] == "birth")) {
#           performance$birth[i] <- "incorrect"
#         } else {performance$birth[i] <- "correct"}
#       } # End nested if-else statement
#     } # End if-else statement
# # Deaths
#       if (isTRUE(growth_criteria$state[i] == cluster_dynamics$state[j] & # For any cluster state...
#                  cluster_dynamics$death[j] == "death" &
#                  growth_criteria$death[i] == "death")) {# If true cluster dynamic classification is contained within the reported cluster dynamic,
#         performance$death[i] <- "correct"
#       } else {
#         if (isTRUE(growth_criteria$state[i] == cluster_dynamics$state[j] & # For any cluster state...
#                    cluster_dynamics$death[j] == "death" &
#                    growth_criteria$death[i] == "NA")) {
#           performance$death[i] <- "incorrect"
#         } else{
#           if (isTRUE(growth_criteria$state[i] == cluster_dynamics$state[j] & # For any cluster state...
#                      cluster_dynamics$death[j] == "NA" &
#                      growth_criteria$death[i] == "death")) {
#             performance$death[i] <- "incorrect"
#           } else {performance$death[i] <- "correct"}
#         } # End nested if-else statement
#       } # End if-else statement
  } # End loop along cluster_dynamics
} # End loop along growth_criteria

## Save performance metrics for each run ################################################################
#saveRDS(clusters, file=paste0("clusters_", sim_index, ".rds"))
#saveRDS(true_clusters, file=paste0("true_clusters_", sim_index, ".rds"))

growth_criteria$sim <- sim_index
performance$sim <- sim_index

setwd("/blue/salemi/brittany.rife/dynamite/simulations/perc_0.15")
write.table(states_present, file=paste0('num_true_clusters_', sim_index, ".tab"), sep='\t', quote=F, row.names = F)
write.table(growth_criteria, file=paste0('sim_growth_stats_', sim_index, ".tab"), sep='\t', quote=F, row.names = F) 
write.table(performance, file=paste0('sim_performance_', sim_index, ".tab"), sep='\t', quote=F, row.names=F) 

#########################################################################################


