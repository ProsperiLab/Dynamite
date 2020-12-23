
## Initialize stable working environment and store time of intiation
rm(list=ls())



# List of packages for session
.packages <-  c("optparse", "remotes", "phytools", "treeio", "dplyr",  
                "plyr", "tidyr", "tidytree", "data.table",
                "parallel", "stringr", "rlist",  "phylodyn",
                "ggtree", "ggplot2", "TreeTools", "geiger", "adephylo", "inflection") # May need to incorporate code for familyR (https://rdrr.io/github/emillykkejensen/familyR/src/R/get_children.R) i fno longer supported.
.github_packages <- c("emillykkejensen/familyR") # "mrc-ide/skygrowth"

# # Install CRAN packages (if not already installed)
# .inst <- .packages %in% installed.packages()
# if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# .inst_github <- .packages %in% installed.packages()
# ## Install GitHub packages(if not already installed)
# if(length(.github_packages[!.inst_github]) > 0) try(remotes::install_github(.github_packages[!.inst_github]))
# if(length(.github_packages[!.inst_github]) > 0) try(devtools::install_github(.github_packages[!.inst_github]))

# Load packages into session 
lapply(.packages, require, character.only=TRUE)
lapply(gsub(".+\\/(.+)", "\\1", .github_packages), require, character.only=TRUE)
numCores <- detectCores()

option_list = list(
  make_option(c("-s", "--sim_index"), type="numeric", default=666),
  make_option(c("-a", "--cluster"), type="character", default="b", 
              help="choice of cluster algorithm from c (Phylopart's cladewise) or b (DYNAMITE's branchwise) [default= phylopart]", metavar="character"),
  make_option(c("-l", "--leaves"), type="character", default="", 
              help="choice of transformation to tree from bifurcating or addLeaves [default=empty]", metavar="character"),
  make_option(c("-t", "--threshold"), type="numeric", default=0.075, 
              help="branch length threshold [default= 0.05]", metavar="numeric")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$s)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

print(opt)

### Functions ##############################################################################################################

`%notin%` <- Negate(`%in%`) # Just plain useful
`%!=na%` <- function(e1, e2) (is.na(e1) & is.na(e2)) # Also useful for finding nas in two dataframes

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
add.test.clust <- function() {
  dirtrans_clusters <- readRDS(file=paste0("dirtrans_clusters_", opt$sim_index, ".rds"))
#  n <- length(true_clusters) # Can be replaced with below since all true clusters now (but if switching back remember to change i to n+i below)
  true_clusters <- list()
  state <- list()
  for (i in seq_along(dirtrans_clusters)){
  state[[i]] <- dirtrans_clusters[[i]]$state[1]
  true_clusters[[i]] <- data.frame(cluster_id = paste0(state[[i]], "_direct"), # This will be some arbitrary name you assign to the cluster, such as C1, C2, C3, etc.
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
    nlist <- geiger::tips(tree,node)
    foo <- distmat[nlist,nlist]
    return(median(foo[upper.tri(foo,diag=FALSE)]))
  } ## Given a node, tree, and distance matrix, return median   pairwise patristic distance (MPPD) of its leaves
  get.node.full.MPPD <- function(node,tree,distmat){
    nlist <- geiger::tips(tree, node)
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
    merge.overlap.clust() %>% # Will not work well with support in labels
    merge.nested.clust() %>%
    list.filter(length(label) >= 5)
  return(clusters)
}
benchmark <- function(clusters, true_clusters) {
# First benchmark true clusters
compareClusters <- function(x,y) { #true clusters=x, clusters=y
  state <- NA
  proportion <- NA
  additional <- NA
  true_taxa <- NA
  identified_taxa <- NA
  
  if (isTRUE(
     names(which.max(table(gsub(".+\\_([A-Z])", "\\1", unlist(str_split(x$taxa, ',')))))) %in% 
     names(table(gsub(".+\\_([A-Z])", "\\1", y$label))) &
     sum(unlist(str_split(x$taxa, ',')) %in% y$label) >= 0.70*length(unlist(str_split(x$taxa, ','))))) {
     #length(y$label))) {
     state <- names(which.max(table(gsub(".+\\_([A-Z])", "\\1", unlist(str_split(x$taxa, ','))))))
     proportion <- sum(y$label %in% unlist(str_split(x$taxa, ',')))/
       length(unlist(str_split(x$taxa, ',')))
     additional <- (length(y$label)-sum(y$label %in% unlist(str_split(x$taxa, ','))))/
       length(y$label)
     true_taxa <- x$taxa
     identified_taxa <- paste(y$label, sep=",", collapse=",")
  } else {
      state1 <-  names(which.max(table(gsub(".+\\_([A-Z])", "\\1", unlist(str_split(x$taxa, ','))))))
      proportion1=0
      additional1=NA
      true_taxa1 <- x$taxa
      identified_taxa1 <- NA
      
      state2 <-   names(which.max(table(gsub(".+\\_([A-Z])", "\\1", y$label))))
      proportion2=NA
      additional2 = 1
      true_taxa2 <- NA
      identified_taxa2 <- paste(y$label, sep=",", collapse=",")
      
      state <- rbind(state1, state2)
      proportion <- rbind(proportion1, proportion2)
      additional <- rbind(additional1, additional2)
      true_taxa <- rbind(true_taxa1, true_taxa2)
      identified_taxa <- rbind(identified_taxa1, identified_taxa2)
      
    }
  performance <- data.frame(state=state, proportion=proportion, additional=additional, 
                    true_taxa=true_taxa, identified_taxa=identified_taxa, 
                    stringsAsFactors = F)
  return(performance)
}
  
performance <- mclapply(true_clusters, function(x) lapply(clusters, function(y) compareClusters(x,y)), mc.cores=numCores)
performance <- mclapply(performance, rbindlist, mc.cores = numCores)
performance <- rbindlist(performance) %>% dplyr::distinct()

 
 ## If multiple clusters as part of the same true cluster are identified, add proportions together (because distinct)
 # performance <- lapply(performance, function(x) {
 #   state=x$state[1]
 #   taxa=x$taxa
 #   proportion=sum(x$proportion)
 #   
 #   return(data.frame(state=state, proportion, taxa=taxa, stringsAsFactors = F))
 # })
 

## Get rid of copies of 1 and 0 proportions
copy <- performance
unwanted <- NULL
i <-1
while(i<=nrow(copy)) {
  for (j in 1:nrow(performance)) {
    if(isTRUE(copy$proportion[i]==0 & performance$proportion[j]>0 &
              copy$state[i]==performance$state[j] &
              copy$true_taxa[i]==performance$true_taxa[j])) {
      unwanted <- rbind(unwanted, copy[i,])
      }
  } # End loop along performance
i=i+1
}

if(!is.null(unwanted)) {
performance <- setdiff(performance, unwanted)
} else {performance <- performance }

## Get rid of copies of 1 and 0 proportions
copy <- performance
unwanted <- NULL
i <-1
while(i<=nrow(copy)) {
  for (j in 1:nrow(performance)) {
    if(isTRUE(is.na(copy$proportion[i]) & performance$proportion[j]>0 &
              copy$state[i]==performance$state[j] &
              copy$identified_taxa[i]==performance$identified[j])) {
      unwanted <- rbind(unwanted, copy[i,])
    }
  } # End loop along performance
  i=i+1
}
if(!is.null(unwanted)) {
  performance <- setdiff(performance, unwanted)
} else {performance <- performance }
return(performance)
}
calculateNe <- function(cluster, tree) {
  if (isTRUE("phylo" %in% class(cluster))) {  
    fit <- ltt(tree, plot=F)
    } else {
      if (isTRUE("tbl" %in% class(cluster) | "data.table" %in% class(cluster))) {
      x <- keep.tip(tree, cluster$label[cluster$label %in% tree$tip.label])
     #  res=round(max(nodeHeights(x))/days) # Total Ne set to every two weeks, but setting clusters up for every week
      #  fit <- tryCatch(skygrowth.map(x, res=res, tau0=0.1), error=function(e) NULL)
      #  p <- data.frame(time=fit$time, nemed=fit$ne, ne_ci=fit$ne_ci)
      fit <- tryCatch(ltt(x, plot=F), error=function(e) NULL)
      }
    }
  return(fit)
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
      Re.dist = sample(1+psi*(Ne$ne[i]-Ne$ne[i-1])/((Ne$time[i]-Ne$time[i-1])*Ne$ne[i-1]),
                       size=s, replace=F)
      logRe = log(Re.dist)
      SElogRe = sd(logRe)/sqrt(s) # standard deviation of 2, sample size of 10
      LCL = exp(mean(logRe) - Z*SElogRe) 
      UCL = exp(mean(logRe) + Z*SElogRe) 
      Re[[i-1]] <- data.frame(time = time, mean_Re=mean(Re.dist), conf.int=paste0("(", LCL, "," ,UCL, ")"))
      } 
      } else {Re[[1]] <- NULL}
  Re <- do.call("rbind",Re)
}
# calculateOster <- function(cluster, tree) {
#   
#   if (isTRUE("phylo" %in% class(cluster))) {  
#     cluster_size <- length(cluster$tip.label)-1
#     sum_heights <- sum(nodeHeights(cluster))
#     longest <- max(nodeHeights(cluster))
#     tr <- cluster_size/sum_heights + longest
#   } else {
#     if (isTRUE("tbl" %in% class(cluster) | "data.table" %in% class(cluster))) {
#       x <- keep.tip(tree, cluster$label[cluster$label %in% tree$tip.label])
#       cluster_size <- length(x$tip.label)-1
#       sum_heights <- sum(nodeHeights(x))
#       longest <- max(nodeHeights(x))
#       tr <- cluster_size/sum_heights + longest
#     }
#   }
#   return(tr)
# }

#tr <- lapply(clusters, function(x) calculateOster(x, tree@phylo))
  
  
modelEffective <- function() {
  
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
  
  total_Ne_df <- data.frame(ne=total_Ne$ltt, time=total_Ne$times) %>% distinct()
  
  clusters_Ne <- lapply(clusters, function(x) calculateNe(x, tree@phylo))
  clusters_gamma <- mclapply(clusters_Ne, function(x) return(x$gamma))
  assign("clusters_gamma", clusters_gamma, envir = globalenv())
  
  clusters_Ne_df <- mclapply(clusters_Ne, function(x) {
    distinct(data.frame(ne=x$ltt, time=x$times))
  }, mc.cores=numCores)
  
  # modelNe_log <- function(Ne_df) {
  #     df <- data.frame(x=Ne_df$time, y=Ne_df$ne)
  #     growth_fit <- nls(y~I(K / (1 + ((K - N0) / N0) * exp(-r * x))), data=df,
  #                       start=list(N0=df$y[1], K=100, r=0.10), algorithm = "port",
  #                       lower = list(N0 = 1, K = 5, r=0),
  #                       upper = list(N0 = 5, K=1000, r=10))
  # 
  #     s_growth_fit <- summary(growth_fit)
  #     model.coef <- NULL
  #     model.coef <- rbind.data.frame(model.coef,
  #                                         data.frame(N0=s_growth_fit$coefficients[1],
  #                                                    K=s_growth_fit$coefficients[2],
  #                                                    r=s_growth_fit$coefficients[3]))
  # 
  #     growth_fx <- function(x, N0, K, r) {
  #        y = K / (1 + ((K - N0) / N0) * exp(-r * x))
  #        if (isTRUE(y<0)) {
  #          y = 0
  #        } else {y = y}
  #        return(y)}
  #     growth_fx1 <- function(x) growth_fx(x, model.coef$N0, model.coef$K,
  #                                         model.coef$r)
  #     tmp <- data.frame(y=growth_fx1(df$x), x=df$x)
  #     r2 <- summary(lm(df$y~tmp$y))$adj.r.squared
  #     
  #     new_times <- seq(min(df$x), max(df$x), 1)
  #     modeled_growth <- data.frame(time=new_times, Ne=sapply(new_times, growth_fx1), model="logistic")
  #     modeled_growth <- cbind(modeled_growth, model.coef)
  #     return(modeled_growth)
  # } # End function
  # modelNe_decay <- function(Ne_df) {
  #   df <- data.frame(x=Ne_df$time, y=Ne_df$ne)
  #   growth_fit <- lm(df$y ~ sin(2*pi/df$x)+cos(2*pi/df$x))
  #   model.coef <- NULL
  #   model.coef <- rbind.data.frame(model.coef,
  #                                  data.frame(Beta=growth_fit$coefficients[1], 
  #                                             k1=growth_fit$coefficients[2],
  #                                             k2=growth_fit$coefficients[3]))
  #   growth_fx <- function(x, k1, k2, Beta){
  #     y = k1*sin(2*pi/x)+ k2*cos(2*pi/x) + Beta
  #     if (isTRUE(y<0) ){
  #       y = 0
  #     } else {y = y}
  #     return(y)}
  #   growth_fx1 <- function(x) growth_fx(x, model.coef$k1, model.coef$k2, model.coef$Beta)
  #   new_times <- seq(min(df$x), max(df$x), 1)
  #   modeled_growth <- data.frame(time=new_times, Ne=sapply(new_times, growth_fx1), model="decay")
  #   if (is.nan(modeled_growth$Ne[1])) {
  #     modeled_growth <- modeled_growth[-1,]
  #   }
  #   max_ne <- max(modeled_growth$Ne)
  #   min_ne <- min(modeled_growth$Ne)
  #   max_t <- modeled_growth$time[modeled_growth$Ne==max_ne]
  #   min_t <- min(modeled_growth$time)
  #   r = (max_ne-min_ne)/(max_t-min_t)
  #   r2 <- summary(growth_fit)$adj.r.squared
  #   
  #   modeled_growth <- cbind(modeled_growth, model.coef, r, r2)
  #   return(modeled_growth)
  # } # End function
  # modelNe_constant <- function(Ne_df) {
  #   df <- data.frame(x=Ne_df$time, y=Ne_df$ne)
  #   growth_fit <- nls(y ~ I(a * x^b), data=df, start=list(b=0.1, a=1), algorithm = "port",
  #                     lower=list(b=0, a=0.1),
  #                     upper=list(b=1, a=10), nls.control(warnOnly=T))
  #   s_growth_fit <- summary(growth_fit)
  #   model.coef <- NULL
  #   model.coef <- rbind.data.frame(model.coef,
  #                                  data.frame(b=s_growth_fit$coefficients[1], 
  #                                             a=s_growth_fit$coefficients[2]))
  #   
  #   growth_fx <- function(x, b, a) {
  #     y = a * x^b
  #     if (isTRUE(y<0)) {
  #       y = 0
  #     } else {y = y}
  #     return(y)}
  #   growth_fx1 <- function(x) growth_fx(x, model.coef$b, model.coef$a)
  #   tmp <- data.frame(y=growth_fx1(df$x), x=df$x)
  #   r2 <- summary(lm(df$y~tmp$y))$adj.r.squared
  #   new_times <- seq(min(df$x), max(df$x), 1)
  #   modeled_growth <- data.frame(time=new_times, Ne=sapply(new_times, growth_fx1), model="constant")
  #   r <- list()
  #   for (i in 2:nrow(modeled_growth)) {
  #     r[[i]] <- (modeled_growth$Ne[i]-modeled_growth$Ne[i-1])/(modeled_growth$time[i]-modeled_growth$time[i-1])
  #   }
  #     r <- max(do.call("rbind",r), na.rm=T)
  #   modeled_growth <- cbind(modeled_growth, r, r2)
  #   return(modeled_growth)
  # } # End function
  # modelNe_growth <- function(Ne_df) {
  #   df <- data.frame(x=Ne_df$time, y=Ne_df$ne)
  #   growth_fit <- nls(y ~ I(r * x^b), data=df, start=list(b=2, r=1), algorithm = "port",
  #                     lower=list(b=2, r=0.01),
  #                     upper=list(b=10, r=10), nls.control(warnOnly=T))
  #   s_growth_fit <- summary(growth_fit)
  #   model.coef <- NULL
  #   model.coef <- rbind.data.frame(model.coef,
  #                                  data.frame(b=s_growth_fit$coefficients[1], 
  #                                             r=s_growth_fit$coefficients[2]))
  #   
  #   growth_fx <- function(x, b, r) {
  #     y = r * x^b
  #     if (isTRUE(y<0)) {
  #       y = 0
  #     } else {y = y}
  #     return(y)}
  #   growth_fx1 <- function(x) growth_fx(x, model.coef$b, model.coef$r)
  #   tmp <- data.frame(y=growth_fx1(df$x), x=df$x)
  #   r2 <- summary(lm(df$y~tmp$y))$adj.r.squared
  #   new_times <- seq(min(df$x), max(df$x), 1)
  #   modeled_growth <- data.frame(time=new_times, Ne=sapply(new_times, growth_fx1), model="growth")
  #   r <- list()
  #   for (i in 2:nrow(modeled_growth)) {
  #     r[[i]] <- (modeled_growth$Ne[i]-modeled_growth$Ne[i-1])/(modeled_growth$time[i]-modeled_growth$time[i-1])
  #   }
  #   r <- max(do.call("rbind",r), na.rm=T)
  #   modeled_growth <- cbind(modeled_growth, r, r2)
  #   return(modeled_growth)
  # } # End function
  # 
  # modelFit <- function(x) {
  #   r2_growth <- tryCatch(modelNe_growth(x), error=function(e) {return(data.frame(r2=0))})
  #   r2_constant <- tryCatch(modelNe_constant(x), error=function(e) {return(data.frame(r2=0))})
  #   if(isTRUE(r2_growth$r2[1] >= r2_constant$r2[1])) {
  #     result <- r2_growth
  #   } else {
  #     result <- r2_constant
  #   }
  #   r2_decay <- tryCatch(modelNe_decay(x), error=function(e) {return(data.frame(r2=0))})
  #   if(isTRUE(r2_decay$r2[1] >= r2_constant$r2[1])) {
  #     result2 <- r2_decay
  #   } else {
  #     result2 <- result
  #   }
  #   return(result)
  # }
  # 
  # total_Ne_model <- modelFit(total_Ne_df)
  # assign("total_Ne_model", total_Ne_model, envir = globalenv())
  # 
  # clusters_Ne_model <- lapply(clusters_Ne_df, modelFit)
   
   modelNe_growth <- function(Ne_df) {
     df <- data.frame(x=Ne_df$time, y=Ne_df$ne)
     df <- df[-1,]
     
      repeat {
        i <- 2
       while (i <= nrow(df)) {
       if (isTRUE(df$x[i]==df$x[i-1])){
         min_y <- min(df$y[i-1],df$y[i])
         unwanted <- df[df$y== min_y & (df$x==df$x[i] | df$x==df$x[i-1]),]
       } else{ unwanted <- data.frame(x=NA, y=NA)}
         df <- setdiff(df, unwanted)
       i <- i+1
       }
        if(length(unique(df$x))==length(df$x)) {
          break
        }
      }


     cc <- check_curve(df$x, df$y)$ctype
     if (cc=="convex") {
       model="constant"
     } else {
       if (cc=="concave") {
         model="growth"
       } else { model="mixed"}
     }
     
     rmax <- list()
     for (i in 2:nrow(df)) {
       rmax[[i]] <- (df$y[i]-df$y[i-1])/(df$x[i]-df$x[i-1])
     }
     rmax <- max(do.call("rbind",rmax), na.rm=T)
     model_data <- data.frame(model=model, rmax=1/rmax)
     Ne_df <- data.frame(time=df$x, ne=df$y)
    return(cbind(Ne_df, model_data))
   }
   
   total_Ne_model <- modelNe_growth(total_Ne_df)
   assign("total_Ne_model", total_Ne_model, envir = globalenv())
   clusters_Ne_model <- lapply(clusters_Ne_df, modelNe_growth)
  c_root_age <- list()
  for (i in seq_along(clusters_Ne_model)) {
    c_root_age[[i]] <- unique(subset(metadata$root_age, metadata$cluster_id == names(clusters_Ne_df)[i] & !is.na(metadata$cluster_id)))
    clusters_Ne_model[[i]]$time = clusters_Ne_model[[i]]$time + c_root_age[[i]]
  }
  assign("clusters_Ne_model", clusters_Ne_model, envir = globalenv())
  
  
  clusters_Re <- mclapply(clusters_Ne_model, function(x) {
    calculateRe(x, conf.level = 0.95)}, mc.cores = numCores)
  assign("clusters_Re", clusters_Re, envir = globalenv())
} # End modelNe function

growthMetrics <- function(clusters, clusters_Ne_model, performance) {
  
  growthPeriod <- function(ne_model) {
    max_ne <- max(ne_model$ne)
    max_t <- ne_model$time[ne_model$ne==max_ne][1]
    t_frac <- (max_t-min(ne_model$time))/(max(ne_model$time)-min(ne_model$time))
    rmax <- ne_model$rmax[1]
    
#    t1 <- min(ne_model$time)
#    t2 <- max(ne_model$time)
#    ne1 <- min(ne_model$ne[ne_model$time == t1])
#    ne2 <- max(ne_model$ne[ne_model$time == t2])
#    r <- (ne2-ne1)/(t2-t1)
    model <- ne_model$model[1]
    growth_rate <- data.frame(t_frac=t_frac, rmax=rmax, model=model)
    if (model=="growth") {
      df <- ne_model[ne_model$time >= max_t,]
      lin.fit <- lm(df$ne~df$time)
      if (summary(lin.fit)$coefficients[2] <0 & summary(lin.fit)$coefficients[2,4] < 0.05) {
        model="decay"
      } else {model=model}
    }
    return(data.frame(t_frac=t_frac, rmax=rmax, model=model))
    }
  
  total_Ne_growth <- growthPeriod(total_Ne_model)
  cluster_Ne_growth <- mclapply(clusters_Ne_model, growthPeriod)
  
  

  growth_criteria <- NULL
  ## Calculate slopes for Ne and ltt, and record R0
  for (i in seq_along(clusters)) {
    state <- dplyr::filter(performance, identified_taxa == paste(clusters[[i]]$label, sep=",", collapse=","))$state
    identified_taxa = paste(clusters[[i]]$label, sep=",", collapse=",")
    cluster_size = length(clusters[[i]]$label)
    ## Find amount of time spent in growth phase (based on maximum)
    model <- cluster_Ne_growth[[i]]$model
    max_r <- cluster_Ne_growth[[i]]$rmax
    rel_max_r <- max_r/total_Ne_growth$rmax
    
 #   r <- cluster_Ne_growth[[i]]$r
#    rel_r <- r/total_Ne_growth$r
    
    t_frac <- cluster_Ne_growth[[i]]$t_frac
    ## Report whether growth occurred 
    mean_R0 <- clusters_Re[[i]]$mean_Re[1]
    lower_R0 <- as.numeric(gsub("\\((.+)\\,.+", "\\1", as.character(clusters_Re[[i]]$conf.int[1])[1]))
    upper_R0 <- as.numeric(gsub(".+\\,(.+)\\)", "\\1", as.character(clusters_Re[[i]]$conf.int[1])[1]))
    gamma <- clusters_gamma[[i]]
    growth_criteria <- rbind(growth_criteria, 
                             data.frame(state = state,
                                        identified_taxa=identified_taxa,
                                        size = cluster_size,
                                        model=model,
                                        t_frac=t_frac,
                                        max_r = max_r,
                                        rel_max_r = rel_max_r,
#                                        r = r,
#                                        rel_r = rel_r,
                                        pybus_gamma = gamma,
                                        mean_R0 = mean_R0,
                                        lower_R0 = lower_R0,
                                        upper_R0 = upper_R0, 
                                        dynamic=NA,
                                        priority=NA,
                                        stringsAsFactors = F))
    } # End loop along clusters

  # Label as growing, not growing, or dead according to ltt slope
  for (i in 1:nrow(growth_criteria)) {
    if (isTRUE(growth_criteria$model[i] == "constant")){
      growth_criteria$dynamic[i] <- "static"
    } else {
      if (isTRUE(growth_criteria$model[i] == "mixed" | growth_criteria$model[i] == "growth")){
        growth_criteria$dynamic[i] <- "growth"
    } else {
      if (isTRUE(growth_criteria$model[i] == "decay")){
        growth_criteria$dynamic[i] <- "decay"
    }
    }
    }
    n <- length(unique(metadata$ID))/10000*100
    if (isTRUE(growth_criteria$rel_max_r[i] > n)){
      growth_criteria$priority[i] <- "high priority"
    } else {
      if (isTRUE(growth_criteria$rel_max_r[i] <= n)){
        growth_criteria$priority[i] <- "low priority"
    }
    }
  }
  performance <- merge(performance, growth_criteria, by=c("state", "identified_taxa"), all.x=T)
  return(performance)
}
connectClust <- function(sub_tree, clusters) {
  
  cluster_data <- clusters
  dup_cluster_data <- cluster_data
  full_tree <- as_tibble(sub_tree)
  
  for (i in seq_along(cluster_data)) {
    
    for (j in seq_along(dup_cluster_data[-i])) {
      dup_cluster_data[[j]]$birth_origin <- NA
      #      dup_cluster_data[[j]]$birth_origin <- NA
      ## If the parent of origin node of one cluster (found by referencing full tree) is a child node in one of the other clusters... 
      if (isTRUE(tidytree::parent(full_tree, dup_cluster_data[[j]]$parent[1])$node %in% cluster_data[[i]]$node |
                 ## Or if the parent of the parent of the origin node of one cluster (found by referencing full tree) is a child node in one of the other clusters... (meaning clusters are allowed to be separated by up to two branches)
                 tidytree::parent(full_tree, dup_cluster_data[[j]]$parent[1])$parent %in% cluster_data[[i]]$node)) {
        #        dup_cluster_data[[j]]$birth_origin == paste0("cluster_", cluster_data[[i]]$parent[1] )
        dup_cluster_data[[j]]$birth_origin = paste0(names(cluster_data)[[i]])
      } #End if statement
    } # End for loop along duplicate list
  } # End for loop along original list
  return(dup_cluster_data)
} #End connectClust function
identifyBirths <- function(clusters, sub_tree, performance) {
  clusters <- connectClust(sub_tree, clusters)
  for (i in seq_along(clusters)) { ## Update performance
    if(isTRUE(!is.na(clusters[[i]]$birth_origin[1]))) {
      performance$birth[i] <- "birth"
    } else {performance$birth[i] <-NA}
  }
  return(performance)
} # End of function
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
identifyDeaths <- function(clusters, sub_tree, performance) {
  contemp_taxa <- getContempTaxa(sub_tree)
  
  for (i in seq_along(clusters)) { ## Update performance
    if (isTRUE(any(clusters[[i]]$label %in% contemp_taxa))) {
      performance$death[i] <- NA
    } else{performance$death[i] <- "death"} # End if-else statement
  } # End for loop
  return(performance)
}
quantifyR0perf <- function(performance, cluster_dynamics) {
  performance$R0_diff<- NA
  performance$R0_overlap <- NA
  
  for (i in 1:nrow(performance)) {
    for (j in 1:nrow(cluster_dynamics)) {
      
      # R0_diff
      if (isTRUE(performance$state[i] == cluster_dynamics$state[j])) {
        performance$R0_diff[i] <- abs(performance$mean_R0[i] - cluster_dynamics$mean_R0[j])
      } else {NULL}
      
      # R0 intervals           
      if (isTRUE(performance$state[i] == cluster_dynamics$state[j] &
                 (performance$lower_R0[i] <= cluster_dynamics$lower_R0[j] &
                  performance$upper_R0[i] >= cluster_dynamics$upper_R0[j]) |
                 (performance$lower_R0[i] >= cluster_dynamics$lower_R0[j] &
                  performance$upper_R0[i] <= cluster_dynamics$upper_R0[j]))){
        performance$R0_overlap[i] <- "included"
      } else {
        if (isTRUE(performance$state[i] == cluster_dynamics$state[j] &
                   (performance$lower_R0[i] <= cluster_dynamics$lower_R0[j] & 
                    performance$upper_R0[i] >= cluster_dynamics$lower_R0[j]) |
                   (performance$upper_R0[i] >= cluster_dynamics$upper_R0[j] &
                    performance$lower_R0[i] <= cluster_dynamics$upper_R0[j]))){# If true cluster dynamic classification is contained within the reported cluster dynamic,
          performance$R0_overlap[i] <- "overlapping"
        } else {
          if (isTRUE(is.na(performance$mean_R0[i]))) {
            performance$R0_overlap[i] <- NA
          } else {performance$R0_overlap[i] <- "non-overlapping"}
        }
      }
    }
  }
  return(performance)
} # End function
classify <- function(performance,cluster_dynamics) { #Where x is performance and y is cluster_dynamics
  x <- performance #Just easier for what follows
  z <- cluster_dynamics # Just easier for what follows
  x$dynamic_TF <- NA
  x$birth_TF <- NA
  #    x$death_TF <- NA
  x$priority_TF <- NA
  for (i in 1:nrow(x)) {
    for (k in 1:nrow(z)) {
      if (x$state[i] == z$state[k]) { # For clusters in both tables (identified true clusters)
        x$dynamic_TF[i] <- if_else(isTRUE(x$dynamic[i] == z$dynamic[k] | pmatch(z$dynamic[k], x$dynamic[i]) ==1 |
                                            (x$dynamic[i] == "static" & z$dynamic[k] == "background")), "TRUE", "FALSE")
        x$birth_TF[i] <- if_else(isTRUE(x$birth[i] == z$birth[k] |x$birth[i] %!=na% z$birth[k]), "TRUE", "FALSE")
        #         x$death_TF[i] <- if_else(isTRUE(x$death[i] == z$death[k] |x$death[i] %!=na% z$death[k]), "TRUE", "FALSE")
        x$priority_TF[i] <- if_else(isTRUE((z$dynamic[k] == "background" & x$priority[i] == "low priority") |
                                             (z$dynamic[k] != "background" & x$priority[i] == "high priority")), "TRUE", "FALSE")
      } else {NULL}
    } #End loop along z
    if (isTRUE(is.na(x$identified_taxa[i]) | x$model[i]=="CNBD")) {
      x$dynamic_TF[i] <- NA
      x$birth_TF[i] <- NA
      x$priority_TF[i] <- NA
    }
  } # End loop along x
  return(x)
} # End function


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
# true_clusters <- find.true.clusters(family_tree, clades) %>%
#   merge.overlap.clust.true() %>%
#   list.filter(length(unlist(str_split(taxa, ','))) >= 6)
## Add test "A" clusters, which is a clade of sequences taken directly from the full transmission tree
true_clusters <- add.test.clust()
# states_present <- lapply(true_clusters, function(x) {
#   (gsub(".+_([A-Z]).+", "\\1", x$taxa))
# })
# states_present <- data.frame(sim=opt$sim_index, state = do.call("rbind", states_present), stringsAsFactors = F) %>%
#   group_by(sim, state) %>%
#   dplyr::summarise(state_count = n())

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

## Do not continue if no clusters found ############################################################
if (length(clusters)==0) {
  print("No clusters detected. Saving null output...")
  performance <- lapply(true_clusters, function(x) {
    state <- names(which.max(table(gsub(".+\\_([A-Z])", "\\1", unlist(str_split(x$taxa, ','))))))
    identified_taxa <- NA
    proportion <- 0
    additional <- NA
    true_taxa <- x$taxa
    size <- NA
    model <- NA
    t_frac <- NA
    max_r <- NA
    rel_max_r <- NA
    pybus_gamma <- NA
    mean_R0 <- NA
    lower_R0 <- NA
    upper_R0 <- NA
    dynamic <- NA
    priority <- NA
    birth <- NA
    death <- NA
    R0_overlap <- NA
    R0_diff <- NA
    dynamic_TF <- NA
    birth_TF <- NA
    priority_TF <- NA
    performance <- data.frame(state=state, identified_taxa=identified_taxa, proportion=proportion, additional=additional,
                              true_taxa=true_taxa, size=size, model=model, t_frac=t_frac, max_r=max_r, rel_max_r=rel_max_r,
                              pybus_gamma=pybus_gamma, mean_R0=mean_R0, lower_R0=lower_R0,
                              upper_R0=upper_R0, dynamic=dynamic, priority=priority, birth=birth, death=death,
                              R0_overlap=R0_overlap, R0_diff=R0_diff, dynamic_TF=dynamic_TF, birth_TF=birth_TF, priority_TF=priority_TF)
    return(performance)
  })
  performance <- do.call(rbind, performance)
  metadata <- NULL
} else {
  print("Benchmarking overall performance")
  performance <- benchmark(clusters, true_clusters)
  print("Calculating Ne(t) and Re(t)...")
  modelEffective()
  # ggplot(bind_rows(clusters_Ne_model, .id="df"), aes(x=time, y=Ne, colour=df)) +
  #   geom_line()
  
  # relative_Ne <- mclapply(clusters_Ne_model, function(x) {
  #   Ne_merged <- merge(total_Ne_model, x, by="time", all.x=T)
  #   Ne_merged$rel_Ne <- Ne_merged$Ne.y/Ne_merged$Ne.x
  #   return(dplyr::select(Ne_merged, time, rel_Ne))
  # }, mc.cores=numCores)
  # 
  # ggplot(bind_rows(relative_Ne, .id="df"), aes(x=time, y=rel_Ne, colour=df)) +
  #   geom_line()
  
  print("Calculating growth metrics...")
  performance <- growthMetrics(clusters, clusters_Ne_model, performance)
  
  print("Identifying births...")
  performance <- identifyBirths(clusters, sub_tree, performance)
  
  print("Identifying deaths...")
  performance <- identifyDeaths(clusters, sub_tree, performance)
  
  print("Quantifying R0 estimation performance...")
  performance <- quantifyR0perf(performance, cluster_dynamics)
  
  print("Classifying correct/incorrect inferences...")
  performance <- classify(performance, cluster_dynamics)
  
}# end if-else function for if zero clusters found vs actual clusters found


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

print("Saving results...")
#saveRDS(clusters, file=paste0("clusters_", opt$sim_index, ".rds"))
#saveRDS(true_clusters, file=paste0("true_clusters_", opt$sim_index, ".rds"))


performance$sim <- opt$sim_index
metadata$sim <- opt$sim_index

#write.table(states_present, file=paste0('num_true_clusters_', opt$sim_index, ".tab"), sep='\t', quote=F, row.names = F)
write.table(metadata, file=paste0('metadata_', opt$sim_index, ".tab"), sep='\t', quote=F, row.names = F)
write.table(performance, file=paste0('sim_performance_', opt$sim_index, ".tab"), sep='\t', quote=F, row.names=F) 

#########################################################################################
## Tree export
# gg_tree <- ggtree(tree, mrsd=)
# 
# 
# exportTree <- function(time_tree, cluster_data) {
#   t.tbl <- as_tibble(time_tree)
#   t.tbl <- cbind(t.tbl, gg_tree$data$x)
#   names(t.tbl)[5] = "dates"
#   t.tbl$cluster_id <- ""
#   for (i in seq_along(cluster_data)) {
#     for (j in seq_along(t.tbl$node)) {
#       if (t.tbl$node[j] %in% cluster_data[[i]]$node) {
#         t.tbl$cluster_id[j] = names(cluster_data)[i]
#       } else{t.tbl$node[j] = t.tbl$node[j]}
#     }
#   }
#   
#   t.tbl$dates <- date_decimal(t.tbl$dates)
#   t.tbl$dates <- as.Date(t.tbl$dates)
#   t.tbl$label <- paste(t.tbl$label, t.tbl$dates, sep="|") # Needed for nextstrain, but ignoring for now
#   class(t.tbl) = c("tbl_tree", class(t.tbl))
#   
#   t2 <- as.treedata(t.tbl)
#   return(t2)
# }
# exported_tree <- exportTree(time_tree, cluster_data)
# write.beast(exported_tree, "tree_data.tree")
# 
# rt1 <- Sys.time()
# rt1-rt0


