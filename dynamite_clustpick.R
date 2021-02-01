
##Install and load necessary packages
## Initialize stable working environment and store time of intiation
rm(list=ls())
#setwd(getSrcDirectory()[1]) #When working on the cluster
#dirname(rstudioapi::getActiveDocumentContext()$path) # If working in Rstudio
rt0 <- Sys.time()

# List of packages for session
.packages <-  c("optparse", "remotes", "phytools", "data.tree", 
                "tidytree", "lubridate", "rlist", "familyR", "tidyverse", 
                "ggtree", "parallel", "geiger", "tibble")  # May need to incorporate code for familyR (https://rdrr.io/github/emillykkejensen/familyR/src/R/get_children.R) i fno longer supported.
#github_packages <- c("tothuhien/Rlsd2") # mrc-ide/skygrowth, "mdkarcher/phylodyn" may need to be done if we think our Re values are going to be greater than 5 for any cluster! If the latter, need aslo install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

## Will need to remove install section if using on cluster ###################################
# Install CRAN packages (if not already installed)
#.inst <- .packages %in% installed.packages()
#if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
#.inst_github <- .packages %in% installed.packages()
## Install GitHub packages(if not already installed)
#if(length(github_packages[!.inst_github]) > 0) try(remotes::install_github(github_packages[!.inst_github]))
#if(length(github_packages[!.inst_github]) > 0) try(devtools::install_github(github_packages[!.inst_github]))
## Will need to remove install section if using on cluster ###################################

# Load packages into session 
lapply(.packages, require, character.only=TRUE)
#lapply(gsub(".+\\/(.+)", "\\1", github_packages), require, character.only=TRUE)

option_list = list(
  make_option(c("-t", "--tree"), type="character", default=list.files(pattern = "*.nwk")[[1]], 
              help="tree file name [default= .nwk extension]", metavar="character"),
  # make_option(c("-m", "--metadata"), type="character", default=list.files(pattern="*.csv")[[1]], 
  #             help="metadata file name [default= .csv extension]", metavar="character"),
  # make_option(c("-s", "--seqLen"), type="numeric", default=30000, 
  #             help="sequence length [default= 10000]", metavar="numeric"),
  make_option(c("-c", "--cluster"), type="character", default="b", 
              help="choice of cluster algorithm from c (Phylopart's cladewise) or b (DYNAMITE's branchwise) [default= dynamite]", metavar="character"),
  make_option(c("-l", "--leaves"), type="character", default="addLeaves", 
              help="choice of transformation to tree from bifurcating or addLeaves [default=addLeaves]", metavar="character"),
  make_option(c("-o", "--threshold"), type="numeric", default=0.10, 
              help="branch length threshold [default= 0.10]", metavar="numeric")
  # make_option(c("-a", "--asr"), type="character", default="Y", 
  #             help="option of ancestral state reconstruction for each cluster [default= Y]", metavar="character")
 ); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$tree)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Additional options
`%notin%` <- Negate(`%in%`) # Just plain useful
# The option below is useful when dealing with dates of internal nodes downstream
options(digits=15)


# Supply tree and metadata table as arguments

numCores = try(Sys.getenv("SLURM_CPUS_ON_NODE"))
if (numCores == "") {
  numCores = detectCores()
}

# Script will read in tree from directory based on a number of format
# Function to read Nexus files or Newick files 
message("Checking for tree file...")
checkFortree <- function(tree_file) {
  print("Make sure tree_file is scaled in substitutions/site and not in units of time.")
  # Check tree_file format -- either Newick or Nexus -- using suffices
  if (endsWith(tree_file, "nwk") || endsWith(tree_file,"newick") || endsWith(tree_file,"tree")){
    print(paste("Newick file detected:", tree_file, sep=" "), stderr())
    sub_tree <- read.tree(tree_file)
  } else if (endsWith(tree_file, "nex") || endsWith(tree_file,"nexus") || endsWith(tree_file, "nxs")){
    print("Nexus file detected")
    sub_tree <- read.tree(tree_file) ##Note: will not read in bracketed annotations!
    return(tree_file)
  }  else {
    # Neither Newick nor Nexus identified -- stop
    print("Cannot identify tree_file format. Use nexus (nex,nxs,nexus), newick (nwk,newick), our BEAST output (tre, tree_files) formats.")
  }# end conditional statement
    #assign("sub_tree_file", sub_tree_file, envir=globalenv())
  return(sub_tree)
}# end checkfortree_file function
sub_tree <- checkFortree(opt$tree)

# Find all well-supported clades
message("Defining well-supported clades within the tree using node labels...")
define.clades <- function(sub_tree) {
  ## Grab only subtrees that are supported by bootstrap values >90
  ## We may need to change this in case people have other support values
  ## Note that subtrees() function in ape renumbers nodes, so useless here, since at the end we wish to recombine the information
  message("Creating list of all subtrees...This may take a while.")
  
  if (isTRUE(any(grepl("/", sub_tree$node.label)))) {
    sub_tree$node.label <- gsub("(.+)\\/(.+)", "\\2", sub_tree$node.label)
  }
  
  
  family_tree <- tidytree::as_tibble(sub_tree)
  ## Need to relabel columns so that "parent" and "node" are "from" and "to" for familyR::get_children function
  colnames(family_tree)[1:2] <- c("from", "to")
  ## The dataframe needs to be transformed back into a data.table for future analyses
  family_tree <- data.table::as.data.table(family_tree)
  assign("family_tree", family_tree, envir = globalenv())
  
  message("Extract well-supported subtrees.")
  supported_nodes <- family_tree
  # supported_nodes <- list()
  ## Need to account for multiple support values (defaulting to second)
  # if(max(as.numeric(sub_tree$node.label), na.rm=T) >1) {
  # supported_nodes <- suppressWarnings(filter(family_tree, as.numeric(label) > 90))
  # } else {
  #   supported_nodes <- suppressWarnings(filter(family_tree, as.numeric(label) > 0.90))
  # }
  assign("supported_nodes", supported_nodes, envir = globalenv())  
  
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
  clades <- mclapply(clades, function(x) {
    dplyr::select(x, from, everything()) %>%
      arrange(level) %>%
      filter(!level==0)
  }, mc.cores=numCores) %>%
    list.filter(from[1] != rootnode(sub_tree)) # Whole tree is included because support is alwasy 100% for root, so discard# End loop along clades
  
  # Save to global environment or merge next function.
  assign("clades", clades, envir=globalenv())

} # End defineClades function
define.clades(sub_tree)

# Using a pre-order (root-to-tip) tree traversal, for each node, find subtrees/clades (>=3 nodes) for which cumulative mean branch, or edge, length is within the branch length limit at that level within the tree using the "branch-wise" algorithm, or "b". Note this should only be used with epidemics that are past the exponential phase of growth.

# Alternatively, using a depth-first algorithm developed by Prosperi et al., (2011), find subtrees/clades (>=2 nodes) for which median pairwise patristic distances (MPPD) between nodes is within 1-5% of the distribution of MPPDs of clades within the entire tree. This is referred to as the "clade-wise algorithm", or "c".

## Functions used in DYNAMITE cluster-picking algorithm ############################################################################
merge.nested.clust <- function(clusters) {
  copy <- clusters
  result <- list()
  unwanted <- list()
  for (ct in seq_along(clusters)) {
    for (cc in seq_along(copy)) {
      if (isTRUE(all(clusters[[ct]]$to %in% copy[[cc]]$to) &
                 length(clusters[[ct]]$to) != length(copy[[cc]]$to))) {
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
      if (isTRUE(sum(copy[[cc]]$to %in% clusters[[ct]]$to) > 0.50*length(copy[[cc]]$to)) &
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
  # pdist.clades <- function(clades, tree, distmat, mode=c('leaf', 'all')){
  #   mode <- match.arg(mode)
  #   
  #   if(mode=='leaf'){
  #     distmat <-  p.dist.mat.leaves
  #     mclapply(clades, function(x) {
  #       get.node.leaf.MPPD(x$from[1], tree, distmat)
  #     }, mc.cores=numCores)
  #   } else{
  #     distmat <-  dist.nodes(tree)
  #     mclapply(clades, function(x) {
  #       get.node.full.MPPD(x$from[1], tree, distmat)
  #     }, mc.cores=numCores)
  #   }
  # } ## Determine MPPD for all well-supported clades
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
  # clade_MPPD <- pdist.clades(clades, tree, mode='leaf')
  # assign("clade_MPPD", clade_MPPD, envir=globalenv())
  
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


message("Picking clusters based on user choice of algorithm....")
branch_length_limit <- branchLengthLimit(sub_tree)
if (opt$cluster == "b") {
  clusters <- branchWise(sub_tree, branch_length_limit, make_tree=opt$leaves)
} else {
  if (opt$cluster == "c") {
    clusters <- phylopart(sub_tree, branch_length_limit)
    clusters <- mclapply(clusters, function(x) dplyr::rename(x, parent=from, node=to), mc.cores=numCores)
  } else {
    message("Incorrect cluster_picking algorithm choice. Please choose between 'b' (branch-wise) or 'c' (clade-wise) and run script again.")
  }
}


# Visualization is based on the resulting list (per cluster) of data

message("Data are now being exported as 'cluster_data.RDS' and 'tree_data.tree.'")
#saveRDS(cluster_stats, "cluster_data.RDS")
saveRDS(clusters, paste("clusters", opt$threshold, ".RDS", sep="_"))


# test <- cluster_list[[1]][[1]][[2]]
# test$trait <- as.factor(test$trait)
# ggplot(data=test, aes(x=date, color=trait)) +
#   geom_histogram(aes(y=..count..))

## Tree export

# exportTree <- function(tree, cluster_data) {
#   t.tbl <- as_tibble(tree)
#   t.tbl <- cbind(t.tbl, gg_tree$data$x) %>%
#     dplyr::rename("dates" = "gg_tree$data$x")
#   t.tbl$cluster_id <- NA
#   for (i in seq_along(cluster_data)) {
#     for (j in seq_along(t.tbl$node)) {
#       if (t.tbl$node[j] %in% cluster_data[[i]]$node) {
#         t.tbl$cluster_id[j] = names(cluster_data)[i]
#       } 
#     }
#   }
# 
#   t.tbl$dates <- date_decimal(t.tbl$dates)
#   t.tbl$dates <- as.Date(t.tbl$dates)
#   t.tbl$label <- paste(t.tbl$label, t.tbl$dates, sep="|")
#   t.tbl <- dplyr::select(t.tbl, parent, node, branch.length, label, cluster_id) %>%
#     as_tibble()
#   class(t.tbl) = c("tbl_tree", class(t.tbl))
# 
#   t2 <- as.treedata(t.tbl)
#   text<-message.tree(t2@phylo)
#   strip.nodelabels<-function(text){
#     obj<-strsplit(text,"")[[1]]
#     cp<-grep(")",obj)
#     csc<-c(grep(":",obj),length(obj))
#     exc<-cbind(cp,sapply(cp,function(x,y) y[which(y>x)[1]],y=csc))
#     exc<-exc[(exc[,2]-exc[,1])>1,]
#     inc<-rep(TRUE,length(obj))
#     if(nrow(exc)>0) for(i in 1:nrow(exc)) 
#       inc[(exc[i,1]+1):(exc[i,2]-1)]<-FALSE
#     paste(obj[inc],collapse="")
#   }
#   t2_phylo <- strip.nodelabels(text)
#   t2_phylo <- read.tree(text=t2_phylo)
#   t2@phylo <- t2_phylo
#   return(t2)
# }
# exported_tree <- exportTree(time_tree, cluster_data)
# message.beast(exported_tree, "tree_data.tree")

rt1 <- Sys.time()
rt1-rt0

