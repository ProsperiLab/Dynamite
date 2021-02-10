
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
github_packages <- c("mrc-ide/skygrowth") # "tothuhien/Rlsd2", "mdkarcher/phylodyn" may need to be done if we think our Re values are going to be greater than 5 for any cluster! If the latter, need aslo install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

## Will need to remove install section if using on cluster ###################################
# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
.inst_github <- .packages %in% installed.packages()
## Install GitHub packages(if not already installed)
if(length(github_packages[!.inst_github]) > 0) try(remotes::install_github(github_packages[!.inst_github]))
if(length(github_packages[!.inst_github]) > 0) try(devtools::install_github(github_packages[!.inst_github]))
## Will need to remove install section if using on cluster ###################################

# Load packages into session 
lapply(.packages, require, character.only=TRUE)
lapply(gsub(".+\\/(.+)", "\\1", github_packages), require, character.only=TRUE)

option_list = list(
  make_option(c("-t", "--tree"), type="character", default=list.files(pattern = "*.nwk")[[1]], 
              help="tree file name [default= .nwk extension]", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=list.files(pattern="*.csv")[[1]], 
              help="metadata file name [default= .csv extension]", metavar="character"),
  make_option(c("-s", "--seqLen"), type="numeric", default=30000, 
              help="sequence length [default= 10000]", metavar="numeric"),
  make_option(c("-c", "--cluster"), type="character", default="b", 
              help="choice of cluster algorithm from c (Phylopart's cladewise) or b (DYNAMITE's branchwise) [default= dynamite]", metavar="character"),
  make_option(c("-l", "--leaves"), type="character", default="addLeaves", 
              help="choice of transformation to tree from bifurcating or addLeaves [default=addLeaves]", metavar="character"),
  make_option(c("-r", "--range"), type="numeric", default=30, 
              help="range of branch length thresholds to try [default=30", metavar="numeric"),
  make_option(c("-a", "--asr"), type="character", default="Y", 
              help="option of ancestral state reconstruction for each cluster [default= Y]", metavar="character")
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
write("Checking for tree file...")
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

## Need to force binary tree and to replace zero branch lengths with full bootstrap support


message("Searching for metadata file...")
if(endsWith(opt$metadata, "csv")) {
  metadata <- read.csv(opt$metadata, header=T, stringsAsFactors = F) %>%
    dplyr::filter(ID %in% sub_tree$tip.label)
} else {
  if(endsWith(opt$metadata, "tab") | endsWith(opt$metadata, "txt")) {
    metadata <- read.table(opt$metadata, sep='\t', header=T, stringsAsFactors = F) %>%
      dplyr::filter(ID %in% sub_tree$tip.label)
  } 
}

if(isTRUE("X" %in% colnames(metadata))) { # rownames can be added to metadata if csv file, so need to remove before creating table
    metadata <- dplyr::select(metadata, -"X")
  }    else {metadata <- metadata} 


#Script will first use least squares dating approach developed by To et al. (2016) in the package Rlsd2 (https://github.com/tothuhien/lsd2) to 1) find the optimal root position in the tree and 2) remove outliers with longer-than-expected (penalized) branch lengths, and 3) output both a timed tree and rooted substitutions/site tree to the global environment (updates sub_tree).
# Function to convert sub_tree to time_tree
write("Converting substitution tree into timed tree using lsd2 and dates provided in metadata file. Please make sure the metadata file contains IDs in the first column and dates in the second column. All other columns can be filled as desired...")
checkForDates <- function(metadata) {
  for (i in seq_along(colnames(metadata))) {
    if (tryCatch({
      isTRUE(any(grepl("id", colnames(metadata)[i], ignore.case = T)))}, error = function(e) stderr() )) {
      colnames(metadata)[i] <- "ID"
    } #end first ef statement
    if (tryCatch({
      isTRUE(any(grepl("date", colnames(metadata)[i], ignore.case = T)))}, error = function(e) stderr() )) {
      colnames(metadata)[i] <- "DATE"
    } # End second if statement
  } # End for loop
 
  assign("metadata", metadata, envir = globalenv())
  createSTS <- function(metadata) {
    date_range <- seq(as.Date("01/01/1900", "%d/%m/%Y"), as.Date(Sys.Date(), "%d/%m/%Y"), by="day")
    sts <- metadata$DATE
    if (isTRUE(class(sts) == "integer")) {
      sts <- as.Date(ISOdate(sts))
    } else if (isTRUE(class(sts) == "numeric")) {
      sts <- as.Date(date_decimal(sts))
    } else if (isTRUE(class(sts) == "character")) {
      sts <- as.Date(sts, tryFormats = c("%Y-%m-%d", "%Y/%m/%d"))
      if (isTRUE(min(sts, na.rm=T) < min(date_range))) { # Attempt to correct if in wrong format, which we will know if the first date is before 1900s
        sts <- as.Date(metadata$DATE, tryFormats = c("%m/%d/%Y", "%m-%d-%Y", "%d/%m/%Y", "%d-%m-%Y"))
      }
    } # if-else statement
    metadata$DATE <- sts
    assign("metadata", metadata, envir = globalenv()) # Need dates for metadata table as well for downstream analysis
    sts <- setNames(sts, metadata$ID)
    sts <- sts[names(sts) %in% sub_tree$tip.label]
    return(sts)
  }
  sts <- createSTS(metadata)
   
  ## Checkpoint
#  if (isTRUE(max(as.Date(ISOdate(sts, 1, 1)))> Sys.Date())) {
    if (isTRUE(max(sts) > Sys.Date())) {

    write(paste("The following sequences likely have incorret dates:",
                  return(sts[sts==max(sts)]),
    "Please start DYNAMITE again, placing a .tab metadata file with consistent date information in current folder.",
                  sep=" "))
    stop()
  } # End checkpoint if statement
  ## Will need sts in downstream analyses
  return(sts)
} # End checkForDates function
sts <- checkForDates(metadata)

DateTree <- function(sub_tree, seqLen) {
  # result <- lsd2(inputTree=sub_tree,
  #                inputDate=decimal_date(sts),
  #                estimateRoot="as",
  #                constraint=TRUE,
  #                variance=1,
  #                ZscoreOutlier = 3,
  #                outFile = "lsd2_results",
  #                seqLen = seqLen,
  #                nullblen=-1
  #                )
  # assign("sub_tree", result$newickTree, envir=globalenv())
  # assign("time_tree", result$dateNexusTreeFile@phylo, envir=globalenv())
  # # assign("time_tree_data", result$dateNexusTreeFile, envir=globalenv()) ## Something wrong with node numbers in this tree, so grab from time_tree
  # assign("tmrca", result$tMRCA, envir=globalenv())
  sub_tree <- sub_tree
  result <- dater(sub_tree, decimal_date(sts), s=seqLen, ncpu=numCores, omega0=8E-04)
  assign("sub_tree", result$intree, envir=globalenv())
  assign("time_tree", as.phylo(as_tibble(result)), envir=globalenv())
  # assign("time_tree_data", result$dateNexusTreeFile, envir=globalenv()) ## Something wrong with node numbers in this tree, so grab from time_tree
  assign("tmrca", result$timeOfMRCA, envir=globalenv())
  
}# End DateTree function
DateTree(sub_tree, opt$seqLen)

# Function to specify most recent sampling date (mrsd)
findMRSD <- function(time_tree) {
    date.mrsd <- max(sts[names(sts) %in% time_tree$tip.label])
    mrsd <- decimal_date(date.mrsd)
    
    gg_tree <- ggtree(time_tree, mrsd=date.mrsd) + theme_tree2()
    assign("gg_tree", gg_tree, envir = globalenv())
    node_dates <- gg_tree$data ## Could also do nodeHeights here....
    assign("node_dates", node_dates, envir=globalenv())
    assign("num.mrsd", mrsd, envir=globalenv())
    assign("date.mrsd", date.mrsd, envir=globalenv())
  print(paste0("The updated most recent sampling date is ", date.mrsd))
} # End findMRSD function
findMRSD(time_tree)

# Find all well-supported clades
write("Defining well-supported clades within the tree using node labels...")
define.clades <- function(sub_tree) {
  ## Grab only subtrees that are supported by bootstrap values >90
  ## We may need to change this in case people have other support values
  ## Note that subtrees() function in ape renumbers nodes, so useless here, since at the end we wish to recombine the information
  write("Creating list of all subtrees...This may take a while.")
  family_tree <- tidytree::as_tibble(sub_tree)
  ## Need to relabel columns so that "parent" and "node" are "from" and "to" for familyR::get_children function
  colnames(family_tree)[1:2] <- c("from", "to")
  ## The dataframe needs to be transformed back into a data.table for future analyses
  family_tree <- data.table::as.data.table(family_tree)
  assign("family_tree", family_tree, envir = globalenv())
  
  write("Extract well-supported subtrees.")
  supported_nodes <- list()
  if(max(as.numeric(sub_tree$node.label), na.rm=T) >1) {
  supported_nodes <- suppressWarnings(dplyr::filter(family_tree, as.numeric(label) > 90))
  } else {
    supported_nodes <- suppressWarnings(dplyr::filter(family_tree, as.numeric(label) > 0.90))
  }
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
      if (isTRUE(sum(copy[[cc]]$to %in% clusters[[ct]]$to) > 0.1*length(copy[[cc]]$to)) &
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
branchLengthLimit <- function(tree, range) {
  p.dist.mat.leaves <- cophenetic(tree)
  ## Alternative distance matrix
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
    node_num <- (ntips+2):(ntips+nint) # Return all internal nodes
    if (isTRUE(length(node_num) >= 1E06)) { # If tree too large, sample 10,000 nodes for MPPD calculation
      node_num <- sample(node_num, 10000, replace = F)
    }
    
    if(mode=='leaf'){
      distmat <-  p.dist.mat.leaves
      MPPD <- sapply(node_num,get.node.leaf.MPPD,tree,distmat) # For each node, extract the MPPD (across patristic distances) for the corresponding subtree
      return(data.frame(node_num=node_num, MPPD=MPPD))
    }
    else{
      distmat <-  dist.nodes(tree)
      MPPD <- sapply(node_num,get.node.full.MPPD,tree,distmat) # For each node, extract the MPPD (across all branch lengths) for the corresponding subtree
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
  
  ## Create a vector of MPPDs for plotting and determining branch length limit
  distvec <- pdist.clusttree(tree, mode='leaf')
  hist(distvec$MPPD)
  
  ## Determine MPPDs for all well-supported clades
  clade_MPPD <- pdist.clades(clades, tree, mode='leaf')
  assign("clade_MPPD", clade_MPPD, envir=globalenv())
  
  optimizeThreshold <- function(distvec, range) {
    # thresholds <- 0.01
    # bl <- quantile(distvec$MPPD, thresholds)
    # 
    # while(length(unique(do.call("rbind", branch_length_limit))) == 1){
    n <- range
    y0 <- rexp(n, 10)  # simulate exp. dist.
    bl <- mclapply(y0, function(x) {
      quantile(distvec$MPPD, x)
    }, mc.cores=numCores)
    
    return(bl)
  } # End function
  
  #  branch_length_limit <- quantile(distvec$MPPD, opt$threshold)
  branch_length_limit <- optimizeThreshold(distvec, range)
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
      if (make_tree=="addLeaves") {
        clusters <- addLeaves(tree, clusters)
        clusters <- mclapply(clusters, as_tibble, mc.cores=numCores)
    }
  }
  
  clusters <-  list.filter(clusters, sum(label %in% tree$tip.label) >= 3)
  
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
    #    merge.nested.clust() %>%
    # merge.overlap.clust() %>% # Will not work well with support in labels
    # merge.nested.clust() %>%
    list.filter(length(label) >= 3)
  return(clusters)
}


message("Picking clusters based on user choice of algorithm....")
branch_length_limit <- branchLengthLimit(sub_tree, opt$range)
if (opt$cluster == "b") {
  possible_clusters <- mclapply(branch_length_limit, function(x) branchWise(sub_tree, x, make_tree=opt$leaves), mc.cores = numCores)
  nclust <- mclapply(possible_clusters, length, mc.cores=numCores)
  best <- max(do.call("rbind", nclust))
  clusters <- possible_clusters[which(nclust==best)]
  branch_length_limit <- branch_length_limit[which(nclust==best)]
  
  if (length(clusters)>1) {
    clusters <- last(clusters)[[1]]
  } else {
    clusters <- clusters[[1]]
  }
} else {
  if (opt$cluster == "c") {
    possible_clusters <- mclapply(branch_length_limit, function(x) phylopart(sub_tree, x), mc.cores = numCores)
    nclust <- mclapply(possible_clusters, length, mc.cores=numCores)
    best <- max(do.call("rbind", nclust))
    clusters <- possible_clusters[which(nclust==best)]
    if (length(clusters)>1) {
      clusters <- last(clusters)[[1]]
    } else {
      clusters <- clusters[[1]]
    }
    clusters <- mclapply(clusters, function(x) dplyr::rename(x, parent=from, node=to), mc.cores=numCores)
    branch_length_limit <- branch_length_limit[which(nclust==best)]
    
  } else {
    message("Incorrect cluster_picking algorithm choice. Please choose between 'b' (branch-wise) or 'c' (clade-wise) and run script again.")
  }
}


clusterPhylo <- function(clusters, sub_tree, time_tree) {
  clusters_time_tree <- mclapply(clusters, function(x) {
    tips <- x$label[x$label %in% time_tree$tip.label]
    full_clade <- extract.clade(time_tree, findMRCA(time_tree, tips))
    cluster <- keep.tip(full_clade, tips) # Need to make sure not full clade, but only tips present in cluster
    return(cluster)
    }, mc.cores=numCores)
  assign("clusters_time_tree", clusters_time_tree, envir = globalenv() )
  
  clusters_sub_tree <- mclapply(clusters, function(x) {
    tips <- x$label[x$label %in% sub_tree$tip.label]
    full_clade <- extract.clade(sub_tree, findMRCA(sub_tree, tips))
    cluster <- keep.tip(full_clade, tips) # Need to make sure not full clade, but only tips present in cluster
    return(cluster)
  }, mc.cores=numCores)

#  clusters_sub_tree <- mclapply(clusters, as.phylo, mc.cores=numCores)
  assign("clusters_sub_tree", clusters_sub_tree, envir = globalenv() )
}
clusterPhylo(clusters, sub_tree, time_tree)

#Employ ancestral state reconstruction of traits (optional)

ancStateRecon <- function(tree, metadata) {
## The element lik.anc gives us the marginal ancestral states, also known as the 'empirical Bayesian posterior probabilities.'
write("Estimating likelihood of ancestral states...", stderr())
feed.mode <- list()
fitER <- list()
## set zero-length branches to be 1/1000000 total tree length
dst<-multi2di(tree)
dst$edge.length[dst$edge.length==0]<-max(nodeHeights(tree))*1e-6

rownames(metadata) <- metadata$ID
for (i in 3:ncol(metadata)) { # Requires that taxa id and sampling dates are first two columns
  feed.mode[[i-2]] <- setNames(metadata[,i],rownames(metadata))
  fitER[[i-2]] <- ace(feed.mode[[i-2]], dst, model = "ER", type="discrete")
}
names(fitER) <- colnames(metadata[,3:ncol(metadata)])

## Look here for http://www.phytools.org/Cordoba2017/ex/8/Anc-states-discrete.html for distribution of sampled stochastic character maps when have more computing power

## simulate single stochastic character map using empirical Bayes method
#mtrees<-make.simmap(sub_tree,feed.mode,model="ER",nsim=100)

### Create rate matrix (Q) ################################################
### This will need to be modified to incorporate several columns of traits
### See https://cran.r-project.org/web/packages/filling/filling.pdf for incomplete data ### filling!
# Q <-list()
# for (column in 2:ncol(traits)) {
#   suppressWarnings({ ## We know it fills diagonals with NAs
#   Q[[column]] <- diag(unique(traits[,column]), nrow = length(unique(traits[,column])))
#   })
#   diag(Q[[column]]) = 1-nrow(Q[[column]])
#   Q[[column]][lower.tri(Q[[column]])] <- 1
#   Q[[column]][upper.tri(Q[[column]])] <- 1
#   colnames(Q[[column]]) <- rownames(Q[[column]]) <- unique(traits[,column])
# }
# Q <- plyr::compact(Q)
# names(Q) <- colnames(traits[-1])


anc.states <- fitER
for (i in 1:length(fitER)) {
  anc.states[[i]] <- as.data.frame(fitER[[i]]$lik.anc)
  anc.states[[i]]$node <- as.character(1:dst$Nnode+Ntip(dst))
}
 names(anc.states) <- names(fitER)
 

## Combine node state labels and probabilities with tip state labels (and prob of 1) - may need to go back and use smart bind so that have separate columns for tip.label?

## Should consider transforming data in to states, rather than probabilities, and assign state of "unknown" to nodes with <90% probability of any particular state

tip_labels <- list()
for (i in 3:ncol(metadata)) {
  tip_labels[[i-2]] <- as.data.frame(cbind(rownames(metadata), metadata[,i]), stringsAsFactors = F)
  colnames(tip_labels[[i-2]])[1] <- "node"
  tip_labels[[i-2]] <- dplyr::mutate(tip_labels[[i-2]], value=1, V2) %>%
    spread(V2, value, fill=0)
  ## Reorder tip labels according to order in sub_tree and assign numeric labels given in sub_tree in order to integrate all data into a tibble.
  tip_labels[[i-2]] <- tip_labels[[i-2]][order(match(tip_labels[[i-2]]$node, dst$tip.label)),]
  tip_labels[[i-2]]$node <- as.character(as.numeric(1:length(dst$tip.label)))
}
names(tip_labels) <- colnames(metadata[,3:ncol(metadata)])

for (i in seq_along(anc.states)) {
  anc.states[[i]] <- rbind(tip_labels[[i]], anc.states[[i]])
}

## For now, best to only assign one state (per trait) per node, rather than probabilities, and, in order to incorporate uncertainty, we will assign "NEI" (Not Enough Info) to any node with <90% posterior probability of any given state.

anc.states <- mclapply(anc.states, function(x) {
  gather(x, trait, pr, -node) %>%
    group_by(node) %>%
    dplyr::mutate(max_pr = max(pr)) %>% ## Determine if reliable state found for each node using max probability
    dplyr::mutate(trait = case_when(max_pr < 0.90 ~ "NEI",
                                      max_pr >= 0.90 ~ trait)) %>% ## When a node does not have single, reliable state, assign "NE" for "Not enough info" as the trait to that node
    dplyr::filter(pr>0.90 | trait == "NEI") %>% ## Assign each node a single trait (the one with >90% prob)
    dplyr::select(-pr, -max_pr) %>% ## Get rid of probability columns
    dplyr::arrange(as.numeric(node)) %>% # Order the data by node (not necessary but useful)
    dplyr::distinct() %>%# Left with redundant info for NEI nodes, so remove.
    dplyr::mutate(node = as.integer(node)) %>%
    left_join(., dplyr::select(node_dates, node, x), by="node") %>%
    dplyr::rename(to = "node", DATE = "x")
  #names(anc.states[[i]])[names(anc.states[[i]]) == 'node'] <- 'to'
      ## Need to change node back to integer to integrate with data downstream
#    anc.states[[i]]$to <- as.integer(anc.states[[i]]$to)
#    anc.states[[i]] <- anc.states[[i]][order(anc.states[[i]]$to),]
}, mc.cores=numCores) # End loop among ancestral states

return(anc.states)
} # End asrClusters function on traits

write("Performinc ancestral state reconstruction...")
if (opt$asr == "Y"){
write("ASR performed using the empirical Bayesian method.")
  asr <- ancStateRecon(sub_tree, metadata) #traits
} else {
  if (opt$asr == "N"){
  write("ASR not performed - results will rely on sampled data only.")
} else {
  write("Correct response not detected.")
}
}
  

# Data manipulation for characterization of clusters and export of data for visualization
write("Performing quick data manipulation... Hold on to your butts!")
## Merge metadata with cluster data
dataManip <- function(clusters) {

## Need to convert clusters to list if only a single cluster found (will be in dataframe format, rather than list)

if (class(clusters) == "data.frame") {
  clusters <- list(clusters) 
  } else if (class(clusters) == "list") {
    clusters <- clusters
  } else {warning("class of clusters data unknown, error in dataManip() function or above")}
  
  metadata2 <- gather(metadata, field, trait,  -ID, -DATE, factor_key=TRUE)
  cluster_data <- list()

  if (isTRUE(exists("asr", envir = globalenv()))) { #Only use if ASR incorporated above
    asr2 <- bind_rows(asr, .id = "field")
    asr2$to <- as.integer(asr2$to)

#   cluster_data <- mclapply(clusters, function(c) {
#     merge(c, asr2, by.x="node" , by.y="to") %>%
# #      dplyr::rename(parent = "from", node = "to") %>%
#       dplyr::select(parent, everything()) %>%
#       merge(., dplyr::select(node_dates, node, x), by.x=c("node", "DATE"), by.y=c("node", "x"))
#   }, mc.cores=numCores)
#   
#   } else {
    
    
    anc.state.data <- mclapply(clusters, function(x) {
      anc.field <- dplyr::select(asr2, field, to) %>%
        dplyr::filter(to == x$parent[1])
      anc.trait <- dplyr::select(asr2,trait, to) %>%
        dplyr::filter(to == x$parent[1])
      x <- merge(anc.field, anc.trait, by="to")
      return(x)
    }, mc.cores = numCores)
    
    
   
    for (i in seq_along(clusters)) {
      cluster_data[[i]] <- merge(clusters[[i]], metadata2, by.x = "label", by.y="ID", all=T) %>%
        dplyr::select(parent, everything())
      cluster_data[[i]] <- merge(cluster_data[[i]], anc.state.data[[i]], by.x=c("parent", "field", "trait"), by.y=c("to", "field", "trait"), all=T)
    } 
    } else {
      
      for (i in seq_along(clusters)) {
        cluster_data[[i]] <- merge(clusters[[i]], metadata2, by.x = "label", by.y="ID", all.x=T) %>%
          dplyr::select(parent, node, field, trait, DATE)
      }
    

}
names(cluster_data) <- names(clusters)
return(cluster_data)

} # End dataManip() function
cluster_data <- dataManip(clusters)

## Now DYNAMITE determines if clusters are related by connecting the children of each cluster to the root of remaining clusters
connectClust <- function(sub_tree, cluster_data) {
  
  cluster_data <- cluster_data
  dup_cluster_data <- cluster_data
  full_tree <- as_tibble(sub_tree)
  
  for (i in seq_along(cluster_data)) {
    
    for (j in seq_along(dup_cluster_data)) {
      dup_cluster_data[[j]]$birth_origin <- NA
      #      dup_cluster_data[[j]]$birth_origin <- NA
      ## If the parent of origin node of one cluster (found by referencing full tree) is a child node in one of the other clusters... 
      if (isTRUE(names(cluster_data)[[i]] != names(dup_cluster_data)[[j]] &
                 (tidytree::parent(full_tree, dup_cluster_data[[j]]$parent[1])$node %in% cluster_data[[i]]$node |
                 ## Or if the parent of the parent of the origin node of one cluster (found by referencing full tree) is a child node in one of the other clusters... (meaning clusters are allowed to be separated by up to two branches)
                 tidytree::parent(full_tree, dup_cluster_data[[j]]$parent[1])$parent %in% cluster_data[[i]]$node))) {
        #        dup_cluster_data[[j]]$birth_origin == paste0("cluster_", cluster_data[[i]]$parent[1] )
        dup_cluster_data[[j]]$birth_origin = paste0(names(cluster_data)[[i]])
      } #End if statement
    } # End for loop along duplicate list
  } # End for loop along original list
  return(dup_cluster_data)
} #End connectClust function
cluster_data <- connectClust(sub_tree, cluster_data) 

# Added tree statistics, but currently using machine learning to determine which of the following
# are actually better predictors of cluster dynamics
calculateNe <- function(tree) {
  tree <- multi2di(tree)
  res=round(max(nodeHeights(tree))/2) # Daily
  fit <- tryCatch(skygrowth.map(tree, res=res, tau0=0.1), error=function(e) NULL)
  p <- data.frame(time=max(nodeHeights(tree))-fit$time, nemed=fit$ne, ne_ci=fit$ne_ci)
  return(p)
}
calculateR0 <- function(Ne, conf.level=0.95) {
  s <- 100 # sample size of 100
  #  psi <- rnorm(s, 0.038, 0.014) #Duration of infection around 14 days
  psi <- rnorm(s, 14, 5) #Duration of infection around 14 days -incubation period of 5 days
  #  psi <- rnorm(s, 9, 5) #Duration of infection around 14 days -incubation period of 5 days
  Z=qnorm(0.5*(1 + conf.level))
  Re <- list()
  if (isTRUE(nrow(Ne) >1)) {
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
  } else {Re <- list(NULL)}
  Re <- do.call("rbind",Re)
  R0 <- Re$mean_Re[1]
  return(R0)
}
calculateLTT <- function(tree) {
  fit <- ltt(tree, plot=F)
  df <- distinct(data.frame(ne=fit$ltt, time=fit$times))
  df <- df[-1,]
  df <- data.frame(x=df$time, y=df$ne)
  
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
  
  
  cc <- check_curve(df$x, log(df$y))$ctype
  # if (cc=="convex") {
  #   model="constant"
  # } else {
  #   if (cc=="concave") {
  #     model="growth"
  #   } else { model="mixed"}
  # }
  result <- data.frame()
  return(cc)
}
calculateGrowth <- function(cluster_Ne, total_r) {# modelNe_log <- function(Ne_df) {
  cluster_df <- data.frame(x=cluster_Ne$time, y=cluster_Ne$nemed)
  if (nrow(cluster_df) > 1) {
    max_cluster_ne <- max(cluster_df$y, na.rm=T); min_cluster_ne <- min(cluster_df$y, na.rm=T)
    max_cluster_t <- max(cluster_df$x, na.rm=T); min_cluster_t <- min(cluster_df$x, na.rm=T)
    cluster_r = (max_cluster_ne-min_cluster_ne)/(max_cluster_t-min_cluster_t)
    
    relative_r <- cluster_r/total_r
    
    rmax <- list()
    for (i in 2:nrow(cluster_df)) {
      rmax[[i]] <- (cluster_df$y[i]-cluster_df$y[i-1])/(cluster_df$x[i]-cluster_df$x[i-1])
    }
    rmax <- max(do.call("rbind",rmax), na.rm=T)
    
    #Fraction of time spent in growth
    max_cluster_t.ne <- cluster_df$x[cluster_df$y==max_cluster_ne]
    t_frac <- (max_cluster_t.ne-min_cluster_t)/(max_cluster_t-min_cluster_t)
  } else {
    cluster_r <- NA
    relative_r <- NA
    t_frac <- NA
    rmax <- NA
  }
  
  result <- data.frame(abs_growth_rate=cluster_r, rel_growth_rate=relative_r, t_frac=t_frac, rmax=rmax)
  return(result)
}
calculateGamma <- function(tree) {
  tree <- multi2di(tree)
  gamma <- gammaStat(tree)
  return(gamma)
}
calculateOster <- function(tree) {
  tree <- multi2di(tree)  
  cluster_size <- length(tree$tip.label)-1
  sum_heights <- sum(nodeHeights(tree))
  longest <- max(nodeHeights(tree))
  co <- cluster_size/sum_heights + longest
  return(co)
}
timeData <- function(tree, clusters) {
  sts <- distRoot(tree@phylo, tips = "all", method="patristic")
  sts <- data.frame(time = sts, ID=names(sts))
  assign("sts", setNames(sts$time, sts$ID), envir = globalenv())
  
  tbl_list <- mclapply(clusters, function(x) {
    mrsd <- max(sts$time[sts$ID %in% x$tip.label], na.rm=T)
    cluster_size <- length(x$tip.label)
    root_age <- mrsd-max(nodeHeights(x))
    x <- cbind(as_tibble(x), mrsd = mrsd,
               cluster_size = cluster_size,
               root_age = root_age 
    )}, mc.cores=numCores) # End creation of tbl_list
  
  assign("tbl_list", tbl_list, envir = globalenv())
  # clust_metadata <- rbindlist(tbl_list, fill=T) %>%
  #   dplyr::select(-parent, -node, -branch.length) # No longer need parenta and node information, since not unique
  # names(clust_metadata)[1] <- "ID"
  
  # sts_list <- mclapply(clusters, function(c) {
  #   c <- sts[sts$ID %in% c$tip.label,]
  #   c <- setNames(c$time, c$ID)
  # }, mc.cores=numCores)
  # 
  # return(sts_list)
}
calculateNeeBD <- function(tree) {
  tree <- multi2di(tree)
  #tree = prune.extinct.taxa(tree)
  fit.bd <- tryCatch(birthdeath(tree), error=function(e) NULL)
  if (!is.null(fit.bd)) {
    b <- bd(fit.bd)[1]
  } else {
    b <- NULL
  }
  return(b)
}
calculatePD <- function(tree) {
  tree <- multi2di(tree)
  pd <- sum(tree$edge.length)
  return(pd)
}
calculateCC <- function(tree) {
  tree <- multi2di(tree)
  out <- capture.output(
    tryCatch(cherry(tree), error=function(e) NULL)
  )
  ntips <- as.numeric(gsub("Number of tips\\: ([0-9]+) ", "\\1", out[5]))
  ncherries <- as.numeric(gsub("Number of cherries\\: ([0-9]+) ", "\\1", out[6]))
  cpn <- ncherries/ntips
  return(cpn)
}



gatherStats <- function(time_tree, clusters_time_tree, clusters_sub_tree) {
  total_Ne <- calculateNe(time_tree)
#  total_gamma <- calculateGamma(time_tree)
#  total_NeeBD <- calculateNeeBD(time_tree)
#  total_oster <- calculateOster(sub_tree)
#  total_PD <- calculatePD(sub_tree)
  
  total_df <- data.frame(x=total_Ne$time, y=total_Ne$nemed)
  max_total_ne <- max(total_df$y); min_total_ne <- min(total_df$y)
  max_total_t <- total_df$x[total_df$y==max_total_ne]; min_total_t <- min(total_df$x)
  total_r = (max_total_ne-min_total_ne)/(max_total_t-min_total_t)
  
  cluster_Ne <- mclapply(clusters_time_tree, calculateNe, mc.cores=numCores)
  cluster_R0 <- mclapply(cluster_Ne, calculateR0, conf.level=0.95, mc.cores=numCores)
  cluster_gamma <- mclapply(clusters_time_tree, calculateGamma, mc.cores=numCores)
  cluster_growth_rates <- mclapply(cluster_Ne, function(x) {
    calculateGrowth(x, total_r)}, mc.cores=numCores)
  cluster_ltt_shape <- mclapply(clusters_time_tree, calculateLTT, mc.cores=numCores)
  cluster_NeeBD <- mclapply(clusters_time_tree, calculateNeeBD, mc.cores=numCores)
  cluster_cherries <- mclapply(clusters_sub_tree, calculateCC, mc.cores=numCores)
  cluster_oster <- mclapply(clusters_sub_tree, calculateOster, mc.cores=numCores)
  cluster_PD <- mclapply(clusters_sub_tree, calculatePD, mc.cores=numCores)
  
  ## Populate cluster_dynamics table
  cluster_dynamics <- rep(list(data.frame(gamma = as.numeric(NA),
                                 R0 = as.numeric(NA),
                                 birth_rate = as.numeric(NA),
                                 #rel_birth_rate = as.numeric(NA),
                                 abs_growth_rate = as.numeric(NA),
                                 #rel_growth_rate = as.numeric(NA),
                                 r_max = as.numeric(NA),
                                 fraction_time_growth = as.numeric(NA),
                                 oster = as.numeric(NA),
                                 PD = as.numeric(NA),
                                 #rel_PD = as.numeric(NA),
                                 ltt_shape = as.character(NA),
                                 cherries = as.character(NA))),
                          length(clusters))
  
  for (i in seq_along(cluster_dynamics)) {
    cluster_dynamics[[i]]$gamma <- cluster_gamma[[i]]
    cluster_dynamics[[i]]$oster <- cluster_oster[[i]]
    cluster_dynamics[[i]]$birth_rate <- cluster_NeeBD[[i]]
    #cluster_dynamics[[i]]$rel_birth_rate <- cluster_NeeBD[[i]]/total_NeeBD
    cluster_dynamics[[i]]$PD <- cluster_PD[[i]]
    #cluster_dynamics[[i]]$rel_PD <- cluster_PD[[i]]/total_PD
    cluster_dynamics[[i]]$R0 <- cluster_R0[[i]]
    cluster_dynamics[[i]]$abs_growth_rate <- cluster_growth_rates[[i]]$abs_growth_rate
    #cluster_dynamics[[i]]$rel_growth_rate <- cluster_growth_rates[[i]]$rel_growth_rate
    cluster_dynamics[[i]]$fraction_time_growth <- cluster_growth_rates[[i]]$t_frac
    cluster_dynamics[[i]]$r_max <- cluster_growth_rates[[i]]$rmax
    cluster_dynamics[[i]]$ltt_shape <- cluster_ltt_shape[[i]]
    cluster_dynamics[[i]]$cherries <- as.numeric(cluster_cherries[[i]])
  }
  
  names(cluster_dynamics) <- names(clusters)
  cluster_dynamics <- dplyr::bind_rows(cluster_dynamics, .id = "cluster_id")
  
  return(cluster_dynamics)

}

cluster_tree_stats <- gatherStats(time_tree, clusters_time_tree, clusters_sub_tree)

# Need to merge tree_stats with ASR data

cluster_data <- dplyr::bind_rows(cluster_data, .id = "cluster_id") %>%
  dplyr::filter(!is.na(DATE))

cluster_info <- list(trait_distributions=cluster_data,
                     tree_stats=cluster_tree_stats)

# Transform sampled tree into a tbl object and assign cluster IDs to internal and external nodes found in list of clusters

clusters.tbl <- dplyr::bind_rows(clusters, .id="cluster_id")
tree.tbl <- as_tibble(time_tree)
tree.tbl$cluster_id <- "Background"
for (i in 1:nrow(clusters.tbl)) {
  for (j in 1:nrow(tree.tbl)) {
    if (isTRUE(tree.tbl$parent[j] == clusters.tbl$parent[i])) {
      tree.tbl$cluster_id[j] = clusters.tbl$cluster_id[i]
    } # End if statement
  } # end loop along tree
} # End loop along cluster_data


#sampled.tree$label <- paste(sampled.tree$label, sampled.tree$state, sampled.tree$date, sep="|")
# sampled.tree <- dplyr::select(sampled.tree, parent, node, branch.length, label, cluster_id) %>%
#   as_tibble()
class(tree.tbl) = c("tbl_tree", class(tree.tbl))
t2 <- as.treedata(tree.tbl)

write("Data are now being exported as 'cluster_info_<tree>.RDS' and 'dynamite_<tree>.tree.'")
saveRDS(cluster_info, paste0("cluster_info_", opt$tree, ".RDS"))
write.table(branch_length_limit, paste0("branch_length_limit_", opt$tree, ".txt"))
write.beast(t2, paste0("dynamite_", opt$tree, ".tree")) #### NEED THIS OUTPUT####################################

## Combined results


rt1 <- Sys.time()
rt1-rt0

