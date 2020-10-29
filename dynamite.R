
##Install and load necessary packages
## Initialize stable working environment and store time of intiation
rm(list=ls())
setwd(getSrcDirectory()[1]) #When working on the cluster
#dirname(rstudioapi::getActiveDocumentContext()$path) # If working in Rstudio
rt0 <- Sys.time()

# List of packages for session
.packages <-  c("optparse", "remotes", "phytools", "data.tree", 
                "tidytree", "rlist", "familyR", "tidyverse", 
                "ggtree", "parallel", "geiger", "tibble")  # May need to incorporate code for familyR (https://rdrr.io/github/emillykkejensen/familyR/src/R/get_children.R) i fno longer supported.
github_packages <- c("mrc-ide/skygrowth", "tothuhien/Rlsd2") # mrc-ide/skygrowth, "mdkarcher/phylodyn" may need to be done if we think our Re values are going to be greater than 5 for any cluster! If the latter, need aslo install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

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
  make_option(c("-m", "--metadata"), type="character", default=list.files(pattern="*.tab")[[1]], 
              help="metadata file name [default= .tab extension]", metavar="character"),
  make_option(c("-s", "--seqLen"), type="numeric", default=30000, 
              help="sequence length [default= 30000]", metavar="numeric"),
  make_option(c("-a", "--cluster"), type="character", default=list.files(pattern="*results.csv")[[1]], 
              help="choice of cluster algorithm from c (Phylopart's cladewise) or b (DYNAMITE's branchwise) [default= phylopart]", metavar="character"),
  make_option(c("-asr", "--asr"), type="character", default="N", 
              help="option of ancestral state reconstruction for each cluster [default= N]", metavar="character"),
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Additional options
`%notin%` <- Negate(`%in%`) # Just plain useful
# The option below is useful when dealing with dates of internal nodes downstream
options(digits=15)


# Supply tree and metadata table as arguments
threshold=as.numeric(0.10)

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
if (isFALSE(is.binary(sub_tree))) {
  sub_tree <- multi2di(sub_tree)
  
  sub_tree$node.label[sub_tree$node.label==""] <- "100"
} else {
  sub_tree <- sub_tree
}



message("Searching for metadata file...")
metadata <- read.table(opt$metadata, sep = '\t', header=T, stringsAsFactors = F) %>%
  dplyr::filter(ID %in% tree$tip.label)

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
  createSTS <- function(metadata_mod) {
    date_range <- seq(as.Date("12/01/2019", "%d/%m/%Y"), as.Date(Sys.Date(), "%d/%m/%Y"), by="day")
    sts <- metadata_mod$DATE
    if (class(sts) == "integer") {
      sts <- as.Date(ISOdate(sts))
    } else if (class(sts) == "numeric") {
      sts <- as.Date(date_decimal(sts))
    } else if (class(sts) == "character") {
      sts <- as.Date(sts, tryFormats = c("%Y-%m-%d", "%Y/%m/%d", "%m/%d/%Y", "%m-%d-%Y", "%d/%m/%Y", "%d-%m-%Y"))
      if (min(sts) < min(date_range)) { # Attempt to correct if in wrong format, which we will know if the first date is before Dec 2019
        sts <- as.Date(sts, "%d/%m/%Y")
      }
    } # if-else statement
    sts <- setNames(sts, metadata$ID)
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
sts <- checkForDates(sub_tree)

DateTree <- function(sub_tree, seqLen) {
#   result <- dater(sub_tree, decimal_date(sts), ncpu=4, parallel_foreach=TRUE)
  result <- lsd2(inputTree=sub_tree,
                 inputDate=decimal_date(sts),
                 estimateRoot="as",
                 constraint=TRUE,
                 variance=1,
                 ZscoreOutlier = 3,
                 outFile = "lsd2_results",
                 seqLen = seqLen,
                 nullblen=-1
                 )
  assign("sub_tree", result$newickTree, envir=globalenv())
  assign("time_tree", result$dateNexusTreeFile@phylo, envir=globalenv())
  # assign("time_tree_data", result$dateNexusTreeFile, envir=globalenv()) ## Something wrong with node numbers in this tree, so grab from time_tree
  assign("tmrca", result$tMRCA, envir=globalenv())
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
  supported_nodes <- suppressWarnings(filter(family_tree, as.numeric(label) > 90))
  } else {
    supported_nodes <- suppressWarnings(filter(family_tree, as.numeric(label) > 0.90))
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


# We can now use a coalescent model, which considers the genealogical process of our sample of taxa taken from an assumed large population that changes in time deterministically. The population size is assumed to be homogeneous and under neutral evolution; although in practice these assumptions are violated, the ‘effective population size’, Ne, can often still be derived, which gives the same coalescence rate as an idealized population of size N. During exponential growth, there is a linear relationship between the prevalence and the incidence, and hence the coalescence rate is directly proportional to the number of infected individuals (Frost and Volz, 2010). Assuming our sampled sequences exhibit the behavior described above, we will use a non-parametric coalescent model to estimate Ne and a polynomial model fit to the Ne to determine the time window during which the exponential growth phase exists. 

# Since genealogies representative of exponentially increasing populations often provide very little information about effective population size near the present (or most recent sample) (de Silva et al. 2012), original "skyline" (Pybus et al., 2000; Strimmer and Pybus, 2001) estimators with Brownian motion priors (Minin et al. 2008; Palacios and Minin 2013) on Ne may produce estimates which stabilize at a constant level even when the true size is increasing or decreasing exponentially. Volz and Didelot's "skygrowth" model offers a more realistic prior that is defined in terms of the growth rate of Ne. 

# Following Ne estimation using skygrowth and definition of the exponential growth period, the median branch length (scaled in substitutions/site) during this period of time can be used as the cutoff to determine the incorporation of internal and external nodes into a "transmission cluster" throughout the remainder of the tree. The reasoning behind this is that, during this time, the internal nodes have a higher probability of being considered individuals involved in direct transmission than the remainder of the tree and enable us to define the median branch length as the median genetic distance separating direct transmission events. 

# Using a pre-order (root-to-tip) tree traversal, for each node, find subtrees/clades (>=3 nodes) for which cumulative mean branch, or edge, length is within the branch length limit at that level within the tree using the "branch-wise" algorithm, or "b". Note this should only be used with epidemics that are past the exponential phase of growth.

# Alternatively, using a depth-first algorithm developed by Prosperi et al., (2011), find subtrees/clades (>=2 nodes) for which median pairwise patristic distances (MPPD) between nodes is within 1-5% of the distribution of MPPDs of clades within the entire tree. This is referred to as the "clade-wise algorithm", or "c".

## Functions used in DYNAMITE cluster-picking algorithm ############################################################################
branchWise <- function(tree) {
  findThreshold <- function(tree) {
  sub_tree <- as.phylo(tree)
  b0 <- BNPR(multi2di(tree))
  b0$data$E_log[b0$data$E_log<0] = 0
  plot_BNPR( b0 )
  p <- data.frame(time=rev(b0$data$time), Elog=b0$data$E_log)
  mL <- try(drm(data=p, Elog~time, fct = LL.3(), type = "continuous"), silent=T)
  p2 <- plot(mL)
  
  dY <- diff(p2$`1`)/diff(p2$time)  # the derivative of your function
  dX <- rowMeans(embed(p2$time,2)) # centers the X values for plotting
  plot(dX,dY,type="l",main="Derivative") #check
  d1 <- data.frame(x=dX, y=dY)
  d1$x[which.max(d1$y)] # This is slightly less than the inflection point, which is 50% of time
  
  time_epi_peak <- max(nodeHeights(tree, root.edge=TRUE)) - d1$x[which.max(d1$y)]
  
  # fit <- fit_easylinear(p2$`rev(b0$x)`, p2$`1`, quota=0.99) #Foound this to be the best, except when prob_exit and p_trans are both high
  # plot(fit)
  # time_epi_peak <- max(nodeHeights(sub_tree, root.edge=TRUE)) - min(fit@fit$model$x)
  
  time.tree.exp <- paleotree::timeSliceTree(sub_tree, time_epi_peak)
  #  bl_rescaled <- time.tree.exp$edge.length*(time_tree_data$mean.rate)
  bl_rescaled <- time.tree.exp$edge.length
  plot(density(bl_rescaled))
  
  branch_length_limit <- median(bl_rescaled) ## Ideally will need to multiply branch lengths by estimated rate.
  #branch_length_limit <- mean(bl_rescaled) ## Ideally will need to multiply branch lengths by estimated rate. # Sometimes median works best
  ## Maybe we should take both and whichever is larger?
  #branch_length_limit <- mean(time.tree.exp$edge.length) + qnorm(.95)*(sd(time.tree.exp$edge.length)/sqrt(Ntip(time.tree.exp)))
  
  return(branch_length_limit)
} # End findThreshold()
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
    results[[c]] <- extract.clade(sub_tree, clusters[[c]]$from[1])
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

  branch_length_limit <- findThreshold(sub_tree)

  clusters <- mclapply(clades, pickClust, mc.cores=numCores) %>%
    compact() %>%
    merge.nested.clust() %>%
    merge.overlap.clust() %>%
    merge.nested.clust() %>% ## Why am I having to run this again????
    list.filter(length(label) >= 5)
    #list.filter(sum(label %in% sub_tree$tip.label) >= 5)
    ## Remove singleton nodes (non-bifurcating branches) by using bifurcate() function or extracting entire clade.
    #clusters <- bifurcate(sub_tree, clusters)
    #clusters <- pullClade(sub_tree, clusters)
  clusters <- addLeaves(tree, clusters)
  
  return(clusters)
}
phylopart <- function(tree) {
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
    mclapply(., function(x){
      dplyr::select(x, from, to, branch.length, label) %>%
        dplyr::arrange(from,to) 
    }, mc.cores=numCores) %>%
    unique()
  
  return(result)
}  # In case you want to remove this and consider only fully nested clusters
  ### Create matrix of each pairwise patristic distance for external leaves using the following
  leaves <- sample(tree$tip.label, 0.50*length(tree$tip.label))
  leaves <- expand.grid(leaves,leaves)
  p.dist.leaves <- sapply(seq_len(nrow(leaves)), ## Create list of all pairwise combinations of IDs using expand.grid()
                        function(k) { #future_sapply actually slower here!
                          i <- leaves[k,1]
                          j <- leaves[k,2]
                          fastDist(tree, i,j)
                        })
  p.dist.mat.leaves <- matrix(p.dist.leaves,
                            nrow=Ntip(tree), ncol=Ntip(tree),
                            dimnames=list(tree$tip.label,tree$tip.label))


  ## Create a vector of MPPDs for plotting and determining branch length limit
  distvec <- pdist.clusttree(tree, mode='all')
  hist(distvec$MPPD)

  ## Determine MPPDs for all well-supported clades
  clade_MPPD <- pdist.clades(clades, tree, mode='all')

  phylopart.threshold <- threshold
  branch_length_limit <- quantile(distvec$MPPD, phylopart.threshold)

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

write("Picking clusters based on user choice of algorithm....")
if (opt$cluster == "b") {
  clusters <- branchWise(sub_tree)
} else {
  if (opt$cluster == "c") {
    clusters <- phylopart(sub_tree)
  } else {
    write("Incorrect cluster_picking algorithm choice. Please choose between 'b' (branch-wise) or 'c' (clade-wise) and run script again.")
  }
}


#Employ ancestral state reconstruction of traits (optional)

ancStateRecon <- function(tree, metadata) {
## The element lik.anc gives us the marginal ancestral states, also known as the 'empirical Bayesian posterior probabilities.'
write("Estimating likelihood of ancestral states...", stderr())
feed.mode <- list()
fitER <- list()
## set zero-length branches to be 1/1000000 total tree length
dst<-tree
dst$edge.length[dst$edge.length==0]<-max(nodeHeights(tree))*1e-6

rownames(metadata) <- metadata[,1]
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


## Now need to figure out if node numbers match in lik.anc with those of sub_tree
# library(RColorBrewer)
# #display.brewer.all(colorblindFriendly = TRUE)
# cols<-brewer.pal(length(unique(traits[,1])), "Set2")
# 
# plotTree(sub_tree,lwd=1)
# nodelabels(node=1:sub_tree$Nnode+Ntip(sub_tree),
#            pie=fitER[[1]]$lik.anc,piecol=cols,cex=0.4)
# add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
#     y=0.8*par()$usr[3],fsize=0.8)

### Still not sure yet, so need to make sure in the end...
#write("Manipulating data for easy output...", stderr())
anc.states <- fitER
for (i in 1:length(fitER)) {
  anc.states[[i]] <- as.data.frame(fitER[[i]]$lik.anc)
  anc.states[[i]]$node <- as.character(1:tree$Nnode+Ntip(tree))
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
  tip_labels[[i-2]] <- tip_labels[[i-2]][order(match(tip_labels[[i-2]]$node, tree$tip.label)),]
  tip_labels[[i-2]]$node <- as.character(as.numeric(1:length(tree$tip.label)))
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
    dplyr::rename(to = "node", date = "x")
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
  write("ASR not performed - results will rely on sampling dates.")
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

  if (exists("asr", envir = globalenv())) { #Only use if ASR incorporated above
    asr2 <- bind_rows(asr, .id = "field")
    asr2$to <- as.integer(asr2$to)

  cluster_data <- mclapply(clusters, function(x) {
    merge(x, asr2, by = "to") %>%
      dplyr::rename(parent = "from", node = "to") %>%
      dplyr::select(parent, everything())
  }, mc.cores=numCores)
  } else {
    
    ## Need to create simething similar to asr table for tips only
    
     cluster_data <- mclapply(clusters, function(y) {
       dplyr::rename(y, parent = "from", node = "to") %>%
         dplyr::select(parent, everything()) %>%
         left_join(., dplyr::select(node_dates, node, x), by="node") %>%
         dplyr::rename(date=x)
  }, mc.cores=numCores)
}

return(cluster_data)

} # End dataManip() function

cluster_data <- dataManip(clusters)

## Now DYNAMITE determines if clusters are related by connecting the children of each cluster to the root of remaining clusters
connectClust <- function(cluster_data) {

cluster_data <- mclapply(cluster_data, as_tibble, mc.cores=numCores)
dup_cluster_data <- cluster_data
full_tree <- as_tibble(sub_tree)

for (i in seq_along(cluster_data)) {
  cluster_data[[i]]$birth_origin <- NA
  for (j in seq_along(dup_cluster_data[-i])) {
    ## If the parent of origin node of one cluster (found by referencing full tree) is a child node in one of the other clusters... 
    if (tidytree::parent(full_tree, dup_cluster_data[[j]]$parent[1])$parent %in% 
      cluster_data[[i]]$node |
      ## Or if the parent of the parent of the origin node of one cluster (found by referencing full tree) is a child node in one of the other clusters... (meaning clusters are allowed to be separated by up to two branches)
      tidytree::parent(full_tree, 
                  parent(full_tree, dup_cluster_data[[j]]$parent[1])$parent)$parent %in%
      cluster_data[[i]]$node) {
    cluster_data[[j]]$birth_origin == paste0("cluster_", cluster_data[[i]]$parent[1] )
    } #End if statement
  } # End for loop along duplicate list
} # End for loop along original list
return(cluster_data)
} #End connectClust function

cluster_data <- connectClust(cluster_data)


# Added statistics, such as effective population size, basic reproductive number, Pybus's gamma, and others have been included below. Please note that Skygrowth-based estimates of Re and/or R0 are not accurate for clusters with true R0>5 or highly-sampled transmission clusters.

write("Claculating some tree statistics....Please note that Skygrowth-based estimates of Re and/or R0 are not accurate for clusters with true R0>5 or highly-sampled transmission clusters.")
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
  } else {Re <- list(NULL)}
  Re <- do.call("rbind",Re)
}

Ne_fulltree <- calculateNe(time_tree)
Re_fulltree <- calculateRe(Ne_fulltree)

gatherStats <- function(cluster) {
  tmp <- keep.tip(time_tree, cluster$node)
  speciation_rate <-yule(as.phylo(tmp))$lambda
  Pybus_gamma <- ltt(as.phylo(tmp), plot=FALSE, gamma=TRUE)$gamma ## Maybe also consider plottng full ltt as well?
  Ne <- calculateNe(as.phylo(tmp))
  Re <- calculateRe(Ne, conf.level = 0.95)
  
  TMRCA <- min(cluster$date)
  timespan <- max(cluster$date) - min(cluster$date)
  cluster_origin <- cluster$birth_origin[1]
  
  ## Put distribution data into list of tables
  distdata <- dplyr::select(cluster, field, trait, date) %>%
    dplyr::group_by(field) %>%
    dplyr::group_split()

  fields <- do.call(rbind, mclapply(distdata, function(x) {
    x$field[1]
  }, mc.cores=numCores)) ## For some reason order is not the same as order of column names in metadata

  names(distdata) <- fields

  ## Dynamics
  dynamics <- list()
  dynamics$birth <- NA
  dynamics$death <- NA
  dynamics$rapid <- NA
  if (!is.na(cluster$birth_origin[1])) {
    dynamics$birth <- "birth"
  } else {dynamics$birth <- NA}
  if (max(cluster$date) < num.mrsd) {
    dynamics$death <- "death"
  } else {dynamics$death <- NA}
  if (isTRUE(Re$mean_Re[1] > Re_fulltree$mean_Re[1])) {
    dynamics$rapid <- "rapid"
  } else {dynamics$rapid <- NA}
  return(tibble::lst(distdata, TMRCA, timespan, cluster_origin, speciation_rate, Pybus_gamma, Ne, Re, dynamics))
}

cluster_stats <- mclapply(cluster_data, gatherStats, mc.cores=numCores)

# Visualization is based on the resulting list (per cluster) of data

write("Data are now being exported as 'cluster_data.RDS' and 'tree_data.tree.'")
saveRDS(cluster_stats, "cluster_data.RDS")

# test <- cluster_list[[1]][[1]][[2]]
# test$trait <- as.factor(test$trait)
# ggplot(data=test, aes(x=date, color=trait)) +
#   geom_histogram(aes(y=..count..))

## Tree export
exportTree <- function(time_tree, cluster_data) {
  t.tbl <- as_tibble(time_tree)
  t.tbl <- cbind(t.tbl, gg_tree$data$x)
  names(t.tbl)[5] = "dates"
  t.tbl$cluster_id <- ""
  for (i in seq_along(cluster_data)) {
    for (j in seq_along(t.tbl$node)) {
      if (t.tbl$node[j] %in% cluster_data[[i]]$node) {
        t.tbl$cluster_id[j] = names(cluster_data)[i]
      } else{t.tbl$node[j] = t.tbl$node[j]}
    }
  }

  t.tbl$dates <- date_decimal(t.tbl$dates)
  t.tbl$dates <- as.Date(t.tbl$dates)
  t.tbl$label <- paste(t.tbl$label, t.tbl$dates, sep="|") # Needed for nextstrain, but ignoring for now
  class(t.tbl) = c("tbl_tree", class(t.tbl))

  t2 <- as.treedata(t.tbl)
  return(t2)
}
exported_tree <- exportTree(time_tree, cluster_data)
write.beast(exported_tree, "tree_data.tree")

rt1 <- Sys.time()
rt1-rt0

