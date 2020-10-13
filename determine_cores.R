##Install and load necessary packages
## Initialize stable working environment and store time of intiation
rm(list=ls())
setwd(getSrcDirectory()[1]) #When working on the cluster
#dirname(rstudioapi::getActiveDocumentContext()$path) # If working in Rstudio

# List of packages for session
.packages <-  c("ape")  # May need to incorporate code for familyR (https://rdrr.io/github/emillykkejensen/familyR/src/R/get_children.R) i fno longer supported.

## Will need to remove install section if using on cluster ###################################
# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session 
lapply(.packages, require, character.only=TRUE)

args = commandArgs(trailingOnly=TRUE)
tree_file = args[1]

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
sub_tree <- checkFortree(tree_file)

## Need to force binary tree and to replace zero branch lengths with full bootstrap support

numClades <- function(tree) {
  bs <- as.numeric(tree$node.label)
  if (max(bs, na.rm=T) >1){
    num <- length(bs >= 90)
  } else {
    num <- length(bs[bs >= 0.90])
  }
  write(paste0("Total number of well-supported clades is ", num), stderr())
}

numClades(sub_tree)