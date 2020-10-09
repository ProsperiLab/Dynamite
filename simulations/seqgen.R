## The following uses Seq-Gen R version to evolve sequences along trees

## Initialize stable working environment and store time of intiation
rm(list=ls())


.packages <-  c("ggtree", "treeio", "phyclust") # May need to incorporate code for familyR (https://rdrr.io/github/emillykkejensen/familyR/src/R/get_children.R) i fno longer supported.

# Load packages into session 
lapply(.packages, require, character.only=TRUE)

## Apply benchmarking to each simulation
args = commandArgs(trailingOnly=TRUE)
sim_index = as.numeric(args[1]) # arg 1 fraction new
#sim_index <- "test" # test


tree_list = list.files(pattern=paste0("sim_", sim_index, ".+\\.tree$"))
print("reading in tree...")
tree = lapply(tree_list, read.beast)[[1]]

text <- write.tree(tree@phylo)
strip.nodelabels<-function(text){
  obj<-strsplit(text,"")[[1]]
  cp<-grep(")",obj)
  csc<-c(grep(":",obj),length(obj))
  exc<-cbind(cp,sapply(cp,function(x,y) y[which(y>x)[1]],y=csc))
  exc<-exc[(exc[,2]-exc[,1])>1,]
  inc<-rep(TRUE,length(obj))
  if(nrow(exc)>0) for(i in 1:nrow(exc)) 
    inc[(exc[i,1]+1):(exc[i,2]-1)]<-FALSE
  paste(obj[inc],collapse="")
}

text <- strip.nodelabels(text)
text <- read.tree(text=text)
text <- multi2di(text)
seqdata <- seqgen(opts="-s8.219178e-05 -mGTR -i0.601 -a2.35 -r0.32512,1.07402,0.26711,0.25277,2.89976,1.00000 -f0.299,0.183,0.196,0.322", rooted.tree=text)
seqdata <- as.vector(seqdata)

write.table(seqdata, paste0("seqdata_", sim_index, ".phy"), 
            quote=F, row.names = F, col.names = F)
