## Initialize stable working environment and store time of intiation
rm(list=ls())

#setwd("/Users/macbook/Dropbox (UFL)/DYNAMITE/HIVdynamite/nosoi_simulations")
# List of packages for session
.packages <-  c("phytools", "ape", "parallel", "ggplot2", "viridis", "igraph", "ggnetwork", "ggpubr", "ggtree", "treeio", "ape", "remotes", "dplyr", "plyr") # May need to incorporate code for familyR (https://rdrr.io/github/emillykkejensen/familyR/src/R/get_children.R) i fno longer supported.
github_packages <- c("slequime/nosoi", "emillykkejensen/familyR") 

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
.inst_github <- .packages %in% installed.packages()
## Install GitHub packages(if not already installed)
if(length(github_packages[!.inst_github]) > 0) try(remotes::install_github(github_packages[!.inst_github]))
if(length(github_packages[!.inst_github]) > 0) try(devtools::install_github(github_packages[!.inst_github]))

# Load packages into session 
lapply(.packages, require, character.only=TRUE)
lapply(gsub(".+\\/(.+)", "\\1", github_packages), require, character.only=TRUE)

### Set seeds #############################################################################
# get simulation parameters

#args = commandArgs(trailingOnly=TRUE)
sim_index = as.numeric(args[1]) # arg 1 fraction new
seeds = readLines('seeds.txt')
set.seed(seeds[sim_index])
#sim_index="test"
numCores = detectCores()


## Matrix generation #######################################################################
traits <- data.frame(location=rbind('A','B','C', 'D', 'E'))
Q <-list()
for (column in 1:ncol(traits)) {
  suppressWarnings({ ## We know it fills diagonals with NAs
    Q[[column]] <- diag(unique(traits[,column]), nrow = length(unique(traits[,column])))
  })
  #    diag(Q[[column]]) = 1-nrow(Q[[column]])
  diag(Q[[column]]) = 0
  non.diag <- 1/(nrow(Q[[column]])-1)
  Q[[column]][lower.tri(Q[[column]])] <- non.diag
  Q[[column]][upper.tri(Q[[column]])] <- non.diag
  colnames(Q[[column]]) <- rownames(Q[[column]]) <- unique(traits[,column])
  ## Modify Q matrix to allow for birth of cluster F from cluster E
#  Q[[column]][4,] <- c(rep(0.0, 4),1.0) # E only gives rise to F
#  Q[[column]][1,] <- c(0.0,rep(1/3,3),0.0) # A cannot give rise to F; the remaining clusters can stay at 0.2 because they have a 0 probability of leaving (see below)
}


Q <- plyr::compact(Q)
names(Q) <- colnames(traits[-1])


# OR
#create matrix layout
# df.matrix = df %>% spread(To, value=N, fill= 0)
# df.matrix = as.matrix(df.matrix[-1])
# rownames(df.matrix) = colnames(df.matrix)
# 
# #Get probabilities (beware, rows should sum up to 1)
# df_transpose = t(df.matrix)
# probabilities <- apply(df_transpose, 1, function(i) i/sum(i))
# transition.matrix = t(probabilities)



# #pExit daily probability for a host to leave the simulation (either cured, died, etc.).
p_Exit_fct  <- function(t, t_sero){ 
  if(t <= t_sero){p=0}
  if(t > t_sero){p=0.95} 
  return(p)
}
t_sero_fct <- function(x){rnorm(x,14,3)} #Approximately 14 days

#pMove probability (per unit of time) for a host to do move, i.e. to leave its current state (for example, leaving state “A”). 
#It should not be confused with the probabilities extracted from the structure.matrix, which represent the probability to go 
#to a specific location once a movement is ongoing (for example, going to “B” or “C” while coming from “A”).
p_Move_fct  <- function(t, current.in, host.count){
  
  if(current.in=="A" & host.count < 2){return(0)}
  if(current.in=="A" & host.count >= 2){return(0.00075)} # Wait a couple of weeks before initiation of clusters
  if(current.in=="B"){return(0)}
  if(current.in=="C"){return(0)}
#  if(current.in=="D"){return(0.0015)}
  if(current.in=="D"){return(0)}
  if(current.in=="E"){return(0)}
}
#  t_clust_fct <- function(x){rnorm(x,mean = 3,sd=1)}

n_contact_fct = function(t, current.in, host.count){
    
  if(current.in=="A") {p=abs(round(rnorm(1, 20, 1), 0))}
  if(current.in=="B") {p=abs(round(rnorm(1, 4, 1), 0))}
  if(current.in=="C") {p=abs(round(rnorm(1, 4, 1), 0))}
  if(current.in=="D") {p=abs(round(rnorm(1, 8, 1), 0))}
  if(current.in=="E") {p=abs(round(rnorm(1, 4, 1), 0))}
  
#  if(current.in=="C"){p=4*10*exp(0.0005*host.count)/((10-4)+4*exp(0.0005*host.count))}

  return(p)
}

#pTrans represents the probability of transmission over time (when a contact occurs).
# in the form of a threshold function: before a certain amount of time since initial infection, the host does not transmit (incubation time, which we call t_incub), and after that time, it will transmit with a certain (constant) probability (which we call p_max). This function is dependent of the time since the host’s infection t.
p_Trans_fct <- function(t, current.in, host.count, t_incub){
  if(t < t_incub){p=0}
  if(t >= t_incub & current.in=="A"){p=0.015}
  if(t >= t_incub & current.in=="B"){p=0.125}
  if(t >= t_incub & current.in=="C"){p=0.175}
  if(t >= t_incub & current.in=="D"){p=0.0875}
  if(t >= t_incub & current.in=="E"){p=0.275}
  return(p)
}

t_incub_fct <- function(x){rnorm(x,mean = 5,sd=2)} #Approximately 4.2 days
#p_max_fct <- function(x){rbeta(x,shape1 = 1,shape2=8)}# Mean of roughly 0.10

# Starting the simulation ------------------------------------

#set.seed(101)

SimulationSingle <- nosoiSim(type="single", # Number of hosts
                             popStructure="discrete", #discrete or continuous
                             structure.matrix = Q[[1]], # prob matrix defined above (row sums to 1, with diags as zero)
                             length.sim = 365, # Max number of time units (can be days, months, weeks, etc.)
                             max.infected = 10000, #maximum number of individuals that can be infected during the simulation.
                             init.individuals = 1, #number of individuals (an integer above 1) that will start a transmission chain. Keep in mind that you will have as many transmission chains as initial individuals, which is equivalent as launching a number of independent nosoi simulations.
                             init.structure = "A",
                             
                             pExit = p_Exit_fct,
                             param.pExit=list(t_sero=t_sero_fct),
                             timeDep.pExit=FALSE,
                             diff.pExit=FALSE,
                             
                             pMove = p_Move_fct,
                             hostCount.pMove=TRUE,
                             param.pMove=NA,
                             timeDep.pMove=FALSE,
                             diff.pMove=TRUE,
                             
                             nContact=n_contact_fct,
                             hostCount.nContact=TRUE,
                             param.nContact = NA,
                             timeDep.nContact=FALSE,
                             diff.nContact=TRUE,
                             
                             # nContact=time_contact,
                             # hostCount.nContact=TRUE,
                             # param.nContact = NA,
                             # timeDep.nContact=FALSE,
                             # diff.nContact=TRUE,
                             
                             pTrans = p_Trans_fct,
                             hostCount.pTrans=TRUE,
                             param.pTrans = list(t_incub=t_incub_fct),
                             timeDep.pTrans=FALSE,
                             diff.pTrans=TRUE,
                             
                             prefix.host="S",
                             print.progress=TRUE,
                             print.step=100)

sum_sim <- summary(SimulationSingle)
cumulative.table <- getCumulative(SimulationSingle)
dynamics.table <- getDynamic(SimulationSingle)
cum.p <- ggplot(data=cumulative.table, aes(x=t, y=Count)) + geom_line() + theme_minimal() + 
  labs(x="Time (t)",y="Cumulative count of infected hosts") + scale_y_log10()
cum.p.c <- ggplot(data=dynamics.table, aes(x=t, y=Count, color=state)) + geom_line() + theme_minimal() + 
  labs(x="Time (t)",y="Number of active infected hosts") + scale_y_log10()

#ggpubr::ggarrange(cum.p, cum.p.c, widths = 2, heights = 1, legend="right")
#ggsave("cov_network_cuminf_max10000.png", plot=last_plot())


## Grab tree #########################################################################################
save.tree <- function(){
  sim.tree <- getTransmissionTree(SimulationSingle)
  # ggtree(test.nosoiA.tree, color = "gray30") + geom_nodepoint(aes(color=state)) + geom_tippoint(aes(color=state)) + 
  #   theme_tree2() + xlab("Time (t)") + theme(legend.position = c(0,0.8), 
  #                                            legend.title = element_blank(),
  #                                            legend.key = element_blank()) 
  #set.seed(5905950) 
  
  get_sample <- function(SimulationSingle) {
    table.hosts <- getTableHosts(SimulationSingle, pop="A")
    sampled.hosts <- sample(table.hosts$hosts.ID, round(2.5*last(cum.p$data$t)), replace=F)
    int.nodes <- sample((Ntip(sim.tree@phylo)+2):(Ntip(sim.tree@phylo)+sim.tree@phylo$Nnode)) #randomize order of internal nodes
    tcs.list <- rep(list(NULL),3)
    ## move along internal nodes and stop when you have three clades with 10-20 individuals.
    s <- round(rnorm(10,15,3))
    for (i in seq_along(tcs.list)){
      while(is.null(tcs.list[[i]])){
        true.cluster <- extract.clade(sim.tree@phylo, sample(int.nodes, 1))
        true.cluster <- subset(sim.tree@data, sim.tree@data$host %in% true.cluster$tip.label)
        if(length(true.cluster$host) %in% s &
           length(grep("A", true.cluster$state)) >= 0.95*length(true.cluster$state)) {
           tcs.list[[i]] <- true.cluster
        }
      }
    }
    ## Add these individuals to list of randomly sampled invididuals
    sampled.hosts <- unique(c(sampled.hosts, unname(unlist(lapply(tcs.list, "[", 'host')))))
    ## Extract tree for list of individuals from the full simulation tree
    sampled.tree <- sampleTransmissionTreeFromExiting(sim.tree, sampled.hosts)
    assign("tcs.list", tcs.list, envir = globalenv()) # Remember now a tibble
    return(sampled.tree)
  }
  
  sampled.tree <- get_sample(SimulationSingle)
  
  ## Save sampled trees as RDS for generation of ctl file for indelible
  #saveRDS(file="sampled_trees.rds", sampled.trees)
  
  sampled.tree <- as_tibble(sampled.tree)
  # dirtrans_clusters <- mclapply(tcs.list, function(x) {
  #   subset(sampled.tree, sampled.tree$label %in% x$host)}, 
  #   mc.cores=numCores)
  
  for (i in seq_along(tcs.list)) {
    tcs.list[[i]]$host <- paste(tcs.list[[i]]$host, tcs.list[[i]]$state, sep="_") }
  
  sampled.tree$label <- paste(sampled.tree$label, sampled.tree$state,sep="_")
  sampled.tree <- as.treedata(sampled.tree)
  
  
  write.beast(sampled.tree, paste0('sim_', sim_index, "_sampled_10000.tree")) #### NEED THIS OUTPUT####################################
  
  ## Export mean R0 for entire simulation for benchmarking
  write.table(data.frame(mean_R0=sum_sim$R0$R0.mean,
                         upper = quantile(sum_sim$R0$R0.dist, 0.95),
                         lower = quantile(sum_sim$R0$R0.dist, 0.05)), 
              paste0("R0_", sim_index, ".tab"), sep='\t', quote=F, row.names = F)
  
  saveRDS(tcs.list, paste0("dirtrans_clusters_", sim_index, ".rds"))
  
  text <- write.tree(sampled.trees@phylo)
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
  seqdata <- seqgen(opts="-s8.219178e-04 -mGTR -i0.601 -a2.35 -r0.32512,1.07402,0.26711,0.25277,2.89976,1.00000 -f0.299,0.183,0.196,0.322", rooted.tree=text)
  seqdata <- as.vector(seqdata)
  
  write.table(seqdata, paste0("seqdata_", sim_index, ".phy"), 
              quote=F, row.names = F, col.names = F)
    }


## If any cluster exceeds background population in number, do not output results; however, if not, proceed with tree extraction.
max.t <- max(dynamics.table$t)
if(isTRUE(all(dynamics.table$Count[dynamics.table$state == 'A' & 
                                   dynamics.table$t == max.t] > 
              dynamics.table$Count[dynamics.table$state %in% c('B', 'C', 'D', 'E') & 
                                   dynamics.table$t == max.t]))) {
  save.tree()
} else {NULL}
## End ###################################################################################

require(phyclust)


## Network ##############
# data.sim <- getTableHosts(SimulationSingle, "A")
# saveRDS(data.sim, "data.sim.rds")
# graph.simA <- graph.data.frame(data.sim[-1,c("inf.by","hosts.ID")],directed=T,vertices = data.sim)
# 
# 
# graph.simA.network <- ggnetwork(graph.simA, layout = with_kk()) #using ggnetwork to provide the layout
# 
# #plotting the network (color is time of infection)
# plot1 <- ggplot(graph.simA.network, aes(x = x, y = y, xend = xend, yend = yend)) +
#   geom_edges(color = "grey70",arrow = arrow(length = unit(0.3, "lines"), type = "open")) +
#   geom_nodes(aes(color=inf.time)) + scale_color_viridis(name="Time of infection",option = "plasma") +
#   theme_blank()
# 
# #plotting the network (color is state the host was infected in)
# plot2 <- ggplot(graph.simA.network, aes(x = x, y = y, xend = xend, yend = yend)) +
#   geom_edges(color = "grey70",arrow = arrow(length = unit(0.3, "lines"), type = "open")) +
#   geom_nodes(aes(color=inf.in)) + scale_color_discrete(name="State of infection") +
#   theme_blank()
# 
# ggpubr::ggarrange(plot1,plot2,widths = 2, heights = 1, legend="bottom")
# 
# ggsave("cov_network_sim3_max10000.png", plot=last_plot())




#### Visualization and animations #######################################################
# library(networkDynamic)
# library(ndtv)
# library(dplyr)
# 
# SimulationA <- getTableHosts(test.nosoiA, "A")
# SimulationA$numID <- as.numeric(as.factor(SimulationA$hosts.ID)) #gets a unique numeric ID for each host
# SimulationA[is.na(out.time),out.time:=(test.nosoiA$total.time)] #active hosts get a "end time" (i.e. the end of the simulation)
# 
# #Create "edges" table (line between nodes)
# edges = SimulationA[,c("hosts.ID","inf.by","inf.time","out.time")]
# colnames(edges) = c("to","from","infected","recovered")
# edges$end = rep(test.nosoiA$total.time+1,nrow(edges))
# 
# #Create "vertex" table (nodes)
# vertex2 = edges %>% group_by(to) %>% mutate(vertex=to)
# vertex2 = vertex2[,c("vertex","infected","recovered","end")]
# vertex2 = as.data.frame(vertex2)
# vertex2$numericID = c(1:nrow(vertex2))
# vertex2$infected = as.numeric(vertex2$infected)
# vertex2$recovered = as.numeric(vertex2$recovered)
# vertex2$end = as.numeric(vertex2$end)
# 
# vertex3=vertex2[,c(1,5)]
# colnames(vertex3) = c("to","numeric.ID")
# 
# test=left_join(edges, vertex3, by=c("to"))
# colnames(vertex3) = c("from","numeric.ID")
# 
# edges2=left_join(test, vertex3, by=c("from"))
# colnames(edges2) = c("to","from","infected","recovered","end","to.num","from.num")
# 
# edges2[1,7] = 1
# 
# #Create main network using networkDynamic
# test.network = networkDynamic(vertex.spells=vertex2[,c(2,4,5)],edge.spells =edges2[,c(3,5,7,6)])
# 
# #Animated plot where "active" hosts are red, inactive gray
# activate.vertex.attribute(test.network,'color','gray',onset=0,terminus=47)
# activate.vertex.attribute(test.network,'color','red',onset=vertex2$infected,terminus=vertex2$recovered)
# 
# #Render animation
# library(tsna)
# saveVideo(render.animation(test.network,vertex.col = "color",
#                            render.par=list(show.time=TRUE),
#                            edge.col = "gray60",
#                            displaylabels=FALSE, render.cache='none',mode="kamadakawai"),video.name="transmission-nosoi.mp4")
# 
# 
# ## Visualize summary of transitions between states
# #Loading data from the simulation
# 
# test.nosoiA.data <- getTableState(test.nosoiA)
# 
# #Loop to get the number of infected in each state per unit of time
# results <- data.frame()
# for (i in 1:max(test.nosoiA.data$time.from)) {
#   temp <- subset(test.nosoiA.data, time.from <= i & c(time.to > i |is.na(time.to)))[,c("hosts.ID","state")]
#   temp$time <- i
#   results <- data.table::rbindlist(c(list(temp),list(results)))
# }
# 
# test2 <- results %>% group_by(time,state) %>% summarise(N=length(hosts.ID)) %>% left_join(layout,by=c("state"="name"))
# 
# #Loop to get the transitions between states and their time
# results2=data.frame()
# for (i in unique(test.nosoiA.data$hosts.ID)) {
#   subset.current <- test.nosoiA.data[i]
#   if(nrow(subset.current) > 1){
#     for (j in 1:(nrow(subset.current)-1)){
#       temp <- data.table(hosts.ID=i,time=as.integer(subset.current[j]$time.to),from=subset.current[j]$state,to=subset.current[j+1]$state)
#       results2 <- data.table::rbindlist(c(list(temp),list(results2)))
#     }}}
# 
# test3 <- results2 %>% group_by(time,from,to) %>% summarise(N=length(hosts.ID)) %>% left_join(layout,by=c("from"="name")) %>% left_join(layout,by=c("to"="name"),suffix=c("from","to"))
# 
# #Animated plot (using gganimate):
# animated.plot <- ggplot() +
#   geom_point(data=test2,aes(x,y, color=state,size=N)) +
#   geom_curve(data=test3, aes(x=xfrom,y=yfrom,xend=xto,yend=yto),arrow = arrow(length = unit(0.03, "npc"),type = "closed"),curvature = 0.2,color="gray50") +
#   scale_color_viridis(guide=FALSE,discrete=TRUE) + theme_blank() + ylim(-0.5,1.2) + xlim(-0.5,1.2) + geom_text(data=layout,aes(x=x,y=y-0.2,label = name), size = 6,color="black") +
#   geom_text(data=test2,aes(x=x,y=y+0.2,label = N), size = 6,color="black") +
#   transition_states(time) + scale_size_continuous(guide=FALSE,range = c(5, 30)) +
#   labs(title = "Time: {closest_state}")
# 
# animate(animated.plot, nframes=test.nosoiA$total.time*2+10,duration=40,end_pause=10)




# How frequently are direct transmission events sampled?
# host.pairs <- expand.grid(sampled.tree@data$host, sampled.tree@data$host, stringsAsFactors = F)
# node.pairs <- expand.grid(as.integer(sampled.tree@data$node), as.integer(sampled.tree@data$node))
# true.host.pairs <- cbind(table.hosts[,2], table.hosts[,1]) # Needs to be switched because "inf.by" is currently in second column
# for (i in 1:nrow(host.pairs)) {
#   for (j in 1:nrow(sampled.tree@phylo$edge)) {
#     for (k in 1:nrow(true.host.pairs)) {
#       if (isTRUE(node.pairs[i,] == sampled.tree@phylo$edge[j,]) & #For node pairs connected by edges in tree,
#           isTRUE(host.pairs[i,] %in% true.host.pairs[k,])) { # If the host pairs represented by those nodes are connected by infection in the simulation table,
#         host.pairs$transmission[i] <- "direct"
#       } else {
#         host.pairs$transmission[i] <- "indirect"
#          } # End if-else statement
#     } # End loop along table.hosts
#   } # End loop along edges
# } # End loop along host and node pairs
