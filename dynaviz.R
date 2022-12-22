
.packages <- c("dplyr", "tidyr", "treeio", "ggtree", "shiny", "cowplot", "ggplot2", "RColorBrewer", "gridExtra")
.github_packages=c("cutr")
               
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) BiocManager::install(.packages[!.inst], type = "source", checkBuilt = TRUE)
.inst_git <- .github_packages %in% installed.packages()
if(length(.github_packages[!.inst_git]) > 0) BiocManager::install(.github_packages[!.inst_git], type = "source", checkBuilt = TRUE)

invisible(lapply(.packages, library, warn.conflicts=FALSE, character.only=TRUE))
invisible(lapply(.github_packages, library, warn.conflicts=FALSE, character.only=TRUE))


`%notin%` <- Negate(`%in%`)

## Read in trait distributions results from DYNAMITE
traits <- read.csv(list.files(pattern="dynamite_trait_*")) %>%
  mutate(DATE=as.Date(DATE)) %>%
  mutate(trait = if_else(trait=="", "NA", trait)) %>%
  mutate(trait = if_else(field=="Age" & trait<=18, "1-18", 
                         if_else(field=="Age" & trait>18 & trait<=35, "18-35",
                                 if_else(field=="Age" & trait>35 & trait<=55, "35-55",
                                         if_else(field=="Age" & trait>55, ">55", trait))))) ## Group age into discrete categories
traits = filter(traits, field %notin% c("genbank_accession", "accession","epiweek", "author")) %>%
  select(-birth_origin)

traits = pivot_wider(traits, names_from=field, values_from=trait)

for (c in 4:ncol(traits)) {
  if(isTRUE(!is.na(as.numeric(unlist(traits[,c])[1])))) {
    traits[,c] = as.numeric(unlist(traits[,c]))
    traits[,c] = smart_cut(unlist(traits[,c]),
                                                               quantile(unlist(traits[,c]), na.rm=T), 
                                                               labels = ~paste0(.y[1],'-',.y[2]-1), 
                                                               simplify = FALSE)
  }
}

num_data = cbind(select(traits, cluster_id, ID), traits[, !sapply(traits, is.character)]) %>%
  gather(field, trait, -c(DATE, cluster_id, ID))

chr_data = cbind(select(traits, DATE), traits[, sapply(traits, is.character)]) %>%
  gather(field, trait, -c(DATE, cluster_id, ID))

stats <- read.csv(list.files(pattern="stats")) %>%
  gather(field, trait, PD:Oster) 
  



## Read in timed tree from DYNAMITE
tree_file = list.files(pattern="dynamite_timetree*")
if(isTRUE(length(tree_file)!=0)) {
  tree = read.beast(tree_file)
} else {
  tree = read.beast(list.files(pattern="dynamite_subtree*"))
  
}


## Option to save pdfs of the distribution plots

p1 = ggplot(chr_data) +
  geom_density(aes(x=DATE, fill=trait), alpha=0.8) +
  facet_wrap(cluster_id~field, scales="free_y", nrow=length(unique(chr_data$cluster_id)), ncol=length(unique(chr_data$field))) +
  theme(text = element_text(size=20)) +
  theme_minimal()

p2 = ggplot(num_data) +
  geom_density(aes(x=DATE, fill=trait), alpha=0.8) +
  facet_wrap(cluster_id~field, scales="free_y", nrow=length(unique(num_data$cluster_id)), ncol=length(unique(num_data$field))) +
  theme(text = element_text(size=20), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20) ) +
  theme_minimal()

ggsave(grid.arrange(p1, p2, ncol=2), "trait_distributions.pdf", width=11, height=8.5, units="in")

p3 = ggplot(stats, aes(x=trait, fill=field)) +
  geom_histogram(alpha=0.8) +
  facet_wrap(cluster_id~field, scales="free_y", nrow=length(unique(stats$cluster_id)), ncol=length(unique(stats$field))) +
  theme(text = element_text(size=20), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20) ) +
  theme_minimal()

ggsave(p3, "tree_stats_by_cluster.pdf", width=11, height=8.5, units="in")



## Plotting of the clusters on the trees
## Find minimum and maximum dates for R Shiny plotting
min_date <- min(traits$DATE)
max_date <- max(traits$DATE)

grouping <- data.frame(Cluster=unique(traits$cluster_id))

rownames(grouping) <- unique(traits$ID)

n1 <- length(unique(grouping$Cluster))

df <- grouping %>% mutate(label=rownames(.))

tree2 <- full_join(tree, df, by = "label")

col <- setNames(brewer.pal(n1, "Set3"), unique(as.character(grouping$Cluster)))

t <- ggtree(tree2, mrsd=max_date) +
  theme_tree2()
#  geom_tiplab(aes(label=VAX), align=TRUE, linetype="blank") 
#  geom_point2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 90), color="dimgrey") +
#  geom_treescale(linesize = 1.5, width=5E-04, x=0.003, y=4700, offset=50)

if(isTRUE(length(tree_file)!=0)) {
  p_full <- gheatmap(t, grouping, colnames=F,
                   colnames_position = "top",
                   offset = 0.01, width=1.25,
                   color=rgb(1, 1, 1, alpha=0.0001)) +
  scale_fill_manual(values=col, na.translate=F) +
  scale_x_ggtree() +
  scale_x_date(date_breaks="1 month", date_labels = "%b") +
  scale_y_continuous(expand=c(0, 0.3)) +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.position = "none",
        text=element_text(size=20)) +
  ylim(NA, 500)

ggsave(plot=p_full, "cluster_tree.pdf", width = 10, height = 11, units = "in", limitsize = FALSE)

leg1 <- cowplot::get_legend(p_full)
plot_grid(leg1)
ggsave(plot=last_plot(), "legend.pdf")

} else {
  p_full <- gheatmap(t, grouping, colnames=F,
                     colnames_position = "top",
                     offset = 0.01, width=1.25,
                     color=rgb(1, 1, 1, alpha=0.0001)) +
    scale_fill_manual(values=col, na.translate=F) +
    scale_x_ggtree() +
    scale_y_continuous(expand=c(0, 0.3)) +
    theme(panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          legend.position = "none",
          text=element_text(size=20)) +
    ylim(NA, 500)
  
  ggsave(plot=p_full, "cluster_tree.pdf", width = 10, height = 11, units = "in", limitsize = FALSE)
  
  leg1 <- cowplot::get_legend(p_full)
  plot_grid(leg1)
  ggsave(plot=last_plot(), "legend.pdf")
}





## Use R Shiny package for interactive plotting of distribution data
ui <- fluidPage(
  titlePanel(title=h4("Cluster trait distributions", align="center")),
  sidebarPanel( 
    sliderInput("date_spec", "Date:", min = min_date, max = max_date, step=7, value=c(min_date, max_date))),
  mainPanel(plotOutput("plot2")))

server <- function(input,output){
  
   dat <- reactive({
     test <- traits[traits$DATE %in% seq(from=min(input$date_spec),to=max(input$date_spec),by=1),]
#     print(test)
    test
  })
  
  output$plot2<-renderPlot({
    grid.arrange(p1,p2,p3, ncol=3)
  }, height = 8000,width = 4000)
  }


shinyApp(ui, server)



