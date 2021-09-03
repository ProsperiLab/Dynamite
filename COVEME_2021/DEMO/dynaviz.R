
.packages <- c("dplyr", "tidyr", "treeio", "ggtree", "shiny", "cowplot", "ggplot2", "RColorBrewer")
               
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) BiocManager::install(.packages[!.inst], type = "source", checkBuilt = TRUE)
invisible(lapply(.packages, library, warn.conflicts=FALSE, character.only=TRUE))


## Read in trait distributions results from DYNAMITE
traits <- read.csv(list.files(pattern="dynamite_trait_*")[[1]]) %>%
  mutate(DATE=as.Date(DATE)) %>%
  mutate(trait = if_else(trait=="", "NA", trait)) %>%
  mutate(trait = if_else(field=="Age" & trait<=18, "1-18", 
                         if_else(field=="Age" & trait>18 & trait<=35, "18-35",
                                 if_else(field=="Age" & trait>35 & trait<=55, "35-55",
                                         if_else(field=="Age" & trait>55, ">55", trait))))) %>% ## Group age into discrete categories
  mutate(trait=if_else(trait=="Gainesville ", "Gainesville", trait)) # Mistake in reporting location by lab

## Read in timed tree from DYNAMITE
tree <- read.beast(list.files(pattern="dynamite_timetree*")[[1]])


## Option to save pdfs of the distribution plots

ggplot(traits) +
  geom_density(aes(x=DATE, fill=trait), alpha=0.8) +
  facet_wrap(cluster_id~field, scales="free_y", nrow=length(unique(traits$cluster_id)), ncol=length(unique(traits$field))) +
  theme(text = element_text(size=20), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20) ) +
  theme_minimal()

ggsave(plot=last_plot(), "trait_distributions.pdf", width=8, height=8, units="in")


## Plotting of the clusters on the trees
## Find minimum and maximum dates for R Shiny plotting
min_date <- min(traits$DATE)
max_date <- max(traits$DATE)


traits_wide  <- spread(traits, field, trait) %>%
  select(ID, cluster_id, Location, Vaccinated)


grouping <- data.frame(Cluster=traits_wide$cluster_id, Origin=traits_wide$Location, VAX=traits_wide$Vaccinated)

rownames(grouping) <- traits_wide$ID


n1 <- length(unique(grouping$Cluster))
n2 <- length(unique(grouping$Origin))
n3 <- length(unique(grouping$VAX))

df <- grouping %>% mutate(label=rownames(.)) %>%
  mutate(VAX = if_else(VAX=="Y",
                       gsub(".+(FL[A-Za-z\\-]+\\d+).+", "\\1", label),
                       ""),
         VAX = gsub("\\-VAC", "", VAX))

tree2 <- full_join(tree, df, by = "label")

col <- setNames(c(brewer.pal(n1, "Pastel1"), brewer.pal(n2, "Paired"), brewer.pal(n3, "Reds")),
                c(unique(as.character(grouping$Cluster)), unique(grouping$Origin), unique(grouping$VAX)))

t <- ggtree(tree2, mrsd=max_date) +
  theme_tree2()
#  geom_tiplab(aes(label=VAX), align=TRUE, linetype="blank") 
#  geom_point2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 90), color="dimgrey") +
#  geom_treescale(linesize = 1.5, width=5E-04, x=0.003, y=4700, offset=50)

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


p_full <- gheatmap(t, grouping, colnames=F,
                   colnames_position = "top",
                   offset = 0.01, width=1.25,
                   color=rgb(1, 1, 1, alpha=0.0001)) +
  scale_fill_manual(values=col, na.translate=F) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0, 0.3)) +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        #        legend.position = "none",
        text=element_text(size=20)) +
  ylim(NA, 500)

leg1 <- cowplot::get_legend(p_full)
plot_grid(leg1)
ggsave(plot=last_plot(), "legend.pdf")





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
    ggplot(dat()) +
      geom_density(aes(x=DATE, fill=trait), alpha=0.8) +
      facet_wrap(cluster_id~field, scales="free_y", nrow=length(unique(traits$cluster_id)), ncol=length(unique(traits$field))) +
      theme(text = element_text(size=20), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20) ) +
      theme_minimal()
  }, height = 2000,width = 2000)
  }


shinyApp(ui, server)



