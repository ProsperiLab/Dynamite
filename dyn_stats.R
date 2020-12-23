rm(list=ls())

# List of packages for session
.packages <-  c("dplyr", "plyr", "ggplot2", "gridExtra", "data.table", "parallel", "pbmcapply") # May need to incorporate code for familyR (https://rdrr.io/github/emillykkejensen/familyR/src/R/get_children.R) i fno longer supported.
# Load packages into session 
lapply(.packages, require, character.only=TRUE)
numCores <- detectCores()

## Set up directories ####
main_dir <- paste0("/Volumes/research/PhylodynamicLab/Projects/DYNAMITE/nosoi_sims/Final_all_direct/")
#main_dir <- paste0("/blue/salemi/brittany.rife/dynamite/simulations/")
setwd(dir = main_dir)

getAllDirs <- function() {
  sub_dirs_thresholds <- expand.grid("perc_", c("0.01", "0.05", "0.10", "0.15", "0.20"))
  sub_dirs_thresholds <- paste0(sub_dirs_thresholds$Var1, sub_dirs_thresholds$Var2)

  sub_dirs_branchwise <- expand.grid("b/addLeaves/", sub_dirs_thresholds)
  sub_dirs_branchwise <- paste0(sub_dirs_branchwise$Var1, sub_dirs_branchwise$Var2)
  
  sub_dirs_cladewise <- expand.grid("c/", sub_dirs_thresholds)
  sub_dirs_cladewise <- paste0(sub_dirs_cladewise$Var1, sub_dirs_cladewise$Var2)
  
  sub_dirs <- c(sub_dirs_branchwise, sub_dirs_cladewise)
  return(sub_dirs)
}

sub_dirs <- getAllDirs()

## Extract data ####
performance_data <- list()
performance_data <- mclapply(sub_dirs, function(x) {
  current_dir <- paste0(main_dir, x)
  setwd(dir = current_dir)
  tables_list <- list.files(pattern="*.tab")
  performance <- rbindlist(mclapply(tables_list[grep(".+performance_.+tab", tables_list)], function(y) {
  read.table(y, header=T, stringsAsFactors = F, sep='\t')}, mc.cores=numCores), fill=T)
  performance$threshold <- as.numeric(gsub(".+(0\\.\\d+)", "\\1", x))
  performance$algorithm <- gsub("([a-z]{1}).+", "\\1", x)
  return(performance)
  }, mc.cores=numCores)


performance_data <- rbindlist(performance_data, use.names=TRUE, fill=T)

##################################################################

## Need to report if multiple clusters reported as a single cluster
performance_data <- performance_data %>%
  dplyr::group_by(algorithm, threshold, sim) %>%
  dplyr::mutate(dup = ifelse(duplicated(identified_taxa, incomparables=NA), 1, 0),
                perc_dup = sum(dup)/n()) %>%
  dplyr::ungroup()
  

labeler <- function(x) {
  x$algorithm <- factor(x$algorithm, 
                        levels=c("b", "c"),
                        labels = c("branchwise", "cladewise"))
  
  x$state <- factor(x$state, 
                    levels=c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                    labels = c("A (R0~2.2)", "B (R0~3.5)", "C (R0~5.5)", "D (R0~5.5)", "E (R0~7.2)",
                               "F (R0~2.7)", "G (R0~8)", "H (R0~5.5)", "I (R0~5.5)"))
  return(x)
}
performance_data <- labeler(performance_data)


overallRates <- function(x) {
  TPR_FPR <- dplyr::select(x, algorithm, threshold, sim, proportion, perc_dup)%>%
    dplyr::group_by(algorithm, threshold, sim, perc_dup) %>%
    dplyr::summarize(num_TP_found=length(proportion[proportion>0]),
                     num_TP_not_found=length(proportion[proportion==0]),
                     total_true=num_TP_found + num_TP_not_found,
                     TPR=num_TP_found/total_true*100,
                     num_FP_found=length(proportion[is.na(proportion)]),
                     total_num_found = num_TP_found + num_FP_found,
                     FPR=num_FP_found/total_num_found) %>%
    mutate(FPR=if_else(is.nan(FPR), 0, FPR)) %>%
    mutate(FPR=FPR*100)
  return(TPR_FPR)
} 
TPR_FPR_overall <- overallRates(performance_data)
median_overall <- dplyr::group_by(TPR_FPR_overall, algorithm, threshold) %>%
  dplyr::summarize(med_TPR = median(TPR),
            med_FPR = median(FPR),
            med_perc_dup = median(perc_dup))
  
detectionRatesA <- function(x) {
TPR_A <- dplyr::select(x, algorithm, threshold, sim, state, proportion)%>%
  filter(state=='A (R0~2.2)') %>%
  dplyr::group_by(algorithm, threshold, sim) %>%
  dplyr::summarize(num_TP_found=length(proportion[proportion>0]),
                num_TP_not_found=length(proportion[proportion==0]),
                total_true=num_TP_found+num_TP_not_found,
                TPR=num_TP_found/total_true*100,
                num_FP=length(proportion[is.na(proportion)]),
                total_num_found = num_TP_found + num_FP,
                FPR=num_FP/total_num_found) %>%
  mutate(FPR=if_else(is.nan(FPR), 0, FPR)) %>%
  mutate(FPR=FPR*100)
return(TPR_A)
}
TPR_FPR_A <- detectionRatesA(performance_data)
median_A <- dplyr::group_by(TPR_FPR_A, algorithm, threshold) %>%
  dplyr::summarize(med_TPR_A = median(TPR),
                   med_FPR_A = median(FPR))


detectionRatesNonA <- function(x) {
  relative <- dplyr::select(x, algorithm, threshold, sim, state, proportion)%>%
    filter(state!='A (R0~2.2)') %>%
    dplyr::group_by(algorithm, threshold, sim, state) %>%
    dplyr::summarize(total_num_found=n(),
                     num_TP_found=length(proportion[proportion>0]),
                     num_TP_not_found=length(proportion[proportion==0]),
                     total_true=num_TP_found+num_TP_not_found,
                     TPR=num_TP_found/total_true*100,
                     num_NTPR=length(proportion[is.na(proportion)]),
                     total_num_found = num_TP_found + num_NTPR,
                     NTPR=num_NTPR/total_num_found) %>%
    mutate(NTPR=if_else(is.nan(NTPR), 0, NTPR)) %>%
    mutate(NTPR=NTPR*100)
  return(relative)
}
TPR_FPR_nonA <- detectionRatesNonA(performance_data)
median_nonA <- dplyr::group_by(TPR_FPR_nonA, algorithm, threshold) %>%
  dplyr::summarize(med_TPR_nonA = median(TPR),
                   med_NTPR_nonA = median(NTPR))


# Make sure all states present
length(unique(TPR_FPR_nonA$state))==8


# Plot overall TPR and FPR ###############################################################################    
# New facet label names for dose variable


setwd(dir=paste0(main_dir, "a_figures"))

overall.PR <- function(detection_rates, mode) {
  #   colors <- factor(c("% Merged" = "lightgreen", "FPR" = "red", "TPR" = "blue"),
  #                    levels=c("red", "blue", "lightgreen"))
  #   TPR_FPR <- dplyr::select(detection_rates, algorithm, threshold, sim, TPR, FPR, perc_dup) %>%
  #     ggplot() +
  #     geom_jitter(aes(x=as.factor(threshold),y=TPR, color="TPR"), shape=0, size=1.5, width = 0.25, height = 1) +
  #     geom_jitter(aes(x=as.factor(threshold),y=FPR, color="FPR"), shape=0, size=1.5, width = 0.25, height = 1) +
  #     geom_jitter(aes(x=as.factor(threshold),y=perc_dup, color="% Merged"), shape=0, size=1.5, width = 0.25, height = 1) +
  #     geom_boxplot(aes(x=as.factor(threshold),y=TPR), color="black", outlier.shape = NA, alpha=0) +
  #     geom_boxplot(aes(x=as.factor(threshold),y=FPR), color="grey", outlier.shape = NA, alpha=0) +
  #     stat_summary(aes(x=as.factor(threshold),y=perc_dup), color="darkgrey", fun = median, fun.min = median, fun.max = median,
  #                  geom = "crossbar", width = 0.75) +
  #     
  #     geom_hline(yintercept=50, linetype="dashed") +
  #     labs(x="Threshold",
  #          y="Rate",
  #          color = "Legend") +
  #     scale_color_manual(values = colors) +
  #     theme_bw() +
  #     theme(legend.title = element_blank(),
  #           panel.grid.major = element_blank(),
  #           panel.grid.minor = element_blank(),
  #           text=element_text(size=12),
  #           axis.text.x = element_text(angle=45, hjust=1)) +
  #     scale_y_continuous(limits=c(0,100), labels = function(x) paste0(sprintf("%.0f", x),"%")) +
  #     facet_grid(algorithm~.)
  #     # stat_summary(aes(x=as.factor(threshold),y=TPR), color="black", fun = median, fun.min = median, fun.max = median,
  #     #              geom = "crossbar", width = 0.75) +
  #     # stat_summary(aes(x=as.factor(threshold),y=FPR), color="grey", fun = median, fun.min = median, fun.max = median,
  #     #              geom = "crossbar", width = 0.5) +
  #     colors <- c("#66a182", "#2e4057")
  
  TPR_FPR <- dplyr::select(detection_rates, algorithm, threshold, sim, TPR, FPR, perc_dup) %>%
    dplyr::group_by(algorithm, threshold, sim, perc_dup) %>%
    dplyr::mutate(TPR.FPR = TPR/FPR) %>%
    dplyr::mutate(TPR.FPR=ifelse(is.nan(TPR.FPR), 0, TPR.FPR)) %>%
    dplyr::group_by(algorithm, threshold) %>%
    dplyr::summarize(med_TPR.FPR=median(TPR.FPR),
                  upper_TPR.FPR = quantile(TPR.FPR, 0.95),
                  lower_TPR.FPR = quantile(TPR.FPR, 0.05),
                  med_perc_dup = median(perc_dup),
                  upper_perc_dup = quantile(perc_dup, 0.95),
                  lower_perc_dup = quantile(perc_dup, 0.05)) %>%
    ggplot() +
    geom_line(aes(x=threshold, y=med_TPR.FPR, group=algorithm, color=algorithm), size=1.5) +
    geom_ribbon(aes(x=threshold, ymin = lower_TPR.FPR, ymax = upper_TPR.FPR, fill = algorithm), alpha=0.5) +
    geom_bar(fill=NA, stat="identity", position="dodge", size=1.5,
             aes(x=threshold, y=med_perc_dup/10, group=algorithm, color=algorithm)) +
    labs(x="Threshold",
         y="TPR/FPR",
         color = "Legend") +
   scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme_bw() +
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          axis.text.x = element_text(angle=45, hjust=1)) +
    scale_y_continuous(breaks=seq(0,10,2),
                       sec.axis = sec_axis(~.*10, name="% Merged", breaks=seq(0,100,10)))

  return(TPR_FPR)
}
TPR_FPR_plot <- overall.PR(TPR_FPR_overall)
TPR_FPR_plot
ggsave(plot=last_plot(), filename="overall_TPR_FPR.pdf", width=6, height=4, units="in")

PR.A <- function(detection_rates, mode) {
#   colors <- factor(c("% Merged" = "lightgreen", "FPR" = "red", "TPR" = "blue"),
#                    levels=c("red", "blue", "lightgreen"))
# TPR_FPR <- dplyr::select(detection_rates, algorithm, threshold, sim, TPR, FPR) %>%
#   ggplot() +
#   # geom_boxplot(aes(x=as.factor(threshold),y=TPR), color="black", outlier.shape = NA, alpha=0) +
#   # geom_boxplot(aes(x=as.factor(threshold),y=FPR), color="grey", outlier.shape = NA, alpha=0) +
#   geom_jitter(aes(x=as.factor(threshold),y=TPR, color="TPR"), shape=0, size=1.5, width = 0.25, height = 1) +
#   geom_jitter(aes(x=as.factor(threshold),y=FPR, color="FPR"), shape=0, size=1.5, width = 0.25, height = 1) +
#   geom_hline(yintercept=50, linetype="dashed") +
#   labs(x="Threshold",
#        y="Rate",
#        color = "Legend") +
#   scale_color_manual(values = colors) +
#   theme_bw() +
#   theme(legend.title = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         text=element_text(size=12),
#         axis.text.x = element_text(angle=45, hjust=1)) +
#   scale_y_continuous(limits=c(0,100), labels = function(x) paste0(sprintf("%.0f", x),"%")) +
#   facet_grid(algorithm~.) +
#   stat_summary(aes(x=as.factor(threshold),y=TPR), color="black", fun = median, fun.min = median, fun.max = median,
#                geom = "crossbar", width = 0.75) +
#   stat_summary(aes(x=as.factor(threshold),y=FPR), color="grey", fun = median, fun.min = median, fun.max = median,
#                geom = "crossbar", width = 0.5)
  colors <- c("#D16103", "#293352")
  
  TPR_FPR <- dplyr::select(detection_rates, algorithm, threshold, sim, TPR, FPR) %>%
    dplyr::group_by(algorithm, threshold, sim) %>%
    dplyr::mutate(TPR.FPR = TPR/FPR) %>%
    dplyr::mutate(TPR.FPR=ifelse(is.nan(TPR.FPR), 0, TPR.FPR)) %>%
    dplyr::group_by(algorithm, threshold) %>%
    dplyr::summarize(med_TPR.FPR=median(TPR.FPR),
                     upper_TPR.FPR = quantile(TPR.FPR, 0.95),
                     lower_TPR.FPR = quantile(TPR.FPR, 0.05)) %>%
    ggplot() +
    geom_line(aes(x=threshold, y=med_TPR.FPR, group=algorithm, color=algorithm), size=1.5) +
    geom_ribbon(aes(x=threshold, ymin = lower_TPR.FPR, ymax = upper_TPR.FPR, fill = algorithm), alpha=0.5) +
    labs(x="Threshold",
         y="TPR/FPR for group A",
         color = "Legend") +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme_bw() +
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          axis.text.x = element_text(angle=45, hjust=1)) +
    scale_y_continuous(breaks=seq(0,10,2))
  

return(TPR_FPR)
}
TPR_FPR_A_plot <- PR.A(TPR_FPR_A)
TPR_FPR_A_plot
ggsave(plot=last_plot(), filename="TPR_FPR_A.pdf", width=6, height=4, units="in")

PR.nonA <- function(detection_rates, mode) {
  colors <- factor(c("% Merged" = "lightgreen", "FPR" = "red", "TPR" = "blue", "NTPR"="purple"),
                   levels=c("red", "blue", "lightgreen", "purple"))
  TPR_FPR <- dplyr::select(detection_rates, algorithm, threshold, sim, state, TPR, NTPR) %>%
    ggplot() +
    # geom_boxplot(aes(x=as.factor(threshold),y=TPR), color="black", outlier.shape = NA, alpha=0) +
    # geom_boxplot(aes(x=as.factor(threshold),y=NTPR), color="grey", outlier.shape = NA, alpha=0) +
    geom_jitter(aes(x=as.factor(threshold),y=TPR, color="TPR"), shape=0, size=1.5, width=0.25, height=1) +
    geom_jitter(aes(x=as.factor(threshold),y=NTPR, color="NTPR"), shape=0, size=1.5, width=0.25, height=1) +
    geom_hline(yintercept=50, linetype="dashed") +
    labs(x="Threshold",
         y="Rate",
         color = "Legend") +
    scale_color_manual(values = colors) +
    theme_bw() +
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          axis.text.x = element_text(angle=45, hjust=1)) +
    scale_y_continuous(limits=c(0,100), labels = function(x) paste0(sprintf("%.0f", x),"%")) +
    facet_grid(algorithm~state) +
    stat_summary(aes(x=as.factor(threshold),y=TPR), color='black', fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", width = 0.75) +
    # stat_summary(aes(x=as.factor(threshold),y=NTPR), color="grey", fun = median, fun.min = median, fun.max = median,
    #              geom = "crossbar", width = 0.5)
    geom_boxplot(aes(x=as.factor(threshold),y=NTPR), color="darkgrey", outlier.shape = NA, alpha=0)
    
    
  
  return(TPR_FPR)
}

TPR_FPR_nonA_plot <- PR.nonA(TPR_FPR_nonA)
TPR_FPR_nonA_plot
ggsave(plot=last_plot(), filename="TPR_FPR_nonA.pdf", width=10, height=4, units="in")



# Plot ############################################################################################
per.state.proportion <- function(performance_data) {
 
  colors <- c("p" = "#FFD524", "a" = "#55185D") 
  prop <- dplyr::select(performance_data, algorithm, threshold, sim, state, proportion, additional) %>%
    dplyr::filter(proportion>0) %>%
    tidyr::drop_na() %>%
    dplyr::group_by(algorithm, threshold, sim, state) %>%
    ggplot() +
    geom_jitter(aes(x=as.factor(threshold),y=additional*100, color="a"), shape=0, size=1.5, width=0.25, height=0.25) +
    geom_jitter(aes(x=as.factor(threshold),y=proportion*100, color="p"), shape=0, size=1.5, width=0.25, height=0.25) +
    geom_hline(yintercept=50, linetype="dashed") +
    labs(x="Threshold",
         y="Proportion of taxa") +
    theme_bw() +
    scale_color_manual(values = colors, 
                       labels=c("% taxa not true",
                                "% true taxa identified")) +
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          axis.text.x = element_text(angle=45, hjust=1)) +
    scale_y_continuous(limits=c(0,100), labels = function(x) paste0(sprintf("%.0f", x),"%")) +
    facet_grid(algorithm~state) +
    stat_summary(aes(x=as.factor(threshold),y=proportion*100), color='black', fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", width = 0.75) +
    stat_summary(aes(x=as.factor(threshold),y=additional*100), color='grey', fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", width = 0.75)
  

  
  return(prop)
}
proportion_plot <- per.state.proportion(performance_data)
proportion_plot
ggsave(plot=last_plot(), filename="proportions2.pdf", width=10, height=4, units="in")

## Subset performance_data to only tue positives
  

#birthPerformance <- function(a) {
  TP_performance <- a[!is.na(a$dynamic),] %>%
    select(algorithm, threshold, sim, state, birth) %>%
    filter(state=="I (R0~5.5)") %>% 
    dplyr::group_by(algorithm, threshold, sim, birth) %>%
    dplyr::summarise(perc=n()) %>% 
    dplyr::group_by(algorithm, threshold, sim) %>% 
    dplyr::mutate(perc=perc/sum(perc)*100) %>%
    dplyr::filter(birth==TRUE) %>% dplyr::select(-birth) %>%
    ggplot() +
    geom_col(aes(x=as.factor(threshold), y=perc)) +
    labs(x="Threshold",
       y=paste0("% correctly identified births")) +
    theme_bw() +
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          axis.text.x = element_text(angle=45, hjust=1)) +
    scale_y_continuous(limits=c(0,100), labels = function(x) paste0(sprintf("%.0f", x),"%")) +
    facet_grid(~algorithm) 
    # stat_summary(aes(x=as.factor(threshold), y=as.numeric(b)), color='black', fun = median, fun.min = median, fun.max = median,
    #              geom = "crossbar", width = 0.75)
  return(TP_performance)
}

# birth_plot <- birthPerformance(performance_data2)
# birth_plot

dynamicPerformance <- function(a) {
  TP_performance <- a[!is.na(a$dynamic_TF),] %>%
    select(algorithm, threshold, state, dynamic_TF) %>%
    dplyr::filter(algorithm=="branchwise") %>%
    dplyr::group_by(threshold, state, dynamic_TF) %>%
    dplyr::summarise(perc=n()) %>% 
    dplyr::group_by(threshold, state) %>% 
    dplyr::mutate(perc=perc/sum(perc)*100)
  i<-1  
  while(i < nrow(TP_performance)) { 
    if(TP_performance$dynamic_TF[i]=="FALSE" & TP_performance$perc[i]==100) {
      new_row <- TP_performance[i,] %>% 
        dplyr::mutate(dynamic_TF=TRUE, perc=0)
      TP_performance <- rbind(TP_performance, new_row)
    }
    i=i+1
  }
  TP_p <- TP_performance %>%
    dplyr::filter(dynamic_TF==TRUE) %>%
    tidyr::drop_na() %>%
    ggplot(aes(x=as.factor(threshold), y=perc, fill=state)) +
    geom_bar(colour = "black", position="dodge", stat="identity") +
    labs(x="Threshold",
         y=paste0("% correctly identified dynamic status")) +
    theme_bw() +
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          axis.text.x = element_text(angle=45, hjust=1)) +
    scale_y_continuous(limits=c(0,100), labels = function(x) paste0(sprintf("%.0f", x),"%"))
  # stat_summary(aes(x=as.factor(threshold), y=as.numeric(b)), color='black', fun = median, fun.min = median, fun.max = median,
  #              geom = "crossbar", width = 0.75)
  return(TP_p)
}

dynamic_plot <- dynamicPerformance(performance_data)
dynamic_plot
ggsave(plot=last_plot(), filename="dynamics.pdf", width=8, height=4, units="in")


priorityPerformance <- function(a) {
  TP_performance <- a[!is.na(a$priority_TF),] %>%
    select(algorithm, threshold, state, priority_TF) %>%
    dplyr::filter(algorithm=="branchwise") %>%
    dplyr::group_by(threshold, state, priority_TF) %>%
    dplyr::summarise(perc=n()) %>% 
    dplyr::group_by(threshold, state) %>% 
    dplyr::mutate(perc=perc/sum(perc)*100)
    i<-1  
    while(i < nrow(TP_performance)) { 
      if(TP_performance$priority_TF[i]=="FALSE" & TP_performance$perc[i]==100) {
        new_row <- TP_performance[i,] %>% mutate(priority_TF=TRUE, perc=0)
        TP_performance <- rbind(TP_performance, new_row)
      }
      i=i+1
    }
    TP_p <- TP_performance %>%
      dplyr::filter(priority_TF==TRUE) %>%
      tidyr::drop_na() %>%
    ggplot(aes(x=as.factor(threshold), y=perc, fill=state)) +
    geom_bar(colour = "black", position="dodge", stat="identity") +
    labs(x="Threshold",
         y=paste0("% correctly identified priority status")) +
    theme_bw() +
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          axis.text.x = element_text(angle=45, hjust=1)) +
    scale_y_continuous(limits=c(0,100), labels = function(x) paste0(sprintf("%.0f", x),"%"))
  # stat_summary(aes(x=as.factor(threshold), y=as.numeric(b)), color='black', fun = median, fun.min = median, fun.max = median,
  #              geom = "crossbar", width = 0.75)
  return(TP_p)
}

priority_plot <- priorityPerformance(performance_data)
priority_plot
ggsave(plot=last_plot(), filename="priority.pdf", width=8, height=4, units="in")

stats <- merge(merge(median_overall, median_A, by=c("algorithm", "threshold")),
               median_nonA, by=c("algorithm", "threshold"))

write.csv(stats, "dyn_stats.csv", quote=F, row.names = F, col.names = T)


