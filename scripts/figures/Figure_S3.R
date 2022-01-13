###################
### Figure_S3.R ###
###################

suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(stringr) )
suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(tidyr) )
suppressPackageStartupMessages( library(ggplot2) )
suppressPackageStartupMessages( library(ggpubr) )
suppressPackageStartupMessages( library(cowplot) )
suppressPackageStartupMessages( library(RColorBrewer) )
suppressPackageStartupMessages( library(extrafont) )

colPal <- list(Blue = "#003366", 
               Red = "#E31B23", 
               MediumBlue = "#005CAB", 
               LightBlue = "#DCEEF3", 
               AccentYellow = "#FFC325", 
               Grey = hsv(0,0,0.5),
               CoolGrey = "#E6F1EE")

theme_set(theme_classic(base_family = "Arial"))

### LOAD DATA ############

high_conf_events <- fread("data/high_confidence_events.csv")
psi_data_filtered <- fread("data/psi_data_filtered_annotated.csv") 

group_levels <- c("Male-Prenatal", "Female-Prenatal", "Male-Postnatal", "Female-Postnatal", "Male-DM1", "Female-DM1", "Male-Unaffected", "Female-Unaffected")

BrainSpan_psi <- psi_data_filtered %>% 
  dplyr::filter(dataset == "BrainSpan") %>%
  dplyr::mutate(
    group = paste(sex, group, sep = "-"),
    group = factor(group, levels = group_levels),
    age = age_in_days)

Otero_psi <- psi_data_filtered %>%
  dplyr::filter(dataset == "Otero et al. (2021)") %>%
  dplyr::mutate(
    group = paste(sex, group, sep = "-"),
    group = factor(group, levels = group_levels))

### STATISTICS #####

# compute correlation between exon inclusion and donor age
corrResults <- data.table(matrix(ncol=7,nrow=0))
for(sex_ in unique(BrainSpan_psi$sex)) {
  for(event_ in high_conf_events$event_id){
    data <- filter(BrainSpan_psi, eventID == event_, sex == sex_)
    results <- cor.test(data$exonInclusion, data$age, method="spearman", exact=F)
    corrResults <- rbind(corrResults, data.table(event_, 
                                                 sex_,
                                                 results$statistic,
                                                 results$p.value,
                                                 NA, 
                                                 results$estimate,
                                                 results$method), use.names=F)
  }
}

colnames(corrResults) <- c("eventID", "sex", "statistic", "p.value", "adj.p.value", "estimate", "method")
corrResults$adj.p.value <- p.adjust(corrResults$p.value, method = "fdr")

# get highlighted events
high_conf_events <- arrange(high_conf_events, dPSI_DM1_CTRL)
highlight_events <- rbind(head(high_conf_events, 4), tail(high_conf_events, 4))

# define order of events based on DM1 inclusion and event name
order <- highlight_events$SE_event_name[ order(highlight_events$DM1_inclusion, highlight_events$SE_event_name, decreasing = F)]
highlight_events <- highlight_events[match(order, highlight_events$SE_event_name),]

# filter for highlight events
BrainSpan_psi <- dplyr::filter(BrainSpan_psi, SE_event_name %in% highlight_events$SE_event_name)
Otero_psi <- dplyr::filter(Otero_psi, SE_event_name %in% highlight_events$SE_event_name)

# order events in BrainSpan data
BrainSpan_psi$SE_event_name <- factor(BrainSpan_psi$SE_event_name, levels = order)

# create table for annotation of correlation coefficient
anno <- data.frame(
  eventID = highlight_events$event_id,
  SE_event_name = highlight_events$SE_event_name,
  sex = c(rep("Male", 8), rep("Female", 8))) %>% 
  merge(corrResults[, c("eventID", "sex", "estimate")]) %>%
  dplyr::mutate(label = paste0("R = ", sprintf("%0.3f", round(estimate, digits = 3)))) %>% 
  dplyr::select(-eventID, -estimate)

### SPLICING OVER AGE #####

exonAgePlot <- function(data){
  data %>%
    ggplot(aes(y = exonInclusion, x = log(age))) +
    geom_smooth(method = "auto", se = FALSE, aes(color = sex)) +
    theme_classic() +
    theme(text = element_text(size = 20), 
          strip.text.x = element_blank(),
          strip.background =element_rect(color="white", fill="white"),
          axis.text.x = element_text(angle = -60, hjust = 0),
          axis.title.y = element_text(angle = 0, vjust = 0.5),
          legend.key.size = unit(1.5, "cm")) +
    labs(color = "") +
    ylab(expression(psi)) + xlab("log10(Age)") + 
    scale_color_manual(values = c(colPal$MediumBlue, colPal$Grey)) +
    scale_x_continuous(breaks = c(log(90), log(280), log(3 * 365 + 280), 
                                  log(10 * 365 + 280), log(20 * 365 + 280)), 
                       labels = c("3 Mos", "Birth", "3 Yrs", "10 Yrs", "20 Yrs")) +
    facet_wrap(~ SE_event_name, ncol = 4, scales = "fixed")
}

BrainSpan.exon.age.increased <- BrainSpan_psi %>%
  filter(SE_event_name %in% highlight_events$SE_event_name[highlight_events$DM1_inclusion == "+"]) %>%
  exonAgePlot()

BrainSpan.exon.age.decreased <- BrainSpan_psi %>%
  filter(SE_event_name %in% highlight_events$SE_event_name[highlight_events$DM1_inclusion == "-"]) %>%
  exonAgePlot()

#### LOAD OTERO DATA ####

# order events in Otero data
Otero_psi$SE_event_name <- factor(Otero_psi$SE_event_name, levels = order)

### MORE STATISTICS ####

df <- bind_rows(BrainSpan_psi, Otero_psi)
stat_results <- compare_means(exonInclusion ~ group, df, group.by = c("SE_event_name", "dataset"), p.adjust.method = "fdr")
stat_results$p.adj.sig <- gtools::stars.pval(stat_results$p.adj)
stat_results$p.adj.sig[!str_detect(stat_results$p.adj.sig, "\\*")] <- "ns"

# filter for comparisons of interest
stat_results <- dplyr::filter(stat_results, 
                              group1 == "Male-Prenatal" & group2 == "Female-Prenatal" |
                                group1 == "Male-Postnatal" & group2 == "Female-Postnatal" |
                                group1 == "Male-DM1" & group2 == "Female-DM1" |
                                group1 == "Male-Unaffected" & group2 == "Female-Unaffected")

### PSI BOXPLOTS ##############

# function to plot boxplots
exonBoxplot <- function(data, title = "", colors = NULL, dataset_ = NULL){
  data %>%
    ggboxplot(y = "exonInclusion", x = "group", fill = "group") +
    xlab("") + ylab(expression(psi)) +
    ggtitle(title) +
    scale_y_continuous(limits=c(0.0,1.2), breaks = seq(0,1,by=0.25)) +
    scale_fill_manual(values = colors) +
    theme_classic() +
    theme(
      text = element_text(size = 20),
      legend.title = element_blank(),
      legend.key.size = unit(1, "cm"),
      legend.text = element_text(size = 18),
      axis.title.y = element_text(angle = 0, vjust = 0.5),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      strip.background = element_rect(color="white", fill="white"),
      strip.placement = "outside",
      panel.spacing.y = unit(2, "lines")) +
    facet_wrap(~ SE_event_name, ncol = 4, strip.position = "top", scales = "free_x") +
    stat_pvalue_manual(
      filter(stat_results,
        dataset == dataset_,
        SE_event_name %in% data$SE_event_name),
      label = "p.adj.sig", y.position = 1.1)
}

#---BRAINSPAN---#
BrainSpan.exon.bp.increased <- BrainSpan_psi %>%
  filter(SE_event_name %in% highlight_events$SE_event_name[highlight_events$DM1_inclusion == "+"]) %>%
  exonBoxplot(
    colors = c("#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"), 
    dataset = "BrainSpan"
  )

BrainSpan.exon.bp.decreased <- BrainSpan_psi %>%
  filter(SE_event_name %in% highlight_events$SE_event_name[highlight_events$DM1_inclusion == "-"]) %>%
  exonBoxplot(
    colors = c("#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"), 
    dataset = "BrainSpan"
  )

#---Otero---#
Otero.exon.bp.increased <- Otero_psi %>%
  filter(SE_event_name %in% highlight_events$SE_event_name[highlight_events$DM1_inclusion == "+"]) %>%
  exonBoxplot(
    colors = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"), 
    dataset = "Otero et al. (2021)"
    )

Otero.exon.bp.decreased <- Otero_psi %>% 
  filter(SE_event_name %in% highlight_events$SE_event_name[highlight_events$DM1_inclusion == "-"]) %>%
  exonBoxplot(
    colors = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"),
    dataset = "Otero et al. (2021)"
  )

### ARRANGE FIGURE ###########

l_Otero <- get_legend(Otero.exon.bp.increased)
l_BrainSpan <- get_legend(BrainSpan.exon.bp.increased)
l_age <- get_legend(BrainSpan.exon.age.decreased)

plot_grid(
  plot_grid(Otero.exon.bp.decreased + theme(legend.position = "none"), l_Otero, 
            BrainSpan.exon.bp.decreased + theme(legend.position = "none"), l_BrainSpan, 
            BrainSpan.exon.age.decreased + theme(legend.position = "none"), l_age, 
            rel_widths = c(0.8,0.2), ncol = 2),
  plot_grid(Otero.exon.bp.increased + theme(legend.position = "none"), l_Otero, 
            BrainSpan.exon.bp.increased + theme(legend.position = "none"), l_BrainSpan, 
            BrainSpan.exon.age.increased + theme(legend.position = "none"), l_age, 
            rel_widths = c(0.8,0.2), ncol = 2),
  ncol = 1, labels = "AUTO", label_size = 20, scale = c(0.95,0.95))

ggsave2(filename = "results/Figure_S3.pdf", width = 12, height = 14, dpi = 300)
