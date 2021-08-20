##################
### Figure_2.R ###
##################

suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(stringr) )
suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(tidyr) )
suppressPackageStartupMessages( library(ggplot2) )
suppressPackageStartupMessages( library(ggpubr) )
suppressPackageStartupMessages( library(cowplot) )
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

BrainSpan_psi <- psi_data_filtered %>% 
  dplyr::filter(dataset == "BrainSpan") %>%
  dplyr::mutate(
    group = factor(group, levels = c("Prenatal", "Postnatal")),
    age = age_in_days)

Otero_psi <- psi_data_filtered %>%
  dplyr::filter(dataset == "Otero et al. (2021)") %>%
  dplyr::mutate(group = factor(group, levels = c("DM1", "Unaffected")))

### STATISTICS #####

# compute correlation between exon inclusion and donor age
corrResults <- data.table(matrix(ncol=6,nrow=0))
for(event_ in high_conf_events$event_id){
  data <- filter(BrainSpan_psi, eventID == event_)
  results <- cor.test(data$exonInclusion, data$age, method="spearman", exact=F)
  corrResults <- rbind(corrResults, data.table(event_, 
                                               results$statistic,
                                               results$p.value,
                                               NA, 
                                               results$estimate,
                                               results$method), use.names=F)
}
colnames(corrResults) <- c("eventID", "statistic", "p.value", "adj.p.value", "estimate", "method")
corrResults$adj.p.value <- p.adjust(corrResults$p.value, method = "fdr")

# save table with correlation results (required for script to create Table S2/S3/S4)
fwrite(corrResults, "data/BrainSpan_correlation_events_age.txt")

# get highlighted events
high_conf_events <- arrange(high_conf_events, dPSI_DM1_CTRL)
highlight_events <- rbind(head(high_conf_events, 4), tail(high_conf_events, 4))

# define order of events based on DM1 inclusion and event name
order <- highlight_events$SE_event_name[ order(highlight_events$DM1_inclusion, highlight_events$SE_event_name, decreasing = F)]
highlight_events <- highlight_events[match(order, highlight_events$SE_event_name),]

# order events in BrainSpan data
BrainSpan_psi$SE_event_name <- factor(BrainSpan_psi$SE_event_name, levels = order)

# create table for annotation of correlation coefficient
anno <- data.frame(
  SE_event_name = highlight_events$SE_event_name,
  label = paste0("R = ",
                 sprintf("%0.3f",
                         round(corrResults$estimate[match(highlight_events$event_id, corrResults$eventID)],
                               digits = 3)))
)

anno <- anno %>%
  mutate(x = c(rep(log(280), 3), log(18 * 365 + 280),
               rep(log(12 * 365 + 280),4)),
         y = c(0.8,0.8,0.8,0.25,
               0.95,0.95,0.95,0.95))

### SPLICING OVER AGE #####

exonAgePlot <- function(data){
  data %>%
    ggplot(aes(y = exonInclusion, x = log(age))) +
    geom_point(aes(color = group), size = 0.5) + 
    geom_smooth(method = "auto", se = FALSE, color= "Black") +
    theme_classic() +
    theme(text = element_text(size = 20), 
          strip.text.x = element_blank(),
          strip.background =element_rect(color="white", fill="white"),
          #strip.placement = "outside",
          legend.position = "none",
          axis.text.x = element_text(angle = -60, hjust = 0),
          axis.title.y = element_text(angle = 0, vjust = 0.5)
          #panel.spacing.y = unit(2, "lines")
    ) +
    ylab(expression(psi)) + xlab("log10(Age)") + 
    scale_color_manual(values = c(colPal$MediumBlue, colPal$Grey)) +
    scale_x_continuous(breaks = c(log(90), log(280), log(3 * 365 + 280), 
                                  log(10 * 365 + 280), log(20 * 365 + 280)), 
                       labels = c("3 Mos", "Birth", "3 Yrs", "10 Yrs", "20 Yrs")) +
    geom_text(data = filter(anno, SE_event_name %in% data$SE_event_name),
              aes(x = x, y = y, label = label), size = 4) +
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

### PSI BOXPLOTS ##############

# function to plot boxplots
exonBoxplot <- function(data, title = "", colors = NULL, dataset_ = NULL){
  data %>%
    ggboxplot(y = "exonInclusion", x = "group", fill = "group") +
    xlab("") + ylab(expression(psi)) +
    ggtitle(title) +
    scale_y_continuous(limits=c(0.0,1.2), breaks = seq(0,1,by=0.25)) +
    scale_fill_manual(values = colors ) +
    #scale_fill_grey(start = 0.5) +
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
    colors = c(colPal$MediumBlue, colPal$Grey), 
    dataset = "BrainSpan"
  )

BrainSpan.exon.bp.decreased <- BrainSpan_psi %>%
  filter(SE_event_name %in% highlight_events$SE_event_name[highlight_events$DM1_inclusion == "-"]) %>%
  exonBoxplot(
    colors = c(colPal$MediumBlue, colPal$Grey), 
    dataset = "BrainSpan"
  )

#---Otero---#
Otero.exon.bp.increased <- Otero_psi %>%
  filter(SE_event_name %in% highlight_events$SE_event_name[highlight_events$DM1_inclusion == "+"]) %>%
  exonBoxplot(
    colors = c(colPal$Red, colPal$Grey), 
    dataset = "Otero et al. (2021)"
    )

Otero.exon.bp.decreased <- Otero_psi %>% 
  filter(SE_event_name %in% highlight_events$SE_event_name[highlight_events$DM1_inclusion == "-"]) %>%
  exonBoxplot(
    colors = c(colPal$Red, colPal$Grey),
    dataset = "Otero et al. (2021)"
  )

### ARRANGE FIGURE ###########

l_Otero <- get_legend(Otero.exon.bp.increased)
l_BrainSpan <- get_legend(BrainSpan.exon.bp.increased)

plot_grid(
  plot_grid(Otero.exon.bp.increased + theme(legend.position = "none"), l_Otero, 
            BrainSpan.exon.bp.increased + theme(
              legend.position = "none" #, 
              #strip.text = element_blank()
              ), l_BrainSpan, 
            BrainSpan.exon.age.increased, NULL, 
            rel_widths = c(0.8,0.2), ncol = 2),
  plot_grid(Otero.exon.bp.decreased + theme(legend.position = "none"), l_Otero, 
            BrainSpan.exon.bp.decreased + theme(
              legend.position = "none" #, 
              #strip.text = element_blank()
              ), l_BrainSpan, 
            BrainSpan.exon.age.decreased, NULL, 
            rel_widths = c(0.8,0.2), ncol = 2),
  ncol = 1, labels = "AUTO", label_size = 20, scale = c(0.95,0.95))

ggsave2(filename = "results/Figure_2.pdf", width = 12, height = 14, dpi = 300)
