#####################
### Figures_S11.R ###
#####################

suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(stringr) )
suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(tidyr) )
suppressPackageStartupMessages( library(ggplot2) )
suppressPackageStartupMessages( library(cowplot) )
suppressPackageStartupMessages( library(corrplot) )
suppressPackageStartupMessages( library(ppcor) )
suppressPackageStartupMessages( library(RColorBrewer) )
suppressPackageStartupMessages( library(ggpubr) )
suppressPackageStartupMessages( library(extrafonts) )

colPal <- list(Blue = "#003366", 
               Red = "#E31B23", 
               MediumBlue = "#005CAB", 
               LightBlue = "#DCEEF3", 
               AccentYellow = "#FFC325", 
               Grey = hsv(0,0,0.5),
               CoolGrey = "#E6F1EE")

theme_set(theme_classic(base_family = "Arial"))

### LOAD DATA ##########

splicingFactors <- fread("lib/splicing_factors.txt") %>% 
  dplyr::filter(str_detect(geneSymbol, "MBNL|CELF"))

log2CPM_data <- fread("data/log2CPM_data_filtered_annotated.csv") %>% 
  dplyr::filter(geneID %in% splicingFactors$geneID)

group_levels <- c("Male-Prenatal", "Female-Prenatal", "Male-Postnatal", "Female-Postnatal", "Male-DM1", "Female-DM1", "Male-Unaffected", "Female-Unaffected")

BrainSpan_log2CPM <- log2CPM_data %>%
  dplyr::filter(dataset == "BrainSpan") %>%
  dplyr::mutate(
    group = paste(sex, group, sep = "-"),
    group = factor(group, levels = group_levels),    
    age = age_in_days)

Otero_log2CPM <- log2CPM_data %>%
  dplyr::filter(dataset == "Otero et al. (2021)") %>%
  dplyr::mutate(
    group = paste(sex, group, sep = "-"),
    group = factor(group, levels = group_levels)
  )

### SPLICING FACTORS OVER AGE ###########

agePlot <- function(data, y_max = 0){
  data %>%
    ggplot(aes(y = log2CPM, x = log(age))) +
    geom_smooth(method = "auto", formula = y ~ x, se = FALSE, aes(color = sex)) +
    ylim(0, y_max) +
    theme_classic() +
    theme(text = element_text(size = 12), 
          strip.text.x = element_blank(),
          strip.background =element_rect(color="white", fill="white"),
          axis.text.x = element_text(angle = -60, hjust = 0),
          legend.key.size = unit(1.5, "cm")) +
    labs(color = "") +
    ylab("log2(CPM)") + xlab("log10(Age)") + 
    scale_color_manual(values = c(colPal$MediumBlue, colPal$Grey)) +
    scale_x_continuous(breaks = c(log(90), log(280), log(3 * 365 + 280), 
                                  log(10 * 365 + 280), log(20 * 365 + 280)), 
                       labels = c("3 Mos", "Birth", "3 Yrs", "10 Yrs", "20 Yrs")) +
    facet_wrap(~ geneSymbol, ncol = 4, scales = "free_y")
}

ageplots <- list()

for(gene in c("CELF1", "MBNL1", "MBNL2")){
  
  p <- BrainSpan_log2CPM %>%
    filter(geneSymbol == gene) %>%
    agePlot(y_max = 10)
  ageplots <- append(ageplots, list(p))
}

### MORE STATISTICS ######

df <- bind_rows(BrainSpan_log2CPM, Otero_log2CPM)

stat_results <- compare_means(log2CPM ~ group, df, group.by = c("geneSymbol", "dataset"), p.adjust.method = "fdr")
stat_results$p.adj.sig <- gtools::stars.pval(stat_results$p.adj)
stat_results$p.adj.sig[!str_detect(stat_results$p.adj.sig, "\\*")] <- "ns"

# filter for comparisons of interest
stat_results <- dplyr::filter(stat_results, 
                              group1 == "Male-Prenatal" & group2 == "Female-Prenatal" |
                                group1 == "Male-Postnatal" & group2 == "Female-Postnatal" |
                                group1 == "Male-DM1" & group2 == "Female-DM1" |
                                group1 == "Male-Unaffected" & group2 == "Female-Unaffected")

### BOXPLOTS #####

# function to plot boxplots
exonBoxplot <- function(data, title = "", colors = NULL, dataset_ = NULL, sig_y_pos = 300){
  data %>%
    ggboxplot(y = "log2CPM", x = "group", fill = "group") +
    xlab("") + ylab("log2(CPM)") +
    ggtitle(title) +
    scale_y_continuous(breaks = c(0, 2.5, 5, 7.5, 10), limits = c(0, 10)) +
    scale_fill_manual(values = colors) +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      legend.title = element_blank(),
      legend.key.size = unit(1, "cm"),
      legend.text = element_text(size = 10),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      strip.background = element_rect(color="white", fill="white"),
      strip.placement = "outside",
      panel.spacing.y = unit(2, "lines")) +
    facet_wrap(~ geneSymbol, ncol = 4, strip.position = "top", scales = "free_x") +
    stat_pvalue_manual(
      filter(stat_results, 
             dataset == dataset_,
             geneSymbol %in% data$geneSymbol), 
      label = "p.adj.sig", y.position = sig_y_pos)
}

df <- bind_rows(BrainSpan_log2CPM, Otero_log2CPM)

plots <- data.frame(
  dataset = rep(c("Otero et al. (2021)", "BrainSpan"), 3), 
  gene = rep(c("CELF1", "MBNL1", "MBNL2"), each = 2),
  ylim = log2(rep(c(470, 350, 450), each = 2)),
  sig_pos = c(9, 9.5, 9, 7.5, 9.5, 9)
)

boxplots <- list()

for(i in 1:nrow(plots)){
  palette <- switch(as.character(plots$dataset[i]),
                          "Otero et al. (2021)" = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"),
                          "BrainSpan" = c("#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
  )
  
  p <- df %>%
    filter(
      dataset == plots$dataset[i],
      geneSymbol == plots$gene[i]
    ) %>%
    exonBoxplot(
      dataset = plots$dataset[i],
      colors = palette,
      sig_y_pos = plots$sig_pos[i]
    )
  
  boxplots <- append(boxplots, list(p))
}

### ARRANGE FIGURE ###########

l_Otero <- get_legend(boxplots[[1]])
l_BrainSpan <- get_legend(boxplots[[2]])
l_age <- get_legend(ageplots[[1]])

plot_grid(
  plot_grid(
    boxplots[[1]] + theme(legend.position = "none"),
    boxplots[[3]] + theme(legend.position = "none"),
    boxplots[[5]] + theme(legend.position = "none"),
    l_Otero, ncol = 4
  ), plot_grid(
    boxplots[[2]] + theme(legend.position = "none"),
    boxplots[[4]] + theme(legend.position = "none"),
    boxplots[[6]] + theme(legend.position = "none"),
    l_BrainSpan, ncol = 4
  ), plot_grid(
    ageplots[[1]] + theme(legend.position = "none"),
    ageplots[[2]] + theme(legend.position = "none"),
    ageplots[[3]] + theme(legend.position = "none"),
    l_age, ncol = 4
  ),
  ncol = 1, labels = "", label_size = 20, rel_heights = c(0.325,0.325,0.35)
)

ggsave2(filename = "results/Figure_S11.pdf", width = 8, height = 8, dpi = 300)

### COMPUTE LOGFC FOR MBNL1, MBNL2 AND CELF1 EXPRESSION ###

Otero_logFC <- Otero_log2CPM %>%
  dplyr::filter(geneSymbol %in% c("MBNL1", "MBNL2", "MBNL3", "CELF1")) %>%
  dplyr::group_by(geneSymbol, group) %>%
  dplyr::summarise(mean = mean(log2CPM)) %>%
  tidyr::spread(key = "group", value = "mean") %>%
  dplyr::mutate(
    logFC_DM1 = log2(2^`Male-DM1`/2^`Female-DM1`),
    logFC_Unaffected = log2(2^`Male-Unaffected`/2^`Female-Unaffected`))

BrainSpan_logFC <- BrainSpan_log2CPM %>%
  dplyr::filter(geneSymbol %in% c("MBNL1", "MBNL2", "MBNL3","CELF1")) %>%
  dplyr::group_by(geneSymbol, group) %>%
  dplyr::summarise(mean = mean(log2CPM)) %>%
  tidyr::spread(key = "group", value = "mean") %>%
  dplyr::mutate(
    logFC_Prenatal = log2(2^`Male-Prenatal`/2^`Female-Prenatal`),
    logFC_Postnatal = log2(2^`Male-Postnatal`/2^`Female-Postnatal`))

Otero_logFC[, c("geneSymbol", "logFC_DM1", "logFC_Unaffected")]
BrainSpan_logFC[, c("geneSymbol", "logFC_Prenatal", "logFC_Postnatal")]
