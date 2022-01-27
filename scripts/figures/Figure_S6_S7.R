#######################
### Figures_S6_S7.R ###
#######################

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

BrainSpan_log2CPM <- log2CPM_data %>%
  dplyr::filter(dataset == "BrainSpan") %>%
  dplyr::mutate(
    group = factor(group, levels = c("Prenatal", "Postnatal")),
    age = age_in_days)

Otero_log2CPM <- log2CPM_data %>%
  dplyr::filter(dataset == "Otero et al. (2021)") %>%
  dplyr::mutate(group = factor(group, levels = c("DM1", "Unaffected")))


### SPLICING FACTORS OVER AGE ###########

# compute correlation between splicing factor expression and donor age
corrResults <- data.table(matrix(ncol=6,nrow=0))
for(gene_ in splicingFactors$geneSymbol){
  data <- filter(BrainSpan_log2CPM, geneSymbol == gene_)
  results <- cor.test(data$log2CPM, data$age, method="spearman", exact=F)
  corrResults <- rbind(corrResults, data.table(gene_, 
                                               results$statistic,
                                               results$p.value,
                                               NA, 
                                               results$estimate,
                                               results$method), use.names=F)
}
colnames(corrResults) <- c("geneSymbol", "statistic", "p.value", "adj.p.value", "estimate", "method")
corrResults$adj.p.value <- p.adjust(corrResults$p.value, method = "fdr")

# create table for annotation of correlation coefficient
anno <- data.frame(
  geneSymbol = corrResults$geneSymbol,
  label = paste0("R = ", sprintf("%0.3f", round(corrResults$estimate, digits = 3))),
  x = log(10 * 365 + 280),
  y = c(rep(2, 8), 5)
  #y = c(rep(400,6),rep(250,3))
  #y = sapply(c(250, rep(400,5), rep(250,3)), log)
)

agePlot <- function(data, y_max = 0, ncol){
  data %>%
    ggplot(aes(y = log2CPM, x = log(age))) +
    geom_point(aes(color = group), size = 0.5) + 
    geom_smooth(method = "auto", formula = y ~ x, se = FALSE, color= "Black") +
    ylim(-1, y_max) +
    theme_classic() +
    theme(text = element_text(size = 12), 
          strip.text.x = element_blank(),
          strip.background =element_rect(color="white", fill="white"),
          #strip.placement = "outside",
          legend.position = "none",
          axis.text.x = element_text(angle = -60, hjust = 0)
          #panel.spacing.y = unit(2, "lines")
    ) +
    ylab("log2(CPM)") + xlab("log10(Age)") + 
    scale_color_manual(values = c(colPal$MediumBlue, colPal$Grey)) +
    scale_x_continuous(breaks = c(log(90), log(280), log(3 * 365 + 280), 
                                  log(10 * 365 + 280), log(20 * 365 + 280)), 
                       labels = c("3 Mos", "Birth", "3 Yrs", "10 Yrs", "20 Yrs")) +
    geom_text(data = filter(anno, geneSymbol %in% data$geneSymbol),
              aes(x = x, y = y, label = label), size = 3) +
    facet_wrap(~ geneSymbol, ncol = ncol, scales = "fixed")
}

CELF.ageplots <- BrainSpan_log2CPM %>%
  filter(str_detect(geneSymbol, "CELF")) %>%
  agePlot(y_max = 10, ncol = 6)

MBNL.ageplots <- BrainSpan_log2CPM %>%
  filter(str_detect(geneSymbol, "MBNL")) %>%
  agePlot(y_max = 10, ncol = 3)

### MORE STATISTICS ######

df <- bind_rows(BrainSpan_log2CPM, Otero_log2CPM)

stat_results <- compare_means(log2CPM ~ group, df, group.by = c("geneSymbol", "dataset"), p.adjust.method = "fdr") %>% arrange(geneSymbol)
stat_results$p.adj.sig <- gtools::stars.pval(stat_results$p.adj)
stat_results$p.adj.sig[!str_detect(stat_results$p.adj.sig, "\\*")] <- "ns"
stat_results$y.pos <- c(9.5, 9, 10, 10, 9.5, 8, 8.5, 8.5, 9, 8.5, 6.5, 6.5, 8, 9, 9, 10, 5.5, 5)
  
### BOXPLOTS #####

# function to plot boxplots
makeBoxplot <- function(data, title = "", colors = NULL, stat_results = NULL, ncol){
  data %>%
    ggplot() + 
    geom_boxplot(mapping = aes(y = log2CPM, x = group, fill = group), fatten = 1.5) +
    #ylim(2.5, 10.0) +
    xlab("") + ylab("log2(CPM)") +
    ggtitle(title) +
    scale_fill_manual(values = colors ) +
    scale_y_continuous(breaks = c(0, 2.5, 5, 7.5, 10), limits = c(-1, 11)) +
    #scale_fill_grey(start = 0.5) +
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
    facet_wrap(~ geneSymbol, ncol = ncol, strip.position = "top", scales = "free_x") +
    stat_pvalue_manual(stat_results, label = "p.adj.sig", y.position = "y.pos")
}

CELF.BrainSpan.boxplots <- BrainSpan_log2CPM %>%
  filter(str_detect(geneSymbol, "CELF")) %>%
  makeBoxplot(colors = c(colPal$MediumBlue, colPal$Grey),
              ncol = 6,
              stat_results = filter(stat_results, dataset == "BrainSpan", str_detect(geneSymbol, "CELF"))
  )

CELF.Otero.boxplots <- Otero_log2CPM %>%
  filter(str_detect(geneSymbol, "CELF")) %>%
  makeBoxplot(colors = c(colPal$Red, colPal$Grey),
              ncol = 6,
              stat_results = filter(stat_results, dataset == "Otero et al. (2021)", str_detect(geneSymbol, "CELF"))
  )

MBNL.BrainSpan.boxplots <- BrainSpan_log2CPM %>%
  filter(str_detect(geneSymbol, "MBNL")) %>%
  makeBoxplot(colors = c(colPal$MediumBlue, colPal$Grey),
              ncol = 3,
              stat_results = filter(stat_results, dataset == "BrainSpan", str_detect(geneSymbol, "MBNL"))
  )

MBNL.Otero.boxplots <- Otero_log2CPM %>%
  filter(str_detect(geneSymbol, "MBNL")) %>%
  makeBoxplot(colors = c(colPal$Red, colPal$Grey),
              ncol = 3,
              stat_results = filter(stat_results, dataset == "Otero et al. (2021)", str_detect(geneSymbol, "MBNL"))
  )

### ARRANGE FIGURE ###########

l_Otero <- get_legend(MBNL.Otero.boxplots)
l_BrainSpan <- get_legend(MBNL.BrainSpan.boxplots)

plot_grid(
  plot_grid(
    CELF.Otero.boxplots + theme(legend.position = "none"),
    CELF.BrainSpan.boxplots + theme(legend.position = "none"),
    CELF.ageplots + theme(legend.position = "none"),
    ncol = 1
  ),
  plot_grid(l_Otero, l_BrainSpan, ncol = 1),
  ncol = 2, labels = "", label_size = 20, rel_widths = c(0.8, 0.2), rel_heights = c(0.325,0.325,0.35)
)

ggsave2(filename = "results/Figure_S6.pdf", width = 10, height = 6, dpi = 300)

plot_grid(
  plot_grid(
    MBNL.Otero.boxplots + theme(legend.position = "none"),
    MBNL.BrainSpan.boxplots + theme(legend.position = "none"),
    MBNL.ageplots + theme(legend.position = "none"),
    ncol = 1
  ),
  plot_grid(l_Otero, l_BrainSpan, ncol = 1),
  ncol = 2, labels = "", label_size = 20, rel_widths = c(0.8, 0.2), rel_heights = c(0.325,0.325,0.35)
)

ggsave2(filename = "results/Figure_S7.pdf", width = 6, height = 6, dpi = 300)

