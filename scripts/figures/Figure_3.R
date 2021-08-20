###################
### Figures_3.R ###
###################

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

# save table with correlation results
fwrite(corrResults, "data/BrainSpan_correlation_splicingFactors_age.txt")

# create table for annotation of correlation coefficient
anno <- data.frame(
  geneSymbol = corrResults$geneSymbol,
  label = paste0("R = ", sprintf("%0.3f", round(corrResults$estimate, digits = 3))),
  x = log(10 * 365 + 280),
  y = rep(2.5, 9)
  #y = c(rep(400,6),rep(250,3))
  #y = sapply(c(250, rep(400,5), rep(250,3)), log)
)

agePlot <- function(data, y_max = 0){
  data %>%
    ggplot(aes(y = log2CPM, x = log(age))) +
    geom_point(aes(color = group), size = 0.5) + 
    geom_smooth(method = "auto", formula = y ~ x, se = FALSE, color= "Black") +
    ylim(0, y_max) +
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

### BOXPLOTS #####

# function to plot boxplots
exonBoxplot <- function(data, title = "", colors = NULL, dataset_ = NULL, sig_y_pos = 300, y_max){
  data %>%
    ggboxplot(y = "log2CPM", x = "group", fill = "group") +
    #ylim(2.5, 10.0) +
    xlab("") + ylab("log2(CPM)") +
    ggtitle(title) +
    scale_fill_manual(values = colors ) +
    scale_y_continuous(breaks = c(0, 2.5, 5, 7.5, 10), limits = c(0, 10)) +
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
  #sig_pos = c(350, 460, 330, 150, 400, 300)
  sig_pos = c(9, 9.5, 9, 7.5, 9.5, 9)
  
)

boxplots <- list()

for(i in 1:nrow(plots)){
  
  color <- if (plots$dataset[i] == "Otero et al. (2021)") c(colPal$Red, colPal$Grey) else c(colPal$MediumBlue, colPal$Grey)
  
  p <- df %>%
    filter(
      dataset == plots$dataset[i],
      geneSymbol == plots$gene[i]
    ) %>%
    exonBoxplot(
      colors = color,
      dataset = plots$dataset[i],
      y_max = plots$ylim[i],
      sig_y_pos = plots$sig_pos[i]
    )
  
  boxplots <- append(boxplots, list(p))
}

### ARRANGE FIGURE ###########

l_Otero <- get_legend(boxplots[[1]])
l_BrainSpan <- get_legend(boxplots[[2]])

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
    ageplots[[1]],
    ageplots[[2]],
    ageplots[[3]],
    NULL, ncol = 4
  ),
  ncol = 1, labels = "", label_size = 20, rel_heights = c(0.325,0.325,0.35)
)

ggsave2(filename = "results/Figure_3.pdf", width = 8, height = 8, dpi = 300)


### COMPUTE LOGFC FOR MBNL1, MBNL2 AND CELF1 EXPRESSION ###

Otero_logFC <- Otero_log2CPM %>%
  dplyr::filter(geneSymbol %in% c("MBNL1", "MBNL2", "MBNL3", "CELF1")) %>%
  dplyr::group_by(geneSymbol, group) %>%
  dplyr::summarise(mean = mean(log2CPM)) %>%
  tidyr::spread(key = "group", value = "mean") %>%
  dplyr::mutate(logFC = log2(2^Unaffected/2^DM1))

BrainSpan_logFC <- BrainSpan_log2CPM %>%
  dplyr::filter(geneSymbol %in% c("MBNL1", "MBNL2", "MBNL3","CELF1")) %>%
  dplyr::group_by(geneSymbol, group) %>%
  dplyr::summarise(mean = mean(log2CPM)) %>%
  tidyr::spread(key = "group", value = "mean") %>%
  dplyr::mutate(logFC = log2(2^Postnatal/2^Prenatal))

Otero_logFC
BrainSpan_logFC

