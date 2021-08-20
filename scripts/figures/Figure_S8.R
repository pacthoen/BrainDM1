###################
### Figure_S8.R ###
###################

suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(stringr) )
suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(tidyr) )
suppressPackageStartupMessages( library(ggplot2) )
suppressPackageStartupMessages( library(cowplot) )
suppressPackageStartupMessages( library(ggcorrplot) )
suppressPackageStartupMessages( library(ppcor) )
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

# load external functions
source("scripts/partialCorr.R")

### LOAD DATA ###########

high_conf_events <- fread("data/high_confidence_events.csv")

splicingFactors <- fread("lib/splicing_factors.txt") %>% 
  dplyr::filter(geneSymbol %in% c("CELF1", "MBNL1", "MBNL2"))

Otero_metadata <- fread("data/sample_metadata.csv") %>%
  dplyr::filter(dataset == "Otero et al. (2021)")

Otero_psi <- fread("data/psi_data_filtered_annotated.csv") %>%
  dplyr::filter(dataset == "Otero et al. (2021)") %>%
  dplyr::mutate(group = factor(group, levels = c("DM1", "Unaffected")))

Otero_log2CPM <- fread("data/log2CPM_data_filtered_annotated.csv") %>% 
  dplyr::filter(geneID %in% splicingFactors$geneID) %>%
  dplyr::filter(dataset == "Otero et al. (2021)") %>%
  dplyr::mutate(group = factor(group, levels = c("DM1", "Unaffected")))

makeCorrPlots <- function(psi, log2CPM, sampleIDs, title){
  
  # get logCPM values of splicing factor genes
  factorExpr <- log2CPM %>% 
    filter(sampleID %in% sampleIDs) %>%
    dplyr::select(geneSymbol, log2CPM, sampleID) %>% 
    spread(geneSymbol, log2CPM) %>% dplyr::select(-sampleID) %>% as.matrix()
  
  # get exon inclusion for DM1/brain related exons
  exonExpr <- psi %>% 
    filter(sampleID %in% sampleIDs) %>%
    dplyr::select(SE_event_name, exonInclusion, sampleID) %>% 
    spread(SE_event_name, exonInclusion) %>% dplyr::select(-sampleID) %>% as.matrix()
  
  # exons with increased inclusion in DM1 followed by exons with decreased inclusion
  order <- high_conf_events$SE_event_name[order(high_conf_events$order)]
  exonExpr <- exonExpr[ , match(order, colnames(exonExpr))]
  
  # correlate logCPM values of splicing factors with exon fractions of DM1 exons
  corr <- psych::corr.test(factorExpr, exonExpr, method = "spearman", adjust = "fdr")
  
  # little hack to allow labelling of sig. correlations (ggcorrplot's default can only label insig. correlations)
  sigMat <- ifelse(corr$p < 0.05, 1, 0)
  
  highlightColIndex <- which(colnames(corr$r) %in% high_conf_events$SE_event_name[high_conf_events$DM1_inclusion == "+"])
  
  fontColor <- rep("Black", ncol(corr$r))
  fontColor[highlightColIndex] <- "Black" #colPal$AccentYellow
  
  corrPlot <- corr$r %>% 
    ggcorrplot(
      method = "square", type = "full", ggtheme = theme_classic(),
      p.mat = sigMat, sig.level = 0.05, insig = "pch", pch = "*", pch.cex = 3, 
      outline.color = "white", colors = c(colPal$Red, "White", colPal$Blue), #brewer.pal(n = 3, name = "PuOr"),
      tl.srt = 60, title = title) +
    theme_classic() +
    theme(
      title = element_text(size = 14),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 60),
      axis.text = element_text(hjust = 1, size = 12, color = fontColor),
      # color = c(rep("Black", sum(high_conf_events$DM1_inclusion == "-")),
      #           rep(colPal$MediumBlue, sum(high_conf_events$DM1_inclusion == "+")))),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 15),
      legend.key.size = unit(1, "cm"))  +
    scale_fill_gradientn(name = "Corr", 
                         colors = colorRampPalette(c(colPal$Red, "White", colPal$Blue))(100),
                         limits = c(-1.0, 1.0))
  
  return(corrPlot)
}

Otero.all.plots <- makeCorrPlots(Otero_psi, Otero_log2CPM, 
                                 Otero_metadata$sampleID, 
                                 title = "All samples")

Otero.DM1.plots <- makeCorrPlots(Otero_psi, Otero_log2CPM, 
                                 Otero_metadata$sampleID[str_which(Otero_metadata$group, "DM1")],
                                 title = "Only DM1")

Otero.CTRL.plots <- makeCorrPlots(Otero_psi, Otero_log2CPM, 
                                  Otero_metadata$sampleID[str_which(Otero_metadata$group, "Unaffected")],
                                  title = "Only unaffected")

### ARRANGE FIGURE ###########

l1 <- get_legend(Otero.all.plots)

plot_grid(Otero.all.plots + theme(legend.position = "none"),
          Otero.DM1.plots + theme(legend.position = "none", axis.text.y = element_blank()),
          Otero.CTRL.plots + theme(legend.position = "none", axis.text.y = element_blank()),
          l1,
          labels = c(""), label_size = 20,
          ncol = 4, nrow = 1, rel_widths = c(1.65,1.025,1.025,0.4), 
          scale = c(0.9,0.9,0.9,1))

ggsave2(filename = "results/Figure_S8.pdf", width = 12, height = 8, dpi = 300)

