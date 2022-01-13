####################
### Figure_S12.R ###
####################

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
source(paste0(repo_dir, "/analysis/functions/partialCorr.R"))

### LOAD DATA ##########

high_conf_events <- fread("data/high_confidence_events.csv")

# define order of events based on DM1 inclusion and gene name
eventOrder <- high_conf_events$SE_event_name[ order(high_conf_events$DM1_inclusion, high_conf_events$SE_event_name, decreasing = T)]

splicingFactors <- fread("lib/splicing_factors.txt") %>% 
  dplyr::filter(geneSymbol %in% c("CELF1", "MBNL1", "MBNL2"))

BrainSpan_psi <- fread("data/psi_data_filtered_annotated.csv") %>% 
  dplyr::filter(dataset == "BrainSpan") %>%
  dplyr::mutate(
    group = factor(group, levels = c("Prenatal", "Postnatal")),
    age = age_in_days)

BrainSpan_log2CPM <- fread("data/log2CPM_data_filtered_annotated.csv") %>% 
  dplyr::filter(geneID %in% splicingFactors$geneID) %>%
  dplyr::filter(dataset == "BrainSpan") %>%
  dplyr::mutate(
    group = factor(group, levels = c("Prenatal", "Postnatal")),
    age = age_in_days)

BrainSpan_metadata <- fread("data/sample_metadata.csv") %>%
  dplyr::filter(dataset == "BrainSpan")

makeCorrPlots <- function(psi, log2CPM, sampleIDs, title){
  
  # get logCPM values of splicing factor genes
  factorExpr <- log2CPM %>% 
    filter(sampleID %in% sampleIDs) %>%
    dplyr::select(geneSymbol, log2CPM, sampleID) %>% 
    spread(geneSymbol, log2CPM) %>% dplyr::select(-sampleID) %>% as.matrix()
  
  # get exon inclusion for DM1/brain related events
  exonExpr <- psi %>% 
    filter(sampleID %in% sampleIDs) %>%
    dplyr::select(SE_event_name, exonInclusion, sampleID) %>% 
    spread(SE_event_name, exonInclusion) %>% dplyr::select(-sampleID) %>% as.matrix()
  
  # events with increased inclusion in DM1 followed by events with decreased inclusion
  order <- high_conf_events$SE_event_name[order(high_conf_events$order)]
  exonExpr <- exonExpr[ , match(order, colnames(exonExpr))]
  
  # correlate logCPM values of splicing factors with exon fractions of DM1 events
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

BrainSpan.all.plots <- makeCorrPlots(BrainSpan_psi, BrainSpan_log2CPM, 
                                     BrainSpan_metadata$sampleID, 
                                     title = "All samples")

BrainSpan.prenatal.plots <- makeCorrPlots(BrainSpan_psi, BrainSpan_log2CPM, 
                                          BrainSpan_metadata$sampleID[str_which(BrainSpan_metadata$group, "Prenatal")],
                                          title = "Only prenatal")

BrainSpan.postnatal.plots <- makeCorrPlots(BrainSpan_psi, BrainSpan_log2CPM, 
                                           BrainSpan_metadata$sampleID[str_which(BrainSpan_metadata$group, "Postnatal")],
                                           title = "Only postnatal")

### ARRANGE FIGURE ###########

l1 <- get_legend(BrainSpan.all.plots)

plot_grid(BrainSpan.all.plots + theme(legend.position = "none"),
          BrainSpan.prenatal.plots + theme(legend.position = "none", axis.text.y = element_blank()),
          BrainSpan.postnatal.plots + theme(legend.position = "none", axis.text.y = element_blank()),
          l1,
          labels = c(""), label_size = 20,
          ncol = 4, nrow = 1, rel_widths = c(1.65,1.025,1.025,0.4), 
          scale = c(0.9,0.9,0.9,1))

ggsave2(filename = "results/Figure_S12.pdf", width = 12, height = 8, dpi = 300)

