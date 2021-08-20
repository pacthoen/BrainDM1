#####################
### Figure_4_S6.R ###
#####################

suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(stringr) )
suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(tidyr) )
suppressPackageStartupMessages( library(ggplot2) )
suppressPackageStartupMessages( library(cowplot) )
suppressPackageStartupMessages( library(ggcorrplot) )
suppressPackageStartupMessages( library(RColorBrewer) )
suppressPackageStartupMessages( library(extrafonts) )

colPal <- list(Blue = "#003366", 
               Red = "#E31B23", 
               MediumBlue = "#005CAB", 
               LightBlue = "#DCEEF3", 
               AccentYellow = "#FFC325", 
               Grey = hsv(0,0,0.5),
               CoolGrey = "#E6F1EE")

theme_set(theme_classic(base_family = "Arial"))

# plotting function
makeCorrplot <- function(corrMat, pMat, title = "", highlight = NULL, wCoeff=F){
  
  if( wCoeff == TRUE ){
    pch <- ""
    
  } else {
    # little hack to allow labeling of sig. correlations (ggcorrplot's default can only label insig. correlations)
    pMat <- ifelse(pMat < 0.05, 1, 0)
    pch <- "*"
  }

  fontColor <- rep("Black", ncol(corrMat))
  fontColor[highlight] <- "Black" # colPal$AccentYellow (do not highlight these splicing events)
  
  corrPlot <- corrMat %>% 
    ggcorrplot(
      method = "square", type = "full", ggtheme = theme_classic(),
      p.mat = pMat, sig.level = 0.05, lab = wCoeff, lab_size = 2, pch = pch, pch.cex = 3, 
      outline.color = "white", colors = c(colPal$Red, "White", colPal$Blue),
      tl.srt = 60) +
    ggtitle(title) +
    theme_classic() +
    theme(
      title = element_text(size = 12),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 60),
      axis.text = element_text(hjust = 1, size = 12, color = fontColor),
                               # color = c(rep("Black", sum(high_conf_events$DM1_inclusion == "-")),
                               #           rep(colPal$MediumBlue, sum(high_conf_events$DM1_inclusion == "+")))),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.key.size = unit(1, "cm")) +
    scale_fill_gradientn(name = "Corr", 
                         colors = colorRampPalette(c(colPal$Red, "White", colPal$Blue))(100),
                         limits = c(-1.0, 1.0))
  
  return(corrPlot)
} 

### LOAD DATA ###########

splicingFactors <- fread("lib/splicing_factors.txt") %>% 
  dplyr::filter(geneSymbol %in% c("CELF1", "MBNL1", "MBNL2"))

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

GTEx_psi <- psi_data_filtered %>%
  dplyr::filter(
    dataset == "GTEx",
    tissue == "Frontal Cortex (BA9)") 

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

GTEx_log2CPM <- log2CPM_data %>%
  dplyr::filter(
    dataset == "GTEx",
    tissue == "Frontal Cortex (BA9)")

### CORRELATION BETWEEN EXONS AND SPLICING FACTORS ##############

# get logCPM values of splicing factor genes
Otero_factorExpr <- Otero_log2CPM %>% 
  dplyr::select(geneSymbol, log2CPM, sampleID) %>% 
  spread(geneSymbol, log2CPM) %>% dplyr::select(-sampleID) %>% as.matrix()

# get exon inclusion for DM1/brain related exons
Otero_exonExpr <- Otero_psi %>% 
  dplyr::select(SE_event_name, exonInclusion, sampleID) %>% 
  spread(SE_event_name, exonInclusion) %>% dplyr::select(-sampleID) %>% as.matrix()

# correlate logCPM values of splicing factors with exon fractions of DM1 exons
Otero_corr <- psych::corr.test(Otero_factorExpr, Otero_exonExpr, method = "spearman", adjust = "fdr")

# hierarchical clustering
clustered <- hclust(dist(1 - t(Otero_corr$r), diag = T, upper = T))

# change order based on clustering
eventOrder <- colnames(Otero_corr$r)[clustered$order]
Otero_corr$r <- Otero_corr$r[, eventOrder]
Otero_corr$p <- Otero_corr$p[, eventOrder]

# save order to high_conf_events
high_conf_events$order <- match(high_conf_events$SE_event_name, eventOrder)
fwrite(high_conf_events, "data/high_confidence_events.csv")

highlightColIndex <- which(colnames(Otero_corr$r) %in% high_conf_events$SE_event_name[high_conf_events$DM1_inclusion == "+"])

Otero.corrPlot <- makeCorrplot(Otero_corr$r, Otero_corr$p, "DM1/Unaffected", highlight = highlightColIndex)

#---BRAINSPAN---#

# get logCPM values of splicing factor genes
BrainSpan_factorExpr <- BrainSpan_log2CPM %>% 
  dplyr::select(geneSymbol, log2CPM, sampleID) %>% 
  spread(geneSymbol, log2CPM) %>% dplyr::select(-sampleID) %>% as.matrix()

# get exon inclusion for DM1/brain related exons
BrainSpan_exonExpr <- BrainSpan_psi %>% 
  dplyr::select(SE_event_name, exonInclusion, sampleID) %>% 
  spread(SE_event_name, exonInclusion) %>% dplyr::select(-sampleID) %>% as.matrix()

# exons with increased inclusion in DM1 followed by exons with decreased inclusion
BrainSpan_exonExpr <- BrainSpan_exonExpr[ , match(eventOrder, colnames(BrainSpan_exonExpr))]

# correlate logCPM values of splicing factors with exon fractions of DM1 exons
BrainSpan_corr <- psych::corr.test(BrainSpan_factorExpr, BrainSpan_exonExpr, method = "spearman", adjust = "fdr")

BrainSpan.corrPlot <- makeCorrplot(BrainSpan_corr$r, BrainSpan_corr$p, "Developing brain", highlight = highlightColIndex)

#---GTEx---#

# get logCPM values of splicing factor genes
GTEx_factorExpr <- GTEx_log2CPM %>% 
  dplyr::select(geneSymbol, log2CPM, sampleID) %>% 
  spread(geneSymbol, log2CPM) %>% dplyr::select(-sampleID) %>% as.matrix()

# get exon inclusion for DM1/brain related exons
GTEx_exonExpr <- GTEx_psi %>% 
  dplyr::select(SE_event_name, exonInclusion, sampleID) %>% 
  spread(SE_event_name, exonInclusion) %>% dplyr::select(-sampleID) %>% as.matrix()

# exons with increased inclusion in DM1 followed by exons with decreased inclusion
GTEx_exonExpr <- GTEx_exonExpr[ , match(eventOrder, colnames(GTEx_exonExpr))]

# correlate logCPM values of splicing factors with exon fractions of DM1 exons
GTEx_corr <- psych::corr.test(GTEx_factorExpr, GTEx_exonExpr, method = "spearman", adjust = "fdr")

GTEx.corrPlot <- makeCorrplot(GTEx_corr$r, GTEx_corr$p, "Adult brain", highlight = highlightColIndex)

### COLOR LEGEND #########

# pdf(paste0(output_dir, "/Figure_5_Legend.pdf"), width = 10, height = 4)
# par(mar = rep(0,4))
# plot(0,xlim = c(0,1), ylim = c(0,1), type = "n")
# corrplot::colorlegend(colorRampPalette(c(colPal$Red, "White", colPal$Blue))(100),
#             #colorRampPalette(brewer.pal(n = 8, name = "PuOr"))(100), 
#             seq(-1.0,1.0,0.2), xlim = c(0,0.1))
# dev.off()

### ARRANGE FIGURE ###########

legend <- get_legend(BrainSpan.corrPlot)

plot_grid(BrainSpan.corrPlot + theme(legend.position = "none"),
          Otero.corrPlot + theme(legend.position = "none", axis.text.y = element_blank()),
          GTEx.corrPlot + theme(legend.position = "none", axis.text.y = element_blank()),
          legend,
          ncol = 4, rel_widths = c(1.65,1.025,1.025,0.4), labels = "", label_size = 20, scale = c(0.9,0.9,0.9,1))

ggsave2(filename = "results/Figure_4.pdf", width = 8, height = 10, dpi = 300)

### FIGURE WITH CORRELATION COEFFICIENTS ##############

BrainSpan.corrPlot <- makeCorrplot(BrainSpan_corr$r, BrainSpan_corr$p, "Developing brain", wCoef = T, highlight = highlightColIndex)

Otero.corrPlot <- makeCorrplot(Otero_corr$r, Otero_corr$p, "DM1/Unaffected", wCoef = T, highlight = highlightColIndex)

GTEx.corrPlot <- makeCorrplot(GTEx_corr$r, GTEx_corr$p, "Adult brain", wCoef = T, highlight = highlightColIndex)

legend <- get_legend(BrainSpan.corrPlot)

plot_grid(BrainSpan.corrPlot + theme(legend.position = "none"),
          Otero.corrPlot + theme(legend.position = "none", axis.text.y = element_blank()),
          GTEx.corrPlot + theme(legend.position = "none", axis.text.y = element_blank()),
          legend,
          ncol = 4, rel_widths = c(1.65,1.025,1.025,0.4), labels = "", label_size = 20, scale = c(0.9,0.9,0.9,1))

ggsave2(filename = "results/Figure_S6.pdf", width = 8, height = 10, dpi = 300)

