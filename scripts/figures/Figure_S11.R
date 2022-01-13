####################
### Figure_S11.R ###
####################

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

comparison_results <- fread("data/BrainSpan_Otero_PSI_comparison_results.csv")

psi_data <- fread("data/psi_data.csv")

splicingFactors <- fread("lib/splicing_factors.txt") %>% 
  dplyr::filter(geneSymbol %in% c("CELF1", "MBNL1", "MBNL2"))

high_conf_events <- fread("data/high_confidence_events.csv")

sample_metadata <- fread("data/sample_metadata.csv")

log2CPM_data <- fread("data/log2CPM_data_filtered_annotated.csv") %>% 
  dplyr::filter(geneID %in% splicingFactors$geneID)

### DATA WRANGLING ##########

# create table with filtered results
results_filtered <- comparison_results %>% 
  dplyr::filter(
    p.val_DM1_CTRL < 0.01,
    abs(dPSI_DM1_CTRL) > 0.2,
    !(event_id %in% high_conf_events$event_id)
  )

# use gene name as SE_event_name 
results_filtered$SE_event_name <- results_filtered$gene_name

# add letters to duplicate genes
duplicates <- unique(results_filtered$SE_event_name[duplicated(results_filtered$SE_event_name)])
for(duplicate in duplicates){
  i <- which(results_filtered$SE_event_name == duplicate)
  results_filtered$SE_event_name[i] <- paste0(results_filtered$SE_event_name[i], letters[seq(i)])
}

# filter complete psi data for high confidence events
psi_data_filtered <- psi_data %>%
  dplyr::filter(eventID %in% results_filtered$event_id) %>%
  data.table::melt(variable.name = "sampleID", value.name = "exonInclusion", id.vars = "eventID") %>%
  merge(y = sample_metadata, by = "sampleID") %>%
  dplyr::mutate(SE_event_name = results_filtered$SE_event_name[ match(eventID, results_filtered$event_id)])

Otero_psi <- psi_data_filtered %>%
  dplyr::filter(dataset == "Otero et al. (2021)") %>%
  dplyr::mutate(group = factor(group, levels = c("DM1", "Unaffected")))

Otero_log2CPM <- log2CPM_data %>%
  dplyr::filter(dataset == "Otero et al. (2021)") %>%
  dplyr::mutate(group = factor(group, levels = c("DM1", "Unaffected")))

### CORRELATION BETWEEN EXONS AND SPLICING FACTORS ##############

# get logCPM values of splicing factor genes
Otero_factorExpr <- Otero_log2CPM %>% 
  dplyr::select(geneSymbol, log2CPM, sampleID) %>% 
  spread(geneSymbol, log2CPM) %>% dplyr::select(-sampleID) %>% as.matrix()

# get exon inclusion for DM1/brain related exons
Otero_exonExpr <- Otero_psi %>% 
  dplyr::select(SE_event_name, exonInclusion, sampleID) %>% 
  spread(SE_event_name, exonInclusion) %>% dplyr::select(-sampleID) %>% as.matrix()

# get events with increased and decreased exon inclusion
events_increased <- results_filtered$SE_event_name[results_filtered$DM1_inclusion == "+"]
events_decreased <- results_filtered$SE_event_name[results_filtered$DM1_inclusion == "-"]

# correlate logCPM values of splicing factors with exon fractions of DM1 exons
Otero_corr_increased <- psych::corr.test(Otero_factorExpr, Otero_exonExpr[, events_increased], method = "spearman", adjust = "fdr")
Otero_corr_decreased <- psych::corr.test(Otero_factorExpr, Otero_exonExpr[, events_decreased], method = "spearman", adjust = "fdr")

# hierarchical clustering
increased_clustered <- hclust(dist(1 - t(Otero_corr_increased$r), diag = T, upper = T))
decreased_clustered <- hclust(dist(1 - t(Otero_corr_decreased$r), diag = T, upper = T))

# change order based on clustering
Otero_corr_increased$r <- Otero_corr_increased$r[, colnames(Otero_corr_increased$r)[increased_clustered$order]]
Otero_corr_decreased$r <- Otero_corr_decreased$r[, colnames(Otero_corr_decreased$r)[decreased_clustered$order]]

### PLOTTING ##############

Otero.corrPlot.increased <- makeCorrplot(Otero_corr_increased$r, Otero_corr_increased$p, "Increased in DM1", wCoef = F)
Otero.corrPlot.decreased <- makeCorrplot(Otero_corr_decreased$r, Otero_corr_decreased$p, "Decreased in DM1", wCoef = F)

Otero.corrPlot.increased
ggsave2(filename = "results/Figure_S11a.pdf", width = 8, height = 10, dpi = 300)

Otero.corrPlot.decreased
ggsave2(filename = "results/Figure_S11b.pdf", width = 8, height = 12, dpi = 300)

