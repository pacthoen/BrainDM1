###################
### Figure_S9.R ###
###################

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

### LOAD DATA ########

sample_metadata <- fread("data/sample_metadata.csv")

high_conf_events <- fread("data/high_confidence_events.csv")

# define order of events based on DM1 inclusion and gene name
eventOrder <- high_conf_events$SE_event_name[ order(high_conf_events$DM1_inclusion, high_conf_events$SE_event_name, decreasing = T)]

splicingFactors <- fread("lib/splicing_factors.txt") %>% 
  dplyr::filter(geneSymbol %in% c("CELF1", "MBNL1", "MBNL2"))

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

### CORRELATION BETWEEN EVENTS AND SPLICING FACTORS ##############

# get logCPM values of splicing factor genes
Otero_factorExpr <- Otero_log2CPM %>% 
  dplyr::select(geneSymbol, log2CPM, sampleID) %>% 
  spread(geneSymbol, log2CPM) %>% dplyr::select(-sampleID) %>% as.matrix()

# get exon inclusion for DM1/brain related events
Otero_exonExpr <- Otero_psi %>% 
  dplyr::select(SE_event_name, exonInclusion, sampleID) %>% 
  spread(SE_event_name, exonInclusion) %>% dplyr::select(-sampleID) %>% as.matrix()

# correlate logCPM values of splicing factors with psi of DM1 events
Otero_corr <- psych::corr.test(Otero_factorExpr, Otero_exonExpr, method = "spearman", adjust = "fdr")

# hierarchical clustering
clustered <- hclust(dist(1 - t(Otero_corr$r), diag = T, upper = T))

# change order based on clustering
eventOrder <- colnames(Otero_corr$r)[clustered$order]
Otero_corr$r <- Otero_corr$r[, eventOrder]
Otero_corr$p <- Otero_corr$p[, eventOrder]

events <- c("ADD1", "DMD")
genes <- c("CELF1","MBNL1", "MBNL2")

anno <- data.frame(SE_event_name = character(), geneSymbol = character(), dataset = character(), label = character())

for(event in events) {
  for(gene in genes) {
    anno <- rbind(anno,
                  data.frame(
                    SE_event_name = event,
                    geneSymbol = gene,
                    dataset = "Otero et al. (2021)",
                    label = paste0("R = ", sprintf("%0.3f", round(Otero_corr$r[gene, event], digits = 3)))
                  )
    )
  }
}

# get logCPM values of splicing factor genes
BrainSpan_factorExpr <- BrainSpan_log2CPM %>% 
  dplyr::select(geneSymbol, log2CPM, sampleID) %>% 
  spread(geneSymbol, log2CPM) %>% dplyr::select(-sampleID) %>% as.matrix()

# get exon inclusion for DM1/brain related events
BrainSpan_exonExpr <- BrainSpan_psi %>% 
  dplyr::select(SE_event_name, exonInclusion, sampleID) %>% 
  spread(SE_event_name, exonInclusion) %>% dplyr::select(-sampleID) %>% as.matrix()

# events with increased inclusion in DM1 followed by events with decreased inclusion
BrainSpan_exonExpr <- BrainSpan_exonExpr[ , match(eventOrder, colnames(BrainSpan_exonExpr))]

# correlate logCPM values of splicing factors with psi of DM1 events
BrainSpan_corr <- psych::corr.test(BrainSpan_factorExpr, BrainSpan_exonExpr, method = "spearman", adjust = "fdr")

for(event in events) {
  for(gene in genes) {
    anno <- rbind(anno,
                  data.frame(
                    SE_event_name = event,
                    geneSymbol = gene,
                    dataset = "BrainSpan",
                    label = paste0("R = ", sprintf("%0.3f", round(BrainSpan_corr$r[gene, event], digits = 3)))
                  )
    )
  }
}    

# get logCPM values of splicing factor genes
GTEx_factorExpr <- GTEx_log2CPM %>% 
  dplyr::select(geneSymbol, log2CPM, sampleID) %>% 
  spread(geneSymbol, log2CPM) %>% dplyr::select(-sampleID) %>% as.matrix()

# get exon inclusion for DM1/brain related events
GTEx_exonExpr <- GTEx_psi %>% 
  dplyr::select(SE_event_name, exonInclusion, sampleID) %>% 
  spread(SE_event_name, exonInclusion) %>% dplyr::select(-sampleID) %>% as.matrix()

# events with increased inclusion in DM1 followed by events with decreased inclusion
GTEx_exonExpr <- GTEx_exonExpr[ , match(eventOrder, colnames(GTEx_exonExpr))]

# correlate logCPM values of splicing factors with psi of DM1 events
GTEx_corr <- psych::corr.test(GTEx_factorExpr, GTEx_exonExpr, method = "spearman", adjust = "fdr")

for(event in events) {
  for(gene in genes) {
    anno <- rbind(anno,
                  data.frame(
                    SE_event_name = event,
                    geneSymbol = gene,
                    dataset = "GTEx",
                    label = paste0("R = ", sprintf("%0.3f", round(GTEx_corr$r[gene, event], digits = 3)))
                  )
    )
  }
} 

### SCATTERPLOT ####

anno$x <- 6
anno$y[anno$SE_event_name == "DMD"] <- 0.95
anno$y[anno$SE_event_name == "ADD1"] <- 0.02


BrainSpan_plotData <- merge(x = BrainSpan_psi %>% 
                              dplyr::select(sampleID, exonInclusion, SE_event_name, group) %>% 
                              dplyr::mutate(sampleID = as.character(sampleID)),
                            y = BrainSpan_log2CPM %>% dplyr::select(sampleID, log2CPM, geneSymbol),
                            by = "sampleID"
)
BrainSpan_plotData$dataset <- "BrainSpan"

Otero_plotData <- merge(x = Otero_psi %>% 
        dplyr::select(sampleID, exonInclusion, SE_event_name, group) %>% 
        dplyr::mutate(sampleID = as.character(sampleID)),
      y = Otero_log2CPM %>% dplyr::select(sampleID, log2CPM, geneSymbol),
      by = "sampleID"
)
Otero_plotData$dataset <- "Otero et al. (2021)"

plotData <- rbind(BrainSpan_plotData, Otero_plotData)

makeScatterplot <- function(data, event, gene, dataset_, colors) {
  data %>%
    dplyr::filter(
      dataset %in% dataset_,
      geneSymbol %in% gene,
      SE_event_name %in% event) %>%
    ggplot() +
    geom_point(aes(x = log2CPM, y = exonInclusion, colour = group), size = 0.8) +
    scale_x_continuous(breaks = c(2, 6, 10), limits = c(2,10)) +
    scale_colour_manual(values = colors) +
    coord_fixed(ratio = 1) +
    xlab("log2(CPM)") + ylab(expression(psi)) +
    ylim(0.0,1.0) + 
    guides(colour = guide_legend(override.aes = list(size=5))) +
    theme_classic() +
    theme(text = element_text(size = 20), 
          aspect.ratio = 1,
          legend.title = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5),
          #strip.background = element_rect(color="white", fill="white"),
          strip.background = element_blank(),
          strip.text = element_blank(),
          strip.placement = "outside",
          panel.spacing = unit(2, "lines")) +
    geom_text(data = dplyr::filter(anno, dataset %in% dataset_), aes(x = x, y = y, label = label), size = 4) +
    facet_wrap(SE_event_name ~ geneSymbol, ncol = 3) 
}

plot1 <- makeScatterplot(plotData, 
                         event = c("ADD1", "DMD"), 
                         gene = c("CELF1", "MBNL1", "MBNL2"), 
                         dataset_ = "BrainSpan",
                         colors = c(colPal$MediumBlue, colPal$Grey)
)
plot2 <- makeScatterplot(plotData, 
                         event = c("ADD1", "DMD"), 
                         gene = c("CELF1", "MBNL1", "MBNL2"), 
                         dataset_ = "Otero et al. (2021)",
                         colors = c(colPal$Red, colPal$Grey)
)

l_BrainSpan <- get_legend(plot1)
l_Otero <- get_legend(plot2)

plot_grid(
  plot1 + theme(legend.position = "none"),
  plot2 + theme(legend.position = "none"),
  ncol = 2, labels = "", label_size = 20, scale = c(0.9,0.9)
)

ggsave2(filename = "results/Figure_S9.pdf", width = 10, height = 6, dpi = 300)

