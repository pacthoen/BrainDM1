######################
### Figure_S1_S4.R ###
######################

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
suppressPackageStartupMessages( library(ggpubr) )

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

### LOAD DATA ############

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

GTEx_psi_all <- psi_data_filtered %>% 
  dplyr::filter(dataset == "GTEx")

GTEx_psi_FC <- psi_data_filtered %>%
  dplyr::filter(
    dataset == "GTEx",
    tissue == "Frontal Cortex (BA9)") 

### STATISTICS #####

### Otero's unaffected samples vs GTEx frontal cortex

data <- data.frame(sampleID = character(), SE_event_name = character(), exonInclusion = numeric(), dataset = character())

Otero_psi_unaff <- Otero_psi[Otero_psi$group == "Unaffected",]

data <- rbind(data,
              data.frame(sampleID = Otero_psi_unaff$sampleID, SE_event_name = Otero_psi_unaff$SE_event_name, exonInclusion= Otero_psi_unaff$exonInclusion,  dataset = "Otero et al. (2021)"),
              data.frame(sampleID = GTEx_psi_FC$sampleID, SE_event_name = GTEx_psi_FC$SE_event_name, exonInclusion = GTEx_psi_FC$exonInclusion, dataset = "GTEx")
)

# check for normality
normalityResults <- data.frame(event = character(), tissue = character(), group = character(), p.value =numeric())
for(event in unique(data$SE_event_name)){
  eventData <- data[data$SE_event_name == event,]
  
  Otero_res <- shapiro.test(eventData$exonInclusion[eventData$dataset == "Otero et al. (2021)"])
  GTEx_res <- shapiro.test(eventData$exonInclusion[eventData$dataset == "GTEx"])
  
  normalityResults <- rbind(normalityResults,
                            data.frame(event, 
                                       group = c("Otero et al. (2021)", "GTEx"), 
                                       p.value = c(Otero_res$p.value, GTEx_res$p.value)
                            )
  )
}
normalityResults$p.adj <- p.adjust(normalityResults$p.value, method = "fdr")

cat(paste0("Shapiro-Wilk test results:\n", sum(normalityResults$p.adj < 0.05), " out of ", nrow(normalityResults), " observations are normally distributed.\n"))

# non-parametric kruskal wallis test
kruskalResults <- data.frame(group = character(), p.value =numeric())
for(event in unique(data$SE_event_name)){
  
  eventData <- data[data$SE_event_name == event,]
  
  res <- kruskal.test(eventData$exonInclusion, eventData$dataset)
  
  dPSI <- abs(mean(eventData$exonInclusion[eventData$dataset == "Otero et al. (2021)"]) -
                mean(eventData$exonInclusion[eventData$dataset == "GTEx"])
  )
  
  kruskalResults <- rbind(kruskalResults,
                          data.frame(event, p.value = res$p.value, dPSI = dPSI
                          )
  )
}
kruskalResults$p.adj <- p.adjust(kruskalResults$p.value, method = "fdr")

cat(paste0("Kruskal-Wallis test results:\nPSI of ", sum(kruskalResults$p.adj < 0.05), " out of ", nrow(kruskalResults), " events is significantly different between unaffected controls (Otero et al., 2021) and GTEx subjects.\n"))
cat(paste0("Significant events:\n", paste(kruskalResults$event[kruskalResults$p.adj < 0.05], collapse = ", "), "\n"))
cat(paste0("Non-significant events:\n", paste(kruskalResults$event[kruskalResults$p.adj > 0.05], collapse = ", "), "\n"))
cat(paste0("Events with |dPSI| > 0.2:\n", paste(kruskalResults$event[kruskalResults$dPSI > 0.2], collapse = ", "), "\n"))

nsEvents <- kruskalResults$event[kruskalResults$p.adj > 0.05]

### COMPARE PSI FOR NON-SIGNIFICANT EVENTS ###

BrainSpan_psi$dataset <- BrainSpan_psi$group

data <- rbind(
  dplyr::select(Otero_psi_unaff, sampleID, SE_event_name, exonInclusion, tissue, dataset),
  dplyr::select(BrainSpan_psi, sampleID, SE_event_name, exonInclusion, tissue, dataset),
  dplyr::select(GTEx_psi_all, sampleID, SE_event_name, exonInclusion, tissue, dataset)
)

# remove everything within parentheses from tissue description and abbreviate tissue type
data$tissue <- abbreviate(str_replace(data$tissue, " \\(.*\\)", ""))

data$tissue <- str_replace(data$tissue, "Mdpc", "Medial")
data$tissue <- str_replace(data$tissue, "Orpc", "Orbital")
data$tissue <- str_replace(data$tissue, "Vnpc", "Ventrolateral")
data$tissue <- str_replace(data$tissue, "Drpc", "Dorsolateral")

data$dataset <- factor(data$dataset, levels = c("Otero et al. (2021)", "Prenatal", "Postnatal", "GTEx"))

incEvents <- c("MBNL2a", "MBNL2c", "PLA2G6", "TACC2")
decEvents <- c("ARHGAP44", "ATP1B3", "MAPT", "SORBS1b")

# function to plot boxplots
makeBoxplot <- function(data, events, fill = "white", color) {
  data %>% 
    dplyr::filter(SE_event_name %in% events) %>%
    dplyr::mutate(SE_event_name = factor(SE_event_name, levels = events)) %>%
    ggboxplot(y = "exonInclusion", x = "tissue", fill = fill) +
    xlab("") + ylab(expression(psi)) +
    scale_y_continuous(limits=c(0.0,1.0), breaks = seq(0,1,by=0.25)) +
    scale_fill_manual(values = color) +
    theme_classic() +
    theme(
      text = element_text(size = 20),
      legend.position = "top",
      legend.title = element_blank(),
      legend.key.size = unit(1, "cm"),
      legend.text = element_text(size = 18),
      axis.title.y = element_text(angle = 0, vjust = 0.5),
      axis.text.x = element_text(angle = -60, hjust = 0),
      strip.background = element_rect(color="white", fill="white"),
      strip.placement = "outside",
      panel.spacing = unit(2, "lines")) +
    facet_wrap(~ SE_event_name, ncol = 4, strip.position = "top", scales = "free_x") 
  
  # stat_pvalue_manual(
  #   filter(stat_results,
  #          dataset == dataset_,
  #          SE_event_name %in% data$SE_event_name),
  #   label = "p.adj.sig", y.position = 1.1)
}

boxplots_BrainSpan <- makeBoxplot(
  data = filter(data, dataset %in% c("Prenatal", "Postnatal")), 
  events = c(incEvents, decEvents),
  fill = "dataset",
  color = c(colPal$MediumBlue, colPal$Grey)
)
boxplots_BrainSpan

ggsave2(filename = "results/Figure_S1.pdf", width = 10, height = 10, dpi = 300)

boxplots_GTEx <- makeBoxplot(
  data = filter(data, dataset %in% c("Otero et al. (2021)", "GTEx")), 
  events = c(incEvents, decEvents), 
  fill = "dataset",
  color = c(colPal$Red, colPal$Grey)
)
boxplots_GTEx

ggsave2(filename = "results/Figure_S4.pdf", width = 10, height = 8, dpi = 300)


