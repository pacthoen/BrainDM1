##################
### Table_S1.R ###
##################

suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(stringr) )
suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(tidyr) )

### LOAD DATA ############

high_conf_events <- fread("data/high_confidence_events.csv")

BrainSpan_psi_data <- fread("data/psi_data_filtered_annotated.csv") %>%
  dplyr::filter(dataset == "BrainSpan") %>%
  dplyr::mutate(group = factor(group, levels = c("Prenatal", "Postnatal")))

BrainSpan_metadata <- fread("data/sample_metadata.csv") %>%
  dplyr::filter(dataset == "BrainSpan") %>%
  dplyr::mutate(group = factor(group, levels = c("Prenatal", "Postnatal")))


### STATISTICS #####

normality_results <- data.frame(event = character(), tissue = character(), group = character(), p.value =numeric())
for(event in unique(BrainSpan_psi_data$SE_event_name)){
  event_data <- BrainSpan_psi_data[BrainSpan_psi_data$SE_event_name %in% event]
  
  for(tissue in unique(BrainSpan_psi_data$tissue)){
    tissue_data <- event_data[event_data$tissue %in% tissue]
    
    prenatalResults <- shapiro.test(tissue_data$exonInclusion[tissue_data$group == "Prenatal"])
    postnatalResults <- shapiro.test(tissue_data$exonInclusion[tissue_data$group == "Postnatal"])
    
    normality_results <- rbind(normality_results,
                              data.frame(event, tissue, 
                                         group = c("Prenatal", "Postnatal"), 
                                         p.value = c(prenatalResults$p.value, postnatalResults$p.value)
                              )
    )
  }
}
normality_results$p.adj <- p.adjust(normality_results$p.value, method = "fdr")

cat(paste0("Shapiro-Wilk test results:\n", sum(normality_results$p.adj < 0.05), " out of ", nrow(normality_results), " observations are normally distributed.\n"))

kruskal_results <- data.frame()
for(event in unique(BrainSpan_psi_data$SE_event_name)){
  event_data <- BrainSpan_psi_data[BrainSpan_psi_data$SE_event_name %in% event, ]
  
  for(group_ in levels(BrainSpan_psi_data$group)){
    group_data <- event_data %>% dplyr::filter(as.character(group) %in% group_)
    
    res <- kruskal.test(group_data$exonInclusion, group_data$tissue)
    
    psiMean <- group_data %>% 
      dplyr::group_by(tissue) %>% 
      dplyr::summarise(mean = mean(exonInclusion)) %>%
      tidyr::spread(key = tissue, value = mean)
    
    eventID <- high_conf_events$event_id[match(event, high_conf_events$SE_event_name)]
    
    kruskal_results <- rbind(kruskal_results,
                            data.frame(event, MISO.eventID = eventID, group = group_, psiMean, p.value = res$p.value)
    )
  }
}
kruskal_results$p.value.FDR <- p.adjust(p = kruskal_results$p.value, method = "fdr")

kruskal_results <- dplyr::rename(kruskal_results, 
                                SE.event = event,
                                mean.psi.medial = "Medial.prefrontal.cortex",
                                mean.psi.orbital = "Orbital.prefrontal.cortex",
                                mean.psi.dorsolateral = "Dorsolateral.prefrontal.cortex",
                                mean.psi.ventrolateral = "Ventrolateral.prefrontal.cortex")

cat(paste0("Kruskal-Wallis test results:\n", sum(kruskal_results$p.adj < 0.05), " out of ", nrow(kruskal_results), " events are significantly different between frontal cortex subregions.\n"))

fwrite(kruskal_results, "results/Table_S1.csv")
