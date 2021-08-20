####################
### comparePSI.R ###
#####################

# This script compares the transcriptome-wide PSI of exon-skipping events between
# prenatal and postnatal samples from the BrainSpan dataset and
# DM1 and unaffected samples from the Otero et al. (2021) dataset
# by performing a Wilcoxon rank-sum test.
#
# It will create important data files for further downstream analysis:
# - high_confidenve_events.txt
# - BrainSpan_Otero_PSI_comparison_results.csv
# - psi_data_filtered_annotated.csv

# load libraries
suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(tidyr) )
suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(matrixTests) )

# load files 
sample_metadata <- data.table::fread("data/sample_metadata.csv")
psi_data <- data.table::fread("data/psi_data.csv")
miso_annotation <- rtracklayer::import.gff("lib/SE.hg38.annotated.gff3")

# get vector of sampleIDs for each dataset
BrainSpan_samples <- sample_metadata$sampleID[sample_metadata$dataset == "BrainSpan"]
Otero_samples <- sample_metadata$sampleID[sample_metadata$dataset == "Otero et al. (2021)"]
GTEx_samples <- sample_metadata$sampleID[sample_metadata$dataset == "GTEx"]

# only consider events without any missing values (consider each dataset separately)
keep <- complete.cases(psi_data[, BrainSpan_samples, with = F]) | 
  complete.cases(psi_data[, Otero_samples, with = F]) |
  complete.cases(psi_data[, GTEx_samples, with = F])

psi_data <- psi_data[keep, ]

# gather sample ids into one column (instead of one column for every sample id)
psi_data <- data.table::melt(psi_data, 
  variable.name = "sampleID", 
  value.name = "exonInclusion", 
  id.vars = c("eventID"))

psi_data <- psi_data %>% 
  merge(y = sample_metadata[, c("sampleID", "dataset", "group")],  
    by = "sampleID") %>%
  dplyr::mutate(group = dplyr::recode(group, 
    "Prenatal" = "BrainSpan (Prenatal)",
    "Postnatal" = "BrainSpan (Postnatal)",
    "DM1" = "Otero2021 (DM1)",
    "Unaffected" = "Otero2021 (Unaffected)")) 

# calculate mean psi per group
psi_mean <- psi_data %>% 
  dplyr::group_by(group, eventID) %>% 
  dplyr::summarise(mean = mean(exonInclusion)) %>% # here one should probably set na.rm = TRUE to avoid lots of NA values
  tidyr::spread(group, mean) 

# calculate delta psi
psi_mean$dPSI_pre_post <- psi_mean$`BrainSpan (Prenatal)` - psi_mean$`BrainSpan (Postnatal)`
psi_mean$dPSI_DM1_CTRL <- psi_mean$`Otero2021 (DM1)` - psi_mean$`Otero2021 (Unaffected)`

# prepare data  for wilcoxon test 
grouping <- unique(data.table::data.table(psi_data[, c("sampleID", "group")]))
psi_data <- psi_data %>% 
  dplyr::select(sampleID, eventID, exonInclusion) %>%
  tidyr::spread(sampleID, exonInclusion) %>%
  as.data.table()

### Wilcoxon Test ##################

## BrainSpan
prenatal <- psi_data[complete.cases(psi_data[, BrainSpan_samples, with = F]), 
  c("eventID", as.character(grouping$sampleID[grouping$group == "BrainSpan (Prenatal)"])), with = F]

postnatal <- psi_data[complete.cases(psi_data[, BrainSpan_samples, with = F]), 
  c("eventID", as.character(grouping$sampleID[grouping$group == "BrainSpan (Postnatal)"])), with = F]

wilcox_BrainSpan <- matrixTests::row_wilcoxon_twosample(x = prenatal[,-"eventID"], 
  y = postnatal[,-"eventID"], 
  alternative = "two.sided", 
  exact = F) 

wilcox_BrainSpan$pvalue[is.na(wilcox_BrainSpan$pvalue)] <- 1.0
wilcox_BrainSpan$p.adj <- p.adjust(wilcox_BrainSpan$pvalue, method = "fdr")
stopifnot(all(prenatal$eventID == postnatal$eventID))
wilcox_BrainSpan$eventID <- prenatal$eventID
wilcox_BrainSpan <- dplyr::rename(wilcox_BrainSpan, p.val_pre_post = pvalue, p.adj_pre_post = p.adj)

## Otero et al. (2021)
DM1 <- psi_data[complete.cases(psi_data[, Otero_samples, with = F]) , 
  c("eventID", as.character(grouping$sampleID[grouping$group == "Otero2021 (DM1)"])), with = F]

unaffected <- psi_data[complete.cases(psi_data[, Otero_samples, with = F]) , 
  c("eventID", as.character(grouping$sampleID[grouping$group == "Otero2021 (Unaffected)"])), with = F]

wilcox_Otero <- matrixTests::row_wilcoxon_twosample(x = DM1[,-"eventID"], 
  y = unaffected[,-"eventID"], 
  alternative = "two.sided", 
  exact = F)

wilcox_Otero$pvalue[is.na(wilcox_Otero$pvalue)] <- 1.0
wilcox_Otero$p.adj <- p.adjust(wilcox_Otero$pvalue, method = "fdr")
stopifnot(all(DM1$eventID == unaffected$eventID))
wilcox_Otero$eventID <- DM1$eventID
wilcox_Otero <- dplyr::rename(wilcox_Otero, p.val_DM1_CTRL = pvalue, p.adj_DM1_CTRL = p.adj)

# create final results table
results <- merge(
  x = wilcox_BrainSpan[c("eventID", "p.val_pre_post", "p.adj_pre_post")],
  y = wilcox_Otero[c("eventID", "p.val_DM1_CTRL", "p.adj_DM1_CTRL")],
  by = "eventID",
  all = TRUE) %>%
  merge(psi_mean, by = "eventID") %>%
  dplyr::mutate(DM1_inclusion = ifelse(dPSI_DM1_CTRL > 0, "+", "-"),
    gene_id = miso_annotation$geneID[match(eventID, miso_annotation$eventID)],
    gene_name = miso_annotation$geneName[match(eventID, miso_annotation$eventID)],
    SE_event_name = miso_annotation$eventName[match(eventID, miso_annotation$eventID)]) %>%
  dplyr::rename(event_id = eventID) %>%
  dplyr::filter(!is.na(SE_event_name))

# drop SE_event_name from results 
results$SE_event_name <- NULL

# create table with filtered results
results_filtered <- results %>% dplyr::filter(
  p.val_pre_post < 0.01, 
  p.val_DM1_CTRL < 0.01,
  abs(dPSI_pre_post) > 0.2,
  abs(dPSI_DM1_CTRL) > 0.2)

# use gene name as SE_event_name 
results_filtered$SE_event_name <- results_filtered$gene_name

# add letters to duplicate genes
duplicates <- unique(results_filtered$SE_event_name[duplicated(results_filtered$SE_event_name)])
for(duplicate in duplicates){
  i <- which(results_filtered$SE_event_name == duplicate)
  results_filtered$SE_event_name[i] <- paste0(results_filtered$SE_event_name[i], letters[seq(i)])
}

# define order of exons based on DM1 inclusion and exon name
order <- results_filtered$SE_event_name[ order(results_filtered$DM1_inclusion, results_filtered$SE_event_name, decreasing = F)]
results_filtered <- results_filtered[match(order, results_filtered$SE_event_name),]

# save complete and filtered results
data.table::fwrite(results_filtered, "data/high_confidence_events.csv")
data.table::fwrite(results, "data/BrainSpan_Otero_PSI_comparison_results.csv")

# filter complete psi data for high confidence events
psi_data_filtered <- psi_data %>%
  dplyr::filter(eventID %in% results_filtered$event_id) %>%
  data.table::melt(variable.name = "sampleID", value.name = "exonInclusion", id.vars = "eventID") %>%
  merge(y = sample_metadata, by = "sampleID") %>%
  dplyr::mutate(
    SE_event_name = results_filtered$SE_event_name[ match(eventID, results_filtered$event_id)],
    SE_event_name = factor(SE_event_name, levels = order))

# save filtered and annotated psi data
data.table::fwrite(psi_data_filtered, "data/psi_data_filtered_annotated.csv")

### Estimation of false positives #########

Otero_shuffled <- t(apply(psi_data[psi_data$eventID %in% results$event_id[ abs(results$dPSI_DM1_CTRL) > 0.2 & results$p.val_DM1_CTRL < 0.01], 
  as.character(grouping$sampleID[grouping$group %in% c("Otero2021 (DM1)",
    "Otero2021 (Unaffected)")]), 
  with = F],
  1, sample))

colnames(Otero_shuffled) <- as.character(grouping$sampleID[grouping$group %in% c("Otero2021 (DM1)", "Otero2021 (Unaffected)")])

wilcox_Otero_shuffled <- row_wilcoxon_twosample(Otero_shuffled[,as.character(grouping$sampleID[grouping$group == "Otero2021 (DM1)"])], 
  Otero_shuffled[,as.character(grouping$sampleID[grouping$group == "Otero2021 (Unaffected)"])], 
  "two.sided", exact = F)

print(paste0("Otero2021 - Sig. exons after shuffling: ", sum(wilcox_Otero_shuffled$pvalue < 0.01)))
print(paste0("Percentage of false positives: ", round(sum(wilcox_Otero_shuffled$pvalue < 0.01) / nrow(wilcox_Otero_shuffled), 4)))

BrainSpan_shuffled <- t(apply(psi_data[psi_data$eventID %in% results$event_id[ abs(results$dPSI_pre_post) > 0.2 & results$p.val_pre_post < 0.01], 
  as.character(grouping$sampleID[grouping$group %in% c("BrainSpan (Prenatal)",
    "BrainSpan (Postnatal)")]), 
  with = F],
  1, sample))

colnames(BrainSpan_shuffled) <- as.character(grouping$sampleID[grouping$group %in% c("BrainSpan (Prenatal)", "BrainSpan (Postnatal)")])

wilcox_BrainSpan_shuffled <- row_wilcoxon_twosample(BrainSpan_shuffled[,as.character(grouping$sampleID[grouping$group == "BrainSpan (Prenatal)"])], 
  BrainSpan_shuffled[,as.character(grouping$sampleID[grouping$group == "BrainSpan (Postnatal)"])], 
  "two.sided", exact = F)

print(paste0("BrainSpan - Sig. exons after shuffling: ", sum(wilcox_BrainSpan_shuffled$pvalue < 0.01)))
print(paste0("Percentage of false positives: ", round(sum(wilcox_BrainSpan_shuffled$pvalue < 0.01) / nrow(wilcox_BrainSpan_shuffled), 4)))

