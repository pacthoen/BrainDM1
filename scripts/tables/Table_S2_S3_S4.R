########################
### Table_S2_S3_S4.R ###
########################

suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(stringr) )
suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(tidyr) )
suppressPackageStartupMessages( library(rtracklayer) )
suppressPackageStartupMessages( library(purrr) )

# load external function
source("scripts/liftOverGRange.R")

### load files #############

high_conf_events <- fread("data/high_confidence_events.csv")
miso_annotation <- import.gff("lib/SE.hg38.annotated.gff3")
BrainSpan_corr_events_age <- fread("data/BrainSpan_correlation_events_age.txt")
psi_comparison_results <- fread("data/BrainSpan_Otero_PSI_comparison_results.csv")

### create annotation tables #############

exonAnn <- miso_annotation[ match(paste0(high_conf_events$event_id, ".A.se"), miso_annotation$ID) ]

# table with skipped exon event coordinates
coorTable <- data.table(
  SE.event = high_conf_events$SE_event_name,
  `DM1 Inclusion` = ifelse(high_conf_events$DM1_inclusion == "+", "Increased", "Decreased"),
  `Developmental Transition` = ifelse(
    BrainSpan_corr_events_age$estimate[match(exonAnn$eventID, BrainSpan_corr_events_age$eventID)] > 0, "Increase", "Decrease"),
  Chr = as.character(seqnames(exonAnn)),
  Start = start(exonAnn),
  End = end(exonAnn),
  Strand = as.character(strand(exonAnn)),
  Length = width(exonAnn),
  geneID = exonAnn$geneID,
  MISO.eventID = exonAnn$eventID
)

coorTable <- arrange(coorTable, desc(`DM1 Inclusion`), SE.event)

# split event column to get access to individual exon coordinates
eventCoor <- str_split(psi_comparison_results$event_id, ":|@")

# note: event column also includes coordinates for flanking exons
# e.g. chr17:44039687:44039836:+@chr17:44049225:44049311:+@chr17:44055741:44055806:+

psi_comparison_results <- psi_comparison_results %>%
  mutate(Chr = map_chr(eventCoor, 1),
         Strand = map_chr(eventCoor, 4),
         Start_up = as.numeric(map_chr(eventCoor, 2)),
         End_up = as.numeric(map_chr(eventCoor, 3)),
         Start_se = as.numeric(map_chr(eventCoor, 6)),
         End_se = as.numeric(map_chr(eventCoor, 7)),
         Start_dn = as.numeric(map_chr(eventCoor, 10)),
         End_dn = as.numeric(map_chr(eventCoor, 11))
  ) %>%
  dplyr::rename(MISO.eventID = event_id,
                HGNC.genesymbol = gene_name,
                Ensembl.geneID = gene_id)

psi_comparison_results$Length_up <- psi_comparison_results$End_up - psi_comparison_results$Start_up + 1
psi_comparison_results$Length_se <- psi_comparison_results$End_se - psi_comparison_results$Start_se + 1
psi_comparison_results$Length_dn <- psi_comparison_results$End_dn - psi_comparison_results$Start_dn + 1

# BrainSpan significant events with dPSI and exon coordinates
sigEvents_BrainSpan <- psi_comparison_results %>%
  dplyr::select(-`Otero2021 (DM1)`, -`Otero2021 (Unaffected)`, 
                -dPSI_DM1_CTRL, -p.val_DM1_CTRL, -p.adj_DM1_CTRL, -p.adj_pre_post) %>%
  na.omit() %>%
  filter(p.val_pre_post < 0.01,
         abs(dPSI_pre_post) > 0.2)

fwrite(sigEvents_BrainSpan, "results/Table_S2.csv")

# Otero significant events with dPSI and exon coordinates
sigEvents_Otero <- psi_comparison_results %>%
  dplyr::select(-`BrainSpan (Prenatal)`, -`BrainSpan (Postnatal)`, 
                -dPSI_pre_post, -p.val_pre_post, -p.adj_pre_post, -p.adj_DM1_CTRL) %>%
  na.omit() %>%
  filter(p.val_DM1_CTRL < 0.01,
         abs(dPSI_DM1_CTRL) > 0.2)

fwrite(sigEvents_Otero, "results/Table_S3.csv")

# overlapping significant events with dPSI and exon coordinates
sigEvents_overlap <- psi_comparison_results %>%
  na.omit() %>%
  filter(p.val_DM1_CTRL < 0.01 & p.val_pre_post < 0.01,
         abs(dPSI_DM1_CTRL) > 0.2 & abs(dPSI_pre_post) > 0.2) %>%
  dplyr::select(-p.adj_pre_post, -p.adj_DM1_CTRL)

sigEvents_overlap <- cbind(
  SE.event = high_conf_events$SE_event_name[match(sigEvents_overlap$MISO.eventID, high_conf_events$event_id)],
  sigEvents_overlap)

sigEvents_overlap$Length_se_divisible_by_three <- ifelse((sigEvents_overlap$Length_se %% 3) == 0, "Yes", "No")

fwrite(sigEvents_overlap, "results/Table_S4.csv")

# ### create tables with hg19 coordinates ########
# 
# # BrainSpan
# events_hg38 <- miso_annotation[match(paste0(sigEvents_BrainSpan$MISO.eventID, ".A.se"), miso_annotation$ID)]
# events_hg19 <- liftOverGRange(events_hg38, "hg38", "hg19", libPath = "lib")
# events_hg19 <- events_hg19 %>%
#   as.data.frame() %>%
#   dplyr::select(Assembly, seqnames, start, end, width, strand, source, geneID, geneName, eventID) %>%
#   dplyr::rename(MISO.eventID.hg38 = eventID,
#     assembly = Assembly,
#     chr = seqnames,
#     geneSymbol = geneName) %>%
#   dplyr::group_by(assembly, chr, start, end, width, strand, source, geneID, geneSymbol) %>%
#   dplyr::summarise(MISO.eventID.hg38 = stringr::str_c(MISO.eventID.hg38, collapse = ";")) %>%
#   dplyr::ungroup()
# 
# fwrite(events_hg19, "results/BrainSpan_sigEvents_hg19.csv")
# 
# # Otero et al. (2021)
# events_hg38 <- miso_annotation[match(paste0(sigEvents_Otero$MISO.eventID, ".A.se"), miso_annotation$ID)]
# events_hg19 <- liftOverGRange(events_hg38, "hg38", "hg19", libPath = "lib")
# events_hg19 <- events_hg19 %>%
#   as.data.frame() %>%
#   dplyr::select(Assembly, seqnames, start, end, width, strand, source, geneID, geneName, eventID) %>%
#   dplyr::rename(MISO.eventID.hg38 = eventID,
#     assembly = Assembly,
#     chr = seqnames,
#     geneSymbol = geneName) %>%
#   dplyr::group_by(assembly, chr, start, end, width, strand, source, geneID, geneSymbol) %>%
#   dplyr::summarise(MISO.eventID.hg38 = stringr::str_c(MISO.eventID.hg38, collapse = ";")) %>%
#   dplyr::ungroup()
# 
# fwrite(events_hg19, "results/Otero_sigEvents_hg19.csv")

