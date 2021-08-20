##################
### Figure_5.R ###
##################

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

### PARTIAL CORRELATION #######

#---BRAINSPAN---#

# get exon inclusion for DM1/brain related exons
BrainSpan_exonExpr <- BrainSpan_psi %>% 
  dplyr::select(SE_event_name, exonInclusion, sampleID) %>% 
  spread(SE_event_name, exonInclusion) %>%
  as.data.table()

# exons with increased inclusion in DM1 followed by exons with decreased inclusion
setcolorder(BrainSpan_exonExpr, c("sampleID", eventOrder))

# compute zero-order correlation
corr <- psych::corr.test(as.matrix(BrainSpan_exonExpr[,-"sampleID"]), 
                         method = "spearman", adjust = "fdr")

# splicing factor covariate: logCPM values of splicing factor genes
BrainSpan_factorExpr <- BrainSpan_log2CPM %>% 
  dplyr::select(geneSymbol, log2CPM, sampleID) %>% 
  spread(geneSymbol, log2CPM) %>% 
  arrange(match(sampleID, BrainSpan_exonExpr$sampleID)) %>%
  dplyr::select(-sampleID) 

colnames(BrainSpan_factorExpr) <- paste0(colnames(BrainSpan_factorExpr), "_covariate")

# age covariate
BrainSpan_age <- sample_metadata$age_in_days[match(BrainSpan_exonExpr$sampleID, sample_metadata$sampleID)]

# compute partial correlations
pcorr_splicingFactors <- partialCorr(as.data.table(cbind(BrainSpan_exonExpr[,-"sampleID"], BrainSpan_factorExpr)), 
                                     exons = eventOrder, 
                                     covariates = colnames(BrainSpan_factorExpr))
pcorr_age <- partialCorr(as.data.table(cbind(BrainSpan_exonExpr[,-"sampleID"], age = BrainSpan_age)), 
                                     exons = eventOrder, 
                                     covariates = "age")

df <- rbind(data.frame(corr = as.vector(abs(corr$r[upper.tri(corr$r, diag=F)])),
                       covariates = "Zero-order", dataset = "Developing brain"),
            data.frame(corr = as.vector(abs(pcorr_splicingFactors$r[upper.tri(pcorr_splicingFactors$r, diag=F)])),
                       covariates = "Controlled for splicing factors", dataset = "Developing brain"),
            data.frame(corr = as.vector(abs(pcorr_age$r[upper.tri(pcorr_age$r, diag=F)])),
                       covariates = "Controlled for age", dataset = "Developing brain")
)

# create nested list with comparison data (for significance testing)
comparisons <- list(
  list(
    comparison = "Zero-order vs. age",
    dataset = "Developing brain",
    corr1 = corr$r,
    corr2 = pcorr_age$r,
    n = length(unique(BrainSpan_psi$sampleID))
  ), list(
    comparison = "Zero-order vs. splicing factors",
    dataset = "Developing brain",
    corr1 = corr$r,
    corr2 = pcorr_splicingFactors$r,
    n = length(unique(BrainSpan_psi$sampleID))
  )
)

#---OTERO---#

# get exon inclusion for DM1/brain related exons
Otero_exonExpr <- Otero_psi %>% 
  dplyr::select(SE_event_name, exonInclusion, sampleID) %>% 
  spread(SE_event_name, exonInclusion) %>%
  as.data.table()

# exons with increased inclusion in DM1 followed by exons with decreased inclusion
setcolorder(Otero_exonExpr, c("sampleID", eventOrder))

# compute zero-order correlation
corr <- psych::corr.test(as.matrix(Otero_exonExpr[,-"sampleID"]), 
                         method = "spearman", adjust = "fdr")

# splicing factor covariate: logCPM values of splicing factor genes
Otero_factorExpr <- Otero_log2CPM %>% 
  dplyr::select(geneSymbol, log2CPM, sampleID) %>% 
  spread(geneSymbol, log2CPM) %>% 
  arrange(match(sampleID, Otero_exonExpr$sampleID)) %>%
  dplyr::select(-sampleID) 

colnames(Otero_factorExpr) <- paste0(colnames(Otero_factorExpr), "_covariate")

# age covariate
Otero_age <- sample_metadata$age[match(Otero_exonExpr$sampleID, sample_metadata$sampleID)]

# group covariate
Otero_group <- data.frame(group = ifelse(
  Otero_psi$group[match(Otero_exonExpr$sampleID, Otero_psi$sampleID)] == "DM1", 1, 0))

# compute partial correlations
pcorr_splicingFactors <- partialCorr(as.data.table(cbind(Otero_exonExpr[,-"sampleID"], Otero_factorExpr)), 
                                     exons = eventOrder, 
                                     covariates = colnames(Otero_factorExpr))
pcorr_age <- partialCorr(as.data.table(cbind(Otero_exonExpr[,-"sampleID"], age = Otero_age)), 
                         exons = eventOrder, 
                         covariates = "age")
pcorr_group <- partialCorr(as.data.table(cbind(Otero_exonExpr[,-"sampleID"], Otero_group)), 
                         exons = eventOrder, 
                         covariates = colnames(Otero_group))

df <- rbind(df,
            data.frame(corr = as.vector(abs(corr$r[upper.tri(corr$r, diag=F)])),
                       covariates = "Zero-order", dataset = "DM1/Unaffected"),
            data.frame(corr = as.vector(abs(pcorr_splicingFactors$r[upper.tri(pcorr_splicingFactors$r, diag=F)])),
                       covariates = "Controlled for splicing factors", dataset = "DM1/Unaffected"),
            data.frame(corr = as.vector(abs(pcorr_age$r[upper.tri(pcorr_age$r, diag=F)])),
                       covariates = "Controlled for age", dataset = "DM1/Unaffected"),
            data.frame(corr = as.vector(abs(pcorr_group$r[upper.tri(pcorr_group$r, diag=F)])),
                       covariates = "Controlled for DM1/Unaffected", dataset = "DM1/Unaffected")
)

# append nested list with comparison data (for significance testing)
comparisons <- append(comparisons,
                      list(list(
                        comparison = "Zero-order vs. age",
                        dataset = "DM1/Unaffected",
                        corr1 = corr$r,
                        corr2 = pcorr_age$r,
                        n = length(unique(Otero_psi$sampleID))
                      ), list(
                        comparison = "Zero-order vs. splicing factors",
                        dataset = "DM1/Unaffected",
                        corr1 = corr$r,
                        corr2 = pcorr_splicingFactors$r,
                        n = length(unique(Otero_psi$sampleID))
                      ), list(
                        comparison = "Zero-order vs. DM1/Unaffected",
                        dataset = "DM1/Unaffected",
                        corr1 = corr$r,
                        corr2 = pcorr_group$r,
                        n = length(unique(Otero_psi$sampleID))
                      ))
)

#---GTEx---#

# get exon inclusion for DM1/brain related exons
GTEx_exonExpr <- GTEx_psi %>% 
  dplyr::select(SE_event_name, exonInclusion, sampleID) %>% 
  spread(SE_event_name, exonInclusion) %>%
  as.data.table()

# exons with increased inclusion in DM1 followed by exons with decreased inclusion
setcolorder(GTEx_exonExpr, c("sampleID", eventOrder))

# compute zero-order correlation
corr <- psych::corr.test(as.matrix(GTEx_exonExpr[,-"sampleID"]), 
                         method = "spearman", adjust = "fdr")

# splicing factor covariate: logCPM values of splicing factor genes
GTEx_factorExpr <- GTEx_log2CPM %>% 
  dplyr::select(geneSymbol, log2CPM, sampleID) %>% 
  spread(geneSymbol, log2CPM) %>% 
  arrange(match(sampleID, GTEx_exonExpr$sampleID)) %>%
  dplyr::select(-sampleID) 

colnames(GTEx_factorExpr) <- paste0(colnames(GTEx_factorExpr), "_covariate")

# age covariate
GTEx_age <- sample_metadata$age[match(GTEx_exonExpr$sampleID, sample_metadata$sampleID)]

# compute partial correlations
pcorr_splicingFactors <- partialCorr(as.data.table(cbind(GTEx_exonExpr[,-"sampleID"], GTEx_factorExpr)), 
                                     exons = eventOrder, 
                                     covariates = colnames(GTEx_factorExpr))
pcorr_age <- partialCorr(as.data.table(cbind(GTEx_exonExpr[,-"sampleID"], age = GTEx_age)), 
                         exons = eventOrder, 
                         covariates = "age")

df <- rbind(df,
            data.frame(corr = as.vector(abs(corr$r[upper.tri(corr$r, diag=F)])),
                       covariates = "Zero-order", dataset = "Adult brain"),
            data.frame(corr = as.vector(abs(pcorr_splicingFactors$r[upper.tri(pcorr_splicingFactors$r, diag=F)])),
                       covariates = "Controlled for splicing factors", dataset = "Adult brain"),
            data.frame(corr = as.vector(abs(pcorr_age$r[upper.tri(pcorr_age$r, diag=F)])),
                       covariates = "Controlled for age", dataset = "Adult brain")
)

df$dataset <- factor(df$dataset, levels = c("Developing brain", "DM1/Unaffected", "Adult brain"))

# append nested list with comparison data (for significance testing)
comparisons <- append(comparisons, 
                      list(list(
                        comparison = "Zero-order vs. age",
                        dataset = "Adult brain",
                        corr1 = corr$r,
                        corr2 = pcorr_age$r,
                        n = length(unique(GTEx_psi$sampleID))
                      ), list(
                        comparison = "Zero-order vs. splicing factors",
                        dataset = "Adult brain",
                        corr1 = corr$r,
                        corr2 = pcorr_splicingFactors$r,
                        n = length(unique(GTEx_psi$sampleID))
                      ))
)

df1 <- filter(df, covariates != "Controlled for splicing factors")
df2 <- filter(df, covariates %in% c("Zero-order", "Controlled for splicing factors"))

vline1 <- summarise(group_by(df1,dataset, covariates), median = median(corr))
vline2 <- summarise(group_by(df2,dataset, covariates), median = median(corr))


### SIGNIFICANCE TESTING ############

# using a paired t-test with FDR-correction

testResults <- data.frame(comparison = character(), dataset = character(), statistic = numeric(), p.value = numeric())

for(comparison in comparisons){
  # res <- psych::cortest.normal(R1 = comparison$corr1, R2 = comparison$corr2, 
  #                                n1 = comparison$n, n2 = comparison$n)
  
  # comparison <- comparisons[[1]]
  
  res <- wilcox.test(comparison$corr1, comparison$corr2, paired = TRUE)
  
  testResults <- rbind(testResults, data.frame(
    comparison = comparison$comparison,
    dataset = comparison$dataset,
    statistic = res$statistic, #res$chi2,
    p.value = res$p.value #pchisq(res$chi2, res$df)
  ))
}

testResults$p.adj <- p.adjust(testResults$p.value, method = "fdr")


### PLOTTING ##################

partial.plot.1 <- ggplot(
  df1, aes(x = corr, color = covariates, fill = covariates)) +
  geom_density(position = "identity", alpha = 0.3) +
  #geom_freqpoly(position = "identity", binwidth = 0.05, size = 1.5) +
  xlab("|R|") + ylab("Density") +
  scale_x_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0)) +
  scale_color_manual("Event cross-correlations", values = c(colPal$Grey, colPal$MediumBlue, colPal$AccentYellow)) +
  scale_fill_manual("Event cross-correlations", values = c(colPal$Grey, colPal$MediumBlue, colPal$AccentYellow)) +
  ylim(0, 4) +
  theme_classic() +
  theme(text = element_text(size = 20), 
        axis.text = element_text(color = "black"),
        strip.text = element_text(size = 20),
        strip.background = element_rect(color="white", fill="white"),
        panel.spacing = unit(2, "lines")) +
  facet_wrap(~ dataset, ncol = 3) +
  geom_vline(data = vline1, aes(xintercept=median, color = covariates), linetype="dashed")

partial.plot.2 <- ggplot(
  df2, aes(x = corr, color = covariates, fill = covariates)) +
  geom_density(position = "identity", alpha = 0.3) +
  #geom_freqpoly(position = "identity", binwidth = 0.05, size = 1.5) +
  xlab("|R|") + ylab("Density") +
  scale_x_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0)) +
  scale_fill_manual("Event cross-correlations", values = c(colPal$Grey, colPal$Red)) +
  scale_color_manual("Event cross-correlations", values = c(colPal$Grey, colPal$Red)) +
  ylim(0, 4) +
  theme_classic() +
  theme(text = element_text(size = 20), 
        axis.text = element_text(color = "black"),
        strip.text = element_text(size = 20),
        strip.background = element_rect(color="white", fill="white"),
        panel.spacing = unit(2, "lines")) +
  facet_wrap(~ dataset, ncol = 3) +
  geom_vline(data = vline2, aes(xintercept=median, color = covariates), linetype="dashed")

### ARRANGE FIGURE ############

legend1 <- get_legend(partial.plot.1)
legend2 <- get_legend(partial.plot.2)

plot_grid(partial.plot.1 + theme(legend.position = "none"), legend1,
          partial.plot.2 + theme(legend.position = "none"), legend2,
          ncol = 2, scale = rep(0.95,4), rel_widths = c(0.7,0.3), 
          labels = "", label_size = 20)

ggsave2(filename = "results/Figure_5.pdf", 
        device = cairo_pdf, width = 12, height = 6, dpi = 300)

