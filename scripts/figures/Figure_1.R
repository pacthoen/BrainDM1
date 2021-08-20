##################
### Figure_1.R ###
##################

suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(stringr) )
suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(tidyr) )
suppressPackageStartupMessages( library(ggplot2) )
suppressPackageStartupMessages( library(VennDiagram) )
suppressPackageStartupMessages( library(cowplot) )
suppressPackageStartupMessages( library(gghighlight) )
suppressPackageStartupMessages( library(RColorBrewer) )
suppressPackageStartupMessages( library(ggrepel) )
suppressPackageStartupMessages( library(extrafont) )

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

colPal <- list(Blue = "#003366", 
               Red = "#E31B23", 
               MediumBlue = "#005CAB", 
               LightBlue = "#DCEEF3", 
               AccentYellow = "#FFC325", 
               Grey = hsv(0,0,0.5),
               CoolGrey = "#E6F1EE")

theme_set(theme_classic(base_family = "Arial"))

### LOAD DATA ############

sample_metadata <- fread("data/sample_metadata.csv")
high_conf_events <- fread("data/high_confidence_events.csv")
psi_comparison_results <- fread("data/BrainSpan_Otero_PSI_comparison_results.csv") 

### DATASET OVERVIEW #########

#---BRAINSPAN---#
BrainSpan.overview.plot <- sample_metadata %>%
  dplyr::filter(dataset == "BrainSpan") %>%
  dplyr::mutate(group = factor(group, levels = c("Prenatal", "Postnatal"))) %>%
  dplyr::select(subjectID, age_in_days, group, sex) %>% distinct() %>% 
  ggplot(aes(y = log(age_in_days), x = group, fill = sex)) +
  geom_dotplot(binaxis = "y", dotsize = 1.5,
               stackdir = "center", position = position_dodge(0.5)) +
  scale_y_continuous(breaks = c(log(90), log(280), log(3 * 365 + 280), 
                                log(10 * 365 + 280), log(20 * 365 + 280)), 
                     labels = c("3 Mos", "Birth", "3 Yrs", "10 Yrs", "20 Yrs")) +
  scale_fill_manual(values = c(colPal$MediumBlue, colPal$Grey)) +
  ggtitle("Developing brain (BrainSpan)") +
  xlab("") + ylab("log10(Age)") +
  theme_classic() +
  theme(text = element_text(size = 20),
        title = element_text(size = 15),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"))

#---OTERO---#
Otero.overview.plot <- sample_metadata %>%
  dplyr::filter(dataset == "Otero et al. (2021)") %>%
  dplyr::select(sampleID, sex, age, group) %>% distinct() %>%
  ggplot(aes(x = group, y = age, fill = sex)) +
  geom_dotplot(binaxis = "y", dotsize = 1.5,
               stackdir = "center", position = position_dodge(0.5)) +
  scale_fill_manual(values = c(colPal$MediumBlue, colPal$Grey)) +
  ggtitle("DM1/Unaffected (Otero et al., 2021)") +
  xlab("") + ylab("Age") +
  ylim(38, 80) +
  theme_classic() +
  theme(text = element_text(size = 20),
        title = element_text(size = 15),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"))

### MEAN PSI PLOTS ##############

highlight_events <- arrange(high_conf_events, dPSI_DM1_CTRL)
highlight_events <- rbind(head(highlight_events, 4), tail(highlight_events, 4))

#---BRAINSPAN---#
BrainSpan.psi.plot <- psi_comparison_results %>%
  filter(!is.na(`BrainSpan (Postnatal)`) & !is.na(`BrainSpan (Prenatal)`)) %>%
  ggplot(aes(y = `BrainSpan (Prenatal)`, x = `BrainSpan (Postnatal)`)) +
  geom_point(size = 0.8) +
  geom_point(data = filter(psi_comparison_results, event_id %in% high_conf_events$event_id), 
             shape = 15, size = 3, color = colPal$MediumBlue) +
  gghighlight(p.val_pre_post < 0.01 & abs(dPSI_pre_post) > 0.2) +
  geom_label_repel(aes(label = high_conf_events$SE_event_name[match(event_id, high_conf_events$event_id)]),
                   filter(psi_comparison_results, event_id %in% highlight_events$event_id),
                   point.padding = 0.5, size = 3) +
  ylab(expression("Prenatal mean" ~ psi)) + 
  xlab(expression("Postnatal mean" ~ psi)) + 
  ylim(0,1) + xlim(0,1) +
  theme_classic() +
  theme(text = element_text(size = 15)) 

#---Otero---#
Otero.psi.plot <- ggplot(
  filter(psi_comparison_results, !is.na(`Otero2021 (DM1)`) & !is.na(`Otero2021 (Unaffected)`)), 
  aes(y = `Otero2021 (DM1)`, x = `Otero2021 (Unaffected)`)) +
  geom_point(size = 0.8) +
  geom_point(data = filter(psi_comparison_results, event_id %in% high_conf_events$event_id), 
             shape = 15, size = 3, color = colPal$MediumBlue) +
  gghighlight(p.val_DM1_CTRL < 0.01 & abs(dPSI_DM1_CTRL) > 0.2) +
  geom_label_repel(aes(label = high_conf_events$SE_event_name[match(event_id, high_conf_events$event_id)]),
                   filter(psi_comparison_results, event_id %in% highlight_events$event_id),
                   point.padding = 0.5, size = 3) +
  ylab(expression("DM1 mean" ~ psi)) + 
  xlab(expression("Unaffected mean" ~ psi)) + 
  ylim(0,1) + xlim(0,1) +
  theme_classic() +
  theme(text = element_text(size = 15))

### BAR PLOT ##########

sigResults <- psi_comparison_results %>%
  dplyr::select(event_id, p.val_pre_post, p.val_DM1_CTRL, 
                dPSI_pre_post, dPSI_DM1_CTRL) %>%
  mutate(sig_in_both = p.val_pre_post < 0.01 & p.val_DM1_CTRL < 0.01 &
           abs(dPSI_pre_post) > 0.2 & abs(dPSI_DM1_CTRL) > 0.2,
         sig_in_one = (p.val_pre_post < 0.01 & abs(dPSI_pre_post) > 0.2) | 
           (p.val_DM1_CTRL < 0.01 & abs(dPSI_DM1_CTRL) > 0.2),
         direction = ifelse(sign(dPSI_pre_post) == sign(dPSI_DM1_CTRL), "Same", "Opposite")) %>%
  filter(sig_in_one == TRUE) %>%
  mutate(significant = ifelse(sig_in_both & !is.na(sig_in_both), "In both datasets", "In one dataset")) %>%
  dplyr::select(-sig_in_both, -sig_in_one)

sigResults$direction[is.na(sigResults$direction)] <- "Unknown"
sigResults$direction <- factor(sigResults$direction, levels = c("Same", "Opposite", "Unknown"))

bar.plot.withNA <- ggplot(
  sigResults, aes(x = significant, fill = direction)) +
  geom_bar(width = 0.5) +
  geom_text(stat = "count", 
            aes(label=..count..), 
            position = position_stack(vjust = 0.5), 
            color="white") +
  scale_fill_manual(values = c(colPal$MediumBlue, colPal$Red, colPal$Grey)) +
  ylab("Number of events") + xlab("Significant") + labs(fill = "Splicing direction") +
  theme_classic() +
  theme(text = element_text(size = 15),
        strip.background =element_rect(color="white", fill="white"),
        strip.placement = "outside")

sigResults <- sigResults[!sigResults$direction == "Unknown", ]

bar.plot.withoutNA <- ggplot(
  sigResults, aes(x = significant, fill = direction)) +
  geom_bar(width = 0.5) +
  geom_text(stat = "count", 
            aes(label=..count..), 
            position = position_stack(vjust = 0.5), 
            color="white") +
  scale_fill_manual(values = c(colPal$MediumBlue, colPal$Grey)) +
  ylab("Number of events") + xlab("Significant") + labs(fill = "Splicing direction") +
  theme_classic() +
  theme(text = element_text(size = 15),
        strip.background =element_rect(color="white", fill="white"),
        strip.placement = "outside")

### VENN DIAGRAM ########

# create venn diagram
venn.plot.withNA <- venn.diagram(x = list(na.omit(psi_comparison_results$event_id[psi_comparison_results$p.val_pre_post < 0.01 & 
                                                                            abs(psi_comparison_results$dPSI_pre_post) > 0.2]),
                                   na.omit(psi_comparison_results$event_id[psi_comparison_results$p.val_DM1_CTRL < 0.01 & 
                                                                            abs(psi_comparison_results$dPSI_DM1_CTRL) > 0.2])),
                          category.names = c("Events related to brain development", "Events related to DM1"),
                          fill = c(colPal$MediumBlue, colPal$Red), #RColorBrewer::brewer.pal(3, "RdYlBu")[c(1,3)],
                          cat.cex = 1.5, cat.pos = c(0,170), cat.dist = c(0.1,0.1), lwd = 1, #lty = "blank", 
                          fontfamily = "serif", cat.fontfamily = "serif",
                          filename = NULL)

# create venn diagram but remove all missing values
psi_comparison_results <- na.omit(psi_comparison_results)
venn.plot.withoutNA <- venn.diagram(x = list(psi_comparison_results$event_id[psi_comparison_results$p.val_pre_post < 0.01 & 
                                                                                       abs(psi_comparison_results$dPSI_pre_post) > 0.2],
                                              psi_comparison_results$event_id[psi_comparison_results$p.val_DM1_CTRL < 0.01 & 
                                                                                       abs(psi_comparison_results$dPSI_DM1_CTRL) > 0.2]),
                                     category.names = c("Events related to brain development", "Events related to DM1"),
                                     fill = c(colPal$MediumBlue, colPal$Red), #RColorBrewer::brewer.pal(3, "RdYlBu")[c(1,3)],
                                     cat.cex = 1.5, cat.pos = c(0,170), cat.dist = c(0.1,0.1), lwd = 1, #lty = "blank", 
                                     fontfamily = "serif", cat.fontfamily = "serif",
                                     filename = NULL)

### ARRANGE FIGURE ###########

plot_grid(NULL, 
          plot_grid(BrainSpan.overview.plot,
                    Otero.overview.plot, 
                    ncol = 2, labels = c("A", "B"), label_size = 20, scale = c(0.9,0.9)),
          plot_grid(BrainSpan.psi.plot,
                    Otero.psi.plot,
                    ncol = 2, labels = c("C", "D"), label_size = 20, scale = c(0.9,0.9)),
          plot_grid(venn.plot.withNA,
                    bar.plot.withNA,
                    ncol = 2, labels = c("E", "F"), label_size = 20, scale = c(0.9,0.9,0.9)),
          ncol = 1, labels = "Figure 1", label_size = 20, rel_heights = c(0.1, 0.5, 1, 0.50))

ggsave2(filename = "results/Figure_1.pdf", width = 12, height = 12, dpi = 300)

