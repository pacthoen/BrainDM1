###################
### Figure_S3.R ###
###################

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

### BARPLOT RANKING SIG EVENTS ##############

# prepare data for plotting
data <- psi_comparison_results %>%
  dplyr::select(event_id, dPSI_pre_post, p.val_pre_post, dPSI_DM1_CTRL, p.val_DM1_CTRL) %>%
  dplyr::mutate(
    BrainSpan_sigEvent = ifelse(abs(dPSI_pre_post) > 0.2 & p.val_pre_post < 0.01, TRUE, FALSE),
    Otero_sigEvent = ifelse(abs(dPSI_DM1_CTRL) > 0.2 & p.val_DM1_CTRL < 0.01, TRUE, FALSE),
    SE_event_name = high_conf_events$SE_event_name[match(event_id, high_conf_events$event_id)]
  ) %>%
  dplyr::filter(BrainSpan_sigEvent | Otero_sigEvent)

data$BrainSpan_sigEvent[is.na(data$BrainSpan_sigEvent)] <- FALSE
data$Otero_sigEvent[is.na(data$Otero_sigEvent)] <- FALSE
data$eventOfInterest <- ifelse(is.na(data$SE_event_name), FALSE, TRUE) 
data$SE_event_name[is.na(data$SE_event_name)] <- ""

#---BRAINSPAN---#
data %>%
  dplyr::filter(!is.na(dPSI_pre_post), BrainSpan_sigEvent) %>%
  dplyr::mutate(abs_dPSI_pre_post = abs(dPSI_pre_post)) %>%
  ggplot() +
  geom_bar(
    aes(x = reorder(event_id, abs_dPSI_pre_post), 
        y = abs_dPSI_pre_post,
        fill = eventOfInterest), 
    stat = "identity", width = 0.7) +
  scale_fill_manual(values = c(colPal$Grey, colPal$Red)) +
  scale_y_continuous(limits=c(0.0,1.0), breaks = seq(0,1,by=0.25)) +
  ylab("| dPSI |") + xlab("Significant events") +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 18),
    #axis.title.y = element_text(angle = 0, vjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )

#---Otero---#
p <- data %>%
  dplyr::filter(!is.na(dPSI_DM1_CTRL), Otero_sigEvent) %>%
  dplyr::mutate(abs_dPSI_DM1_CTRL = abs(dPSI_DM1_CTRL),
                eventOfInterest = factor(ifelse(eventOfInterest == TRUE, "Yes", "No"), levels = c("Yes", "No"))) %>%  
  ggplot(aes(x = reorder(event_id, abs_dPSI_DM1_CTRL),
             y = abs_dPSI_DM1_CTRL)) +
  geom_bar(
    aes(fill = eventOfInterest), 
    stat = "identity", width = 0.6) +
  geom_text(aes(label = ifelse(eventOfInterest == "Yes", "x", "")),
    position=position_dodge(width=0.9), vjust=-0.25) +
  scale_fill_manual(values = c(colPal$MediumBlue, colPal$Grey)) +
  scale_y_continuous(limits=c(0.0,1.0), breaks = seq(0,1,by=0.25)) +
  ylab("| dPSI |") + xlab("Sig. events in DM1/Unaffected") +
  guides(fill=guide_legend(title="Developmentally regulated")) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 18),
    #axis.title.y = element_text(angle = 0, vjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )  

l <- get_legend(p)

plot_grid(
  l,
  p + theme(legend.position = "none"),
  scale = c(0.9, 0.9), rel_heights = c(0.2,0.8), labels = "", nrow = 2
)

ggsave2(filename = "results/Figure_S3.pdf", width = 10, height = 4, dpi = 300)

### significance testing

x <- data %>%
  dplyr::filter(!is.na(dPSI_DM1_CTRL), Otero_sigEvent) %>%
  dplyr::mutate(abs_dPSI_DM1_CTRL = abs(dPSI_DM1_CTRL))

testResults <- wilcox.test(x$abs_dPSI_DM1_CTRL[x$eventOfInterest == TRUE], 
                           x$abs_dPSI_DM1_CTRL[x$eventOfInterest == FALSE],
                           alternative = "greater")

cat("Wilcoxon test - p-value:", testResults$p.value, "\n")
