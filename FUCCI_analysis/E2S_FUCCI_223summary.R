### R script for data analysis of live-cell imaging FUCCI data. 
### Input should be a csv file with single cell data created by running 
### CellProfiler with TrackObjects module and CPtrackR to fix the output.
### 
### The input data should at least contain the following columns:
### - Nuclei_TrackObjects_Label_30
### - Nuclei_TrackObjects_ParentObjectNumber_30
### - Nuclei_Number_Object_Number
### - plateID 
### - locationID
### - treatment
### - dose_uM
### - timeID
### - timeAfterExposure
### - 
### 
### Author: Muriel Heldring
### Email: m.m.heldring@lacdr.leidenuniv.nl
### Date created: 2021/10/14
### 

## Clear environment
rm(list = ls())
gc()
getwd()

### USER DEFINED VARIABLES ###

DATE <- "20240322"

PATH_TO_INPUT <- "/data/muriel/Projects/E2S/Analysis/00_Imaging/Exp007_FUCCI_Britt/Input/"
OUTPUT_PATH <- "/data/muriel/Projects/E2S/Analysis/00_Imaging/Exp007_FUCCI_Britt/Output/"
READOUT <- "hoechst" # "hoechst" "H2B"

FIRST_TIMEID <- 2 # TimeID for the first measurement. Usually, due to frame shifts, I tell CellProfiler to 
# skip the first image, so then the first timeID in the data table should be 2.
TIME_STEP <- 0.5 # Time interval between measurements in hours


### END USER DEFINED VARIABLES ###

## Package installation (if necessary) and loading
# install.packages("tidyverse")
library(pacman)
p_load(pheatmap)
p_load(ggplot2)
p_load(grid)
p_load("viridis")
p_load(DBI)
p_load(tidyverse)
p_load(ggnewscale)
p_load(multidplyr)
p_load(ggpubr)
p_load(rstatix)

# Functions
getPCN <- function(x) {
  # Get the unique values and how often they are repeated
  y <- rle(x)
  # Bind the values in blocks of three with a dash in between
  sequences <- paste0(y$values,"-",c(y$values[-1],NA),"-",c(y$values[-1:-2],NA,NA))
  # Repeat every sequence of 3 according to the length of the chunk
  out <- rep(c(NA,sequences[-length(sequences)]),y$length)
  out
}

# Create output folder
dirName <- paste0(OUTPUT_PATH,READOUT,"_segmentation/",DATE,"_Figures")
if (!dir.exists(dirName)){
  dir.create(dirName)
  print("Created directory!")}


### Start analysis

# Load data experiment 223
load(file = paste0(OUTPUT_PATH,"/",READOUT,"_segmentation/RData/","223","_",READOUT,"_tracks_unique.R"))
tracks_unique <- tracks_unique %>% add_column(Experiment = "223") %>% relocate(Experiment)

# Refactor condition
tracks_unique <- tracks_unique %>% 
  mutate(condition = ifelse(condition == "10pM", "10 pM",
                            ifelse(condition == "1nM", "1 nM",
                                   ifelse(condition == "100nM", "100 nM","DMSO"))),
         condition = factor(condition, levels = c("DMSO", "10 pM", "1 nM","100 nM")))

# Change the time unit to hr^-1 instead of 0.5 hr^-1
tracks_unique <- tracks_unique %>% group_by(Experiment,well_name_p, cid, TrackID) %>%
  mutate(PhaseLength = TIME_STEP * PhaseLength)

# Filter only the full length phases
tracks_unique <- tracks_unique %>% 
  filter(PrevCurNext %in% c("G1-G1/S-G2", "G1/S-G2-G1","G2-G1-G1/S"))

# Select only the FUCCI data and remove the FUCCI/WT mix
tracks_unique <- tracks_unique %>% filter(grepl("02", well_name_p))
unique(tracks_unique %>% pull(well_name_p))

# Get phase statistics
phase_stats <- tracks_unique %>% ungroup() %>% 
  group_by(Experiment,cid, well_name_p, well_name, Phase_corrected, 
           treatment,condition,dose_uM,cell_line, TrackID) %>%
  summarise(PhaseLength = unique(PhaseLength)) %>%
  ungroup() 


ggplot(phase_stats %>% filter(cell_line == "FUCCI"),
       aes(x = condition, y = PhaseLength, fill = Phase_corrected)) +
  geom_point(shape = 21, size=1, position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), alpha = 0.25, color = "black") +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  ylab("Phase duration (h)") + xlab("E2 concentration") +
  scale_fill_manual(values = c("G1" = "red","G1/S" = "gold", "S-G2-M" = "green"),
                    name = "Phase") +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 12),
        legend.position = "top") +
  facet_grid(. ~ Phase_corrected)
ggsave(paste0(OUTPUT_PATH,READOUT,"_segmentation/",DATE,"_Figures/Fig4F_phase_length_stats_FUCCI223data.pdf"), width = 8.5, height = 3.5)

dataStats <- phase_stats %>% filter(cell_line == "FUCCI",
                       Phase_corrected == "G1") %>% group_by(condition) %>%
  summarise(medianPhaseDur = median(PhaseLength))

# Check for normality ## VIOLATED!
d <- phase_stats %>% filter(cell_line == "FUCCI") %>% mutate(phase_condition = paste0(Phase_corrected,"_",condition))
# Build the linear model
model  <- lm(PhaseLength ~ phase_condition, data = d)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))

# Check normality assumption by group
d %>%
  group_by(phase_condition) %>%
  shapiro_test(PhaseLength)

# Non-parametric test
# Welch One way ANOVA test
res.aov2 <- d %>% welch_anova_test(PhaseLength ~ phase_condition)
# Pairwise comparisons (Games-Howell)
pwc2 <- d %>% games_howell_test(PhaseLength ~ phase_condition)
# Visualization: box plots with p-values
pwc2 <- pwc2 %>% add_xy_position(x = "phase_condition", step.increase = 1)
# Filter out the relevant comparisons (within phases)
pwc2 <- pwc2 %>% mutate(group1_cond = str_split(group1,"_", simplify = T)[,1],
                        group2_cond = str_split(group2,"_", simplify = T)[,1]) %>%
  filter(group1_cond == group2_cond)
ggboxplot(d, x = "phase_condition", y = "PhaseLength") +
  stat_pvalue_manual(pwc2, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov2, detailed = TRUE),
    caption = get_pwc_label(pwc2)
  )

