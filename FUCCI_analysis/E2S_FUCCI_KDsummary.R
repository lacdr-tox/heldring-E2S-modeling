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

DATE <- "20221216"

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

# Make a cluster
if(!exists("cluster")) {
  cluster <- new_cluster(35)
  cluster_library(cluster,"tidyverse")
  cluster_send(cluster,
               ma <- function(x, n = 5, side = 2){
                 stats::filter(x, rep(1 / n, n), sides = side)
               })
  cluster_send(cluster,
               getPCN <- function(x) {
                 # Get the unique values and how often they are repeated
                 y <- rle(x)
                 # Bind the values in blocks of three with a dash in between
                 sequences <- paste0(y$values,"-",c(y$values[-1],NA),"-",c(y$values[-1:-2],NA,NA))
                 # Repeat every sequence of 3 according to the length of the chunk
                 out <- rep(c(NA,sequences[-length(sequences)]),y$length)
                 out
               })
  cluster_send(cluster,
               getPrev <- function(x) {
                 # Get the unique values and how often they are repeated
                 y <- rle(x)
                 # Bind the values in blocks of three with a dash in between
                 prev <- y$values
                 # Repeat every sequence of 3 according to the length of the chunk
                 out <- rep(c(NA,prev[-length(prev)]),y$length)
                 out
               })}

### Start analysis

# Load data experiment 230
load(file = paste0(OUTPUT_PATH,"/",READOUT,"_segmentation/RData/","230","_",READOUT,"_tracks_unique.R"))
tracks_230 <- tracks_unique
tracks_230 <- tracks_230 %>% add_column(Experiment = "230") %>% relocate(Experiment)

# Load data experiment 233
load(file = paste0(OUTPUT_PATH,"/",READOUT,"_segmentation/RData/","233","_",READOUT,"_tracks_unique.R"))
tracks_233 <- tracks_unique
tracks_233 <- tracks_233 %>% add_column(Experiment = "233") %>% relocate(Experiment)

tracks_unique <- bind_rows(tracks_230, tracks_233)

# Rename treatment
tracks_unique <- tracks_unique %>% 
  mutate(treatment = ifelse(treatment == "siKP", "siControl",treatment))

# # Update PrevCurNext column
# tracks_unique <- tracks_unique %>% ungroup() %>% 
#   group_by(TrackNameUnique, well_name_p) %>%
#   partition(cluster) %>% 
#   mutate(PrevCurNext = getPCN(Phase_corrected)) %>%
#   collect()

# Rename phases
tracks_unique <- tracks_unique %>% 
  mutate(Phase_corrected = ifelse(Phase_corrected == "G2", "S-G2-M", Phase_corrected))

# # From the original tracks, delete duplicates but do not split mother and daughter cells
# # Collapse duplicate rows into one
# tracks_unique <- track_filtered_pks %>% ungroup() %>%
#   filter(!duplicated(track_filtered_pks %>% ungroup() %>%
#                        select(Experiment,well_name_p,cid,timeID,Cdt1,Geminin, meanHoechst, nucleiSize)))
# # Remove tracks of length less than 5
# tracks_unique <- tracks_unique %>%
#   group_by(Experiment,well_name_p, cid, TrackID) %>%
#   mutate(track_length = n()) %>%
#   filter(track_length > 5) %>% ungroup()
# # # Recalculate phase lengths
# tracks_unique <- tracks_unique %>% group_by(Experiment,well_name_p, cid, TrackID) %>%
#   mutate(PhaseLength = unlist(sapply(rle(Phase_corrected)$length,function(x){rep(x,x)}, simplify = T)))
# # Change the time unit to hr^-1 instead of 0.5 hr^-1
# tracks_unique <- tracks_unique %>% group_by(Experiment,well_name_p, cid, TrackID) %>%
#   mutate(PhaseLength = TIME_STEP * PhaseLength)
# #tracks_unique <- tracks_unique %>% filter(well_name_p == "C02_1")

# # Filter only the full length phases
tracks_unique <- tracks_unique %>% 
  filter(PrevCurNext %in% c("G1-G1/S-G2", "G1/S-G2-G1","G2-G1-G1/S"))

# Select only the FUCCI data and remove the FUCCI/WT mix
tracks_unique <- tracks_unique %>% filter(grepl("B", well_name_p) | grepl("C", well_name_p))
unique(tracks_unique %>% pull(well_name_p))

# Get phase statistics
phase_stats <- tracks_unique %>% ungroup() %>% 
  group_by(Experiment,cid, well_name_p, well_name, Phase_corrected, 
           treatment,condition,dose_uM,cell_line, TrackID) %>%
  summarise(PhaseLength = unique(PhaseLength)) %>%
  ungroup() 

ggplot(phase_stats %>% filter(cell_line == "FUCCI",
                              treatment %in% c("mock", "siGREB1", "siPGR", "siTFF1"),
                              condition == "100nM"),
       aes(x = treatment, y = PhaseLength, fill = Phase_corrected)) + 
  geom_point(shape = 21, size=1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), alpha = 0.25, color = "black") +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  ylab("Phase duration (h)") + xlab("Condition") +
  scale_fill_manual(values = c("G1" = "red","G1/S" = "gold", "S-G2-M" = "green"), 
                    name = "Phase") + 
  theme_classic() + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 12))
ggsave(paste0(OUTPUT_PATH,READOUT,"_segmentation/",DATE,"_Figures/Fig5A_phase_length_stats_KDdata.pdf"), width = 8, height = 3)


# # Do statistical test
# res.aov2 <- aov(PhaseLength ~ treatment + Phase_corrected, data = phase_stats %>% filter(cell_line == "FUCCI",
#                                                                                          condition == "100nM"))
# print(summary(res.aov2))
# print(TukeyHSD(res.aov2, which = "treatment"))
# 
# res.aov3 <- aov(PhaseLength ~ treatment + Phase_corrected, data = phase_stats %>% filter(cell_line == "FUCCI",
#                                                                                          condition == "100nM"))
# print(summary(res.aov3))
# print(TukeyHSD(res.aov3, which = "treatment"))
# 

# Do pairwise t-test to check if phase lengths differ between treatment conditions
# Do for G1
#tmp <- phase_stats %>% filter(cell_line == "FUCCI")
G1_data <- phase_stats %>% filter(cell_line == "FUCCI", 
                                  Phase_corrected == "G1",
                                  treatment %in% c("mock", "siGREB1", "siPGR", "siTFF1"),
                                  condition == "100nM")
print(pairwise.t.test(G1_data %>% pull (PhaseLength), G1_data %>% pull (treatment),
                      p.adjust.method = "bonferroni"))

# Do for G1/S
G1S_data <- phase_stats %>% filter(cell_line == "FUCCI", 
                                   Phase_corrected == "G1/S",
                                   treatment %in% c("mock", "siGREB1", "siPGR", "siTFF1"),
                                   condition == "100nM")
print(pairwise.t.test(G1S_data %>% pull (PhaseLength), G1S_data %>% pull (treatment),
                      p.adjust.method = "bonferroni"))

# Do for G2
G2_data <- phase_stats %>% filter(cell_line == "FUCCI", 
                                  Phase_corrected == "S-G2-M",
                                  treatment %in% c("mock", "siGREB1", "siPGR", "siTFF1"),
                                  condition == "100nM")
print(pairwise.t.test(G2_data %>% pull (PhaseLength), G2_data %>% pull (treatment),
                      p.adjust.method = "bonferroni"))

