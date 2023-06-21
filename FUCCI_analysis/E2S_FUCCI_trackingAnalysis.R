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
EXP <- "230"
READOUT <- "hoechst" # "hoechst" "H2B"

FIRST_TIMEID <- 2 # TimeID for the first measurement. Usually, due to frame shifts, I tell CellProfiler to 
# skip the first image, so then the first timeID in the data table should be 2.
TIME_STEP <- 0.5 # Time interval between measurements in hours


### END USER DEFINED VARIABLES ###

## Package installation (if necessary) and loading
# install.packages("tidyverse")
library(pacman)
p_load(multidplyr)
p_load(pheatmap)
p_load(ggplot2)
p_load(grid)
p_load("viridis")
p_load(DBI)
p_load(tidyverse)
p_load(ggnewscale)
p_load(readr)

# FUNCTIONS

make_full_tracks <- function(data, ftid){
  ### Fill timepoint t1
  tracks_t1 <- data %>% filter(timeID == ftid) %>% pull(alt_uid)
  df <- tibble(TrackName = tracks_t1,
               !!paste0("t",ftid) := tracks_t1)
  tracks_prev <- tracks_t1
  prev_col_name <- paste0("t",ftid)
  
  ### Fill following timepoints
  # Get the number of timepoints
  tps <- unique(data %>% pull(timeID))
  for (i in tps[2:length(tps)]) {
    tracks_next <- data %>% filter(timeID == i) %>% pull(alt_uid)
    next_col_name <- paste0("t",i)
    existing <- tracks_next[which(tracks_next %in% tracks_prev)]
    branch <- tracks_next[which(!tracks_next %in% tracks_prev & grepl("\\.",tracks_next))]
    new <- tracks_next[which(!tracks_next %in% tracks_prev & !grepl("\\.",tracks_next))]
    
    # Check whether all tracks are used and unique
    if (!all(sort(c(new,existing,branch)) == sort(tracks_next))) {warning("Not all new track points are taken along or some are double.")}
    
    # Make tibbles from existing, branch and new tracks
    existing <- tibble(!!prev_col_name := existing,
                       !!next_col_name := existing)
    new <- tibble(TrackName = new,
                  !!next_col_name := new)
    branch <- tibble(!!prev_col_name := sapply(branch,function(x){ 
      paste0(str_split(x,"\\.")[[1]][1:length(str_split(x,"\\.")[[1]])-1],collapse = ".")}),
      !!next_col_name := branch)
    
    # Add the existing track points to t2
    if (!dim(branch)[1] == 0) {
      df <- left_join(df, bind_rows(branch,existing), by = prev_col_name)
    } else {
      df <- left_join(df, existing, by = prev_col_name)
    }
    df <- bind_rows(df, new)
    
    # Rename previous colname and previous track
    tracks_prev <- tracks_next
    prev_col_name <- next_col_name
  }
  
  # Change double track names
  df <- df %>% mutate(TrackNameUnique = seq(1:nrow(df)))
  
  return(df)
}

# Calculate moving average
ma <- function(x, n = 5, side = 2){
  stats::filter(x, rep(1 / n, n), sides = side)
}


# Make track df
get_tracks <- function(well_id, data, ftid){
  # Get subset of data for one well location
  data_tmp <- data %>% filter(well_name_p == well_id)
  
  # Initiate temporary output df
  make_full_tracks(data_tmp, ftid = ftid) %>% mutate(well_name_p = well_id) %>% relocate(c("TrackName","TrackNameUnique","well_name_p"))
  #return(data_tmp)
}

getPCN <- function(x) {
  # Get the unique values and how often they are repeated
  y <- rle(x)
  # Bind the values in blocks of three with a dash in between
  sequences <- paste0(y$values,"-",c(y$values[-1],NA),"-",c(y$values[-1:-2],NA,NA))
  # Repeat every sequence of 3 according to the length of the chunk
  out <- rep(c(NA,sequences[-length(sequences)]),y$length)
  out
}

getPrev <- function(x) {
  # Get the unique values and how often they are repeated
  y <- rle(x)
  # Bind the values in blocks of three with a dash in between
  prev <- y$values
  # Repeat every sequence of 3 according to the length of the chunk
  out <- rep(c(NA,prev[-length(prev)]),y$length)
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

# Get data frames
files <- dir(PATH_TO_INPUT)
track_data_files <- files[grepl("fixed_tracks",files)]
layout_metadata_files <- files[grepl("layout",files)]

input_file <- track_data_files[grepl(EXP,track_data_files) &
                                 grepl(READOUT,track_data_files)]
layout <- layout_metadata_files[grepl(EXP,layout_metadata_files)]

## Start analysis ##

# Load the data
data_orig <- read_csv(paste0(PATH_TO_INPUT,input_file), 
                      col_types = cols(alt_uid = col_character()))
layout_orig <- read_csv(paste0(PATH_TO_INPUT,layout))

# Remove ...1 column
layout <- layout_orig %>% select(-c(X1))

# Select and rename columns to be kept
colnames(data_orig)
data <- data_orig %>% dplyr::rename(Cdt1 = nuclei_Intensity_MeanIntensity_image_cdt1,
                                    Geminin = nuclei_Intensity_MeanIntensity_image_geminin,
                                    meanHoechst = nuclei_Intensity_MeanIntensity_image_hoechst,
                                    integratedHoechst = nuclei_Intensity_IntegratedIntensity_image_hoechst,
                                    meanH2B = nuclei_Intensity_MeanIntensity_image_h2b,
                                    integratedH2B = nuclei_Intensity_IntegratedIntensity_image_h2b,
                                    nucleiSize = nuclei_AreaShape_Area,
                                    plateID = Image_Metadata_PlateID,
                                    well = Image_Metadata_Well,
                                    well_name = Image_Metadata_well_name,
                                    well_name_p = Image_Metadata_well_name_p,
                                    position = Image_Metadata_imageNr) %>%
  select(c("position","ImageNumber","cid","uid","alt_uid",
           "plateID","well","well_name_p","well_name",
           "nuclei_Number_Object_Number", 
           "Cdt1","Geminin","meanHoechst","integratedHoechst",
           "meanH2B","integratedH2B","nucleiSize"))

# Merge data with layout metadata
data <- left_join(data,layout, by = "well_name") %>% 
  relocate(c("treatment","condition","dose_uM","cell_line",
             "cid","uid","alt_uid"))

# Do min-max normalisation on data
data <- data %>% mutate(Cdt1 = ((Cdt1 - min(Cdt1, na.rm = T))/
                                  (max(Cdt1, na.rm = T) - min(Cdt1, na.rm = T))),
                        Geminin = ((Geminin - min(Geminin, na.rm = T))/
                                     (max(Geminin, na.rm = T) - min(Geminin, na.rm = T))),
                        meanHoechst = ((meanHoechst - min(meanHoechst, na.rm = T))/
                                         (max(meanHoechst, na.rm = T) - min(meanHoechst, na.rm = T))),
                        integratedHoechst = ((integratedHoechst - min(integratedHoechst, na.rm = T))/
                                               (max(integratedHoechst, na.rm = T) - min(integratedHoechst, na.rm = T))),
                        meanH2B = ((meanH2B - min(meanH2B, na.rm = T))/
                                     (max(meanH2B, na.rm = T) - min(meanH2B, na.rm = T))),
                        integratedH2B = ((integratedH2B - min(integratedH2B, na.rm = T))/
                                           (max(integratedH2B, na.rm = T) - min(integratedH2B, na.rm = T))),
                        nucleiSize = ((nucleiSize - min(nucleiSize, na.rm = T))/
                                        (max(nucleiSize, na.rm = T) - min(nucleiSize, na.rm = T))))

# Reset ImageNumber to timeID
data <- data %>% group_by(well_name_p,treatment,condition,dose_uM,cell_line) %>% 
  mutate(timeID = ImageNumber - min(ImageNumber)+FIRST_TIMEID) %>% 
  ungroup %>% relocate(c("well_name","well_name_p","timeID"))

# Remove quickly appearing and disappearing cells, based on the uid:
# If the uid only exists for 5 frames or less, this is removed
# Only remove track dead ends
data <- data %>% group_by(well_name_p, uid) %>% 
  mutate(fragment_length = length(uid)) %>% ungroup() %>%
  group_by(well_name_p) %>% 
  partition(cluster) %>% 
  mutate(has_daughter = sapply(alt_uid, 
                               function(x){ifelse(paste0(x,".1") %in% alt_uid |
                                                    paste0(x,".2") %in% alt_uid,T,F)})) %>% 
  filter(!(fragment_length < 5 & !has_daughter)) %>%
  collect()

# Make full tracks per well_id
wellIDs <- unique(data %>% pull(well_name_p))
names(wellIDs) <- wellIDs
tracks_df_list <- lapply(wellIDs, get_tracks, data, ftid = FIRST_TIMEID)
names(tracks_df_list)

# Bind the track data frames
tracks_df <- tracks_df_list[[1]]
for (i in seq(2,length(tracks_df_list))) {
  tracks_df <- bind_rows(tracks_df,tracks_df_list[[i]])
}

# Collect the GFP expression info per track
df <- pivot_longer(tracks_df,cols = colnames(tracks_df)[grepl("t[1-9]",colnames(tracks_df))], 
                   names_to = "time",values_to = "alt_uid") %>% mutate(timeID = as.numeric(substring(time, 2)))
track_info <- left_join(df, data, by = c("timeID","alt_uid","well_name_p"))
track_info <- track_info %>% mutate(time_h = (timeID * TIME_STEP) - TIME_STEP)

# Calculate moving average
track_info <- track_info %>% group_by(well_name_p, TrackNameUnique) %>% 
  partition(cluster) %>% 
  mutate(ma_Cdt1 = ma(Cdt1, n = 10),
         ma_Geminin = ma(Geminin, n = 10),
         ma_meanHoechst = ma(meanHoechst, n = 10),
         ma_integratedHoechst = ma(integratedHoechst, n = 10),
         ma_meanH2B = ma(meanH2B, n = 10),
         ma_integratedH2B = ma(integratedH2B, n = 10),
         ma_nucleiSize = ma(nucleiSize, n = 10),
         track_break = ifelse(c(FALSE, alt_uid[-1L]!= alt_uid[-length(alt_uid)]),time_h,NA)) %>% 
  collect()

save(track_info, file = paste0(OUTPUT_PATH,"/",READOUT,"_segmentation/RData/",EXP,"_",READOUT,"_track_info.R"))
load(file = paste0(OUTPUT_PATH,"/",READOUT,"_segmentation/RData/",EXP,"_",READOUT,"_track_info.R"))

# Filter the data to get most relevant tracks
# Tracks of relevance are: 
# - the ones that are at least 30 hours long
track_filtered <- track_info %>% group_by(TrackNameUnique, well_name_p) %>% 
  mutate(track_length = sum(!is.na(meanHoechst))) %>% 
  ungroup %>% 
  filter(track_length > 30,
         !is.na(cid))

# Get the mitosis moments based on:
# - cell size
# - Hoechst intensity
# - Geminin drop

# Define peaks, dips and transitions
track_filtered_pks <- track_filtered %>% group_by(TrackNameUnique, well_name_p) %>% 
  partition(cluster) %>% 
  mutate(diff_hoechst = ma(abs(meanHoechst - ma_meanHoechst), n = 2),
         diff_size = c(NA,diff(nucleiSize)),
         diff_geminin = c(NA,diff(ma(Geminin, n = 2))),
         diff_cdt1 = c(NA,ma(diff(ma_Cdt1)))) %>%
  mutate(gem_over_cdt1 = ifelse(ma_Geminin == 0, -100,
                                ifelse(ma_Cdt1 < 0.0001 & ma_Geminin < 0.0001, -100,
                                       ifelse(ma_Cdt1 == 0 & ma_Geminin > 0.0001, 100,
                                              log(ma_Geminin/(ma_Cdt1))))),
         peaks_Hoechst_binary = diff_hoechst > 0.03,
         peaks_Hoechst_time = ifelse(peaks_Hoechst_binary,time_h,NA),
         peak_size_binary = diff_size > 0.05,
         peak_size_binshift = c(NA,peak_size_binary[-length(peak_size_binary)]),
         dip_size_binary = diff_size < -0.05,
         peaks_size_time = ifelse(dip_size_binary & (!peak_size_binshift),time_h,NA),
         dip_geminin_binary = diff_geminin < -0.5*sd(Geminin, na.rm = T),
         peaks_geminin_time = ifelse(dip_geminin_binary,time_h,NA),
         cdt1_decrease = diff_cdt1 < -0.15*sd(Cdt1, na.rm = T)
         ) %>%
  collect()

track_filtered_pks <- track_filtered_pks %>% group_by(TrackNameUnique, well_name_p, time_h) %>%
  partition(cluster) %>%
  mutate(Consistency = sum(c(dip_geminin_binary,(dip_size_binary & !peak_size_binshift),peaks_Hoechst_binary), na.rm = T),
         Break = dip_geminin_binary & Consistency > 0,
         Division = ifelse(dip_geminin_binary & Consistency == 3,"High",
                           ifelse(dip_geminin_binary & Consistency == 2,"Medium",
                                  ifelse(dip_geminin_binary & Consistency == 1,"Low",NA)))) %>%
  collect()

track_filtered_pks <- track_filtered_pks %>% ungroup() %>% group_by(TrackNameUnique, well_name_p) %>%
  partition(cluster) %>% 
  mutate(
    # If multiple breaks are identified, only take the last break as breaking point
    Break_filtered = ifelse(c(diff(Break),NA) == -1,T,F),
    Division_time = ifelse(Break_filtered, time_h,NA),
    Cycle_length = replace(Division_time,which(!is.na(Division_time)),c(NA,diff(na.omit(Division_time)))),
    Phase = ifelse(gem_over_cdt1 == -100, "G1",
                        ifelse(gem_over_cdt1 == 100, "G2",
                               ifelse(abs(gem_over_cdt1) < 2,"G1/S",
                                      ifelse(ma_Geminin > ma_Cdt1, "G2", "G1"))))) %>%
  select(-c(peaks_Hoechst_binary,peak_size_binary,peak_size_binshift,
            diff_hoechst,diff_size,diff_geminin,diff_cdt1)) %>%
  collect()

# Filter only rows with less than 5 breaks
track_filtered_pks <- track_filtered_pks %>% 
  filter(sum(Break_filtered, na.rm = T) < 5) %>% 
  ungroup() %>% group_by(cid, well_name_p) %>%
  mutate(TrackID = TrackNameUnique - min(TrackNameUnique),
         TrackCount = length(unique(TrackID))) %>%
  arrange(well_name_p, TrackNameUnique, timeID) %>% ungroup()

# Rename some phases
track_filtered_pks <- track_filtered_pks %>% ungroup() %>% group_by(TrackNameUnique, well_name_p) %>%
  partition(cluster) %>% 
  mutate(PrevCurNext = getPCN(Phase),
         Prev = getPrev(Phase),
         PhaseLength = rep(rle(Phase)$length,rle(Phase)$length)) %>%
  mutate(Phase_corrected = ifelse(PrevCurNext == "G2-G1/S-G1", "G1", 
                                  ifelse(PrevCurNext %in% c("G2-G1-G2",
                                                            "G2-G1/S-G2",
                                                            "G1-G2-G1",
                                                            "G1-G1/S-G1",
                                                            "G1/S-G1-G1/S",
                                                            "G1/S-G2-G1/S") & PhaseLength < 5,
                                         Prev, Phase))) %>%
  collect()

phase_stats <- track_filtered_pks %>% group_by(TrackNameUnique, well_name_p) %>%
  mutate(phaseLength = unlist(sapply(rle(Phase_corrected)$length,function(x){seq(1,x)}, simplify = T)),
         phaseChange = c(!Phase_corrected[-1] == Phase_corrected[-length(Phase_corrected)],NA),
         timeUntilPhaseChange = ifelse(phaseChange,time_h,NA),
         timeUntilPhaseChange = replace(timeUntilPhaseChange,which(!is.na(timeUntilPhaseChange)),
                                        c(na.omit(timeUntilPhaseChange)[1],diff(na.omit(timeUntilPhaseChange))))) %>%
  filter(phaseChange, timeUntilPhaseChange > 5) %>% 
  mutate(phaseID = row_number()) %>% 
  ungroup() 
phase_stats <- phase_stats %>% group_by(cid, well_name_p, well_name, Phase_corrected, phaseID,
                                         treatment,condition,dose_uM,cell_line) %>%
  summarise(phaseLength = mean(timeUntilPhaseChange, na.rm = T)) %>%
  ungroup() 

if (EXP == 223) {
  phase_stats <- phase_stats %>% mutate(E2_conc = ifelse(condition == "100nM", "100 nM",
                                                         ifelse(condition == "10pM", "0.01 nM",
                                                                ifelse(condition == "1nM", "1 nM","0 nM"))))
  ggplot(phase_stats %>% filter(cell_line == "FUCCI")) + #, condition == "100nM")) + #filter(well_name_p %in% c("B04_1"))) + #
    #geom_boxplot(aes(x = treatment, y = phaseLength, fill = Phase_corrected)) +
    geom_boxplot(aes(x = E2_conc, y = phaseLength, fill = Phase_corrected)) +
    ylab("Phase length (h)") + xlab("E2 concentration") +
    scale_fill_manual(values = c("G1" = "red", "G2" = "green","G1/S" = "gold")) + 
    theme_classic()
} else if (EXP %in% c(230, 233)) {
  ggplot(phase_stats %>% filter(cell_line == "FUCCI")) + #, condition == "100nM")) + #filter(well_name_p %in% c("B04_1"))) + #
    geom_boxplot(aes(x = treatment, y = phaseLength, fill = Phase_corrected)) +
    ylab("Phase length (h)") + xlab("Condition") +
    scale_fill_manual(values = c("G1" = "red", "G2" = "green","G1/S" = "gold")) + 
    theme_classic()
  }
#ggsave(paste0(OUTPUT_PATH,READOUT,"_segmentation/",DATE,"_Figures/",EXP,"_phase_length_unfiltered.pdf"), width = 12, height = 3)

# Update PrevCurNext column
track_filtered_pks <- track_filtered_pks %>% ungroup() %>% 
  group_by(TrackNameUnique, well_name_p) %>%
  partition(cluster) %>% 
  mutate(PrevCurNext = getPCN(Phase_corrected)) %>%
  collect()

## Filter the tracks into single tracks
# Make a plot for the filtered unique tracks
POSITION = "B02_1"
numofplots <- dim(distinct(track_filtered_pks %>% filter(well_name_p == POSITION) %>% ungroup %>%
                             select(cid)))[1]
gplot <- ggplot(track_filtered_pks %>% filter(well_name_p == POSITION)) +
  #geom_vline(aes(xintercept = Division_time, color = Division)) +
  geom_line(aes(x = time_h, y = nucleiSize), color = 'grey') +
  geom_line(aes(x = time_h, y = meanHoechst), color = 'blue') +
  geom_line(aes(x = time_h, y = Cdt1), color = 'red') +
  geom_line(aes(x = time_h, y = Geminin), color = 'green') +
  new_scale_color() +
  #geom_point(aes(x = time_h, y = 0.8, color = Phase), pch = "|") +
  geom_point(aes(x = time_h, y = -0.15, color = Phase_corrected), pch = "|") +
  scale_colour_manual(values = c("G1" = "red", "G1/S" = "gold", "G2" = "green")) +
  theme_classic() +
  ggtitle(paste(EXP,", ", POSITION)) +
  #facet_wrap(facets = ~ TrackNameUnique, ncol = 10)
  facet_grid(cid ~ TrackID, scales = "free_y")
ggsave(paste0(OUTPUT_PATH,READOUT,"_segmentation/",DATE,"_Figures/",EXP,"_Complete_Tracks_", POSITION,".pdf"), plot = gplot, 
       width = 10, height = ceiling(numofplots), limitsize = F)

# From the original tracks, delete duplicates but do not split mother and daughter cells
# Collapse duplicate rows into one
tracks_unique <- track_filtered_pks %>% ungroup() %>%
  filter(!duplicated(track_filtered_pks %>% ungroup() %>%
                       select(well_name_p,cid,timeID,Cdt1,Geminin, meanHoechst, nucleiSize)))
# Remove tracks of length less than 6
tracks_unique <- tracks_unique %>%
  group_by(well_name_p, cid, TrackID) %>%
  mutate(track_length = n()) %>%
  filter(track_length > 5) %>% 
  ungroup()

# # Recalculate phase lengths
tracks_unique <- tracks_unique %>% group_by(well_name_p, cid, TrackID) %>%
  mutate(PhaseLength = unlist(sapply(rle(Phase_corrected)$length,function(x){rep(x,x)}, simplify = T)),
         Phase_corrected = ifelse(Phase_corrected == "G2", "S-G2-M", Phase_corrected))

# Make a plot for the filtered unique tracks
POSITION = "B02_1"
numofplots <- dim(distinct(tracks_unique %>% filter(well_name_p == POSITION) %>% ungroup %>%
                             select(cid)))[1]
gplot <- ggplot(tracks_unique %>% filter(well_name_p == POSITION)) +
  #geom_vline(aes(xintercept = Division_time, color = Division)) +
  geom_line(aes(x = time_h, y = nucleiSize), color = 'grey') +
  geom_line(aes(x = time_h, y = meanHoechst), color = 'blue') +
  geom_line(aes(x = time_h, y = Cdt1), color = 'red') +
  geom_line(aes(x = time_h, y = Geminin), color = 'green') +
  new_scale_color() +
  #geom_point(aes(x = time_h, y = 0.8, color = Phase), pch = "|") +
  geom_point(aes(x = time_h, y = -0.15, color = Phase_corrected), pch = "|") +
  scale_colour_manual(values = c("G1" = "red", "G1/S" = "gold", "S-G2-M" = "green")) +
  theme_classic() +
  ggtitle(paste(EXP,", ", POSITION)) +
  #facet_wrap(facets = ~ TrackNameUnique, ncol = 10)
  facet_grid(cid ~ TrackID, scales = "free_y")
ggsave(paste0(OUTPUT_PATH,READOUT,"_segmentation/",DATE,"_Figures/",EXP,"_Unique_Tracks_", POSITION,".pdf"), plot = gplot, 
       width = 10, height = ceiling(numofplots), limitsize = F)

if (EXP == "223"){
  # Make a plot for one unique track
  POSITION = "C02_3" # C02_3 is 10 pM FUCCI
  # Make a plot for one track (in Experiment 223, C02_3, 98)
  cols <- c("Nuclei size"="grey","Hoechst intensity"="blue","Cdt1 intensity"="red","Geminin intensity"="green")
  ggplot(tracks_unique %>% filter(well_name_p == POSITION, cid == 98, TrackID == 1) %>%
           filter(!is.na(Phase_corrected))) +
    geom_line(aes(x = time_h, y = nucleiSize, color = 'Nuclei size'), size = 1.5) +
    geom_line(aes(x = time_h, y = meanHoechst, color = 'Hoechst intensity'), size = 1.5) +
    geom_line(aes(x = time_h, y = Cdt1, color = 'Cdt1 intensity'), size = 1.5) +
    geom_line(aes(x = time_h, y = Geminin, color = 'Geminin intensity'), size = 1.5) +
    scale_colour_manual(values = cols,
                        name = "Readout") +
    new_scale_color() +
    xlab("Time (h)") + ylab("Intensity (a.u.)") +
    geom_point(aes(x = time_h, y = -0.15, color = Phase_corrected), pch = "|", size = 5) +
    scale_colour_manual(values = c("G1" = "red", "G1/S" = "gold", "S-G2-M" = "green"),
                        name = "Phase") +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 12)) #+
  # ggtitle(paste(EXP,", ", POSITION, ", cid = 98, TrackID = 1"))
  ggsave(paste0(OUTPUT_PATH,READOUT,"_segmentation/",DATE,"_Figures/Fig4D_",EXP,"_", POSITION,"_Track98",".pdf"), 
         width = 6, height = 4)

  # Make a heatmap-like plot for DMSO conditions
  POSITION = "B02_3" # B02_3 is DMSO
  # Remove initial NAs
  hmdata <- tracks_unique %>%
    select(well_name_p, cid, TrackID, time_h, TrackNameUnique, Phase_corrected) %>%
    filter(well_name_p == POSITION)
  hmdata <- hmdata %>%
    group_by(well_name_p, cid, TrackID) %>%
    mutate(time_h = time_h - 2.5)
  track_order <- hmdata %>%
    arrange(TrackNameUnique) %>%
    pull(TrackNameUnique)
  hmdata <- hmdata %>%
    mutate(TrackNameUnique = factor(TrackNameUnique, levels = rev(unique(track_order))))
  ggplot(hmdata %>% filter(well_name_p == POSITION, time_h >= 0,!is.na(Phase_corrected))) + #!is.na(Phase_corrected))) +
    geom_tile(aes(x = time_h, y = as.factor(TrackNameUnique), fill = Phase_corrected)) +
    scale_fill_manual(values = c("G1" = "red", "G1/S" = "gold", "S-G2-M" = "green"),
                      name = "Phase") +
    xlab("Time (h)") + ylab("Tracks") +
    xlim(0,44) +
    ggtitle("DMSO") +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 20, hjust = 0.5))
  ggsave(paste0(OUTPUT_PATH,READOUT,"_segmentation/",DATE,"_Figures/FigS4A_",EXP,"_", POSITION,"_DMSO_TrackOverview",".pdf"), 
         width = 4, height = 6)
  
  POSITION = "E02_1" # E02_3 is 100 nM
  # Remove initial NAs
  hmdata <- tracks_unique %>%
    select(well_name_p, cid, TrackID, time_h, TrackNameUnique, Phase_corrected) %>%
    filter(well_name_p == POSITION)
  hmdata <- hmdata %>%
    group_by(well_name_p, cid, TrackID) %>%
    mutate(time_h = time_h - 2.5)
  track_order <- hmdata %>%
    arrange(TrackNameUnique) %>%
    pull(TrackNameUnique)
  hmdata <- hmdata %>%
    mutate(TrackNameUnique = factor(TrackNameUnique, levels = rev(unique(track_order))))
  # Make a heatmap-like plot
  ggplot(hmdata %>% filter(well_name_p == POSITION, time_h >= 0,!is.na(Phase_corrected))) + #!is.na(Phase_corrected))) +
    geom_tile(aes(x = time_h, y = as.factor(TrackNameUnique), fill = Phase_corrected)) +
    scale_fill_manual(values = c("G1" = "red", "G1/S" = "gold", "S-G2-M" = "green"),
                      name = "Phase") +
    xlab("Time (h)") + ylab("Tracks") +
    xlim(0,44) +
    ggtitle("100 nM E2") +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 20, hjust = 0.5))
  ggsave(paste0(OUTPUT_PATH,READOUT,"_segmentation/",DATE,"_Figures/Fig4E_",EXP,"_", POSITION,"_100nM_TrackOverview",".pdf"), 
         width = 4, height = 6)
} else if (EXP == "230"){
  POSITION = "C02_1" # C02_1 is 100 nM in Mock
  # Remove initial NAs
  hmdata <- tracks_unique %>%
    select(well_name_p, cid, TrackID, time_h, TrackNameUnique, Phase_corrected) %>%
    filter(well_name_p == POSITION)
  hmdata <- hmdata %>%
    group_by(well_name_p, cid, TrackID) %>%
    mutate(time_h = time_h - 2.5)
  track_order <- hmdata %>%
    arrange(TrackNameUnique) %>%
    pull(TrackNameUnique)
  hmdata <- hmdata %>%
    mutate(TrackNameUnique = factor(TrackNameUnique, levels = rev(unique(track_order))))
  # Make a heatmap-like plot
  ggplot(hmdata %>% filter(well_name_p == POSITION, time_h >= 0,!is.na(Phase_corrected))) + #!is.na(Phase_corrected))) +
    geom_tile(aes(x = time_h, y = as.factor(TrackNameUnique), fill = Phase_corrected)) +
    scale_fill_manual(values = c("G1" = "red", "G1/S" = "gold", "S-G2-M" = "green"),
                      name = "Phase") +
    xlab("Time (h)") + ylab("Tracks") +
    xlim(0,44) +
    theme_classic() +
    ggtitle("Mock") +
    theme(axis.text.y = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 20, hjust = 0.5))
  ggsave(paste0(OUTPUT_PATH,READOUT,"_segmentation/",DATE,"_Figures/FigS5A_",EXP,"_", POSITION,"_Mock_100nM_TrackOverview",".pdf"), 
         width = 4, height = 6)
  
  POSITION = "C05_1" # C02_1 is 100 nM in siPGR
  # Remove initial NAs
  hmdata <- tracks_unique %>%
    select(well_name_p, cid, TrackID, time_h, TrackNameUnique, Phase_corrected) %>%
    filter(well_name_p == POSITION)
  hmdata <- hmdata %>%
    group_by(well_name_p, cid, TrackID) %>%
    mutate(time_h = time_h - 2.5)
  track_order <- hmdata %>%
    arrange(TrackNameUnique) %>%
    pull(TrackNameUnique)
  hmdata <- hmdata %>%
    mutate(TrackNameUnique = factor(TrackNameUnique, levels = rev(unique(track_order))))
  # Make a heatmap-like plot
  ggplot(hmdata %>% filter(well_name_p == POSITION, time_h >= 0,!is.na(Phase_corrected))) + #!is.na(Phase_corrected))) +
    geom_tile(aes(x = time_h, y = as.factor(TrackNameUnique), fill = Phase_corrected)) +
    scale_fill_manual(values = c("G1" = "red", "G1/S" = "gold", "S-G2-M" = "green"),
                      name = "Phase") +
    xlab("Time (h)") + ylab("Tracks") +
    xlim(0,44) +
    ggtitle("siPGR") +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 20, hjust = 0.5))
  ggsave(paste0(OUTPUT_PATH,READOUT,"_segmentation/",DATE,"_Figures/FigS5B_",EXP,"_", POSITION,"_siPGR_100nM_TrackOverview",".pdf"), 
         width = 4, height = 6)
  
  POSITION = "C04_1" # C02_1 is 100 nM in siGREB1
  # Remove initial NAs
  hmdata <- tracks_unique %>%
    select(well_name_p, cid, TrackID, time_h, TrackNameUnique, Phase_corrected) %>%
    filter(well_name_p == POSITION)
  hmdata <- hmdata %>%
    group_by(well_name_p, cid, TrackID) %>%
    mutate(time_h = time_h - 2.5)
  track_order <- hmdata %>%
    arrange(TrackNameUnique) %>%
    pull(TrackNameUnique)
  hmdata <- hmdata %>%
    mutate(TrackNameUnique = factor(TrackNameUnique, levels = rev(unique(track_order))))
  # Make a heatmap-like plot
  ggplot(hmdata %>% filter(well_name_p == POSITION, time_h >= 0,!is.na(Phase_corrected))) + #!is.na(Phase_corrected))) +
    geom_tile(aes(x = time_h, y = as.factor(TrackNameUnique), fill = Phase_corrected)) +
    scale_fill_manual(values = c("G1" = "red", "G1/S" = "gold", "S-G2-M" = "green"),
                      name = "Phase") +
    xlab("Time (h)") + ylab("Tracks") +
    xlim(0,44) +
    ggtitle("siGREB1") +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 20, hjust = 0.5))
  ggsave(paste0(OUTPUT_PATH,READOUT,"_segmentation/",DATE,"_Figures/FigS5C_",EXP,"_", POSITION,"_siGREB1_100nM_TrackOverview",".pdf"), 
         width = 4, height = 6)
  
  POSITION = "C07_1" # C07_1 is 100 nM in siTFF1
  # Remove initial NAs
  hmdata <- tracks_unique %>%
    select(well_name_p, cid, TrackID, time_h, TrackNameUnique, Phase_corrected) %>%
    filter(well_name_p == POSITION)
  hmdata <- hmdata %>%
    group_by(well_name_p, cid, TrackID) %>%
    mutate(time_h = time_h - 2.5)
  track_order <- hmdata %>%
    arrange(TrackNameUnique) %>%
    pull(TrackNameUnique)
  hmdata <- hmdata %>%
    mutate(TrackNameUnique = factor(TrackNameUnique, levels = rev(unique(track_order))))
  # Make a heatmap-like plot
  ggplot(hmdata %>% filter(well_name_p == POSITION, time_h >= 0,!is.na(Phase_corrected))) + #!is.na(Phase_corrected))) +
    geom_tile(aes(x = time_h, y = as.factor(TrackNameUnique), fill = Phase_corrected)) +
    scale_fill_manual(values = c("G1" = "red", "G1/S" = "gold", "S-G2-M" = "green"),
                      name = "Phase") +
    xlab("Time (h)") + ylab("Tracks") +
    xlim(0,44) +
    ggtitle("siTFF1") +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 20, hjust = 0.5))
  ggsave(paste0(OUTPUT_PATH,READOUT,"_segmentation/",DATE,"_Figures/FigS5D_",EXP,"_", POSITION,"_siTFF1_100nM_TrackOverview",".pdf"), 
         width = 4, height = 6)
  }


save(track_filtered_pks, file = paste0(OUTPUT_PATH,"/",READOUT,"_segmentation/RData/",EXP,"_",READOUT,"_track_filtered_pks.R"))
#load(file = paste0(OUTPUT_PATH,"/",READOUT,"_segmentation/RData/",EXP,"_",READOUT,"_track_filtered_pks.R"))

save(tracks_unique, file = paste0(OUTPUT_PATH,"/",READOUT,"_segmentation/RData/",EXP,"_",READOUT,"_tracks_unique.R"))
#load(file = paste0(OUTPUT_PATH,"/",READOUT,"_segmentation/RData/",EXP,"_",READOUT,"_tracks_unique.R"))


