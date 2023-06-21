setwd("/data/muriel/Projects/E2S/Analysis/00_Imaging/ExpXXX_Starvation")

rm(list = ls())

# Load libraries
library(tidyverse)
# install.packages("FME")
# library("FME")
# vignette("FME")
library(reshape2) # library for reshaping data (tall-narrow <-> short-wide)
library(deSolve) # library for solving differential equations
library(minpack.lm) # library for least squares fit using levenberg-marquart algorithm

# Get data files

files <- grep('data.txt', dir("Data_Britt_Starvation"), value=TRUE)

# Load datasets in a list
data_aslist <- list()
for (i in seq_along(files)) {
  # Get the data set
  data_tmp <- read_delim(paste0("Data_Britt_Starvation/",files[i]), delim = "\t")
  # Remove "imageNumber" as column name
  if ("imageNumber" %in% colnames(data_tmp)) {
    data_tmp <- data_tmp %>% select(-c(imageNumber))}
  data_aslist[[i]] <- data_tmp
}

# data_1 <- data_aslist[[1]]
# data_2 <- data_aslist[[2]]
# data_3 <- data_aslist[[3]]
# data_4 <- data_aslist[[4]]
# data_5 <- data_aslist[[5]]
# data_6 <- data_aslist[[6]]
# data_7 <- data_aslist[[7]]
# data_8 <- data_aslist[[8]]
# data_9 <- data_aslist[[9]]

# Bind all data frames into one
all_data <- bind_rows(data_aslist)

# Correct the timeID
all_data <- all_data %>% 
  mutate(timeID_orig = timeID,
         timeID = as.numeric(sapply(plateID, function(id) {str_split(id, "MT_t")[[1]][2]})),
         Protein = sapply(cell_line, function(cl) {str_split(cl, "-")[[1]][2]}),
         well_row = sapply(locationID, function(id) {substr(id,1,1)}))


# Filter only first 4 rows (option 1-4)
all_data <- all_data %>% filter(well_row %in% c("C","D","E","F"))

# Create a column with the options and with the experiment name
all_data <- all_data %>% 
  mutate(option = ifelse(medium == "Comp_medium_NS", 1,
                         ifelse(medium == "Comp_medium_S", 2,
                                ifelse(medium == "Starv_medium_1", 3,
                                       4))),
         experiment = sapply(plateID, function(id) {str_split(id, "_")[[1]][2]})) %>%
  relocate(option)

# Create a column with the time after starvation and the intensity
all_data <- all_data %>% 
  mutate(StarvationTime = ifelse(option == 1, timeID,
                                 ifelse(option == 2, timeID + 24,
                                        ifelse(option == 3, timeID + 48, 
                                               timeID + 72))))

# Pivot to wide format
all_data <- pivot_wider(all_data, id_cols = c("option", "treatment","medium", "dose_uM", "Protein", "StarvationTime",
                                              "plateID", "timeID_orig", "timeID", "cell_line", "locationID","control",
                                              "experiment"),
                        names_from = variable, values_from = value)

# Get one column with intensities
all_data <- all_data %>% 
  mutate(Intensity = ifelse(Protein %in% c("GREB1","TFF1"), 
                            Cytoplasm_N_Intensity_IntegratedIntensity_Image_GFP,
                            Nuclei_Intensity_IntegratedIntensity_Image_GFP))

# Change PGR in PR
all_data <- all_data %>% 
  mutate(Protein = ifelse(Protein == "PGR", "PR", Protein))

# Calculate average of technical replicates
all_data_mean <- all_data %>% group_by(option, treatment, medium, dose_uM, Protein, StarvationTime,
                                       experiment, cell_line, timeID) %>%
  summarise(MeanIntensity = mean(Intensity),
            SDIntensity = sd(Intensity))

# Reproduce figs Britt
ggplot(data = all_data_mean %>% filter(treatment == "DMSO")) +
  geom_point(aes(x = timeID, y = MeanIntensity, color = as.factor(option))) +
  geom_line(aes(x = timeID, y = MeanIntensity, color = as.factor(option)), alpha = 0.25) +
  geom_errorbar(aes(x = timeID, ymin = MeanIntensity - SDIntensity, 
                    ymax =  MeanIntensity + SDIntensity, color = as.factor(option)), width = 2) +
  scale_color_viridis_d(name = "Option") +
  xlab("Time (h)") + ylab("GFP intensity (a.u.)") +
  theme_bw() +
  facet_grid(~Protein, scales="free")

# Get relative changes in intensity
all_data_mean <- all_data_mean %>% group_by(Protein, option) %>% 
  mutate(RelativeIntensity = MeanIntensity / MeanIntensity[1])

# Make summary
data_sum <- all_data_mean %>% group_by(StarvationTime, Protein, treatment) %>%
  summarise(IntensityMean = mean(MeanIntensity),
            IntensitySD = sd(MeanIntensity),
            RelativeIntensityMean = mean(RelativeIntensity),
            RelativeIntensitySD = sd(RelativeIntensity))

# Make plots
# Intensities
ggplot(data = data_sum %>% filter(treatment == "DMSO")) +
  geom_point(aes(x = StarvationTime, y = IntensityMean, color = Protein)) +
  geom_errorbar(aes(x = StarvationTime, ymin = IntensityMean - IntensitySD, 
                    ymax = IntensityMean + IntensitySD, color = Protein), width = 5) +
  xlab("Time in starvation (h)") + 
  ylab("GFP intensity (a.u.)") +
  scale_color_viridis_d() +
  scale_x_continuous(breaks = c(0,24,48,72,96,120)) +
  theme_classic()
ggsave("Figures/FigS3B_MeanExpressionInStarvation.pdf", width = 4, height = 3)

# Relative intensity changes
ggplot(data = data_sum %>% filter(RelativeIntensityMean != 1)) +
  geom_point(aes(x = StarvationTime, y = RelativeIntensityMean, color = Protein)) +
  geom_errorbar(aes(x = StarvationTime, ymin = RelativeIntensityMean - RelativeIntensitySD,
                    ymax = RelativeIntensityMean + RelativeIntensitySD, color = Protein), width = 5) +
  xlab("Time in starvation (h)") + 
  ylab("Relative GFP intensity (a.u.)") +
  xlim(0,120) +
  scale_color_viridis_d() +
  theme_classic()
ggsave("Figures/RelativeExpressionInStarvation.pdf", width = 4, height = 3)

write_csv(data_sum, path = "DataOutput/data_sum.csv")
write_csv(all_data_mean, path = "DataOutput/all_data_mean.csv")
write_csv(all_data, path = "DataOutput/all_data.csv")
