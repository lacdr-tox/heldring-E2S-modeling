# Make figures for paper
# Muriel Heldring
# 2022-08-30

rm(list = ls())

# Load packages
library(tidyverse)
library(scales)

# GFP data and model simulations
setwd("/data/muriel/Projects/E2S/Analysis/00_Imaging/Exp013_GFP_Britt")

# Make figure directory
dir.create("Figures/v2/")

# Load the GFP data
data <- read_csv("/data/muriel/Projects/E2S/Modeling/Input/version2/Exp_E2data.csv")

# Rename dose_uM column
data <- rename(data,dose_nM = dose_uM)
# Make dose_uM column a factor and rename PR
data <- data %>% mutate(dose_nM = as.factor(dose_nM),
                        StateVar = ifelse(StateVar == "PGR", "PR",StateVar))

# Calculate mean and sd
data_summary <- data %>% filter(!dose_nM %in% c(1000, 10000)) %>% #filter(dose_nM %in% c(0.001,0.01,0.1,1,10,100)) %>%
  group_by(StateVar, timeID, dose_nM) %>%
  #mean count and sd
  summarise(Real_GFP_intensity_mean = mean(Real_GFP_intensity, na.rm = T),
            Real_GFP_intensity_sd = sd(Real_GFP_intensity, na.rm = T))

# Make figure 1c
ggplot(data = data_summary, aes(x = timeID, y = Real_GFP_intensity_mean, group = dose_nM)) +
  geom_ribbon(aes(ymin = Real_GFP_intensity_mean-Real_GFP_intensity_sd,
                  ymax = Real_GFP_intensity_mean+Real_GFP_intensity_sd, 
                  fill = dose_nM),
              alpha = 0.2, show.legend = FALSE) +
  geom_point(aes(color = dose_nM), pch = 20) +
  geom_line(aes(color = dose_nM)) +
  ylab("GFP intensity (a.u.)") + xlab("Time (h)") +
  scale_color_viridis_d(name = "Dose (nM)") +
  theme_classic() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.position = "top") +
  scale_fill_viridis_d(name = "Dose (nM)") +
  facet_grid(StateVar~., scales = "free_y")
ggsave("Figures/v2/Fig1C_GFPdata.pdf", width = 4, height = 6)

# Load the simulation data
simu_data_MN03 <- read_csv("/data/muriel/Projects/E2S/Modeling/Output/SimulationData/MN03_simu_data.csv") %>% select(-c(X1))
simu_data_MN08 <- read_csv("/data/muriel/Projects/E2S/Modeling/Output/SimulationData/MN08_simu_data.csv") %>% select(-c(X1))
simu_data_MN10 <- read_csv("/data/muriel/Projects/E2S/Modeling/Output/SimulationData/MN10_simu_data.csv") %>% select(-c(X1))

all_simu_data <- left_join(left_join(simu_data_MN03, simu_data_MN08, 
                                     by = c("dose_nM", "timeID", "StateVar"), suffix = c("_MN03","_MN08")),
                           simu_data_MN10 %>% rename("GFP_MN10" = "GFP"), by = c("dose_nM", "timeID", "StateVar")) %>%
  pivot_longer(cols = c("GFP_MN03","GFP_MN08","GFP_MN10"), names_to = "Model", values_to = "GFP") %>%
  mutate(Model = ifelse(Model == "GFP_MN03","I",
                        ifelse(Model == "GFP_MN08","II","III")))


data_list <- list("MN03" = simu_data_MN03, "MN08" = simu_data_MN08, "MN10" = simu_data_MN10)

for (model in names(data_list)) {
  
  simu_data <- data_list[[model]]
  
  # Make dose_uM column a factor and rename PR
  simu_data <- simu_data %>% mutate(dose_nM = as.factor(dose_nM),
                                    StateVar = ifelse(StateVar == "PGR", "PR",StateVar))
  
  # Make Figs
  facetlabels <- c("0.001" = "0.001 nM",
                   "0.01" = "0.01 nM",
                   "0.1" = "0.1 nM",
                   "1" = "1 nM",
                   "10" = "10 nM",
                   "100" = "100 nM")
  ggplot() +
    geom_point(data = data %>% filter(dose_nM %in% c(0.001,0.01,0.1,1,10,100)), 
               aes(x = timeID, y = Real_GFP_intensity,
                   color = dose_nM), 
               size = 0.2, alpha = 0.2) +
    geom_line(data = simu_data, aes(x = timeID, y = GFP,
                                    color = dose_nM), 
              size = 1) +
    ylab("GFP intensity (a.u.)") + xlab("Time (h)") +
    scale_color_viridis_d(name = "Dose (nM)") +
    theme_classic() +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=16),
          legend.position = 'none') +
    facet_grid(StateVar ~ dose_nM, scales = "free_y", 
               labeller = labeller(.cols = facetlabels))
  if (model == "MN03") {
    ggsave(paste0("Figures/v2/FigSX_GFPdata_and_simulation_",model,".pdf"), width = 7, height = 4)}
  else if (model == "MN08") {
    ggsave(paste0("Figures/v2/FigSX_GFPdata_and_simulation_",model,".pdf"), width = 7, height = 4)}
  else if (model == "MN10") {
    ggsave(paste0("Figures/v2/FigSX_GFPdata_and_simulation_",model,".pdf"), width = 7, height = 4)}
  
  # # Make figure SX
  # ggplot() +
  #   geom_errorbar(data = data_summary %>% filter(dose_nM %in% c(0.001,0.01,0.1,1,10,100)), 
  #                 aes(x = timeID, ymin = Real_GFP_intensity_mean-Real_GFP_intensity_sd, 
  #                     ymax = Real_GFP_intensity_mean+Real_GFP_intensity_sd), 
  #                 color = viridis_pal()(11)[6]) +
  #   #geom_point(aes(x = timeID, y = Real_GFP_intensity_mean), color = "black", pch = 19) +
  #   geom_point(data = data_summary %>% filter(dose_nM %in% c(0.001,0.01,0.1,1,10,100)), 
  #              aes(x = timeID, y = Real_GFP_intensity_mean), 
  #              color = viridis_pal()(11)[5], pch = 20) +
  #   geom_line(data = simu_data, aes(x = timeID, y = GFP), color = "black", size = 1) +
  #   ylab("GFP intensity (a.u.)") + xlab("Time (h)") +
  #   theme_classic() +
  #   theme(axis.text = element_text(size=12),
  #         axis.title = element_text(size=16)) +
  #   facet_grid(StateVar ~ dose_nM)
  # ggsave("Figures/FigSX_GFPdata_and_simulation.pdf", width = 7, height = 5)
}

ggplot() +
  geom_point(data = data %>% filter(dose_nM %in% c(0.001,0.01,0.1,1,10,100)), 
             aes(x = timeID, y = Real_GFP_intensity),
                 color = "black", 
             size = 0.2, alpha = 0.2) +
  geom_line(data = all_simu_data %>% 
              mutate(Model = factor(Model, levels = c("I","II","III"))), 
            aes(x = timeID, y = GFP, color = Model),
            lty = "solid",
            size = 1, alpha = 1) +
  # geom_line(data = all_simu_data %>% filter(Model == "I"), 
  #           aes(x = timeID, y = GFP),
  #           color = "darkred", 
  #           lty = "solid",
  #           size = 1, alpha = 0.5) +
  # geom_line(data = all_simu_data %>% filter(!Model == "III"), 
  #           aes(x = timeID, y = GFP),
  #           color = "black", 
  #           lty = "dashed",
  #           size = 1, alpha = 0.5) +
  ylab("GFP intensity (a.u.)") + xlab("Time (h)") +
  #scale_color_viridis_d(name = "Model", direction = -1, option = "turbo") +
  scale_colour_manual(values = c("#e41a1c", "#4daf4a","#984ea3")) +
  theme_classic() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.position = 'top',
        legend.text = element_text(size=12),
        legend.title = element_text(size=16)) +
  facet_grid(StateVar ~ dose_nM, scales = "free_y", 
             labeller = labeller(.cols = facetlabels))
ggsave(paste0("Figures/v2/FigS1C_GFPdata_and_simulation.pdf"), width = 7, height = 5)

