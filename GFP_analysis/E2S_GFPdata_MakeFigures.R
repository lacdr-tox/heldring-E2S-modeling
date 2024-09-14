# Make figures for paper
# Muriel Heldring
# 2022-08-30

rm(list = ls())

# Load packages
install.packages("pacman")
library(pacman)
# remove.packages("tidyverse")
# remove.packages("glue")
# install.packages("tidyverse")
p_load(glue)
p_load(tidyverse)
p_load(dplyr)
p_load(stringr)
# library(cli)
# library(tidyverse)
# library(glue)
# library(stringr)
p_load(scales)

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
simu_data_MN03 <- read_csv("/data/muriel/Projects/E2S/Modeling/Output/SimulationData/MN03_simu_data.csv") %>% select(-c(...1))
simu_data_MN08 <- read_csv("/data/muriel/Projects/E2S/Modeling/Output/SimulationData/MN08_simu_data.csv") %>% select(-c(...1))
simu_data_MN10 <- read_csv("/data/muriel/Projects/E2S/Modeling/Output/SimulationData/MN10_simu_data.csv") %>% select(-c(...1))

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

# Compare calibration results for Model III
# Load the simulation data
all_files <- dir("/data/muriel/Projects/E2S/Modeling/Output/SimulationData/")
parm_sets <- all_files[grepl("simu_data_parmset",all_files)]

# read in set 1
all_simu_data <- read_csv(paste0("/data/muriel/Projects/E2S/Modeling/Output/SimulationData/",parm_sets[1])) %>% 
  select(-c("X1")) 


for (i in 2:length(parm_sets)) {
  # Get parmset_nr
  setnr <- ifelse(substr(substr(parm_sets[i],23,24),2,2) == "_",substr(parm_sets[i],23,23),substr(parm_sets[i],23,24))

  # Get cost
  cost <- ifelse(substr(substr(parm_sets[i],34,35),2,2) == ".",substr(parm_sets[i],29,34),substr(parm_sets[i],30,35))
  
  # read in set 'i'
  all_simu_data_tmp <- read_csv(paste0("/data/muriel/Projects/E2S/Modeling/Output/SimulationData/",parm_sets[i])) %>% select(-c("X1"))
  
  # Left join sets
  all_simu_data <- left_join(all_simu_data, all_simu_data_tmp, 
                             by = c("dose_nM", "timeID", "StateVar"), 
                             suffix = c("",paste0("_set",setnr,"_cost",cost)))
}

print(colnames(all_simu_data))

# Add _set1_cost294.44 to GFP
setnr_1 <- ifelse(substr(substr(parm_sets[1],23,24),2,2) == "_",substr(parm_sets[1],23,23),substr(parm_sets[1],23,24))
cost_1 <- ifelse(substr(substr(parm_sets[1],34,35),2,2) == ".",substr(parm_sets[1],29,34),substr(parm_sets[1],30,35))
if (setnr_1 == 1){
  all_simu_data <- all_simu_data %>% rename("GFP_set1_cost294.44" = "GFP")} else {
    warning("setnr is not 1!")
  }

colnames(all_simu_data)
sapply(colnames(all_simu_data)[4:78],function(x){str_split(x,"_")[[1]][3]}, USE.NAMES = F)

all_simu_data_pivot <-  all_simu_data %>%
  pivot_longer(cols = -c("dose_nM","timeID","StateVar"), names_to = "Parmset_Cost", values_to = "GFP") %>%
  mutate(Parmset = sapply(Parmset_Cost,function(x){str_split(x,"_")[[1]][2]}, USE.NAMES = F)) %>%
  mutate(Cost = sapply(Parmset_Cost,function(x){str_split(x,"_")[[1]][3]}, USE.NAMES = F)) %>%
  mutate(Optima = ifelse(Cost == "cost287.67","287.67 (n = 8)",
                         ifelse(Cost == "cost294.44", "294.44 (n = 61)",
                                ifelse(Cost == "cost297.17", "297.17 (n = 1)",
                                       ifelse(Cost == "cost307.30", "307.30 (n = 2)",
                                              ifelse(Cost == "cost326.03","326.03 (n = 2)", "369.00 (n=1)"))))))
# all_simu_data <-  all_simu_data %>%
#   pivot_longer(cols = c("GFP_set1","GFP_set2","GFP_set3","GFP_set4","GFP_set5",
#                         "GFP_set6","GFP_set7","GFP_set8","GFP_set9","GFP_set10",
#                         "GFP_set11","GFP_set12","GFP_set13","GFP_set14","GFP_set15"), names_to = "Parmset", values_to = "GFP") %>%
#   mutate(Parmset = substr(Parmset,5,9)) %>%
#   mutate(Optima = ifelse(Parmset == "set6","326.03 (n = 1)",ifelse(Parmset == "set9", "287.67 (n = 1)", "294.44 (n = 13)")))

# Make Figs
facetlabels <- c("0.001" = "0.001 nM",
                 "0.01" = "0.01 nM",
                 "0.1" = "0.1 nM",
                 "1" = "1 nM",
                 "10" = "10 nM",
                 "100" = "100 nM")
ggplot() +
  geom_point(data = data %>% filter(dose_nM %in% c(1,10,100)) %>%
               group_by(dose_nM, StateVar, timeID) %>%
               summarise(MeanGFP = mean(Real_GFP_intensity),
                         sdGFP = sd(Real_GFP_intensity)), 
             aes(x = timeID, y = MeanGFP),
             color = "grey", 
             size = 0.2, alpha = 1) +
  geom_errorbar(data = data %>% filter(dose_nM %in% c(1,10,100)) %>%
                  group_by(dose_nM, StateVar, timeID) %>%
                  summarise(MeanGFP = mean(Real_GFP_intensity),
                            sdGFP = sd(Real_GFP_intensity)),
                aes(x = timeID, ymin=MeanGFP-sdGFP, ymax=MeanGFP+sdGFP),
                color = "grey", 
                size = 0.2, alpha = 1) +
  geom_line(data = all_simu_data_pivot %>% filter(dose_nM %in% c(1,10,100),
                                                  Cost %in% c("cost287.67","cost294.44","cost297.17")) %>%
              mutate(Optima = factor(Optima)),
            aes(x = timeID, y = GFP, color = Optima),
            size = 0.5, alpha = 1) +
  ylab("GFP intensity (a.u.)") + xlab("Time (h)") +
  #scale_color_viridis_d(name = "Optima", direction = -1, option = "turbo") +
  scale_colour_manual(values = c("#e41a1c", "#4daf4a","#984ea3")) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted"))+
  theme_classic() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.position = 'right',
        legend.text = element_text(size=12),
        legend.title = element_text(size=16)) +
  facet_grid(StateVar ~ dose_nM, scales = "free_y", 
             labeller = labeller(.cols = facetlabels))
ggsave(paste0("Figures/v2/FigS2B_GFPdata_and_simulation_new.pdf"), width = 5, height = 5)
