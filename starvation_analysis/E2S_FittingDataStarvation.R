library("FME")
library(tidyverse)

rm(list = ls())

setwd("/data/muriel/Projects/E2S/Analysis/00_Imaging/ExpXXX_Starvation")

# Get the data
all_data <- read_csv("DataOutput/all_data.csv")
data_sum <- read_csv("DataOutput/data_sum.csv")

# Remove everything but DMSO
all_data <- all_data %>% filter(treatment == "DMSO") %>% select(-c(treatment))
data_sum <- data_sum %>% filter(treatment == "DMSO") %>% select(-c(treatment))

# Get the data
data_GREB1 <- as.data.frame(data_sum %>% filter(Protein == "GREB1") %>%
                              pivot_wider(id_cols = StarvationTime, names_from = Protein, values_from = c(IntensityMean, IntensitySD)) %>%
                              dplyr::rename("time" = "StarvationTime",
                                            "GREB1" = "IntensityMean_GREB1",
                                            "sd" = "IntensitySD_GREB1"))

data_PR <- as.data.frame(data_sum %>% filter(Protein == "PR") %>%
                            pivot_wider(id_cols = StarvationTime, names_from = Protein, values_from = c(IntensityMean, IntensitySD)) %>%
                            dplyr::rename("time" = "StarvationTime",
                                          "PR" = "IntensityMean_PR",
                                          "sd" = "IntensitySD_PR"))

data_TFF1 <- as.data.frame(data_sum %>% filter(Protein == "TFF1") %>%
                             pivot_wider(id_cols = StarvationTime, names_from = Protein, values_from = c(IntensityMean, IntensitySD)) %>%
                             dplyr::rename("time" = "StarvationTime",
                                           "TFF1" = "IntensityMean_TFF1",
                                           "sd" = "IntensitySD_TFF1"))

data_long <- as.data.frame(all_data %>% ungroup %>% rename("name" = "Protein", 
                                                           "time" = "StarvationTime",
                                                           "val" = "Intensity") %>% 
                             select(c(name, time, val)) %>%
                             relocate("name"))

data_long_mean <- data_long %>% group_by(name, time) %>%
  summarise(meanval = mean(val))

# Initialize model
DegradModel <- function (pars, G_0, P_0, T_0) {
  
  derivs <- function(time, y, pars) {
    with (as.list(c(pars, y)), {
      dGREB1 <- - dg * GREB1
      dPR <- - dp * PR
      dTFF1 <- - dt * TFF1
      
      return(list(c(dGREB1, dPR, dTFF1)))
    })
  }
  
  # initial conditions
  y <- c(GREB1 = G_0, PR = P_0, TFF1 = T_0)
  
  times <- c(seq(0, 120, 1))
  out <- ode(y = y, parms = pars, times = times, func = derivs)
  
  as.data.frame(out)
}

# Define parameters
pars <- c(dg = 0.01, dp = 0.02, dt = 0.02)

# Simulate
out <- DegradModel(pars = pars, 
                   G_0 = data_long_mean %>% filter(name == "GREB1", time == 0) %>% pull(meanval), 
                   P_0 = data_long_mean %>% filter(name == "PR", time == 0) %>% pull(meanval), 
                   T_0 = data_long_mean %>% filter(name == "TFF1", time == 0) %>% pull(meanval))

# Plot the output
out
plot(out$time, out$GREB1, main = "GREB1", ylab = "Intensity (a.u.)",
     xlab = "Time (hr)", type = "b")

# Make cost function
# ModelCost <- function (pars) {
#   out <- DegradModel(pars, G_0 = 9, P_0 = 3, T_0 = 24)
#   cost0 <- modCost(model = out, obs = data_GREB1, err = "sd")
#   cost1 <- modCost(model = out, obs = data_PR, err = "sd", cost = cost0)
#   return(modCost(model = out, obs = data_TFF1, err = "sd", cost = cost1))
# }
ModelCost <- function (pars) {
  out <- DegradModel(pars, 
                     G_0 = data_long_mean %>% filter(name == "GREB1", time == 0) %>% pull(meanval), 
                     P_0 = data_long_mean %>% filter(name == "PR", time == 0) %>% pull(meanval), 
                     T_0 = data_long_mean %>% filter(name == "TFF1", time == 0) %>% pull(meanval))
  cost0 <- modCost(model = out, obs = data_long, y = "val")
  cost1 <- modCost(model = out, obs = data_long, y = "val", cost = cost0)
  return(modCost(model = out, obs = data_long, y = "val", cost = cost1))
}

# Get cost
# ModelCost(pars)$model
ModelCost(pars)$model

plot(ModelCost(pars), xlab="time")

# Do fitting
ModelCost2 <- function(lpars) {
  ModelCost(c(exp(lpars)))
}
Pars <- pars * 2
Fit <- modFit(f = ModelCost2, p = log(Pars))
exp(coef(Fit))
deviance(Fit)

# Make plot
ini <- DegradModel(pars = c(Pars), 
                   G_0 = data_long_mean %>% filter(name == "GREB1", time == 0) %>% pull(meanval), 
                   P_0 = data_long_mean %>% filter(name == "PR", time == 0) %>% pull(meanval), 
                   T_0 = data_long_mean %>% filter(name == "TFF1", time == 0) %>% pull(meanval))
final <- DegradModel(pars = c(exp(coef(Fit))), 
                     G_0 = data_long_mean %>% filter(name == "GREB1", time == 0) %>% pull(meanval), 
                     P_0 = data_long_mean %>% filter(name == "PR", time == 0) %>% pull(meanval), 
                     T_0 = data_long_mean %>% filter(name == "TFF1", time == 0) %>% pull(meanval))

pdf(file="Figures/FigS3C_FittedStarvationData.pdf", width = 7.3, height = 2.5)
# GREB1
par(mfrow = c(1,3))
plot(x = data_GREB1$time, y = data_GREB1$GREB1, xlab = "Time (h)", ylab = "GFP intensity (a.u.)",
     ylim = c(0,10),main = "GREB1", pch = 19)
#lines(ini$time, ini$GREB1, lty = 2)
lines(final$time, final$GREB1)
legend("topright", c("data", "fitted"),
       lty = c(NA,1), pch = c(19, NA))
# PR
plot(x = data_PR$time, y = data_PR$PR, xlab = "Time (h)", ylab = "GFP intensity (a.u.)",
     ylim = c(0,4),main = "PR", pch = 19)
#lines(ini$time, ini$PR, lty = 2)
lines(final$time, final$PR)
# TFF1
plot(x = data_TFF1$time, y = data_TFF1$TFF1, xlab = "Time (h)", ylab = "GFP intensity (a.u.)",
     ylim = c(0,25),main = "TFF1", pch = 19)
#lines(ini$time, ini$TFF1, lty = 2)
lines(final$time, final$TFF1)
dev.off()

par(mfrow = c(1, 1))

write_csv(final, path = "DegradationSimulation.csv")
