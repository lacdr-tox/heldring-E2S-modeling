# Statistical test on simulation data of phase lengths

rm(list = ls())
gc()
getwd()

# Load the data
simdat_withoutKD <- read_csv("/data/muriel/Projects/E2S/Modeling/Output/SimulationData/MN10_Phase_duration_simu_withoutKD.csv")

ggplot(simdat_withoutKD,
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

# Check for normality ## VIOLATED!
d_without <- simdat_withoutKD %>% mutate(phase_condition = paste0(Phase_corrected,"_",condition))
# Build the linear model
model  <- lm(PhaseLength ~ phase_condition, data = d_without)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))

# Check normality assumption by group
d_without %>%
  group_by(phase_condition) %>%
  shapiro_test(PhaseLength)

# Non-parametric test
# Welch One way ANOVA test
res.aov2 <- d_without %>% welch_anova_test(PhaseLength ~ phase_condition)
# Pairwise comparisons (Games-Howell)
pwc2 <- d_without %>% games_howell_test(PhaseLength ~ phase_condition)
# Visualization: box plots with p-values
pwc2 <- pwc2 %>% add_xy_position(x = "phase_condition", step.increase = 1)
# Filter out the relevant comparisons (within phases)
pwc2 <- pwc2 %>% mutate(group1_cond = str_split(group1,"_", simplify = T)[,1],
                        group2_cond = str_split(group2,"_", simplify = T)[,1]) %>%
  filter(group1_cond == group2_cond)
ggboxplot(d_without, x = "phase_condition", y = "PhaseLength") +
  stat_pvalue_manual(pwc2, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov2, detailed = TRUE),
    caption = get_pwc_label(pwc2)
  )

# Load the data
simdat_withKD <- read_csv("/data/muriel/Projects/E2S/Modeling/Output/SimulationData/MN10_Phase_duration_simu_withKD.csv")

ggplot(simdat_withKD,
       aes(x = treatment, y = PhaseLength, fill = Phase_corrected)) +
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

# Check for normality ## VIOLATED!
d_with <- simdat_withKD %>% mutate(phase_treatment = paste0(Phase_corrected,"_",treatment))
# Build the linear model
model  <- lm(PhaseLength ~ phase_treatment, data = d_with)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))

# Check normality assumption by group
d_with %>%
  group_by(phase_treatment) %>%
  shapiro_test(PhaseLength)

# Non-parametric test
# Welch One way ANOVA test
res.aov2 <- d_with %>% welch_anova_test(PhaseLength ~ phase_treatment)
# Pairwise comparisons (Games-Howell)
pwc2 <- d_with %>% games_howell_test(PhaseLength ~ phase_treatment)
# Visualization: box plots with p-values
pwc2 <- pwc2 %>% add_xy_position(x = "phase_treatment", step.increase = 1)
# Filter out the relevant comparisons (within phases)
pwc2 <- pwc2 %>% mutate(group1_cond = str_split(group1,"_", simplify = T)[,1],
                        group2_cond = str_split(group2,"_", simplify = T)[,1]) %>%
  filter(group1_cond == group2_cond)
ggboxplot(d_with, x = "phase_treatment", y = "PhaseLength") +
  stat_pvalue_manual(pwc2, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov2, detailed = TRUE),
    caption = get_pwc_label(pwc2)
  )

