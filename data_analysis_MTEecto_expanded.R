library(dplyr)
library(ggplot2)
library(lme4)

pairs_data = read.csv("clean_data_MTEecto_expanded.csv", stringsAsFactors = FALSE)

### All plot calculations
pairs_data$needed_mass_change = (pairs_data$constantmetrate_mass - pairs_data$initial_mass) / abs(pairs_data$initial_mass) * 100
pairs_data$actual_mass_change = (pairs_data$final_mass - pairs_data$initial_mass) / abs(pairs_data$initial_mass) * 100
pairs_data$mass_constant_slope = -pairs_data$Ea / .00008617
pairs_data$mass_constant_intercept = 0
pairs_data$initial_temp_K = pairs_data$initial_temp + 273.15
pairs_data$final_temp_K = pairs_data$final_temp + 273.15
pairs_data$temp_axis_num = -(pairs_data$initial_temp_K - pairs_data$final_temp_K)
pairs_data$temp_axis_den = pairs_data$initial_temp_K * pairs_data$final_temp_K
pairs_data$temp_axis = pairs_data$temp_axis_num / pairs_data$temp_axis_den
pairs_data$metab_axis = log(pairs_data$final_metrate) - log(pairs_data$initial_metrate)
pairs_data$expect_metrate_change = pairs_data$temp_axis * pairs_data$mass_constant_slope
pairs_data$observ_metrate_change = pairs_data$metab_axis

### MAIN: Compensation mass plot
ggplot(pairs_data, aes(x = needed_mass_change, y = actual_mass_change)) +
  geom_point() +
  coord_cartesian(xlim = c(-100, 100), ylim = c(-100, 100)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_vline(xintercept = 0, color = "grey") +
  geom_abline(intercept = 0, slope = 1)

### STATS: Mass comparison t-test
t.test(log(pairs_data$final_mass), log(pairs_data$constantmetrate_mass), paired = TRUE)

### MAIN: Variance plot with trend line
pairs_data = pairs_data %>% 
  mutate(point_color = ifelse(Class == "Amphibia", "Amphibians", "Non-Amphibians"))

ggplot(pairs_data, aes(x = expect_metrate_change, y = observ_metrate_change)) +
  geom_point(aes(color = point_color)) +
  scale_color_manual(values = c("red", "black",  "black", "blue")) +
  geom_abline(aes(color = "Model", intercept = 0, slope = 1)) +
  geom_smooth(aes(color = "Trend"), method = "lm", se = FALSE, size = .6) +
  guides(color = guide_legend(override.aes = list(linetype = c("blank", "blank", "solid", "solid"), 
                                                  shape = c(16, 16, NA, NA)))) +
  labs(x = "Expected metabolic rate change", y = "Observed metabolic rate change", color = "") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### SUPPLEMENT: Temp difference mostly affects x-axis
ggplot(pairs_data, aes(x = temp_axis_num, y = temp_axis)) +
  geom_point() +
  facet_wrap(~Class) +
  labs(x = "Temperature difference \n T1 - T2", y = "Temperature axis \n (T1-T2)/(T1*T2)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### STATS: Variance explained by no mass change line
unexplained_variance = var(with(pairs_data, observ_metrate_change - expect_metrate_change))
total_variance = var(with(pairs_data, observ_metrate_change))
r2 = 1 - unexplained_variance / total_variance

### MAIN: With mass percent change lines
pairs_data$mass_intercept_1 = pairs_data$exponent * log(2)
pairs_data$mass_intercept_2 = pairs_data$exponent * log(1.5)
pairs_data$mass_intercept_3 = pairs_data$exponent * log(0.75)
pairs_data$mass_intercept_4 = pairs_data$exponent * log(0.5)

ggplot(pairs_data, aes(x = temp_axis, y = metab_axis)) +
  geom_point() +
  geom_abline(data = pairs_data, aes(slope = mass_constant_slope, intercept = mass_constant_intercept)) +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_intercept_1), color = "grey") +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_intercept_2), color = "grey") +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_intercept_3), color = "grey") +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_intercept_4), color = "grey") +
  facet_wrap(~Class) +
  labs(x = "Temperature axis \n -(T1-T2)/(T1*T2)", y = "Metabolic rate difference \n log(R2) - log(R1)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### SUPPLEMENT: With data trend line
ggplot(pairs_data, aes(x = temp_axis, y = metab_axis)) +
  geom_point() +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_constant_intercept)) +
  geom_smooth(method = "lm", se = FALSE, size = .6) +
  facet_wrap(~Class) +
  labs(x = "Temperature axis \n (T1-T2)/(T1*T2)", y = "Metabolic rate difference \n log(R2)-log(R1)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Full model
#pairs_data$residual_alt = pairs_data$metab_axis - (pairs_data$mass_constant_slope * pairs_data$temp_axis)
pairs_data$residual = pairs_data$exponent * (log(pairs_data$final_mass) - log(pairs_data$initial_mass))

pairs_data$temp_axis_rescale = pairs_data$temp_axis * 100000
full_model = lmer(residual ~ log(initial_mass) + initial_temp + temp_axis_rescale + (1|species) + 
                    (1|Class) + (1|studyID) + (1|study), data = pairs_data)
summary(full_model)

#Likelihood ratio tests
null_mass = lmer(residual ~ initial_temp + temp_axis_rescale + (1|species) + 
                   (1|Class) + (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
model_mass = lmer(residual ~ log(initial_mass) + initial_temp + temp_axis_rescale + (1|species) + 
                    (1|Class) + (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
anova(null_mass, model_mass)

null_temp_init = lmer(residual ~ log(initial_mass) + temp_axis_rescale + (1|species) + 
                        (1|Class) + (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
model_temp_init = lmer(residual ~ log(initial_mass) + initial_temp + temp_axis_rescale + (1|species) + 
                         (1|Class) + (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
anova(null_temp_init, model_temp_init)

### STATS: Linear mixed model
final_model = lmer(residual ~ temp_axis_rescale + (1|species) + (1|Class) + 
                     (1|studyID) + (1|study), data = pairs_data)
summary(final_model)

### SUPPLEMENT: Linear mixed model diagnostics
coef(final_model)
plot(fitted(final_model), residuals(final_model)) #linearity good, heteroscedasticity
abline(h = 0)
hist(residuals(final_model)) # normality looks good in histogram
qqnorm(residuals(final_model))
qqline(residuals(final_model)) #normality of residuals heavy-tailed in qqplot

#Likelihood ratio tests
null_temp_diff = lmer(residual ~ (1|species) + (1|Class) + 
                        (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
model_temp_diff = lmer(residual ~ temp_axis_rescale + (1|species) + (1|Class) + 
                         (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
anova(null_temp_diff, model_temp_diff)

null_sp = lmer(residual ~ temp_axis_rescale + (1|Class) + 
                 (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
model_sp = lmer(residual ~ temp_axis_rescale + (1|species) + (1|Class) + 
                  (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
anova(null_sp, model_sp)

null_class = lmer(residual ~ temp_axis_rescale + (1|species) +  
                    (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
model_class = lmer(residual ~ temp_axis_rescale + (1|species) + (1|Class) + 
                     (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
anova(null_class, model_class)

null_studyID = lmer(residual ~ temp_axis_rescale + (1|species) + (1|Class) + 
                      (1|study), data = pairs_data, REML = FALSE)
model_studyID = lmer(residual ~ temp_axis_rescale + (1|species) + (1|Class) + 
                       (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
anova(null_studyID, model_studyID)

null_study = lmer(residual ~ temp_axis_rescale + (1|species) + (1|Class) + 
                    (1|studyID), data = pairs_data, REML = FALSE)
model_study = lmer(residual ~ temp_axis_rescale + (1|species) + (1|Class) + 
                     (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
anova(null_study, model_study)

### SUPPLEMENT: No systematic variation due to absolute mass or temp
ggplot(pairs_data, aes(x = log(initial_mass), y = residual)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  facet_wrap(~Class, scales = "free_x") +
  ylim(-1.4, 1.4) +
  labs(x = "log(mass)", y = "Mass residual") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(pairs_data, aes(x = initial_temp_K, y = residual)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  facet_wrap(~Class, scales = "free_x") +
  ylim(-1.4, 1.4) +
  labs(x = "Initial temperature (K)", y = "Mass residual") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
