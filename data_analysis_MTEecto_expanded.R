library(dplyr)
library(ggplot2)
library(lme4)

pairs_data = read.csv("clean_data_MTEecto_expanded.csv", stringsAsFactors = FALSE)

pairs_data$mass_constant_slope = -pairs_data$Ea / .00008617
pairs_data$mass_constant_intercept = 0
pairs_data$initial_temp_K = pairs_data$initial_temp + 273.15
pairs_data$final_temp_K = pairs_data$final_temp + 273.15
pairs_data$temp_axis_num = -(pairs_data$initial_temp_K - pairs_data$final_temp_K)
pairs_data$temp_axis_den = pairs_data$initial_temp_K * pairs_data$final_temp_K
pairs_data$temp_axis = pairs_data$temp_axis_num / pairs_data$temp_axis_den
pairs_data$metab_axis = log(pairs_data$final_metrate) - log(pairs_data$initial_metrate)

### All classes together
ggplot(pairs_data, aes(x = temp_axis, y = metab_axis)) +
  geom_point() +
  geom_abline(data = pairs_data, aes(slope = mass_constant_slope, intercept = mass_constant_intercept)) +
  labs(x = "Temperature axis \n -(T1-T2)/(T1*T2)", y = "Metabolic rate difference \n log(R2) - log(R1)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Variance explained by no mass change line
r2 = 1 - var(with(pairs_data, metab_axis - (temp_axis * mass_constant_slope))) / var(with(pairs_data, metab_axis))

### With mass percent change lines
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

### With data trend line
ggplot(pairs_data, aes(x = temp_axis, y = metab_axis)) +
  geom_point() +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_constant_intercept)) +
  geom_smooth(method = "lm", se = FALSE, size = .6) +
  facet_wrap(~Class) +
  labs(x = "Temperature axis \n (T1-T2)/(T1*T2)", y = "Metabolic rate difference \n log(R2)-log(R1)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Temp difference mostly affects x-axis
ggplot(pairs_data, aes(x = temp_axis_num, y = temp_axis)) +
  geom_point() +
  facet_wrap(~Class) +
  labs(x = "Temperature difference \n T1 - T2", y = "Temperature axis \n (T1-T2)/(T1*T2)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### By mass
#pairs_data$residual_alt = pairs_data$metab_axis - (pairs_data$mass_constant_slope * pairs_data$temp_axis)
pairs_data$residual = pairs_data$exponent * (log(pairs_data$final_mass) - log(pairs_data$initial_mass))

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

### Residuals model for species variation
residual_model = lmer(residual ~ log(initial_mass) + initial_temp + temp_axis + (1|species) + (1|Class) + (1|studyID) + (1|study), data = pairs_data)

### Compensation mass plot
pairs_data$needed_mass_change = (pairs_data$constantmetrate_mass - pairs_data$initial_mass) / abs(pairs_data$initial_mass) * 100
pairs_data$actual_mass_change = (pairs_data$final_mass - pairs_data$initial_mass) / abs(pairs_data$initial_mass) * 100

ggplot(pairs_data, aes(x = needed_mass_change, y = actual_mass_change)) +
  geom_point() +
  coord_cartesian(xlim = c(-100, 100), ylim = c(-100, 100)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_vline(xintercept = 0, color = "grey") +
  geom_abline(intercept = 0, slope = 1)

