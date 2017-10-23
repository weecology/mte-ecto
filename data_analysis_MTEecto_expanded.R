library(dplyr)
library(ggplot2)

pairs_data = read.csv("clean_data_MTEecto_expanded.csv", stringsAsFactors = FALSE)

pairs_data$x_axis = (1/(pairs_data$final_temp + 273.15)) - (1/(pairs_data$initial_temp + 273.15))
pairs_data$y_axis = log(pairs_data$final_metrate) - log(pairs_data$initial_metrate)
pairs_data$line_slope = -pairs_data$Ea / .00008617
pairs_data$line_intercept = 0

### With mass percent change lines
pairs_data$mass_intercept_1 = pairs_data$exponent * log(2)
pairs_data$mass_intercept_2 = pairs_data$exponent * log(1.5)
pairs_data$mass_intercept_3 = pairs_data$exponent * log(0.75)
pairs_data$mass_intercept_4 = pairs_data$exponent * log(0.5)

ggplot(pairs_data, aes(x = -x_axis, y = y_axis)) +
  geom_point() +
  geom_abline(data = pairs_data, aes(slope = line_slope, intercept = line_intercept)) +
  geom_abline(aes(slope = line_slope, intercept = mass_intercept_1), color = "grey") +
  geom_abline(aes(slope = line_slope, intercept = mass_intercept_2), color = "grey") +
  geom_abline(aes(slope = line_slope, intercept = mass_intercept_3), color = "grey") +
  geom_abline(aes(slope = line_slope, intercept = mass_intercept_4), color = "grey") +
  facet_wrap(~Class) +
  labs(x = "Temperature axis \n (T2-T1)/(T1*T2)", y = "Metabolic rate difference \n log(R2)-log(R1)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### With data trend line
ggplot(pairs_data, aes(x = -x_axis, y = y_axis)) +
  geom_point() +
  geom_abline(aes(slope = line_slope, intercept = line_intercept)) +
  geom_smooth(method = "lm", se = FALSE, size = .6) +
  facet_wrap(~Class) +
  labs(x = "Temperature axis \n (T2-T1)/(T1*T2)", y = "Metabolic rate difference \n log(R2)-log(R1)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### By temperature difference and absolute temperature
#TODO: cowplot these together? 
pairs_data$temp_diff = pairs_data$final_temp - pairs_data$initial_temp
pairs_data$temp_prod = pairs_data$final_temp * pairs_data$initial_temp

ggplot(pairs_data, aes(x = temp_diff, y = -x_axis)) +
  geom_point() +
  facet_wrap(~Class) +
  labs(x = "Temperature difference \n T2 - T1", y = "Temperature axis \n (T2-T1)/(T1*T2)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(pairs_data, aes(x = temp_prod, y = -x_axis)) +
  geom_point() +
  facet_wrap(~Class) +
  labs(x = "Temperature product \n T2 * T1", y = "Temperature axis \n (T2-T1)/(T1*T2)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### By mass
pairs_data$line_value = pairs_data$line_slope * -pairs_data$x_axis + pairs_data$line_intercept
pairs_data$mass_residual = pairs_data$y_axis - pairs_data$line_value

ggplot(pairs_data, aes(x = log(initial_mass), y = mass_residual)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  facet_wrap(~Class, scales = "free_x") +
  ylim(-1.4, 1.4) +
  labs(x = "log(mass)", y = "Mass residual") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

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
