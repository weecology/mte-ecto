library(dplyr)
library(ggplot2)

pairs_data = read.csv("clean_data_MTEecto_expanded.csv", stringsAsFactors = FALSE)
pairs_data$mass_constant_slope = -pairs_data$Ea / .00008617
pairs_data$mass_constant_intercept = 0
pairs_data$initial_temp_K = pairs_data$initial_temp + 273.15
pairs_data$final_temp_K = pairs_data$final_temp + 273.15
pairs_data$temp_axis_num = -(pairs_data$initial_temp_K - pairs_data$final_temp_K)
pairs_data$temp_axis_den = pairs_data$initial_temp_K * pairs_data$final_temp_K
pairs_data$temp_axis = pairs_data$temp_axis_num / pairs_data$temp_axis_den
pairs_data$metab_axis = log(pairs_data$final_metrate) - log(pairs_data$initial_metrate)

fishes_data = filter(pairs_data, Class == "Actinoperygii")

fishes_data$initial_exponent = 0.0417 * (1/(.0000862*(fishes_data$initial_temp + 273.15))) - 0.8923
fishes_data$final_exponent = 0.0417 * (1/(.0000862*(fishes_data$final_temp + 273.15))) - 0.8923

fishes_data$initial_metrate_adjusted = (fishes_data$initial_mass ^ fishes_data$initial_exponent) * (exp(((fishes_data$Ea) / (.00008617 * (fishes_data$initial_temp + 273.15)))))
fishes_data$final_metrate_adjusted = (fishes_data$initial_mass ^ fishes_data$final_exponent) * (exp(((fishes_data$Ea) / (.00008617 * (fishes_data$final_temp + 273.15)))))
fishes_data$metab_axis_adjusted = fishes_data$final_metrate_adjusted - fishes_data$initial_metrate_adjusted

fishes_data$mass_intercept_1 = fishes_data$exponent * log(2)
fishes_data$mass_intercept_2 = fishes_data$exponent * log(1.5)
fishes_data$mass_intercept_3 = fishes_data$exponent * log(0.75)
fishes_data$mass_intercept_4 = fishes_data$exponent * log(0.5)

ggplot(fishes_data, aes(x = temp_axis, y = metab_axis)) +
  geom_point() +
  geom_point(aes(x = temp_axis, y = metab_axis_adjusted), color = "red", alpha = 0.7) +
  geom_abline(intercept = 0, lty = 3) +
  coord_cartesian(xlim = c(0, 0.0003), ylim = c(-1, 2)) +
  geom_abline(data = fishes_data, aes(slope = mass_constant_slope, intercept = mass_constant_intercept)) +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_intercept_1), color = "grey") +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_intercept_2), color = "grey") +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_intercept_3), color = "grey") +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_intercept_4), color = "grey") +
  labs(x = "Temperature axis \n -(T1-T2)/(T1*T2)", y = "Metabolic rate difference \n log(R2) - log(R1)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
