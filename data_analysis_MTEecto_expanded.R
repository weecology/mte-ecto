library(dplyr)
library(ggplot2)

pairs_data = read.csv("clean_data_MTEecto_expanded.csv", stringsAsFactors = FALSE)

pairs_data$x_axis = (1/(pairs_data$final_temp + 273.15)) - (1/(pairs_data$initial_temp + 273.15))
pairs_data$y_axis = log(pairs_data$final_metrate) - log(pairs_data$initial_metrate)
pairs_data$line_slope = -pairs_data$Ea / .00008617
pairs_data$line_intercept = 0

### Base plot
ggplot(pairs_data, aes(x = -x_axis, y = y_axis)) +
  geom_point() +
  geom_abline(data = pairs_data, aes(slope = line_slope, intercept = line_intercept)) +
  facet_wrap(~Class) +
  labs(x = "Temperature difference \n (T2-T1)/(T1*T2)", y = "Metabolic rate difference \n log(R2)-log(R1)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### With data trend line
ggplot(pairs_data, aes(x = -x_axis, y = y_axis)) +
  geom_point() +
  geom_abline(aes(slope = line_slope, intercept = line_intercept)) +
  geom_smooth(method = "lm", se = FALSE, size = .6) +
  facet_wrap(~Class) +
  labs(x = "Temperature difference \n (T2-T1)/(T1*T2)", y = "Metabolic rate difference \n log(R2)-log(R1)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### By temperature difference and absolute temperature
pairs_data$interval = pairs_data$final_temp - pairs_data$initial_temp
pairs_data = pairs_data %>% mutate(interval_type = ifelse(interval < 16, "small", "large"))
pairs_data$absolute = pairs_data$initial_temp
pairs_data = pairs_data %>% mutate(absolute_type = ifelse(absolute < 14, "small", "large"))
pairs_data$interval_absolute = paste(pairs_data$interval_type, pairs_data$absolute_type, sep = ", ")

ggplot(pairs_data, aes(x = -x_axis, y = y_axis, color = interval_absolute)) +
  geom_point() +
  geom_abline(data = pairs_data, aes(slope = line_slope, intercept = line_intercept)) +
  facet_wrap(~Class) +
  labs(x = "Temperature difference \n (T2-T1)/(T1*T2)", y = "Metabolic rate difference \n log(R2)-log(R1)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### By mass
pairs_data = pairs_data %>% 
  arrange(initial_mass) %>%
  mutate(row_num = 1:nrow(pairs_data)) %>%
  mutate(size_category = case_when(
    row_num < (nrow(pairs_data) / 3) ~ "small", 
    row_num >= (nrow(pairs_data) / 3) & row_num <= ((nrow(pairs_data) / 3) * 2) ~ "medium", 
    row_num > ((nrow(pairs_data) / 3) * 2) ~ "large"))

ggplot(pairs_data, aes(x = -x_axis, y = y_axis, color = size_category)) +
  geom_point() +
  geom_abline(data = pairs_data, aes(slope = line_slope, intercept = line_intercept)) +
  facet_wrap(~Class) +
  labs(x = "Temperature difference \n (T2-T1)/(T1*T2)", y = "Metabolic rate difference \n log(R2)-log(R1)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

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
  labs(x = "Temperature difference \n (T2-T1)/(T1*T2)", y = "Metabolic rate difference \n log(R2)-log(R1)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Averages by species and temperature difference
pairs_data_by_sp_temp = pairs_data %>%
  mutate(temp_diff = final_temp - initial_temp) %>%
  group_by(Class, species, temp_diff, line_slope, line_intercept) %>%
  summarize(count = n(), 
            x_axis_mean = mean(x_axis), 
            y_axis_mean = mean(y_axis), 
            x_axis_sd = sd(x_axis), 
            y_axis_sd = sd(y_axis))

ggplot(pairs_data_by_sp_temp, aes(x = -x_axis_mean, y = y_axis_mean)) +
  geom_point() +
  geom_abline(data = pairs_data_by_sp_temp, aes(slope = line_slope, intercept = line_intercept)) +
  facet_wrap(~Class) +
  labs(x = "Temperature difference \n (T2-T1)/(T1*T2)", y = "Metabolic rate difference \n log(R2)-log(R1)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### With error bars
ggplot(pairs_data_by_sp_temp, aes(x = -x_axis_mean, y = y_axis_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = y_axis_mean - y_axis_sd, ymax = y_axis_mean + y_axis_sd), width = .00001, color = "red") +
  #geom_errorbarh(aes(xmin = -(x_axis_mean - x_axis_sd), xmax = -(x_axis_mean + x_axis_sd)),  height = .5, color = "orange") +
  geom_abline(data = pairs_data_by_sp_temp, aes(slope = line_slope, intercept = line_intercept)) +
  facet_wrap(~Class) +
  labs(x = "Temperature difference \n (T2-T1)/(T1*T2)", y = "Metabolic rate difference \n log(R2)-log(R1)") +
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

pairs_data = pairs_data %>% 
  mutate(equal_interval = case_when(
    interval < 10 ~ "small", 
    interval >= 10 & interval < 20 ~ "medium", 
    interval >= 20 ~ "large"
  ))

ggplot(pairs_data, aes(x = needed_mass_change, y = actual_mass_change, color = equal_interval)) +
  geom_point() +
  coord_cartesian(xlim = c(-100, 100), ylim = c(-100, 100)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_vline(xintercept = 0, color = "grey") +
  geom_abline(intercept = 0, slope = 1)

pairs_data = pairs_data %>%
  arrange(interval) %>% 
  mutate(row_num = 1:nrow(pairs_data)) %>% 
  mutate(equal_group = case_when(
    row_num < (nrow(pairs_data) / 3) ~ "small", 
    row_num >= (nrow(pairs_data) / 3) & row_num < ((nrow(pairs_data) / 3) * 2) ~ "medium", 
    row_num >= (nrow(pairs_data) / 3 * 2) ~ "large"
  ))

ggplot(pairs_data, aes(x = needed_mass_change, y = actual_mass_change, color = equal_group)) +
  geom_point() +
  coord_cartesian(xlim = c(-100, 100), ylim = c(-100, 100)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_vline(xintercept = 0, color = "grey") +
  geom_abline(intercept = 0, slope = 1)

ggplot(pairs_data, aes(interval, fill = equal_group)) +
  geom_histogram()
