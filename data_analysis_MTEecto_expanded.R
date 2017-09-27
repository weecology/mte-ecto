library(dplyr)
library(ggplot2)

pairs_data = read.csv("clean_data_MTEecto_expanded.csv", stringsAsFactors = FALSE)
pairs_data = pairs_data[pairs_data$Class == "Insecta",]
#pairs_data = pairs_data[pairs_data$studyID == "FAL02A",]

pairs_data$x_axis = (1/(pairs_data$final_temp + 273.15)) - (1/(pairs_data$initial_temp + 273.15))
pairs_data$y_axis = log(pairs_data$final_metrate) - log(pairs_data$initial_metrate)

#pairs_data$x_axis_alt = ((pairs_data$final_temp + 273.15) - (pairs_data$initial_temp + 273.15)) / ((pairs_data$final_temp + 273.15) * (pairs_data$initial_temp + 273.15))
#pairs_data$x_axis_num = ((pairs_data$final_temp + 273.15) - (pairs_data$initial_temp + 273.15))
#pairs_data$x_axis_den = ((pairs_data$final_temp + 273.15) * (pairs_data$initial_temp + 273.15))

pairs_data$line_slope = -pairs_data$Ea / .00008617
pairs_data$line_intercept = 0

# #T1 = initial_temp; T2 = final_temp
# #pairs_data$x_axis = (pairs_data$final_temp - pairs_data$initial_temp) / (pairs_data$final_temp * pairs_data$initial_temp)
# pairs_data$interval = pairs_data$final_temp - pairs_data$initial_temp
# #interval_type: small is 1-15, large is 16-29
# pairs_data = pairs_data %>% mutate(interval_type = ifelse(interval < 16, "small", "large"))
# pairs_data$absolute = pairs_data$initial_temp
# #absolute_type: small is 2-13, large is 14-25
# pairs_data = pairs_data %>% mutate(absolute_type = ifelse(absolute < 14, "small", "large"))
# pairs_data$interval_absolute = paste(pairs_data$interval_type, pairs_data$absolute_type, sep = ", ")

pairs_data = pairs_data %>% 
  arrange(initial_mass) %>%
  mutate(row_num = 1:nrow(pairs_data)) %>%
  mutate(size_category = case_when(
    row_num < (nrow(pairs_data) / 3) ~ "small", 
    row_num >= (nrow(pairs_data) / 3) & row_num <= ((nrow(pairs_data) / 3) * 2) ~ "medium", 
    row_num > ((nrow(pairs_data) / 3) * 2) ~ "large"))

#TODO: color code according to mass categories
ggplot(pairs_data, aes(x = -x_axis, y = y_axis, color = size_category)) +
  geom_point() +
  geom_abline(data = pairs_data, aes(slope = line_slope, intercept = line_intercept)) 
#+    facet_wrap(~Class)

ggplot(pairs_data, aes(x = -x_axis, y = y_axis)) +
  geom_point() +
  geom_abline(aes(slope = line_slope, intercept = line_intercept)) +
  #geom_smooth(method = "lm", se = FALSE) +
  theme_bw()

#TODO: mass decrease lines
pairs_data = pairs_data %>% 
  mutate(percent_change_mass = ((final_mass - initial_mass) / initial_mass) * 100)
percent_value = -15
initial_over_final = (percent_value / 100) + 1
pairs_data$percent_decrease_5_intercept = pairs_data$exponent * log(initial_over_final)
ggplot(pairs_data, aes(x = -x_axis, y = y_axis)) +
  geom_point() +
  geom_abline(data = pairs_data, aes(slope = line_slope, intercept = line_intercept)) +
  geom_abline(aes(slope = line_slope, intercept = percent_decrease_5_intercept), color = "green") +
  #geom_smooth(method = "lm", se = FALSE) +
  theme_bw()

# example_sp = pairs_data %>% 
#   filter(species == "Drosophila melanogaster  ") %>% 
#   filter(study == "TM61" | study == "N40") %>% 
#   #filter(trial == "A" | trial == "B") %>%
#   mutate(temp_diff = final_temp - initial_temp) %>%
#   filter(temp_diff == 5 | temp_diff == 15)
# example_sp$temp_diff = as.factor(example_sp$temp_diff)
# 
# ggplot(example_sp, aes(x = -x_axis, y = y_axis, color = temp_diff)) +
#   geom_point()
# 
# example_sp_sum = example_sp %>% 
#   group_by(temp_diff) %>% #plus group by species
#   summarize(x_axis_mean = mean(x_axis), 
#             x_axis_sd = sd(x_axis), 
#             y_axis_mean = mean(y_axis), 
#             y_axis_sd = sd(y_axis))
# 
# ggplot(example_sp_sum, aes(x = -x_axis_mean, y = y_axis_mean)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = y_axis_mean - y_axis_sd, ymax = y_axis_mean + y_axis_sd), width = .000001) +
#   geom_errorbarh(aes(xmin = -(x_axis_mean - x_axis_sd), xmax = -(x_axis_mean + x_axis_sd)), height = .01)

pairs_data_test = pairs_data %>%
  mutate(temp_diff = final_temp - initial_temp) %>%
  group_by(species, temp_diff) %>%
  summarize(count = n(), 
            x_axis_mean = mean(x_axis), 
            x_axis_sd = sd(x_axis), 
            y_axis_mean = mean(y_axis), 
            y_axis_sd = sd(y_axis))

ggplot(pairs_data_test, aes(x = -x_axis_mean, y = y_axis_mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = y_axis_mean - y_axis_sd, ymax = y_axis_mean + y_axis_sd), width = .000001, color = "red") +
    geom_errorbarh(aes(xmin = -(x_axis_mean - x_axis_sd), xmax = -(x_axis_mean + x_axis_sd)), height = .05, color = "orange")
  


