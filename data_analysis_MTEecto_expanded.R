library(dplyr)
library(ggplot2)

pairs_data = read.csv("clean_data_MTEecto_expanded.csv", stringsAsFactors = FALSE)
pairs_data = pairs_data[pairs_data$Class == "Insecta",]
#pairs_data = pairs_data[pairs_data$studyID == "FAL02A",]

pairs_data$x_axis = (1/(pairs_data$initial_temp + 273.15)) - (1/(pairs_data$final_temp + 273.15))
pairs_data$y_axis = log(pairs_data$initial_metrate) - log(pairs_data$final_metrate)

pairs_data$line_slope = pairs_data$Ea / .00008617
pairs_data$line_intercept = 0

#T1 = initial_temp; T2 = final_temp
#pairs_data$x_axis = (pairs_data$final_temp - pairs_data$initial_temp) / (pairs_data$final_temp * pairs_data$initial_temp)
pairs_data$interval = pairs_data$final_temp - pairs_data$initial_temp
#interval_type: small is 1-15, large is 16-29
pairs_data = pairs_data %>% mutate(interval_type = ifelse(interval < 16, "small", "large"))
pairs_data$absolute = pairs_data$initial_temp
#absolute_type: small is 2-13, large is 14-25
pairs_data = pairs_data %>% mutate(absolute_type = ifelse(absolute < 14, "small", "large"))
pairs_data$interval_absolute = paste(pairs_data$interval_type, pairs_data$absolute_type, sep = ", ")

#pairs_data = pairs_data %>% filter(interval_type == "large")
#pairs_data = pairs_data %>% filter(absolute_type == "large")

#TODO: color code according to mass categories
ggplot(pairs_data, aes(x = x_axis, y = y_axis)) +
  geom_point(aes(colour = interval_absolute)) +
  geom_abline(data = pairs_data, aes(slope = line_slope, intercept = line_intercept))


all_studyIDs_plots = 
ggplot() +
  geom_point(data = pairs_data, aes(x = log(initial_mass), y = log(constantmass_metrate)), color = "orange", size = 2) +
  geom_vline(data = pairs_data, aes(xintercept = log(initial_mass)), color = "orange", size = 0.2) +
  geom_point(data = pairs_data, aes(x = log(constantmetrate_mass), y = log(initial_metrate)), color = "red", size = 2) +
  stat_smooth(data = pairs_data, aes(x = log(constantmetrate_mass), y = log(initial_metrate)), color = "red", method = "lm", se = FALSE, size = 0.2) +
  geom_point(data = pairs_data, aes(x = log(final_mass), y = log(final_metrate)), color = "blue", size = 2) +
  stat_smooth(data = pairs_data, aes(x = log(final_mass), y = log(final_metrate)), color = "blue", method = "lm", se = FALSE, size = 0.2) +
  geom_point(data = pairs_data, aes(x = log(initial_mass), y = log(initial_metrate)), color = "black", size = 2) +
  xlab("Log mass") +
  ylab("Log metabolic rate") +
  facet_wrap(~studyID, scales = "free")

ggsave("all_pairs_studyID.jpg", plot = all_studyIDs_plots, width = 20, height = 20)

#added legend manually

#TODO: connecting same temp points and remove colored lines
ggplot() +
  geom_point(data = pairs_data, aes(x = log(initial_mass), y = log(constantmass_metrate)), color = "orange", size = 2) +
  #geom_point(data = pairs_data, aes(x = log(constantmetrate_mass), y = log(initial_metrate)), color = "red", size = 2) +
  geom_point(data = pairs_data, aes(x = log(final_mass), y = log(final_metrate)), color = "blue", size = 2) +
  geom_line(data = pairs_data, aes(x = log(final_mass), y = log(final_metrate)))
  #geom_point(data = pairs_data, aes(x = log(initial_mass), y = log(initial_metrate)), color = "black", size = 2) +
  xlab("Log mass") +
  ylab("Log metabolic rate") +
  facet_wrap(~studyID, scales = "free")
