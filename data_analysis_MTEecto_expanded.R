library(dplyr)
library(ggplot2)

pairs_data = read.csv("clean_data_MTEecto_expanded.csv", stringsAsFactors = FALSE)

all_studyIDs_plots = ggplot() +
  geom_point(data = pairs_data, aes(x = log(initial_mass), y = log(constantmass_metrate)), color = "orange", size = 2) +
  geom_vline(data = pairs_data, aes(xintercept = log(initial_mass)), color = "orange", size = 1) +
  geom_point(data = pairs_data, aes(x = log(constantmetrate_mass), y = log(initial_metrate)), color = "red", size = 2) +
  stat_smooth(data = pairs_data, aes(x = log(constantmetrate_mass), y = log(initial_metrate)), color = "red", method = "lm", se = FALSE, fullrange = TRUE) +
  geom_point(data = pairs_data, aes(x = log(final_mass), y = log(final_metrate)), color = "blue", size = 2) +
  stat_smooth(data = pairs_data, aes(x = log(final_mass), y = log(final_metrate)), color = "blue", method = "lm", se = FALSE, fullrange = TRUE) +
  geom_point(data = pairs_data, aes(x = log(initial_mass), y = log(initial_metrate)), color = "black", size = 2) +
  xlab("Mass") +
  ylab("Metabolic rate") +
  facet_wrap(~studyID, scales = "free")

ggsave("all_pairs_studyID.jpg", plot = all_studyIDs_plots, width = 20, height = 20)

#added legend manually
