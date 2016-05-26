# Read in clean data
replicates_data = read.csv("clean_data_MTEecto.csv", stringsAsFactors = FALSE)

#---------------METABOLIC RATE ANALYSIS & VISUALIZATION----------

### Metabolic rate percent differences
replicates_data$PD_MR_CS = ((replicates_data$constantmass_metrate / replicates_data$initial_metrate) - 1) * 100
replicates_data$PD_MR_VS = ((replicates_data$final_metrate / replicates_data$initial_metrate) - 1) * 100

library(dplyr)
species_data = replicates_data %>% 
  group_by(species) %>% 
  summarise_each(funs(mean, sd), PD_MR_CS, PD_MR_VS)

occurrences = replicates_data %>% 
  group_by(species) %>% 
  summarise(count = n())
species_data$occurrences = occurrences$count

# Looking at variation in percent differences for each species w/ sufficient pairs
species_data = species_data[order(species_data$occurrences, decreasing = TRUE),]
sub_species_data = species_data[1:7,]
bp_VS = barplot(sub_species_data$PD_MR_VS_mean)
arrows(bp_VS, sub_species_data$PD_MR_VS_mean + sub_species_data$PD_MR_VS_sd, bp_VS, sub_species_data$PD_MR_VS_mean, angle = 90, code = 0)
arrows(bp_VS, sub_species_data$PD_MR_VS_mean - sub_species_data$PD_MR_VS_sd, bp_VS, sub_species_data$PD_MR_VS_mean, angle = 90, code = 0)
bp_CS = barplot(sub_species_data$PD_MR_CS_mean)
arrows(bp_CS, sub_species_data$PD_MR_CS_mean + sub_species_data$PD_MR_CS_sd, bp_CS, sub_species_data$PD_MR_CS_mean, angle = 90, code = 0)
arrows(bp_CS, sub_species_data$PD_MR_CS_mean - sub_species_data$PD_MR_CS_sd, bp_CS, sub_species_data$PD_MR_CS_mean, angle = 90, code = 0)

# PD distributions by species
for (species_name in unique(sub_species_data$species)){
  sub = filter(replicates_data, species == species_name)
  #plot(density(sub$PD_MR_CS), main = c(species_name, "CS"))
  hist(sub$PD_MR_CS, main = c(species_name, "CS"))
  abline(v = mean(sub$PD_MR_CS))
  #plot(density(sub$PD_MR_VS), main = c(species_name, "VS"))
  hist(sub$PD_MR_VS, main = c(species_name, "VS"))
  abline(v = mean(sub$PD_MR_VS))
}

# Looking at one species
single = replicates_data[replicates_data$species == "Rana sylvatica",]
barplot(single$PD_MR_CS, col = "white")
barplot(single$PD_MR_VS, col = rgb(0, 0, 0, alpha = 0.3), add = TRUE)
barplot(single$constantmass_metrate, col = rgb(0, 0, 1, alpha = 0.3))
barplot(single$final_metrate, col = rgb(0, 1, 0, alpha = 0.3), add = TRUE)
barplot(single$initial_metrate, col = rgb(1, 0, 0, alpha = 0.3), add = TRUE)

species_data$PD_MR_diff = species_data$PD_MR_CS_mean - species_data$PD_MR_VS_mean
species_data = species_data[order(species_data$PD_MR_diff, decreasing = TRUE),]

# Barplot
barplot(species_data$PD_MR_VS_mean, col = "green", xlab = "Species", ylab = "MR change due to temp increase")
barplot(species_data$PD_MR_CS_mean, col = rgb(0, 0, 0, alpha = 0.3), add = TRUE)
legend("topright", c("constant mass/no STR", "varying mass/STR"), fill = c("white", "green"))

# T-test or something? 
#TODO: paired t-test comparing constant mass and varying mass MR PDs

#------------------MASS ANALYSIS & VISUALIZATION---------------

### Mass percent difference
replicates_data$PD_mass = ((replicates_data$constantmetrate_mass / replicates_data$initial_mass) - 1) * 100

mass_comparison = replicates_data %>%
  group_by(species) %>%
  summarise_each(funs(mean), PD_mass)
species_data$PD_mass = mass_comparison$PD_mass

# Plot
plot(density(species_data$PD_mass), main = "Needed decrease in mass", xlim = c(-40, 10))
species_data$PD_mass_null = rnorm(45, mean = 0, sd = sd(species_data$PD_mass))
lines(density(species_data$PD_mass_null), col = "red")
legend("topright", c("red = null", "black = actual"))
