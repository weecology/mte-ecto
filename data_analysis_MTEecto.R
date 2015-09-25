# Read in clean data
replicates_data = read.csv("clean_data_MTEecto.csv")

# Look at distributions of all mass and metabolic rate values
par(mfrow = c(3, 2))
hist(replicates_data$initial_mass)
hist(replicates_data$initial_metrate)
hist(replicates_data$final_mass)
hist(replicates_data$final_metrate)
hist(replicates_data$constantmass_metrate)
hist(replicates_data$constantmetrate_mass)

# Do normality test on all mass and metabolic rate values

pvalues_list = list(non_trans = numeric())

IM_norm = shapiro.test(replicates_data$initial_mass)
pvalues_list[1,1] = IM_norm$p.value

IME_norm = shapiro.test(replicates_data$initial_metrate)
pvalues_list[[2]] = IME_norm$p.value

shapiro.test(replicates_data$final_mass)
shapiro.test(replicates_data$final_metrate)
shapiro.test(replicates_data$constantmass_metrate)
shapiro.test(replicates_data$constantmetrate_mass)


#-----------------------TRANSFORM MASS & MR VALUES-----------------------------

# Log transform all metabolic rate & mass values
replicates_data$log_initial_mass = log(replicates_data$initial_mass)
replicates_data$log_initial_metrate = log(replicates_data$initial_metrate)
replicates_data$log_final_mass = log(replicates_data$final_temp)
replicates_data$log_final_metrate = log(replicates_data$final_metrate)
replicates_data$log_constantmass_metrate = log(replicates_data$constantmass_metrate)
replicates_data$log_constantmetrate_mass = log(replicates_data$constantmetrate_mass)

# Histograms of transformed values
par(mfrow = c(3, 2))
hist(replicates_data$log_initial_mass)
hist(replicates_data$log_initial_metrate)
hist(replicates_data$log_final_mass)
hist(replicates_data$log_final_metrate)
hist(replicates_data$log_constantmass_metrate)
hist(replicates_data$log_constantmetrate_mass)

# Comparison of final and constant mass transformed metabolic rates
replicates_data$metrate_comparison = replicates_data$log_constantmass_metrate - replicates_data$log_final_metrate
mean(replicates_data$metrate_comparison)
hist(replicates_data$metrate_comparison)
t.test(replicates_data$log_final_metrate, replicates_data$log_constantmass_metrate, paired = TRUE)
boxplot(replicates_data$log_final_metrate, replicates_data$log_constantmass_metrate, main = "Mass change mediates MR increase", xlab = "Mass", ylab = "Log metabolic rate", names = c("actual (i.e., varying)", "constant"))

# Comparison of final and constant MR transformed masses
replicates_data$mass_comparison = replicates_data$log_constantmetrate_mass - replicates_data$log_final_mass
mean(replicates_data$mass_comparison)
hist(replicates_data$mass_comparison)
t.test(replicates_data$log_final_mass, replicates_data$log_constantmetrate_mass, paired = TRUE)
boxplot(replicates_data$log_final_mass, replicates_data$log_constantmetrate_mass, main = "Required mass change to maintain MR", xlab = "Metabolic rate", ylab = "Log mass", names = c("actual", "constant (i.e., compensation mass)"))

#-------------------PERCENT DIFFERENCES OF MASS & MR VALUES--------------------

# Percent differences for metabolic rates
# 1: initial to final (empirical), expect increase
# 2: initial to constant mass (theoretical), expect increase
# 3: final to constant mass, expect increase
replicates_data$PD_1 = ((replicates_data$final_metrate / replicates_data$initial_metrate) - 1) * 100
replicates_data$PD_2 = ((replicates_data$constantmass_metrate / replicates_data$initial_metrate) - 1) * 100
replicates_data$PD_3 = ((replicates_data$constantmass_metrate / replicates_data$final_metrate) - 1) * 100

# Percent differences for masses
# 4: initial to final (empirical), expect decrease
# 5: initial to constant metabolic rate (theoretical), expect decrease
# 6: final to constant metabolic rate, expect decrease
replicates_data$PD_4 = ((replicates_data$final_mass / replicates_data$initial_mass) - 1) * 100
replicates_data$PD_5 = ((replicates_data$constantmetrate_mass / replicates_data$initial_mass) - 1) * 100
replicates_data$PD_6 = ((replicates_data$constantmetrate_mass / replicates_data$final_mass) - 1) * 100

# Metabolic rate PD comparison for each observationn
barplot(replicates_data$PD_1, col = "green", ylab = "MR change from initial")
barplot(replicates_data$PD_2, col = rgb(0, 0, 0, alpha = 0.3), add = TRUE)
legend("topright", c("Actual", "Constant mass"), fill = c("green", "white"))

# Metabolic rate PD comparison for actual and constant mass conditions
boxplot(replicates_data$PD_1, replicates_data$PD_2, main = "Mass change mediates MR increase", xlab = "Mass", ylab = "Metabolic rate PD", names = c("actual (i.e., varying)", "constant"))

# Mass PD comparison for each observationn
barplot(replicates_data$PD_4, col = "orange", ylab = "Mass change from initial")
barplot(replicates_data$PD_5, col = rgb(0, 0, 0, alpha = 0.3), add = TRUE)
legend("topright", c("Actual", "Compensation"), fill = c("orange", "white"))

# Mass PD comparison for actual and constant metabolic rate conditions
boxplot(replicates_data$PD_4, replicates_data$PD_5, main = "Required mass change to maintain MR", xlab = "Metabolic rate", ylab = "Mass PD", names = c("actual", "constant (i.e., compensation mass)"))

# ------------------Mixed model----------

library(lme4)

replicates_data$releveled_class = relevel(factor(replicates_data$Class), "Insecta")

replicates_data$paper_id = gsub("[A-Z]$", "", as.character(replicates_data$studyID))

summary(lmer(mass_comparison ~ (1|studyID), data = replicates_data))
summary(lmer(mass_comparison ~ (1|species), data = replicates_data))
summary(lmer(mass_comparison ~ releveled_class + (1|species), data = replicates_data))

summary(lmer(mass_comparison ~ releveled_class + (1|species) + (1|paper_id), data = replicates_data))

# Don't publish based on these numbers
round(ranef(lmer(mass_comparison ~ releveled_class + (1|species), data = replicates_data))[[1]], 3)

#---------------------SPECIES VALUES------------------------

# Create species dataset
# This code is way awful
library(dplyr)

# Get columns for averages
species_data = replicates_data %>%
  group_by(species) %>%
  #summarise(count = n())
  summarise_each(funs(mean), initial_mass, final_mass, initial_metrate, final_metrate, 
                 constantmass_metrate, constantmetrate_mass, PD_1, PD_2, PD_3, PD_4, PD_5, PD_6) 

# Add in species occurrence values
occurrences = replicates_data %>%
  group_by(species) %>%
  summarise(count = n())

species_data$occurrences = occurrences$count

# Density plots for t-test for metabolic rate + normality test + t-test
plot(density(species_data$PD_1), col = "black", ylim = c(0, 0.19), main = "Comparison of MR PDs", xlab = "Percent Difference", ylab = "PD Density")
legend("topright", c("1: initial to final", "2: initial to constant mass"), fill = c("black", "green"))
lines(density(species_data$PD_2), col = "green")
abline(v = 0, lty = 2, col = "red")
shapiro.test(species_data$PD_1)
shapiro.test(species_data$PD_2)
t.test(species_data$PD_1, species_data$PD_2, paired = TRUE)

# Density plots for t-test for mass + normality test + t-test
plot(density(species_data$PD_4), col = "black", xlim = c(-40, 40), ylim = c(0, 0.25), main = "Comparison of Mass PDs", xlab = "Percent Difference", ylab = "PD Density")
lines(density(species_data$PD_5), col = "orange")
legend("topright", c("4: initial to final", "5: inital to constant MR"), fill = c("black", "orange"))
abline(v = 0, lty = 2, col = "red")
shapiro.test(species_data$PD_4)
shapiro.test(species_data$PD_5)
t.test(species_data$PD_4, species_data$PD_5, paired = TRUE)

# Add in class column
# FIXME: lookup is a dataframe in the cleaning R script
#lookup_test = lookup[match(species_data$species, replicates_data$species),]
#species_data$class = lookup_test$Class
