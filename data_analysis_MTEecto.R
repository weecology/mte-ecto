#-----------------------TRANSFORM MASS & MR VALUES-----------------------------

# Log transform all metabolic rate & mass values
three_degree_pairs$log_initial_mass = log(three_degree_pairs$initial_mass)
three_degree_pairs$log_initial_metrate = log(three_degree_pairs$initial_metrate)
three_degree_pairs$log_final_mass = log(three_degree_pairs$final_temp)
three_degree_pairs$log_final_metrate = log(three_degree_pairs$final_metrate)
three_degree_pairs$log_constantmass_metrate = log(three_degree_pairs$constantmass_metrate)
three_degree_pairs$log_constantmetrate_mass = log(three_degree_pairs$constantmetrate_mass)

# Histograms of transformed values
hist(three_degree_pairs$log_initial_mass)
hist(three_degree_pairs$log_initial_metrate)
hist(three_degree_pairs$log_final_mass)
hist(three_degree_pairs$log_final_metrate)
hist(three_degree_pairs$log_constantmass_metrate)
hist(three_degree_pairs$log_constantmetrate_mass)

# Comparison of final and constant mass transformed metabolic rates
three_degree_pairs$metrate_comparison = three_degree_pairs$log_constantmass_metrate - three_degree_pairs$log_final_metrate
mean(three_degree_pairs$metrate_comparison)
hist(three_degree_pairs$metrate_comparison)
t.test(three_degree_pairs$log_final_metrate, three_degree_pairs$log_constantmass_metrate, paired = TRUE)
boxplot(three_degree_pairs$log_final_metrate, three_degree_pairs$log_constantmass_metrate, main = "Mass change mediates MR increase", xlab = "Mass", ylab = "Log metabolic rate", names = c("actual (i.e., varying)", "constant"))

# Comparison of final and constant MR transformed masses
three_degree_pairs$mass_comparison = three_degree_pairs$log_constantmetrate_mass - three_degree_pairs$log_final_mass
mean(three_degree_pairs$mass_comparison)
hist(three_degree_pairs$mass_comparison)
t.test(three_degree_pairs$log_final_mass, three_degree_pairs$log_constantmetrate_mass, paired = TRUE)
boxplot(three_degree_pairs$log_final_mass, three_degree_pairs$log_constantmetrate_mass, main = "Required mass change to maintain MR", xlab = "Metabolic rate", ylab = "Log mass", names = c("actual", "constant (i.e., compensation mass)"))

#-------------------PERCENT DIFFERENCES OF MASS & MR VALUES--------------------

# Percent differences for metabolic rates
# 1: initial to final (empirical), expect increase
# 2: initial to constant mass (theoretical), expect increase
# 3: final to constant mass, expect increase
three_degree_pairs$PD_1 = ((three_degree_pairs$final_metrate / three_degree_pairs$initial_metrate) - 1) * 100
three_degree_pairs$PD_2 = ((three_degree_pairs$constantmass_metrate / three_degree_pairs$initial_metrate) - 1) * 100
three_degree_pairs$PD_3 = ((three_degree_pairs$constantmass_metrate / three_degree_pairs$final_metrate) - 1) * 100

# Percent differences for masses
# 4: initial to final (empirical), expect decrease
# 5: initial to constant metabolic rate (theoretical), expect decrease
# 6: final to constant metabolic rate, expect decrease
three_degree_pairs$PD_4 = ((three_degree_pairs$final_mass / three_degree_pairs$initial_mass) - 1) * 100
three_degree_pairs$PD_5 = ((three_degree_pairs$constantmetrate_mass / three_degree_pairs$initial_mass) - 1) * 100
three_degree_pairs$PD_6 = ((three_degree_pairs$constantmetrate_mass / three_degree_pairs$final_mass) - 1) * 100

# Metabolic rate PD comparison for each observationn
barplot(three_degree_pairs$PD_1, col = "green", ylab = "MR change from initial")
barplot(three_degree_pairs$PD_2, col = rgb(0, 0, 0, alpha = 0.3), add = TRUE)
legend("topright", c("Actual", "Constant mass"), fill = c("green", "white"))

# Metabolic rate PD comparison for actual and constant mass conditions
boxplot(three_degree_pairs$PD_1, three_degree_pairs$PD_2, main = "Mass change mediates MR increase", xlab = "Mass", ylab = "Metabolic rate PD", names = c("actual (i.e., varying)", "constant"))

# Mass PD comparison for each observationn
barplot(three_degree_pairs$PD_4, col = "orange", ylab = "Mass change from initial")
barplot(three_degree_pairs$PD_5, col = rgb(0, 0, 0, alpha = 0.3), add = TRUE)
legend("topright", c("Actual", "Compensation"), fill = c("orange", "white"))

# Mass PD comparison for actual and constant metabolic rate conditions
boxplot(three_degree_pairs$PD_4, three_degree_pairs$PD_5, main = "Required mass change to maintain MR", xlab = "Metabolic rate", ylab = "Mass PD", names = c("actual", "constant (i.e., compensation mass)"))

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

# ------------------Mixed model----------

library(lme4)

three_degree_pairs$releveled_class = relevel(factor(three_degree_pairs$Class), "Insecta")

three_degree_pairs$paper_id = gsub("[A-Z]$", "", as.character(three_degree_pairs$studyID))

summary(lmer(mass_comparison ~ (1|studyID), data = three_degree_pairs))
summary(lmer(mass_comparison ~ (1|species), data = three_degree_pairs))
summary(lmer(mass_comparison ~ releveled_class + (1|species), data = three_degree_pairs))

summary(lmer(mass_comparison ~ releveled_class + (1|species) + (1|paper_id), data = three_degree_pairs))

# Don't publish based on these numbers
round(ranef(lmer(mass_comparison ~ releveled_class + (1|species), data = three_degree_pairs))[[1]], 3)

#---------------------SPECIES VALUES------------------------

# Create species dataset
# This code is way awful
library(dplyr)

# Get columns for averages
species_data = three_degree_pairs %>%
  group_by(species) %>%
  #summarise(count = n())
  summarise_each(funs(mean), initial_mass, final_mass, initial_metrate, final_metrate, 
                 constantmass_metrate, constantmetrate_mass, PD_1, PD_2, PD_3, PD_4, PD_5, PD_6) 

# Add in species occurrence values
occurrences = three_degree_pairs %>%
  group_by(species) %>%
  summarise(count = n())

species_data$occurrences = occurrences$count

# Add in class column
lookup_test = lookup[match(species_data$species, three_degree_pairs$species),]
species_data$class = lookup_test$Class
