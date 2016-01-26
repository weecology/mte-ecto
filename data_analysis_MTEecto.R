# Read in clean data
replicates_data = read.csv("clean_data_MTEecto.csv")

#--------------------NORMALITY OF MASS & MR VALUES-------------------------

# Look at distributions of all mass and metabolic rate values
par(mfrow = c(3, 2))
hist(replicates_data$initial_mass)
hist(replicates_data$initial_metrate)
hist(replicates_data$final_mass)
hist(replicates_data$final_metrate)
hist(replicates_data$constantmass_metrate)
hist(replicates_data$constantmetrate_mass)

# Generate list of normality test p-values for all mass and MR values
IM_norm = shapiro.test(replicates_data$initial_mass)
IMR_norm = shapiro.test(replicates_data$initial_metrate)
FM_norm = shapiro.test(replicates_data$final_mass)
FMR_norm = shapiro.test(replicates_data$final_metrate)
CMR_norm = shapiro.test(replicates_data$constantmass_metrate)
CM_norm = shapiro.test(replicates_data$constantmetrate_mass)

non_trans = c(IM_norm$p.value, IMR_norm$p.value, FM_norm$p.value, FMR_norm$p.value, CMR_norm$p.value, CM_norm$p.value)

#-----------------NORMALITY OF TRANSFORMED MASS & MR VALUES--------------------

# Log transform all metabolic rate & mass values
replicates_data$log_initial_mass = log(replicates_data$initial_mass)
replicates_data$log_initial_metrate = log(replicates_data$initial_metrate)
replicates_data$log_final_mass = log(replicates_data$final_mass)
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

# Add log transformed values to normality test p-values list
LIM_norm = shapiro.test(replicates_data$log_initial_mass)
LIMR_norm = shapiro.test(replicates_data$log_initial_metrate)
LFM_norm = shapiro.test(replicates_data$log_final_mass)
LFMR_norm = shapiro.test(replicates_data$log_final_metrate)
LCMR_norm = shapiro.test(replicates_data$log_constantmass_metrate)
LCM_norm = shapiro.test(replicates_data$log_constantmetrate_mass)

trans = c(LIM_norm$p.value, LIM_norm$p.value, LFM_norm$p.value, LFMR_norm$p.value, LCMR_norm$p.value, LCM_norm$p.value)
norm_pvalues = data.frame(non_trans, trans)

#------------------DIFFERENCES IN TRANSFORMED FINAL & CONSTANT VALUES--------------

par(mfrow = c(2, 2))
# Comparison of final and constant mass transformed metabolic rates
# Positive values mean MR is greater when mass does not vary
# Bootstrap gets distribution of average log difference to check that value
t.test(replicates_data$log_final_metrate, replicates_data$log_constantmass_metrate, paired = TRUE)
boxplot(replicates_data$log_final_metrate, replicates_data$log_constantmass_metrate, main = "Mass change mediates MR increase", xlab = "Mass", ylab = "Log metabolic rate", names = c("actual (i.e., varying)", "constant"))
replicates_data$metrate_comparison = replicates_data$log_constantmass_metrate - replicates_data$log_final_metrate
mean(replicates_data$metrate_comparison)
hist(replicates_data$metrate_comparison)
MR_bootstrap = replicate(1000, {
  boot_rows = sample.int(nrow(replicates_data), replace = TRUE)
  boot_dataset = replicates_data[boot_rows, ]
  mean(log(boot_dataset$constantmass_metrate) - log(boot_dataset$final_metrate))
}
)
hist(MR_bootstrap)
abline(v = mean(replicates_data$metrate_comparison), lwd = 3, col = "red")

par(mfrow = c(2, 2))
# Comparison of final and constant MR transformed masses
# Negative values mean mass would have to be much smaller for MR to be maintained
# Bootstrap gets distribution of average log difference to check that value
t.test(replicates_data$log_final_mass, replicates_data$log_constantmetrate_mass, paired = TRUE)
boxplot(replicates_data$log_final_mass, replicates_data$log_constantmetrate_mass, main = "Required mass change to maintain MR", xlab = "Metabolic rate", ylab = "Log mass", names = c("actual", "constant (i.e., compensation mass)"))
replicates_data$mass_comparison = replicates_data$log_constantmetrate_mass - replicates_data$log_final_mass
mean(replicates_data$mass_comparison)
hist(replicates_data$mass_comparison)
mass_bootstrap = replicate(1000, {
  boot_rows = sample.int(nrow(replicates_data), replace = TRUE)
  boot_dataset = replicates_data[boot_rows, ]
  mean(log(boot_dataset$constantmetrate_mass) - log(boot_dataset$final_mass))
  }
)
hist(mass_bootstrap)
abline(v = mean(replicates_data$mass_comparison), lwd = 3, col = "red")

#-----------PERCENT DIFFERENCES OF MASS & MR VALUES TO INCLUDE INITIAL----------

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
par(mfrow = c(1, 1))
barplot(replicates_data$PD_1, col = "green", ylab = "MR change from initial")
barplot(replicates_data$PD_2, col = rgb(0, 0, 0, alpha = 0.3), add = TRUE)
legend("topright", c("Actual", "Constant mass"), fill = c("green", "white"))

# Metabolic rate PD comparison for actual and constant mass conditions
boxplot(replicates_data$PD_2, replicates_data$PD_1, main = "Mass change mediates MR increase", ylab = "Percent change in MR", names = c("direct effect of temp", "direct & indirect effects of temp"))
abline(h = 0, col = "red", lty = 2)

# Mass PD comparison for each observationn
barplot(replicates_data$PD_4, col = "orange", ylab = "Mass change from initial")
barplot(replicates_data$PD_5, col = rgb(0, 0, 0, alpha = 0.3), add = TRUE)
legend("topright", c("Actual", "Compensation"), fill = c("orange", "white"))

# Mass PD comparison for actual and constant metabolic rate conditions
boxplot(replicates_data$PD_5, replicates_data$PD_4, main = "Required mass change to maintain MR", ylab = "Percent change in mass", names = c("compensation", "empirical"))
abline(h = 0, col = "red", lty = 2)

# ------------------Mixed model----------

library(lme4)

# LMM including all fixed and random effects of interest
LMM = lmer(mass_comparison ~ initial_temp + initial_mass_grams + Class + (1 | species) + (1 | studyID), data = replicates_data)
summary(LMM)

# # From meeting with Dave
# replicates_data$releveled_class = relevel(factor(replicates_data$Class), "Insecta")
# replicates_data$paper_id = gsub("[A-Z]$", "", as.character(replicates_data$studyID))
# summary(lmer(mass_comparison ~ (1|studyID), data = replicates_data))
# summary(lmer(mass_comparison ~ (1|species), data = replicates_data))
# summary(lmer(mass_comparison ~ releveled_class + (1|species), data = replicates_data))
# summary(lmer(mass_comparison ~ releveled_class + (1|species) + (1|paper_id), data = replicates_data))
# # Don't publish based on these numbers
# round(ranef(lmer(mass_comparison ~ releveled_class + (1|species), data = replicates_data))[[1]], 3)

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

# Add log transformed values to species dataset
species_data$log_initial_mass = log(species_data$initial_mass)
species_data$log_initial_metrate = log(species_data$initial_metrate)
species_data$log_final_mass = log(species_data$final_mass)
species_data$log_final_metrate = log(species_data$final_metrate)
species_data$log_constantmass_metrate = log(species_data$constantmass_metrate)
species_data$log_constantmetrate_mass = log(species_data$constantmetrate_mass)

#---------------------CLASS VALUES------------------------

# Get columns for averages
classes_data = replicates_data %>%
  group_by(Class) %>%
  #summarise(count = n())
  summarise_each(funs(mean), initial_mass, final_mass, initial_metrate, final_metrate, 
                 constantmass_metrate, constantmetrate_mass, PD_1, PD_2, PD_3, PD_4, PD_5, PD_6) 

# Add in species occurrence values
occurrences = replicates_data %>%
  group_by(Class) %>%
  summarise(count = n())

classes_data$occurrences = occurrences$count

replicates_data %>%
  group_by(Class) %>%
  print(Class)

boxplot(replicates_data$PD_2, replicates_data$PD_1, main = "Mass change mediates MR increase", ylab = "Percent change in MR", names = c("direct effect of temp", "direct & indirect effects of temp"))
abline(h = 0, col = "red", lty = 2)

