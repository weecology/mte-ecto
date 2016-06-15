# Read in clean data
pairs_data = read.csv("clean_data_MTEecto.csv", stringsAsFactors = FALSE)

#---------------------------PAIRS VISUALIZATION & ANALYSIS-------------------

# Metabolic rate
pairs_data$PD_MR_CS = ((pairs_data$constantmass_metrate / pairs_data$initial_metrate) - 1) * 100
pairs_data$PD_MR_VS = ((pairs_data$final_metrate / pairs_data$initial_metrate) - 1) * 100
pairs_data$MR_diff = pairs_data$PD_MR_CS - pairs_data$PD_MR_VS
pairs_data$PD_mass = ((pairs_data$final_mass / pairs_data$initial_mass) - 1) * 100
pairs_data = pairs_data[order(pairs_data$MR_diff),]
barplot(pairs_data$PD_MR_VS, col = "green", xlab = "Pairs", ylab = "MR change due to temp increase")
barplot(pairs_data$PD_MR_CS, col = rgb(0, 0, 0, alpha = 0.3), add = TRUE)
legend("topright", c("constant mass/no STR", "varying mass/STR"), fill = c("white", "green"))

# Linear mixed model for difference in metabolic rates
pairs_data$log_MR_diff = log(pairs_data$constantmass_metrate) - log(pairs_data$final_metrate)
library(lme4)
initial_model = lmer(log_MR_diff ~ initial_temp + (1|species) + (1|Class) + (1|studyID) + (1|study), data = pairs_data)
summary(initial_model)
#removed random effects study and studyID because they explained no variability
full_model = lmer(log_MR_diff ~ initial_temp + (1|species) + (1|Class), data = pairs_data)
summary(full_model)
coef(full_model)
hist(resid(full_model))
plot(fitted(full_model), residuals(full_model)) #linearity good, heteroscedasticity
abline(h = 0)
hist(residuals(full_model))
qqnorm(residuals(full_model))
qqline(residuals(full_model)) #normality of residuals not good, but that's okay? 

# Likelihood ratio test for temp w/ linear regression plot
null_temp = lmer(log_MR_diff ~ (1|species) + (1|Class), data = pairs_data, REML = FALSE)
model_temp = lmer(log_MR_diff ~ initial_temp + (1|species) + (1|Class), data = pairs_data, REML = FALSE)
anova(null_temp, model_temp)
par(mfrow = c(1, 1))
lr = lm(pairs_data$log_MR_diff ~ pairs_data$initial_temp)
plot(pairs_data$initial_temp, pairs_data$log_MR_diff)
abline(lr)

# Likelihood ratio test for species
null_sp = lmer(log_MR_diff ~ initial_temp + (1|Class), data = pairs_data, REML = FALSE)
model_sp = lmer(log_MR_diff ~ initial_temp + (1|species) + (1|Class), data = pairs_data, REML = FALSE)
anova(null_sp, model_sp)

# Likelihood ratio test for class
null_class = lmer(log_MR_diff ~ initial_temp + (1|species), data = pairs_data, REML = FALSE)
model_class = lmer(log_MR_diff ~ initial_temp + (1|species) + (1|Class), data = pairs_data, REML = FALSE)
anova(null_class, model_class)

#---------------METABOLIC RATE ANALYSIS & VISUALIZATION----------

### Metabolic rate percent differences

#MR PDs for all pairs grouped by species
sp_names = data.frame(table(pairs_data$species))
sp_names = sp_names[order(sp_names$Freq, decreasing = TRUE),]

par(mfrow = c(3, 5))
for(name in sp_names$Var1){
  sp_df = pairs_data[pairs_data$species == name,]
  barplot(sp_df$PD_MR_VS, col = "green", main = name, ylim = c(15, 90))
  barplot(sp_df$PD_MR_CS, col = rgb(0, 0, 0, alpha = 0.3), add = TRUE)
  abline(h = 0)
}

# Average pairs' PD values to get species PD values
library(dplyr)
species_data = pairs_data %>% 
  group_by(species) %>% 
  summarise_each(funs(mean, sd), PD_MR_CS, PD_MR_VS, initial_metrate, final_metrate, constantmass_metrate)

occurrences = pairs_data %>% 
  group_by(species) %>% 
  summarise(count = n())
species_data$occurrences = occurrences$count

# Barplot
par(mfrow = c(1, 1))
species_data$PD_MR_diff = species_data$PD_MR_CS_mean - species_data$PD_MR_VS_mean
species_data = species_data[order(species_data$PD_MR_diff),]
barplot(species_data$PD_MR_VS_mean, col = "green", xlab = "Species", ylab = "MR change due to temp increase")
barplot(species_data$PD_MR_CS_mean, col = rgb(0, 0, 0, alpha = 0.3), add = TRUE)
legend("topright", c("constant mass/no STR", "varying mass/STR"), fill = c("white", "green"))

# # Average pairs' metabolic rates to get species MR values, then do PD
# species_data$PD_MR_CS_calc = ((species_data$constantmass_metrate_mean / species_data$initial_metrate_mean) - 1) * 100
# species_data$PD_MR_VS_calc = ((species_data$final_metrate_mean / species_data$initial_metrate_mean) - 1) * 100
# barplot(species_data$PD_MR_VS_calc, col = rgb(1, 0, 0, alpha = 0.2), xlab = "Species", ylab = "MR change due to temp increase", add = TRUE)
# barplot(species_data$PD_MR_CS_calc, col = rgb(0, 0, 1, alpha = 0.2), add = TRUE)
# legend("topright", c("constant mass/no STR", "varying mass/STR"), fill = c("white", "green"))

#------------------MASS ANALYSIS & VISUALIZATION---------------

#TODO: change PD_mass variable name

### Mass percent difference
pairs_data$PD_mass = ((pairs_data$constantmetrate_mass / pairs_data$initial_mass) - 1) * 100

### Mass differences for 1:1 plot
pairs_data$needed_mass_diff = pairs_data$initial_mass - pairs_data$constantmetrate_mass
pairs_data$actual_mass_diff = pairs_data$initial_mass - pairs_data$final_mass

### Mass ratio for 1:1 plot
pairs_data$needed_mass_ratio = pairs_data$constantmetrate_mass / pairs_data$initial_mass
pairs_data$actual_mass_ratio = pairs_data$final_mass / pairs_data$initial_mass

mass_comparison = replicates_data %>%
  group_by(species) %>%
  summarise_each(funs(mean), PD_mass)
species_data$PD_mass = mass_comparison$PD_mass

diffs = replicates_data %>%
  group_by(species) %>%
  summarise_each(funs(mean), needed_mass_diff, actual_mass_diff)
species_data$needed_mass_diff = diffs$needed_mass_diff
species_data$actual_mass_diff = diffs$actual_mass_diff

ratios = replicates_data %>%
  group_by(species) %>%
  summarise_each(funs(mean), needed_mass_ratio, actual_mass_ratio)
species_data$needed_mass_ratio = ratios$needed_mass_ratio
species_data$actual_mass_ratio = ratios$actual_mass_ratio

# Plot
plot(density(species_data$PD_mass), main = "Needed decrease in mass", xlim = c(-40, 10))
species_data$PD_mass_null = rnorm(45, mean = 0, sd = sd(species_data$PD_mass))
lines(density(species_data$PD_mass_null), col = "red")
legend("topright", c("red = null", "black = actual"))

plot(species_data$needed_mass_diff, species_data$actual_mass_diff, xlim = c(-86000000, 209581910), ylim = c(-86000000, 209581910))
lines(x = c(-86000000, 209581910), y = c(-86000000, 209581910))
abline(h = 0, v = 0, col = "red")

species_data$needed_mass_diff_log = log(species_data$needed_mass_diff)
species_data$actual_mass_diff_log = log(species_data$actual_mass_diff)
plot(species_data$needed_mass_diff_log, species_data$actual_mass_diff_log, xlim = c(5, 20), ylim = c(5, 20))
lines(x = c(5, 20), y = c(5, 20))
abline(h = 0, v = 0, col = "red")

plot(species_data$needed_mass_ratio, species_data$actual_mass_ratio, xlim = c(0.6, 1.1), ylim = c(0.7, 1.35))
lines(x = c(0.6, 1.35), y = c(0.6, 1.35))
abline(h = 1, v = 1, col = "red", lty = 2)
