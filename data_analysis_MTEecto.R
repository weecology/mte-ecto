# Read in clean data
pairs_data = read.csv("clean_data_MTEecto.csv", stringsAsFactors = FALSE)

#----------------METABOLIC RATE VISUALIZATION & ANALYSIS FOR PAIRS-------------------

# Metabolic rate
pairs_data$PD_MR_CS = ((pairs_data$constantmass_metrate / pairs_data$initial_metrate) - 1) * 100
pairs_data$PD_MR_VS = ((pairs_data$final_metrate / pairs_data$initial_metrate) - 1) * 100
pairs_data$MR_diff = pairs_data$PD_MR_CS - pairs_data$PD_MR_VS
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

# Full model but w/ size as fixed effect
# initial_mass vs. initial_mass_mg ?
full_model_size = lmer(log_MR_diff ~ initial_temp + initial_mass_mg + (1|species) + (1|Class), data = pairs_data)
summary(full_model_size)

# Likelihood ratio test for temp w/ linear regression plot & by species variation
null_temp = lmer(log_MR_diff ~ (1|species) + (1|Class), data = pairs_data, REML = FALSE)
model_temp = lmer(log_MR_diff ~ initial_temp + (1|species) + (1|Class), data = pairs_data, REML = FALSE)
anova(null_temp, model_temp)
par(mfrow = c(1, 1))
lr = lm(pairs_data$log_MR_diff ~ pairs_data$initial_temp)
plot(pairs_data$initial_temp, pairs_data$log_MR_diff)
abline(lr)

pairs_data = pairs_data[order(pairs_data$species),]
library(dplyr)
by_species = group_by(pairs_data, species)
temperatures = summarise(by_species, 
                         pairs_count = n(),
                         minimum = min(initial_temp), 
                         maximum = max(initial_temp), 
                         average = mean(initial_temp), 
                         standard_dev = sd(initial_temp))
barplot(temperatures$maximum)
barplot(temperatures$minimum, col = "red", add = TRUE)

# Likelihood ratio test for size w/ lin reg plot & by species variation
null_size = lmer(log_MR_diff ~ initial_temp + (1|species) + (1|Class), data = pairs_data, REML = FALSE)
model_size = lmer(log_MR_diff ~ initial_temp + initial_mass_mg + (1|species) + (1|Class), data = pairs_data, REML = FALSE)
anova(null_size, model_size)
lr_size = lm(pairs_data$log_MR_diff ~ pairs_data$initial_mass_mg)
plot(pairs_data$initial_mass_mg, pairs_data$log_MR_diff)
abline(lr_size)

sizes = summarise(by_species, 
                  minimum = min(initial_mass), 
                  maximum = max(initial_mass), 
                  average = mean(initial_mass), 
                  standard_dev = sd(initial_mass))
barplot(log(sizes$maximum))
barplot(log(sizes$minimum), col = "red", add = TRUE)

# Likelihood ratio test for species
null_sp = lmer(log_MR_diff ~ initial_temp + (1|Class), data = pairs_data, REML = FALSE)
model_sp = lmer(log_MR_diff ~ initial_temp + (1|species) + (1|Class), data = pairs_data, REML = FALSE)
anova(null_sp, model_sp)

# Likelihood ratio test for class
null_class = lmer(log_MR_diff ~ initial_temp + (1|species), data = pairs_data, REML = FALSE)
model_class = lmer(log_MR_diff ~ initial_temp + (1|species) + (1|Class), data = pairs_data, REML = FALSE)
anova(null_class, model_class)

#---------------METABOLIC RATE VISUALIZATION & ANALYSIS FOR SPECIES----------

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

# Average metabolic rates by species
by_species = group_by(pairs_data, species, Class)
species_data = summarise(by_species, 
                         pairs_count = n(), 
                         initial_metrate = mean(initial_metrate), 
                         final_metrate = mean(final_metrate), 
                         constantmass_metrate = mean(constantmass_metrate))

# Calculations and percent difference plots
species_data$log_MR_diff = log(species_data$constantmass_metrate) - log(species_data$final_metrate)
species_data$PD_MR_CS = (species_data$constantmass_metrate - species_data$initial_metrate) / abs(species_data$initial_metrate) * 100
species_data$PD_MR_VS = (species_data$final_metrate - species_data$initial_metrate) / abs(species_data$initial_metrate) * 100

library(tidyr)
library(ggplot2)
species_data_plot = species_data[c("species", "Class", "PD_MR_CS", "PD_MR_VS")]
species_data_plot = gather(species_data_plot, key, value, -species, -Class)
species_data_plot$class_abbr = as.numeric(factor(species_data_plot$Class))
ggplot(species_data_plot, aes(x = species, y = value, fill = factor(key))) +
  geom_bar(stat = "identity", position = "identity") + 
  facet_grid(~class_abbr, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  scale_fill_manual(values = c(rgb(0, 0, 1, alpha = 0.5), rgb(1, 0, 0, alpha = 0.5)), 
                    labels = c("Constant-mass metabolic rate", "Varying-mass metabolic rate")) +
  guides(fill = guide_legend(title = NULL)) + 
  theme(legend.position = "top") + 
  scale_y_continuous(name = "Percent change from initial mass (%)") + scale_x_discrete(name = "Species")

# T-test for comparing means of varying size and constant size metabolic rates
par(mfrow = c(1, 1))
hist(species_data$log_MR_diff) #normality of difference between pairs
t.test(log(species_data$constantmass_metrate), log(species_data$final_metrate), paired = TRUE)

#------------------MASS ANALYSIS & VISUALIZATION---------------

### Mass ratio for 1:1 plot
mass_averages = pairs_data %>%
  group_by(species) %>%
  summarise_each(funs(mean), initial_mass, final_mass, constantmetrate_mass)

species_data = merge(x = species_data, y = mass_averages, by = "species", all.x = TRUE)
species_data$needed_mass_ratio = species_data$constantmetrate_mass / species_data$initial_mass
species_data$actual_mass_ratio = species_data$final_mass / species_data$initial_mass

# Plot
plot(species_data$needed_mass_ratio, species_data$actual_mass_ratio, xlim = c(0.63, 1.37), ylim = c(0.63, 1.37), pch = 20, xlab = "Needed mass ratio", ylab = "Actual mass ratio")
lines(x = c(0, 5), y = c(0, 5))
abline(h = 1, v = 1, col = "red", lty = 2)

# Temporary plot
plot(density(log(species_data$constantmetrate_mass)))
lines(density(log(species_data$final_mass)), col = "red")

# T test for comparing means of actual and needed mass
t.test(log(species_data$final_mass), log(species_data$constantmetrate_mass), paired = TRUE)

