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

asterisks = species_data %>%
  group_by(species) %>%
  select(species, PD_MR_CS, PD_MR_VS) %>%
  mutate(value_real = ifelse(PD_MR_CS > PD_MR_VS, PD_MR_CS, PD_MR_VS), 
         color = ifelse(PD_MR_CS > PD_MR_VS, "PD_MR_CS", "PD_MR_VS"))
asterisks = merge(x = asterisks, y = species_data_plot[,c("class_abbr", "species")], by = "species", all.x = TRUE)
asterisks$value = rep(-3, nrow(asterisks))

ggplot(species_data_plot, aes(x = species, y = value)) +
  geom_bar(aes(fill = key), position = "dodge", stat = "identity") +
  theme_bw() +
  facet_grid(~class_abbr, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c(rgb(0, 0, 0), rgb(0, 1, 0)), 
                    labels = c("Constant-mass metabolic rate", "Varying-mass metabolic rate")) +
  scale_y_continuous(name = "Non-STR & STR metabolic rate changes (%)") + scale_x_discrete(name = "Species") +
  geom_text(aes(colour = factor(color)), data = asterisks, label = "*", size = 10) +
  scale_colour_manual(values = c(rgb(0, 0, 0), rgb(0, 1, 0))) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 65)) +
  theme(strip.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 13), 
        axis.title = element_text(size = 16))

class_abbreviation = unique(species_data_plot[,c("Class", "class_abbr")])
class_abbreviation = class_abbreviation[order(class_abbreviation$class_abbr),]

# T-test for comparing means of varying size and constant size metabolic rates
par(mfrow = c(1, 1))
hist(species_data$log_MR_diff) #normality of difference between pairs
t.test(log(species_data$constantmass_metrate), log(species_data$final_metrate), paired = TRUE)

#------------------MASS ANALYSIS & VISUALIZATION---------------

### Mass percentage for 1:1 plot
mass_averages = pairs_data %>%
  group_by(species) %>%
  summarise_each(funs(mean), initial_mass, final_mass, constantmetrate_mass)

species_data = merge(x = species_data, y = mass_averages, by = "species", all.x = TRUE)
species_data$needed_mass_percent = species_data$constantmetrate_mass / species_data$initial_mass * 100
species_data$actual_mass_percent = species_data$final_mass / species_data$initial_mass * 100

# Plot
plot(50, 50, xlim = c(60, 137), ylim = c(60, 137), xlab = "", ylab = "")
mtext("Compensation mass change", side = 1, line = 2.5, font = 2, cex = 1.5)
mtext("(needed mass/initial mass) * 100", side = 1, line = 3.5, cex = 0.75)
mtext("Observed mass change", side = 2, line = 2.8, font = 2, cex = 1.5)
mtext("(actual mass/initial mass) * 100", side = 2, line = 2, cex = 0.75)
polygon(c(50, 50, 100), c(50, 100, 100), col = rgb(0, 1, 0, alpha = 0.3), border = NA)
polygon(c(50, 100, 100, 50), c(100, 100, 140, 140), col = rgb(0, 1, 0), border = NA)
points(species_data$needed_mass_percent, species_data$actual_mass_percent, pch = 20)
lines(x = c(0, 500), c(0, 500))
abline(h = 100, v = 100, lty = 3)
text(113.9082, 128, "mass change for temp-\n erature to have no effect\n on metabolic rate", cex = 0.9, xpd = TRUE)
arrows(x0 = 113.9818, y0 = 122.8864, x1 = 113.7489, y1 = 115.0543, length = 0.08)
text(130.9794, 91, "mass does\n not change")
arrows(x0 = 130.6321, y0 = 94.34712, x1 = 130.6321, y1 = 99.01680, length = 0.08)
text(112, 64.21285, "mass does not\n need to change")
arrows(x0 = 105.2136, y0 = 64.41466, x1 = 101.9626, y1 = 64.41466, length = 0.08)

# T test for comparing means of actual and needed mass
t.test(log(species_data$final_mass), log(species_data$constantmetrate_mass), paired = TRUE)

