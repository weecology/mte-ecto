# Read in clean data
pairs_data = read.csv("clean_data_MTEecto.csv", stringsAsFactors = FALSE)

#---------------------------PAIRS-------------------

# Metabolic rate
pairs_data$PD_MR_CS = ((pairs_data$constantmass_metrate / pairs_data$initial_metrate) - 1) * 100
pairs_data$PD_MR_VS = ((pairs_data$final_metrate / pairs_data$initial_metrate) - 1) * 100
pairs_data$MR_diff = pairs_data$PD_MR_CS - pairs_data$PD_MR_VS
pairs_data$PD_mass = ((pairs_data$final_mass / pairs_data$initial_mass) - 1) * 100
pairs_data = pairs_data[order(pairs_data$MR_diff),]
barplot(pairs_data$PD_MR_VS, col = "green", xlab = "Pairs", ylab = "MR change due to temp increase")
barplot(pairs_data$PD_MR_CS, col = rgb(0, 0, 0, alpha = 0.3), add = TRUE)
legend("topright", c("constant mass/no STR", "varying mass/STR"), fill = c("white", "green"))

pairs_data = pairs_data[order(pairs_data$initial_temp),]
barplot(pairs_data$PD_MR_VS, col = "green", xlab = "Pairs", ylab = "MR change due to temp increase")
barplot(pairs_data$PD_MR_CS, col = rgb(0, 0, 0, alpha = 0.3), add = TRUE)
legend("topright", c("constant mass/no STR", "varying mass/STR"), fill = c("white", "green"))

pairs_data$log_MR_diff = log(pairs_data$constantmass_metrate) - log(pairs_data$final_metrate)

#TODO: MR PDs for all pairs grouped by species
sp_names = data.frame(table(pairs_data$species))
sp_names = sp_names[order(sp_names$Freq, decreasing = TRUE),]

par(mfrow = c(3, 5))
for(name in sp_names$Var1){
  sp_df = pairs_data[pairs_data$species == name,]
  barplot(sp_df$PD_MR_VS, col = "green", main = name, ylim = c(15, 90))
  barplot(sp_df$PD_MR_CS, col = rgb(0, 0, 0, alpha = 0.3), add = TRUE)
  abline(h = 0)
}

library(lme4)
#removed random effects study and studyID because they explained no variability
full_model = lmer(log_MR_diff ~ initial_temp + (1|species) + (1|Class), data = pairs_data)
summary(full_model)

# Likelihood Ratio Test for temp
null_temp = lmer(log_MR_diff ~ (1|species) + (1|Class), data = pairs_data, REML = FALSE)
model_temp = lmer(log_MR_diff ~ initial_temp + (1|species) + (1|Class), data = pairs_data, REML = FALSE)
anova(null_temp, model_temp)

par(mfrow = c(1, 1))
lr = lm(pairs_data$log_MR_diff ~ pairs_data$initial_temp)
plot(pairs_data$initial_temp, pairs_data$log_MR_diff)
abline(lr)

# Likelihood Ratio Test for species
null_sp = lmer(log_MR_diff ~ initial_temp + (1|Class), data = pairs_data, REML = FALSE)
model_sp = lmer(log_MR_diff ~ initial_temp + (1|species) + (1|Class), data = pairs_data, REML = FALSE)
anova(null_sp, model_sp)

# Likelihood Ratio Test for class
null_class = lmer(log_MR_diff ~ initial_temp + (1|species), data = pairs_data, REML = FALSE)
model_class = lmer(log_MR_diff ~ initial_temp + (1|species) + (1|Class), data = pairs_data, REML = FALSE)
anova(null_class, model_class)

#---------------METABOLIC RATE ANALYSIS & VISUALIZATION----------

### Metabolic rate percent differences

library(dplyr)
species_data = pairs_data %>% 
  group_by(species) %>% 
  summarise_each(funs(mean, sd), PD_MR_CS, PD_MR_VS)

occurrences = pairs_data %>% 
  group_by(species) %>% 
  summarise(count = n())
species_data$occurrences = occurrences$count

# Looking at variation in percent differences for each species w/ sufficient pairs
species_data = species_data[order(species_data$occurrences, decreasing = TRUE),]
sub_species_data = species_data[1:4,]

par(mfrow = c(2, 2))
sp1 = pairs_data[pairs_data$species == "Drosophila melanogaster  ",]
sp1$PD_mass = ((sp1$final_mass / sp1$initial_mass) - 1) * 100
sp1 = sp1[order(sp1$PD_mass, decreasing = TRUE),]
sp1$MR_diff = sp1$PD_MR_VS - sp1$PD_MR_CS
barplot(sp1$PD_MR_VS, col = "green", main = "Drosophila melanogaster")
barplot(sp1$PD_MR_CS, col = rgb(0, 0, 0, alpha = 0.3), add = TRUE)
abline(h = mean(sp1$PD_MR_VS), col = "green")
abline(h = mean(sp1$PD_MR_CS))

sp2 = pairs_data[pairs_data$species == "Drosophila simulans  ",]
sp2$PD_mass = ((sp2$final_mass / sp2$initial_mass) - 1) * 100
sp2 = sp2[order(sp2$PD_mass, decreasing = TRUE),]
sp2$MR_diff = sp2$PD_MR_VS - sp2$PD_MR_CS
barplot(sp2$PD_MR_VS, col = "green", main = "Drosophila simulans")
barplot(sp2$PD_MR_CS, col = rgb(0, 0, 0, alpha = 0.3), add = TRUE)
abline(h = mean(sp2$PD_MR_VS), col = "green")
abline(h = mean(sp2$PD_MR_CS))

sp3 = pairs_data[pairs_data$species == "Rana sylvatica",]
sp3$PD_mass = ((sp3$final_mass / sp3$initial_mass) - 1) * 100
sp3 = sp3[order(sp3$PD_mass, decreasing = TRUE),]
sp3$MR_diff = sp3$PD_MR_VS - sp3$PD_MR_CS
barplot(sp3$PD_MR_VS, col = "green", main = "Rana sylvatica")
barplot(sp3$PD_MR_CS, col = rgb(0, 0, 0, alpha = 0.3), add = TRUE)
abline(h = mean(sp3$PD_MR_VS), col = "green")
abline(h = mean(sp3$PD_MR_CS))

sp4 = pairs_data[pairs_data$species == "Acyrthosiphon pisum ",]
sp4$PD_mass = ((sp4$final_mass / sp4$initial_mass) - 1) * 100
sp4 = sp4[order(sp4$PD_mass, decreasing = TRUE),]
sp4$MR_diff = sp4$PD_MR_VS - sp4$PD_MR_CS
barplot(sp4$PD_MR_VS, col = "green", main = "Acyrthosiphon pisum")
barplot(sp4$PD_MR_CS, col = rgb(0, 0, 0, alpha = 0.3), add = TRUE)
abline(h = mean(sp4$PD_MR_VS), col = "green")
abline(h = mean(sp4$PD_MR_CS))

par(mfrow = c(2, 2))
for (species_name in unique(sub_species_data$species)){
  sub = filter(pairs_data, species == species_name)
  hist(sub$PD_MR_VS, main = c(species_name), col = "green")
  abline(v = mean(sub$PD_MR_VS), col = "green")
  hist(sub$PD_MR_CS, col = rgb(0, 0, 0, alpha = 0.3), add = TRUE)
  abline(v = mean(sub$PD_MR_CS))
}

par(mfrow = c(1, 2))
bp_VS = barplot(sub_species_data$PD_MR_VS_mean)
arrows(bp_VS, sub_species_data$PD_MR_VS_mean + sub_species_data$PD_MR_VS_sd, bp_VS, sub_species_data$PD_MR_VS_mean, angle = 90, code = 0)
arrows(bp_VS, sub_species_data$PD_MR_VS_mean - sub_species_data$PD_MR_VS_sd, bp_VS, sub_species_data$PD_MR_VS_mean, angle = 90, code = 0)
bp_CS = barplot(sub_species_data$PD_MR_CS_mean)
arrows(bp_CS, sub_species_data$PD_MR_CS_mean + sub_species_data$PD_MR_CS_sd, bp_CS, sub_species_data$PD_MR_CS_mean, angle = 90, code = 0)
arrows(bp_CS, sub_species_data$PD_MR_CS_mean - sub_species_data$PD_MR_CS_sd, bp_CS, sub_species_data$PD_MR_CS_mean, angle = 90, code = 0)

# Barplot
par(mfrow = c(1, 1))
species_data$PD_MR_diff = species_data$PD_MR_CS_mean - species_data$PD_MR_VS_mean
species_data = species_data[order(species_data$PD_MR_diff),]
barplot(species_data$PD_MR_VS_mean, col = "green", xlab = "Species", ylab = "MR change due to temp increase")
barplot(species_data$PD_MR_CS_mean, col = rgb(0, 0, 0, alpha = 0.3), add = TRUE)
legend("topright", c("constant mass/no STR", "varying mass/STR"), fill = c("white", "green"))

# T-test or something? 
#TODO: paired t-test comparing constant mass and varying mass MR PDs


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
