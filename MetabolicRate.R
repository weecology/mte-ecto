
#--------------CALCULATING EXPONENTS AND ACTIVATION ENERGY FOR METABOLIC EQUATIONS-------------------------------

### For the metabolic rate equation Metabolic Rate = constant * (Mass ^exponent) * e ^ Activation Energy/kt
### This section uses metabolic rate, temperature, and size data to calculated fitted exponent and activation
### energy for the taxonomic classes for which we could find data. 

# Classes: Branchiopoda, Malacostraca, Maxillopoda
# Data source: Makrieva
MR_data = read.csv("Class_metabolicrates_Makrievadata.csv", stringsAsFactors = FALSE)
Class_list=unique(MR_data$Class)
InvK=1/(273.15+MR_data$TC) #getting all the units correct
Mkg=MR_data$Mg/1000
W=MR_data$qWkg*Mkg
MR_data = cbind(MR_data, InvK, W)

#obtaining class-specific exponents for metabolic rate scaling for Branchiopoda Malacostraca and Maxillopoda
class_slopes=c()
for (current_class in Class_list){
  Class_data=subset(MR_data, MR_data$Class == current_class)
  basic=lm(log(Class_data$W) ~ log(Class_data$Mg) + Class_data$InvK)
  Mslope=basic$coefficients[2]
  Tslope=basic$coefficients[3]
  Ea=Tslope*0.00008617
  class_slopes=rbind(class_slopes, c(current_class, Mslope, Tslope, Ea))
}

# Classes: Amphibia
# Data source: White
MR_Amphi = read.csv("Whiteetal_Amphibiandata.csv", stringsAsFactors = FALSE)
A_InvK=1/(273.15+MR_Amphi$TC)
A_Mkg=MR_Amphi$Mg/1000
MR_Amphi = cbind(MR_Amphi, A_InvK, A_Mkg)
A_basic=lm(log(MR_Amphi$Watts) ~ log(MR_Amphi$A_Mkg) + MR_Amphi$A_InvK)
Mslope=A_basic$coefficients[2]
Tslope=A_basic$coefficients[3]
Ea=Tslope*0.00008617
class_slopes=rbind(class_slopes, c("Amphibia", Mslope, Tslope, Ea))
class_slopes=as.data.frame(class_slopes, stringsAsFactors = FALSE)

# Classes: Actinoperygii
# Data source: Gillooly
#obtaining super-class-specific exponents for metabolic rate scaling in fish
MR_fish = read.csv("gillooly_fish.csv", stringsAsFactors = FALSE)
f_Mkg=MR_fish$Mg/1000
MR_fish = cbind(MR_fish, f_Mkg)
F_basic=lm(log(MR_fish$W) ~ log(MR_fish$f_Mkg) + MR_fish$invK)
Mslope=F_basic$coefficients[2]
Tslope=F_basic$coefficients[3]
Ea=Tslope*0.00008617
class_slopes=rbind(class_slopes, c("Actinoperygii", Mslope, Tslope, Ea))
class_slopes=as.data.frame(class_slopes, stringsAsFactors = FALSE)

#removing unneccesary column and renaming columns
class_values=class_slopes[c(-3)]
names(class_values)[names(class_values) == "V1"] = "Class"
names(class_values)[names(class_values) == "log(Class_data$Mg)"] = "exponent"
names(class_values)[names(class_values) == "Class_data$InvK"] = "Ea"

# Classes: Insecta, Gastropoda, Eurotatoria, Entognatha
# Data source (Insecta only): Chown
# Remaining classes' values from theory
I_class="Insecta"
other_class=c("Gastropoda","Eurotatoria", "Entognatha")
I_slope=0.75
I_Ea=-0.62
Other_Ea=-0.63
other_values=cbind(I_class,I_slope,I_Ea)
for (current_class in other_class){
  current_values=cbind(current_class,I_slope,Other_Ea)
  other_values=rbind(other_values, current_values)
}

other_values = as.data.frame(other_values)
names(other_values)[names(other_values) == "I_class"] = "Class"
names(other_values)[names(other_values) == "I_slope"] = "exponent"
names(other_values)[names(other_values) == "I_Ea"] = "Ea"
class_values = rbind(class_values, other_values)
class_values = transform(class_values, exponent = as.numeric(exponent), Ea = as.numeric(Ea))

#---------------------CREATE THREE DEGREE DATASET----------------------------
### Within a scientific paper, the body size-temperature study may have replicated
### experiments or have slightly different experiments that were conducted. Each experiment
### is called a replicate and was given a unique identifier in "MTEEcto_data.csv". 

# Read in compiled experimental dataset, has one temperature and replicate per row
original_data = read.csv("MTEEcto_data.csv", stringsAsFactors = FALSE, fileEncoding = "latin1")

# Turn unique row names into column in original dataset
original_data = data.frame(as.numeric(rownames(original_data)), original_data)
colnames(original_data)[1] = "row_names"

# Create list of experiment names and respective species
replicate_list_original=unique(original_data$studyID)

# Create dataset containing all temperature pairs within each replicate
all_pairs = c()
for(current_replicate in replicate_list_original){
  single_replicate = subset(original_data, original_data$studyID == current_replicate)
  single_replicate_rownames_combos = combn(single_replicate$row_names, 2)
  single_replicate_rownames_combos = t(single_replicate_rownames_combos)
  single_replicate_temps_combos = combn(single_replicate$temp, 2)
  single_replicate_temps_combos = t(single_replicate_temps_combos)
  single_replicate_masses_combos = combn(single_replicate$mass, 2)
  single_replicate_masses_combos = t(single_replicate_masses_combos)
  single_replicate_combos = cbind(single_replicate_rownames_combos, single_replicate_temps_combos, single_replicate_masses_combos)
  single_replicate_combos = data.frame(single_replicate_combos)
  single_replicate_combos = cbind(current_replicate, single_replicate_combos)
  all_pairs = rbind(all_pairs, single_replicate_combos)
}
names(all_pairs) = c("studyID", "initial_temp_row", "final_temp_row", "initial_temp", "final_temp", "initial_mass", "final_mass")

# Subset dataset to contain only temperature pairs that differ by 3*
all_pairs$temp_diff = all_pairs$final_temp - all_pairs$initial_temp
three_degree_pairs = subset(all_pairs, all_pairs$temp_diff == 3)

# Add species to three degree difference pairs
lookup = original_data[match(three_degree_pairs$studyID, original_data$studyID),]
three_degree_pairs$species = lookup$Species
three_degree_pairs$Class = lookup$Class

# List of species and number of pairs per species
species_duplicates = sort(table(three_degree_pairs$species), decreasing=TRUE)

# Add class-specific values to pairs dataset
store_class_values = class_values[match(three_degree_pairs$Class, class_values$Class),]
three_degree_pairs = cbind(three_degree_pairs, store_class_values$exponent, store_class_values$Ea)
colnames(three_degree_pairs)[11] = "exponent"
colnames(three_degree_pairs)[12] = "Ea"

# List of classes and number of pairs per class
class_distribution = table(three_degree_pairs$Class)

#------------------CALCULATE METABOLIC RATES/MASSES--------------------------

# Use MTE equation to calculate metabolic rates (initial, final, and w/ constant/initial mass)
three_degree_pairs$initial_metrate = (three_degree_pairs$initial_mass ^ three_degree_pairs$exponent) * (exp(((three_degree_pairs$Ea) / (.00008617 * (three_degree_pairs$initial_temp + 273.15)))))
three_degree_pairs$final_metrate = (three_degree_pairs$final_mass ^ three_degree_pairs$exponent) * (exp(((three_degree_pairs$Ea) / (.00008617 * (three_degree_pairs$final_temp + 273.15)))))
three_degree_pairs$constantmass_metrate = (three_degree_pairs$initial_mass ^ three_degree_pairs$exponent) * (exp(((three_degree_pairs$Ea) / (.00008617 * (three_degree_pairs$final_temp + 273.15)))))

# Use rearranged MTE equation to calculate mass w/ constant/initial metabolic rate
three_degree_pairs$constantmetrate_mass = (three_degree_pairs$initial_metrate / (exp(three_degree_pairs$Ea / (.00008617 * (three_degree_pairs$final_temp + 273.15))))) ^ (1 / three_degree_pairs$exponent)

#----------------CALCULATE LOG TRANSFORMED VALUES-----------------------------

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

#--------------------EXPLORATORY STATS AND FIGURES-----------------------


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

### Summary values

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

### Figures

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

### Statistical tests
