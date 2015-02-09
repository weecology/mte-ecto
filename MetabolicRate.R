
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

#--------------------EXPLORATORY STATS AND FIGURES-----------------------

### Summary values

# Percent differences for metabolic rates
  # 1: initial to final (empirical), expect increase
  # 2: initial to constant mass (theoretical), expect increase
  # 3: final to constant mass, expect increase
percent_diff1 = ((three_degree_pairs$final_metrate / three_degree_pairs$initial_metrate) - 1) * 100
percent_diff2 = ((three_degree_pairs$constantmass_metrate / three_degree_pairs$initial_metrate) - 1) * 100
percent_diff3 = ((three_degree_pairs$constantmass_metrate / three_degree_pairs$final_metrate) - 1) * 100

# Percent differences for masses
  # 4: initial to final (empirical), expect decrease
  # 5: initial to constant metabolic rate (theoretical), expect decrease
  # 6: final to constant metabolic rate, expect decrease
percent_diff4 = ((three_degree_pairs$final_mass / three_degree_pairs$initial_mass) - 1) * 100
percent_diff5 = ((three_degree_pairs$constantmetrate_mass / three_degree_pairs$initial_mass) - 1) * 100
percent_diff6 = ((three_degree_pairs$constantmetrate_mass / three_degree_pairs$final_mass) - 1) * 100

# Average percent differences
avg_pd_1 = mean(percent_diff1)
avg_pd_2 = mean(percent_diff2)
avg_pd_3 = mean(percent_diff3)
avg_pd_4 = mean(percent_diff4)
avg_pd_5 = mean(percent_diff5)
avg_pd_6 = mean(percent_diff6)

# Standard deviations of percent differences
sd_pd_1 = sd(percent_diff1)
sd_pd_2 = sd(percent_diff2)
sd_pd_3 = sd(percent_diff3)
sd_pd_4 = sd(percent_diff4)
sd_pd_5 = sd(percent_diff5)
sd_pd_6 = sd(percent_diff6)

### Statistical tests

# T-tests for metabolic rates
  # 1: initial to final (empirical)
  # 2: initial to constant mass (theoretical)
  # 3: final to constant mass
ttest1 = t.test(three_degree_pairs$initial_metrate, three_degree_pairs$final_metrate, paired = TRUE)
ttest2 = t.test(three_degree_pairs$initial_metrate, three_degree_pairs$constantmass_metrate, paired = TRUE)
ttest3 = t.test(three_degree_pairs$final_metrate, three_degree_pairs$constantmass_metrate, paired = TRUE)

# T-tests for masses
  # 4: initial to final (empirical)
  # 5: initial to constant metabolic rate (theoretical)
  # 6: final to constant metabolic rate
ttest4 = t.test(three_degree_pairs$initial_mass, three_degree_pairs$final_mass, paired = TRUE)
ttest5 = t.test(three_degree_pairs$initial_mass, three_degree_pairs$constantmetrate_mass, paired = TRUE)
ttest6 = t.test(three_degree_pairs$final_mass, three_degree_pairs$constantmetrate_mass, paired = TRUE)

### Figures

# Histograms not useful
# # Histogram of all metabolic rates
# hist(log(three_degree_pairs$initial_metrate))
# hist(log(three_degree_pairs$final_metrate), col = rgb(0, 0, 1, 0.5), add = T)
# hist(log(three_degree_pairs$constantmass_metrate), col = rgb(1, 0, 0, 0.5), add = T)
# 
# # Histogram of all mass values
# hist(log(three_degree_pairs$initial_mass))
# hist(log(three_degree_pairs$final_mass), col = rgb(0, 0, 1, 0.5), add = T)
# hist(log(three_degree_pairs$constantmetrate_mass), col = rgb(1, 0, 0, 0.5), add = T)

# Density plot of metabolic rates, log-transformed
plot(density(log(three_degree_pairs$initial_metrate)), col = "orange", main = "", 
     xlab = "log(metabolic rate value)", ylab = "metabolic rate density")
lines(density(log(three_degree_pairs$final_metrate)), col = "green")
lines(density(log(three_degree_pairs$constantmass_metrate)), col = "purple")
legend("topright", c("initial", "final", "constant mass"), title = "Metabolic rates:",
       fill = c("orange", "green", "purple"))

# Density plot of mass values, log-transformed
plot(density(log(three_degree_pairs$initial_mass)), col = "orange", main = "", 
     xlab = "log(mass value)", ylab = "mass density")
lines(density(log(three_degree_pairs$final_mass)), col = "green")
lines(density(log(three_degree_pairs$constantmetrate_mass)), col = "purple")
legend("topright", c("initial", "final", "constant met rate"), title = "Masses:",
       fill = c("orange", "green", "purple"))

# 1:1 plots not useful, though mass values one is slightly clearer
# # 1:1 plot of metabolic rates
# plot(log(three_degree_pairs$initial_metrate), log(three_degree_pairs$final_metrate), col = "green")
# points(log(three_degree_pairs$initial_metrate), log(three_degree_pairs$constantmass_metrate), col = "purple")
# lines(x = c(-40, 0), y = c(-40, 0))
# 
# # 1:1 plot of mass values
# plot(log(three_degree_pairs$initial_mass), log(three_degree_pairs$final_mass), col = "green")
# points(log(three_degree_pairs$initial_mass), log(three_degree_pairs$constantmetrate_mass), col = "purple")
# lines(x = c(-7, 7), y = c(-7, 7))

# Barplot of percent differences for metabolic rates; expect 0 < diff3 < diff1 < diff2
barplot(percent_diff2, col = rgb(0, 0, 1, 0.5), ylab = "Metabolic rate percent diff")
barplot(percent_diff1, add = T)
barplot(percent_diff3, col = rgb(1, 0, 0, 0.5), add = T)
legend("topright", c("actual diff (1)", "no mass change diff (2)", "(3)"), cex = 0.8, 
       title = "Met rate comparison", fill = c("white", "blue", "red"))

# Barplot of percent differences for masses; expect 0 > diff4 > diff 5, uncertain diff6
barplot(percent_diff5, col = rgb(0, 0, 1, 0.5), ylab = "Mass percent diff")
barplot(percent_diff6, col = rgb(1, 0, 0, 0.5), add = T)
barplot(percent_diff4, add = T)
legend("topright", c("actual diff (4)", "needed diff (5)", "(6)"), cex = 0.8, title = "Mass comparison", 
       fill = c("white", "blue", "red"))

# Barplot of average percent difference for metabolic rates
# Possibly add SD in future
barplot(c(avg_pd_1, avg_pd_2, avg_pd_3), col = c("white", "blue", "red"), 
        names.arg = c("actual diff (1)", "no mass change diff (2)", "(3)"), 
        xlab = "Met rate comparison", ylab = "Avg percent diff")

# Barplot of average percent difference for masses
# Possibly add SD in future
barplot(c(avg_pd_4, avg_pd_5, avg_pd_6), col = c("white", "blue", "red"), 
        names.arg = c("actual diff (4)", "needed diff (5)", "(6)"), 
        xlab = "Mass comparison", ylab = "Avg percent diff")

#-----------------------------ABSOLUTE TEMP-------------------------------

# Arrange dataset according to class
#three_degree_pairs = three_degree_pairs[with(three_degree_pairs, order(Class)), ]

# Arrange dataset according to initial temperature
#three_degree_pairs = three_degree_pairs[with(three_degree_pairs, order(initial_temp)), ]

