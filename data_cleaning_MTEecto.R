library(dplyr)

#-------------------------CALCULATE CLASS VALUES----------------------------

### Use existing physiological values, literature estimates, or theoretical values
### to obtain scaling exponents and activation energies for all classes

# Classes: Branchiopoda, Malacostraca, Maxillopoda
# Data source: Makrieva
MR_data = read.csv("Class_metabolicrates_Makrievadata.csv", stringsAsFactors = FALSE)
Class_list=unique(MR_data$Class)
InvK=1/(273.15+MR_data$TC) #getting all the units correct
Mkg=MR_data$Mg/1000
W=MR_data$qWkg*Mkg
MR_data = cbind(MR_data, InvK, W)

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

#---------------------REARRANGE EXPERIMENTAL DATASET----------------------------

### Get dataset with mass-temperature data paired for 3* temperature difference
### Each experiment has unique identifier

# Read in original compiled experimental dataset
original_data = read.csv("MTEEcto_data.csv", stringsAsFactors = FALSE, fileEncoding = "latin1")

# Turn unique row names into column in original dataset
original_data = data.frame(as.numeric(rownames(original_data)), original_data)
colnames(original_data)[1] = "row_names"

# Put temperatures into ascending order
original_data = original_data %>%
  arrange(temp)

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
names(all_pairs) = c("studyID", "initial_temp_row", "final_temp_row", "initial_temp", "final_temp", "initial_mass_mg", "final_mass_mg")

# Subset dataset to contain only the lowest temperature with all its possible pairs
lowest_deg_pairs = all_pairs %>%
  group_by(studyID) %>%
  filter(initial_temp == min(initial_temp))

# Add species to three degree difference pairs
lookup = original_data[match(lowest_deg_pairs$studyID, original_data$studyID),]
lowest_deg_pairs$species = lookup$Species
lowest_deg_pairs$Class = lookup$Class
lowest_deg_pairs$ref = lookup$ref

# List of species and number of pairs per species
species_duplicates = sort(table(lowest_deg_pairs$species), decreasing=TRUE)

# Add class-specific values to pairs dataset
store_class_values = class_values[match(lowest_deg_pairs$Class, class_values$Class),]
lowest_deg_pairs = as.data.frame(lowest_deg_pairs)
lowest_deg_pairs = cbind(lowest_deg_pairs, store_class_values$exponent, store_class_values$Ea)
colnames(lowest_deg_pairs)[11] = "exponent"
colnames(lowest_deg_pairs)[12] = "Ea"

# List of classes and number of pairs per class
class_distribution = table(lowest_deg_pairs$Class)

# Change mass units from milligrams to nanograms
lowest_deg_pairs$initial_mass = lowest_deg_pairs$initial_mass_mg * 1000000
lowest_deg_pairs$final_mass = lowest_deg_pairs$final_mass_mg * 1000000

# Separate studyID into study and trial columns
library(stringr)
study = str_sub(lowest_deg_pairs$studyID, 1, -2)
trial = str_sub(lowest_deg_pairs$studyID, -1)
lowest_deg_pairs = cbind(study, trial, lowest_deg_pairs)

#------------------CALCULATE METABOLIC RATES/MASSES--------------------------

### Use MTE equation: Metabolic Rate = constant * (Mass ^exponent) * e ^ Activation Energy/kt

# Use MTE equation to calculate metabolic rates (initial, final, and w/ constant/initial mass)
lowest_deg_pairs$initial_metrate = (lowest_deg_pairs$initial_mass ^ lowest_deg_pairs$exponent) * (exp(((lowest_deg_pairs$Ea) / (.00008617 * (lowest_deg_pairs$initial_temp + 273.15)))))
lowest_deg_pairs$final_metrate = (lowest_deg_pairs$final_mass ^ lowest_deg_pairs$exponent) * (exp(((lowest_deg_pairs$Ea) / (.00008617 * (lowest_deg_pairs$final_temp + 273.15)))))
lowest_deg_pairs$constantmass_metrate = (lowest_deg_pairs$initial_mass ^ lowest_deg_pairs$exponent) * (exp(((lowest_deg_pairs$Ea) / (.00008617 * (lowest_deg_pairs$final_temp + 273.15)))))

# Use rearranged MTE equation to calculate mass w/ constant/initial metabolic rate
lowest_deg_pairs$constantmetrate_mass = (lowest_deg_pairs$initial_metrate / (exp(lowest_deg_pairs$Ea / (.00008617 * (lowest_deg_pairs$final_temp + 273.15))))) ^ (1 / lowest_deg_pairs$exponent)

# Create CSV of cleaned and prepared data
write.csv(lowest_deg_pairs, file = "clean_data_MTEecto_expanded.csv")
