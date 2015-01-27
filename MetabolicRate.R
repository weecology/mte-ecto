
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

other_values=as.data.frame(other_values)
names(other_values)[names(other_values) == "I_class"] = "Class"
names(other_values)[names(other_values) == "I_slope"] = "exponent"
names(other_values)[names(other_values) == "I_Ea"] = "Ea"
class_values=rbind(class_values, other_values)
class_values1= transform(class_values, exponent = as.numeric(exponent), Ea = as.numeric(Ea))

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
names(all_pairs) = c("studyID", "lower_temp_row", "higher_temp_row", "lower_temp", "higher_temp", "lower_temp_mass", "higher_temp_mass")

# Subset dataset to contain only temperature pairs that differ by 3*
all_pairs$temp_diff = all_pairs$higher_temp - all_pairs$lower_temp
three_degree_pairs = subset(all_pairs, all_pairs$temp_diff == 3)

# Add species to three degree difference pairs
lookup = original_data[match(three_degree_pairs$studyID, original_data$studyID),]
three_degree_pairs$species = lookup$Species

# List of species and number of pairs per species
species_duplicates = sort(table(three_degree_pairs$species), decreasing=TRUE)

# Use row names for 3* difference pairs to create dataset where each row has 
# one temperature and replicate
three_degree_rows = c(three_degree_pairs$lower_temp_row, three_degree_pairs$higher_temp_row)
three_degree_rows = sort(unique(three_degree_rows))

three_degree_data = c()
for (current_row in three_degree_rows){
  row_data = subset(original_data, original_data$row_names == current_row)
  three_degree_data = rbind(three_degree_data, row_data)
}

# Create list of experiment names and respective species
replicate_list = unique(three_degree_data$studyID)

#List of replicates with respective species and classes
replicate_info = c()
for (current_replicate in replicate_list){
  replicate_data = subset(three_degree_data, three_degree_data$studyID == current_replicate)
  replicate_species = unique(replicate_data$Species)
  replicate_class = unique(replicate_data$Class)
  replicate_info = rbind(replicate_info, c(current_replicate, replicate_species, replicate_class))
}

##-------------------------CALCULATING METABOLIC RATES------------------------

### This section pulls out the size-temperature data and metabolic rate parameters for a specific taxonomic class
### It extracts the data from a unique replicate. FOr that replicate it uses the size-temperature data to calculates 
### metabolic rates for the following cases: 1) no size change as temperature increases, 2) size and temperature increase
### as observed in the data. 

MTE_allreps = c()
#unused empty vector?
MTE_allcompare =c()
Classes_all = unique(three_degree_data$Class)

for (index_class in Classes_all){
  class_TSD = subset(three_degree_data, three_degree_data$Class == index_class)
  class_MTE = subset(class_values1, class_values1$Class == index_class)
  class_replicate = unique(class_TSD$studyID)
    
    for (index_replicate in class_replicate){
      MTE_repdata = c()
      current_data = subset (class_TSD, class_TSD$studyID == index_replicate)
      temps = unique(current_data$temp)
 
        for (current_temp in temps){
          min_mass=current_data$mass[current_data$temp==current_temp]
          repeater=length(temps)
          MTE_initial=(min_mass^(class_MTE$exponent))*(exp(((class_MTE$Ea)/(.00008617*(current_temp+273.15)))))
          MTE_nochange=min_mass^(class_MTE$exponent)*exp(class_MTE$Ea/(.00008617*(temps+273.15)))
          MTE_change=(current_data$mass^(class_MTE$exponent))*(exp(((class_MTE$Ea)/(.00008617*(temps+273.15)))))
          rep_vector=c(rep(index_replicate,repeater))
          Tmin_vector=c(rep(current_temp, repeater))
          MTE_repdata=cbind(temps,Tmin_vector, MTE_initial, MTE_nochange, MTE_change)
          MTE_repdata=as.data.frame(MTE_repdata)
          MTE_rep=data.frame(cbind(MTE_repdata, rep_vector), stringsAsFactors=FALSE)
          MTE_allreps=rbind(MTE_allreps, MTE_rep)     
          }
     }
}

# Limit to 3* difference pairs
MTE_allreps$temp_diff = MTE_allreps$temps-MTE_allreps$Tmin_vector
metabolic_rates = subset(MTE_allreps, MTE_allreps$temp_diff == 3) 

##------------------------Analysis 1: Q10 CALCULATIONS--------------------------

### Use temperatures and metabolic rates to determine Q10 values for body size change
### and no body size change for each temperature/replicate combination, then replicate
### average and species average

# Individual Q10s
Q_nochange=(metabolic_rates$MTE_nochange/metabolic_rates$MTE_initial)^(3/metabolic_rates$temp_diff)
Q_change=(metabolic_rates$MTE_change/metabolic_rates$MTE_initial)^(3/metabolic_rates$temp_diff) 
Q_diff=100*((Q_change/Q_nochange)-1)

Q_individual=c()
Q_individual=cbind(metabolic_rates, Q_nochange, Q_change, Q_diff)

# Replicate average Q10
Q_replicate=c()
for (unique_rep in replicate_list){
  rep_subset = subset(Q_individual, Q_individual$rep_vector == unique_rep) 
  sp_name = subset(replicate_info, replicate_info[,1] == unique_rep,)
  avg_Q_nochange = mean(rep_subset$Q_nochange)
  avg_Q_change = mean(rep_subset$Q_change)
  Q_replicate = rbind(Q_replicate, c(avg_Q_nochange,avg_Q_change, sp_name))
}

Q_replicate = as.data.frame(Q_replicate, stringsAsFactors = FALSE)
colnames(Q_replicate) = c("Q_nochange", "Q_change", "studyID", "species", "class")
Q_replicate = transform(Q_replicate, Q_nochange = as.numeric(Q_nochange), Q_change = as.numeric(Q_change))

# Species average Q10
Q_species=c()
species=unique(Q_replicate$species)
for (current_sp in species){
  sp_subset = subset(Q_replicate, Q_replicate$species == current_sp) 
  sp_class = unique(sp_subset$class)
  sp_Q_nochange = mean(sp_subset$Q_nochange)
  sp_Q_change = mean(sp_subset$Q_change)
  Q_species = rbind(Q_species, c(sp_Q_nochange, sp_Q_change, current_sp, sp_class))
}

Q_species = as.data.frame(Q_species, stringsAsFactors = FALSE)
colnames(Q_species) = c("Q_nochange", "Q_change", "species", "class")
Q_species = transform(Q_species, Q_nochange = as.numeric(Q_nochange), Q_change = as.numeric(Q_change))

Q_species$Q_diff = ((Q_species$Q_change / Q_species$Q_nochange) - 1) * 100
Q_species = transform(Q_species, Q_diff = as.numeric(Q_diff))
Q_species = Q_species[order(Q_species$Q_change, Q_species$Q_nochange),]

### Summary stats and figures

# Average Q10 values
mean_Q_nochange = mean(Q_species$Q_nochange)
mean_Q_change = mean(Q_species$Q_change)
mean_Q_diff = mean(Q_species$Q_diff)

# Q10 t-test across classes
ttest_Q = t.test(Q_species$Q_nochange, Q_species$Q_change, paired = TRUE)

# Q10 change and no change distributions across classes
Q_comparison = plot(density(Q_species$Q_nochange), main = "", xlab = "Q10 values",
                    ylab = "Q10 density")
lines(density(Q_species$Q_change), lty = 2)
legend("topleft", inset=.05, c("Q3 no change", "Q3 change"), lty = c(1, 2))

# Q10 change and no change percent difference figure
Q_percent_diff = plot(density(Q_species$Q_diff), main = "", xlab = "Q10 percent differences")

# Subset by class
Actinoperygii = subset(Q_species, Q_species$class == "Actinoperygii")
Amphibia = subset(Q_species, Q_species$class == "Amphibia")
Branchiopoda = subset(Q_species, Q_species$class == "Branchiopoda")
Insecta = subset(Q_species, Q_species$class == "Insecta")
Malacostraca = subset(Q_species, Q_species$class == "Malacostraca")
Maxillopoda = subset(Q_species, Q_species$class == "Maxillopoda")
Entognatha = subset(Q_species, Q_species$class == "Entognatha")

Qdiff_Amphibia=mean(Amphibia$Q_diff)
Qdiff_Branchiopoda=mean(Branchiopoda$Q_diff)
Qdiff_Malacostraca=mean(Malacostraca$Q_diff)
Qdiff_Maxillopoda=mean(Maxillopoda$Q_diff)
Qdiff_Insecta=mean(Insecta$Q_diff)
Qdiff_Actinoperygii=mean(Actinoperygii$Q_diff)
Qdiff_Entognatha=mean(Entognatha$Q_diff)

# Class Q10 percent difference distributions
Q_diff_range=c(-40, 100)
y_range=c(0,.15)
plot(density(Insecta$Q_diff), xlim=Q_diff_range, ylim=y_range, lwd=2, main="", col="goldenrod1",  
    xlab="Percent difference in Q3 values")
points(Branchiopoda$Q_diff, .1, lwd=2, col="coral4", pch=19)
lines(density(Malacostraca$Q_diff),lwd=2, col="blue" )
lines(density(Maxillopoda$Q_diff), lwd=2, col="grey19")
points(Amphibia$Q_diff, .1, lwd=2, col="magenta", pch=19)
points(Actinoperygii$Q_diff, .1, lwd=2, col="green", pch=19)
points(Entognatha$Q_diff, .1, lwd=2, col="tomato2", pch=19)

legend("topright", inset=.05, title="Line (l) or Point (p) colors", c("Insecta (l)", "Branchiopoda (p)",  "Malacostraca (l)", "Maxillopoda (l)", "Amphibia (p)", "Actinoperygii (p)", "Entognatha (p)"),
  fill=c("goldenrod1", "coral4", "blue", "grey19", "magenta", "green", "tomato2" ),)
       
##-----------------CALCULATE COMPENSATION BODY MASS----------------------------------

### Using metabolic rates for initial mass and temperature, which have already ("MTE_initial") 
### been calculated, to determine what mass each species would need to be to compensate for a
### 3*C increase in temperature. 


#create dataframe containing initial temp & mass rows for each replicate with calculated
# metabolic rate from warming data frame
compensation_mass_calcs = unique(warming[,c("rep_vector", "Tmin_vector", "MTE_initial")])
colnames(compensation_mass_calcs) = c("rep_vector", "initial_temp", "initial_MetRate")

#add original data to new dataframe
compensation_mass_data = merge(compensation_mass_calcs, TSD_data, by.x=c("rep_vector", "initial_temp"), 
                               by.y=c("studyID", "temp"))
names(compensation_mass_data)[names(compensation_mass_data) == "mass"] = "initial_mass"
compensation_mass_data = compensation_mass_data[order(compensation_mass_data$Species),]

#add increased temperature to dataframe
compensation_mass_data$comp_temp = compensation_mass_data$initial_temp + 3

#add class-specific exponent and activation energy (from class_values1) to dataframe
getting_class_values = class_values1[match(compensation_mass_data$Class, class_values1$Class),]
compensation_mass_data = cbind(compensation_mass_data, getting_class_values$exponent, getting_class_values$Ea)
names(compensation_mass_data)[names(compensation_mass_data) == "getting_class_values$exponent"] = "class_exponent"
names(compensation_mass_data)[names(compensation_mass_data) == "getting_class_values$Ea"] = "class_Ea"

#calculate compensation mass using rearranged MTE equation
compensation_mass_data$comp_mass = (compensation_mass_data$initial_MetRate / 
                                  (exp(compensation_mass_data$class_Ea / (.00008617 * 
                                  (compensation_mass_data$comp_temp + 273.15))))) ^ (1 / 
                                  compensation_mass_data$class_exponent)

#determine difference between initial mass and compensation mass (should be positive)
compensation_mass_data$mass_diff = compensation_mass_data$initial_mass - compensation_mass_data$comp_mass

#determine what percent of initial mass the compensation mass is
#see CompMassTest.xls for test of first replicate
#why are these all basically the same number? in 60%ish range; maybe driven by class exponent & Ea?
compensation_mass_data$mass_reduction = (compensation_mass_data$comp_mass / 
                                           compensation_mass_data$initial_mass) * 100
hist(compensation_mass_data$mass_reduction, breaks=2, xlim=c(0,100))

#attempt to compare mass reduction percents for each class
insect_class = compensation_mass_data[compensation_mass_data$Class == "Insecta",]

#histograms for each class' mass reduction
pdf("ClassMassReduction.pdf")
par(mfrow=c(3,3))

for (current_class in class_values1$Class){
  class_subset = subset(compensation_mass_data, compensation_mass_data$Class == current_class)
  hist(class_subset$mass_reduction, xlim=c(50,90), main=NULL)
  mtext(paste("Class:", class_subset$Class))
}

dev.off()

# #figures to compare initial and compensation masses
# #plot(rownames(compensation_mass_data), compensation_mass_data$initial_mass)
# initial_mass_hist = hist(compensation_mass_data$initial_mass)
# comp_mass_hist = hist(compensation_mass_data$comp_mass)
# plot(initial_mass_hist, breaks=seq(0,1000,by=20), col=rgb(0,0,1,1/4))
# plot(comp_mass_hist, breaks=seq(0,1000,by=20), col=rgb(1,0,0,1/4), add=T)

#averaging of compensation masses and comparison metrics

##-------------EMPIRICAL MASS CHANGES FOR EXTREME TEMP CHANGES------------------

# Calculate mass reduction for most extreme temperature changes in each of the 
# replicates to determine what magnitude of mass change is biologically feasible


# Get lowest and highest temperatures (and then corresponding masses) for each replicate
empirical_mass_change_extreme = c()
for(current_replicate in replicate){
  single_replicate = subset(TSD_data, TSD_data$studyID == current_replicate)
  highest_temperature = subset(single_replicate, single_replicate$temp == max(single_replicate$temp))
  lowest_temperature = subset(single_replicate, single_replicate$temp == min(single_replicate$temp))
  empirical_mass_change_extreme = rbind(empirical_mass_change_extreme, c(current_replicate, highest_temperature$temp, highest_temperature$mass, lowest_temperature$temp, lowest_temperature$mass))
}

# Necessary changes to make dataframe useful
empirical_mass_change_extreme = data.frame(empirical_mass_change_extreme)
names(empirical_mass_change_extreme) = c("studyID", "highest_temp", "highest_temp_mass", "lowest_temp", "lowest_temp_mass")
empirical_mass_change_extreme[,c("highest_temp", "highest_temp_mass", "lowest_temp", "lowest_temp_mass")] = as.numeric(as.character(unlist(empirical_mass_change_extreme[,c("highest_temp", "highest_temp_mass", "lowest_temp", "lowest_temp_mass")])))

# Get temperature difference between high and low temps for each replicate
empirical_mass_change_extreme$temp_difference = empirical_mass_change_extreme$highest_temp - empirical_mass_change_extreme$lowest_temp

# Get mass reduction value for each replicate, with any increase in body size resulting in null
empirical_mass_change_extreme$mass_reduction = (empirical_mass_change_extreme$highest_temp_mass / empirical_mass_change_extreme$lowest_temp_mass) * 100
empirical_mass_change_extreme$mass_reduction[empirical_mass_change_extreme$mass_reduction >= 100] = NA
hist(empirical_mass_change_extreme$mass_reduction)

# Create histogram for each of 3 temp difference categories:
#   Realistic range: 0-4 degrees
#   Low range: 0-10 degrees
#   Medium range: 10-20 degrees
#   High range: 20-30 degrees
realistic_temp_diff = subset(empirical_mass_change_extreme, empirical_mass_change_extreme$temp_difference <= 4)
low_temp_diff = subset(empirical_mass_change_extreme, empirical_mass_change_extreme$temp_difference < 10)
med_temp_diff = subset(empirical_mass_change_extreme, (10 <= empirical_mass_change_extreme$temp_difference) & (empirical_mass_change_extreme$temp_difference < 20))
high_temp_diff = subset(empirical_mass_change_extreme, empirical_mass_change_extreme$temp_difference >= 20)

hist(realistic_temp_diff$mass_reduction, xlim = c(0,100))
hist(low_temp_diff$mass_reduction, xlim = c(0,100))
hist(med_temp_diff$mass_reduction, xlim = c(0,100))
hist(high_temp_diff$mass_reduction, xlim = c(0,100))


##-------------EMPIRICAL MASS CHANGES FOR REALISTIC TEMP CHANGES--------------

# Calculate mass reduction for realistic temperature change (i.e., 3*), replacing values over 100 w/ NA
empirical_mass_change_realistic$mass_reduction = (empirical_mass_change_realistic$higher_temp_mass / empirical_mass_change_realistic$lower_temp_mass) * 100
empirical_mass_change_realistic$mass_reduction[empirical_mass_change_realistic$mass_reduction >= 100] = NA
realistic_histogram = hist(empirical_mass_change_realistic$mass_reduction, breaks=3, xlim = c(0,100))
realistic_histogram$breaks
realistic_histogram$counts
(realistic_histogram$counts / sum(realistic_histogram$counts)) * 100

# Compare empirical realistic temperature change mass reduction to theoretical
theoretical_histogram = hist(compensation_mass_data$mass_reduction, breaks=2, xlim=c(0,100))
theoretical_histogram$breaks
theoretical_histogram$counts
(theoretical_histogram$counts / sum(theoretical_histogram$counts)) * 100

# Plot together
realistic_histogram = hist(empirical_mass_change_realistic$mass_reduction, breaks=3, xlim = c(0,100), ylim = c(0,200), col = "red")
theoretical_histogram = hist(compensation_mass_data$mass_reduction, breaks=2, xlim=c(0,100), add = T, col = "blue")

# Comparing empirical mass change to theoretical mass changed needed to maintain metabolic rate (i.e., compensation mass)

# Took this function from another project
get_lookup_matches = function(lookup, data_col1, data_col2, lookup_col1, lookup_col2){
  # Get info from lookup for two-column matches between lookup and data
  #
  # Args: 
  #   lookup: Dataframe that contains desired information
  #   data_col1: First column from dataset to match
  #   data_col2: Second column from dataset to match
  #   lookup_col1: First column from lookup to match to data_col1
  #   lookup_col2: Second column from lookup to match to data_col2
  #
  # Returns: 
  #   Lookup information that matches dataset columns
  dataset_info = interaction(data_col1, data_col2)
  lookup_info = interaction(lookup_col1, lookup_col2)
  lookup_matches = lookup[match(dataset_info, lookup_info),]
}

getting_comp = get_lookup_matches(compensation_mass_data, empirical_mass_change_realistic$studyID, empirical_mass_change_realistic$lower_temp, compensation_mass_data$rep_vector, compensation_mass_data$initial_temp)
empirical_mass_change_realistic$mass_reduction_theory = getting_comp$mass_reduction
empirical_mass_change_realistic$mass_red_diff = empirical_mass_change_realistic$mass_reduction - empirical_mass_change_realistic$mass_reduction_theory
mass_reduction_comparison = hist(empirical_mass_change_realistic$mass_red_diff)
plot(empirical_mass_change_realistic$higher_temp, empirical_mass_change_realistic$mass_red_diff)
abline(lm(empirical_mass_change_realistic$mass_red_diff ~ empirical_mass_change_realistic$higher_temp))
abline(h=0, col="red")

##-----------------------ABSOLUTE TEMPERATURE EFFECTS--------------------

# Plot absolute temp and Q10 value for each individual in entire dataset and 3* subset

# Did for entire dataset, but doesn't really work because many temp differences
# plot(Q_file$temps, Q_file$Q_nochange)
# abline(lm(Q_file$Q_nochange ~ Q_file$temps))
# plot(Q_file$temps, Q_file$Q_change)
# abline(lm(Q_file$Q_change ~ Q_file$temps))
# 
# min(Q_file$Q_change)
# max(Q_file$Q_change)

plot(Q_file_subset$temps_subset, Q_file_subset$Q_nochange_subset, ylim=c(1.0, 1.4))
abline(lm(Q_file_subset$Q_nochange_subset ~ Q_file_subset$temps_subset))
abline(h=1, col="red")
plot(Q_file_subset$temps_subset, Q_file_subset$Q_change_subset)
abline(lm(Q_file_subset$Q_change_subset ~ Q_file_subset$temps_subset))
abline(h=1, col="red")

min(Q_file_subset$Q_change_subset)
max(Q_file_subset$Q_change_subset)


