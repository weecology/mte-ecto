
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

#---------------------THREE DEGREE SUBSET DATASET----------------------------
# Read in compiled experimental dataset, has one temperature and replicate per row
original_data = read.csv("MTEEcto_data.csv", stringsAsFactors = FALSE, fileEncoding = "latin1")

# Turn unique row names into column in original dataset
original_data = data.frame(as.numeric(rownames(original_data)), original_data)
colnames(original_data)[1] = "row_names"

# Create list of experiment names and respective species
replicate_list=unique(original_data$studyID)

# Create dataset containing all temperature pairs within each replicate
all_pairs = c()
for(current_replicate in replicate_list){
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

##-------------------EXTRACTING SPECIES AND REPLICATE INFO--------------------

### Within a scientific paper, the body size-temperature study may have replicated
### experiments or have slightly different experiments that were conducted. Each experiment
### is called a replicate and was given a unique identifier in "MTEEcto_data.csv". This 
### section of code extracts information so that each replicate can be analyzed separately

rep_linker=c()

#extracting species and replicate info
for (current_replicate in replicate){
  replicate_data = subset(TSD_data, TSD_data$studyID == current_replicate)
  replicate_species= unique(replicate_data$Species)
  replicate_class=unique(replicate_data$Class)
  rep_linker=rbind(rep_linker, c(current_replicate, replicate_species, replicate_class))
}
rep_linker_sort= rep_linker[order(rep_linker[,2], rep_linker[,1]),]
#write.table(rep_linker_sort, file = "rep_linker.csv", sep = ",", col.names = TRUE)
class_uniquereplicate=unique(rep_linker_sort[,1])

##---------------------------CALCULATING METABOLIC RATES FOR TEMPERATURE/SIZE COMBOS FROM DATA----------------------

### This section pulls out the size-temperature data and metabolic rate parameters for a specific taxonomic class
### It extracts the data from a unique replicate. FOr that replicate it uses the size-temperature data to calculates 
### metabolic rates for the following cases: 1) no size change as temperature increases, 2) size and temperature increase
### as observed in the data. 

MTE_allreps=c()
#unused empty vector?
MTE_allcompare=c()
Classes_all=unique(TSD_data$Class)

for (index_class in Classes_all){
  class_TSD= subset(TSD_data, TSD_data$Class == index_class)
  class_MTE= subset(class_values1, class_values1$Class == index_class)
  class_replicate=unique(class_TSD$studyID)
    
    for (index_replicate in class_replicate){
      MTE_repdata=c()
      current_data= subset (class_TSD, class_TSD$studyID == index_replicate)
      temps=unique(current_data$temp)
 
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

temp_diff=MTE_allreps$temps-MTE_allreps$Tmin_vector
MTE_allreps=cbind(MTE_allreps, temp_diff)
warming= subset (MTE_allreps, MTE_allreps$temp_diff > 0) #cleans out decreasing temperature situations

##-----------------CALCULATING AVG REPLICATE Q3 FOR THE GENERATED TEMPERATURE-METABOLIC RATE DATA------------------

### Using the temperature-metabolic rate data generated above and stored in the dataframe 'warming'
### to calculate an average Q3 value for each replicate

#change 4*C change to 3*C change

Q_nochange=(warming$MTE_nochange/warming$MTE_initial)^(3/warming$temp_diff)
Q_change=(warming$MTE_change/warming$MTE_initial)^(3/warming$temp_diff) 
Q_diff=100*((Q_change/Q_nochange)-1)

Q_file=c()
Q_file=cbind(warming, Q_nochange, Q_change, Q_diff)
#write.table(Q_file, file = "ClassQ_file.csv", sep = ",", col.names = TRUE)
#END RAW DATAFILE CREATION

#Average Q for each study
rep_avgs=c()
for (unique_rep in class_uniquereplicate){
  rep_subset= subset (Q_file, Q_file$rep_vector == unique_rep) 
  #sp_name= rep_linker_sort[,2][ which (rep_linker_sort[,1] == unique_rep),]
  sp_name= subset(rep_linker_sort, rep_linker_sort[,1] == unique_rep,)
  avg_Q_nochange=mean(rep_subset$Q_nochange)
  avg_Q_change=mean(rep_subset$Q_change)
  rep_avgs=rbind(rep_avgs, c(unique_rep, avg_Q_nochange,avg_Q_change, sp_name))
}
 rep_avgs_1=  as.data.frame(rep_avgs, stringsAsFactors = FALSE)
 rep_avgs_2= transform(rep_avgs_1, V2 = as.numeric(V2), V3 = as.numeric(V3))
write.table(rep_avgs_2, file = "class_rep_avgs_2.csv", sep = ",", col.names = TRUE)


##---------------------------AVERAGE Q3 FOR EACH SPECIES-------------------------------------------

### If more than one replicate for a species, calculates an average across the replciates so that
### each species only has one resulting datapoint

Q_1=c()
species=unique(rep_avgs_2$V5)
for (current_sp in species){
  sp_subset= subset (rep_avgs_2, rep_avgs_2$V5 == current_sp) 
  sp_class=unique(sp_subset[,6])
  sp_Q_nochange=mean(sp_subset[,2])
  sp_Q_change=mean(sp_subset[,3])
  Q_1=rbind(Q_1, c(current_sp, sp_class, sp_Q_nochange, sp_Q_change))
}

Q_2=  as.data.frame(Q_1, stringsAsFactors = FALSE)
Q_sp= transform(Q_2, V3 = as.numeric(V3), V4 = as.numeric(V4))
Q_diff=100*((Q_sp$V4/Q_sp$V3)-1)
ClassQ_sp=cbind(Q_sp, Q_diff)
ClassQ_sp=ClassQ_sp[order(ClassQ_sp[,2], ClassQ_sp[,1]),]
ClassQ_sp=  as.data.frame(ClassQ_sp, stringsAsFactors = FALSE)

names(ClassQ_sp)[names(ClassQ_sp) == "V1"] = "Species"
names(ClassQ_sp)[names(ClassQ_sp) == "V2"] = "Class"
names(ClassQ_sp)[names(ClassQ_sp) == "V3"] = "Q3_noTSR"
names(ClassQ_sp)[names(ClassQ_sp) == "V4"] = "Q3_TSR"
ClassQ_sp= transform(ClassQ_sp, Q_diff = as.numeric(Q_diff))
write.table(ClassQ_sp, file = "ClassQ_sp.csv", sep = ",", col.names = TRUE)

#-------------------------SUMMARY STATS AND GRAPHS -------------------------
mean_Q3noTSR=mean(ClassQ_sp$Q3_noTSR)
mean_Q3TSR=mean(ClassQ_sp$Q3_TSR)
t.test(ClassQ_sp$Q3_noTSR,ClassQ_sp$Q3_TSR, paired=TRUE )
mean_Qdiff=mean(ClassQ_sp$Q_diff)

Qdiff.gt10=subset(ClassQ_sp, ClassQ_sp$Q_diff > 10)
Qdiff.lt10=subset(ClassQ_sp, ClassQ_sp$Q_diff < -10)
Qdiff.lt15=subset(ClassQ_sp, ClassQ_sp$Q_diff < -15)
Qdiff.lt20=subset(ClassQ_sp, ClassQ_sp$Q_diff < -20)
num.sp.gt.10.Qdiff=length(Qdiff.gt10$Q_diff)
num.sp.lt.10.Qdiff=length(Qdiff.lt10$Q_diff)
num.sp.lt.15.Qdiff=length(Qdiff.lt15$Q_diff)
num.sp.lt205.Qdiff=length(Qdiff.lt20$Q_diff)
num.sp=length(ClassQ_sp$Q_diff)

#graph of TSR and no TSR distributions - all classes pooled
yes_range=range(ClassQ_sp$Q3_noTSR)
plot(density(ClassQ_sp$Q3_noTSR), xlim = yes_range, main="", col="dark grey", lwd=2, 
    xlab="Q3")
lines(density(ClassQ_sp$Q3_TSR), lwd=2)
legend("topleft", inset=.05, c("Q3-no size", "Q3-size"), fill=c("grey", "black"))

#graph of percent change in Q3 - all classes pooled
yes_range=range(ClassQ_sp$Q_diff)
plot(density(ClassQ_sp$Q_diff), xlim = yes_range, main="", col="black", lwd=2, 
    xlab="Percent difference between Q3-size response and Q3-no size")

#Class-specific graph and results
Actinoperygii=subset(ClassQ_sp, ClassQ_sp$Class == "Actinoperygii")
Amphibia=subset(ClassQ_sp, ClassQ_sp$Class == "Amphibia")
Branchiopoda =subset(ClassQ_sp, ClassQ_sp$Class == "Branchiopoda")
Insecta =subset(ClassQ_sp, ClassQ_sp$Class == "Insecta")
Malacostraca=subset(ClassQ_sp, ClassQ_sp$Class == "Malacostraca")
Maxillopoda=subset(ClassQ_sp, ClassQ_sp$Class == "Maxillopoda")
Gastropoda=subset(ClassQ_sp, ClassQ_sp$Class == "Gastropoda")
Eurotatoria=subset(ClassQ_sp, ClassQ_sp$Class == "Eurotatoria")
Entognatha=subset(ClassQ_sp, ClassQ_sp$Class == "Entognatha")

Q_diff_range=c(-40,40)
y_range=c(0,.12)
plot(density(Amphibia$Q_diff), xlim=Q_diff_range, ylim=y_range, lwd=2, main="", col="goldenrod1",  
    xlab="Percent difference in Q3 values")
lines(density(Branchiopoda$Q_diff), lwd=2, col="coral4")
lines(density(Malacostraca$Q_diff),lwd=2, col="blue" )
lines(density(Maxillopoda$Q_diff), lwd=2, col="grey19")
lines(density(Insecta$Q_diff), lwd=2, col="magenta")
points(Gastropoda$Q_diff, .1, lwd=2, col="skyblue", pch=19)
points(Actinoperygii$Q_diff, .1, lwd=2, col="green", pch=19)
points(Eurotatoria$Q_diff, c(.1, .1), lwd=2, col="purple", pch=19)
points(Entognatha$Q_diff, .1, lwd=2, col="tomato2", pch=19)

legend("topright", inset=.05, title="Line (l) or Point colors (p)",c("Actinoperygii (p)", "Amphibia (l)", 
  "Branchiopoda (l)","Entognatha (p)", "Eurotatoria (p)","Gastropoda (p)","Insecta (l)", "Malacostraca (l)", "Maxillopoda (l)"), 
  fill=c("green", "goldenrod1", "coral4", "tomato2", "purple", "skyblue", "magenta", "blue", "grey19" ),)

Qdiff_Amphibia=mean(Amphibia$Q_diff)
Qdiff_Branchiopoda=mean(Branchiopoda$Q_diff)
Qdiff_Malacostraca=mean(Malacostraca$Q_diff)
Qdiff_Maxillopoda=mean(Maxillopoda$Q_diff)
Qdiff_Insecta=mean(Insecta$Q_diff)
Qdiff_Gastropoda=mean(Gastropoda$Q_diff)
Qdiff_Actinoperygii=mean(Actinoperygii$Q_diff)
Qdiff_Eurotatoria=mean(Eurotatoria$Q_diff)
Qdiff_Entognatha=mean(Entognatha$Q_diff)

#Exemplar graph for metabolic rate differences using Drosphila
Dmelano=subset(warming, warming$rep_vector == "TM61D")
DM_tmin=min(Dmelano$Tmin_vector)
DM_initialR1_1=Dmelano$MTE_initial[Dmelano$Tmin_vector == DM_tmin]
DM_initialR1=unique(DM_initialR1_1)

x=seq(0,4, 1)
plot(x, DM_initialR1*(1.387^(x/4)), type="l", lwd=2, xlab="Temperature change (degrees Celsius)",
     ylab="Metabolic Rate (watts)")
lines(x, DM_initialR1*(1.296^(x/4)), col="gray",lwd=2)
legend("topleft", inset=.05, c("No Size Change (Q3=1.387)", "Size Change (Q3=1.296"), fill=c("black", "gray"))

R2_TSR=DM_initialR1*(1.296^(x/4))
R2_noTSR=DM_initialR1*(1.387^(x/4))
Initial_N=100
Initial_E=Initial_N*DM_initialR1

N_TSR=Initial_E/R2_TSR
N_noTSR=Initial_E/R2_noTSR
N.data=cbind(x, R2_TSR, R2_noTSR, N_TSR, N_noTSR)
N.data= transform(N.data, x = as.numeric(x))

plot(N.data$x, N.data$N_noTSR, type="o", lwd=2, xlab="Temperature change (degrees Celsius)",
     ylab="Abundance")
lines(N.data$x, N.data$N_TSR, type="o", lwd=2, col="gray")
legend("topright", inset=.05, c("No Size Change (Q3=1.387)", "Size Change (Q3=1.296"), fill=c("black", "gray"))

                      
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


##-------------------Q10 CALCULATIONS ON 3* TEMP DIFF SUBSET--------------

# Repeating previous analysis to get Q10 values for replicate pairs, except only 
# for pairs that have a 3 degree difference in temperature; new sections indicated
# by four pound symbols

#### Get all columns from TSD_data into the subset dataset
all_subset_rownames = c(empirical_mass_change_realistic$row_name_1, empirical_mass_change_realistic$row_name_2)
all_subset_rownames = sort(unique(all_subset_rownames))

three_degree_subset_data = c()
for (current_row_name in all_subset_rownames){
  subset_data = subset(TSD_data, TSD_data$row_names == current_row_name)
  three_degree_subset_data = rbind(three_degree_subset_data, subset_data)
}

# Code from previous Q10 analysis, just replaced TSD_data with three_degree_subset_data
# What is the purpose of all this code?
replicate_subset=unique(three_degree_subset_data$studyID)
rep_linker_subset=c()

#extracting species and replicate info
for (current_replicate in replicate_subset){
  replicate_data_subset = subset(three_degree_subset_data, three_degree_subset_data$studyID == current_replicate)
  replicate_species_subset = unique(replicate_data_subset$Species)
  replicate_class_subset = unique(replicate_data_subset$Class)
  rep_linker_subset = rbind(rep_linker_subset, c(current_replicate, replicate_species_subset, replicate_class_subset))
}
rep_linker_sort_subset = rep_linker_subset[order(rep_linker_subset[,2], rep_linker_subset[,1]),]
#write.table(rep_linker_sort, file = "rep_linker.csv", sep = ",", col.names = TRUE)
class_uniquereplicate_subset=unique(rep_linker_sort_subset[,1])


MTE_allreps_subset=c()
#unused empty vector?
MTE_allcompare_subset=c()
Classes_all_subset=unique(three_degree_subset_data$Class)

for (index_class_subset in Classes_all_subset){
  class_TSD_subset= subset(three_degree_subset_data, three_degree_subset_data$Class == index_class_subset)
  class_MTE_subset= subset(class_values1, class_values1$Class == index_class_subset)
  class_replicate_subset=unique(class_TSD_subset$studyID)

  for (index_replicate_subset in class_replicate_subset){
    MTE_repdata_subset=c()
    current_data_subset= subset(class_TSD_subset, class_TSD_subset$studyID == index_replicate_subset)
    temps_subset=unique(current_data_subset$temp)

    for (current_temp_subset in temps_subset){
      min_mass_subset=current_data_subset$mass[current_data_subset$temp==current_temp_subset]
      repeater_subset=length(temps_subset)
      MTE_initial_subset=(min_mass_subset^(class_MTE_subset$exponent))*(exp(((class_MTE_subset$Ea)/(.00008617*(current_temp_subset+273.15)))))
      MTE_nochange_subset=min_mass_subset^(class_MTE_subset$exponent)*exp(class_MTE_subset$Ea/(.00008617*(temps_subset+273.15)))
      MTE_change_subset=(current_data_subset$mass^(class_MTE_subset$exponent))*(exp(((class_MTE_subset$Ea)/(.00008617*(temps_subset+273.15)))))
      rep_vector_subset=c(rep(index_replicate_subset,repeater_subset))
      Tmin_vector_subset=c(rep(current_temp_subset, repeater_subset))
      MTE_repdata_subset=cbind(temps_subset,Tmin_vector_subset, MTE_initial_subset, MTE_nochange_subset, MTE_change_subset)
      MTE_repdata_subset=as.data.frame(MTE_repdata_subset)
      MTE_rep_subset=data.frame(cbind(MTE_repdata_subset, rep_vector_subset), stringsAsFactors=FALSE)
      MTE_allreps_subset=rbind(MTE_allreps_subset, MTE_rep_subset)     
    }
  }
}

temp_diff_subset=MTE_allreps_subset$temps_subset-MTE_allreps_subset$Tmin_vector_subset
MTE_allreps_subset=cbind(MTE_allreps_subset, temp_diff_subset)
warming_subset= subset(MTE_allreps_subset, MTE_allreps_subset$temp_diff_subset > 0) #cleans out decreasing temperature situations

#### Want just 3 degree difference pairs
warming_subset = subset(warming_subset, warming_subset$temp_diff_subset == 3)


# Replicates' average Q10s
Q_nochange_subset=(warming_subset$MTE_nochange_subset/warming_subset$MTE_initial_subset)^(3/warming_subset$temp_diff_subset)
Q_change_subset=(warming_subset$MTE_change_subset/warming_subset$MTE_initial_subset)^(3/warming_subset$temp_diff_subset) 
Q_diff_subset=100*((Q_change_subset/Q_nochange_subset)-1)

Q_file_subset=c()
Q_file_subset=cbind(warming_subset, Q_nochange_subset, Q_change_subset, Q_diff_subset)
#write.table(Q_file, file = "ClassQ_file.csv", sep = ",", col.names = TRUE)
#END RAW DATAFILE CREATION

#Average Q for each study
rep_avgs_subset=c()
for (unique_rep_subset in class_uniquereplicate_subset){
  rep_subset_subset = subset(Q_file_subset, Q_file_subset$rep_vector_subset == unique_rep_subset) 
  #sp_name= rep_linker_sort[,2][ which (rep_linker_sort[,1] == unique_rep),]
  sp_name_subset= subset(rep_linker_sort_subset, rep_linker_sort_subset[,1] == unique_rep_subset,)
  avg_Q_nochange_subset=mean(rep_subset_subset$Q_nochange_subset)
  avg_Q_change_subset=mean(rep_subset_subset$Q_change_subset)
  rep_avgs_subset=rbind(rep_avgs_subset, c(unique_rep_subset, avg_Q_nochange_subset,avg_Q_change_subset, sp_name_subset))
}
rep_avgs_1_subset=  as.data.frame(rep_avgs_subset, stringsAsFactors = FALSE)
rep_avgs_2_subset= transform(rep_avgs_1_subset, V2 = as.numeric(V2), V3 = as.numeric(V3))
write.table(rep_avgs_2_subset, file = "class_rep_avgs_2_subset.csv", sep = ",", col.names = TRUE)


# Species' average Q10s
Q_1_subset=c()
species_subset=unique(rep_avgs_2_subset$V5)
for (current_sp_subset in species_subset){
  sp_subset_subset= subset(rep_avgs_2_subset, rep_avgs_2_subset$V5 == current_sp_subset) 
  sp_class_subset=unique(sp_subset_subset[,6])
  sp_Q_nochange_subset=mean(sp_subset_subset[,2])
  sp_Q_change_subset=mean(sp_subset_subset[,3])
  Q_1_subset=rbind(Q_1_subset, c(current_sp_subset, sp_class_subset, sp_Q_nochange_subset, sp_Q_change_subset))
}

Q_2_subset=  as.data.frame(Q_1_subset, stringsAsFactors = FALSE)
Q_sp_subset= transform(Q_2_subset, V3 = as.numeric(V3), V4 = as.numeric(V4))
Q_diff_subset=100*((Q_sp_subset$V4/Q_sp_subset$V3)-1)
ClassQ_sp_subset=cbind(Q_sp_subset, Q_diff_subset)
ClassQ_sp_subset=ClassQ_sp_subset[order(ClassQ_sp_subset[,2], ClassQ_sp_subset[,1]),]
ClassQ_sp_subset=as.data.frame(ClassQ_sp_subset, stringsAsFactors = FALSE)

names(ClassQ_sp_subset)[names(ClassQ_sp_subset) == "V1"] = "Species"
names(ClassQ_sp_subset)[names(ClassQ_sp_subset) == "V2"] = "Class"
names(ClassQ_sp_subset)[names(ClassQ_sp_subset) == "V3"] = "Q3_noTSR"
names(ClassQ_sp_subset)[names(ClassQ_sp_subset) == "V4"] = "Q3_TSR"
ClassQ_sp_subset= transform(ClassQ_sp_subset, Q_diff_subset = as.numeric(Q_diff_subset))
write.table(ClassQ_sp_subset, file = "ClassQ_sp_subset.csv", sep = ",", col.names = TRUE)


# Summary stats and graphs
mean_Q3noTSR_subset=mean(ClassQ_sp_subset$Q3_noTSR)
mean_Q3TSR_subset=mean(ClassQ_sp_subset$Q3_TSR)
t.test(ClassQ_sp_subset$Q3_noTSR,ClassQ_sp_subset$Q3_TSR, paired=TRUE )
mean_Qdiff_subset=mean(ClassQ_sp_subset$Q_diff_subset)

Qdiff.gt10_subset=subset(ClassQ_sp_subset, ClassQ_sp_subset$Q_diff_subset > 10)
Qdiff.lt10_subset=subset(ClassQ_sp_subset, ClassQ_sp_subset$Q_diff_subset < -10)
Qdiff.lt15_subset=subset(ClassQ_sp_subset, ClassQ_sp_subset$Q_diff_subset < -15)
Qdiff.lt20_subset=subset(ClassQ_sp_subset, ClassQ_sp_subset$Q_diff_subset < -20)
num.sp.gt.10.Qdiff_subset=length(Qdiff.gt10_subset$Q_diff_subset)
num.sp.lt.10.Qdiff_subset=length(Qdiff.lt10_subset$Q_diff_subset)
num.sp.lt.15.Qdiff_subset=length(Qdiff.lt15_subset$Q_diff_subset)
num.sp.lt205.Qdiff_subset=length(Qdiff.lt20_subset$Q_diff_subset)
num.sp_subset=length(ClassQ_sp_subset$Q_diff_subset)

# #graph of TSR and no TSR distributions - all classes pooled
# yes_range_subset=range(ClassQ_sp_subset$Q3_noTSR)
# plot(density(ClassQ_sp_subset$Q3_noTSR), xlim = yes_range_subset, main="", col="dark grey", lwd=2, 
#      xlab="Q3")
# lines(density(ClassQ_sp_subset$Q3_TSR), lwd=2)
# legend("topleft", inset=.05, c("Q3-no size", "Q3-size"), fill=c("grey", "black"))
# 
# #graph of percent change in Q3 - all classes pooled
# yes_range_subset=range(ClassQ_sp_subset$Q_diff_subset)
# plot(density(ClassQ_sp_subset$Q_diff_subset), xlim = yes_range_subset, main="", col="black", lwd=2, 
#      xlab="Percent difference between Q3-size response and Q3-no size")

#Class-specific graph and results
Actinoperygii_subset=subset(ClassQ_sp_subset, ClassQ_sp_subset$Class == "Actinoperygii")
Amphibia_subset=subset(ClassQ_sp_subset, ClassQ_sp_subset$Class == "Amphibia")
Branchiopoda_subset =subset(ClassQ_sp_subset, ClassQ_sp_subset$Class == "Branchiopoda")
Insecta_subset =subset(ClassQ_sp_subset, ClassQ_sp_subset$Class == "Insecta")
Malacostraca_subset=subset(ClassQ_sp_subset, ClassQ_sp_subset$Class == "Malacostraca")
Maxillopoda_subset=subset(ClassQ_sp_subset, ClassQ_sp_subset$Class == "Maxillopoda")
Gastropoda_subset=subset(ClassQ_sp_subset, ClassQ_sp_subset$Class == "Gastropoda")
Eurotatoria_subset=subset(ClassQ_sp_subset, ClassQ_sp_subset$Class == "Eurotatoria")
Entognatha_subset=subset(ClassQ_sp_subset, ClassQ_sp_subset$Class == "Entognatha")

# Q_diff_range_subset=c(-40,40)
# y_range_subset=c(0,.18)
# plot(density(Amphibia_subset$Q_diff_subset), xlim=Q_diff_range_subset, ylim=y_range_subset, lwd=2, main="", col="goldenrod1",  
#      xlab="Percent difference in Q3 values")
# lines(density(Branchiopoda_subset$Q_diff_subset), lwd=2, col="coral4")
# lines(density(Malacostraca_subset$Q_diff_subset),lwd=2, col="blue" )
# lines(density(Maxillopoda_subset$Q_diff_subset), lwd=2, col="grey19")
# lines(density(Insecta_subset$Q_diff_subset), lwd=2, col="magenta")
# points(Gastropoda_subset$Q_diff_subset, .1, lwd=2, col="skyblue", pch=19)
# points(Actinoperygii_subset$Q_diff_subset, .1, lwd=2, col="green", pch=19)
# points(Eurotatoria_subset$Q_diff_subset, c(.1, .1), lwd=2, col="purple", pch=19)
# points(Entognatha_subset$Q_diff_subset, .1, lwd=2, col="tomato2", pch=19)
# 
# legend("topright", inset=.05, title="Line (l) or Point colors (p)",c("Actinoperygii (p)", "Amphibia (l)", 
#                                                                      "Branchiopoda (l)","Entognatha (p)", "Eurotatoria (p)","Gastropoda (p)","Insecta (l)", "Malacostraca (l)", "Maxillopoda (l)"), 
#        fill=c("green", "goldenrod1", "coral4", "tomato2", "purple", "skyblue", "magenta", "blue", "grey19" ),)
# 
Qdiff_Amphibia_subset=mean(Amphibia_subset$Q_diff_subset)
Qdiff_Branchiopoda_subset=mean(Branchiopoda_subset$Q_diff_subset)
Qdiff_Malacostraca_subset=mean(Malacostraca_subset$Q_diff_subset)
Qdiff_Maxillopoda_subset=mean(Maxillopoda_subset$Q_diff_subset)
Qdiff_Insecta_subset=mean(Insecta_subset$Q_diff_subset)
Qdiff_Gastropoda_subset=mean(Gastropoda_subset$Q_diff_subset)
Qdiff_Actinoperygii_subset=mean(Actinoperygii_subset$Q_diff_subset)
Qdiff_Eurotatoria_subset=mean(Eurotatoria_subset$Q_diff_subset)
Qdiff_Entognatha_subset=mean(Entognatha_subset$Q_diff_subset)


# Compare entire data and 3* subset results figures
compare_plot_range = c(0.9, 1.7)
plot(density(ClassQ_sp_subset$Q3_noTSR), xlim = compare_plot_range, main="", col="dark grey", lwd=2, 
     xlab="Q3")
lines(density(ClassQ_sp_subset$Q3_TSR), lwd=2)
lines(density(ClassQ_sp$Q3_noTSR, lwd=2))
lines(density(ClassQ_sp$Q3_TSR), lwd=2)
legend("topleft", inset=.05, c("Q3-no size", "Q3-size"), fill=c("grey", "black"))


# #### Direct (non-Q10) comparison of initial and body size change metabolic rates
# warming_subset$MTE_change_diff = 100 * (warming_subset$MTE_change_subset / warming_subset$MTE_initial_subset)
# plot(density(warming_subset$MTE_change_diff), lwd=2)


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


