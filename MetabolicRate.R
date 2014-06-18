
#--------------CALCULATING EXPONENTS AND ACTIVATION ENERGY FOR METABOLIC EQUATIONS-------------------------------

### For the metabolic rate equation Metabolic Rate = constant * (Mass ^exponent) * e ^ Activation Energy/kt
### This section uses metabolic rate, temperature, and size data to calculated fitted exponent and activation
### energy for the taxonomic classes for which we could find data. 

##-----------------------------ANALYSES WITH MAKRIEVA DATA FOR Branchiopoda Malacostraca and Maxillopoda---------------
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
  Ea=Tslope*0.000086
  class_slopes=rbind(class_slopes, c(current_class, Mslope, Tslope, Ea))
}

##-----------------------------ANALYSES WITH WHITE DATA FOR AMPHIBIANS--------------------------------------------
MR_Amphi = read.csv("Whiteetal_Amphibiandata.csv", stringsAsFactors = FALSE)
A_InvK=1/(273.15+MR_Amphi$TC)
A_Mkg=MR_Amphi$Mg/1000
MR_Amphi = cbind(MR_Amphi, A_InvK, A_Mkg)
A_basic=lm(log(MR_Amphi$Watts) ~ log(MR_Amphi$A_Mkg) + MR_Amphi$A_InvK)
Mslope=A_basic$coefficients[2]
Tslope=A_basic$coefficients[3]
Ea=Tslope*0.000086
class_slopes=rbind(class_slopes, c("Amphibia", Mslope, Tslope, Ea))
class_slopes=as.data.frame(class_slopes, stringsAsFactors = FALSE)

##-----------------------------ANALYSES WITH GILLOOLY DATA FOR FISH-----------------------------------------------

#obtaining super-class-specific exponents for metabolic rate scaling in fish
MR_fish = read.csv("gillooly_fish.csv", stringsAsFactors = FALSE)
f_Mkg=MR_fish$Mg/1000
MR_fish = cbind(MR_fish, f_Mkg)
F_basic=lm(log(MR_fish$W) ~ log(MR_fish$f_Mkg) + MR_fish$invK)
Mslope=F_basic$coefficients[2]
Tslope=F_basic$coefficients[3]
Ea=Tslope*0.000086
class_slopes=rbind(class_slopes, c("Actinoperygii", Mslope, Tslope, Ea))
class_slopes=as.data.frame(class_slopes, stringsAsFactors = FALSE)

##----------------CLEANING UP DATAFRAME WITH CLASS-SPECIFIC PARAMETERS FOR METABOLIC EQUATION------------------------

#removing unneccesary column and renaming columns
class_values=class_slopes[c(-3)]
names(class_values)[names(class_values) == "V1"] = "Class"
names(class_values)[names(class_values) == "log(Class_data$Mg)"] = "exponent"
names(class_values)[names(class_values) == "Class_data$InvK"] = "Ea"

##---------------ADDING Insecta VALUES FROM CHOWN AND THEORY VALUES FOR REMAINING CLASSES------------------------------
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

#----------------------Q10 CALCULATIONS---------------------------------------

TSD_data = read.csv("MTEEcto_data.csv", stringsAsFactors = FALSE, fileEncoding = "latin1")
#stringsAsFactors prevents R from reading strings as factors, which mucks up further analysis
#must re-encode data to properly read in data file

##-------------------EXTRACTING SPECIES AND REPLICATE INFO--------------------

### Within a scientific paper, the body size-temperature study may have replicated
### experiments or have slightly different experiments that were conducted. Each experiment
### is called a replicate and was given a unique identifier in "MTEEcto_data.csv". This 
### section of code extracts information so that each replicate can be analyzed separately

replicate=unique(TSD_data$studyID)
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
          MTE_initial=(min_mass^(class_MTE$exponent))*(exp(((class_MTE$Ea)/(.000086*(current_temp+273.15)))))
          MTE_nochange=min_mass^(class_MTE$exponent)*exp(class_MTE$Ea/(.000086*(temps+273.15)))
          MTE_change=(current_data$mass^(class_MTE$exponent))*(exp(((class_MTE$Ea)/(.000086*(temps+273.15)))))
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

##-----------------CALCULATING AVG REPLICATE Q4 FOR THE GENERATED TEMPERATURE-METABOLIC RATE DATA------------------

### Using the temperature-metabolic rate data generated above and stored in the dataframe 'warming'
### to calculate an average Q4 value for each replicate

Q_nochange=(warming$MTE_nochange/warming$MTE_initial)^(4/warming$temp_diff)
Q_change=(warming$MTE_change/warming$MTE_initial)^(4/warming$temp_diff) 
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


##---------------------------AVERAGE Q4 FOR EACH SPECIES-------------------------------------------

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
names(ClassQ_sp)[names(ClassQ_sp) == "V3"] = "Q4_noTSR"
names(ClassQ_sp)[names(ClassQ_sp) == "V4"] = "Q4_TSR"
ClassQ_sp= transform(ClassQ_sp, Q_diff = as.numeric(Q_diff))
write.table(ClassQ_sp, file = "ClassQ_sp.csv", sep = ",", col.names = TRUE)

#-------------------------SUMMARY STATS AND GRAPHS -------------------------
mean_Q4noTSR=mean(ClassQ_sp$Q4_noTSR)
mean_Q4TSR=mean(ClassQ_sp$Q4_TSR)
t.test(ClassQ_sp$Q4_noTSR,ClassQ_sp$Q4_TSR, paired=TRUE )
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
yes_range=range(ClassQ_sp$Q4_noTSR)
plot(density(ClassQ_sp$Q4_noTSR), xlim = yes_range, main="", col="dark grey", lwd=2, 
    xlab="Q4")
lines(density(ClassQ_sp$Q4_TSR), lwd=2)
legend("topleft", inset=.05, c("Q4-no size", "Q4-size"), fill=c("grey", "black"))

#graph of percent change in Q4 - all classes pooled
yes_range=range(ClassQ_sp$Q_diff)
plot(density(ClassQ_sp$Q_diff), xlim = yes_range, main="", col="black", lwd=2, 
    xlab="Percent difference between Q4-size response and Q4-no size")

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
    xlab="Percent difference in Q4 values")
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
legend("topleft", inset=.05, c("No Size Change (Q4=1.387)", "Size Change (Q4=1.296"), fill=c("black", "gray"))

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
legend("topright", inset=.05, c("No Size Change (Q4=1.387)", "Size Change (Q4=1.296"), fill=c("black", "gray"))

                      

