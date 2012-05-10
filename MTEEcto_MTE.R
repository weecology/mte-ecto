#File for calculating Q4s for each species and for generating population-level comparisons

data = read.csv("MTEEcto_data.csv", stringsAsFactors = FALSE)
#stringsAsFactors prevents R from reading strings as factors, which mucks up further analysis

replicate=unique(data$replicate)
rep_linker=c()
#extracting species and replicate info
for (current_replicate in replicate){
  replicate_data = subset(data, data$replicate == current_replicate)
  replicate_species= unique(replicate_data$Species)
  replicate_class=unique(replicate_data$Class)
  rep_linker=rbind(rep_linker, c(current_replicate, replicate_species, replicate_class))
}
rep_linker_sort= rep_linker[order(rep_linker[,2], rep_linker[,1]),]
write.table(rep_linker_sort, file = "rep_linker.csv", sep = ",", col.names = TRUE)

#Creating RAW datafile
MTE_allreps=c()
MTE_allcompare=c()

for (current_replicate in replicate){
  MTE_repdata=c()
  current_data= subset (data, data$replicate == current_replicate)
  temps=unique(current_data$temp)
  
  for (current_temp in temps){
    min_mass=current_data$mass[current_data$temp==current_temp]
    repeater=length(temps)
  
    MTE_nochange=(min_mass^0.75)*(exp((-0.63/(.000086*(temps+273.15)))))
    MTE_change=(current_data$mass^0.75)*(exp((-0.63/(.000086*(temps+273.15)))))
    MTE_initial=(min_mass^0.75)*(exp((-0.63/(.000086*(current_temp+273.15)))))
        
    
    rep_vector=c(rep(current_replicate,repeater))
    Tmin_vector=c(rep(current_temp, repeater))
    MTE_repdata=cbind(temps,Tmin_vector, MTE_initial, MTE_nochange, MTE_change)
    MTE_repdata=as.data.frame(MTE_repdata)
    MTE_rep=data.frame(cbind(MTE_repdata, rep_vector), stringsAsFactors=FALSE)
    MTE_allreps=rbind(MTE_allreps, MTE_rep)
      
 }}
temp_diff=MTE_allreps$temps-MTE_allreps$Tmin_vector
MTE_allreps=cbind(MTE_allreps, temp_diff)
warming=warming= subset (MTE_allreps, MTE_allreps$temp_diff > 0)

Q_nochange=(warming$MTE_nochange/warming$MTE_initial)^(4/warming$temp_diff)
Q_change=(warming$MTE_change/warming$MTE_initial)^(4/warming$temp_diff) 
Q_diff=100*((Q_nochange/Q_change)-1)

Q_file=c()
Q_file=cbind(warming, Q_nochange, Q_change, Q_diff)
write.table(Q_file, file = "Q_file.csv", sep = ",", col.names = TRUE)
#END RAW DATAFILE CREATION

#Average Q for each study
rep_avgs=c()
for (unique_rep in replicate){
  rep_subset= subset (Q_file, Q_file$rep_vector == unique_rep) 
  #sp_name= rep_linker_sort[,2][ which (rep_linker_sort[,1] == unique_rep),]
  sp_name= subset(rep_linker_sort, rep_linker_sort[,1] == unique_rep,)
  avg_Q_nochange=mean(rep_subset$Q_nochange)
  avg_Q_change=mean(rep_subset$Q_change)
  rep_avgs=rbind(rep_avgs, c(unique_rep, avg_Q_nochange,avg_Q_change, sp_name))
  #rep_avgs=data.frame(rbind(rep_avgs, c(unique_rep, avg_Q_nochange,avg_Q_change, sp_name)), stringsAsFactors=FALSE)
}
 rep_avgs_1=  as.data.frame(rep_avgs, stringsAsFactors = FALSE)
 rep_avgs_2= transform(rep_avgs_1, V2 = as.numeric(V2), V3 = as.numeric(V3))
write.table(rep_avgs_2, file = "rep_avgs_2.csv", sep = ",", col.names = TRUE)

#END STUDY AVERAGES

#Average Q for species
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
Q_diff=100*((Q_sp$V3/Q_sp$V4)-1)
Q_sp=cbind(Q_sp, Q_diff)
Q_diff_alt=100*((Q_sp$V4/Q_sp$V3)-1)
Q_sp=cbind(Q_sp, Q_diff_alt)
colnames(Q_sp)=c("Species", "Class", "Q4", "Q4TSR", "Q4Q4TSR_diff", "Q4TSRQ4_diff")
write.table(Q_sp, file = "Q_sp.csv", sep = ",", col.names = TRUE)
# END AVGERAGE Q FOR SPECIES

#Pulling out population-level data
pop.study = read.csv("Population_data.csv", stringsAsFactors = FALSE)
pop.list=unique(pop.study$Study)
pop.data=subset(rep_avgs_2, rep_avgs_2$V4 %in% pop.list)
pop.data=as.data.frame(pop.data)
colnames(pop.data)=c("Study", "Q_nochange", "Q_TSD", "Study_copy", "Species", "Class")
pop.data=pop.data[order(pop.data$Study),]
pop.study=pop.study[order(pop.study$Study),]
pop.data=cbind(pop.data, pop.study$temperature)

#export pop.data table to add unique paper column
write.table(pop.data, file = "pop.data.csv", sep = ",", col.names = TRUE)

#import modified pop.data file
pop.data.1 = read.csv("pop.data.1.csv", stringsAsFactors = FALSE)
unique_papers=unique(pop.data.1$Papercode)
pop.results=c()
for (current_paper in unique_papers){
  paper.data=subset(pop.data.1, pop.data.1$Papercode == current_paper)
  QTSDregress=lm(paper.data$Q_TSD ~ paper.data$temperature)
  QTSDsummary=summary(QTSDregress)
  paper.results=cbind(current_paper, QTSDsummary$coefficients[2,1],QTSDsummary$coefficients[2,4])
  pop.results=rbind(pop.results, paper.results)
}
colnames(pop.results)=c("Papercode", "slope", "pvalue")

################ANALYSES#######################
average_Qdiff=mean(Q_sp$Q4TSRQ4_diff)
median_Qdiff=median(Q_sp$Q4TSRQ4_diff)
minmax=range(Q_sp$Q4TSRQ4_diff)

############PLOTS#####################
jpeg("C:/Documents and Settings/skme/My Dropbox/Morgane/MTEEcto/TSD_percentdiff.jpeg")
plot(density(Q_sp$Q4TSRQ4_diff),  col="red",main=" ",
xlab="Percent difference between TSD and non_TSD Q4s",)
title(main="How much would metabolic response to temperature be overestimated", line=3)
title(main="by ignoring TSD?", line=2)
dev.off()

jpeg("C:/Documents and Settings/skme/My Dropbox/Morgane/MTEEcto/TSD_MTE_density.jpeg")
#density functions of Qs
yes_range=range(Q_sp$Q4TSR)
plot(density(Q_sp$Q4), xlim = yes_range, main="", col="red",  
    xlab="x-fold increase in metabolic rate")
lines(density(Q_sp$Q4TSR))
legend("topleft", inset=.05, title="line colors",c("no TSD", "TSD"), fill=c("red", "black"))
dev.off()

Actinoperygii=subset(Q_sp, Q_sp$Class == "Actinoperygii")
Amphibia=subset(Q_sp, Q_sp$Class == "Amphibia")
Branchiopoda =subset(Q_sp, Q_sp$Class == "Branchiopoda")
Entognatha =subset(Q_sp, Q_sp$Class == "Entognatha")
Eurotatoria =subset(Q_sp, Q_sp$Class == "Eurotatoria")
Insecta =subset(Q_sp, Q_sp$Class == "Insecta")
Gastropoda =subset(Q_sp, Q_sp$Class == "Gastropoda")
Malacostraca=subset(Q_sp, Q_sp$Class == "Malacostraca")
Maxillopoda=subset(Q_sp, Q_sp$Class == "Maxillopoda")

jpeg("C:/Documents and Settings/skme/My Dropbox/Morgane/MTEEcto/Q_diff_classes.jpeg")
Q_diff_range=c(-40,40)
y_range=c(0,.1)
plot(density(Amphibia$Q4TSRQ4_diff), xlim=Q_diff_range, ylim=y_range, main="", col="goldenrod1",  
    xlab="Percent difference in Q4 values")
lines(density(Branchiopoda$Q4TSRQ4_diff), col="coral4")
lines(density(Malacostraca$Q4TSRQ4_diff),col="lightblue" )
lines(density(Maxillopoda$Q4TSRQ4_diff), col="grey19")
lines(density(Insecta$Q4TSRQ4_diff),col="magenta")
legend("topright", inset=.05, title="Taxonomic Class (# of spp)",c("Amphibia (6)", "Branchiopoda (4)","Insecta (69)", "Malacostraca (9)", "Maxillopoda (16)"), 
  fill=c("goldenrod1", "coral4", "magenta", "lightblue", "grey19" ))
dev.off()
        



  
  