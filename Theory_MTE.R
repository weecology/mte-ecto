setwd("C:/Users/Morgan/Dropbox/Morgane/MTEEcto")

M2_values = c()
start.mass=.001
end.mass=100
interval=.1
M1_masses = seq(start.mass, end.mass, interval)
T1_temps = c(280.95, 282.95, 284.95, 286.95, 288.95)
b=c(.75)
EA=c(0.62)
models=c("Insecta")

#determining required mass change given EA and b for a specified change in T
for (current_mass in M1_masses){
  for (current_T1 in T1_temps){
    T2=seq(current_T1, current_T1+4, 1)
    for(current_T2 in T2){
    logM2=log(current_mass)+(EA/(.000086*b))*((1/current_T2) - (1/current_T1))
  M2_values=rbind(M2_values, c(current_mass, current_T1, current_T2, logM2 ))}}}
M2_values=as.matrix(M2_values)
M2=exp(M2_values[,4])
M2M1=M2/M2_values[,1]
T_diff=M2_values[,3]-M2_values[,2]
M2_values=cbind(M2_values, M2, M2M1, T_diff)
colnames(M2_values) = c("M1", "T1", "T2", "LogM2", "M2", "M2M1", "T_diff")
M2_values=as.data.frame(M2_values)

plot_data=subset(M2_values, M2_values$T_diff > 0)
plot_data=plot_data[order(plot_data$T_diff),]
T1_280=subset(plot_data, plot_data$T1 == 280.95)
T1_282=subset(plot_data, plot_data$T1 == 282.95)
T1_284=subset(plot_data, plot_data$T1 == 284.95)
T1_286=subset(plot_data, plot_data$T1 == 286.95)
T1_288=subset(plot_data, plot_data$T1 == 288.95)

subtitle_b=paste("b= ", b)
subtitle_EA=paste("Ea= ", EA)
subtitle_b_EA=paste(subtitle_b, subtitle_EA, sep=", ")
subtitle=paste(models, subtitle_b_EA, sep=", ")

plot(T1_280$M2M1 ~ T1_280$T_diff, type="b", col="blue", xlab="T2-T1 (deg K)", ylab=
  "M2/M1", )
title(main="Required size shift to maintain metabolic rate as temperature increases",)
text(1.7,.65, subtitle)
points(T1_282$M2M1 ~ T1_282$T_diff, type="b", col="darkturquoise")
points(T1_284$M2M1 ~ T1_284$T_diff, type="b", col="orange")
points(T1_286$M2M1 ~ T1_286$T_diff, type="b", col="red")
points(T1_288$M2M1 ~ T1_288$T_diff, type="b", col="darkmagenta")
legend(3.4, 0.9, c("280.95K","282.95K","284.95K", "286.95K","288.95K"), pch=1, col=c("blue", "darkturquoise","orange","red",
        "darkmagenta"), title="Initial Temperature")

#regression to determine if T1, M1, and T_diff are important for M2M1
multiple_regression= lm(M2_values$M2M1 ~ M2_values$M1 + M2_values$T1 + M2_values$T_diff)
summary(multiple_regression)


#examining Eab ratio
TSD_data = read.csv("MTEEcto_data.csv", stringsAsFactors = FALSE)

replicate=unique(TSD_data$replicate)
rep_freq=c()
for (unique_replicate in replicate){
  rep_number=length(TSD_data$replicate[TSD_data$replicate==unique_replicate])
  unique_rep_number=cbind(unique_replicate, rep_number)
  rep_freq=rbind(rep_freq, c(unique_rep_number))
}
  
frequent_replicates=subset(rep_freq, rep_freq[,2] > 3)
#frequent_replicates=c("KMD99A")
lnM=log(TSD_data$mass)
TSD_data=cbind(TSD_data, lnM)

LM_results=c()
for (unique_replicate in frequent_replicates){
  LM_subset=subset(TSD_data, TSD_data$replicate==unique_replicate)
  LM_fit=lm(LM_subset$lnM ~ LM_subset$temp)
  LM_coeff=coef(LM_fit)
  LM_results=rbind(LM_results, c(unique_replicate, LM_coeff, summary(LM_fit)$r.squared))
}
  
