library(dplyr)
library(ggplot2)
library(lme4)
library(grid)
library(cowplot)

pairs_data = read.csv("clean_data_MTEecto_expanded.csv", stringsAsFactors = FALSE)

### Dataset references
forster_list = c("Abdullahi & Parry (1985) ", " Abdullahi & Parry (1985  )", 
                 "Leandro Tiselius P & Queiroga H(2006) ", "Hansen et al.(2009) ", 
                 "Booth & Kiddell(2007) ", "Lamb & MacKay(1987)", " MacKay & Lamb(1988)  ", 
                 " Jalali Tirry & De Clercq(2009_  ", " Van Den Heuvel(1963) ", 
                 " Parker(1982)  ", " Bryant   Thomas & Bale (1997  )", 
                 " Bryant Thomas & Bale (1997  )", " Bryant Thomas & Bale(1997)  ", 
                 " Carroll & Hoyt(1986)  ", " Forster et al (2011)", " Fischer et al. (2003)  ", 
                 " Fischer Brakefield & Zwaan( 2003)  ", " Fischer et al. 2004  ", 
                 " Campbell et al. 2001  ", " Uye (1988  )", " Stillwell et al. (2007  )", 
                 " Sulkin & McKeen (1994) ", " Frouz Ali &  Lobinske (2002  )", 
                 " Konstantinov (1958  )", " Strixino & Strixino (1985  )", 
                 " Reigada & Godoy 2006  ", " Sweeney & Vannote 1984 ", " Pechenik 1984  ", 
                 " Pechenik & Lima (1984  )", " Mullens & Rutz (1983  )", 
                 " Akey Potter &Jones (1978  )", " Tsitsipis (1980)  ", 
                 " Dufresne & Hebert 1998  ", " Brambilla 1982  ", " Walls & Ventelä (1998) ", 
                 " Alvarez & Nicieza(2002)  ", " de Jong(2010) ", " David et al.( 1997)  ", 
                 " David et al.(1997  )", " David et al(1997)  ", " Ray (1960  )", 
                 " Ray 1960  ", " Ray(1960)  ", " Crill Huey & Gilchrist 1996 ", 
                 "Karan etal (1998)", " Karan et al. (1998)  ", " Petavy et al. (2001)  ", 
                 " Petavy et al.(2001  )", " Tantawy & Mallah (1961 )", " Neel (1940  )", 
                 " David et al. (1994  )", " David Legout& Moreteau (2006  )", 
                 " David & Clavel (1967  )", " Druger (1962 )", " Villarreal & Hernandez-Llamas (2005  )", 
                 " Stam van de Leemkule & Ernsting(1996 )", " Larsen Madden & Nault (1990  )", 
                 " Tsukada Asai & Higuchi(2005  )", " Kunisch & Anger (1984  )", " Anger (1984  )", 
                 " Bryant   Thomas & Bale (1997  )", " Bryant Thomas & Bale (1997  )", 
                 " Bryant Thomas & Bale(1997)  ", " De Block & Stoks (2003  )", 
                 " Fischer & Fiedler (2002  )", " Karl Janowitz& Fischer (2008  )", 
                 " Tependino & Parker (1986  )", " Radmacher &Strohm (2009 )", 
                 " Ouellet & Chabot (2005 )", " Uye (1991  )", " Shuster & Guthrie (1999  )", 
                 " Seikai (1986  )", " Soderstrom (1988  )", " Seko et al. 2005  ", 
                 " Lebedeva & Gerasimova (1985  )", " Rahman (1969 )", " Jones Hart & Bull (1982)  ", 
                 " Miller 1977  ", " Harkey & Semlitsch (1988  )", " Lee et al. (2003  )", 
                 " Li et al. (2009 )", " Allsopp (1981 )", " Orizaola & Laurila( 2009  )", 
                 " Smith-Gill & Berven (1979  )", " Berven & Gill(1983  )", 
                 " Laugen et al. (2005  )", " Blanckenhorn(2000  )", " Lonsdale & Levinton (1985  )", 
                 " Perrin (1988  )", " Kimoto Uye &Onbe (1986  )", " Stillwell & Fox (2005  )", 
                 " Stillwell Moya-Larano& Fox (2008  )", " Stelzer 2002 ", 
                 " Whitman (1986  )", " Karan Moreteau & David (1999  )")

pairs_data = pairs_data %>% 
  mutate(forster = ifelse(ref %in% forster_list, "Y", "N"))

addtl_refs = pairs_data %>% 
  filter(forster == "N") %>% 
  group_by(ref) %>% 
  summarize(counts = n())

pairs_data$Class[pairs_data$Class == "Actinoperygii"] = "Actinopterygii"

### Dataset metrics
range(pairs_data$initial_temp)
range(pairs_data$final_temp)
nrow(pairs_data)
length(unique(pairs_data$species))
length(unique(pairs_data$study))
length(unique(pairs_data$Class))
range(pairs_data$initial_mass_mg)
min(pairs_data$initial_mass_mg) * 1000 #ug
max(pairs_data$initial_mass_mg) / 1000 #g

### All plot calculations
pairs_data$needed_mass_change = (pairs_data$constantmetrate_mass - pairs_data$initial_mass) / abs(pairs_data$initial_mass) * 100
pairs_data$actual_mass_change = (pairs_data$final_mass - pairs_data$initial_mass) / abs(pairs_data$initial_mass) * 100
pairs_data$mass_constant_slope = -pairs_data$Ea / .00008617
pairs_data$mass_constant_intercept = 0
pairs_data$initial_temp_K = pairs_data$initial_temp + 273.15
pairs_data$final_temp_K = pairs_data$final_temp + 273.15
pairs_data$temp_axis_num = -(pairs_data$initial_temp_K - pairs_data$final_temp_K)
pairs_data$temp_axis_den = pairs_data$initial_temp_K * pairs_data$final_temp_K
pairs_data$temp_axis = pairs_data$temp_axis_num / pairs_data$temp_axis_den
pairs_data$metab_axis = log(pairs_data$final_metrate) - log(pairs_data$initial_metrate)
pairs_data$expect_metrate_change = pairs_data$temp_axis * pairs_data$mass_constant_slope
pairs_data$observ_metrate_change = pairs_data$metab_axis

### MAIN: Compensation mass plot
x_text1 = textGrob("Decreasing needed mass", gp = gpar(cex = .8))
x_text2 = textGrob("Increasing needed mass", gp = gpar(cex = .8))
y_text1 = textGrob("Increasing actual mass", gp = gpar(cex = .8), rot = 90)
y_text2 = textGrob("Decreasing actual mass", gp = gpar(cex = .8), rot = 90)

plot_comp_mass = ggplot(pairs_data, aes(x = needed_mass_change, y = actual_mass_change)) +
  geom_point(size = 0.8) +
  coord_cartesian(xlim = c(-100, 100), ylim = c(-100, 100)) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_vline(xintercept = 0, color = "grey") +
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Compensation mass change (%)", y = "Observed mass change (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x = element_text(face = "bold", margin = unit(c(1.65, 0, 0, 0), "cm")), 
        axis.title.y = element_text(face = "bold", margin = unit(c(0, 1.55, 0, 0), "cm"))) +
  annotation_custom(grob = x_text1, xmin = -50, xmax = -50, ymin = -130, ymax = -130) +
  annotation_custom(grob = x_text2, xmin = 50, xmax = 50, ymin = -130, ymax = -130) +
  annotation_custom(grob = y_text1, xmin = -135, xmax = -135, ymin = 50, ymax = 50) +
  annotation_custom(grob = y_text2, xmin = -135, xmax = -135, ymin = -50, ymax = -50) +
  annotation_custom(grob = linesGrob(arrow = arrow(angle = 20, length = unit(0.15, "inches"), type = "open", ends = "first")), xmin = -100, xmax = -10, ymin = -135, ymax = -135) +
  annotation_custom(grob = linesGrob(arrow = arrow(angle = 20, length = unit(0.15, "inches"), type = "open", ends = "last")), xmin = 100, xmax = 10, ymin = -135, ymax = -135) +
  annotation_custom(grob = linesGrob(arrow = arrow(angle = 20, length = unit(0.15, "inches"), type = "open", ends = "last")), xmin = -140, xmax = -140, ymin = 10, ymax = 100) +
  annotation_custom(grob = linesGrob(arrow = arrow(angle = 20, length = unit(0.15, "inches"), type = "open", ends = "first")), xmin = -140, xmax = -140, ymin = -100, ymax = -10) +
  annotate("text", x = 75, y = 25, label = "Actual mass equals\n needed mass") +
  annotation_custom(grob = linesGrob(arrow = arrow(angle = 20, length = unit(0.15, "inches"), type = "open", ends = "first")), xmin = 25, xmax = 50, ymin = 20, ymax = 20)

plot_comp_mass <- ggplotGrob(plot_comp_mass)
plot_comp_mass$layout$clip[plot_comp_mass$layout$name=="panel"] <- "off"
grid.draw(plot_comp_mass)
ggsave("figures/fig2.jpg", plot = plot_comp_mass, width = 6, height = 6)

### STATS: Mass comparison t-test
mean(pairs_data$actual_mass_change)
sd(pairs_data$actual_mass_change)

comp_mass_df = data.frame(decrease_insuff = length(which((pairs_data$actual_mass_change > pairs_data$needed_mass_change) & pairs_data$actual_mass_change < 0)), 
                          decrease_exceed = length(which(pairs_data$actual_mass_change < pairs_data$needed_mass_change)), 
                          increase = length(which(pairs_data$actual_mass_change > 0)), 
                          same = length(which(pairs_data$actual_mass_change == 0)))
comp_mass_df = data.frame(t(comp_mass_df))
comp_mass_df = tibble::rownames_to_column(comp_mass_df)
comp_mass_df = comp_mass_df %>% 
  rename(mass = rowname, 
         number = t.comp_mass_df.) %>% 
  mutate(percent = number / nrow(pairs_data) * 100)

t.test(log(pairs_data$final_mass), log(pairs_data$constantmetrate_mass), paired = TRUE)

pairs_data_subset = pairs_data %>% 
  group_by(studyID) %>% 
  filter(temp_axis_num == min(temp_axis_num))
t.test(log(pairs_data_subset$final_mass), log(pairs_data_subset$constantmetrate_mass), paired = TRUE)

### MAIN: Variance plot with trend line
pairs_data = pairs_data %>% 
  mutate(point_color = ifelse(Class == "Amphibia", "Amphibians", "Non-amphibians"))

plot_metrates = ggplot(pairs_data, aes(x = expect_metrate_change, y = observ_metrate_change)) +
  geom_point(aes(color = point_color)) +
  scale_color_manual(values = c("dark orange", "blue",  "black")) +
  geom_abline(aes(color = "Size does not change", intercept = 0, slope = 1)) +
  geom_hline(yintercept = 0, color = "grey") +
  guides(color = guide_legend(override.aes = list(linetype = c("blank", "blank", "solid"), 
                                                  shape = c(16, 16, NA)))) +
  labs(x = "Metabolic rate change with constant size", y = "Metabolic rate change with varying size", color = "") +
  xlim(min(pairs_data$observ_metrate_change), max(pairs_data$expect_metrate_change)) +
  ylim(min(pairs_data$observ_metrate_change), max(pairs_data$expect_metrate_change)) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("figures/fig3.jpg", plot = plot_metrates)

### STATS: Variance explained by no mass change line
metrate_df = data.frame(decrease_exceed = length(which(round(pairs_data$observ_metrate_change, 7) < 0)), 
                        increase = length(which(round(pairs_data$observ_metrate_change, 7) > round(pairs_data$expect_metrate_change, 7))), 
                        decrease_insuff = length(which((round(pairs_data$observ_metrate_change, 7) < round(pairs_data$expect_metrate_change, 7)) & round(pairs_data$observ_metrate_change, 7) > 0)), 
                        same = length(which(round(pairs_data$observ_metrate_change, 7) == round(pairs_data$expect_metrate_change, 7))))

metrate_df = data.frame(t(metrate_df))
metrate_df = tibble::rownames_to_column(metrate_df)
metrate_df = metrate_df %>% 
  rename(metrate = rowname, 
         number = t.metrate_df.) %>% 
  mutate(percent = number / nrow(pairs_data) * 100)

unexplained_variance = var(with(pairs_data, observ_metrate_change - expect_metrate_change))
total_variance = var(with(pairs_data, observ_metrate_change))
r2 = 1 - unexplained_variance / total_variance

### MAIN: With mass percent change lines
pairs_data$mass_intercept_1 = pairs_data$exponent * log(2)
pairs_data$mass_intercept_2 = pairs_data$exponent * log(1.5)
pairs_data$mass_intercept_3 = pairs_data$exponent * log(0.75)
pairs_data$mass_intercept_4 = pairs_data$exponent * log(0.5)

add_grob_outside <- function(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data){
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = FALSE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax))
}

plot_classes = ggplot(data = pairs_data, aes(x = temp_axis, y = metab_axis)) +
  geom_point(size = 0.8) +
  geom_abline(data = pairs_data, aes(slope = mass_constant_slope, intercept = mass_constant_intercept)) +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_intercept_1), color = "grey") +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_intercept_2), color = "grey") +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_intercept_3), color = "grey") +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_intercept_4), color = "grey") +
  facet_wrap(~Class) +
  labs(x = expression(atop("Temperature axis", -(T[1]-T[2])/(T[1]*T[2]))), 
       y = expression(atop("Metabolic rate change with varying size", log(R[2]) - log(R[1])))) +
  add_grob_outside(grob = textGrob("Increasing temperature difference", gp = gpar(cex = .6)), 
                   xmin = 0.00016, xmax = 0.00016, ymin = -1.85, ymax = -1.85, data = data.frame(Class = "Malacostraca")) +
  add_grob_outside(grob = linesGrob(arrow = arrow(angle = 20, length = unit(0.15, "inches"), type = "open", ends = "last")), 
                   xmin = 0.00001183752, xmax = 0.00033713062, ymin = -2.2, ymax = -2.2, 
                   data = data.frame(Class = "Malacostraca")) +
  theme_bw() +
  theme(axis.title.x = element_text(margin = unit(c(1, 0, 0, 0), "cm")),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_classes <- ggplotGrob(plot_classes)
plot_classes$layout[grepl("panel", plot_classes$layout$name), ]$clip <- "off"
grid.draw(plot_classes)

ggsave("figures/fig4.jpg", plot = plot_classes, width = 7, height = 5)

### SUPPLEMENT: Temp difference mostly affects x-axis
plot_temp_axis = ggplot(pairs_data, aes(x = temp_axis_num, y = temp_axis)) +
  geom_point(size = 0.8) +
  facet_wrap(~Class) +
  labs(x = expression(atop("Temperature difference", T[1]-T[2])), 
       y = expression(atop("Temperature axis", -(T[1]-T[2])/(T[1]*T[2])))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("figures/supp8.jpg", plot = plot_temp_axis, width = 7, height = 5)

### SUPPLEMENT: With data trend line
plot_class_trend = ggplot(pairs_data, aes(x = temp_axis, y = metab_axis)) +
  geom_point(size = 0.8) +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_constant_intercept)) +
  geom_smooth(method = "lm", se = FALSE, size = .6, color = "red") +
  facet_wrap(~Class) +
  labs(x = expression(atop("Temperature axis", -(T[1]-T[2])/(T[1]*T[2]))), 
       y = expression(atop("Metabolic rate change with varying size", log(R[2]) - log(R[1])))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("figures/supp6.jpg", plot = plot_class_trend, width = 7, height = 5)

### Full model
pairs_data$residual = pairs_data$metab_axis - (pairs_data$mass_constant_slope * pairs_data$temp_axis)

pairs_data$temp_axis_rescale = pairs_data$temp_axis * 100000
full_model = lmer(residual ~ log(initial_mass) + initial_temp + temp_axis_rescale + (1|species) + 
                    (1|Class) + (1|studyID) + (1|study), data = pairs_data)
summary(full_model)

#Likelihood ratio tests
null_mass = lmer(residual ~ initial_temp + temp_axis_rescale + (1|species) + 
                   (1|Class) + (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
model_mass = lmer(residual ~ log(initial_mass) + initial_temp + temp_axis_rescale + (1|species) + 
                    (1|Class) + (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
anova(null_mass, model_mass)

null_temp_init = lmer(residual ~ log(initial_mass) + temp_axis_rescale + (1|species) + 
                        (1|Class) + (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
model_temp_init = lmer(residual ~ log(initial_mass) + initial_temp + temp_axis_rescale + (1|species) + 
                         (1|Class) + (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
anova(null_temp_init, model_temp_init)

### STATS: Linear mixed model
final_model = lmer(residual ~ temp_axis_rescale + (1|species) + (1|Class) + 
                     (1|studyID) + (1|study), data = pairs_data)
summary(final_model)

random_effects = as.data.frame(VarCorr(final_model))
random_effects_vars = random_effects$vcov
residuals_var = var(pairs_data$residual)
random_effects_vars / residuals_var * 100

### SUPPLEMENT: Linear mixed model diagnostics
diag_df = data.frame(Residuals = residuals(final_model), Fitted = fitted(final_model))

diag_linear_hetero = ggplot(diag_df, aes(x = Fitted, y = Residuals)) +
  geom_point() +
  geom_hline(yintercept = 0)

diag_normal1 = ggplot(diag_df, aes(x = Residuals)) +
  geom_histogram(bins = 50) +
  xlim(-.5, .5)

diag_normal2 = ggplot(diag_df, aes(sample = Residuals)) +
  stat_qq()

diags = plot_grid(diag_linear_hetero, diag_normal1, diag_normal2, labels = c("A", "B", "C"))
ggsave("figures/supp5.jpg", plot = diags, width = 7, height = 5)

#Likelihood ratio tests
null_temp_diff = lmer(residual ~ (1|species) + (1|Class) + 
                        (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
model_temp_diff = lmer(residual ~ temp_axis_rescale + (1|species) + (1|Class) + 
                         (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
anova(null_temp_diff, model_temp_diff)

null_sp = lmer(residual ~ temp_axis_rescale + (1|Class) + 
                 (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
model_sp = lmer(residual ~ temp_axis_rescale + (1|species) + (1|Class) + 
                  (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
anova(null_sp, model_sp)

null_class = lmer(residual ~ temp_axis_rescale + (1|species) +  
                    (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
model_class = lmer(residual ~ temp_axis_rescale + (1|species) + (1|Class) + 
                     (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
anova(null_class, model_class)

null_studyID = lmer(residual ~ temp_axis_rescale + (1|species) + (1|Class) + 
                      (1|study), data = pairs_data, REML = FALSE)
model_studyID = lmer(residual ~ temp_axis_rescale + (1|species) + (1|Class) + 
                       (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
anova(null_studyID, model_studyID)

null_study = lmer(residual ~ temp_axis_rescale + (1|species) + (1|Class) + 
                    (1|studyID), data = pairs_data, REML = FALSE)
model_study = lmer(residual ~ temp_axis_rescale + (1|species) + (1|Class) + 
                     (1|studyID) + (1|study), data = pairs_data, REML = FALSE)
anova(null_study, model_study)

### SUPPLEMENT: No systematic variation due to absolute mass or temp
residuals_mass = ggplot(pairs_data, aes(x = log(initial_mass), y = residual)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  facet_wrap(~Class, scales = "free_x") +
  ylim(-1.4, 1.4) +
  labs(x = "log(mass)", y = "Mass residual") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

residuals_temp = ggplot(pairs_data, aes(x = initial_temp_K, y = residual)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  facet_wrap(~Class, scales = "free_x") +
  ylim(-1.4, 1.4) +
  labs(x = "Initial temperature (K)", y = "Mass residual") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

residuals = plot_grid(residuals_mass, residuals_temp, labels = c("A", "B"), ncol = 1)
ggsave("figures/supp4.jpg", plot = residuals, width = 8, height = 10)

### SUPPLEMENT: temperature-dependent exponent for fishes
fishes_data = filter(pairs_data, Class == "Actinopterygii")

fishes_data$initial_exponent = 0.0417 * (1/(.0000862*(fishes_data$initial_temp + 273.15))) - 0.8923
fishes_data$final_exponent = 0.0417 * (1/(.0000862*(fishes_data$final_temp + 273.15))) - 0.8923

fishes_data$initial_metrate_adjusted = (fishes_data$initial_mass ^ fishes_data$initial_exponent) * (exp(((fishes_data$Ea) / (.00008617 * (fishes_data$initial_temp + 273.15)))))
fishes_data$final_metrate_adjusted = (fishes_data$initial_mass ^ fishes_data$final_exponent) * (exp(((fishes_data$Ea) / (.00008617 * (fishes_data$final_temp + 273.15)))))
fishes_data$metab_axis_adjusted = fishes_data$final_metrate_adjusted - fishes_data$initial_metrate_adjusted

fishes_data$mass_intercept_1 = fishes_data$exponent * log(2)
fishes_data$mass_intercept_2 = fishes_data$exponent * log(1.5)
fishes_data$mass_intercept_3 = fishes_data$exponent * log(0.75)
fishes_data$mass_intercept_4 = fishes_data$exponent * log(0.5)

ggplot(fishes_data, aes(x = temp_axis, y = metab_axis)) +
  geom_point() +
  geom_point(aes(x = temp_axis, y = metab_axis_adjusted), color = "red", alpha = 0.7) +
  geom_abline(intercept = 0, lty = 3) +
  coord_cartesian(xlim = c(0, 0.0003), ylim = c(-1, 2)) +
  geom_abline(data = fishes_data, aes(slope = mass_constant_slope, intercept = mass_constant_intercept)) +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_intercept_1), color = "grey") +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_intercept_2), color = "grey") +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_intercept_3), color = "grey") +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_intercept_4), color = "grey") +
  labs(x = "Temperature axis \n -(T1-T2)/(T1*T2)", y = "Metabolic rate difference \n log(R2) - log(R1)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("figures/supp2.jpg")

