library(dplyr)
library(ggplot2)
library(lme4)
library(grid)

pairs_data = read.csv("clean_data_MTEecto_expanded.csv", stringsAsFactors = FALSE)

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
  geom_point() +
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
comp_mass_df = comp_mass_df %>% 
  rename(number = t.comp_mass_df.) %>% 
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
  guides(color = guide_legend(override.aes = list(linetype = c("blank", "blank", "solid"), 
                                                  shape = c(16, 16, NA)))) +
  labs(x = "Metabolic rate change with constant size", y = "Metabolic rate change with varying size", color = "") +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("figures/fig3.jpg", plot = plot_metrates)

### STATS: Variance explained by no mass change line
unexplained_variance = var(with(pairs_data, observ_metrate_change - expect_metrate_change))
total_variance = var(with(pairs_data, observ_metrate_change))
r2 = 1 - unexplained_variance / total_variance

### MAIN: With mass percent change lines
pairs_data$mass_intercept_1 = pairs_data$exponent * log(2)
pairs_data$mass_intercept_2 = pairs_data$exponent * log(1.5)
pairs_data$mass_intercept_3 = pairs_data$exponent * log(0.75)
pairs_data$mass_intercept_4 = pairs_data$exponent * log(0.5)

annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data){
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax))
}

data_1 = data.frame(Class = "Malacostraca")

plot_classes = ggplot(pairs_data, aes(x = temp_axis, y = metab_axis)) +
  geom_point() +
  geom_abline(data = pairs_data, aes(slope = mass_constant_slope, intercept = mass_constant_intercept)) +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_intercept_1), color = "grey") +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_intercept_2), color = "grey") +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_intercept_3), color = "grey") +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_intercept_4), color = "grey") +
  facet_wrap(~Class) +
  labs(x = expression(atop("Increasing temperature difference", -(T[1]-T[2])/(T[1]*T[2]))), 
       y = expression(atop("Metabolic rate change with varying size", log(R[2]) - log(R[1])))) +
  #annotation_custom2(grob = linesGrob(arrow = arrow(type = "open", ends = "last")), xmin = 0, xmax = .0003, ymin = -2, ymax = -2, data = data_1) +
  annotation_custom2(grob = linesGrob(), xmin = 0, xmax = .0003, ymin = -2, ymax = -2, data = data_1) +
  theme_bw() +
  theme(axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "cm")),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_classes <- ggplotGrob(plot_classes)
plot_classes$layout$clip[plot_classes$layout$name=="panel"] <- "off"
grid.draw(plot_classes)

#Code from here: https://stackoverflow.com/questions/10525957/how-to-draw-lines-outside-of-plot-area-in-ggplot2
test = data.frame(
  group=c(rep(1,6), rep(2,6)),
  subgroup=c( 1,1,1,2,2,2,1,1,1,2,2,2),
  category=c( rep(1:3, 4)),
  count=c( 10,80,10,5,90,5,  10,80,10,5,90,5   )
)
p <- ggplot(test) +
  geom_bar(aes(subgroup, count, fill = category), stat = "identity") +
  facet_grid(. ~ group) +
  theme(legend.position = "none",  
        plot.margin = unit(c(1,5,1,1), "lines"))
data = data.frame(group=2)
annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data){
    layer(data = data, stat = StatIdentity, position = PositionIdentity, 
          geom = ggplot2:::GeomCustomAnn,
          inherit.aes = TRUE, params = list(grob = grob, 
                                            xmin = xmin, xmax = xmax, 
                                            ymin = ymin, ymax = ymax))
  }
p1 = p + annotation_custom2(linesGrob(), xmin = 2.75, xmax = 2.75, ymin = 85, ymax = 100, data = data)
gt <- ggplotGrob(p1)
gt$layout[grepl("panel", gt$layout$name), ]$clip <- "off"
grid.draw(gt)



### SUPPLEMENT: Temp difference mostly affects x-axis
ggplot(pairs_data, aes(x = temp_axis_num, y = temp_axis)) +
  geom_point() +
  facet_wrap(~Class) +
  labs(x = "Temperature difference \n T1 - T2", y = "Temperature axis \n (T1-T2)/(T1*T2)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### SUPPLEMENT: With data trend line
ggplot(pairs_data, aes(x = temp_axis, y = metab_axis)) +
  geom_point() +
  geom_abline(aes(slope = mass_constant_slope, intercept = mass_constant_intercept)) +
  geom_smooth(method = "lm", se = FALSE, size = .6) +
  facet_wrap(~Class) +
  labs(x = "Temperature axis \n (T1-T2)/(T1*T2)", y = "Metabolic rate difference \n log(R2)-log(R1)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Full model
#pairs_data$residual_alt = pairs_data$metab_axis - (pairs_data$mass_constant_slope * pairs_data$temp_axis)
pairs_data$residual = pairs_data$exponent * (log(pairs_data$final_mass) - log(pairs_data$initial_mass))

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

### SUPPLEMENT: Linear mixed model diagnostics
coef(final_model)
plot(fitted(final_model), residuals(final_model)) #linearity good, heteroscedasticity
abline(h = 0)
hist(residuals(final_model)) # normality looks good in histogram
qqnorm(residuals(final_model))
qqline(residuals(final_model)) #normality of residuals heavy-tailed in qqplot

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
ggplot(pairs_data, aes(x = log(initial_mass), y = residual)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  facet_wrap(~Class, scales = "free_x") +
  ylim(-1.4, 1.4) +
  labs(x = "log(mass)", y = "Mass residual") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(pairs_data, aes(x = initial_temp_K, y = residual)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  facet_wrap(~Class, scales = "free_x") +
  ylim(-1.4, 1.4) +
  labs(x = "Initial temperature (K)", y = "Mass residual") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
