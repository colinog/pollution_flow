#gitcreds::gitcreds_set()
#theme_set(theme_bw())
#renv::snapshot()
library(tidyverse)
library(glmmTMB)  
library(openxlsx)
library(DHARMa)
library(car)
library(emmeans)
library(rstatix)
pj = position_dodge(width = 0.3)

##### Data visualization ########
#I included the daily pesticide class fluctuation
#for the entry into the flumes
##### Pesticide concentration in water ########
pest_water = readRDS("All_code_data_visualization/Data/pesticide_water.rds")

#####Summary of pesticide concentration in water######

water_sum = pest_water |> 
  filter(!Concentration_µg_L %in% c("<LOQ", NA)) |> #exclude where conc. is <LOQ and NA
  mutate(across(c(Concentration_µg_L, LOQ), as.numeric)) |> # convert the column to numeric variable
  filter(!Class %in% "Metabolite") |> # exclude were class was metabolite
  rstatix::group_by(Analyte, Class) |> # group by date and pesticide class
  rstatix::get_summary_stats(Concentration_µg_L, type = "full") |>  #summarise by concentration
  arrange(Class)
  
openxlsx::write.xlsx(water_sum, "All_code_data_visualization/95_pesticide/water_pest_sum.xlsx", rowNames = F)

####Time series plot of daily pesticide concentration in water ######



pest_water |> 
  filter(!Concentration_µg_L %in% c("<LOQ", NA)) |> #exclude where conc. is <LOQ and NA
  mutate(across(c(Concentration_µg_L, LOQ), as.numeric)) |> # convert the column to numeric variable
  filter(!Class %in% "Metabolite") |> # exclude were class was metabolite
  rstatix::group_by(date, Class) |> # group by date and pesticide class
  rstatix::get_summary_stats(Concentration_µg_L, type = "full") |> #summarise by concentration
  ggplot(aes(date, median,  shape = Class, group = Class))+
  geom_line()+
  geom_point()+
  scale_shape_manual(values = c(5,6,7))+
  labs(y = "Median pesticide concentrations (µg/L)", x = "Sampling days")+
  facet_grid(Class~.)+ #Separated by pesticide class with same y axis value
  scale_x_date(breaks = seq(min(pest_water$date), max(pest_water$date), by = "6 days"), #date sequence showing every six days
               labels = function(x) format(x, "%d-%b")  # Show only month and day
  )+
  geom_vline(
    xintercept = as.Date("2021-06-15"),  # Add a vertical line at 15 June
    linetype = "dashed", color = "black"   # Customize the line appearance
  ) +
  annotate(
    "text", 
    x = as.Date("2021-06-15"), y = 0.08,  # Position the text on the plot
    label = "Start of Low-flow", 
    angle = 90, hjust = 1.2, vjust = -0.5, size = 4, color = "black"
  )+
  theme(panel.background = element_blank(),
        strip.text.y =  element_text(angle = 90, size = 11, face = "bold"),
        axis.line = element_line(color = "black"), 
        axis.title.y = element_text(face = "bold"),
        axis.title = element_text(size = 12, color = "black",
                                  face = "bold"),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        legend.position = "bottom",
        legend.text = element_text(size = 11, face = "bold", color = "black"),
        legend.title = element_text(size = 11),
        plot.tag = element_text(face = "bold"),
  )

ggsave("All_code_data_visualization/plot/pest_water.png", dpi = 300, width = 22, height = 18, units = "cm")

####### Pesticide concentration in sediment under control and low-flow treatments

sed_pest = readRDS("All_code_data_visualization/Data/sediment_pesticide.rds")
#####Summary of sediment pesticide ###########
sed_sum = sed_pest |> 
  dplyr::group_by(Treatment,Class,Pesticides  ) |> 
  rstatix::get_summary_stats(Concentration_µg_kg, type = "full") |> 
  mutate(Treatment = case_when(
    Treatment == "C" ~ "Control",
    Treatment == "D" ~ "Low-flow",
    TRUE ~ Treatment  # Keeps other values unchanged
  ))

openxlsx::write.xlsx(sed_sum, "All_code_data_visualization/95_pesticide/sed_pest_sum.xlsx", rowNames = F)

#####Visualize by mean pesticide concentration under control and low-flow


#The treatment column is a categorical variable with control (C) and low-flow (D)
sed_pest |> 
  dplyr::group_by(Treatment, Week) |> 
  rstatix::get_summary_stats(Concentration_µg_kg, type = "full")|> 
  ggplot(aes(Treatment, mean, fill = Week, shape = Week, linetype = Week))+
  geom_point(position = pj, size = 3, color = "black")+
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), 
                width = 0.1, position = pj, color = "black") +
  scale_shape_manual(values = c(5,6))+ #shape by treatment in week 4 and 6
  ylim(0, 2.5)+ #y-axis limit
  scale_x_discrete(labels = c("Control", "Low-flow"))+
  labs(y = expression(bold("Mean pesticide concentration (µg kg"^-1*")")), 
       x = "Treatment")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.title.y = element_text(face = "bold"),
    axis.title = element_text(size = 12, color = "black",
                              face = "bold"),
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 12, face = "bold", color = "black"),
    legend.title = element_text(size = 12),
    plot.tag = element_text(face = "bold"),
  )

ggsave("All_code_data_visualization/plot/pest_sed.png", dpi = 300, width = 15, height = 12, units = "cm")

####### The mean sum concentration to estimate percentage difference ######
#aggregate(Concentration_µg_kg ~ Treatment + Week, data = sed_pest, FUN = mean)

avg_sum = sed_pest |> 
  dplyr::group_by(Treatment, Week) |> 
  rstatix::get_summary_stats(Concentration_µg_kg, type = "mean_ci")|> 
  mutate(Treatment = case_when(
    Treatment == "C" ~ "Control",
    Treatment == "D" ~ "Low-flow"
  )) |> 
  mutate(Week = case_when(
    Week == 4 ~ "Wk4",
    Week == 6 ~ "Wk6"
  ))

#### Compute the percentage difference for both treatment by week ##### 

p_diff = avg_sum |> 
  select(Treatment, Week, mean) |> 
  pivot_wider(names_from = Week, values_from = mean) |> 
  mutate(per_diff = (Wk6 - Wk4) / Wk4 * 100)

print(p_diff) #



# library(ggridges)
# ggplot(sed_pest, aes(x = Concentration_µg_kg, y = Class, fill = Class)) + 
#   geom_density_ridges(rel_min_height = 0.01)+
#   scale_fill_viridis_d()+
#   theme(legend.position = "none")

#####Wet emergence sample collected over five weeks #########

abu_bio = readRDS("All_code_data_visualization/Data/wet_sample.rds") #read in emergence and biomass data

######### standardise abundance visualization ##############
abu_bio |> 
  dplyr::group_by(Treatment, Week) |>  #group by Treatment and Week
  rstatix::get_summary_stats(CPUE, type = "full")|> # summarise abundance using catch per unit effort (CPUE) z
  ggplot(aes(Week,mean, fill = Treatment, color = Treatment, 
             shape = Treatment, linetype = Treatment, group = Treatment))+
  geom_point( position = pd, size = 3, color = "black")+
  geom_line(position = pd, color = "black", linewidth = 0.3)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.1, position = pd, color = "black") +
  scale_fill_viridis_d(option = "plasma")+
  scale_linetype_manual(values = c("solid", "dashed"))+
  labs(y = expression(bold("Abundance (ind. m"^-2*" day"^-1*")")), 
       x = expression(bold("Sampling Week")))+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.title.y = element_text(face = "bold"),
    axis.title = element_text(size = 13, color = "black",
                              face = "bold"),
    axis.text = element_text(size = 13, face = "bold", color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 12, face = "bold", color = "black"),
    legend.title = element_text(size = 12),
    strip.text = element_text(size = 13, face = "bold"),
    plot.tag = element_text(face = "bold")
  )

ggsave("All_code_data_visualization/plot/abundance.png", dpi = 300, width = 15, height = 12, units = "cm") 
########## standardize biomass visualization #####
abu_bio |> 
  dplyr::group_by(Treatment, Week) |>  #group by Treatment and Week
  rstatix::get_summary_stats(mass_flux, type = "full")|> # summarise biomass using mass flux
  ggplot(aes(Week,mean, fill = Treatment, color = Treatment, 
             shape = Treatment, linetype = Treatment, group = Treatment))+
  geom_point( position = pd, size = 3, color = "black")+
  geom_line(position = pd, color = "black", linewidth = 0.3)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.1, position = pd, color = "black") +
  scale_fill_viridis_d(option = "plasma")+
  scale_linetype_manual(values = c("solid", "dashed"))+
  labs(y = expression(bold("Biomass (mg m"^-2*" day"^-1*")")), 
       x = expression(bold("Sampling Week")))+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.title.y = element_text(face = "bold"),
    axis.title = element_text(size = 13, color = "black",
                              face = "bold"),
    axis.text = element_text(size = 13, face = "bold", color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 12, face = "bold", color = "black"),
    legend.title = element_text(size = 12),
    strip.text = element_text(size = 13, face = "bold"),
    plot.tag = element_text(face = "bold")
  )

ggsave("All_code_data_visualization/plot/biomass.png", dpi = 300, width = 15, height = 12, units = "cm")
###### spider abundance #########

spider_dat = readRDS("All_code_data_visualization/Data/spider_data.rds")

spider_dat |> 
  dplyr::group_by(Treatment,Family) |> # group by Treatment and Family
  rstatix::get_summary_stats(count, type = "full")|> # summarize by the spider count
  ggplot(aes(Treatment, mean, fill = Family, color = Family, 
             shape = Family, linetype = Family))+
  geom_point(position = pj, size = 3, color = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.1, position = pj, color = "black") +
  ylim(0, 50)+
  scale_shape_manual(values = c(0,8))+
  labs(y = "Mean abundance per riparian area", x = "Treatment")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.title.y = element_text(face = "bold"),
    axis.title = element_text(size = 12, color = "black",
                              face = "bold"),
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 12, face = "bold", color = "black"),
    legend.title = element_text(size = 12),
    plot.tag = element_text(face = "bold"),
  )

ggsave("All_code_data_visualization/plot/spider.png", dpi = 300, width = 15, height = 12, units = "cm")


######Statistical analysis of all tested hypothesis ######
####################Sediment Pesticide Concentration########################
sed_pest$Concentration <- 1/sqrt(sed_pest$Concentration_µg_kg) 

#inverse square root transformation of response variable concentration
#because initial model with the variable was not normally distributed

fp = glmmTMB(Concentration ~  Treatment + Week
             + (1|Flume),
             data = sed_pest, lognormal(link = "log")) 
summary(fp)

##initial model to check for spatial autocorrelation
res_pest <- residuals(fp) #extract residual from model
###calculate and plot the acf residuals
auto_cor_pest <- acf(res_pest, lag.max=40, plot=FALSE)
plot(auto_cor_pest, main="ACF of Residuals", xlab="Lag", ylab="ACF")

# Perform Ljung-Box test on residuals to confirm when necessary
#Box.test(residuals, lag = 20, type = "Ljung-Box")


#Treatment factors are: C = Control and D = Low-flow
fit_sed = glmmTMB(Concentration ~  Treatment + Week
                  + (1|Flume) + ar1(Week+0|Flume),
                  dispformula = ~Week, ##Week was added as dispersion factor
                  data = sed_pest, family = Gamma(link = "log"))

##The interaction terms were not significant in the initial model



summary(fit_sed)
diagnose(fit_sed)
testDispersion(fit_sed)
res_sim = simulateResiduals(fit_sed)
plot(res_sim)
print(res_sim)
testResiduals(res_sim)
plot(res_sim)
plotQQunif(res_sim)
plotResiduals(res_sim)


#Fit type 2 Anova
Anova(fit_sed, type="II") #The factors were not significant 

#I can continue with the pairwise comparison or report the summary of the glmm model
mt = emmeans(fit_sed, specs = pairwise ~ Treatment,
             adjust = "bonferroni", type = "response")

mt$contrasts %>%
  rbind()

mtw = emmeans(fit_sed, specs = pairwise ~ Week|Treatment,
              adjust = "bonferroni", type = "response")


mtw$contrasts %>%
  rbind()



####Emergence data for both mass flux and emergence rate
f = glmmTMB(CPUE ~ Treatment + Week +
              (1 | Flume),
            data = abu_bio,
            family = Gamma(link = "log")) ##fit glmm model

res_ab <- residuals(f) # extract residuals
#calculate and plot the act
auto_cor_ab <- acf(res_ab, lag.max=40, plot=FALSE)
plot(auto_cor_ab, main="ACF of Residuals", xlab="Lag", ylab="ACF")

#The model residuals shows serial autocorrelation

## Perform Ljung-Box test on residuals to confirm when necessary
#Box.test(residuals, lag = 20, type = "Ljung-Box")

##run the model with sampling time and ar1, the interaction terms were not significant

fit_ = glmmTMB(CPUE ~ Treatment + Week +
                 (1 | Flume) + ar1(Week+0|Flume), #sparial autocorrelation structure
               data = abu_bio,
               family = Gamma(link = "log")) ##We use Gamma with log link function 

summary(fit_) ##model summary
diagnose(fit_)
testDispersion(fit_) ##no dispersion found
res_sim = simulateResiduals(fit_, plot = F)
plot(res_sim) ##The residuals are normally distributed with homogenous variance
plotQQunif(res_sim) ##Plot the qqplot
testResiduals(res_sim) ##This generate all plot and values of all diagnosis test
plotResiduals(res_sim)
plot(res_sim)

##Type 2 Anova
Anova(fit_, type="II") #Only the fixed factor week was significant

emmeans(fit_, specs = pairwise ~ Treatment,
        adjust = "bonferroni", type = "response") 
#Pairwise comparison of the treatment

m_c = emmeans(fit_, specs = pairwise ~ Week|Treatment,
              adjust = "bonferroni", type = "response")

m_c$contrasts %>%
  rbind() 

write.xlsx(m_c$contrasts %>%
             rbind(),"All_code_data_visualization/95_pesticide/cpue.xlsx" )




#emmeans(fit_, ~Treatment * Week, component = "response")


##The pairwise comparison within specific treatment across weeks 

#########Biomass flux#####
#fit inital glmm model

f1 = glmmTMB(mass_flux ~ Treatment + Week +
               (1 | Flume),
             data = abu_bio,
             family = Gamma(link = "log"))
#extract model residuals
res_bio <- residuals(f1)
#calculate and plot the acf of residuals
auto_cor_bio <- acf(res_bio, lag.max=40, plot=FALSE)
plot(auto_cor_bio, main="ACF of Residuals", xlab="Lag", ylab="ACF")
#The model residual shows serial autocorrelation

# Perform Ljung-Box test on residuals to confirm when necessary
#Box.test(residuals, lag = 20, type = "Ljung-Box")

fit_1 = glmmTMB(mass_flux ~ Treatment + Week +
                  (1 | Flume) + ar1(Week+0|Flume),
                data = abu_bio,
                family = Gamma(link = "log")) ##fit the final model with Gamma log link function


summary(fit_1)
diagnose(fit_1)
testDispersion(fit_1)
testUniformity(fit_1)
res_sim = simulateResiduals(fit_1, plot = F)
plot(res_sim)
plotQQunif(res_sim)
plotResiduals(res_sim)
plot(res_sim)

##Type 2 Anova 
Anova(fit_1, type="II")

##Pairwise comparison because the factor #Week is significant

emmeans(fit_1, specs = pairwise ~ Treatment,
        adjust = "bonferroni", type = "response")

m_f = emmeans(fit_1, specs = pairwise ~ Week|Treatment,
               adjust = "bonferroni", type = "response")

m_f$contrasts %>%
  rbind() 

write.xlsx(m_f$contrasts %>%
             rbind(),"All_code_data_visualization/95_pesticide/mass.xlsx" )
########Spider count model#######

fit_sp = glmmTMB(count ~  Treatment * Family + (1|Flume), 
                 data = spider_dat, family = poisson(link = "log"))
#The interaction terms were significant
summary(fit_sp)
diagnose(fit_sp)
testDispersion(fit_sp)
res_sim = simulateResiduals(fit_sp, plot = F)
plot(res_sim)
plotQQunif(res_sim)
plotResiduals(res_sim)
plot(res_sim)


#Fit type 2 Anova
Anova(fit_sp, type="II") #Significant fixed effects and their interaction 

##We therefore conduction a pairwise comparison
#########Compute pairwise comparison##############
sp_t = emmeans(fit_sp,specs = pairwise ~ Treatment,
               adjust = "bonferroni", type = "response") 
sp_t$contrasts %>%
  rbind() 
#Treatment level comparison, there was no significant difference between control and low-flow

sp_e = emmeans(fit_sp, specs = pairwise ~ Treatment|Family,
                     adjust = "bonferroni", type = "response")
sp_e$contrasts %>%
  rbind() 
#Tetragnatha shows a significant difference between Control and Low flow, with the Control condition 
#showing a much higher response (ratio = 2.577, p = 0.0158). with control2.58 times higher than under low-flow. 
#Lycosidae shows No significant difference between Control and Low flow treatments



# print(e_results)
# 
# pairwise_comparisons <- pairs(e_results)
# summary(pairwise_comparisons)
# result = print(pairwise_comparisons)




#####Physical chemical properties measured by Verena and Gemma########
# phch_1 = openxlsx::read.xlsx("../Verena_Schreiner/RSM_2021_phch.xlsx", 
#                              sheet = "macroinv", detectDates = TRUE) |> 
#   filter(!Flume %in% c(3,8,12,15) ) 
# 
# saveRDS(phch_1,"Data/RSM_2021_phch.rds" )

phch_1 = readRDS("All_code_data_visualization/Data/RSM_2021_phch.rds")
##I removed flume 3,8,12, and 15 which was not included in low-flow experiment

# wtw_1 = readxl::read_excel("../Verena_Schreiner/RSM_2021_phch.xlsx", 
#                            sheet = "WTW", 
#                            col_types = c("date", "text", 
#                                          "numeric","numeric",
#                                          "numeric","numeric")) |>  
#   mutate(month = format(date_time, "%m")) |> 
#   filter(!month %in% c("05","09")) 
# 
# saveRDS(wtw_1,"Data/wtw.rds")

wtw_1 = readRDS("All_code_data_visualization/Data/wtw.rds")

##I excluded the month of May and September data,Low-flow was June to July



avg_wtw = wtw_1 |> ##Below, I rename the row to correct overlapping name for analysis
  #dplyr::mutate(across(`sampling point`, as.factor)) |> 
  mutate(`sampling point` = case_when(
    `sampling point` == "Flume01_intlet_I" ~ "Flume01_inlet",
    `sampling point` == "Flume04_intlet_I" ~ "Flume04_inlet",
    `sampling point` == "Flume07_intlet_I" ~ "Flume07_inlet",
    `sampling point` == "Flume10_intlet_I" ~ "Flume10_inlet",
    `sampling point` == "Flume13_intlet_I" ~ "Flume13_inlet",
    `sampling point` == "Flume16_intlet_I" ~ "Flume16_inlet",
    `sampling point` == "Flume07_inlet_II" ~ "Flume7_inlet",
    `sampling point` == "Flume07_outlet_I" ~ "Flume07_outlet",
    `sampling point` == "Flume07_outlet_II" ~ "Flume07_outlet",
    `sampling point` == "Flume13_outlet_I" ~ "Flume13_outlet",
    `sampling point` == "Flume13_inlet_I" ~ "Flume13_inlet",
    `sampling point` == "Flume13_intlet_II" ~ "Flume13_inlet",
    `sampling point` == "Flume13_outlet_II" ~ "Flume13_outlet",
    `sampling point` == "Flume01_inlet_I" ~ "Flume01_inlet",
    `sampling point` == "Flume04_inlet_I" ~ "Flume04_inlet",
    `sampling point` == "Flume07_inlet_I" ~ "Flume07_inlet",
    `sampling point` == "Flume10_inlet_I" ~ "Flume10_inlet",
    `sampling point` == "Flume13_inlet_I" ~ "Flume13_inlet",
    `sampling point` == "Flume16_inlet_I" ~ "Flume16_inlet",
    `sampling point` == "Flume01_outlet_I" ~ "Flume01_outlet",
    `sampling point` == "Flume04_outlet_I" ~ "Flume04_outlet",
    `sampling point` == "Flume10_outlet_I" ~ "Flume10_outlet",
    `sampling point` == "Flume16_outlet_I" ~ "Flume16_outlet",
    
    TRUE ~ `sampling point`
  )) |>
  mutate(
    Treatment = case_when(
      `sampling point` == "Flume01_inlet" ~ "Control",
      `sampling point` == "Flume01_outlet" ~ "Control",
      `sampling point` == "Flume04_inlet" ~ "Low flow",
      `sampling point` == "Flume04_outlet" ~ "Low flow",
      `sampling point` == "Flume07_inlet" ~ "Control",
      `sampling point` == "Flume07_outlet" ~ "Control",
      `sampling point` == "Flume10_inlet" ~ "Low flow",
      `sampling point` == "Flume10_outlet" ~ "Low flow",
      `sampling point` == "Flume13_inlet" ~ "Control",
      `sampling point` == "Flume13_outlet" ~ "Control",
      `sampling point` == "Flume16_inlet" ~ "Low flow",
      `sampling point` == "Flume16_outlet" ~ "Low flow",
      TRUE ~ `sampling point`
    )
  ) |>  
  mutate(
    location = case_when(
      `sampling point` == "Flume01_inlet" ~ "inlet",
      `sampling point` == "Flume01_outlet" ~ "outlet",
      `sampling point` == "Flume04_inlet" ~ "inlet",
      `sampling point` == "Flume04_outlet" ~ "outlet",
      `sampling point` == "Flume07_inlet" ~ "inlet",
      `sampling point` == "Flume07_outlet" ~ "outlet",
      `sampling point` == "Flume10_inlet" ~ "inlet",
      `sampling point` == "Flume10_outlet" ~ "outlet",
      `sampling point` == "Flume13_inlet" ~ "inlet",
      `sampling point` == "Flume13_outlet" ~ "outlet",
      `sampling point` == "Flume16_inlet" ~ "inlet",
      `sampling point` == "Flume16_outlet" ~ "outlet",
      TRUE ~ `sampling point`
    )
  ) |> 
  dplyr::group_by(Treatment, location) |>
  #mutate_if(is.numeric, round, digits = 3) |> 
  rstatix::get_summary_stats(type = "full") ##mean_se for mean and standard error

openxlsx::write.xlsx(avg_wtw,"All_code_data_visualization/95_pesticide/avg_wtw.xlsx" )


sum_phch = phch_1 |>
  filter(phase %in% "treatment") |> ##I use treatment phase as it was during the experimental phase
  dplyr::mutate(across(Flume, factor)) |>
  group_by(location, treatment) |> 
  rstatix::get_summary_stats(type = "full") |>  #estimation of the max velocity and depth in inlet and outlet stretch
  mutate(treatment = case_when(
    treatment == "C" ~ "Control",
    treatment == "D" ~ "Low-flow",
    TRUE ~ treatment  # Keeps other values unchanged
  )) |> 
  arrange(treatment, location)
openxlsx::write.xlsx(sum_phch,"All_code_data_visualization/95_pesticide/avg_phch.xlsx" )

####For estimation of maximum flow velocity in the stream before the start of the experiment####
phch = phch_1 |>
  filter(phase %in% "colonization") |>
  dplyr::mutate(across(c(Flume), as.factor)) |>
  dplyr::group_by(location,treatment) |>
  rstatix::get_summary_stats(type = "full")


##write.table(phch, file = "phch.txt", row.names = FALSE, sep = "\t", quote = FALSE)





