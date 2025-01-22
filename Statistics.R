library(tidyverse)
library(glmmTMB)  
library(bbmle)
library(openxlsx)
library(lme4)
library(DHARMa)
library(car)
library(multcomp)
library(stargazer)
library(broom)
library(MuMIn)
library(car)
library(performance)
library(emmeans)



packageVersion("tidyverse")

citation("emmeans")

passive_stat = emerg_dat |> 
  #filter(Sampling %in% "Passive") |> 
  dplyr::mutate(Days = 7) |> 
  dplyr::mutate(Area = 1) |> 
  dplyr::mutate(CPUE = Calculated.No..Ind/(Area * Days)) |> 
  dplyr::mutate(mass_flux = Emergence.mass..mg./(Area * Days)) |> 
  mutate(Treatment = if_else(Treatment == "Low-flow", "Low flow", Treatment)) |> 
  dplyr::mutate_at(vars(c(Flume,Sampling,Position,Treatment,
                          Size.group, Experiment.Week)), as.factor) |>  
  rename(Sampling_week = Experiment.Week) |> 
  rename(mass = Emergence.mass..mg.)
  #mutate(sub_tag = factor(Position, labels = c("I: inlet", "II: middle", "III: outlet")))
  #dplyr::mutate(Week = Sampling_week)
  #dplyr::mutate(treat_week = paste(Sampling_week, Treatment, sep = "-"))

str(passive_stat)

# Adding blocks to the treatment ------------------------------------------
# passive_stat
# Block 1 (1,2)
# Block 2 (4,5)
# Block 3 (6,7)
# Block 4 (9,10)
# Block 5 (11,14)
# Block 6 (13,16)

passive_block = passive_stat |>
  mutate(
    Block = case_when(
      Flume == "1" ~ "1",
      Flume == "2" ~ "1",
      Flume == "5" ~ "2",
      Flume == "4" ~ "2",
      Flume == "7" ~ "3",
      Flume == "6" ~ "3",
      Flume == "9" ~ "4",
      Flume == "10" ~ "4",
      Flume == "11" ~ "5",
      Flume == "14" ~ "5",
      Flume == "13" ~ "6",
      Flume == "16" ~ "6",
      TRUE ~ Flume
    ))

str(passive_block)
passive_block = passive_block |> 
  dplyr::mutate(across(Block, as.factor))
# Abundance flux ----------------------------------------------------------
options(na.action = "na.fail")

wk2 = passive_stat |> 
  filter(Sampling_week %in% 2)

fit_1 = glmmTMB(CPUE ~ Treatment * Position + 
                  (1 | Flume) + ar1(Week + 0 | Flume), 
                data = wk2, 
                family = Gamma(link = "log"))

it_1 = lm(CPUE ~
                 Treatment * Position,
                  data = wk1,
               #family = Gamma(link = "log")
)
summary(fit_1)
Anova(fit_1, type="II")

fit_1 = glmmTMB(CPUE ~ Treatment * Position + Week +
               (1|Sampling_week) + (1|Flume),
               family = gaussian(link = "log"),
               data = passive_stat)



summary(fit_1)
Anova(fit_1, type = 2)


#######################
e_results <- emmeans(fit_1, ~ Treatment * Position)
pairwise_comparisons <- pairs(e_results)
summary(pairwise_comparisons)
################

str(passive_stat)

p_up = passive_stat |> 
  dplyr::filter(Position %in% "inlet")

p_mid = passive_stat |> 
  dplyr::filter(Position %in% "middle")

p_out = passive_stat |> 
  dplyr::filter(Position %in% "outlet")

passive_stat |> 
  ggplot(aes(Week, mass_flux)) +
  geom_line()

str(passive_stat)

summary(passive_stat$mass_flux)

####Emergence data for both mass flux and emergence rate
p_sample = passive_sample |> 
  filter(Sampling %in% "Passive")

a_sample = passive_sample |> 
  filter(Sampling %in% "Active")


f = glmmTMB(CPUE ~ Treatment * Week +
              (1 | Flume),
            data = a_sample,
            family = Gamma(link = "log"))

residuals <- residuals(f)
autocorrelation <- acf(residuals, lag.max=40, plot=FALSE)
plot(autocorrelation, main="ACF of Residuals", xlab="Lag", ylab="ACF")

Box.test(residuals, lag = 20, type = "Ljung-Box")


fit_ = glmmTMB(CPUE ~ Treatment * Week +
                  (1 | Flume) + ar1(Week+0|Flume),
                data = p_sample,
              family = Gamma(link = "log"))
summary(fit_)

diagnose(fit_)
testDispersion(fit_)
res_sim = simulateResiduals(fit_, plot = F)
plot(res_sim)
plotQQunif(res_sim)
plotResiduals(res_sim)
plot(res_sim)


res2 <- recalculateResiduals(res_sim, group=passive_stat$Week, rotation=estimated) 

plot(res2)


Anova(fit_, type="II")

#########Compute pairwise comparison##############
em_results <- emmeans(fit_, ~ Treatment * Sampling_week,
                      adjust = "bonferroni")
pairwise_comparisons <- pairs(em_results)
summary(pairwise_comparisons)
print(pairwise_comparisons)


#########

f1 = glmmTMB(mass_flux ~ Treatment + Sampling_week +
              (1 | Flume),
            data = passive_stat,
            family = Gamma(link = "log"))

residuals <- residuals(f1)
autocorrelation <- acf(residuals, lag.max=40, plot=FALSE)
plot(autocorrelation, main="ACF of Residuals", xlab="Lag", ylab="ACF")

Box.test(residuals, lag = 20, type = "Ljung-Box")

fit_1 = glmmTMB(mass_flux ~ Treatment * Sampling_week +
                  (1 | Flume) + ar1(Sampling_week+0|Flume),
                data = passive_stat,
                family = Gamma(link = "log"))

m_sum = summary(fit_1)
fixed_effects <- as.data.frame(m_sum$coefficients$cond)
print(fixed_effects)
write.xlsx(fixed_effects, "./plot/mass.xlsx")


diagnose(fit_1)
testDispersion(fit_1)
res_sim = simulateResiduals(fit_1, plot = F)
plot(res_sim)
plotQQunif(res_sim)
plotResiduals(res_sim)
plot(res_sim)



Anova(fit_1, type="II")

#########Compute pairwise comparison##############
e_results <- emmeans(fit_1, ~ Treatment * Sampling_week,
                     adjust = "bonferroni")
pairwise_comparisons <- pairs(e_results)
summary(pairwise_comparisons)
result = print(pairwise_comparisons)

write.xlsx(result, "./plot/emm_mass.xlsx")

str(spider_dat)
#########Spider###
tet_sp = spider_dat |> 
  filter(Family %in% "Tetragnatha") |> 
  #dplyr::select(-c(Family, Treat_fam)) |> 
  dplyr::mutate(across(c(Flume, Treatment, Family), as.factor))

lyc_sp = spider_dat |> 
  filter(Family %in% "lycosidae") |> 
  #dplyr::select(-c(Family, Treat_fam)) |> 
  dplyr::mutate(across(c(Flume, Treatment, Family), as.factor))

######Tetragnathidae####
library(nlme)
t = nlme::gls(count ~ Treatment, tet_sp,
              weights = varIdent(form = ~1|Flume))

summary(t)

kruskal.test(count ~ Treatment, data = tet_sp)
wilcox.test(count ~ Treatment, data=tet_sp) 

plot(t)
anova(t)

fit_tet = glmmTMB(count ~  Treatment + (1|Flume), 
                  data = spider_dat, family = poisson(link = "log"))

summary(fit_tet)
diagnose(fit_tet)
testDispersion(fit_tet)
res_sim = simulateResiduals(fit_tet, plot = F)
plot(res_sim)
plotQQunif(res_sim)
plotResiduals(res_sim)
plot(res_sim)



Anova(fit_1, type="II")


######Lycosidae#####

fit_lyc = glmer(count ~  Treatment + (1|Flume), 
                  data = lyc_sp,
                family = poisson)

summary(fit_lyc)
diagnose(fit_lyc)
testDispersion(fit_lyc)
res_sim = simulateResiduals(fit_lyc, plot = F)
plot(res_sim)
plotQQunif(res_sim)
plotResiduals(res_sim)
plot(res_sim)



Anova(fit_1, type="II")

#############################





ft1 = dredge(fit_1, beta = c("sd"), 
                   evaluate = T, rank = AICc, fixed = T)
print(ft1)

best_model <- get.models(ft1, delta == 0)[[1]]
summary(best_model)

Anova(best_model, type = 2)

summary(glht(best_model))
#########Compute pairwise comparism##############
e_results <- emmeans(best_model, ~ Treatment)
pairwise_comparisons <- pairs(e_results)
summary(pairwise_comparisons)
print(pairwise_comparisons)
# Compute profile likelihood confidence interval -------------------------------------

confint(fit_1, method="profile")

# Perform a likelihood ratio test to assess the significance of predictors
drop1(fit_1, test="Chisq")

#############################


summary(best_model)
diagnose(best_model)
testDispersion(best_model)

res_sim = simulateResiduals(best_model, plot = F)
plot(res_sim)
plotQQunif(res_sim)
plotResiduals(res_sim)
plot(res_sim)
model_performance(fit_1)



# Biomass flux ------------------------------------------------------------

fit_2 = glmmTMB(mass_flux ~ Sampling_week + Treatment * Position +
                  (1|Flume), data = passive_stat,
                family = Gamma(link = "log"))


summary(fit_2)
ft2 = dredge(fit_2, beta = c("sd"), 
                   evaluate = T, rank = AICc, fixed = T)
print(ft2)

best_model2 <- get.models(ft2, delta == 0)[[1]]

Anova(best_model2, type = 2)
summary(glht(best_model2))

# Print the summary of the best model
#summary(best_model2)

summary(best_model2)
diagnose(best_model2)
testDispersion(best_model2)

res_sim2 = simulateResiduals(best_model2, plot = F)
plotQQunif(res_sim2)
plotResiduals(res_sim2)
plot(res_sim2)
r_2 = performance::r2(best_model2)
print(r_2). ##Conditionl and marginal r square
model_performance(best_model2)



g_mod = mgcv::gam(CPUE ~ as.numeric(Sampling_week), data = passive_stat)


plot(CPUE ~ Sampling_week, data = passive_stat)

summary(g_mod)

options(na.action = "na.fail")
# fit_comp <- glht(fit_1, linfct = mcp(Experiment.Week = "Tukey",
#                                      Treatment = "Tukey",
#                                      Position = "Tukey"
#                                      ))
# summary(fit_comp)
#Best fit model
# ft1_model = dredge(fit_1, beta = c("none", "sd", "partial.sd"), 
#                    evaluate = T, rank = AICc)
# print(ft1_model)
# 
# best_model2 <- get.models(ft1_model, delta == 0)[[1]]
# 
# # Print the summary of the best model
# summary(best_model2)
# 
# c_int = confint(best_model2, method = "profile")
# 
# str(c_int)
# 
# # Convert to a data frame for easier manipulation with dplyr
# profile_intervals_df <- as.data.frame(c_int) %>%
#   rownames_to_column(var = "parameter")
# 
# 
# # Identify and transform random effect parameters using dplyr
# transformed_intervals_df <- profile_intervals_df %>%
#   mutate(
#     `2.5 %` = exp(`2.5 %`),
#     `97.5 %` = exp(`97.5 %`)
#   )
# 
# 
# best_model2$sdr$pdHess
# 
# diagnose(best_model2)
# testDispersion(best_model2)
# 
# res_sim = simulateResiduals(best_model2, plot = F)
# plotQQunif(res_sim)
# plotResiduals(res_sim)
# plot(res_sim)
# #######
# diagnose(
#   best_model2,
#   eval_eps = 1e-05,
#   evec_eps = 0.01,
#   big_coef = 10,
#   big_sd_log10 = 3,
#   big_zstat = 5,
#   check_coefs = TRUE,
#   check_zstats = TRUE,
#   check_hessian = TRUE,
#   check_scales = TRUE,
#   explain = TRUE
# )
# 
# 
# ########
# confint(best_model)
# exp_coef = exp(coef(summary(fit_2))$cond[, "Estimate"])
# print(exp_coef)
# # Multiple comparism ------------------------------------------------------
# 
# fit_comp <- glht(best_model2, linfct = mcp(Experiment.Week = "Tukey"))
# summary(fit_comp)
# 
# sum_fit = summary(best_model2)
# 
# # Extract the fixed effects from the model summary
# fixed_effects2 <- sum_fit$coefficients$cond
# 
# # Convert fixed effects to a data frame if necessary
# fixed_effect2_df <- as.data.frame(fixed_effects)
# 
# write.xlsx(fixed_effect2_df, "CPUE.xlsx")
# 
# broom::glance(summary(best_model2))


# Biomass data model ------------------------------------------------------



# fit_comp <- glht(fit_2, linfct = mcp(Experiment.Week = "Tukey",
#                                      Treatment = "Tukey",
#                                      Position = "Tukey"))
# summary(fit_comp)


############
# ft2_model = dredge(fit_2, beta = c("none", "sd", "partial.sd"), 
#                    evaluate = T, rank = AICc)
# print(ft2_model)
# 
# # Best fit model ----------------------------------------------------------
# # We use the get.model function in mumin package to select the best fit model 
# #following the dredge function 
# #the dredge function based on the lowest AICc value.
# best_model <- get.models(ft2_model, delta == 0)[[1]]
# 
# # Print the summary of the best model
# summary(best_model)
# 
# #####################
# diagnose(best_model)
# testDispersion(best_model)
# 
# res_sim = simulateResiduals(best_model, plot = TRUE)
# 
# plotQQunif(res_sim)
# plotResiduals(res_sim)
# plot(res_sim)

#####Multiple comparison
#fit_comp <- glht(fit_1, linfct = mcp(Position = "Tukey"))
#summary(fit_comp)
###################



#vif(fit_1)
ft1_model = dredge(fit_1, beta = c("none", "sd", "partial.sd"), 
                   evaluate = T, rank = AICc, fixed = T)
print(ft1_model)

md_avg = model.avg(ft1_model, fit = T, subset = delta < 4)
summary(md_avg)



# Tetragnatha ----------------------------------------------------------

tet_ = tet_ |>
  mutate(
    rep = case_when(
      Flume == "1" ~ "Flume1",
      Flume == "2" ~ "Flume2",
      Flume == "5" ~ "Flume3",
      Flume == "4" ~ "Flume4",
      Flume == "7" ~ "Flume7",
      Flume == "6" ~ "Flume6",
      Flume == "9" ~ "Flume9",
      Flume == "10" ~ "Flume10",
      Flume == "11" ~ "Flume11",
      Flume == "14" ~ "Flume14",
      Flume == "13" ~ "Flume13",
      Flume == "16" ~ "Flume16",
      TRUE ~ Flume
    ))

tet_ <- tet_count |>
  mutate(
    Flume = case_when(
      Flume == "1" ~ "Flume1",
      Flume == "2" ~ "Flume2",
      Flume == "5" ~ "Flume5",
      Flume == "4" ~ "Flume4",
      Flume == "7" ~ "Flume7",
      Flume == "6" ~ "Flume6",
      Flume == "9" ~ "Flume9",
      Flume == "10" ~ "Flume10",
      Flume == "11" ~ "Flume11",
      Flume == "14" ~ "Flume14",
      Flume == "13" ~ "Flume13",
      Flume == "16" ~ "Flume16",
      TRUE ~ as.character(Flume)  # Convert Flume to character
    )
  )


lyco_ = lyco_count |> 
  rename(Family = family)

tet_ = tet_ |> 
  rename(Family = family)

spider_stat = rbind(tet_,lyco_)

write.xlsx(spider_stat, "spider.xlsx")
write.xlsx(wk5, "abu_bio.xlsx")

sp_em = openxlsx::read.xlsx("spider_emergence.xlsx", sheet = "sp_em")

spider_stat$Flume
spider_stat = spider_stat |> 
  dplyr::mutate(across(c(Treatment, family, Flume), as.factor)) 
#dplyr::mutate(Treat_fam = paste(Treatment, Family, sep = "-"))
sp_em = sp_em |> 
  dplyr::mutate(across(c(Treatment, Family, Flume), as.factor)) |> 

sp_tet = sp_em |> 
  dplyr::filter(!Family %in% c("Lycosidae", "biomass_flux")) |> 
  droplevels() |> 
  mutate(Family = relevel(Family, ref = "Tetragnathidae")) 

sp_lyc = sp_em |> 
  dplyr::filter(!Family %in% c("Tetragnathidae", "biomass_flux")) |> 
  droplevels() |> 
  mutate(Family = relevel(Family, ref = "Lycosidae")) |> 
  
  
summary(sp_lyc)

lyc = spider_stat |> 
  dplyr::filter(Family %in% "Lycosidae") |> 
  rstatix::group_by(Treatment) |> 
  rstatix::get_summary_stats(count, type = "mean_se")

tet = spider_stat |> 
  dplyr::filter(Family %in% "Tetragnathidae") |> 
  rstatix::group_by(Treatment) |> 
  rstatix::get_summary_stats(count, type = "mean_se")


fit_sp = glmmTMB(count ~  Treatment * Family + (1|Flume), 
                  data = spider_dat, family = poisson(link = "log"))
summary(fit_sp)
diagnose(fit_sp)
testDispersion(fit_sp)
simulateResiduals(fit_sp, plot = T)
Anova(fit_sp, type = 2)


#########Pairwise Comparison########
e_tet <- emmeans(fit_sp, ~ Treatment * Family,
                 adjust = "bonferroni")
pairwise_comparisons <- pairs(e_tet)
summary(pairwise_comparisons)
tet_res = print(pairwise_comparisons)

write.xlsx(tet_res, "./plot/spider_stat.xlsx")

##################

em_bio = passive_sample |> 
  dplyr::filter(Experiment.Week %in% 5)

dev.off()
# plot(all_dat$Tetragnathidae,all_dat$emergence)
# abline(lm(emergence ~ Tetragnathidae, data = all_dat), col = "red", lwd = 2)
# cor_coeff <- cor(all_dat$Tetragnathidae, all_dat$emergence, use = "complete.obs")
# print(cor_coeff)

fit_em = glmmTMB(Tetragnathidae ~  Treatment * emergence + (1|Flume), 
                  data = all_dat, family = poisson(link = "log"))

summary(fit_em)
diagnose(fit_em)
testDispersion(fit_em)
simulateResiduals(fit_em, plot = T)
Anova(fit_em, type = 2)

summary(all_dat$emergence)
e_res = emmeans(fit_em, ~ Treatment * emergence, at = list(emergence = c(5, 15, 25, 35, 50, 70)))
pairwise_comparisons <- contrast(e_res, method = "pairwise")

# View results with p-values
summary(pairwise_comparisons, infer = c(TRUE, TRUE))  # infer = TRUE provides CIs and p-values

#########Pairwise Comparison########
e_tet <- emmeans(fit_em,~ Treatment * emergence,
                 adjust = "bonferroni")
pairwise_comparisons <- pairs(e_tet)
summary(pairwise_comparisons)
print(pairwise_comparisons)
##################

fit_bio = glmmTMB(mass_flux ~  Treatment + (1|Flume), 
                  data = em_bio, family = poisson(link = "log"))
summary(fit_bio)
diagnose(fit_bio)
testDispersion(fit_bio)
simulateResiduals(fit_bio, plot = T)
Anova(fit_bio, type = 2)

#########Pairwise Comparison########
e_tet <- emmeans(fit_sp, ~ Treatment * Family,
                 adjust = "bonferroni")
pairwise_comparisons <- pairs(e_tet)
summary(pairwise_comparisons)
print(pairwise_comparisons)

# Print the summary of the best model
summary(fit_tet)

diagnose(fit_tet)
testDispersion(fit_tet)
simulateResiduals(fit_tet, plot = T)

Anova(fit_tet, type = 2)

#########Pairwise Comparison########
e_tet <- emmeans(fit_tet, ~ Treatment * Family,
                 adjust = "bonferroni")
pairwise_comparisons <- pairs(e_tet)
summary(pairwise_comparisons)
print(pairwise_comparisons)

########################################
#####Correlation#####
bio <- sp_em |>
  dplyr::filter(Family == "biomass_flux") |>
  dplyr::group_by(Flume,Treatment) |>
  dplyr::summarise(biomass = mean(count, na.rm = TRUE)) |>
  dplyr::arrange(Treatment,desc(Flume)) |> 
  ungroup() |> 
  dplyr::select(biomass)

emerg = sp_em |> 
  dplyr::filter(Family == "emergence_rate") |>
  dplyr::group_by(Flume,Treatment) |>
  dplyr::summarise(emergence = mean(count, na.rm = TRUE)) |> 
  dplyr::arrange(Treatment,desc(Flume)) |> 
  ungroup() |> 
  dplyr::select(emergence)

lyco = sp_em |> 
  dplyr::filter(Family %in% "Lycosidae") |>
  dplyr::arrange(Treatment,desc(Flume)) |> 
  dplyr::select(count) |> 
  rename(Lycosidae = count)

Tet = sp_em |> 
  dplyr::filter(Family %in% "Tetragnathidae") |> 
  dplyr::arrange(Treatment,desc(Flume)) |> 
  dplyr::select(Flume, Treatment,count) |> 
  rename(Tetragnathidae = count)

all_dat = cbind(bio,emerg,lyco, Tet)

ggplot(all_dat, aes(x = biomass, y = Tetragnathidae)) +
  geom_point() +
  labs(x = "Biomass", y = "Lycosidae", title = "Biomass vs Lycosidae")


plot(all_dat$biomass, all_dat$emergence, all_dat$Lycosidae, all_dat$Tetragnathidae,
     xlab = "Emergence Biomass",
     ylab = "Spider Count",
     main = "Scatter Plot of Emergence Biomass vs Spider Count")

correlation <- cor(all_dat$emergence, all_dat$Tetragnathidae)
correlation_test <- cor.test(all_dat$emergence, all_dat$Tetragnathidae)

# Print results
print(correlation)        # Correlation coefficient
print(correlation_test)    # Correlation coefficient + p-value + confidence interval

##########


ft_tet = dredge(fit_tet,beta = "sd",
                   evaluate = T, rank = AICc)
print(ft_tet)


# select best fit model ---------------------------------------------------

best_tet <- get.models(ft_tet, delta == 0)[[1]]

Anova(best_tet, type = 2)



# Print the summary of the best model
summary(best_tet)

diagnose(best_tet)
testDispersion(best_tet)

res_sim = simulateResiduals(best_tet)
plot(res_sim)
model_performance(best_tet)

#########Pairwise Comparison########
e_tet <- emmeans(best_tet, ~ Treatment * Family)
pairwise_comparisons <- pairs(e_tet)
summary(pairwise_comparisons)
print(pairwise_comparisons)

########################################
######Sediment concentration#####

# sed_stat = dat_sed |> 
#   mutate(Flume = str_extract(Sample, "F\\d+")) |> 
#   mutate(across(Flume, as.factor))
# 
# summary(sed_stat$Concentration_µg_kg)
# 
# str(sed_stat)
# 
f = lmer(Concentration ~ Treatment * Week +
              (1 + Week | Flume), data = sed_pest)

summary(f)
hist(residuals(f), main = "Histogram of Residuals", xlab = "Residuals")
qqnorm(residuals(f))
qqline(residuals(f), col = "red")
shapiro.test(residuals(f))
# library(nlme)

# sed_stat$Concentration <- log(sed_stat$Concentration_µg_kg+1)
# fit_sed <- lmer(Concentration ~ Treatment + Week +
#                   (1|Flume),
#                data = sed_stat)
# summary(fit_sed)
library(mgcv)
model <- gam(Concentration ~ Week * Treatment,
             data = sed_pest, family = Gamma())

# Residual diagnostics
mgcv::gam.check(model_gamm$gam)

model_gamm <- gamm(Concentration ~ Week * Treatment,
                   random = list(Flume = ~1),
                   data = sed_pest, family = Gamma())
# Summary of the fixed effects
summary(model_gamm$gam)

# Summary of the random effects
summary(model_gamm$lme)

plot(model_gamm$lme, shade = T)

model_gammpql <-MASS::gammpql(
  fixed = Concentration ~ s(Week) + Treatment,
  random = ~1 | Flume,
  family = Gamma,
  data = sed_pest
)




residuals <- residuals(model_gamm$lme)
autocorrelation <- acf(residuals, lag.max=40, plot=T)
plot(autocorrelation, main="ACF of Residuals", xlab="Lag", ylab="ACF")

Box.test(residuals, lag = 20, type = "Ljung-Box")

model_bam <- bam(Concentration ~ Week + Treatment + s(Flume, bs = "re"),
                 data = sed_pest, family = gaussian())

############################################
sed_pest$Concentration <- log(sed_pest$Concentration_µg_kg+1)

str(sed_pest)

str(sediment_final)

fp = glmmTMB(Concentration ~  Treatment * Week
             + (1|Flume),
             data = sed_pest, family = Gamma(link = "log"))

summary(fp)
diagnose(fp)
simulateResiduals(fp, plot = T)

residuals <- residuals(fp)
autocorrelation <- acf(residuals, lag.max=40, plot=T)
plot(autocorrelation, main="ACF of Residuals", xlab="Lag", ylab="ACF")

Box.test(residuals, lag = 20, type = "Ljung-Box")

ab = sed_pest |> 
  dplyr::filter(!Class %in% "Herbicide")

hist(sed_pest$Concentration_µg_kg)
str(sediment_final)


fit_sed = glmmTMB(Concentration ~   Treatment * Week
                  + (1|Flume) + ar1(Week + 0 | Flume),
                  data = sed_pest,
                  family = Gamma(link = "log"))


summary(fit_sed)
diagnose(fit_sed)
testDispersion(fit_sed)
test
res_sim = simulateResiduals(fit_sed, plot = T)
#plotResiduals(res_sim, sed_stat$Pesticide_Class)
#testOutliers(res_sim)
#testQuantiles(res_sim)
print(res_sim)
plot(res_sim)
plotQQunif(res_sim)
plotResiduals(res_sim)
Anova(fit_sed, type = 2)

# emm <- emmeans(fit_sed, pairwise ~ Treatment | Week, type = "response")  # Pairwise comparisons
# summary(emm)

#########Pairwise Comparison########
e_sed <- emmeans(fit_sed, ~ Treatment + Week,
                 adjust = "bonferroni", detailed = T)
pairwise_comparisons <- pairs(e_sed)
summary(pairwise_comparisons, infer = TRUE)
p_result = print(pairwise_comparisons)



#+ ar1(Week + 0 | Flume)
# fit_sed = glm(Concentratinbinom12()# fit_sed = glm(Concentration_µg_kg ~  Treatment * Week,
#                   data = sed_stat, family = Gamma(link = "inverse"))
# 
# 
# residuals_glmm <- residuals(fit_sed, type = "deviance")
# qqnorm(residuals_glmm)
# qqline(residuals_glmm, col = "red")
# 
# shapiro.test(residuals_glmm)
# ks.test(residuals_glmm, "pnorm", mean = mean(residuals_glmm), sd = sd(residuals_glmm))
# hist(residuals_glmm)
# hist(sed_stat$Concentration_µg_kg)

residuals <- residuals(fit_sed)
shapiro.test(residuals)

summary(fit_sed)
diagnose(fit_sed)
testDispersion(fit_sed)
res_sim = simulateResiduals(fit_sed, plot = T)
#plotResiduals(res_sim, sed_stat$Pesticide_Class)
#testOutliers(res_sim)
#testQuantiles(res_sim)
print(res_sim)
plot(res_sim)
plotQQunif(res_sim)
plotResiduals(res_sim)
Anova(fit_sed, type = 2)

#########Pairwise Comparison########
e_sed <- emmeans(fit_sed, ~ Treatment * Week,
                 adjust = "bonferroni")
pairwise_comparisons <- pairs(e_sed)
summary(pairwise_comparisons, infer = TRUE)
p_result = print(pairwise_comparisons)


# # Extract scaled residuals
# residuals <- res_sim$scaledResiduals
# 
# # Plot density of residuals
# plot(density(residuals),
#      main = "Density of Residuals",
#      xlab = "Residuals",
#      ylab = "Density",
#      col = "blue",
#      lwd = 2)
# 
# # Add a normal distribution curve
# curve(dnorm(x, mean = mean(residuals), sd = sd(residuals)), 
#       col = "red", 
#       lwd = 2, 
#       add = TRUE)
# 
# # Optional: Add a uniform distribution curve (for DHARMa residuals)
# curve(dunif(x, min = 0, max = 1), 
#       col = "green", 
#       lwd = 2, 
#       add = TRUE)
# 
# # Add legend
# legend("topright", legend = c("Residual Density", "Normal Curve", "Uniform Curve"),
#        col = c("blue", "red", "green"), lwd = 2)
# 
# 
# # Calculate residual variance by class
# class_var <- sed_stat %>%
#   mutate(residual = res_sim$scaledResiduals) %>%
#   group_by(Pesticide_Class) %>%
#   summarise(var_residual = var(residual))
# 
# # Compute weights inversely proportional to variance
# class_var <- class_var %>%
#   mutate(weight = 1 / var_residual)
# 
# # Join weights back to the dataset
# sed_stat <- sed_stat %>%
#   left_join(class_var, by = "Pesticide_Class")
# 
# sed_stat$weight <- sed_stat$weight / max(sed_stat$weight)
# 
# model <- glmmTMB(Concentration ~ Treatment * Week * Pesticide_Class
#                  + (1 | Flume) , 
#                  data = sed_stat, 
#                  family = Gamma(link = "log"))
# 
# summary(model)
# diagnose(model)
# testDispersion(model)
# res_sim = simulateResiduals(model, plot = T)
# plotResiduals(res_sim, sed_stat$Pesticide_Class)
# #testOutliers(res_sim)
# #testQuantiles(res_sim)
# print(res_sim)
# plot(res_sim)
# plotQQunif(res_sim)
# plotResiduals(res_sim)
# Anova(fit_sed, type = 2)
#ft_tet = dredge(fit_sed,beta = "sd",
#                 evaluate = T, rank = AICc)
# print(ft_tet)
# 
# 
# # select best fit model ---------------------------------------------------
# 
#best_tet <- get.models(ft_tet, delta == 0)[[1]]
# 
# Anova(best_tet, type = 2)

# Print the summary of the best model




#########Pairwise Comparison########
e_sed <- emmeans(fit_sed, ~ Treatment * Week,
                 adjust = "bonferroni")
pairwise_comparisons <- pairs(e_sed)
summary(pairwise_comparisons, infer = TRUE)
p_result = print(pairwise_comparisons)

write.xlsx(p_result, "./plot/pesticide.xlsx")

#citation("performance")

#Marginal R² (R²_m = 0.385): This indicates that 39% of the variance in the response variable is explained by the fixed effects alone.
#Conditional R² (R²_c = 0.884): This indicates that 88% of the variance is explained by both fixed and random effects combined.

##############################
exp_coef = exp(coef(summary(fit_2))$cond[, "Estimate"])
print(exp_coef)

citation ("glmmTMB")
citation ("DHARMa")
citation("MuMIn")
citation("tidyverse")
packageVersion("DHARMa")
citation("performance")
citation("tidyverse")
citation("rstatix")

#Gamma(link = "log")

summary(fit_1)
diagnose(best_model)
testDispersion(best_model)

res_sim = simulateResiduals(best_model, plot = TRUE)
plotQQunif(res_sim)
plotResiduals(res_sim)
plot(res_sim)

fit_comp <- glht(fit_1, linfct = mcp(Position = "Tukey"))
summary(fit_comp)

sum_fit = summary(fit1)

stargazer(sum_fit$coefficients$cond,
          title = "Seasonal variability model",
          type = "text",
          keep.stat = c("all"))

#rikz <- read.delim("./plot/RIKZ.txt", header = TRUE, sep = "\t")
remove(rikz)

