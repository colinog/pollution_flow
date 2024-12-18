library(openxlsx)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(gridExtra)
#source("./theme_pub.R")
library(stargazer)
library(openxlsx)
# Pesticide data (Sediment and Water) -------------------------------------

dat_water = read.csv("data_water.csv")|>
  filter(!Concentration %in% c("<LOQ", NA)) |> 
  mutate(across(Concentration, as.numeric))

head(sediment_pesticide)

str(sediment_pesticide)

dat_sed = read.csv("./Corrected_Alex/Sediment_Pest_Concentration.csv") |> 
  rename(Sample = Sample.ID) |> 
  dplyr::mutate(across(c(Week,Treatment,Treatment_Week,Pesticide_Class), as.factor))

str(dat_water)
water_sum = dat_water %>%
  dplyr::mutate(across(c(Analyte,Type,Sample), as.factor)) |> 
  filter(!Type %in% "Metabolite") |> 
  group_by(Type, Analyte) %>%
  rstatix::get_summary_stats(Concentration, type = "full")
  #dplyr::filter(Type %in% "Fungicide") |> 
  #arrange(desc(n)) |>
  #slice(1:5)

mean(sediment_sum$n)

write.csv(water_sum, "water_con.csv")
write.csv(sediment_sum, "sediment_con.csv")
str(dat_sed)

sediment_sum = dat_sed %>%
  group_by(Treatment, Week, Pesticide_Class) %>%
  rstatix::get_summary_stats(Concentration_µg_kg, type = "full") 
  arrange(desc(n)) |> 
  slice(1:5)

write.xlsx(sediment_sum, "sedinment_class.xlsx")

sediment_sum |> 
  summarise(mean(n))

# Emergence and Spider data -----------------------------------------------

emerg_insect = read.csv("./original_spider_emergence/Emergence_RSM_Mass.csv", 
                        header = T, sep = ";") |> 
  mutate(Diptera = rowSums(across(c(no..Chironomids, No..Other.Nematocera)), na.rm = T)) |> 
  mutate(Diptera = if_else(Diptera == 0, NA_real_, Diptera))

#######Percentage estimation######
emerg_p = emerg_insect |> 
  filter(Sampling %in% "Passive")
str(emerg_p)

emerg_tab = emerg_p |> 
  dplyr::select(-c(1:8,13,19)) |> 
  tidyr::pivot_longer(
    5:10,
    names_to = "Order_family",
    values_to = "Count")

emerg_sum = emerg_tab |> 
  dplyr::mutate(Days = 7) |> 
  dplyr::mutate(Area = 1) |> 
  dplyr::filter(!Count %in% NA) |> 
  dplyr::mutate(CPUE = Count/(Area * Days)) |> 
  group_by(Order_family, Position, Treatment) |> 
  rstatix::get_summary_stats(CPUE, type = "mean_se") 

write.xlsx(emerg_sum, "./plot/emerg_sum_SI.xlsx")
  



sum(sum(emerg_p$no..Chironomids, na.rm = T),
sum(emerg_p$No..Other.Nematocera, na.rm = T),
sum(emerg_p$No..mayflies, na.rm = T),
sum(emerg_p$No..Plecoptera, na.rm = T),
sum(emerg_p$no..Trichoptera, na.rm = T))

sum(emerg_insect$Chironomids, na.rm = T)

(sum(emerg_p$no..Trichoptera, na.rm = T)/12611)*100



#Select the column necessary for analysis
emerg_dat = emerg_insect |> 
  dplyr::select(Date.2021,Flume,Sampling,Position,Treatment,Experiment.Week,Size.group,
         Emergence.mass..mg.,Calculated.No..Ind) 
 

# Convert the wide format to long format using the count value ------------
#Create a new column "family" with only the order or family name
 
no_chiro = emerg_insect |> 
  select(8:18) |> 
  tidyr::pivot_longer(
    7:11,
    names_to = "No_individual",
    values_to = "Count") |> 
   mutate(family = str_remove(No_individual,"no\\..|No.."))

passive_sample = emerg_dat |> 
  filter(Sampling %in% "Passive") |> 
  dplyr::mutate(Days = 7) |> 
  dplyr::mutate(Area = 1) |> 
  dplyr::mutate(CPUE = Calculated.No..Ind/(Area * Days)) |> 
  dplyr::mutate(mass_flux = Emergence.mass..mg./(Area * Days)) |> 
  dplyr::mutate_at(vars(c(Flume,Sampling,Position,Treatment,
                        Size.group)), as.factor) 
  #mutate(Treatment = if_else(Treatment == "Low-flow", "Low flow", Treatment))


passive_sample |> 
  ggplot(aes(Treatment,mass_flux, fill = Treatment))+
  geom_boxplot(position = "dodge2") +
  theme_Publication()

a =passive_sample |> 
  dplyr::select(Treatment, Position,Experiment.Week, mass_flux) |> 
  rstatix::group_by(Treatment, Experiment.Week) |> 
  rstatix::get_summary_stats(type = "full")



pd = position_dodge(width = 0.1)

passive_sample |> 
  dplyr::group_by(Treatment, Experiment.Week) |> 
  rstatix::get_summary_stats(CPUE, type = "full")|> 
  ggplot(aes(Experiment.Week,mean, fill = Treatment, color = Treatment, 
             shape = Treatment, linetype = Treatment))+
  geom_point( position = pd, size = 3, color = "black")+
  geom_line(position = pd, color = "black", linewidth = 0.3)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.1, position = pd, color = "black") +
  scale_fill_viridis_d(option = "plasma")+
  scale_linetype_manual(values = c("solid", "dashed"))+
  labs(y = expression(bold("Emergence rate (ind. m"^-2*" day"^-1*")")), 
       x = "Sampling Week")+
  #facet_wrap(~ Size.group, ncol = 1) 
  #labs(tag = "a)")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.title.y = element_text(face = "bold"),
    axis.title = element_text(size = 15, color = "black",
                              face = "bold"),
    axis.text = element_text(size = 15, face = "bold", color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 15, face = "bold", color = "black"),
    legend.title = element_text(size = 17),
    strip.text = element_text(size = 15, face = "bold"),
    plot.tag = element_text(face = "bold")
  )
  
ggsave("./plot/cpue.png", dpi = 300, width = 15, height = 12, units = "cm")  




combined_plot <- a + b + plot_layout(ncol = 2)

patchwork::wrap_plots(a, b, ncol = 2, widths = 1)

print(combined_plot)
################
passive_sample |> 
  dplyr::group_by(Treatment, Position) |> 
  rstatix::get_summary_stats(CPUE, type = "mean_se")|> 
  ggplot(aes(Treatment,mean, fill = Treatment, color = Treatment, 
             shape = Treatment))+
  geom_point(position = pd, size = 3, color = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.1, position = pd, color = "black") +
  #scale_fill_viridis_d(option = "plasma")+
  facet_wrap(~Position, ncol = 1)+
  #scale_x_discrete(name ="Size class", 
                   #limits=c("Small", "Medium", "Large"))+
  theme_Publication(base_size = 14) +
  labs(y = "Mean spider count")


passive_sample |> 
  dplyr::group_by(Treatment) |> 
  rstatix::get_summary_stats(mass_flux, type = "mean_se")|> 
  ggplot(aes(x = Treatment, y = mean, fill = Treatment, color = Treatment, shape = Treatment)) +
  geom_point(size = 3, color = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1, color = "black") +
  #facet_wrap(~Position, ncol = 1) +
  #ylim(0, 90)+
  #scale_x_discrete(expand = expansion(mult = c(3, 3))) +  # Adjust spacing
  #theme_Publication(base_size = 14) +
  labs(y = "")+
  # theme(
  #   axis.text.x = element_text(angle = 90, hjust = 1),# Rotate x-axis labels if needed
  #   legend.position = "none"
  # )+
  
  
  labs(tag = "a)")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.title.y = element_text(face = "bold"),
    axis.title = element_text(size = 15, color = "black",
                              face = "bold"),
    axis.text = element_text(size = 15, face = "bold", color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 15, face = "bold", color = "black"),
    legend.title = element_text(size = 17),
    strip.text = element_text(size = 15, face = "bold"),
    plot.tag = element_text(face = "bold")
  )

ggsave("./plot/cpue_t.tiff",dpi = 600, width = 12, height = 20, units = "cm")

b

m_b = passive_sample |> 
  dplyr::group_by(Treatment, Position) |> 
  rstatix::get_summary_stats(mass_flux, type = "mean_se")
# library(magick)
# 
# ggsave("mass_try.png",plot = b, dpi = 600, width = 10, height = 25, units = "cm")
# 
# img_plot <- cowplot::ggdraw() + cowplot::draw_image("mass_try.png")
# print(img_plot)
# 
# combined_plot <- a + img_plot + plot_layout(ncol = 2)
# print
############

ggsave("mass.png",dpi = 600, width = 20, height = 25, units = "cm")

ggboxplot(data = passive_sample, x = "Treatment",y = "CPUE",
          color = "Position", 
          add = c("jitter", "mean_se"), 
          #shape = "Treatment"
          width = 0.05
          )+
  facet_wrap(~Position, ncol = 1)+
  theme_Publication(base_size = 14)+
  scale_colour_Publication1()

ggboxplot(data = passive_sample, x = "Treatment", y = "mass_flux",
          color = "Position", 
          add = c("jitter", "mean_se"), 
          width = 0.7) +
  facet_wrap(~Position, ncol = 1)   # Facet into one column
  theme_minimal(base_size = 14) +    # Optional: adjust theme
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels if needed
    plot.margin = margin(10, 10, 10, 10))  # Adjust plot margins



passive_sample |> 
  dplyr::group_by(Position, Treatment) |> 
  #dplyr::summarise(mass_flux = n(), .groups = "keep") |> 
  ggplot(aes(Treatment,mass_flux, colour = Position))+
  geom_point(alpha = 0.1)+
  stat_summary(
    fun.data = mean_cl_boot, #boot
    geom = "errorbar",
    width = 0.1,
    position = position_dodge2(0.6))+
    scale_x_discrete(labels = c("Control", "Low Flow"))+
  facet_wrap(~Position,ncol = 1)



ggplot(passive_sample, aes(x = Treatment, y = CPUE, fill = Position)) +
  geom_col() + 
  theme_minimal() +
  facet_wrap(~Position, ncol = 1)
  labs(title = "Boxplot of Emergence Rate by Treatment and Position")

###########
# Create the first plot (adjust as per your actual plot)
plot1 <- ggplot(passive_sample, aes(x = Treatment, y = CPUE, color = Position)) +
  geom_point() + 
  theme_minimal() + 
  labs(title = "Emergence Rate vs Treatment (By Position)")

# Create the boxplot for emergence rate vs treatment, separated by position
boxplot <- ggplot(passive_sample, aes(x = Treatment, y = CPUE, fill = Position)) +
  geom_boxplot() + 
  theme_minimal() +
  labs(title = "Boxplot of Emergence Rate by Treatment and Position")

library(patchwork)
# Combine the two plots with patchwork
combined_plot <- plot1 + boxplot + plot_layout(ncol = 2)

# Display the combined plot
print(combined_plot)

# Abundance proportion per order per position -----------------------------------------
str(emerg_insect)

insect_long = emerg_insect |> 
  dplyr::filter(!Sampling %in% "Active") |> #Only the passive sampling samples was used
  select(Flume,Sampling, Position,Treatment,
         No..mayflies, No..Other.Nematocera, 
         No..Plecoptera, no..Chironomids, no..Trichoptera, Diptera,-Comments) |> 
  pivot_longer(5:10,
               names_to = "No.Order",
               values_to = "count") |> 
  mutate(Order = str_remove(No.Order, "No\\.\\.|no\\.\\.")) |> 
  dplyr::mutate(Days = 7) |> 
  dplyr::mutate(Area = 1) |> 
  dplyr::mutate(CPUE = count/(Area * Days))

insect_biomass = emerg_insect |> 
  select(Date.2021,Flume,Sampling,Treatment, Size.group,
         Emergence.mass..mg.) |> 
  dplyr::mutate(Days = 7) |> 
  dplyr::mutate(Area = 1) |> 
  dplyr::mutate(biomass = Emergence.mass..mg./(Area * Days)) 


insect_avg = insect_long |> 
  rstatix::group_by(Treatment,Position, Order) |> 
  rstatix::get_summary_stats(CPUE, type = "full")


# percentage composition --------------------------------------------------

total = sum(insect_long$count, na.rm = T)

em_t = insect_long |> 
  dplyr::group_by(Order) |> 
  dplyr::summarise(perc = sum(count, na.rm = T)) |> 
  dplyr::mutate(n = (perc/total)* 100)


# end percentage composition ----------------------------------------------


insect_perc <- insect_long |> 
  dplyr::group_by(Treatment, Position, Order) |> 
  dplyr::summarise(Total = sum(count, na.rm = TRUE), .groups = "keep") |> 
  dplyr::group_by(Treatment,Position) |> 
  dplyr::mutate(Total_Position = sum(Total, na.rm = TRUE)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(Proportion =round((Total / Total_Position) * 100,3)) 

  
openxlsx::write.xlsx(insect_avg, "insect_avg.xlsx")


# Total spider count ------------------------------------------------------
tot_spider = openxlsx::read.xlsx("./seafile/Spider Count October 2021.xlsx", sheet = "Sheet1") |> 
  filter(!Treatment %in% "X") |> 
  select(-Total_Number) |> 
  pivot_longer(3:5,
               names_to = "size",
               values_to = "count") |>
  dplyr::mutate(across(size, as.factor))

str(tot_spider)

tot_spider |> 
  dplyr::group_by(Treatment, size) |> 
  rstatix::get_summary_stats(count, type = "mean_se")|> 
  ggplot(aes(size,mean, fill = Treatment, color = Treatment, 
             shape = Treatment))+
  geom_point(position = pd, size = 3, color = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.3, position = pd, color = "black") +
  scale_fill_viridis_d(option = "plasma")+
  scale_x_discrete(name ="Size class", 
                   limits=c("Small", "Medium", "Large"))+
  theme_Publication(base_size = 14) +
  labs(y = "Mean spider count")

# Tetragnatha data --------------------------------------------------------

tet_dat = read.csv("./original_spider_emergence/Tetragnatha_Sampling.csv", 
                        header = T, sep = ";") |> 
  dplyr::mutate(Treatment = if_else(Treatment == "drought", 
                                    "Low flow", Treatment)) 
  #dplyr::mutate_at(vars(c(Treatment, Flume)), as.factor)

str(tet_count)

tet_count = tet_dat |> 
  dplyr::group_by(Flume,Treatment) |> 
  dplyr::summarise(count = n()) |> 
  dplyr::mutate(Family = "Tetragnathidae")
  #dplyr::mutate(F_T = str_c(Flume,Treatment, sep = "-")) |> 
  #dplyr::mutate_at(vars(c(Treatment)), as.factor)


tet_count = tet_dat |> 
  dplyr::group_by(Flume,Treatment) |> 
  dplyr::summarise(count = n(), .groups = "keep")

tet_dat |> 
  dplyr::group_by(Flume,Treatment) |> 
  dplyr::summarise(count = n(), .groups = "keep") |> 
  ggplot(aes(Treatment,count))+
  geom_point(alpha = 0.4)+
  stat_summary(
    fun.data = mean_cl_boot, #boot
    geom = "errorbar",
    width = 0.1,
    position = position_dodge(0.6)) +
  scale_x_discrete(labels = c("Control", "Low Flow"))

# Tetragnathidae and Lycosidae --------------------------------------------

lyco_count = lyco_data |> 
  dplyr::mutate(Family = "Lycosidae") |> 
  select(-Female) |> 
  dplyr::mutate(Treatment = case_when(Treatment == "C" ~ "Control",
                                      Treatment == "D"~ "Low flow",
                                      TRUE ~ Treatment)) |> 
  rename(count = Total)
tet_count = tet_count 
spider_long = rbind(tet_count,lyco_count)

pj = position_dodge(width = 0.3)

wk5 = passive_sample |> 
  dplyr::filter(Experiment.Week %in% 5)

spider_ = spider_stat |> 
  mutate(
  Family = case_when(
    Family == "lycosidae" ~ "Lycosidae",
    TRUE ~ "Tetragnathidae"
  )
)



spider_ |> 
  dplyr::group_by(Treatment, Family) |> 
  rstatix::get_summary_stats(count, type = "mean_se")|> 
  ggplot(aes(Treatment, mean, fill = Family, color = Family, 
             shape = Family, linetype = Family))+
  geom_point(position = pj, size = 3, color = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.1, position = pj, color = "black") +
  ylim(0, 50)+
  scale_shape_manual(values = c(0,8))+
  scale_x_discrete(labels = c("Control", "Low-flow"))+
  #scale_fill_viridis_d(option = "plasma")+
  #theme_Publication(base_size = 14) +
  labs(y = "Mean number of spiders", x = "Treatment")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.title.y = element_text(face = "bold"),
    axis.title = element_text(size = 12, color = "black",
                              face = "bold"),
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 12, face = "bold", color = "black"),
    legend.title = element_text(size = 14),
    #strip.text = element_text(size = 15, face = "bold"),
    plot.tag = element_text(face = "bold"),
    #plot.margin = margin(t = 0, r = 0, b = 0, l = 10)
  )

ggsave("./plot/spider.png",dpi = 300, width = 14, height = 12, units = "cm")

############################

sp_lyc |> 
  dplyr::group_by(Treatment, Family) |> 
  rstatix::get_summary_stats(count, type = "mean_se")|> 
  ggplot(aes(Treatment, mean, fill = Family, shape = Family, linetype = Family))+
  geom_point(position = pj, size = 3, color = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.1, position = pj, color = "black") +
  scale_shape_manual(values = c(0,8))+
  #ylim(0, 50)+
  scale_x_discrete(labels = c("Control", "Low-flow"))+
  #scale_fill_viridis_d(option = "plasma")+
  #theme_Publication(base_size = 14) +
  #labs(y = "Mean number of spiders", x = "Treatment")+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.title.y = element_text(face = "bold"),
    axis.title = element_text(size = 12, color = "black",
                              face = "bold"),
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 12, face = "bold", color = "black"),
    legend.title = element_text(size = 14),
    #strip.text = element_text(size = 15, face = "bold"),
    plot.tag = element_text(face = "bold"),
    #plot.margin = margin(t = 0, r = 0, b = 0, l = 10)
  )

##########Sediment concentration ######
dat_sed |> 
  dplyr::group_by(Treatment, Week) |> 
  rstatix::get_summary_stats(Concentration_µg_kg, type = "mean_se")|> 
  ggplot(aes(Treatment, mean, fill = Week, shape = Week, linetype = Week))+
  geom_point(position = pj, size = 3, color = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.1, position = pj, color = "black") +
  scale_shape_manual(values = c(5,6))+
  #ylim(0, 0.2)+
  scale_x_discrete(labels = c("Control", "Low-flow"))+
  #scale_fill_viridis_d(option = "plasma")+
  #theme_Publication(base_size = 14) +
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
    legend.title = element_text(size = 14),
    #strip.text = element_text(size = 15, face = "bold"),
    plot.tag = element_text(face = "bold"),
    #plot.margin = margin(t = 0, r = 0, b = 0, l = 10)
  )

ggsave("./plot/pesticide.png",dpi = 300, width = 16, height = 12, units = "cm")

        

# Lycosidae data ----------------------------------------------------------

n = lyco_data |> 
  slice(1)

lyco_data = openxlsx::read.xlsx("./original_spider_emergence/Ale to Collins/Original files/RSM Spider Samples Overview 18.02.22.xlsx",
                                sheet = "Sheet1") |> 
  select(-5:-13) |> 
  slice(-1:-4) |> 
  setNames(n) |> 
  rename(Total = `Total lycosidae estimate`) |> 
  rename(Female = `Wespenspinnen (female)`) |> 
  mutate(across(Treatment, as.factor)) |> 
  dplyr::mutate(across(c(Total, Female), as.numeric))

str(lyco_data)

lyco_data |> 
  dplyr::group_by(Treatment) |> 
  rstatix::get_summary_stats(Total, type = "mean_se") |> 
  ggplot(aes(Treatment,mean)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.1, position = pd, color = "black")+
  geom_point(position = pd, size = 3, color = "black")+
  scale_x_discrete(labels = c("Control", "Low Flow"))

ggplot(lyco_data, aes(x = Treatment, y = Total)) +
  geom_point(alpha = 0.4)+
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    position = position_jitter(0.01)
  )+
  scale_x_discrete(labels = c("Control", "Low Flow"))

tet_count = as.data.frame(tet_count)
str(tet_count)
l_test = pairwise_t_test(Total ~ Treatment, paired = T, 
                         p.adjust.method = "bonferroni",
                         data = lyco_data)



ch_d = cohens_d(count ~ Family, data = spider_stat)

spider_stat = data.frame(spider_stat)

ch_d = lyco_count |> 
  #mutate(across(Treatment, Position, sep = "_")) |> 
  cohens_d(count ~ Treatment, var.equal = F, conf.level = 0.95,
           ci.type = "perc", ci = T)

group1 <- tet_count %>% filter(Treatment == "Control") %>% pull(count)
group2 <- tet_count %>% filter(Treatment == "Low flow") %>% pull(count)
cohen_d_value <- cohens_d(group1,group2)
# Lycosidae_pitfall_Ale ---------------------------------------------------

# pit_lyco = read.csv("./original_spider_emergence/Ale to Collins/R/spiders_pit_new.csv",
#                     header = T, sep = ",")
# 
# 
# pit_lyco |> 
#   dplyr::group_by(Treatment, Distance) |> 
#   rstatix::get_summary_stats(Abu, type = "mean_se")|> 
#   ggplot(aes(Treatment, mean, fill = Distance, color = Distance, 
#              shape = Distance))+
#   geom_point(position = pd, size = 3, color = "black")+
#   geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
#                 width = 0.3, position = pd, color = "black") +
#   scale_fill_viridis_d(option = "plasma")+
#   theme_Publication() +
#   labs(y = "Mean Lycosidae count", x = "Treatment")


# New_Lycosidae_data_Karol ------------------------------------------------
new_lyco = openxlsx::read.xlsx("./Karol_Lycosidae_2021/Lycosidae-2021_Final.xlsx",
                               sheet = "Lycosidae")

lyco_count = new_lyco |>
  dplyr::filter(family %in% "lycosidae") |>
  dplyr::group_by(Flume,Treatment, family) |>
  dplyr::summarise(count = n(), .groups = "keep")

lyco_count |>
  dplyr::group_by(Treatment) |>
  rstatix::get_summary_stats(count, type = "mean_se") |>
  ggplot(aes(Treatment,mean)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.1, position = pd, color = "black")+
  geom_point(position = pd, size = 3, color = "black")+
  scale_x_discrete(labels = c("Control", "Low-Flow"))

# Effect size -------------------------------------------------------------

eff_size_ab = passive_sample |> 
  #mutate(Treatment_Week = paste(Treatment, Position, sep = "_")) |> 
  cohens_d(CPUE ~ Experiment.Week, var.equal = TRUE, conf.level = 0.95,
           ci.type = "perc", ci = T)

eff_size_bio = passive_sample |> 
  cohens_d(mass_flux ~ Experiment.Week, var.equal = TRUE)

# If there is need to plot  -----------------------------------------------
eff_size_ab |> 
  ggplot( aes(x = interaction(group1, group2),y = effsize)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, size = 1.2)+
  coord_flip()

#####################



# Environmental Variable data ---------------------------------------------
f1_cont = read.csv("../env_data/f1_control.csv")
f16_lowflow= read.csv("../env_data/f16_lowflow.csv")

f16_lf = f16_lowflow |> 
  filter(DO >= 70) |> 
  filter(Conductivity >= 200) |> 
  filter(pH >= 6.5 & pH <= 10)

f1_kt = f1_cont |> 
  filter(DO >= 70) |> 
  filter(Conductivity >= 200) |> 
  filter(pH >= 6.5 & pH <= 10) 

str(f1_kt)


f1_kt |> 
  dplyr::select(Datum, Uhrzeit, DO, pH, Conductivity, Cond_temp) |> 
  dplyr::rename(`Dissolved oxygen (% saturation)` = DO,
         `Conductivity (µS/cm)` = Conductivity,
         `Temperature (°C)` = Cond_temp) |> 
  tidyr::pivot_longer(3:6,
               names_to = "PhCh",
               values_to = "value") |>
  dplyr::mutate(month = str_sub(Datum, 4, 5),
    datetime = as.POSIXct(paste(Datum, Uhrzeit), format = "%d.%m.%y %H:%M:%S"),
    date = as.Date(datetime),
    hour = floor_date(datetime, unit = "hour")) |> 
  dplyr::filter(!month %in% c("08","09")) |> 
  group_by(hour, PhCh) |> 
  rstatix::get_summary_stats(type = "full") |> 
  #filter(!n %in% 1) |> 
  ggplot(aes(x = hour, y = mean, color = PhCh)) +
  #geom_ribbon(aes(ymin = mean - ci, ymax = mean + ci), alpha = 0.1) +
  
  geom_line() +
  facet_wrap(~PhCh, scales = "free_y", ncol = 1) +
  labs(
       x = "Date",
       y = "Measurements") +
  theme_Publication(base_size = 12) +
  scale_colour_Publication1()

ggsave("phch.png",dpi = 600, width = 25, height = 15, units = "cm")


f16_sum = f16_lf |>
  mutate(Datum = as.Date(Datum, format="%d.%m.%y")) |> 
  mutate(monat = format(Datum, "%m")) |> 
  dplyr::group_by(monat) |> 
  rstatix::get_summary_stats(type = "full")

f1_sum = f1_kt |> 
  mutate(Datum = as.Date(Datum, format="%d.%m.%y")) |> 
  mutate(monat = format(Datum, "%m")) |> 
  dplyr::group_by(monat) |> 
  rstatix::get_summary_stats(type = "full")

openxlsx::write.xlsx(f16_sum, "f16_lf_summary.xlsx")
openxlsx::write.xlsx(f1_sum, "f1_kt_summary.xlsx")



stargazer::stargazer(f16_lf, type = "text")
stargazer::stargazer(f1_kt, type = "text")


# Verena_Screiner_phch ----------------------------------------------------

str(phch_1)

phch_1 = openxlsx::read.xlsx("../Verena_Schreiner/RSM_2021_phch.xlsx", 
                             sheet = "macroinv", detectDates = TRUE) |> 
  filter(!Flume %in% c(3,8,12,15) )

wtw_1 = readxl::read_excel("../Verena_Schreiner/RSM_2021_phch.xlsx", 
                            sheet = "WTW", 
                           col_types = c("date", "text", 
                                         "numeric","numeric",
                                         "numeric","numeric")) |>  
   mutate(month = format(date_time, "%m")) |> 
   filter(!month %in% c("05","09"))



avg_wtw = wtw_1 |> 
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
  rstatix::get_summary_stats(type = "mean_se")
  # dplyr::summarise(n = n(),
  #                 temp = mean(`T (°C)`),
  #                 cond = mean(`C (uS/cm)`),
  #                 DOD_perc = mean(`O2 (% sat)`),
  #                 DODµmol = mean(`O2 (umol/L)`),
  #                 .groups = "drop")

write.table(avg_wtw, file = "avg_wtw.txt", row.names = FALSE, sep = "\t", quote = FALSE)

a = phch_1 |> 
  filter(phase %in% "colonization") |> 
  group_by(location) |> 
  rstatix::get_summary_stats(type = "full") #estimation of the max velocity and depth in inlet and outlet stretch

phch = phch_1 |> 
  filter(phase %in% "treatment") |>
  dplyr::mutate(across(c(Flume), as.factor)) |> 
  dplyr::group_by(location,treatment) |> 
  rstatix::get_summary_stats(type = "mean_se")
  # dplyr::summarise(n = n(),
  #                  flow_rate = mean(`flow_vel.(m/sec)`),
  #                  depth = mean(`depth.(cm)`),
  #                  .groups = "drop")
  
write.table(phch, file = "phch.txt", row.names = FALSE, sep = "\t", quote = FALSE)












