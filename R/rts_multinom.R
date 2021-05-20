### ----------------------------------------------------------------
###
### ----------------------------------------------------------------

#### Load required packages  ####

library(ggplot2)
library(dplyr)
library(tableone)
library(cowplot)
library(purrr)
library(tidyr)
library(nimue)
library(furrr)
library(forcats)
library(jcolors)
library(RColorBrewer)
library(reshape2)

## Load functions
source("./R/coverage_matrix.R") 
source("./R/functions.R")
source("./R/extract_data.R")
source("./R/output_analysis_multinom.R")

#### --------------------- Parameters ------------------------ ####

# Model parameters
target_pop <-  50e6
one_year <- 365

# Vaccine parameters
efficacy <- c(0.94,0.63)
efficacy_disease <- c(0.6,0.6)

# Vaccination rate parameters
pop <- squire::get_population(country = "United Kingdom")$n
pop_standardise <- target_pop / sum(pop)
pop <- pop * pop_standardise

daily_vaccine <-  sum(pop[4:17])/304  # Overall all individuals 15+ can be vaccinated in 10 months 
max_vaccine <-  c(0,daily_vaccine,0)


# Time frame parameters
start_vaccination <- 320
end_vaccination <- 320 + 304
tt_vaccine <- c(0, start_vaccination, end_vaccination)

extra_days <- sum(pop[2:3])/daily_vaccine  #  Extra days needed to vaccinate children
end_children <- ceiling(end_vaccination+extra_days)
tt_children <- c(0, start_vaccination, end_children)

final_day <- start_vaccination + (one_year*2) # tow year after vaccination starts 

day0 <- as.Date("2020-02-14")

# Rt sequences
r_seq <- c(3,1,1.5)  # Before vaccination
tt_seq <- c(0,55,200)
r_vacc <- seq(from=1.5, to=3,length.out = 10)  # After vaccination

#### Antivax views ####
cov_matrix <- coverage_draws(by_age)

adul_antivax <- age_antivax(df = cov_matrix,ideal = 0.95,min_age = 15)

adul_antivax$Age <- factor(adul_antivax$Age, levels= c("0-5","5-10","10-15","15-20", "20-25","25-30","30-35", "35-40", "40-45", "45-50",
                                                       "50-55", "55-60" ,"60-65", "65-70", "70-75", "75-80", "80+"))

adul_summary <- melt(adul_antivax[,-101],id.vars = "Age") %>% group_by(Age) %>% summarise_at(vars(value),list(Q1=~quantile(., probs = 0.10),
                                                                                                  median=median, Q3=~quantile(., probs = 0.90))) %>%
  mutate (ideal = adul_antivax$Ideal)  %>%
  relocate(Age, .after = last_col())


child_antivax <- age_antivax(df = cov_matrix,ideal = 0.95,min_age = 5)
child_antivax$Age <- factor(child_antivax$Age, levels= c("0-5","5-10","10-15","15-20", "20-25","25-30","30-35", "35-40", "40-45", "45-50",
                                                         "50-55", "55-60" ,"60-65", "65-70", "70-75", "75-80", "80+"))

child_summary <- melt(child_antivax[,-101],id.vars = "Age") %>% group_by(Age) %>% summarise_at(vars(value),list(Q1=~quantile(., probs = 0.10),
                                                                                                               median=median, Q3=~quantile(., probs = 0.90))) %>%
  mutate (ideal = child_antivax$Ideal)  %>%
  relocate(Age, .after = last_col())


# Rt profiles ####

rt_ad_95 <- rt_profiles(mat_antivax = adul_summary,r_vacc = r_vacc,r_seq = r_seq,tt_seq = tt_seq, 
                        efficacy = efficacy[1],daily_vaccine = daily_vaccine,target_pop = target_pop,
                        start_vaccination = start_vaccination,end_vaccination = end_vaccination,
                        day_0 = day0,pop_byage = pop)

rt_ch_95 <- rt_profiles(mat_antivax = child_summary,r_vacc = r_vacc,r_seq = r_seq,tt_seq = tt_seq, 
                        efficacy = efficacy[1],daily_vaccine = daily_vaccine,target_pop = target_pop,
                        start_vaccination = start_vaccination,end_vaccination = end_children,
                        day_0 = day0,pop_byage = pop)

rt_ad_70 <- rt_profiles(mat_antivax = adul_summary,r_vacc = r_vacc,r_seq = r_seq,tt_seq = tt_seq, 
                        efficacy = efficacy[2],daily_vaccine = daily_vaccine,target_pop = target_pop,
                        start_vaccination = start_vaccination,end_vaccination = end_vaccination,
                        day_0 = day0,pop_byage = pop)

rt_ch_70 <- rt_profiles(mat_antivax = child_summary,r_vacc = r_vacc,r_seq = r_seq,tt_seq = tt_seq, 
                        efficacy = efficacy[2],daily_vaccine = daily_vaccine,target_pop = target_pop,
                        start_vaccination = start_vaccination,end_vaccination = end_children,
                        day_0 = day0,pop_byage = pop)

## Running scenarios #### 

runs_95 <- run_scenarios(mat_antivax = data.frame(adul_summary),rt_list = rt_ad_95,eff = efficacy[1],
                         eff_disease = efficacy_disease[1],tt_vax = tt_vaccine, max_vaccine = max_vaccine)

runs_ch_95 <- run_scenarios(mat_antivax = data.frame(child_summary),rt_list = rt_ch_95,eff = efficacy[1],
                            eff_disease = efficacy_disease[1],tt_vax = tt_children, max_vaccine = max_vaccine)

runs_70  <-  run_scenarios(mat_antivax = data.frame(adul_summary),rt_list = rt_ad_70,eff = efficacy[2],
                           eff_disease = efficacy_disease[2],tt_vax = tt_vaccine, max_vaccine = max_vaccine)

runs_ch_70 <- run_scenarios(mat_antivax = data.frame(child_summary),rt_list = rt_ch_70,eff = efficacy[2],
                            eff_disease = efficacy_disease[2],tt_vax = tt_children, max_vaccine = max_vaccine)
## ----------  Plots      --------- ####

# Figure 3: Rt Profile ####

## Adults 

scenarios <- c("Upper", "Mean","Lower","Ideal")

toplot_rt95 <- bind_profiles(list_rt=rt_ad_95[2:5], scenario_list=scenarios,vacc_start=start_vaccination,
                             vacc_finishes= end_vaccination,day0= day0)
toplot_rt95$Efficacy <- "High"
toplot_rt70 <- bind_profiles(list_rt=rt_ad_70[2:5], scenario_list=scenarios,vacc_start=start_vaccination,
                             vacc_finishes= end_vaccination,day0= day0)
toplot_rt70$Efficacy <- "Moderate"


big_df <- rbind(toplot_rt95,toplot_rt70)
big_df$scenario <- factor(big_df$scenario, levels = c("Lower","Upper", "Mean","Ideal"))

colors <- c ("#ddb6a1","#ddb6a1","#ba6d43","#030303")
plot_adults<- ggplot(big_df, aes(x=asDate))  +
  geom_step(data=big_df, aes(x = asDate, y = Rt, color = scenario,linetype= Efficacy),size=0.8) + theme_bw() + 
  scale_color_manual(values= colors,name= "Profile")+
  geom_vline(xintercept = day0+ start_vaccination, linetype="dotted", color = "#696969", size=0.8)+
  geom_vline(xintercept = day0+ end_vaccination, color = "#696969", linetype="dotted", size=0.8)+
  labs(y = "Rt", x = "Time(days)") +
  scale_x_date(date_breaks = "6 month" , date_labels = "%b-%y",limits=as.Date (c(NA,"2023-01-06")))+
  theme_bw(base_family = 12) +
  theme(legend.position="none",
        strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line())

## Children 

toplot_rt95 <- bind_profiles(list_rt=rt_ch_95[2:5], scenario_list=scenarios,vacc_start=start_vaccination,
                             vacc_finishes= end_vaccination,day0= day0)
toplot_rt95$Efficacy <- "High"

toplot_rt70 <- bind_profiles(list_rt=rt_ch_70[2:5], scenario_list=scenarios,vacc_start=start_vaccination,
                             vacc_finishes= end_vaccination,day0= day0)
toplot_rt70$Efficacy <- "Moderate"


big_df <- rbind(toplot_rt95,toplot_rt70)
big_df$scenario <- factor(big_df$scenario, levels = c("Lower","Upper", "Mean","Ideal"))

colors <- c ("#dabed6","#dabed6","#824979","#030303")
plot_children <- ggplot(big_df, aes(x=asDate))  +
  geom_step(data=big_df, aes(x = asDate, y = Rt, color = scenario,linetype= Efficacy),size=0.8) + theme_bw() + 
  scale_color_manual(values= colors,name= "Profile")+
  geom_vline(xintercept = day0+ start_vaccination, linetype="dotted", color = "#696969", size=0.8)+
  geom_vline(xintercept = day0+ end_children, color = "#696969", linetype="dotted", size=0.8)+
  labs(y = "Rt", x = "Time(days)") +
  scale_x_date(date_breaks = "6 month" , date_labels = "%b-%y",limits=as.Date (c(NA,"2023-01-06")))+
  theme_bw(base_family = 12) +
  theme(legend.position="none",
        strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line())


plot_grid(plot_adults + theme(legend.position = "none"), 
          plot_children+ theme(legend.position = "none"),
          ncol=2,  labels= c("a","b"))



##  Supplementaru Figure: Dynamics ####

## 95% efficiency
cols_names <-c ("Lower" ,"Hesitancy", "Upper", "Ideal") 

to_plot_95 <- to_plot_output(runs_list = runs_95,runs_ch = runs_ch_95,names_vector = cols_names,day0 = day0)

dynamics_95 <- plot_dynamics_short(df = to_plot_95,vacc_start = start_vaccination,vacc_finishes = end_vaccination,vacc_fin_chil = end_children,
                    yscale = 50,ylabel = "Deaths per million",day0 = day0)


## 70% efficiency
to_plot_70 <- to_plot_output(runs_list = runs_70,runs_ch = runs_ch_95,names_vector = cols_names,day0=day0)

dynamics_70 <- plot_dynamics_short(df = to_plot_70, vacc_start = start_vaccination,vacc_finishes = end_vaccination,vacc_fin_chil = end_children,
                                   yscale = 50,ylabel = "Deaths per million",day0 = day0)
leg <- get_legend(dynamics_95)

total <- plot_grid(dynamics_95 + theme(legend.position = "none"),
                   dynamics_70 + theme(legend.position = "none"),
                   leg,#averted95 + theme(legend.position = "none"), 
                   ncol=1,  labels= c("a","b"), rel_heights = c(1,1,0.1))
