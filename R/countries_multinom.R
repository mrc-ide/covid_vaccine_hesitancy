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
source("./R/functions.R")
source("./R/output_analysis_multinom.R")
source("./R/coverage_matrix.R") 


#### ---------------------  GLOBAL PARAMETERS  ------------------------ ####

ideal_coverage <- 0.8
efficacy <- c(0.94,0.63)
efficacy_disease <- c(0.6,0.6)

age_vec <- c("0-5","5-10","10-15","15-20", "20-25","25-30","30-35", "35-40", "40-45", "45-50",
             "50-55", "55-60" ,"60-65", "65-70", "70-75", "75-80", "80+")

#### ---------------------  AntiVax Views  ------------------------ ####

cov_france <- coverage_draws(france)
cov_france_ad <-  summary_coverages(cov_matrix = cov_france)
cov_france_ch <- summary_coverages(cov_matrix = cov_france,adults =  FALSE)

cov_germany <- coverage_draws(germany)
cov_germany_ad <-  summary_coverages(cov_matrix = cov_germany)
cov_germany_ch <-  summary_coverages(cov_matrix = cov_germany,adults =  FALSE)

cov_uk <- coverage_draws(uk)
cov_uk_ad <- summary_coverages(cov_matrix = cov_uk)
cov_uk_ch <- summary_coverages(cov_matrix = cov_uk,adults = FALSE)

#### ---------------------  Running scenarios   ------------------------ ####

# 1. Get country fittings: 
fit_france <- get_rt_fit ("FRA")
fit_germany <- get_rt_fit ("DEU")
fit_uk <- get_rt_fit ("GBR")

# 2. Get parameters:  #####

param_france_95 <- set_parameters(rt_fit = fit_france, country = "France",efficacy = efficacy[1])
param_france_70 <- set_parameters(rt_fit = fit_france, country = "France",efficacy = efficacy[2])

param_germany_95 <- set_parameters(rt_fit = fit_germany, country = "Germany",efficacy = efficacy[1])
param_germany_70 <- set_parameters(rt_fit = fit_germany, country = "Germany",efficacy = efficacy[2])

param_uk_95 <- set_parameters(rt_fit = fit_uk, country = "United Kingdom",efficacy = efficacy[1])
param_uk_70 <- set_parameters(rt_fit = fit_uk, country = "United Kingdom",efficacy = efficacy[2])

# 3. Run the scenarios  ####

out_france <- run_countries(country = "France",mat_cov = cov_france_ad,param_95 = param_france_95,param_70 = param_france_70,
                            efficacy_vec =efficacy,efficacy_disease_vec = efficacy_disease,scenario = "Adults") 

out_germany <- run_countries(country = "Germany",mat_cov = cov_germany_ad,param_95 = param_germany_95,param_70 = param_germany_70,
                             efficacy_vec =efficacy,efficacy_disease_vec = efficacy_disease,scenario = "Adults") 

out_uk <- run_countries(country = "United Kingdom",mat_cov = cov_uk_ad,param_95 = param_uk_95,param_70 = param_uk_70,
                        efficacy_vec =efficacy,efficacy_disease_vec = efficacy_disease,scenario = "Adults") 

# 4. Get cumulative per country  ####


cum_france <- (out_france) %>% filter (t>param_france_95$start_vaccination, t <= (param_france_95$end_vaccination+730)) %>%
  group_by(efficacy,age_group) %>% mutate(cum_hes =cumsum(Hesitancy), 
                                          cum_lower = cumsum(Lower),
                                          cum_upper = cumsum (Upper),
                                          cum_ideal = cumsum(Ideal),
                                          country = "France") %>%
  filter (t == (param_france_95$start_vaccination+730))

cum_germany <- out_germany %>% filter (t>param_germany_95$start_vaccination,t <= (param_germany_95$end_vaccination+730)) %>%
  group_by(efficacy,age_group) %>% mutate(cum_hes =cumsum(Hesitancy), 
                                          cum_lower = cumsum(Lower),
                                          cum_upper = cumsum (Upper),
                                          cum_ideal = cumsum(Ideal),
                                          country = "Germany") %>%
  filter (t == (param_germany_95$start_vaccination+730))

cum_uk <- out_uk %>% filter (t>param_uk_95$start_vaccination,t <= (param_uk_95$end_vaccination+730)) %>%
  group_by(efficacy,age_group) %>% mutate(cum_hes =cumsum(Hesitancy), 
                                          cum_lower = cumsum(Lower),
                                          cum_upper = cumsum (Upper),
                                          cum_ideal = cumsum(Ideal),
                                          country = "U.K.") %>%
  filter (t == (param_uk_95$start_vaccination+730)) #730 or 1095

cum_all <- rbind(cum_france,cum_germany,cum_uk)



# 

# 5. Finding death ratios:

cum_all <- cum_all  %>% group_by(country,age_group,efficacy) %>% 
  mutate (ratio = cum_hes/cum_ideal,
          ratio_lower = cum_lower/cum_ideal,
          ratio_upper = cum_upper/cum_ideal)

# to document

total <- cum_all %>% group_by (country,efficacy) %>% select(t,cum_hes,cum_lower,cum_upper,cum_ideal) %>%summarize_all(sum) %>%
  mutate (ratio = cum_hes/cum_ideal,
          ratio_lower = cum_lower/cum_ideal,
          ratio_upper = cum_upper/cum_ideal)

## ---------------------------   FIGURES ----------------------------- ####

# 1. Death Ratios ####


plot_ratios <- plot_bar_deaths_comp (cum_all)
#plot_bar_deaths_comp (cum_all %>% filter(efficacy =="High"))


# 2. Disease dynamics  ####

dyn_france <- plot_dynamics_countries(big_df= out_france %>% filter(efficacy== "High"), 
                                      vacc_start=param_france_95$start_vaccination, 
                                      vacc_finishes=param_france_95$end_vaccination,yscale=1 , 
                                      ylabel= "Total deaths", day0=param_france_95$day0)+ ylim(0,5000)



dyn_germany <- plot_dynamics_countries(big_df= out_germany %>% filter(efficacy== "High"),
                                       vacc_start=param_germany_95$start_vaccination, 
                                       vacc_finishes=param_germany_95$end_vaccination,yscale=1 , 
                                       ylabel= "Total deaths", day0=param_germany_95$day0) + ylim(0,5000)


dyn_uk <- plot_dynamics_countries(big_df= out_uk %>% filter(efficacy== "High"), 
                                  vacc_start=param_uk_95$start_vaccination, 
                                  vacc_finishes=param_uk_95$end_vaccination,yscale=1 , 
                                  ylabel= "Total deaths", day0=param_uk_95$day0)+ ylim(0,5000)



# 3. Vaccine hesitancy  #### 

plot_hes <- plot_overal_hesitancy (country_list = c("France", "Germany","United Kingdom"),cov_matrix_list = list(cov_france_ad,cov_germany_ad,cov_uk_ad),ini = 1,fini = 4)

# 4. Rt Profiles   ####

last_date <- as.Date("2023-01-01")

rt_france95 <- param_france_95$rt_ad %>% slice(c(1:n(),n())) %>%  mutate (efficacy = "High", country = "France")
rt_france95$asDate[nrow(rt_france95)] <- last_date
rt_france70 <- param_france_70$rt_ad %>% slice(c(1:n(),n())) %>% mutate (efficacy = "Moderate", country = "France")
rt_france70$asDate[nrow(rt_france70)] <- last_date

rt_germany95 <- param_germany_95$rt_ad %>% slice(c(1:n(),n())) %>% mutate (efficacy = "High", country = "Germany")
rt_germany95$asDate[nrow(rt_germany95)] <- last_date
rt_germany70 <- param_germany_70$rt_ad %>% slice(c(1:n(),n())) %>% mutate (efficacy = "Moderate", country = "Germany")
rt_germany70$asDate[nrow(rt_germany70)] <- last_date

rt_uk95 <- param_uk_95$rt_ad %>% slice(c(1:n(),n())) %>% mutate (efficacy = "High", country = "U.K")
rt_uk95$asDate[nrow(rt_uk95)] <- last_date
rt_uk70 <- param_uk_70$rt_ad %>% slice(c(1:n(),n())) %>% mutate (efficacy = "Moderate", country = "U.K")
rt_uk70$asDate[nrow(rt_uk70)] <- last_date

big_df <- rbind(rt_france95, rt_france70,rt_germany95,rt_germany70,rt_uk95,rt_uk70)

plot_rt <-plot_profile_2D (big_df) + theme(legend.position = "bottom")


# rt_france <-  plot_rt(rt = param_france_70$rt_ad,start_vaccination = param_france_95$start_vaccination,end_vaccination = param_france_95$end_vaccination,
#                       end_children = param_france_95$end_vaccination,day0 = param_france_95$day0)
# rt_germany <- plot_rt(rt = param_germany_95$rt_ad,start_vaccination = param_germany_95$start_vaccination,end_vaccination = param_germany_95$end_vaccination,
#                       end_children = param_germany_95$end_vaccination,day0 = param_germany_95$day0)
# rt_uk <- plot_rt(rt = param_uk_70$rt_ad,start_vaccination = param_uk_95$start_vaccination,end_vaccination = param_uk_95$end_vaccination,
#                  end_children = param_uk_95$end_vaccination,day0 = param_uk_95$day0)



# Combining everything:  ####

# Density 
# middle <- plot_grid(hes_france + theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
#                     hes_germany + theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
#                     hes_uk + theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
#                     ncol=3,
#                     scale=0.9) +
#   draw_label("Vaccine uptake", x=0.5, y=  0, vjust=0, angle= 0,size=12) +
#   draw_label("Density ", x=  0, y=0.5, vjust= 1.5, angle=90,size=12)

# Rt 
leg_c <- get_legend(plot_rt)
bottom <-plot_grid(plot_hes,
                   plot_rt + theme(legend.position="none"),
                   ncol=2,
                   labels= c("b","c"))
                   
plot_all <- plot_grid (plot_ratios,bottom,leg_c, ncol=1, rel_heights = c(1,0.7,0.1), labels = c("a",""))

# Supplementary  #### 

s_top <- plot_bar_deaths_all (cum_all %>% filter(efficacy == "High"))

s_bottom <- plot_grid(dyn_france + theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
                      dyn_germany + theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
                      dyn_uk + theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
                      ncol=1,
                      scale=0.9) +
  draw_label("Time(days)", x=0.5, y=  0, vjust=0, angle= 0,size=10) +
  draw_label("Total Deaths ", x=  0, y=0.5, vjust= 1.5, angle=90,size=10)

s_total <- plot_grid (s_bottom,s_top, ncol=2, labels= c("a","b"))
