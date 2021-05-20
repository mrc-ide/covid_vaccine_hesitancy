### ----------------------------------------------------------------
###                 CODE FOR FIGURES 1 AND 2 
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
source ("./R/extract_data.R")
source("./R/functions.R")
source("./R/output_analysis_multinom.R")
source("./R/coverage_matrix.R")     # Code that gets coverage matrix by age group based on samples from survey results 
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

#### Rt profile ####

## 1. Vaccine efficacy 95%

rt95 <- rt_df (r_vacc = r_vacc,r_seq = r_seq,tt_seq = tt_seq,efficacy = efficacy[1],
               daily_vaccine = daily_vaccine,target_pop = target_pop,
               start_vaccination = start_vaccination,day0 = day0,
               end_vaccination = end_vaccination)
rt95$scenario <- "Adults"

 rt_ch_95 <- rt_df(r_vacc = r_vacc,r_seq = r_seq,tt_seq = tt_seq,efficacy = efficacy[1],
                   daily_vaccine = daily_vaccine,target_pop = target_pop,
                   start_vaccination = start_vaccination,day0 = day0, 
                   end_vaccination = end_children)
rt_ch_95$scenario <- "Children"

plot_rt95 <-  plot_rt(rt = rt95,start_vaccination = start_vaccination,end_vaccination = end_vaccination,end_children = end_children,day0 = day0)
  #plot_rt_scenarios (rbind(rt95,rt_ch_95),start_vaccination = start_vaccination,end_vaccination = end_vaccination, 
              #                  end_children=end_children,day0 = day0)

## 2. Vaccine efficacy 70%
rt70 <- rt_df(r_vacc = r_vacc,r_seq = r_seq,tt_seq = tt_seq,efficacy = efficacy[2],
              daily_vaccine = daily_vaccine,target_pop = target_pop,
              start_vaccination = start_vaccination,day0 = day0, 
              end_vaccination = end_vaccination)
rt70$scenario <- "Adults"

rt_ch_70 <- rt_df(r_vacc = r_vacc,r_seq = r_seq,tt_seq = tt_seq,efficacy = efficacy[2],
                   daily_vaccine = daily_vaccine,target_pop = target_pop,
                   start_vaccination = start_vaccination,day0 = day0, 
                   end_vaccination = end_children)
 rt_ch_70$scenario <- "Children"

plot_rt70 <- plot_rt(rt = rt70,start_vaccination = start_vaccination,end_vaccination = end_vaccination,end_children = end_children,day0 = day0)
  #plot_rt_scenarios (rbind(rt70,rt_ch_70),start_vaccination = start_vaccination,end_vaccination = end_vaccination, end_children=end_children,day0 = day0)
#### --------------------- Scenarios  ------------------------####

# 1. Coverage matrix 

cov_matrix <- coverage_draws(by_age)

adul_antivax <- age_antivax(df = cov_matrix,ideal = 0.95,min_age = 15)
adul_antivax$Age <- factor(adul_antivax$Age, levels= c("0-5","5-10","10-15","15-20", "20-25","25-30","30-35", "35-40", "40-45", "45-50",
                                                       "50-55", "55-60" ,"60-65", "65-70", "70-75", "75-80", "80+"))
adul_summary <- melt(adul_antivax[,-101],id.vars = "Age") %>% group_by(Age) %>% summarise_at(vars(value),list(Q1=~quantile(., probs = 0.20),
                                                                                                              median=median, Q3=~quantile(., probs = 0.80))) %>%
  mutate (ideal = adul_antivax$Ideal)  %>%  relocate(Age, .after = last_col())

child_antivax <- age_antivax(df = cov_matrix,ideal = 0.95,min_age = 5)
child_antivax$Age <- factor(child_antivax$Age, levels= c("0-5","5-10","10-15","15-20", "20-25","25-30","30-35", "35-40", "40-45", "45-50",
                                                         "50-55", "55-60" ,"60-65", "65-70", "70-75", "75-80", "80+"))
child_summary <- melt(child_antivax[,-101],id.vars = "Age") %>% group_by(Age) %>% summarise_at(vars(value),list(Q1=~quantile(., probs = 0.20),
                                                                                                            median=median, Q3=~quantile(., probs = 0.80))) %>%
  mutate (ideal = child_antivax$Ideal)%>% relocate(Age, .after = last_col())

# 2. Scenarios 

s_nolift <- list(R0 = r_seq, tt_R0 = tt_seq)
# 95% efficiency UK 
s_novax_95 <- list (R0 = rt95$Rt, tt_R0 = rt95$t)
# 70% efficiency UK 
s_novax_70<-list (R0 = rt70$Rt, tt_R0 = rt70$t)

s_counter <- list(s_nolift,s_novax_95,s_novax_70)

### Running scenarios ####

# No vaccine or interventions 
runs_counter <- lapply (s_counter, function(x){
  
  do.call(run_scenario,x)
  
})

runs_95 <- run_scenarios(mat_antivax = adul_summary,rt_list = list(rt95),eff = efficacy[1],
                         eff_disease = efficacy_disease[1],tt_vax = tt_vaccine, max_vaccine = max_vaccine)

runs_ch_95 <- run_scenarios(mat_antivax = child_summary,rt_list = list(rt95),eff = efficacy[1],
                             eff_disease = efficacy_disease[1],tt_vax = tt_children,max_vaccine = max_vaccine)


runs_70  <-  run_scenarios(mat_antivax = adul_summary,rt_list = list(rt70),eff = efficacy[2],
                           eff_disease = efficacy_disease[2],tt_vax = tt_vaccine,max_vaccine = max_vaccine)


runs_ch_70 <- run_scenarios(mat_antivax = child_summary,rt_list = list(rt70),eff = efficacy[2],
                            eff_disease = efficacy_disease[2],tt_vax = tt_children,max_vaccine = max_vaccine)

### FIGURE 1 ####


no_vax95 <- compartment_out(out_nimue = runs_ch_95[[3]],pcompartment = "deaths")

plot_compartment_dynamics(df=no_vax95, yscale=50,ylab="deaths", day0=day0)+
  xlim(c(as.Date("2020-03-01"),as.Date("2020-12-31")))


## 95% efficiency
cols_names <-c ("Lower" ,"Hesitancy", "Upper", "Ideal") 

to_plot_95 <- to_plot_output(runs_list = runs_95,runs_ch = runs_ch_95,names_vector = cols_names,day0 = day0)

plot_dynamics_95 <- plot_dynamics_short(df = to_plot_95 ,vacc_start = start_vaccination,
                                        vacc_finishes = end_vaccination,vacc_fin_chil = end_children,
                                        yscale = 50,ylabel = "Deaths per million",day0 =  day0) + theme(legend.position = "bottom")


# 70% efficiency

to_plot_70 <- to_plot_output(runs_list = runs_70,runs_ch = runs_ch_70,names_vector = cols_names,day0 = day0)
plot_dynamics_70 <- plot_dynamics_short(df = to_plot_70 ,vacc_start = start_vaccination,
                                        vacc_finishes = end_vaccination,vacc_fin_chil = end_children,
                                        yscale = 50,ylabel = "Deaths per million",day0 =  day0)

top <- plot_grid(plot_rt95+ theme(legend.position = "none"), 
                   plot_dynamics_95 + theme(legend.position = "none"),
                   plot_rt70+ theme(legend.position = "none"), 
                   plot_dynamics_70 +theme(legend.position = "none") ,
                   #averted95 + theme(legend.position = "none"), 
                   ncol=2,  labels= c("a","b","c","d"))

leg <- get_legend(plot_dynamics_95)


total <- plot_grid (top, leg, ncol=1, rel_heights = c(1,0.1))

# To document

x <- to_plot_95 %>% filter(as_date >= as.Date("2022-01-01") & as_date <= as.Date("2023-01-01")& scenario == "Adults")
c <- to_plot_95 %>% filter(as_date >= as.Date("2022-01-01") & as_date <= as.Date("2023-01-01")& scenario == "Children")
max(x$Upper) 
y <- x %>% mutate(cum_hes = cumsum(Hesitancy), 
                                             cum_lower = cumsum(Lower),
                                             cum_upper = cumsum (Upper),
                                             cum_ideal = cumsum(Ideal))

## Figure 2 #### 

# Deaths 

df_95 <- cum_output_vac (runs=runs_95,ini=1,last=4,pcompartment="deaths_cumu",
                         efficacy= "High",start_vaccination=start_vaccination,final_day=final_day,names = cols_names)
#ideal_95 <- cum_output_vac (runs=runs_95,ini=101,last=101,pcompartment="deaths_cumu",
                            #efficacy= "High",scenario= "Ideal",start_vaccination=start_vaccination,final_day=final_day)

df_70 <- cum_output_vac (runs=runs_70,ini=1,last=4,pcompartment="deaths_cumu",
                         efficacy= "Moderate",start_vaccination=start_vaccination,final_day=final_day,names = cols_names)

deaths_df <- rbind(df_95,df_70)

# Hospitalisation 

hf_95 <- cum_output_vac (runs=runs_95,ini=1,last=4,pcompartment="hospitalisations_cumu",
                         efficacy= "High",start_vaccination=start_vaccination,final_day=final_day,names = cols_names)

hf_70 <- cum_output_vac (runs=runs_70,ini=1,last=4,pcompartment="hospitalisations_cumu",
                         efficacy= "Moderate",start_vaccination=start_vaccination,final_day=final_day, names= cols_names)


hosp_df <- rbind(hf_95, hf_70)


# Plot
#Bar plot 

deaths_plot <-plot_bar_vacc_status (df = deaths_df, yscale =  50, ylab ="Cumulative deaths \n per million" ) + 
  theme(legend.position = "bottom")
hosp_plots <- plot_bar_vacc_status (df = hosp_df, yscale =  50, ylab ="Cumulative Hospitalisations \n per million")

legend_bar <- get_legend(deaths_plot)

plot_grid(deaths_plot + theme(legend.position = "none"), 
          hosp_plots+ theme(legend.position = "none"), 
          legend_bar,
          nrow=3,
          rel_heights = c(1,1,0.1),labels = c("a","b"))



## Supplementary Fig 1 ####


coverage_plot <- plot_coverage(cov_matrix)

## Supplementary Fig 2 #### 

# Deaths 

chf_95 <- cum_output_vac (runs=runs_ch_95,ini=1,last=4,pcompartment="deaths_cumu",
                         efficacy= "High",start_vaccination=start_vaccination,final_day=final_day,names = cols_names)

chf_70 <- cum_output_vac (runs=runs_ch_70,ini=1,last=4,pcompartment="deaths_cumu",
                         efficacy= "Moderate",start_vaccination=start_vaccination,final_day=final_day,names = cols_names)

deaths_chf <- rbind(chf_95, chf_70)

# Hospitalisation 


hch_95 <- cum_output_vac (runs=runs_ch_95,ini=1,last=4,pcompartment="hospitalisations_cumu",
                         efficacy= "High",start_vaccination=start_vaccination,final_day=final_day, names= cols_names)

hch_70 <- cum_output_vac (runs=runs_ch_70,ini=1,last=4,pcompartment="hospitalisations_cumu",
                         efficacy= "Moderate",start_vaccination=start_vaccination,final_day=final_day, names = cols_names)


hosp_chf <- rbind(hch_95, hch_70)


# Plot
#Bar plot 

deaths_plot_ch <-plot_bar_vacc_status (df = deaths_chf, yscale =  50, ylab ="Cumulative deaths \n per million" ) + 
  theme(legend.position = "bottom")
hosp_plots_ch <- plot_bar_vacc_status (df = hosp_chf, yscale =  50, ylab ="Cumulative Hospitalisations \n per million")

legend_bar <- get_legend(deaths_plot)

plot_grid(deaths_plot_ch + theme(legend.position = "none"), 
          hosp_plots_ch+ theme(legend.position = "none"), 
          legend_bar,
          nrow=3,
          rel_heights = c(1,1,0.1),labels = c("a","b"))






# ideal_95 <- compartment_out(out_nimue = runs_95[[101]],pcompartment = "deaths") %>% mutate(as_date= as.Date(t, origin = day0))
# ideal_95$scenario <- "Adults"
# ch_95_ideal <-  compartment_out(out_nimue = runs_ch_95[[101]],pcompartment = "deaths") %>% mutate(as_date= as.Date(t, origin = day0))
# ch_95_ideal$scenario <- "Children"
# 
# ideal_95 <- rbind(ideal_95,ch_95_ideal)
# 
# 
# adul_95 <- summary_runs (runs_list= runs_95, ini=1, fini=100, pcompartment= "deaths",day0=day0, scenario = "Adults")
# ch_95 <- summary_runs (runs_list= runs_ch_95, ini=1, fini=100, pcompartment= "deaths",day0=day0,scenario = "Children")
# adul_95 <- rbind (adul_95,ch_95)
# 
# 
# plot_dynamics_95 <- plot_dynamics_multi(hesitancy_df = adul_95,ideal_df = ideal_95,vacc_start = start_vaccination,
#                                         vacc_finishes = end_vaccination,vacc_fin_chil = end_children,
#                                         yscale = 50,ylabel = "Deaths per million",day0 =  day0) + theme(legend.position = "bottom")