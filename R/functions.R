
#### ---------------------  Country simulataions functions  ------------------------ ####

 # Get country fits 
get_rt_fit <- function (country_id, last_date = "2020-12-31"){
  
  # Pull latest country fit
  iso3c <- country_id
  country <- squire::population$country[squire::population$iso3c==iso3c][1]
  json_path <- file.path("https://raw.githubusercontent.com/mrc-ide/global-lmic-reports/master/",iso3c,"input_params.json")
  json <- jsonlite::read_json(json_path)
  
  dates<- c(unlist(lapply(json, "[[", "date")))
  tt_R0 <- c(unlist(lapply(json, "[[", "tt_beta")) + 1)
  rt <- c(unlist(lapply(json, "[[", "Rt")))
  
  out <- data.frame(cbind(dates,tt_R0,rt)) %>% filter(dates <= last_date)
  
  
}

# Get parameters to run the model based on fitting results and vaccine efficacy 
set_parameters <- function(rt_fit,country,efficacy){
  
  vacc_length =304
  
  rt_fit$dates <- as.Date(rt_fit$dates)
  rt_fit$tt_R0 <- as.numeric(rt_fit$tt_R0)
  rt_fit$rt <- as.numeric(rt_fit$rt)
  
  day0 <- rt_fit$dates[1]
  
  # Vaccination rate parameters
  pop <- squire::get_population(country = country)$n
  target_pop <- sum(pop)
  
  daily_vaccine <-  sum(pop[4:17])/vacc_length  # Overall all individuals 15+ can be vaccinated in 10 months 
  max_vaccine <-  c(0,daily_vaccine,0)
  
  
  # Time frame parameters
  start_vaccination <- rt_fit$tt_R0 [nrow(rt_fit)]+1    # Start vaccination January 1st
  end_vaccination <- start_vaccination + vacc_length
  tt_vaccine <- c(0, start_vaccination, end_vaccination)
  
  
  extra_days <- sum(pop[2:3])/daily_vaccine  #  Extra days needed to vaccinate children
  end_children <- ceiling(end_vaccination+extra_days)
  tt_children <- c(0, start_vaccination, end_children)
  
  #rt_fit$rt[nrow(rt_fit)] <- 1.5
  # Rt profile 
  last <-rt_fit$rt[nrow(rt_fit)]
  first <- 4.6 #rt_fit$rt[1]
  
  r_vacc <- seq(from=last, to=first,length.out = 10)  # After vaccination
  
  rt_ad <- rt_df (r_vacc = r_vacc,r_seq = rt_fit$rt,tt_seq = rt_fit$tt_R0,efficacy = efficacy,
                  daily_vaccine = daily_vaccine,target_pop = target_pop,
                  start_vaccination = start_vaccination,day0 = day0,
                  end_vaccination = end_vaccination)
  
  rt_ch <- rt_df (r_vacc = r_vacc,r_seq = rt_fit$rt,tt_seq = rt_fit$tt_R0,efficacy = efficacy,
                  daily_vaccine = daily_vaccine,target_pop = target_pop,
                  start_vaccination = start_vaccination,day0 = day0,
                  end_vaccination = end_children)
  
  out <- list( day0= day0, start_vaccination = start_vaccination, 
               end_vaccination = end_vaccination, rt_ad = rt_ad,
               rt_ch = rt_ch, end_children=end_children,
               daily_vaccine = daily_vaccine,
               tt_ad = tt_vaccine, tt_ch = tt_children)
  
  
  return (out)
  
}

# Function that runs the model for a country 
run_country <- function(country,R0, tt_R0,tt_vaccine, max_vaccine,
                        efficacy=0,target = matrix(0, nrow=1, ncol=17),
                        time_period = (365 * 4),                      
                        immunosenescence = 1, 
                        duration_R = 365, duration_V=Inf,
                        efficacy_disease=0){
  

  efficacy <- set_efficacy(efficacy = efficacy, efficacy_disease= efficacy_disease, immunosenescence = immunosenescence)
  
  # Run
  
  r1 <- nimue::run(
    time_period = time_period,
    R0 = R0, 
    tt_R0 = tt_R0,
    country = country,
    seeding_cases = 5,
    max_vaccine = max_vaccine,
    tt_vaccine = tt_vaccine,
    dur_V = duration_V,
    vaccine_efficacy_infection = efficacy$vaccine_efficacy_infection,
    vaccine_efficacy_disease = efficacy$vaccine_efficacy_disease,
    vaccine_coverage_mat = target,
    dur_R = duration_R
  )
  
  return (r1)
  
}  


# Function that runs the model for a specific country for both high and moderate vaccine efficacy. 
# @param country           Name of country, it will be used by the function that runs the model
# @param mat_cov           Coverages matrices, each column has the final coverage per age group
# @param param_95          Parameters for a high efficacy 
# @param param70           Parameters for a moderate efficacy
# @param efficacy_vec      Vector with infection efficacy
# @param efficacy_disease  Vector with disease efficacies
# @param scenario          Adults or children 
run_countries <- function (country, mat_cov,param_95,param_70,efficacy_vec,efficacy_disease_vec, scenario ){
  
  day0 <- param_95$day0  
  max_vaccine <-  c(0,param_95$daily_vaccine,0)
  
  
  country_95<- run_scenarios(mat_antivax = mat_cov,rt_list = list(param_95$rt_ad),eff = efficacy_vec[1],
                             eff_disease = efficacy_disease_vec[1],tt_vax = param_95$tt_ad,max_vaccine = max_vaccine,
                             country=country,country_model = 1)
  
  country_70 <-  run_scenarios(mat_antivax = mat_cov,rt_list = list(param_70$rt_ad),eff = efficacy_vec[2],
                               eff_disease = efficacy_disease_vec[2],tt_vax = param_70$tt_ad,max_vaccine = max_vaccine,
                               country=country,country_model = 1)
  # Output analyis 
  names_vector<- c ("Lower" ,"Hesitancy", "Upper", "Ideal")
  
  high <- compartment_out(out_nimue = country_95[[1]],pcompartment = "deaths",reduce_age = FALSE) %>% mutate(as_date= as.Date(t, origin = day0))
  colnames(high)[1] <- names_vector[1]
  
  low <-  compartment_out(out_nimue = country_70[[1]],pcompartment = "deaths",reduce_age = FALSE) %>% mutate(as_date= as.Date(t, origin = day0))
  colnames(low)[1] <- names_vector[1]
  
  for ( i in 2:length(names_vector)){
    
    # High efficacy
    temp <- compartment_out(out_nimue = country_95[[i]],pcompartment = "deaths", reduce_age = FALSE) %>% mutate(as_date= as.Date(t, origin = day0))
    colnames(temp)[1] <- names_vector[i]
    
    high <- high %>% left_join(temp)
    
    #Moderate efficacy
    temp2 <- compartment_out(out_nimue = country_70[[i]],pcompartment = "deaths",reduce_age = FALSE) %>% mutate(as_date= as.Date(t, origin = day0))
    colnames(temp2)[1] <- names_vector[i]
    
    low <- low %>% left_join(temp2)
    
  }
  
  
  high$efficacy <- "High"
  low$efficacy <- "Moderate"
  
  return (rbind(high,low))
  
}




#### ---------------------  Rt Profiles functions  ------------------------ ####

# Get times for lifting intervention based on Rt vector, 
# Vaccine coverage and Vaccine effciency 

time_estimation <- function (Rt,efficacy,daily_vaccines,target_pop ){
  
  coverage_needed <- (1-1/Rt)*(1/efficacy)
  
  times <- (coverage_needed*target_pop)/daily_vaccines
  
  return <- times
}


# Rt_ profile based on rt seq. 
rt_df <- function (r_vacc,r_seq, tt_seq, efficacy, daily_vaccine,target_pop, start_vaccination,day0,end_vaccination){
  
  lifting <- time_estimation (Rt = r_vacc ,efficacy = efficacy,daily_vaccines = daily_vaccine,target_pop =  target_pop)
  
  tt_vacc <- lifting + start_vaccination
  
  rt<- data.frame(cbind(c(r_seq,r_vacc),tt_R0 = c(tt_seq,tt_vacc)))
  colnames(rt) <- c("Rt","t")
  rt$asDate <- as.Date(rt$t, origin = day0)
  
  
  pos_min <-  min(which(rt$t > end_vaccination))

  if(pos_min <= length(rt$t)){

    rt$Rt[pos_min:length(rt$Rt)] <- rt$Rt[pos_min-1]

  }


  return (rt)
  
}


## Function that returs Rt profiles for different vaccination coverage profile. 
# The Rt profile is given such that herd immunity thersholds are mantained based on 
# vaccine coverage. 
rt_profiles <- function(mat_antivax, r_vacc, r_seq,tt_seq,efficacy,daily_vaccine,target_pop,start_vaccination,
                        day_0,end_vaccination,pop_byage){
  
  rt_original <- rt_df (r_vacc = r_vacc,r_seq = r_seq,tt_seq = tt_seq,efficacy = efficacy,
                        daily_vaccine = daily_vaccine,target_pop = target_pop,
                        start_vaccination = start_vaccination,day0 = day0,
                        end_vaccination = end_vaccination)
  
  Rts <- list(rt_original)
  
  for (i in 1:(ncol(mat_antivax)-1)){
    
    rt_profile <- rt_original
    
    pop_vacc <- pop_byage * mat_antivax[,i]  # Estimate how many people per age groups gets vaccinated
    
    total_vacc <- sum(pop_vacc)             # Total people vaccinated
    cov <- total_vacc/target_pop
    time_total <- (total_vacc/daily_vaccine) + start_vaccination
    
    last_rt <- match(max(rt_profile$t[rt_profile$t<time_total]), rt_profile$t)
    
    if (last_rt <nrow(rt_profile))
    {
      
      rt_profile$Rt[last_rt+1] <- 1/(1-cov*efficacy)
      rt_profile$t[last_rt+1] <- time_total
      rt_profile$asDate[last_rt+1] <- as.Date(time_total, origin = day0)
      
      rt_profile <- rt_profile[-c((last_rt+2):length(rt_profile$t)), ]
    }
    
    Rts[[i+1]] <- rt_profile
  }
  return (Rts)
  
}

#### ---------------------  Running the model functions  ------------------------ ####

# Function that runs scenarios for each one of the coverages on the mat_antivax
# columns
# @param mat_antivax      Matrix with max. coverage per age group in each column
# @param rt_list          List with Rt profile for each column of mat_antivax or just one general profile
# @param eff              Infection efficacy
# @param eff_disease      Disease efficacy
# @param tt_vax           Vector with times of vaccination. 
run_scenarios <- function (mat_antivax, rt_list,eff, eff_disease, tt_vax, country= "United Kingdom",
                           max_vaccine,country_model=0 ){
  
  mat_antivax <- data.frame(mat_antivax)
  
  cov_list <-  vector(mode = "list", length = ncol(mat_antivax)-1)
  s <- vector(mode = "list", length = ncol(mat_antivax)-1)
  
  
  for (i in 1:length(cov_list)){
    
    # Adults 
    cov_list[[i]] <- priority_matrix_elder(mat_antivax[,i])  
    
    if (length(rt_list)==1){
      
      rt_prof <- rt_list[[1]]  
      
      s[[i]] <- list(R0 = rt_prof$Rt, tt_R0 = rt_prof$t, efficacy =eff,target = cov_list[[i]], 
                     efficacy_disease = eff_disease,
                     max_vaccine = max_vaccine, tt_vaccine = tt_vax, country=country)
    }
    if (length(rt_list) > 1){
      
      rt_prof <- rt_list[[i+1]]     # Scenarios from matrix start in position 2 
      
      s[[i]] <- list(R0 = rt_prof$Rt, tt_R0 = rt_prof$t, efficacy =eff,target = cov_list[[i]], 
                     efficacy_disease = eff_disease,
                     max_vaccine = max_vaccine, tt_vaccine = tt_vax,country=country)
    }
  }
  
  if(country_model ==0){
  
  runs <- lapply (s, function(x){
    
    do.call(run_scenario,x)
  })
  
  } else if (country_model ==1){
   
    runs <- lapply (s, function(x){

      do.call(run_country,x)
    })
  }
  
  
  return(runs)
  
}



# Function that summarises coverage matrix and get median, quantiles and ideal scenario coverage vectors
summary_coverages <- function (cov_matrix, adults= TRUE, ideal=0.95){
  
  age_vec <- c("0-5","5-10","10-15","15-20", "20-25","25-30","30-35", "35-40", "40-45", "45-50",
               "50-55", "55-60" ,"60-65", "65-70", "70-75", "75-80", "80+") 
  
  if (adults){
    
    cov = age_antivax(df = cov_matrix,ideal = ideal ,min_age = 15)
    
  }else {
    cov = age_antivax(df = cov_matrix,ideal = ideal ,min_age = 5)
  }
  cov$Age <- factor(cov$Age, levels= age_vec)
  cov_summary <- melt(cov[,-101],id.vars = "Age") %>% group_by(Age) %>% summarise_at(vars(value),list(Q1=~quantile(., probs = 0.10),
                                                                                                      median=median, Q3=~quantile(., probs = 0.90))) %>%
    mutate (ideal = cov$Ideal)  %>%
    relocate(Age, .after = last_col())
  
  return (cov_summary)
  
  
}



## Function that expands coverage per age group from the survey data to the
#  age groups required in NIMUE
# @param df       dataframe with coverage matrix from data survey sample 
# @param ideal    Ideal coverage used as counterfactual
# @param min_age  Mimimun age for vaccination: Values accepted 15 or 5

age_antivax <- function (df, ideal, min_age){
  
  age_groups <- c("0-5","5-10","10-15","15-20", "20-25","25-30","30-35", "35-40", "40-45", "45-50",
                  "50-55", "55-60" ,"60-65", "65-70", "70-75", "75-80", "80+")
  
  out <- data.frame(matrix(0, nrow= length(age_groups), ncol= ncol(df)-1))
  
  out$Age <- age_groups
  
  out$Ideal<-  ideal
  if(min_age == 15){
    
    out$Ideal[1:3] <- 0
    out[4:5,1:(ncol(df)-1)] <- df[1,2:ncol(df)]
    out[6:7,1:(ncol(df)-1)] <- df[2,2:ncol(df)]
    out[8:9,1:(ncol(df)-1)] <- df[3,2:ncol(df)]
    out[10:11,1:(ncol(df)-1)] <- df[4,2:ncol(df)]
    out[12:13,1:(ncol(df)-1)] <- df[5,2:ncol(df)]
    out[14:17,1:(ncol(df)-1)] <- df[6,2:ncol(df)]
  } else if (min_age == 5){
    out$Ideal[1] <- 0
    out[2:5,1:(ncol(df)-1)] <- df[1,2:ncol(df)]
    out[6:7,1:(ncol(df)-1)] <- df[2,2:ncol(df)]
    out[8:9,1:(ncol(df)-1)] <- df[3,2:ncol(df)]
    out[10:11,1:(ncol(df)-1)] <- df[4,2:ncol(df)]
    out[12:13,1:(ncol(df)-1)] <- df[5,2:ncol(df)]
    out[14:17,1:(ncol(df)-1)] <- df[6,2:ncol(df)]  
  }
  
  out <-out  %>% relocate (Age, .after = last_col())
  
  return(out)
  
}

## Prioritisation matrix 
## @param coverage Vector of 17 positions with the maximum coverage achieved by each age group 
priority_matrix_elder <- function(coverage) {
  
  total <- length(coverage)
  out <- matrix(0, nrow=total, ncol = total)
  
  counting <- 0
  
  for( i in 1:total){
    out[(total-counting),(counting+1):total] <- coverage[(counting+1):total]
    counting <- counting +1
  }
  
  return (out)
  
}

##Run scenarios
run_scenario <- function(R0 = c (3,0.9,2.4), tt_R0 = c(0,120 - 90, 365 - 90), target_pop = 50e6, 
                         income_group = "HIC", hs_constraints = "Present", country="United Kingdom",
                         efficacy = 0.95, efficacy_disease= 0, mode="Infection", target = matrix(0, nrow=1, ncol=17),
                         immunosenescence = 1, duration_R = 365, duration_V=Inf,
                         seeding_cases=20, max_vaccine = 1e8, tt_vaccine=1){
  
  # Population and mixing
  rep_country <- country
  pop <- squire::get_population(country = rep_country)$n
  pop_standardise <- target_pop / sum(pop)
  pop <- pop * pop_standardise
  mm <- squire::get_mixing_matrix(country = rep_country)
  # Hospital capacity
  hc <- get_capacity(country = rep_country, income_group = income_group, pop = pop, hs_constraints = hs_constraints)
  
  # Poorer health outcomes for LMICs and LICs
  pnsdt = get_prob_non_severe_death_treatment(income_group, hs_constraints)
  
  # Vaccine parameters
  # Efficacy
  efficacy <- set_efficacy(efficacy = efficacy, immunosenescence = immunosenescence, efficacy_disease = efficacy_disease )
  
  
  # Run
  #<3
  r1 <- nimue::run(
    time_period = (365 * 4) ,
    R0 = R0, 
    tt_R0 = tt_R0,
    population = pop,
    contact_matrix_set = mm,
    hosp_bed_capacity = hc$hosp_beds,
    ICU_bed_capacity = hc$ICU_beds,
    prob_non_severe_death_treatment = pnsdt,
    seeding_cases = seeding_cases,
    seed = 1,
    max_vaccine = max_vaccine,
    tt_vaccine = tt_vaccine,
    dur_V = duration_V,
    vaccine_efficacy_infection = efficacy$vaccine_efficacy_infection,
    vaccine_efficacy_disease = efficacy$vaccine_efficacy_disease,
    vaccine_coverage_mat = target,
    dur_R = duration_R
  )
  
  return(r1)
  
  
}

# Set hospital and ICU capacity
get_capacity <- function(country, income_group, pop, hs_constraints){
  
  hc <- squire::get_healthcare_capacity(country = country)
  
  # Unconstrained healthcare
  if(hs_constraints == "Absent"){
    hc$hosp_beds <- 1000000
    hc$ICU_beds <- 1000000
  }
  
  if(hs_constraints == "Present"){
    if(income_group %in% c("HIC", "UMIC")){
      hc$hosp_beds <- 1000000
      hc$ICU_beds <- 1000000
    }
    if(income_group %in% c("LMIC", "LIC")){
      hc$ICU_beds <- 0
    }
  }
  
  hc$hosp_beds <- round(hc$hosp_beds * sum(pop) / 1000)
  hc$ICU_beds <- round(hc$ICU_beds * sum(pop) / 1000)
  
  return(hc)
}


# Parameterise poorer health outcomes in LMIC and LIC
get_prob_non_severe_death_treatment <- function(income_group, hs_constraints){
  psdt <- squire:::probs$prob_non_severe_death_treatment
  
  if(income_group  == "LIC" & hs_constraints == "Present"){
    psdt <- c(rep(0.25, 16), 0.5804312)
  }
  return(psdt)
}

# Set vaccine efficacy against infection or disease
set_efficacy <- function(efficacy, immunosenescence, efficacy_disease){
  
  
  out <- list(vaccine_efficacy_infection=rep(efficacy, 17),
              vaccine_efficacy_disease=rep(efficacy_disease,17))
  
  
  out$vaccine_efficacy_infection[14:17] <- out$vaccine_efficacy_infection[14:17] * immunosenescence
  out$vaccine_efficacy_disease[14:17] <- out$vaccine_efficacy_disease[14:17] * immunosenescence
  
  return(out)
}



