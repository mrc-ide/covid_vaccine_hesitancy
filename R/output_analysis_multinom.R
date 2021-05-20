# Function that binds Rt profiles from a list 
bind_profiles <- function(list_rt, scenario_list,vacc_start, vacc_finishes,day0){
  
  big_df<- data.frame(matrix(0,nrow=0, ncol= 4))
  
  for (i in 1:length(list_rt)){
    
    temp <- list_rt[[i]]
    
    end_plot_date <-  day0 + vacc_finishes + 365
    end_plot <-  vacc_finishes + 365
    
    last <- data.frame(temp$Rt[length(temp$Rt)], end_plot, end_plot_date)
    colnames(last) <- colnames(temp)
    temp<- rbind(temp,last)  
    
    temp$scenario <- scenario_list[i]
    
    big_df<- rbind(big_df, temp)
  } 
  
  return (big_df)
}

# Extracting data for one compartment. 
compartment_out <- function (out_nimue, pcompartment, reduce_age= TRUE){
  
  out_format <- nimue::format(out_nimue,reduce_age = reduce_age) %>% filter(compartment == pcompartment )
  return (out_format)
}


# Bind in a big data frame scenarios to compare
bind_to_compare <- function (list_out,compartment, reduce_age=TRUE){
  
  big_df <- list_out[[1]] %>% compartment_out(pcompartment = compartment,reduce_age = reduce_age) 
  
  for (i in 2:length(list_out)){
    
    df <- list_out[[i]] %>% compartment_out(pcompartment = compartment, reduce_age = reduce_age) 
    
    big_df <- rbind (big_df, df)
    
  }
  return (big_df)
  
}

# Disease dynamics summary. Extract dynamics for a compartmet for each run and return summar of stats per time step 
summary_runs <- function (runs_list, ini, fini, pcompartment,day0, scenario,reduce_age = TRUE){
  
  temp <- bind_to_compare(list_out = runs_list[c(ini:fini)], compartment = pcompartment,reduce_age = reduce_age )
  
  if (reduce_age == FALSE) {
    out <- temp %>% group_by(t,age_group) %>% summarize(avg = mean(value), n = n(), sd = sd(value), se = sd/sqrt(n)) %>%
      mutate(as_date= as.Date(t, origin = day0),
             lower_ci= avg - qt(1 - (0.05 / 2), n - 1) * se,
             upper_ci = avg + qt(1 - (0.05 / 2), n - 1) * se,
             scenario= scenario)
    
  }
  else {
    
    out <- temp %>% group_by(t) %>% summarize(avg = mean(value), n = n(), sd = sd(value), se = sd/sqrt(n)) %>%
      mutate(as_date= as.Date(t, origin = day0),
             lower_ci= avg - qt(1 - (0.05 / 2), n - 1) * se,
             upper_ci = avg + qt(1 - (0.05 / 2), n - 1) * se,
             scenario= scenario)
  }
  
  return(out)
}

# Function that gets runs from median, quantiles and ideal scenarios and returns a dataframe to plot dynamics
to_plot_output <- function (runs_list, runs_ch = NULL, names_vector,day0){
  
  adul <- compartment_out(out_nimue = runs_list[[1]],pcompartment = "deaths") %>% mutate(as_date= as.Date(t, origin = day0))
  colnames(adul)[1] <- names_vector[1]
  
  ch <-  compartment_out(out_nimue = runs_ch[[1]],pcompartment = "deaths") %>% mutate(as_date= as.Date(t, origin = day0))
  colnames(ch)[1] <- names_vector[1]
  
  for ( i in 2:length(runs_list)){
    
    #Adults
    temp <- compartment_out(out_nimue = runs_list[[i]],pcompartment = "deaths") %>% mutate(as_date= as.Date(t, origin = day0))
    colnames(temp)[1] <- names_vector[i]
    
    adul <- adul %>% left_join(temp)
    #Children 
    temp2 <- compartment_out(out_nimue = runs_ch[[i]],pcompartment = "deaths") %>% mutate(as_date= as.Date(t, origin = day0))
    colnames(temp2)[1] <- names_vector[i]
    
    ch <- ch %>% left_join(temp2)
    
  }
  
  
  adul$scenario <- "Adults"
  ch$scenario <- "Children"
  adul <- rbind(adul,ch)
  
  return(adul)
  
}


# Function that extracts cumulative values for vacc and anti vacc individuals 
# given start and finish date
extract_cumu <- function(df,start,finish){
  
  out <- data.frame (matrix(0,ncol=0,nrow=2))  
  out$vacc_status <- levels(as.factor(df$vacc_status))
  
  df_v <- df %>% filter (vacc_status == "Vaccinated", t >= start , t <= finish) 
  v <- sum(df_v$value)
  df_uv <- df %>% filter (vacc_status == "Unvaccinated", t >= start , t <= finish)
  uv <- sum(df_uv$value)
  
  out <- out %>% mutate(cum_value = ifelse(vacc_status == "Vaccinated",v,uv))
  
  return(out)
}

# Function that estimate average cumulative values for each vaccination status
# It iterates and summarises data from the different runs
cum_output_vac <- function (runs,ini,last,pcompartment,efficacy,start_vaccination,final_day, names){
  
  df<- data.frame(matrix(0, nrow=0, ncol= 3))
  
  for (i in ini:last){
    temp <- extract_data(x=runs[[i]],reduce_age = TRUE, compartment =  pcompartment) %>% 
      extract_cumu (start =start_vaccination,finish = final_day)
    
    temp$scenario <- names[i]
    
    df <- rbind(df,temp)
    
  }
  
  out <- df %>% mutate(efficacy =efficacy) #df %>% group_by(vacc_status) %>% summarize(avg = mean(cum_value), n = n(), sd = sd(cum_value), se = sd/sqrt(n)) %>%
  # mutate(lower_ci= avg - qt(1 - (0.05 / 2), n - 1) * se,
  #        upper_ci = avg + qt(1 - (0.05 / 2), n - 1) * se,
  #        efficacy = efficacy,
  #        scenario= scenario )
  # 
  return(out)
}





## ---------------------------------------------- ####
##            PLOTTING FUNCTIONS                   ##
## ---------------------------------------------   ##

## 1. 
plot_compartment_dynamics <- function (df, with_age= FALSE,yscale,ylab, day0 ){
  
  
  df$as_date <- as.Date(df$t,origin = day0)
  
  plot <- ggplot() + geom_line(data = df,aes(x= as_date, y = value/yscale), colour="brown4", size=0.8) + ylab(ylab) +
    theme_bw() + xlab("Time")
  
  if(with_age == TRUE){
    
    plot <- plot + facet_wrap(~age_group,scales = "free_y")
  }
  
  return(plot)
}


## 2.  Plot dynamics for different scenarios 
plot_dynamics_multi <- function(hesitancy_df,ideal_df, vacc_start, vacc_finishes, vacc_fin_chil, 
                                yscale, with_age = FALSE, ylabel, day0){
  
  pre_vac <- ideal_df %>% filter (t< vacc_start)
  
  colors <- c("#ba6d43","#824979","#030303")
  
  plot <- ggplot() + 
    geom_line(data = hesitancy_df, aes(x = as_date, y = avg/yscale,color = scenario), size = 1) + 
    scale_color_manual(name="",values=colors[1:2]) +
    geom_ribbon(data = hesitancy_df, aes(x = as_date, ymin = lower_ci/yscale, ymax = upper_ci/yscale, fill=scenario), alpha = .2)+
    scale_fill_manual(name = "", values=colors[1:2]) +
    geom_line (data=ideal_df,aes(x=as_date,y = value/yscale, color=scenario), size = 1, linetype = "dotted")+
    geom_line (data=pre_vac,aes(x=as_date,y = value/yscale), size = 1, color="grey42") +
    geom_vline(xintercept = day0+ vacc_start, linetype="dashed",color = "#696969", size=0.6)+
    geom_vline(xintercept = day0+ vacc_finishes,linetype="dashed", color = colors[1],  size=0.6)+
    geom_vline(xintercept = day0+  vacc_fin_chil , linetype="dashed",color = colors[2], size=0.6)+ 
    theme_bw()  + xlab("Time (days)") + ylab (ylabel)+
    scale_x_date(date_breaks = "6 month" , date_labels = "%b-%y", limits=c(as.Date("2020-03-01"),as.Date("2023-01-06"))) +
    theme(legend.position="none",
          strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line())
  
  if (with_age == TRUE){
    
    plot <- plot +  facet_wrap(~age_group,scales = "free_y")
    
  }  
  
  return(plot) 
}


# 3.  Plor Rt for Figure 1 
plot_rt_scenarios <- function(rt, start_vaccination, end_vaccination, end_children,day0){
  
  colors <- c("#ba6d43","#824979","#030303")
  
  end_plot_date <-  day0 + end_vaccination + 365
  end_plot <-  end_vaccination + 365
  
  x <- rt %>% group_by(scenario) %>%  filter(t == max(t)) 
  x$t <- end_plot
  x$asDate <- end_plot_date
  rt<- rbind(rt,x)  
  
  out<-  ggplot(rt, aes(x = asDate, y = Rt, color=scenario)) +
    geom_step(size=1) + theme_bw(base_family = 14)  + xlab("Time (days)") + ylab("Rt") + 
    scale_color_manual(values=colors)+
    geom_vline(xintercept = day0+ start_vaccination, linetype="dotted", color = "#696969", size=0.8)+
    geom_vline(xintercept = day0+ end_vaccination, color = colors[1], linetype="dotted", size=0.8)+
    geom_vline(xintercept = day0+ end_children, color = colors[2], linetype="dotted", size=0.8)+
    scale_x_date(date_breaks = "6 month" , date_labels = "%b-%y")
  #xlim(day0,end_plot_date)
  
  
  return (out)
  
}


# 3 Bar plots by Vacc status 

plot_bar_vacc_status <- function (df,yscale,ylab){
  
  df$scenario <- as.factor(df$scenario)
  
  x <- df %>% filter(scenario != "Ideal") %>% spread(scenario,cum_value)%>% rename(cum_value = Hesitancy) %>% 
    mutate(scenario = "Hesitancy")
  y <- df %>% filter (scenario == "Ideal") 
  
  toplot <- bind_rows(x,y)
  
  #df_ideal <- df %>% filter (scenario == "Ideal")
  # df_ideal <- df %>% filter (scenario == "Hesitncy")
  
  plot <- ggplot(data = toplot, aes(x = reorder(scenario, cum_value), y = cum_value/yscale, 
                                    fill = vacc_status))+
    geom_col() +
    geom_text(aes(label = round(cum_value/yscale,0)),position = position_stack(vjust = 0.8), size = 3.2) +
    geom_text(aes(label =ifelse(is.na(Lower), " ", paste("\n(",round(Upper/yscale,0), "-",round(Lower/yscale,0),")"))),position = position_stack(vjust = 0.75), size = 3.2)+
    labs(y = ylab, x = "Scenario", fill = "Vaccine status") +
    scale_fill_manual(values = c("#c59e96","#5c8e72")) +
    theme_bw(base_size = 13) +
    theme(#aspect.ratio = 1,
      strip.background = element_rect(fill = NA),
      panel.border = element_blank(),
      axis.line = element_line(),
      legend.text.align = 0) +
    facet_wrap(~efficacy)
  
  return(plot)
  
}


# 4  Plot rt profilie 

plot_rt <- function(rt, start_vaccination, end_vaccination,end_children,day0){
  
  rt <- rt[,1:3]
  colors <- c("#ba6d43","#824979","#030303")
  
  end_plot_date <-  day0 + end_vaccination + 365
  end_plot <-  start_vaccination + 740
  
  last <- data.frame(rt$Rt[length(rt$Rt)], end_plot, end_plot_date)
  colnames(last) <- colnames(rt)
  rt<- rbind(rt,last)  
  
  out<-  ggplot(rt, aes(x = asDate, y = Rt)) +
    geom_step(size=1)  + xlab("Time (days)") + ylab("Rt") + 
    geom_vline(xintercept = day0+ start_vaccination, linetype="dotted", color = "#696969", size=0.8)+
    geom_vline(xintercept = day0+ end_vaccination, color = colors[1], linetype="dotted", size=0.8)+
    geom_vline(xintercept = day0+ end_children, color = colors[2], linetype="dotted", size=0.8)+
    scale_x_date(date_breaks = "6 month" , date_labels = "%b-%y",limits=as.Date (c(NA,"2023-01-06")))+
    theme_bw(base_family = 12) +
    theme(legend.position="none",
          strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line())
  #xlim(day0,end_plot_date)
  
  
  return (out)
  
}

## 5 Plot dynamics for mean, ideal, upper and lower bounds 
plot_dynamics_short <- function(df, vacc_start, vacc_finishes, vacc_fin_chil, 
                                yscale, with_age = FALSE, ylabel, day0){
  library(ggnewscale)
  
  pre_vac <- df %>% filter (t< vacc_start)
  
  colors <- c("#ba6d43","#824979")
  col2 <- c("#2f3237","#a6b380")
  
  plot <- ggplot(df, aes(x= as_date, group=scenario)) + 
    geom_line(aes(y = Hesitancy/yscale,color = scenario), size = 1.1) + 
    scale_color_manual(values=colors, name="Vaccine hesitancy", labels= c("Adults","Children")) +
    geom_line (aes(y=Lower/yscale, color=scenario), alpha=0.7, linetype = "dashed",size=1) +
    geom_line (aes(y=Upper/yscale, color=scenario), alpha=0.7, linetype ="longdash",size=1) +
    new_scale_colour() +
    geom_line (aes(y = Ideal/yscale, color=scenario), size = 1.1)+
    scale_colour_manual(values=col2,name="Ideal Scenario", labels= c("Adults","Children"))+
    geom_vline(xintercept = day0+ vacc_start, linetype="dotted", color = "#696969", size=0.8)+
    geom_vline(xintercept = day0+ vacc_finishes, color = colors[1], linetype="dotted", size=0.8)+
    geom_vline(xintercept = day0+  vacc_fin_chil , color = colors[2], linetype="dotted", size=0.8)+ 
    theme_bw()  + xlab("Time (days)") + ylab (ylabel)+
    scale_x_date(date_breaks = "6 month" , date_labels = "%b-%y", limits=c(as.Date("2020-03-01"),as.Date("2023-01-06"))) +
    theme(legend.position = "bottom",
      strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line()) 
    #+ ylim(0,7)
  
  return(plot)
  
}


#6 
plot_bar_deaths_comp <- function (df){
  
  library(RColorBrewer)
  pal <-  colorRampPalette(brewer.pal(5, "BuPu"))
  
  
  df$country <- as.factor(df$country)
  #df$Country <- fct_rev(df$Country)
  
  df$age_group <- factor(df$age_group, levels = c("0-5","5-10","10-15","15-20","20-25",
                                                  "25-30","30-35","35-40","40-45","45-50",
                                                  "50-55", "55-60","60-65","65-70","70-75","75-80","80+")) 
  
  plot <-  ggplot(data = df, aes(x = age_group, y = ratio, fill=age_group,label = round(ratio,1))) +
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=ratio_upper, ymax=ratio_lower),width=.2) +
    #geom_text(vjust = 3, size = 4) +
    facet_grid(efficacy~country) +
    scale_fill_manual(values=pal(17)) +
    theme_bw(base_family = 12) +
    theme(axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position="none",
          strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line())+
    labs(y = "Cumulative deaths ratio", x = "Age Group")# + ylim (0,8)
  #coord_flip()
  
  return(plot)
}

#7 Plots coverage per age group 
plot_coverage <- function (cov_matrix){
  
  cov_long <- melt(cov_matrix,value_name= "X1")
  
  
  coverage_plot <- ggplot(cov_long, aes(value),group= X1) + #geom_density() + facet_grid(~X1)
    geom_histogram(aes(y=..density..),color= "#9A8B4F", fill="#eeeadd", position="identity") + 
    facet_wrap(~X1, scales="free_y") +
    theme_bw(base_family = 16) +
    theme(axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position="none",
          strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line())+
    labs(y = "Density ", x = "Vaccine uptake") + 
    ylim(c(0,NA))
  
  return (coverage_plot)
  
}


# 8 Plor disease dynamics for country and efficacy scenario  
plot_dynamics_countries <- function(big_df, vacc_start, vacc_finishes,
                                    yscale, with_age = FALSE, ylabel, day0){
  
  total <- big_df %>%  group_by(t) %>% select(t,Hesitancy,Lower,Upper,Ideal) %>%summarize_all(sum)%>% 
    mutate(as_date =as.Date(t, origin = day0))
  
  colors <- c("#6B5876","#030303","#ba6d43","#824979")
  
  plot <- ggplot(data=total) + 
    geom_line(aes(x = as_date, y = Hesitancy/yscale), size = 1, color= colors[1]) + 
    geom_ribbon(aes(x = as_date, ymin = Lower/yscale, ymax = Upper/yscale),fill=colors[1], alpha = .4)+
    geom_line (aes(x=as_date,y = Ideal/yscale), size = 1, color=colors[2])+
    #geom_line (data=pre_vac,aes(x=as_date,y = value/yscale), size = 1, color="grey42") +
    geom_vline(xintercept = day0+ vacc_start, linetype="dashed",color = "#696969", size=0.6)+
    geom_vline(xintercept = day0+ vacc_finishes,linetype="dashed", color = colors[3],  size=0.6)+
    theme_bw()  + xlab("Time (days)") + ylab (ylabel)+
    scale_x_date(date_breaks = "6 month" , date_labels = "%b-%y") +
    theme(legend.position="none",
          strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line())
  
  
  
  if (with_age == TRUE){
    
    plot <- ggplot(data=big_df) + 
      geom_line(aes(x = as_date, y = Hesitancy/yscale), size = 1, color= colors[1]) + 
      geom_ribbon(aes(x = as_date, ymin = Lower/yscale, ymax = Upper/yscale),color=colors[1], alpha = .4)+
      geom_line (aes(x=as_date,y = Ideal/yscale), size = 1, color=colors[2])+
      #geom_line (data=pre_vac,aes(x=as_date,y = value/yscale), size = 1, color="grey42") +
      geom_vline(xintercept = day0+ vacc_start, linetype="dashed",color = "#696969", size=0.6)+
      geom_vline(xintercept = day0+ vacc_finishes,linetype="dashed", color = colors[3],  size=0.6)+
      theme_bw()  + xlab("Time (days)") + ylab (ylabel)+
      scale_x_date(date_breaks = "6 month" , date_labels = "%b-%y")+ facet_wrap(~age_group,scales = "free_y")
    
  }  
  
  return(plot) 
}



# 9. Plot overall vaccine hesitancy given country and cov_matrix
plot_overal_hesitancy <- function(country_list, cov_matrix_list, ini, fini){
  
  mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(15)
  toplot <- data.frame(matrix(0,nrow = 0,ncol=5))
  
  for (i in 1:length(country_list)){
    
    country <- country_list[i]
    cov_matrix <- cov_matrix_list[[i]]
    
    pop <- squire::get_population(country = country)$n
    
    cov_matrix[,ini:fini] <- (cov_matrix[,ini:fini] *pop)/sum(pop)
    
    total <- data.frame(t(colSums(cov_matrix[,ini:fini])))
    total$country <- country
    
    toplot <- rbind(toplot,total)
    
    
  }
  
  toplot$country[which(toplot$country== "United Kingdom")] <- "U.K."
  
  coverage_plot <- ggplot(toplot, aes(x= country, fill=country))+
    geom_bar(aes(y = median),position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=Q1, ymax=Q3),
                  width=.2)+
    scale_fill_manual(values= mycolors, labels= c("France","Germany","U.K")) + 
    theme_bw(base_family = 14) +
    theme(legend.position="none",
          strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line())+
    labs(y = "Vaccine uptake ", x = "Country") + coord_flip() 
  
  return (coverage_plot)
}


# 10. 

plot_bar_deaths_all <- function (df){
  
  library(RColorBrewer)
  pal <-  colorRampPalette(brewer.pal(5, "BuGn"))
  
  
  df$country <- as.factor(df$country)
  #df$Country <- fct_rev(df$Country)
  
  df$age_group <- factor(df$age_group, levels = c("0-5","5-10","10-15","15-20","20-25",
                                                  "25-30","30-35","35-40","40-45","45-50",
                                                  "50-55", "55-60","60-65","65-70","70-75","75-80","80+")) 
  
  plot <-  ggplot(data = df, aes(x = age_group, fill=age_group)) +
    geom_bar(aes(y = cum_hes),position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=cum_lower, ymax=cum_upper),
                  width=.2) +
    scale_fill_manual(values=pal(17)) +
    facet_grid(country~.) +
    theme_bw(base_family = 12) +
    theme(axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position="none",
          strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line())+
    labs(y = "Cumulative deaths", x = "Age Group") 
  #coord_flip()
  
  return(plot)
}

plot_profile_2D <- function(big_df, vacc_start = "2021-01-01"){
  
  
  mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(15)
  plot <- ggplot(big_df, aes(x = asDate, y = Rt, color = country, linetype= efficacy))  +
    geom_step(size=0.8) + theme_bw() + 
    scale_linetype_manual(values=c("solid", "dotted"))+
    geom_vline(xintercept = as.Date(vacc_start), linetype="dotted", color = "#696969", size=0.8) + 
    geom_vline(xintercept = as.Date(vacc_start) + 304, linetype="dotted", color = "#696969", size=0.8)+
    scale_color_manual(values= mycolors) +
    scale_x_date(date_breaks = "6 month" , date_labels = "%b-%y") +
    labs(y = "Rt", x = "Time(days)", color = "Country", linetype = "Efficacy") +
    theme(strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line())
  
  return(plot) 
}


