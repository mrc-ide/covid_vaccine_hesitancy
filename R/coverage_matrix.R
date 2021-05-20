### -------------------------------------------- ####
##             Data importation                   ##
##  -------------------------------------------- ##

# Importing survey cleaned data 
data <- read.csv("./data/imperial_all3.0.csv")
# Make survey response as factor
data$vac_1 <- as.factor(data$vac_1)
data$Vaccine_fct <- factor(data$vac_1, levels = levels(data$vac_1),
                           labels = c(
                             "Strongly agree",
                             "Somewhat agree",
                             "Neutral/no opinion",
                             "Somewhat disagree",
                             "Strongly disagree"))

# Select European countries 
data <- data %>% filter(!is.na(data$agegroup_fct)) %>% 
  filter (!country %in% c("Australia","Canada","Japan","Singapore","South Korea"))   
# Grouping by age 
by_age <- data %>% group_by(agegroup_fct)
by_age <- by_age %>% count(Vaccine_fct)  %>% group_split()


france <- data %>% filter (country == "France")%>% group_by(agegroup_fct)%>% count(Vaccine_fct)  %>% group_split()
germany <- data %>% filter (country == "Germany")%>% group_by(agegroup_fct)%>% count(Vaccine_fct)  %>% group_split()
uk <- data %>% filter (country == "United Kingdom")%>% group_by(agegroup_fct)%>% count(Vaccine_fct)  %>% group_split()

### -------------------------------------------- ####
##      Coverage per age group estimate           ##
##  -------------------------------------------- ##

coverage_draws <- function (by_age){

# Number of uncertainty draws
n_draws <- 100
cov_matrix <- data.frame(matrix(0, ncol =n_draws +1, nrow=length(by_age)))

weighting_alpha <- c(39,30,20,10,1) # ((1 - mu) / var - 1 / mu) * mu ^ 2#
weighting_beta <-  40 - weighting_alpha #weighting_alpha * (1 / mu - 1)

for (i in 1:length(by_age)){

temp <- by_age[[i]]

# !For a single age-group
# Number of survey respondents
pop_size <- sum (temp$n)
# Survey responses
responses <- c("1"=temp$n[1] ,"2"=temp$n[2],"3"=temp$n[3],"4"=temp$n[4],"5"=temp$n[5])
# Draw responses from random multinomial
response_draws <-  rmultinom(n_draws, pop_size, prob = responses / pop_size)
# Weighting = the probability that someone would accept the vaccine given their response
weighting <- c("1" = 1, "2" = 0.25,"3"=0.5,"4"=0.75, "5" = 0)
# Work out the coverage for each draw
coverage <- apply(response_draws, 2, function(response_draw, weighting){
  sum(response_draw * weighting) / sum(response_draw)
}, weighting = weighting)


# Work out the coverage for each draw
#mu <- c(0.0001,0.25,0.5,0.75,1)
#var <- 0.00001
coverage2 <- apply(response_draws, 2, function(response_draw, weighting_alpha, weighting_beta){
  weighting_draw <- rbeta(5, weighting_alpha, weighting_beta)
  sum(response_draw * weighting_draw) / sum(response_draw)
},  weighting_alpha = weighting_alpha, weighting_beta = weighting_beta)
# Visualise the uncertainty
#hist(coverage2)


cov_matrix[i,1] <- temp$agegroup_fct[1]
cov_matrix[i,2:(n_draws+1)] <- coverage2


}

return (cov_matrix)

}


# weighting_alpha <- c(39,30,20,10,1) # ((1 - mu) / var - 1 / mu) * mu ^ 2#
# weighting_beta <-  40 - weighting_alpha #weighting_alpha * (1 / mu - 1)
# 
# x <- seq(0,1,length.out = 100)
# 
# a <- dbeta(x, weighting_alpha[1], weighting_beta[1])
# b <- dbeta(x, weighting_alpha[2], weighting_beta[2])
# c <- dbeta(x, weighting_alpha[3], weighting_beta[3])
# d <- dbeta(x, weighting_alpha[4], weighting_beta[4])
# e <- dbeta(x, weighting_alpha[5], weighting_beta[5])
# 
# df <- data.frame(x,a,b,c,d,e)
# df <- gather(df, func, val, -x)
# 
# gg <- ggplot(df, aes(x=x, y=val, group=func))+ geom_line(aes(color=func),size=1) +
#   scale_color_manual(name="Beta params",
#                      values=c("#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352"),
#                      labels=c("α=39, β=1", "α=30, β=10","α=20, β=20", "α=10, β=30","α=1, β=39"))+theme_bw()



