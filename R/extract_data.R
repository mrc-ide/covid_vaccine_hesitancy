extract_data <- function (x, reduce_age = FALSE, compartment){
  
  
  # Indices of all the model outputs 
   indices <- squire:::odin_index(x$model)   
  # Convert cumulative data to incidence data 
  
  i_convert <-  unlist(indices[compartment]) 
  
  x$output[, i_convert, 1] <- apply(x$output[,i_convert, 1], 2, function(x){
    x - dplyr::lag(x)
  })
 
  index_D <- indices[compartment][[1]]               # Indices for the compartment D - Row: Age groups  Column: Vaccine status
  
  time <- x$output[,1,1]                    # Time vector 
  
  # 1. Extracting indexes row by row of the index matrix. 
  # 2. First collect simulation results of first vaccine status: Unvaccinated
  # 2.1 Then rowSums the other vaccine status to get an overall dynamics over time of: Vaccinated 
  # 3. o is a matrix where Columns: Age groups. nrows= length(time)x2. First half unvaccinated, second half vaccinated. 
  
  o <- apply(index_D, 1, 
             function(a, b){
               c(b[,a[1],1], rowSums(b[,a[2:6],1]))
             }, b = x$output)
  
  
  
  # Length of time vector per simulation 
  vec_t<- length(o[,1])/2
  
  
  # Transforming the first part of the data set: Unvaccinated 
  o1  <- o[c(1:vec_t),]
  if (reduce_age == FALSE) {
    
    o1 <- add_age(o1) 
  }
  else{
    
    o1 <- data.frame(collapse_age(o1))
  }
  
  o1 <- dplyr::bind_rows(o1) %>%  dplyr::mutate(t = rep(time, length(o1[,1])/ length(time)))
  o1$vacc_status <- "Unvaccinated"
  
  # Transforming the second part of the data set: Vaccinated 
  o2 <-  o[c((vec_t+1):(vec_t*2)),]
  
  if (reduce_age == FALSE) {
    
    o2 <- add_age(o2) 
  }
  else{
    
    o2 <- data.frame(collapse_age(o2))
  }
  o2 <-  dplyr::bind_rows(o2) %>%  dplyr::mutate(t = rep(time, length(o2[,1])/ length(time)))
  o2$vacc_status <- "Vaccinated"
  
  # Binding both part for output 
  out <- rbind(o1,o2)
  
  return(out)
  
}

# Function that converts matrix with columns of age groups to matrix of 2 columns: "Age_index", "age_group"
add_age <- function(x){
  m <- matrix(c(rep(1:ncol(x), each = (nrow(x))), as.vector(x)), ncol = 2)
  colnames(m) <- c("age_index", "value")
  
  m <- data.frame(m)
  
  ag <- c(paste0(seq(0, 75, 5), "-", seq(5, 80, 5)), "80+")
  m$age_group = factor(ag[m$age_index], levels = ag)
  m <- m  %>%  dplyr::select(-.data$age_index)
  
  return(m)
}

collapse_age <- function(x){
  m <- matrix(rowSums(x), ncol = 1)
  colnames(m) <- "value"
  return(m)
}




