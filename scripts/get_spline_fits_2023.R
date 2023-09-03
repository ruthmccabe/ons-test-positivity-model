### get median and uncertainty of all fits into one data frame 

library(tidyverse)
library(readxl)
library(cowplot)
library(matrixStats)

knot_props <- seq(20,40,10)
areas <- c("England","Scotland","Wales","Northern Ireland")
spline_degree <- 3

ons <- read_excel("data/ons_clean.xlsx",
                  sheet="positivity_long") %>%
  mutate(date_low = as.Date(date_low,"%d %B %Y"),
         date_high = as.Date(date_high,"%d %B %Y"),
         date = date_low+((date_high-date_low)/2),
         avg_pos = avg_pos/100,
         lower_pos = lower_pos/100,
         upper_pos = upper_pos/100) %>% 
  filter(date>=as.Date("2022-11-01"))

fits <- c()

for(area in areas){
  for(k in knot_props){
    
    knot_prop <- k/100
    
    fit <- readRDS(paste0("outputs/2023/objects/fit_",k,"_",area,"_2023.RDS")) 
    
    ff <- rstan::extract(fit)
    ff_exp <- exp(ff$Y_hat)/(exp(ff$Y_hat)+1)

    y_out <- t(ff_exp)

    ## try projecting onto whole thing to get more granular estimates 

    ons_area <- ons %>% filter(areaName==area)
    
    # ONS data for model
    X <- as.numeric(ons_area$date)
    
    # knots
    num_data <- length(X)
    knots <- c(min(X)-14,min(X)-7,X[seq(1, length(X), length.out = round(knot_prop*num_data,0))],max(X)+7,max(X)+14)
    num_knots <- length(knots)
    knots_dates <- as.Date(knots,origin="1970-01-01")
    
    ## try projecting onto whole thing to get more granular estimates 
    
    ons_daily <- ons_area %>% complete(date_low = seq.Date(min(date_low),max(date_high),by="day")) %>% 
      dplyr::select(-date,-date_high) %>% 
      fill(time_period,avg_pos,lower_pos,upper_pos,method,areaName,
           .direction="down") %>%
      rename(date = date_low) %>%
      filter(date>=min(ons_area$date),
             date<=max(ons_area$date))
    
    
    X_new <- seq(min(as.numeric(ons_daily$date))-14, max(as.numeric(ons_daily$date))+14, 0.1)
    X_daily <- as.numeric(ons_daily$date)
    B_true <- splines::bs(X_new, knots = knots, df = num_knots+spline_degree-1,
                          degree=spline_degree, intercept = FALSE)
    B_true <- t(predict(B_true, X_daily))
    Y_array <- array(data=NA, dim=c(nrow(ff$a), length(X_daily)))
    

    #a0<-mean(ff$a0)
    for(i in seq_len(nrow(ff$a))){
      a <- c(ff$a[i,],0)
      Y_array[i,] <- as.vector(a%*%B_true)
    }  

    Y_array_exp <- exp(Y_array)/(exp(Y_array)+1)
    y_out <- t(Y_array_exp)


    # output the spline fit
    df_out <- data.frame("date"=as.Date(ons_daily$date),"Y"=colQuantiles(Y_array_exp,probs=c(0.025,0.5,0.975))) %>% 
      mutate(knot_prop=k,areaName=area)

    fits <- rbind(fits,df_out)
    
  }
  
}


saveRDS(fits,"outputs/2023/objects/spline_fits_combined.RDS")




