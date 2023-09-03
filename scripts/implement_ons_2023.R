## write function that takes country as input and then runs to get spline fit, r and R

library(readxl)
library(tidyverse)
library(rstan)
library(splines)
library(matrixStats)
library(cowplot)
library(Hmisc)
library(parallel)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

ons <- read_excel("data/ons_clean.xlsx",
                  sheet="positivity_long") %>%
  mutate(date_low = as.Date(date_low,"%d %B %Y"),
         date_high = as.Date(date_high,"%d %B %Y"),
         date = date_low+((date_high-date_low)/2),
         avg_pos = avg_pos/100,
         lower_pos = lower_pos/100,
         upper_pos = upper_pos/100) # take date as mid period for plotting purposes


## filter to only be from the beginning of November 2022
ons <- ons %>% filter(date>=as.Date("2022-11-01"))

# ggplot(ons,aes(x=date,y=avg_pos))+
#   geom_point()+
#   facet_wrap(~areaName)+
#   scale_x_date(limits=c(as.Date("2022-11-01"),NA))

#stan model
spline_prev_model <-stan_model("stan_models/b_spline_ons_beta.stan")

gen_time <- readRDS("data/generation_time.RDS")

# ## testing purposes
# data <- ons
# area <- "Scotland"
# knot_prop <- 0.2
# spline_degree <- 3
# stan_model <- spline_prev_model
# iterations <- 20000
# warm <- 2000
# no_chains <- 4
# generation_time <- gen_time
# tau_max <- 14


get_estimates_2023 <- function(data,area,
                               knot_prop,
                               spline_degree = 3,stan_model,
                               iterations,warm,no_chains,
                               generation_time,
                               tau_max,
                               cl){
  #### for parallel computing
  library(readxl); library(tidyverse); library(rstan); library(splines); library(matrixStats); library(cowplot)
  
  # subset ONS data per country
  ons_area <- data %>% filter(areaName==area) %>%
    mutate(var = ((upper_pos - lower_pos)/3.92)^2,
           ons_sd = sqrt(var),
           alpha_ons = avg_pos*(((avg_pos * (1- avg_pos))/var) - 1),
           beta_ons = alpha_ons*((1 - avg_pos)/avg_pos))
  
  ons_area <- merge(ons_area,generation_time %>% rename(sd_gen = sd),by="date",all.x=TRUE)
  
  # ONS data for model
  X <- as.numeric(ons_area$date)
  Y_var <- ons_area$var
  alpha_ons = ons_area$alpha_ons
  beta_ons = ons_area$beta_ons
  
  # knots
  num_data <- length(X)
  knots <- c(min(X)-14,min(X)-7,X[seq(1, length(X), length.out = round(knot_prop*num_data,0))],max(X)+7,max(X)+14)
  num_knots <- length(knots)
  knots_dates <- as.Date(knots,origin="1970-01-01")
  
  fit<-sampling(stan_model,iter=iterations,
                control = list(adapt_delta=0.95,
                               max_treedepth = 10),
                warmup = warm,
                chains = no_chains,
                data = list(num_data = num_data,
                            num_knots = num_knots,
                            knots = knots,
                            spline_degree = spline_degree,
                            alpha_ons = alpha_ons,
                            beta_ons = beta_ons,
                            X = X))
  
  ff <- rstan::extract(fit)
  ff_exp <- exp(ff$Y_hat)/(exp(ff$Y_hat)+1)
  
  y_out <- t(ff_exp)
  
  ## try projecting onto whole thing to get more granular estimates 
  
  ons_daily <- ons_area %>% complete(date_low = seq.Date(min(date_low),max(date_high),by="day")) %>% 
    dplyr::select(-date,-date_high) %>% 
    fill(time_period,avg_pos,lower_pos,upper_pos,method,areaName,
         var,ons_sd,alpha_ons,beta_ons,central,sd_gen,gen_alpha,gen_beta,distribution,
         .direction="down") %>%
    rename(date = date_low)
  
  
  X_new <- seq(min(as.numeric(ons_daily$date))-14, max(as.numeric(ons_daily$date))+14, 0.1)
  X_daily <- ons_daily$date
  B_true <- splines::bs(X_new, knots = knots, df = num_knots+spline_degree-1,
                        degree=spline_degree, intercept = FALSE)
  B_true <- t(predict(B_true, X_daily))
  Y_array <- array(data=NA, dim=c(nrow(ff$a), length(X_daily)))
  
  #a0<-mean(ff$a0)
  for(i in seq_len(nrow(ff$a))){
    #a <- array(NA, num_knots+spline_degree-1)
    a <- c(ff$a[i,],0)
    #a0 <- ff$a0[i]
    # for(j in seq_len(length(a))){
    #   a[j] <- ff$a[i,j]
    #   #a[j] <- mean(ff$a[,j])
    # }
    Y_array[i,] <- as.vector(a%*%B_true)
  }  
  
  Y_array_exp <- exp(Y_array)/(exp(Y_array)+1)
  y_out <- t(Y_array_exp)
  
  
  # output the spline fit
  df_out <- data.frame("date"=as.Date(ons_daily$date),"Y"=colQuantiles(Y_array_exp,probs=c(0.025,0.5,0.975)))
  #df_out2 <- data.frame("date"=as.Date(ons_daily$date),"Y"=colQuantiles(Y_array_exp,probs=c(0.025,0.5,0.975)))
  
  lags <- seq(1:tau_max)
  lag_names <- paste("lag", formatC(lags, width = nchar(max(lags)), flag = "0"),
                     sep = "_")
  lag_functions <- setNames(paste("dplyr::lag(., ", lags, ")"), lag_names)
  
  
  make_df <- function(x){
    
    library(tidyverse)
    
    # create data frame with the spline fit from each column of the stan output
    data_it <- data.frame("date"=as.Date(ons_daily$date),"spline"=x,"distribution"=ons_daily$distribution,
                          "gen_alpha" = ons_daily$gen_alpha,"gen_beta"=ons_daily$gen_beta,
                          "gen_mean" = ons_daily$central,"gen_sd"=ons_daily$sd_gen)
    
    colnames(data_it)[2] <- "spline"
    
    data_it <- data_it %>%
      # mutate to the growth rate
      mutate(deriv = c(NA,diff(spline)/diff(as.numeric(date))),
             growth = deriv/spline) %>%
      # now we make columns for the two parameters depending on the relevant distribution
      mutate(param1 = ifelse(is.na(gen_alpha),gen_mean,gen_alpha),
             param2 = ifelse(is.na(gen_beta),gen_sd,gen_beta)) %>% #tibble() %>%
      group_by(date) %>%
      # get the relevant serial interval distribution
      mutate(dist = ifelse(distribution=="gamma",
                           list(dgamma(1:tau_max,shape=param1,rate=param2)),
                           ifelse(distribution=="lognormal",
                                  list(dlnorm(1:tau_max,meanlog=param1,sdlog=param2)),NA)),
             g_a = ifelse(distribution=="gamma",
                          sum(dgamma(1:tau_max,shape=param1,rate=param2)),
                          ifelse(distribution=="lognormal",
                                 sum(dlnorm(1:tau_max,meanlog=param1,sdlog=param2)),NA)))
    
    
    # now add the lags into the data
    data_it <- data_it %>% ungroup() %>% mutate_at(vars(spline), funs_(lag_functions)) %>%
      group_by(date) %>% #mutate(list(c(c_across(lag_01:lag_14)))) %>%
      #group_by(date) %>%
      mutate(values_to_mult = list(c(c_across(first(lag_names):last(lag_names))))) %>%
      #dplyr::select(-starts_with(lag_names)) %>%
      mutate(result = list(dist[[1]]*values_to_mult[[1]]),
             integral = sum(dist[[1]]*values_to_mult[[1]]),
             Rt = spline/(integral/g_a))
    
    return(data_it %>% dplyr::select(date,growth,Rt))
    
  }
  
  ## make y_out into a list for use with lapply
  y_out_list <- apply(y_out,2,list)
  
  
  # ####### add the cluster call thing here
  #library(parallel)
  # cl <- makeCluster(6)
  # clusterExport(cl = cl, c("make_df","y_out_list","ons_daily",
  #                          "lag_functions","lag_names","tau_max"))
  
  
  ## for each data frame put into a list and then bind together
  output_it <- bind_rows(parLapply(cl = cl,X = y_out_list, fun = make_df))
  
  epi_ests <- output_it %>% group_by(date) %>% summarise(date = mean(date),
                                                         r.2.5. = quantile(growth,0.025,na.rm=TRUE),
                                                         r.50. = quantile(growth,0.5,na.rm=TRUE),
                                                         r.97.5. = quantile(growth,0.975,na.rm=TRUE),
                                                         R.2.5. = quantile(Rt,0.025,na.rm=TRUE),
                                                         R.50. = quantile(Rt,0.5,na.rm=TRUE),
                                                         R.97.5. = quantile(Rt,0.975,na.rm=TRUE))
  
  saveRDS(epi_ests,paste0("outputs/2023/objects/r_R_summary_",
                          knot_prop*100,"_",area,"_2023.RDS"))
  
  #### add in the plot here
  panel_A <- ggplot()+
    theme_bw()+
    geom_ribbon(data=df_out,aes(x=date,ymin=Y.2.5.*100,ymax=Y.97.5.*100),alpha=0.8,fill="#077b8a")+
    geom_line(data=df_out,aes(x=date,y=Y.50.*100),lwd=1,col="#077b8a")+
    geom_point(data=ons_area,aes(x=date,y=avg_pos*100))+
    geom_errorbar(data=ons_area,aes(x=date,ymin=lower_pos*100,ymax=upper_pos*100))+
    labs(x="Date",y="ONS CIS Test \nPositivity (%)",tag="A")+
    scale_x_date(date_breaks = "1 months",date_labels = "%b-%Y")+
    scale_y_continuous(breaks=seq(0,12,2))
  
  panel_B <- ggplot()+
    theme_bw()+
    geom_ribbon(data=epi_ests,aes(x=date,ymin=r.2.5.,ymax=r.97.5.),alpha=0.8,fill="#077b8a")+
    geom_line(data=epi_ests,aes(x=date,y=r.50.),lwd=1,col="#077b8a")+
    geom_abline(slope=0,intercept=0,linetype="dashed")+
    labs(x="Date",y="Growth Rate",tag="B")+
    scale_x_date(date_breaks = "1 months",date_labels = "%b-%Y")
  
  panel_D <- ggplot()+
    theme_bw()+
    geom_ribbon(data=epi_ests,aes(x=date,ymin=R.2.5.,ymax=R.97.5.),alpha=0.8,fill="#077b8a")+
    geom_line(data=epi_ests,aes(x=date,y=R.50.),lwd=1,col="#077b8a")+
    geom_abline(slope=0,intercept=1,linetype="dashed")+
    labs(x="Date",y="Effective reproduction number",tag="C")+
    scale_x_date(date_breaks = "1 months",date_labels = "%b-%Y")
  
  plot_grid(panel_A,panel_B,panel_D,
            ncol=1,align="v",axis="lr")
  ggsave(paste0("outputs/2023/",area,"_",knot_prop*100,"_2023.png"),height=6,width=4)
  
  saveRDS(fit,paste0("outputs/2023/objects/fit_",
                     knot_prop*100,"_",area,"_2023.RDS"))
  
  
  return(list(spline_fit = df_out,
              epi_ests = epi_ests))
  
}



#library(parallel)
cl <- makeCluster(6)
clusterExport(cl = cl, c("get_estimates_2023","ons","spline_prev_model",
                         "gen_time"))

#### Scotland
scotland_0.2 <- get_estimates_2023(data = ons,area = "Scotland",
                                   knot_prop = 0.2,
                                   spline_degree = 3,stan_model = spline_prev_model,
                                   iterations = 20000,warm = 2000,
                                   no_chains = 4,
                                   tau_max = 14,
                                   generation_time = gen_time,
                                   #official_data = official,
                                   cl = cl)

scotland_0.3 <- get_estimates_2023(data = ons,area = "Scotland",
                                   knot_prop = 0.3,
                                   spline_degree = 3,stan_model = spline_prev_model,
                                   iterations = 20000,warm = 2000,
                                   no_chains = 4,
                                   tau_max = 14,
                                   generation_time = gen_time,
                                   #official_data = official,
                                   cl = cl)

scotland_0.4 <- get_estimates_2023(data = ons,area = "Scotland",
                                   knot_prop = 0.4,
                                   spline_degree = 3,stan_model = spline_prev_model,
                                   iterations = 20000,warm = 2000,
                                   no_chains = 4,
                                   tau_max = 14,
                                   generation_time = gen_time,
                                   #official_data = official,
                                   cl = cl)





#### England
england_0.2 <- get_estimates_2023(data = ons,area = "England",
                                  knot_prop = 0.2,
                                  spline_degree = 3,stan_model = spline_prev_model,
                                  iterations = 20000,warm = 2000,
                                  no_chains = 4,
                                  tau_max = 14,
                                  generation_time = gen_time,
                                  #official_data = official,
                                  cl = cl)

england_0.3 <- get_estimates_2023(data = ons,area = "England",
                                  knot_prop = 0.3,
                                  spline_degree = 3,stan_model = spline_prev_model,
                                  iterations = 20000,warm = 2000,
                                  no_chains = 4,
                                  tau_max = 14,
                                  generation_time = gen_time,
                                  #official_data = official,
                                  cl = cl)

england_0.4 <- get_estimates_2023(data = ons,area = "England",
                                  knot_prop = 0.4,
                                  spline_degree = 3,stan_model = spline_prev_model,
                                  iterations = 20000,warm = 2000,
                                  no_chains = 4,
                                  tau_max = 14,
                                  generation_time = gen_time,
                                  #official_data = official,
                                  cl = cl)





#### Wales
wales_0.2 <- get_estimates_2023(data = ons,area = "Wales",
                                knot_prop = 0.2,
                                spline_degree = 3,stan_model = spline_prev_model,
                                iterations = 20000,warm = 2000,
                                no_chains = 4,
                                tau_max = 14,
                                generation_time = gen_time,
                                #official_data = official,
                                cl = cl)

wales_0.3 <- get_estimates_2023(data = ons,area = "Wales",
                                knot_prop = 0.3,
                                spline_degree = 3,stan_model = spline_prev_model,
                                iterations = 20000,warm = 2000,
                                no_chains = 4,
                                tau_max = 14,
                                generation_time = gen_time,
                                #official_data = official,
                                cl = cl)

wales_0.4 <- get_estimates_2023(data = ons,area = "Wales",
                                knot_prop = 0.4,
                                spline_degree = 3,stan_model = spline_prev_model,
                                iterations = 20000,warm = 2000,
                                no_chains = 4,
                                tau_max = 14,
                                generation_time = gen_time,
                                #official_data = official,
                                cl = cl)




#### Northern Ireland
NorthernIreland_0.2 <- get_estimates_2023(data = ons,area = "Northern Ireland",
                                          knot_prop = 0.2,
                                          spline_degree = 3,stan_model = spline_prev_model,
                                          iterations = 20000,warm = 2000,
                                          no_chains = 4,
                                          tau_max = 14,
                                          generation_time = gen_time,
                                          #official_data = official,
                                          cl = cl)

NorthernIreland_0.3 <- get_estimates_2023(data = ons,area = "Northern Ireland",
                                          knot_prop = 0.3,
                                          spline_degree = 3,stan_model = spline_prev_model,
                                          iterations = 20000,warm = 2000,
                                          no_chains = 4,
                                          tau_max = 14,
                                          generation_time = gen_time,
                                          #official_data = official,
                                          cl = cl)

NorthernIreland_0.4 <- get_estimates_2023(data = ons,area = "Northern Ireland",
                                          knot_prop = 0.4,
                                          spline_degree = 3,stan_model = spline_prev_model,
                                          iterations = 20000,warm = 2000,
                                          no_chains = 4,
                                          tau_max = 14,
                                          generation_time = gen_time,
                                          #official_data = official,
                                          cl = cl)



stopCluster(cl)
