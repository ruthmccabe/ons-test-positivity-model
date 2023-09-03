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

ons <- ons %>% filter(date<=as.Date("2022-12-31"))


#stan models 
spline_prev_model <-stan_model("stan_models/b_spline_ons_beta.stan")


official <- read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&metric=transmissionRateGrowthRateMax&metric=transmissionRateMax&metric=transmissionRateMin&metric=transmissionRateGrowthRateMin&format=csv") %>%
  mutate(date = as.Date(date)) %>% 
  mutate(transmissionRateGrowthRateMax = ifelse(is.na(transmissionRateGrowthRateMax),100,
                                                transmissionRateGrowthRateMax),
         transmissionRateGrowthRateMin = ifelse(is.na(transmissionRateGrowthRateMin),100,
                                                transmissionRateGrowthRateMin),
         transmissionRateMax = ifelse(is.na(transmissionRateMax),100,transmissionRateMax),
         transmissionRateMin = ifelse(is.na(transmissionRateMin),100,transmissionRateMin)) %>%
  group_by(areaName) %>% 
  mutate(reporting_group = seq(1,length(areaName),1)) %>%
  complete(date = seq.Date(min(date),max(date),by="day")) %>%
  fill(areaCode,areaName,areaType,
       transmissionRateGrowthRateMax,transmissionRateGrowthRateMin,
       transmissionRateMax,transmissionRateMin,
       reporting_group,
       .direction = "up") %>% 
  data.frame() %>%
  mutate(transmissionRateGrowthRateMax = ifelse(transmissionRateGrowthRateMax==100,NA,
                                                transmissionRateGrowthRateMax),
         transmissionRateGrowthRateMin = ifelse(transmissionRateGrowthRateMin==100,NA,
                                                transmissionRateGrowthRateMin),
         transmissionRateMax = ifelse(transmissionRateMax==100,NA,transmissionRateMax),
         transmissionRateMin = ifelse(transmissionRateMin==100,NA,transmissionRateMin)) %>%
  mutate(transmissionRateGrowthRateMid = (transmissionRateGrowthRateMax+transmissionRateGrowthRateMin)/2,
         transmissionRateMid = (transmissionRateMax+transmissionRateMin)/2) 



gen_time <- readRDS("data/generation_time.RDS")

# ## testing purposes
# data <- ons
# area <- "England"
# knot_prop <- 0.2
# spline_degree <- 3
# stan_model <- spline_prev_model
# iterations <- 5000
# warm <- 2000
# no_chains <- 4
# generation_time <- gen_time
# tau_max <- 14


get_estimates_sens <- function(data,area,
                               knot_prop,
                               spline_degree = 3,stan_model,
                               iterations,warm,no_chains,
                               generation_time,
                               tau_max,
                               official_data,
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
  
  ########## spline fit ##########
  
  
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
  
  saveRDS(fit,paste0("outputs/objects/fit_",
                     knot_prop*100,"_",area,".RDS"))
  
  fit_extract <- rstan::extract(fit)
  Yhat_exp <- exp(fit_extract$Y_hat)/(exp(fit_extract$Y_hat)+1)
  
  y_out <- t(Yhat_exp)
    
  ### save model convergence metrics
  fit_convergence <- data.frame("parameter"=rownames(summary(fit)$summary),summary(fit)$summary) 
  rownames(fit_convergence) <- c()
  saveRDS(fit_convergence,paste0("outputs/objects/fit_convergence_",
                        knot_prop*100,"_",area,".RDS"))
  
  
  ########## projection for daily granularity ##########
  
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
  B_true_proj <- t(predict(B_true, X_daily))
  B_true_proj <- B_true_proj[1:ncol(fit_extract$a),]
  Y_daily <- array(data=NA, dim=c(nrow(fit_extract$a), length(X_daily)))
  
  for(i in seq_len(nrow(fit_extract$a))){
    a <- fit_extract$a[i,]
    Y_daily[i,] <- as.vector(a%*%B_true_proj)
  }  
  
  Y_daily_exp <- exp(Y_daily)/(exp(Y_daily)+1)
  y_out <- t(Y_daily_exp)
  
  
  # output the spline fit
  df_out <- data.frame("date"=as.Date(ons_daily$date),
                       "Y"=colQuantiles(Y_daily_exp,probs=c(0.025,0.05,0.5,0.95,0.975))) %>%
    filter(date>=min(ons_area$date),date<=max(ons_area$date))
  
  saveRDS(df_out,paste0("outputs/objects/spline_fit_daily_",
                          knot_prop*100,"_",area,".RDS"))
  
  
  
  ## check that the daily and original fits agree
  df_original <- data.frame("date"=as.Date(ons_area$date),"Y"=colQuantiles(Yhat_exp,
                                                                       probs=c(0.025,0.05,0.5,0.95,0.975)))
  
  ggplot()+
    theme_bw()+
    geom_line(data=df_original, aes(x=date,y=Y.50.,col="Original"),lwd=1)+
    geom_line(data=df_out, aes(x=date,y=Y.50.,col="Daily"))
  ggsave(paste0("outputs/check_daily/",area,"_",knot_prop*100,".png"))
  
  
  ########## estimating R(t) and r(t) ##########
  
  
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
  

  # # ###### add the cluster call thing here
  # library(parallel)
  # cl <- makeCluster(6)
  # clusterExport(cl = cl, c("make_df","y_out_list","ons_daily",
  #                          "lag_functions","lag_names","tau_max"))

  
  ## for each data frame put into a list and then bind together
  output_it <- bind_rows(parLapply(cl = cl,X = y_out_list, fun = make_df))
  
  epi_ests <- output_it %>% group_by(date) %>% summarise(date = mean(date),
                                                         r.2.5. = quantile(growth,0.025,na.rm=TRUE),
                                                         r.5. = quantile(growth,0.05,na.rm=TRUE),
                                                         r.50. = quantile(growth,0.5,na.rm=TRUE),
                                                         r.95. = quantile(growth,0.95,na.rm=TRUE),
                                                         r.97.5. = quantile(growth,0.975,na.rm=TRUE),
                                                         R.2.5. = quantile(Rt,0.025,na.rm=TRUE),
                                                         R.5. = quantile(Rt,0.05,na.rm=TRUE),
                                                         R.50. = quantile(Rt,0.5,na.rm=TRUE),
                                                         R.95. = quantile(Rt,0.95,na.rm=TRUE),
                                                         R.97.5. = quantile(Rt,0.975,na.rm=TRUE)) %>%
    filter(date>=min(ons_area$date),date<=max(ons_area$date))
  
  saveRDS(epi_ests,paste0("outputs/objects/r_R_summary_",
                          knot_prop*100,"_",area,".RDS"))
  
  
  ########## comparison with official estimates ##########
  
  official_area <- official %>% filter(areaName==area) 
  
  official_comb <- merge(epi_ests,official_area,all.x=TRUE)
  
  #### add in the plot here
  panel_A <- ggplot()+
    theme_bw()+
    geom_ribbon(data=df_out,aes(x=date,ymin=Y.2.5.*100,ymax=Y.97.5.*100),alpha=0.8,fill="#077b8a")+
    geom_line(data=df_out,aes(x=date,y=Y.50.*100),lwd=1,col="#077b8a")+
    geom_point(data=ons_area,aes(x=date,y=avg_pos*100))+
    geom_errorbar(data=ons_area,aes(x=date,ymin=lower_pos*100,ymax=upper_pos*100))+
    labs(x="Date",y="ONS CIS Test \nPositivity (%)",tag="A")+
    scale_x_date(date_breaks = "3 months",date_labels = "%b-%Y")+
    scale_y_continuous(breaks=seq(0,12,2))
  
  panel_B <- ggplot()+
    theme_bw()+
    geom_ribbon(data=official_comb,
                aes(x=date,ymin=transmissionRateGrowthRateMin/100,ymax=transmissionRateGrowthRateMax/100),
                alpha=0.55,
                #fill="#b1f3fb"
                fill="#d72631")+
    geom_ribbon(data=official_comb,aes(x=date,ymin=r.2.5.,ymax=r.97.5.),alpha=0.8,fill="#077b8a")+
    geom_line(data=official_comb,aes(x=date,y=r.50.),lwd=1,col="#077b8a")+
    geom_abline(slope=0,intercept=0,linetype="dashed")+
    labs(x="Date",y="Growth Rate",tag="B")+
    scale_x_date(date_breaks = "5 months",date_labels = "%b-%Y")
  
  panel_C <- ggplot(official_comb)+
    theme_bw()+
    geom_point(aes(x=transmissionRateGrowthRateMid/100,y=r.50.),col="#a2d5c6",pch=1)+
    #geom_point(aes(x=transmissionRateGrowthRateMax/100,y=r.97.5.,col="Upper"))+
    geom_abline(slope=1,intercept=0,linetype="dashed")+
    geom_smooth(aes(x=transmissionRateGrowthRateMid/100,y=r.50.),method="lm",fill="#a2d5c6",col="#a2d5c6")+
    #geom_smooth(aes(x=transmissionRateGrowthRateMax/100,y=r.97.5.,col="Upper",fill="Upper"),method="lm")+
    # scale_color_manual(values=c("#edca82","#a9c0a6"))+
    # scale_fill_manual(values=c("#edca82","#a9c0a6"))+
    labs(x="Officially reported estimates",y="Modelled estimates",col="Bound",fill="Bound",tag="C")+
    theme(legend.position = "none")
  
  panel_D <- ggplot()+
    theme_bw()+
    geom_ribbon(data=official_comb,
                aes(x=date,ymin=transmissionRateMin,ymax=transmissionRateMax),alpha=0.55,fill="#d72631")+
    geom_ribbon(data=official_comb,aes(x=date,ymin=R.2.5.,ymax=R.97.5.),alpha=0.8,fill="#077b8a")+
    geom_line(data=official_comb,aes(x=date,y=R.50.),lwd=1,col="#077b8a")+
    geom_abline(slope=0,intercept=1,linetype="dashed")+
    labs(x="Date",y="Effective reproduction number",tag="D")+
    scale_x_date(date_breaks = "5 months",date_labels = "%b-%Y")
  
  panel_E <- ggplot(official_comb)+
    theme_bw()+
    geom_point(aes(x=transmissionRateMin,y=R.50.),col="#a2d5c6",pch=1)+
    #geom_point(aes(x=transmissionRateMax,y=R.97.5.,col="Upper"))+
    geom_abline(slope=1,intercept=0,linetype="dashed")+
    geom_smooth(aes(x=transmissionRateMin,y=R.50.),method="lm",fill="#a2d5c6",col="#a2d5c6")+
    #geom_smooth(aes(x=transmissionRateMax,y=R.97.5.,col="Upper",fill="Upper"),method="lm")+
    # scale_color_manual(values=c("#edca82","#a9c0a6"))+
    # scale_fill_manual(values=c("#edca82","#a9c0a6"))+
    labs(x="Officially reported estimates",y="Modelled estimates",col="Bound",fill="Bound",tag="E")+
    theme(legend.position = "none")
  
  plot_grid(panel_A,
            plot_grid(panel_B,panel_C,align="hv",rel_widths = c(1.5,1)),
            plot_grid(panel_D,panel_E,align="hv",rel_widths = c(1.5,1)),
            ncol=1,rel_heights=c(1,1.5,1.5))
  ggsave(paste0("outputs/",area,"_",knot_prop*100,".png"),height=7,width=8)
  
  
  ########## correlation, proportion agreement and proportion of variation ##########
  
  correlation_lags <- seq(-30,30,1)
  
  cor_outputs <- c()
  prop_agreement_outputs <- c()
  #prop_variation <- c()
  lm_results <- c()
  
  
  for(i in 1:length(correlation_lags)){
    
    ## correlation
    
    cor_outputs_it <- data.frame("correlation" = 
                                   c(round(cor(Lag(official_comb$transmissionRateGrowthRateMid,
                                                   correlation_lags[i]),
                                               official_comb$r.50.,
                                               use="complete",method="pearson"),2),
                                     round(cor(Lag(official_comb$transmissionRateMid,
                                                   correlation_lags[i]),
                                               official_comb$R.50.,
                                               use="complete",method="pearson"),2)),
                                 "metric" = c("r","R")) %>% 
      mutate("knot_prop"=knot_prop,
             "areaName"=area,
             "lag"=correlation_lags[i])
    
    cor_outputs <- rbind(cor_outputs,cor_outputs_it)
    
    ## proportion agreement
    pa_data_it <- rbind(official_comb %>% mutate(Lag = Lag(official_comb$transmissionRateGrowthRateMid,
                                                           correlation_lags[i])) %>%
                          dplyr::select(Lag,r.50.) %>% 
                          rename(official = Lag,
                                 est = r.50.) %>%
                          mutate(agreement_label = factor(ifelse(est>=0&official>=0,"Both positive",
                                                                 ifelse(est<0 & official<0,"Both negative",
                                                                        ifelse(est>=0&official<0,"Should be neg",
                                                                               ifelse(est<0 & official>=0,
                                                                                      "Should be pos",NA))))),
                                 metric = "r"),
                        official_comb %>% mutate(Lag = Lag(official_comb$transmissionRateMid,
                                                           correlation_lags[i])) %>%
                          dplyr::select(Lag,R.50.) %>%
                          rename(official = Lag,
                                 est = R.50.) %>%
                          mutate(agreement_label = factor(ifelse(est>=1&official>=1,"Both positive",
                                                                 ifelse(est<1 & official<1,"Both negative",
                                                                        ifelse(est>=1&official<1,"Should be neg",
                                                                               ifelse(est<1 & official>=1,
                                                                                      "Should be pos",NA))))),
                                 metric = "R")
    ) %>% 
      group_by(metric) %>%
      summarise("both_positive"=length(which(agreement_label=="Both positive")),
                "both_negative"=length(which(agreement_label=="Both negative")),
                "should_be_pos"=length(which(agreement_label=="Should be pos")),
                "should_be_neg"=length(which(agreement_label=="Should be neg")),
                "total_obs"=length(!is.na(agreement_label)),
                "prop_agreement"=(both_positive+both_negative)/total_obs) %>%
      mutate("knot_prop"=knot_prop,
             "areaName"=area,
             "lag"=correlation_lags[i])
    
    prop_agreement_outputs <- rbind(prop_agreement_outputs,pa_data_it)
    
    
    # ## proportion of variation 
    # prop_variation_it <- rbind(
    #   official_comb %>% mutate(Lag = Lag(official_comb$transmissionRateGrowthRateMid,
    #                                      correlation_lags[i]),
    #                            Lag = Lag/100) %>%
    #     dplyr::select(Lag,r.50.) %>% 
    #     rename(official = Lag,
    #            est = r.50.) %>%
    #     filter(!is.na(official),!is.na(est)) %>%
    #     mutate(num = (official - est)^2,
    #            mean = mean(official,na.rm=TRUE),
    #            denom = (official - mean)^2) %>%
    #     summarise("total_numerator" = sum(num,na.rm = TRUE),
    #               "total_denominator" = sum(denom,na.rm = TRUE),
    #               "total" = 1 - (total_numerator/total_denominator),
    #               "total_num/denom" = sum(num/denom,na.rm = TRUE),
    #               "total_1_num/denom" = sum(1-(num/denom),na.rm = TRUE)) %>%
    #     mutate(metric = "r"),
    #   official_comb %>% mutate(Lag = Lag(official_comb$transmissionRateMid,correlation_lags[i])) %>%
    #     dplyr::select(Lag,R.50.) %>% 
    #     rename(official = Lag,
    #            est = R.50.) %>%
    #     filter(!is.na(official),!is.na(est)) %>%
    #     mutate(num = (official - est)^2,
    #            mean = mean(official,na.rm=TRUE),
    #            denom = (official - mean)^2) %>%
    #     summarise("total_numerator" = sum(num,na.rm = TRUE),
    #               "total_denominator" = sum(denom,na.rm = TRUE),
    #               "total" = 1 - (total_numerator/total_denominator),
    #               "total_num/denom" = sum(num/denom,na.rm = TRUE),
    #               "total_1_num/denom" = sum(1-(num/denom),na.rm = TRUE)) %>%
    #     mutate(metric = "R")) %>%
    #   mutate("knot_prop"=knot_prop,
    #          "areaName"=area,
    #          "lag"=correlation_lags[i])
    # 
    # 
    # prop_variation <- rbind(prop_variation,prop_variation_it)
    
    ## linear regression analysis
    
    if(area %in% c("England","Scotland","Wales")){
      
      lm_daily_r_it <- lm(data = official_comb,
                          Lag(transmissionRateGrowthRateMid,correlation_lags[i])/100~r.50.)
      
      lm_daily_r_it_no_int <- lm(data = official_comb,
                          Lag(transmissionRateGrowthRateMid,correlation_lags[i])/100~-1+r.50.) 
      
      lm_daily_R_it <- lm(data = official_comb,
                        Lag(transmissionRateMid,correlation_lags[i])~R.50.)
    
      lm_daily_R_it_no_int <- lm(data = official_comb,
                        Lag(transmissionRateMid,correlation_lags[i])~-1+R.50.)
      
      
      lm_results_it <- data.frame("R2"=c(summary(lm_daily_r_it)$r.squared,
                                         summary(lm_daily_r_it_no_int)$r.squared,
                                         summary(lm_daily_R_it)$r.squared,
                                         summary(lm_daily_R_it_no_int)$r.squared),
                                  "R2_adj"=c(summary(lm_daily_r_it)$adj.r.squared,
                                             summary(lm_daily_r_it_no_int)$adj.r.squared,
                                             summary(lm_daily_R_it)$adj.r.squared,
                                             summary(lm_daily_R_it_no_int)$adj.r.squared),
                                  "coefficient"=c(summary(lm_daily_r_it)$coefficient[2,1],
                                                  summary(lm_daily_r_it_no_int)$coefficient[1],
                                                  summary(lm_daily_R_it)$coefficient[2,1],
                                                  summary(lm_daily_R_it_no_int)$coefficient[1]),
                                  "metric" = c("Instantaeous growth rate","Instantaeous growth rate",
                                               "Effective reproduction number","Effective reproduction number"),
                                  "intercept" = rep(c("Intercept","No Intercept"),2)) %>%
        mutate(knot_prop = knot_prop,
               areaName = area,
               lag = correlation_lags[i])
      
      lm_results <- rbind(lm_results,lm_results_it)
      

    }
    
    if(area=="Northern Ireland"){
      
      # lm_daily_r_it <- lm(data = official_comb %>% filter(areaName==area,knot_prop==k/100),
      #                     Lag(transmissionRateGrowthRateMid,correlation_lags[i])~r.50.)
      # 
      # lm_daily_r_it_no_int <- lm(data = official_comb %>% filter(areaName==area,knot_prop==k/100),
      #                            Lag(transmissionRateGrowthRateMid,correlation_lags[i])~-1+r.50.) 
      
      lm_daily_R_it <- lm(data = official_comb,
                          Lag(transmissionRateMid,correlation_lags[i])~R.50.)
      
      lm_daily_R_it_no_int <- lm(data = official_comb,
                                 Lag(transmissionRateMid,correlation_lags[i])~-1+R.50.)
      
      
      lm_results_it <- data.frame("R2"=c(NA,
                                         NA,
                                         summary(lm_daily_R_it)$r.squared,
                                         summary(lm_daily_R_it_no_int)$r.squared),
                                  "R2_adj"=c(NA,
                                             NA,
                                             summary(lm_daily_R_it)$adj.r.squared,
                                             summary(lm_daily_R_it_no_int)$adj.r.squared),
                                  "coefficient"=c(NA,
                                                  NA,
                                                  summary(lm_daily_R_it)$coefficient[2,1],
                                                  summary(lm_daily_R_it_no_int)$coefficient[1]),
                                  "metric" = c("Instantaeous growth rate","Instantaeous growth rate",
                                               "Effective reproduction number","Effective reproduction number"),
                                  "intercept" = rep(c("Intercept","No Intercept"),2)) %>%
        mutate(knot_prop = knot_prop,
               areaName = area,
               lag = correlation_lags[i])
      
      lm_results <- rbind(lm_results,lm_results_it)
      
      
    }
    
    
  }
  

  saveRDS(cor_outputs,paste0("outputs/objects/correlation_",
                             knot_prop*100,"_",area,".RDS"))
  
  saveRDS(prop_agreement_outputs,paste0("outputs/objects/prop_agreement_",
                             knot_prop*100,"_",area,".RDS"))
  
  # saveRDS(prop_variation,paste0("outputs/objects/prop_variation_",
  #                               knot_prop,"_",area,".RDS"))
  
  saveRDS(lm_results,paste0("outputs/objects/lm_results_",
                                knot_prop*100,"_",area,".RDS"))
  
  
  
  ########## aggregated estimates over time ##########

  output_it_grouped_summary <- merge(output_it %>% filter(date <= max(ons_area$date),
                                                          date >= min(ons_area$date)),
                                     official_comb %>% dplyr::select(date,reporting_group),
                                     by = "date",
                                   all.x=TRUE) %>% 
    filter(!is.na(reporting_group)) %>%
    group_by(reporting_group) %>%
    summarise(date = max(date),
              r.2.5. = quantile(growth,0.025,na.rm=TRUE),
              r.5. = quantile(growth,0.05,na.rm=TRUE),
              r.50. = quantile(growth,0.5,na.rm=TRUE),
              r.95. = quantile(growth,0.95,na.rm=TRUE),
              r.97.5. = quantile(growth,0.975,na.rm=TRUE),
              R.2.5. = quantile(Rt,0.025,na.rm=TRUE),
              R.5. = quantile(Rt,0.05,na.rm=TRUE),
              R.50. = quantile(Rt,0.5,na.rm=TRUE),
              R.95. = quantile(Rt,0.95,na.rm=TRUE),
              R.97.5. = quantile(Rt,0.975,na.rm=TRUE))
  
  saveRDS(output_it_grouped_summary,paste0("outputs/objects/grouped_summary_",
                             knot_prop*100,"_",area,".RDS"))
  
  
  ### plot this and also compare to the official estimates 

 output_summary_comb <- rbind(output_it_grouped_summary %>% 
                                dplyr::select(reporting_group,date,
                                              r.5.,r.50.,r.95.,
                                              R.5.,R.50.,R.95.) %>%
                                rename(r_lower = r.5.,
                                       r_med = r.50.,
                                       r_upper = r.95.,
                                       R_lower = R.5.,
                                       R_med = R.50.,
                                       R_upper = R.95.) %>% 
                                mutate(estimate = "Modelled"),
                               merge(official_area %>% dplyr::select(reporting_group,
                                 transmissionRateGrowthRateMin,
                                 transmissionRateGrowthRateMid,
                                 transmissionRateGrowthRateMax,
                                 transmissionRateMin,
                                 transmissionRateMid,
                                 transmissionRateMax) %>%
   rename(r_lower = transmissionRateGrowthRateMin,
          r_med = transmissionRateGrowthRateMid,
          r_upper = transmissionRateGrowthRateMax,
          R_lower = transmissionRateMin,
          R_med = transmissionRateMid,
          R_upper = transmissionRateMax) %>%
   mutate(r_lower = r_lower/100,
          r_med = NA,
          r_upper = r_upper/100,
          R_med = NA) %>% 
     unique(),
   output_it_grouped_summary %>% dplyr::select(reporting_group,date),by="reporting_group",all.x=TRUE) %>%
     mutate(estimate = "Government published") %>%
     dplyr::select(reporting_group,date,r_lower,r_med,r_upper,R_lower,R_med,R_upper,estimate)) %>%
   mutate(r_midpoint = ifelse(estimate=="Government published",(r_upper + r_lower)/2,r_med),
          R_midpoint = ifelse(estimate=="Government published",(R_upper + R_lower)/2,R_med))
 
 ### combine into one plot
 plot_grid(ggplot(output_summary_comb,
        aes(x=date,ymin=r_lower,ymax=r_upper,y=r_med,group=estimate,colour=estimate))+
   theme_bw()+
   geom_errorbar(position=position_dodge(width=0.25))+
   geom_point(aes(x=date,y=r_med))+
   scale_colour_manual(values = c("#d72631","#077b8a"))+
   theme(legend.position="bottom")+
   geom_hline(yintercept = 0,linetype="dashed")+
   labs(x="Date",y="Instantaneous growth rate",col="",tag="A"),
   ggplot(output_summary_comb,
        aes(x=date,ymin=r_lower,ymax=r_upper,group=estimate,fill=estimate))+
   theme_bw()+
   geom_ribbon(alpha=0.55)+
   scale_fill_manual(values = c("#d72631","#077b8a"))+
   theme(legend.position = "bottom")+
   geom_hline(yintercept = 0,linetype="dashed")+
   labs(x="Date",y="Instantaneous growth rate",fill="",tag="B"),
   ggplot(output_summary_comb,
        aes(x=date,ymin=R_lower,ymax=R_upper,y=R_med,group=estimate,colour=estimate))+
   theme_bw()+
   geom_errorbar(position=position_dodge(width=0.25))+
   geom_point(aes(x=date,y=R_med))+
   scale_colour_manual(values = c("#d72631","#077b8a"))+
   theme(legend.position="bottom")+
   geom_hline(yintercept = 1,linetype="dashed")+
   labs(x="Date",y="Effective reproduction number",col="",tag="C"),
    ggplot(output_summary_comb,
        aes(x=date,ymin=R_lower,ymax=R_upper,group=estimate,fill=estimate))+
   theme_bw()+
   geom_ribbon(alpha=0.55)+
   scale_fill_manual(values = c("#d72631","#077b8a"))+
   theme(legend.position = "bottom")+
   geom_hline(yintercept = 1,linetype="dashed")+
   labs(x="Date",y="Effective reproduction number",fill=""),
   align="hv",axis="tl"
)
 ggsave(paste0("outputs/grouped_output/",area,"_",knot_prop*100,".png"),height=7,width=8)
 
 
 ### compare midpoints with the medians for each of the groups
 
 
 
  # output_summary_comb %>% dplyr::select(reporting_group,date,r_med,r_midpoint,R_med,R_midpoint) %>%
  #   pivot_wider(names_from = variable,
  #               values_from = c(r_med,r_midpoint))

 
 output_summary_for_scatter <- merge(output_it_grouped_summary,
                                     output_summary_comb %>% filter(estimate=="Government published") %>% 
                                       dplyr::select(reporting_group,date,r_midpoint,R_midpoint),
       by=c("date","reporting_group"),
       all.x=TRUE)
 
 r_extrema <- max(c(abs(min(output_summary_for_scatter$r.5.,na.rm=TRUE)),
                abs(min(output_summary_for_scatter$r_midpoint,na.rm=TRUE)),
                abs(max(output_summary_for_scatter$r.5.,na.rm=TRUE)),
                abs(max(output_summary_for_scatter$r_midpoint,na.rm=TRUE))))
 
 R_min <- min(min(output_summary_for_scatter$R.5.,na.rm=TRUE),
              min(output_summary_for_scatter$R_midpoint,na.rm=TRUE))
 
 R_max <- max(max(output_summary_for_scatter$R.5.,na.rm=TRUE),
              max(output_summary_for_scatter$R_midpoint,na.rm=TRUE))
 
 
 plot_grid( ggplot(output_summary_for_scatter,
        aes(x=r_midpoint,y=r.50.))+
   theme_bw()+
   geom_point(col="#a2d5c6")+
   geom_smooth(method="lm",fill="#a2d5c6",col="#a2d5c6")+
   geom_abline(slope=1,intercept=0,linetype="dashed")+
   ylim(c(-r_extrema,r_extrema))+
   xlim(c(-r_extrema,r_extrema))+
     labs(x="Government published estimates",y="Modelled estimates",tag="A"),
   ggplot(output_summary_for_scatter,
        aes(x=R_midpoint,y=R.50.))+
   theme_bw()+
   geom_point(col="#a2d5c6")+
   geom_smooth(method="lm",fill="#a2d5c6",col="#a2d5c6")+
   geom_abline(slope=1,intercept=0,linetype="dashed")+
   ylim(c(R_min,R_max))+
   xlim(c(R_min,R_max))+
     labs(x="Government published estimates",y="Modelled estimates",tag="B"))
 ggsave(paste0("outputs/grouped_output/",area,"_",knot_prop*100,"_agreement.png"),height=4,width=8)
 

}



#library(parallel)
cl <- makeCluster(6)
clusterExport(cl = cl, c("get_estimates_sens","ons","spline_prev_model",
                         "gen_time","official"))

#### Scotland
scotland_0.2 <- get_estimates_sens(data = ons,area = "Scotland",
                                  knot_prop = 0.2,
                                  spline_degree = 3,stan_model = spline_prev_model,
                                  iterations = 20000,warm = 2000,
                                  no_chains = 4,
                                  tau_max = 14,
                                  generation_time = gen_time,
                                  official_data = official,
                                  cl = cl)

scotland_0.3 <- get_estimates_sens(data = ons,area = "Scotland",
                                   knot_prop = 0.3,
                                   spline_degree = 3,stan_model = spline_prev_model,
                                   iterations = 20000,warm = 2000,
                                   no_chains = 4,
                                   tau_max = 14,
                                   generation_time = gen_time,
                                   official_data = official,
                                   cl = cl)

scotland_0.4 <- get_estimates_sens(data = ons,area = "Scotland",
                                   knot_prop = 0.4,
                                   spline_degree = 3,stan_model = spline_prev_model,
                                   iterations = 20000,warm = 2000,
                                   no_chains = 4,
                                   tau_max = 14,
                                   generation_time = gen_time,
                                   official_data = official,
                                   cl = cl)






#### England
england_0.2 <- get_estimates_sens(data = ons,area = "England",
                                   knot_prop = 0.2,
                                   spline_degree = 3,stan_model = spline_prev_model,
                                   iterations = 20000,warm = 2000,
                                   no_chains = 4,
                                   tau_max = 14,
                                   generation_time = gen_time,
                                   official_data = official,
                                   cl = cl)

england_0.3 <- get_estimates_sens(data = ons,area = "England",
                                   knot_prop = 0.3,
                                   spline_degree = 3,stan_model = spline_prev_model,
                                   iterations = 20000,warm = 2000,
                                   no_chains = 4,
                                   tau_max = 14,
                                   generation_time = gen_time,
                                   official_data = official,
                                   cl = cl)

england_0.4 <- get_estimates_sens(data = ons,area = "England",
                                   knot_prop = 0.4,
                                   spline_degree = 3,stan_model = spline_prev_model,
                                   iterations = 20000,warm = 2000,
                                   no_chains = 4,
                                   tau_max = 14,
                                   generation_time = gen_time,
                                   official_data = official,
                                   cl = cl)




#### Wales
wales_0.2 <- get_estimates_sens(data = ons,area = "Wales",
                                  knot_prop = 0.2,
                                  spline_degree = 3,stan_model = spline_prev_model,
                                  iterations = 20000,warm = 2000,
                                  no_chains = 4,
                                  tau_max = 14,
                                  generation_time = gen_time,
                                  official_data = official,
                                  cl = cl)

wales_0.3 <- get_estimates_sens(data = ons,area = "Wales",
                                  knot_prop = 0.3,
                                  spline_degree = 3,stan_model = spline_prev_model,
                                  iterations = 20000,warm = 2000,
                                  no_chains = 4,
                                  tau_max = 14,
                                  generation_time = gen_time,
                                  official_data = official,
                                  cl = cl)

wales_0.4 <- get_estimates_sens(data = ons,area = "Wales",
                                  knot_prop = 0.4,
                                  spline_degree = 3,stan_model = spline_prev_model,
                                  iterations = 20000,warm = 2000,
                                  no_chains = 4,
                                  tau_max = 14,
                                  generation_time = gen_time,
                                  official_data = official,
                                  cl = cl)



#### Northern Ireland
NorthernIreland_0.2 <- get_estimates_sens(data = ons,area = "Northern Ireland",
                                knot_prop = 0.2,
                                spline_degree = 3,stan_model = spline_prev_model,
                                iterations = 20000,warm = 2000,
                                no_chains = 4,
                                tau_max = 14,
                                generation_time = gen_time,
                                official_data = official,
                                cl = cl)

NorthernIreland_0.3 <- get_estimates_sens(data = ons,area = "Northern Ireland",
                                knot_prop = 0.3,
                                spline_degree = 3,stan_model = spline_prev_model,
                                iterations = 20000,warm = 2000,
                                no_chains = 4,
                                tau_max = 14,
                                generation_time = gen_time,
                                official_data = official,
                                cl = cl)

NorthernIreland_0.4 <- get_estimates_sens(data = ons,area = "Northern Ireland",
                                knot_prop = 0.4,
                                spline_degree = 3,stan_model = spline_prev_model,
                                iterations = 20000,warm = 2000,
                                no_chains = 4,
                                tau_max = 14,
                                generation_time = gen_time,
                                official_data = official,
                                cl = cl)





stopCluster(cl)

