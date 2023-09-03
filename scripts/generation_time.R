## generation time update

library(readxl)
library(tidyverse)

gt_ests <- read_excel("documents/generation_time.xlsx",sheet="estimates") %>% 
  mutate(date_begin = as.Date(date_begin),
         date_end = as.Date(date_end))
gt_ests$date_end[nrow(gt_ests)] <- Sys.Date()+300  # until finalise analysis date do this


generation_time <- gt_ests %>% complete(date_begin = seq.Date(min(date_begin),max(date_end),by="day")) %>%
  dplyr::select(-date_end) %>%
  rename(date = date_begin) %>% 
  fill(central,sd,distribution,.direction="down") %>%
  mutate(gen_alpha = ifelse(distribution=="gamma",(central^2)/(sd^2),NA),
         gen_beta = ifelse(distribution=="gamma",central/(sd^2),NA))
saveRDS(generation_time,"data/generation_time.RDS")


ggplot(generation_time %>% filter(date>as.Date("2021-01-01")),aes(x=date))+
  geom_line(aes(y=gen_alpha,col="Shape"))+
  geom_line(aes(y=gen_beta,col="Rate"))+
  theme_bw()+
  labs(x="date",y="Value",col="Parameter")+
  scale_x_date(date_labels="%b-%Y",date_breaks = "3 months")
ggsave("outputs/generation_time/gen_time_overview.png",height=3)


# statistically significant differences between 

which(generation_time$date==as.Date("2020-04-01"))

early_2020 <- rgamma(1000,shape=generation_time$gen_alpha[1],rate=generation_time$gen_beta[1])
mid_2020 <- rlnorm(1000,#shape=generation_time$gen_alpha[which(generation_time$date==as.Date("2020-04-01"))],
                   #rate=generation_time$gen_beta[which(generation_time$date==as.Date("2020-04-01"))]
                   mean=generation_time$central[which(generation_time$date==as.Date("2020-04-01"))],
                   sd=generation_time$sd[which(generation_time$date==as.Date("2020-04-01"))])

hist(early_2020)
hist(mid_2020)

t.test(early_2020,mid_2020) # statistically significant differences 

ks.test(early_2020,mid_2020) # statistically significant differences 



### plot the different generation time distributions 

generation_time_params <- gt_ests %>% mutate(gen_alpha = ifelse(distribution=="gamma",(central^2)/(sd^2),NA),
                                             gen_beta = ifelse(distribution=="gamma",central/(sd^2),NA)) %>%
  mutate(time_label = ifelse(date_begin==as.Date("2020-01-01"),"Wildtype (Jan-Mar)",
                             ifelse(date_begin==as.Date("2020-04-01"),"Wildtype (Apr-Nov)",
                                    ifelse(date_begin==as.Date("2020-12-26"),"Alpha",
                                           ifelse(date_begin==as.Date("2021-05-22"),"Delta",
                                                  ifelse(date_begin==as.Date("2021-12-22"),"Omicron",NA))))))

gen_time_samples <- data.frame()

for(i in 1:nrow(generation_time_params)){
  
  if(generation_time_params$distribution[i]=="gamma"){
      samples <- rgamma(1000,
                    shape=generation_time_params$gen_alpha[i],
                    rate=generation_time_params$gen_beta[i])
  df_it <- data.frame("time_label"=generation_time_params$time_label[i],
                      "samples"=samples)
  } else if(generation_time_params$distribution[i]=="lognormal"){
    samples <- rlnorm(1000,
                      meanlog=generation_time_params$central[i],
                      sdlog=generation_time_params$sd[i])
    df_it <- data.frame("time_label"=generation_time_params$time_label[i],
                        "samples"=samples)
  } else{
    
  }
  
  gen_time_samples <- rbind(gen_time_samples,df_it)
}


gen_time_samples <- gen_time_samples %>% mutate(time_label = factor(time_label,
                                                                    levels=c("Wildtype (Jan-Mar)","Wildtype (Apr-Nov)",
                                                                             "Alpha","Delta","Omicron")))

ggplot(gen_time_samples %>% filter(time_label!="Wildtype (Jan-Mar)"),
       aes(x=samples,fill=time_label,col=time_label))+
  geom_density(alpha=0.5)+
  facet_wrap(~time_label,ncol=1)+
  scale_x_continuous(limits = c(0,100))+
  theme_bw()+
  scale_color_viridis_d(end=0.9)+
  scale_fill_viridis_d(end=0.9)+
  theme(strip.background = element_rect(fill="white"),
        legend.position="bottom")+
  labs(x="Generation time",col="Time period",fill="Time period")







