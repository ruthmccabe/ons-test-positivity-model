### England

library(tidyverse)
library(readxl)
library(matrixStats)
library(rstan)
library(Hmisc)
library(cowplot)

knots <- seq(20,40,10)

correlations <- c()
prop_agreement <- c()
lm_results <- c()

for(k in knots){
  
  correlations_it <- readRDS(paste0("outputs/objects/correlation_",
                                    k,
                                    "_England.RDS"))
  
  correlations <- rbind(correlations,correlations_it)
  
  prop_agreement_it <- readRDS(paste0("outputs/objects/prop_agreement_",
                                      k,
                                      "_England.RDS"))
  
  prop_agreement <- rbind(prop_agreement,prop_agreement_it)
  
  lm_results_it <- readRDS(paste0("outputs/objects/lm_results_",
                                  k,
                                  "_England.RDS"))
  
  lm_results <- rbind(lm_results,lm_results_it)
  
  
}


assessment <- rbind(correlations %>% mutate(assessment="Spearman correlation") %>%
                      rename(value = correlation),
                    prop_agreement %>% mutate(assessment="Agreement proportion") %>% 
                      select(prop_agreement,metric,knot_prop,areaName,lag,assessment) %>%
                      rename(value = prop_agreement),
                    lm_results %>% filter(intercept=="Intercept") %>% mutate(assessment = "R2") %>%
                      select(R2,metric,knot_prop,areaName,lag,assessment) %>%
                      rename(value = R2)) %>%
  mutate(metric = ifelse(metric=="r","Instantaneous growth rate",
                         ifelse(metric=="R","Effective reproduction number",
                                ifelse(metric=="Instanaeous growth rate","Instantaneous growth rate",metric))))

assessment <- assessment %>% group_by(metric,assessment,knot_prop) %>%
  mutate(max_value = max(value),
         opt_lag = ifelse(value==max_value,lag,NA))

assessment %>% group_by(assessment, metric) %>%
  filter(value==max(value))


# plot_grid(ggplot(assessment %>% filter(lag<=10,lag>-20,
#                                        assessment=="Adjusted R2"),
#                  aes(x=knot_prop*100,y=lag,fill=value))+
#             theme_bw()+
#             geom_tile()+
#             facet_grid(metric~assessment)+
#             theme(strip.background = element_rect(fill="white"),
#                   legend.position = "bottom")+
#             labs(x="",
#                  y="Lag applied to government-published estimates (days)",fill="Adjusted R2")+
#             guides(fill=guide_colorbar(ticks=FALSE,title.position = "top"))+
#             scale_fill_viridis_c(end=0.8,breaks=seq(0,1,0.2)),
#           ggplot(assessment %>% filter(lag<=10,lag>-20,
#                                        assessment=="Spearman correlation"),
#                  aes(x=knot_prop*100,y=lag,fill=value))+
#             theme_bw()+
#             geom_tile()+
#             facet_grid(metric~assessment)+
#             theme(strip.background = element_rect(fill="white"),
#                   legend.position="bottom")+
#             labs(x="Knots (% of data points)",
#                  y="",fill="Spearman Correlation")+
#             guides(fill=guide_colorbar(ticks=FALSE,title.position = "top"))+
#             scale_fill_viridis_c(option="inferno",end=0.9,breaks=seq(0,1,0.2)),
#           ggplot(assessment %>% filter(lag<=10,lag>-20,
#                                        assessment=="Agreement proportion"),
#                  aes(x=knot_prop*100,y=lag,fill=value))+
#             theme_bw()+
#             geom_tile()+
#             facet_grid(metric~assessment)+
#             theme(strip.background = element_rect(fill="white"),
#                   legend.position="bottom")+
#             labs(x="",
#                  y="",fill="Agreement proportion")+
#             guides(fill=guide_colorbar(ticks=FALSE,title.position = "top"))+
#             scale_fill_viridis_c(option="rocket",breaks=seq(0,1,0.2)),
#           nrow=1
# )
# ggsave("outputs/analysis/England_heatmap.png",width=7,height=5.5)



# ggplot(assessment %>% filter(lag<=10,lag>-20),aes(x=lag,y=value,col=factor(knot_prop*100)))+
#   theme_bw()+
#   geom_line(lwd=1)+
#   facet_grid(metric~assessment)+
#   theme(strip.background = element_rect(fill="white"),
#         legend.position="bottom")+
#   #geom_vline(aes(xintercept=opt_lag,col=factor(knot_prop)),linetype="dashed")+
#   labs(x="Lag applied to government-published estimates (days)",
#        y="Value of assessment metric",
#        col="Knots (% of data points)")+
#   scale_y_continuous(breaks=seq(0,1,0.2))+
#   scale_colour_viridis_d(option="mako",begin=0.2,end=0.8)
# 
# ggsave("outputs/analysis/England_metrics.png",width=6,height=5)


### maximal lags under the adjusted R2 for 40 knots

assessment %>% filter(knot_prop==0.3,
                      assessment=="R2") %>% group_by(metric) %>%
  filter(value==max(value))


## get all the spline fits as with 2023

spline_fit_daily <- c()

for(k in knots){
  
  spline_fit_daily_it <- readRDS(paste0("outputs/objects/spline_fit_daily_",
                                        k,
                                        "_England.RDS")) %>%
    mutate(knot_prop = k)
  
  spline_fit_daily <- rbind(spline_fit_daily,spline_fit_daily_it)
  
}

ons <- read_excel("data/ons_clean.xlsx",
                  sheet="positivity_long") %>%
  mutate(date_low = as.Date(date_low,"%d %B %Y"),
         date_high = as.Date(date_high,"%d %B %Y"),
         date = date_low+((date_high-date_low)/2),
         avg_pos = avg_pos/100,
         lower_pos = lower_pos/100,
         upper_pos = upper_pos/100) 

ons_area <- ons %>% filter(areaName=="England",
                           date<=as.Date("2022-12-28"))


panel_A <- ggplot()+
  theme_bw()+
  geom_ribbon(data=spline_fit_daily %>% filter(knot_prop==30),
              aes(x=date,ymin=Y.2.5.*100,ymax=Y.97.5.*100),alpha=0.6,fill="#077b8a")+
  geom_line(data=spline_fit_daily %>% filter(knot_prop==30),
            aes(x=date,y=Y.50.*100),lwd=1.5,col="#077b8a")+
  geom_point(data=ons_area,aes(x=date,y=avg_pos*100))+
  geom_errorbar(data=ons_area,aes(x=date,ymin=lower_pos*100,ymax=upper_pos*100))+
  labs(x="Date",y="ONS CIS Test \nPositivity (%)",tag="A")+
  scale_x_date(date_breaks = "3 months",date_labels = "%b-%y")+
  scale_y_continuous(breaks=seq(0,12,2),limits = c(0,9.5))



## get all of the R and r ests

estimates <- c()

for(k in knots){
  
  estimates_it <- readRDS(paste0("outputs/objects/r_R_summary_",
                                 k,
                                 "_England.RDS")) %>%
    mutate(knot_prop=k)
  
  estimates <- rbind(estimates,estimates_it)
  
}


official <- read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&metric=transmissionRateGrowthRateMax&metric=transmissionRateMax&metric=transmissionRateMin&metric=transmissionRateGrowthRateMin&format=csv") %>%
  mutate(date = as.Date(date)) %>% group_by(areaName) %>% 
  complete(date = seq.Date(min(date),max(date),by="day")) %>%
  fill(areaCode,areaName,areaType,
       transmissionRateGrowthRateMax,transmissionRateGrowthRateMin,
       transmissionRateMax,transmissionRateMin,.direction = "up") %>% 
  #mutate(date = date-14) %>% 
  data.frame() %>%
  mutate(transmissionRateGrowthRateMid = (transmissionRateGrowthRateMax+transmissionRateGrowthRateMin)/2,
         transmissionRateMid = (transmissionRateMax+transmissionRateMin)/2)


## merge with official 

official_comb <- merge(estimates %>% filter(knot_prop==30),
                       official %>% filter(areaName=="England")) %>% 
  mutate(transmissionRateGrowthRateMin.lag = Lag(transmissionRateGrowthRateMin,-9),
         transmissionRateGrowthRateMax.lag = Lag(transmissionRateGrowthRateMax,-9),
         transmissionRateGrowthRateMid.lag = Lag(transmissionRateGrowthRateMid,-9),
         transmissionRateMin.lag = Lag(transmissionRateMin,-8),
         transmissionRateMax.lag = Lag(transmissionRateMax,-8),
         transmissionRateMid.lag = Lag(transmissionRateMid,-8)) 

panel_B <- ggplot()+
  theme_bw()+
  geom_ribbon(data=official_comb,
              aes(x=date,ymin=transmissionRateMin.lag,ymax=transmissionRateMax.lag),alpha=0.45,fill="#d72631")+
  geom_ribbon(data=official_comb,aes(x=date,ymin=R.2.5.,ymax=R.97.5.),alpha=0.6,fill="#077b8a")+
  geom_line(data=official_comb,aes(x=date,y=R.50.),lwd=1.5,col="#077b8a")+
  geom_abline(slope=0,intercept=1,linetype="dashed")+
  labs(x="Date",y="Effective reproduction \nnumber R(t)",tag="B")+
  scale_x_date(date_breaks = "5 months",date_labels = "%b-%y")

panel_C <- ggplot(official_comb)+
  theme_bw()+
  geom_point(aes(y=transmissionRateMid.lag,x=R.50.),col="#a2d5c6",pch=1)+
  geom_smooth(aes(y=transmissionRateMid.lag,x=R.50.),
              col="#08605F",fill="#08605F",method="lm",lwd=0.5)+
  geom_abline(slope=1,intercept=0,linetype="dashed")+
  labs(y="Government-published \nestimates",x="ONS-based estimates",tag="C")+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  xlim(c(0.7,1.4))+
  ylim(c(0.7,1.4))


panel_D <- ggplot()+
  theme_bw()+
  geom_ribbon(data=official_comb,
              aes(x=date,ymin=transmissionRateGrowthRateMin.lag/100,
                  ymax=transmissionRateGrowthRateMax.lag/100),alpha=0.45,
              #fill="#b1f3fb"
              fill="#d72631")+
  geom_ribbon(data=official_comb,aes(x=date,ymin=r.2.5.,ymax=r.97.5.),alpha=0.6,fill="#077b8a")+
  geom_line(data=official_comb,aes(x=date,y=r.50.),lwd=1.5,col="#077b8a")+
  geom_abline(slope=0,intercept=0,linetype="dashed")+
  labs(x="Date",y="Instantaneous growth \nrate r(t)",tag="D")+
  scale_x_date(date_breaks = "5 months",date_labels = "%b-%y")

panel_E <- ggplot(official_comb)+
  theme_bw()+
  geom_point(aes(y=transmissionRateGrowthRateMid.lag/100,x=r.50.),col="#a2d5c6",pch=1)+
  geom_smooth(aes(y=transmissionRateGrowthRateMid.lag/100,x=r.50.),
              col="#08605F",fill="#08605F",method="lm",lwd=0.5)+
  geom_abline(slope=1,intercept=0,linetype="dashed")+
  labs(y="Government-published \nestimates",x="ONS-based estimates",tag="E")+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  xlim(c(-0.06,0.07))+
  ylim(c(-0.06,0.07))


plot_grid(panel_A,
          plot_grid(panel_B,panel_C,panel_D,panel_E,align="hv",rel_widths = c(1.5,1)),
          ncol=1,rel_heights=c(1,2))


# panel_F <- ggplot(assessment %>% filter(lag<=10,lag>-20),aes(x=lag,y=value,col=factor(knot_prop*100)))+
#   theme_bw()+
#   geom_line(lwd=1)+
#   facet_grid(assessment~metric)+
#   theme(strip.background = element_rect(fill="white"),
#         legend.position="bottom")+
#   #geom_vline(aes(xintercept=opt_lag,col=factor(knot_prop)),linetype="dashed")+
#   labs(x="Lag applied to government-published \nestimates (days)",
#        y="Value of assessment metric",
#        col="Knots (% of data points)",tag="F")+
#   scale_y_continuous(breaks=seq(0,1,0.2))+
#   guides(col = guide_legend(title.position = "top",nrow=1))+
#   scale_colour_manual(values=c("#45425A","#9EB25D","#077b8a"))

# plot_grid(plot_grid(panel_A,
#           plot_grid(panel_B,panel_C,panel_D,panel_E,align="hv",rel_widths = c(1.5,1)),
#           ncol=1,rel_heights=c(1,2)),
#           panel_F,ncol=2,rel_widths = c(2,1),align="hv")
# ggsave("outputs/analysis/England.png",height=8,width=10)



assessment <- assessment %>% 
  mutate(metric_short = factor(ifelse(metric=="Effective reproduction number","R(t)",
                                      ifelse(metric=="Instantaneous growth rate","r(t)",NA)),
                               levels=c("R(t)","r(t)")),
         opt_lag = ifelse(assessment=="R2"&knot_prop==0.3,opt_lag,NA))

panel_F <- ggplot(assessment %>% filter(lag<=20,lag>-20,assessment=="R2"),
                  aes(x=lag,y=value,col=factor(knot_prop*100)))+
  theme_bw()+
  geom_line(lwd=1)+
  facet_grid(~metric_short)+
  theme(strip.background = element_rect(fill="white"),
        legend.position="none")+
  geom_vline(aes(xintercept=opt_lag,col=factor(knot_prop*100)),linetype="dotted",lwd=1)+
  labs(x="Lag (days)",
       y="R2",
       col="Knots (% of data points)",tag="F")+
  scale_y_continuous(breaks=seq(0,1,0.25),limits = c(0,1))+
  guides(col = guide_legend(title.position = "top",nrow=1))+
  scale_colour_manual(values=c("#588157","#62B6CB","#1B4965"))



panel_G <- ggplot(assessment %>% filter(lag<=20,lag>-20,assessment=="Spearman correlation"),
                  aes(x=lag,y=value,col=factor(knot_prop*100)))+
  theme_bw()+
  geom_line(lwd=1)+
  facet_grid(~metric_short)+
  theme(strip.background = element_rect(fill="white"),
        legend.position="none")+
  #geom_vline(aes(xintercept=opt_lag,col=factor(knot_prop)),linetype="dashed")+
  labs(x="Lag (days)",
       y="Spearman correlation",
       col="Knots (% of data points)",tag="G")+
  scale_y_continuous(breaks=seq(0,1,0.25),limits = c(0,1))+
  guides(col = guide_legend(title.position = "top",nrow=1))+
  scale_colour_manual(values=c("#588157","#62B6CB","#1B4965"))


panel_H <- ggplot(assessment %>% filter(lag<=20,lag>-20,assessment=="Agreement proportion"),
                  aes(x=lag,y=value,col=factor(knot_prop*100)))+
  theme_bw()+
  geom_line(lwd=1)+
  facet_grid(~metric_short)+
  theme(strip.background = element_rect(fill="white"),
        legend.position="bottom")+
  #geom_vline(aes(xintercept=opt_lag,col=factor(knot_prop)),linetype="dashed")+
  labs(x="Lag (days)",
       y="Agreement proportion",
       col="Knots (% of data points)",tag="H")+
  scale_y_continuous(breaks=seq(0,1,0.25),limits = c(0,1))+
  guides(col = guide_legend(title.position = "top",nrow=1))+
  scale_colour_manual(values=c("#588157","#62B6CB","#1B4965"))





# plot_grid(panel_A,
#           plot_grid(panel_B,panel_C,panel_D,panel_E,align="hv",rel_widths = c(1.5,1)),
#           ncol=1,rel_heights=c(1,2))
# 
# 
# plot_grid(panel_A,panel_F,nrow=1,rel_widths=c(2,1))
# plot_grid(panel_B,panel_C,panel_G,nrow=1,rel_widths = c(1.5,1,1),align="hv",axis="tb")
# plot_grid(panel_D,panel_E,panel_H,nrow=1,rel_widths = c(1.5,1,1),align="hv",axis="tb")


plot_grid(plot_grid(panel_A,panel_F,nrow=1,rel_widths=c(2.5,1)),
          plot_grid(plot_grid(panel_B,panel_C,panel_G,nrow=1,rel_widths = c(1.5,1,1),align="hv",axis="tb"),
                    plot_grid(panel_D,panel_E,panel_H,nrow=1,rel_widths = c(1.5,1,1),align="hv",axis="tb"),
                    nrow=2,rel_heights=c(1,1.4),align="hv",axis="lrtb"),
          nrow=3,rel_heights=c(1,2))
ggsave("outputs/analysis/England.png",width=10,height=9)


plot_grid(plot_grid(panel_A + scale_x_date(date_breaks = "4 months",date_labels = "%b-%y"),
                    panel_F+labs(tag="D"),nrow=1,align="hv",axis="ltrb",
                    rel_widths=c(1.75,1)),
          plot_grid(panel_B+scale_x_date(date_breaks = "4 months",date_labels = "%b-%y",
                                         limits = c(as.Date("2020-05-03"),as.Date("2022-12-31"))),
                    panel_G+labs(tag="E"),nrow=1,align="hv",axis="ltrb",
                    rel_widths=c(1.75,1)),
          plot_grid(panel_D+labs(tag="C",y="Instantaneous growth rate r(t)")+
                      scale_x_date(date_breaks = "4 months",date_labels = "%b-%y",
                                   limits = c(as.Date("2020-05-03"),as.Date("2022-12-31"))),
                    panel_H+labs(tag="F"),nrow=1,align="hv",axis="ltrb",
                    rel_widths=c(1.75,1)),
          nrow=3,
          rel_heights=c(1,1,1.4))
ggsave("outputs/analysis/England_without_scatter.png",width=9,height=7)


#### plot without any lags for the supplement 

panel_B_nolag <- ggplot()+
  theme_bw()+
  geom_ribbon(data=official_comb,
              aes(x=date,ymin=transmissionRateMin,ymax=transmissionRateMax),alpha=0.45,fill="#d72631")+
  geom_ribbon(data=official_comb,aes(x=date,ymin=R.2.5.,ymax=R.97.5.),alpha=0.6,fill="#077b8a")+
  geom_line(data=official_comb,aes(x=date,y=R.50.),lwd=1.5,col="#077b8a")+
  geom_abline(slope=0,intercept=1,linetype="dashed")+
  labs(x="Date",y="Effective reproduction \nnumber R(t)",tag="B")+
  scale_x_date(date_breaks = "5 months",date_labels = "%b-%y")

panel_C_nolag <- ggplot(official_comb)+
  theme_bw()+
  geom_point(aes(y=transmissionRateMid,x=R.50.),col="#a2d5c6",pch=1)+
  geom_smooth(aes(y=transmissionRateMid,x=R.50.),col="#08605F",fill="#08605F",method="lm",lwd=0.5)+
  geom_abline(slope=1,intercept=0,linetype="dashed")+
  labs(y="Government-published \nestimates",x="ONS-based estimates",tag="C")+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  xlim(c(0.7,1.4))+
  ylim(c(0.7,1.4))


panel_D_nolag <- ggplot()+
  theme_bw()+
  geom_ribbon(data=official_comb,
              aes(x=date,ymin=transmissionRateGrowthRateMin/100,
                  ymax=transmissionRateGrowthRateMax/100),alpha=0.45,
              #fill="#b1f3fb"
              fill="#d72631")+
  geom_ribbon(data=official_comb,aes(x=date,ymin=r.2.5.,ymax=r.97.5.),alpha=0.6,fill="#077b8a")+
  geom_line(data=official_comb,aes(x=date,y=r.50.),lwd=1.5,col="#077b8a")+
  geom_abline(slope=0,intercept=0,linetype="dashed")+
  labs(x="Date",y="Instantaneous growth \nrate r(t)",tag="D")+
  scale_x_date(date_breaks = "5 months",date_labels = "%b-%y")

panel_E_nolag <- ggplot(official_comb)+
  theme_bw()+
  geom_point(aes(y=transmissionRateGrowthRateMid/100,x=r.50.),col="#a2d5c6",pch=1)+
  geom_smooth(aes(y=transmissionRateGrowthRateMid/100,x=r.50.),col="#08605F",fill="#08605F",method="lm",lwd=0.5)+
  geom_abline(slope=1,intercept=0,linetype="dashed")+
  labs(y="Government-published \nestimates",x="ONS-based estimates",tag="E")+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  xlim(c(-0.06,0.07))+
  ylim(c(-0.06,0.07))


plot_grid(panel_A,
          plot_grid(panel_B_nolag,panel_C_nolag,panel_D_nolag,panel_E_nolag,align="hv",rel_widths = c(1.5,1)),
          ncol=1,rel_heights=c(1,2))
ggsave("outputs/analysis/England_nolag.png",height=8,width=10)



### traceplot
fit <- readRDS("outputs/objects/fit_30_England.RDS")
traceplot(fit,pars=c("a[1]","a[12]","a[24]","a[36]","a[48]",
                     "Y_hat[1]","Y_hat[12]","Y_hat[24]","Y_hat[36]","Y_hat[48]",
                     "tau","gamma"))
ggsave("outputs/analysis/England_traceplots.png")




