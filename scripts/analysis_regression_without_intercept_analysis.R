### regression without intercept analysis

library(tidyverse)
library(readxl)
library(matrixStats)
library(rstan)
library(Hmisc)
library(cowplot)

knots <- seq(20,40,10)
areas <- c("England","Scotland","Wales","Northern Ireland")

lm_results <- c()

for(k in knots){
  for(a in areas){
    
      lm_results_it <- readRDS(paste0("outputs/objects/lm_results_",
                                  k,"_",a,".RDS")) %>%
        mutate(areaName=a)
  
  lm_results <- rbind(lm_results,lm_results_it)
  
  }

}

lm_results <- lm_results %>% 
  mutate(metric = case_when(metric=="Instantaeous growth rate" ~ "Instantaneous growth rate",
                            metric=="Instanaeous growth rate" ~ "Instantaneous growth rate",
                            metric=="Effective reproduction number" ~ "Effective reproduction number"
                            ),
         areaName = factor(areaName,levels=c("England","Scotland","Wales","Northern Ireland"))) 

# ggplot(lm_results %>% filter(intercept=="No Intercept"),
#        aes(x=lag,y=coefficient,col=as.factor(knot_prop)))+
#   theme_bw()+
#   geom_line()+
#   facet_grid(metric~areaName)+
#   geom_hline(yintercept=1,linetype="dashed")+
#   theme(strip.background = element_rect(fill="white"),
#         legend.position="bottom")+
#   labs(x="Lag (days)",y="Linear regression coefficient beta1",col="Knots (as % of data points)")+
#   scale_colour_manual(values=c("#588157","#62B6CB","#1B4965"))+
#   scale_y_continuous(breaks=seq(-0.2,1,0.2))
# ggsave("outputs/analysis/regression_without_intercept")


lm_results_long <- lm_results %>% pivot_longer(cols=c("R2","coefficient"),
                            names_to = "lm_output",
                             values_to = "value") #%>%
  # mutate(lm_output = case_when(lm_output=="coefficient" ~ "Coefficient (Beta1)",
  #                              lm_output=="R2" ~ "R2"))




# ggplot(lm_results_long %>% filter(intercept=="No Intercept",lm_output=="Coefficient (Beta1)"),
#        aes(x=lag,y=value,col=as.factor(knot_prop)))+
#   theme_bw()+
#   geom_line(lwd=1)+
#   facet_grid(areaName~metric)+
#   geom_hline(yintercept=1,linetype="dashed")+
#   theme(strip.background = element_rect(fill="white"),
#         legend.position="bottom")+
#   labs(x="Lag (days)",y="Value of output",col="Knots (as % of data points)")+
#   scale_colour_manual(values=c("#588157","#62B6CB","#1B4965"))+
#   scale_y_continuous(breaks=seq(-0.2,1,0.2))
# ggsave("outputs/analysis/regression_without_intercept.png",width=8,height=6)




lm_results_ci <- lm_results %>%
  mutate(lower_ci = coefficient-(1.96*se),
         upper_ci = coefficient+(1.96*se)) %>%
  filter(lag<=20&lag>=-20)


plot_grid(ggplot(lm_results_ci %>% filter(intercept=="No Intercept",
                                          metric=="Effective reproduction number"))+
            theme_bw()+
            geom_ribbon(aes(x=lag,ymin=lower_ci,ymax=upper_ci,fill=as.factor(knot_prop)),alpha=0.15)+
            geom_line(aes(x=lag,y=coefficient,col=as.factor(knot_prop)),lwd=1)+
            facet_grid(~areaName)+
            geom_hline(yintercept=1,linetype="dashed")+
            theme(strip.background = element_rect(fill="white"),
                  legend.position="none")+
            labs(x="Lag (days)",y="Coefficient (Beta1) \nfor R(t) models",col="Knots (as % of data points)",
                 tag="A")+
            scale_colour_manual(values=c("#588157","#62B6CB","#1B4965"))+
            scale_fill_manual(values=c("#588157","#62B6CB","#1B4965"))+
            scale_y_continuous(breaks=seq(0.9,1.1,0.05)),
          ggplot(lm_results_ci %>% filter(intercept=="No Intercept",
                                          metric=="Instantaneous growth rate"))+
            theme_bw()+
            geom_ribbon(aes(x=lag,ymin=lower_ci,ymax=upper_ci,fill=as.factor(knot_prop)),alpha=0.15)+
            geom_line(aes(x=lag,y=coefficient,col=as.factor(knot_prop)),lwd=1)+
            facet_grid(~areaName)+
            geom_hline(yintercept=1,linetype="dashed")+
            theme(strip.background = element_rect(fill="white"),
                  legend.position="bottom")+
            labs(x="Lag (days)",y="Coefficient (Gamma1) \nfor r(t) models",col="Knots (as % of data points)",
                 fill="Knots (as % of data points)",tag="B")+
            scale_colour_manual(values=c("#588157","#62B6CB","#1B4965"))+
            scale_fill_manual(values=c("#588157","#62B6CB","#1B4965"))+
            scale_y_continuous(breaks=seq(-0.2,1,0.4)),
          nrow=2,rel_heights=c(1,1.1))
ggsave("outputs/analysis/regression_without_intercept.png",width=8,height=6)



# plot_grid(ggplot(lm_results_long %>% 
#                    filter(intercept=="No Intercept",metric=="Effective reproduction number",
#                           lm_output=="Coefficient (Beta1)"),
#                  aes(x=lag,y=value,col=as.factor(knot_prop*100)))+
#             theme_bw()+
#             geom_line(lwd=1)+
#             facet_grid(~areaName)+
#             geom_hline(yintercept=1,linetype="dashed")+
#             theme(strip.background = element_rect(fill="white"),
#                   legend.position="none")+
#             labs(x="Lag (days)",y="Coefficient (Beta1) \nfor R(t) models",col="Knots (as % of data points)",
#                  tag="A")+
#             scale_colour_manual(values=c("#588157","#62B6CB","#1B4965"))+
#             scale_y_continuous(breaks=seq(0.9,1.1,0.05)),
#           ggplot(lm_results_long %>% 
#                    filter(intercept=="No Intercept",metric=="Instantaneous growth rate",
#                           lm_output=="Coefficient (Beta1)"),
#                  aes(x=lag,y=value,col=as.factor(knot_prop*100)))+
#             theme_bw()+
#             geom_line(lwd=1)+
#             facet_grid(~areaName)+
#             geom_hline(yintercept=1,linetype="dashed")+
#             theme(strip.background = element_rect(fill="white"),
#                   legend.position="bottom")+
#             labs(x="Lag (days)",y="Coefficient (Gamma1) \nfor r(t) models",col="Knots (as % of data points)",
#                  tag="B")+
#             scale_colour_manual(values=c("#588157","#62B6CB","#1B4965"))+
#             scale_y_continuous(breaks=seq(-0.2,1,0.5)),
#           nrow=2,rel_heights=c(1,1.1))











