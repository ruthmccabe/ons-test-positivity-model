## plot results for 2023

library(tidyverse)
library(readxl)
library(cowplot)

knots <- seq(20,40,10)
areas <- c("England","Scotland","Wales","Northern Ireland")

ests <- c()

for(a in areas){
  for(k in knots){
    
  ests_it <- readRDS(paste0("outputs/2023/objects/r_R_summary_",k,"_",a,"_2023.RDS")) %>% 
    mutate(knot_prop=k,areaName=a)
  ests <- rbind(ests,ests_it)
    
  }

}

## run script get_spline_fits_2023.R

spline_fits <- readRDS("outputs/2023/objects/spline_fits_combined.RDS") %>%
  filter(knot_prop %in% knots)

ons <- read_excel("data/ons_clean.xlsx",
                  sheet="positivity_long") %>%
  mutate(date_low = as.Date(date_low,"%d %B %Y"),
         date_high = as.Date(date_high,"%d %B %Y"),
         date = date_low+((date_high-date_low)/2),
         avg_pos = avg_pos/100,
         lower_pos = lower_pos/100,
         upper_pos = upper_pos/100) 


full_comb <- merge(spline_fits,ons,all.x=TRUE) %>%
  mutate(areaName = factor(areaName,
                           levels=c("England","Scotland","Wales","Northern Ireland")))

ests <- ests %>%
  mutate(areaName = factor(areaName,
                           levels=c("England","Scotland","Wales","Northern Ireland")))


# ggplot(full_comb %>% filter(date>=as.Date("2023-01-01")&date<max(ons$date)),
#        aes(x=date,y=Y.50.*100,ymin=Y.2.5.*100,ymax=Y.97.5.*100))+
#   geom_ribbon(aes(fill=factor(knot_prop)),alpha=0.2)+
#   geom_line(aes(col=factor(knot_prop)),lwd=1)+
#   geom_point(aes(x=date,y=avg_pos*100))+
#   geom_errorbar(aes(x=date,ymin=lower_pos*100,ymax=upper_pos*100))+
#   theme_bw()+
#   facet_wrap(~areaName,nrow=1)+
#   scale_fill_viridis_d(option="mako",begin=0.2,end=0.8)+
#   scale_colour_viridis_d(option="mako",begin=0.2,end=0.8)+
#   labs(x="Date (2023)",y="Estimated test positivity",
#        col="Knots (as % of data points)",fill="Knots (as % of data points)")+
#   theme(legend.position = "bottom",
#         strip.background = element_rect(fill="white"))+
#   scale_x_date(date_labels = "%d-%b")+
#   guides(col = guide_legend(title.position = "top",nrow=1),
#          fill = guide_legend(title.position = "top",nrow=1))




plot_grid(ggplot(full_comb %>% filter(date>=as.Date("2023-01-01")&date<max(ons$date)),
                 aes(x=date,y=Y.50.*100,ymin=Y.2.5.*100,ymax=Y.97.5.*100))+
            geom_ribbon(aes(fill=factor(knot_prop)),alpha=0.2)+
            geom_line(aes(col=factor(knot_prop)),lwd=1)+
            geom_point(aes(x=date,y=avg_pos*100))+
            geom_errorbar(aes(x=date,ymin=lower_pos*100,ymax=upper_pos*100))+
            theme_bw()+
            facet_wrap(~areaName,nrow=1)+
            # scale_fill_viridis_d(option="mako",begin=0.2,end=0.8)+
            # scale_colour_viridis_d(option="mako",begin=0.2,end=0.8)+
            scale_colour_manual(values=c("#588157","#62B6CB","#1B4965"))+
            scale_fill_manual(values=c("#588157","#62B6CB","#1B4965"))+
            labs(x="Date (2023)",y="Estimated test \npositivity",
                 col="Knots (as % of data points)",fill="Knots (as % of data points)",tag="A")+
            theme(legend.position = "none",
                  strip.background = element_rect(fill="white"))+
            scale_x_date(date_labels = "%d-%b",date_breaks="1 month")+
            guides(col = guide_legend(title.position = "top",nrow=1),
                   fill = guide_legend(title.position = "top",nrow=1)),  
          # rep num
          ggplot(ests %>% filter(date>=as.Date("2023-01-01")&date<max(ons$date)),
       aes(x=date,ymin=R.2.5.,y=R.50.,ymax=R.97.5.))+
  geom_ribbon(aes(fill=factor(knot_prop)),alpha=0.2)+
  geom_line(aes(col=factor(knot_prop)),lwd=1)+
  theme_bw()+
  geom_hline(yintercept=1,linetype="dashed")+
  facet_wrap(~areaName,nrow=1)+
  # scale_fill_viridis_d(option="mako",begin=0.2,end=0.8)+
  # scale_colour_viridis_d(option="mako",begin=0.2,end=0.8)+
    scale_colour_manual(values=c("#588157","#62B6CB","#1B4965"))+
    scale_fill_manual(values=c("#588157","#62B6CB","#1B4965"))+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "none",
        legend.title.align = 0.5)+
  guides(col = guide_legend(title.position = "left",nrow=1),
         fill = guide_legend(title.position = "left",nrow=1))+
  labs(x="Date (2023)",y="Effective reproduction \nnumber R(t)",
       col="Knots (as % of data points)",fill="Knots (as % of data points)",tag="B")+
  scale_x_date(date_labels = "%d-%b",date_breaks="1 month"),
  # growth rate
          ggplot(ests %>% filter(date>=as.Date("2023-01-01")&date<max(ons$date)),
       aes(x=date,ymin=r.2.5.,y=r.50.,ymax=r.97.5.))+
  geom_ribbon(aes(fill=factor(knot_prop)),alpha=0.2)+
  geom_line(aes(col=factor(knot_prop)),lwd=1)+
  theme_bw()+
  geom_hline(yintercept=0,linetype="dashed")+
  facet_wrap(~areaName,nrow=1)+
  # scale_fill_viridis_d(option="mako",begin=0.2,end=0.8)+
  # scale_colour_viridis_d(option="mako",begin=0.2,end=0.8)+
    scale_colour_manual(values=c("#588157","#62B6CB","#1B4965"))+
    scale_fill_manual(values=c("#588157","#62B6CB","#1B4965"))+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "bottom",
        legend.title.align = 0.5)+
    guides(col = guide_legend(title.position = "left",nrow=1),
           fill = guide_legend(title.position = "left",nrow=1))+
  labs(x="Date (2023)",y="Instantaneous \ngrowth rate r(t)",
       col="Knots (as % of data points)",fill="Knots (as % of data points)",tag="C")+
  scale_x_date(date_labels = "%d-%b",date_breaks="1 month"),
  nrow=3,align="v",axis="l",rel_heights = c(1,1,1.3)
  )
ggsave("outputs/2023/overview.png",width=8,height=7)





  

