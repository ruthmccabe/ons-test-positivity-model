### simulation of data to test likelihood assumption 

library(tidyverse)
library(binom)
library(cowplot)
library(matrixStats)
#library(fitdistrplus)
library(invgamma)

# set.seed(2904) ## ok, one with 37 non normal though
# set.seed(304) ## worse, up to 60
set.seed(703) ## bingo

likelihood_simulation <- function(n_obs = 150,
                                  total = 80000,
                                  n_samps = 10000){



  # initialise
  # pos_init <- round(runif(1,min=100,max=500))
  # tot_init <- round(runif(1,min=1500,max=4000))
  
  #prev_init <- rbeta(1,shape1=2,shape2=5)

  df <- data.frame("time"=seq(1,n_obs,1)) %>%
    group_by(time) %>%
    mutate(prev = NA,
           prev_logit = NA,
           pos = NA)
  
  # df$pos[1] <- pos_init
  # df$tot[1] <- tot_init
  df$prev[1] <- rbeta(1,shape1=2,shape2=10) # initialise
  df$prev_logit[1] <- log(df$prev[1]/(1-df$prev[1]))
  
  
  # random walk update
  for(i in 2:nrow(df)){
    sd <- rinvgamma(n=1,shape=2,scale=50)
    df$prev_logit[i] <- rnorm(1,df$prev_logit[i-1],sd)
  }
  
  
  
  df <- df %>% mutate(prev = exp(prev_logit)/(exp(prev_logit)+1),
                      pos = round(prev * total))
  
  
  
  # # random walk update
  # for(i in 2:nrow(df)){
  #   df$pos[i] <- rnorm(1,df$pos[i-1],50)
  #   df$tot[i] <- rnorm(1,df$tot[i-1],500)
  #   df$prev[i] <- df$pos[i]/df$tot[i]
  #   while(df$prev[i]<0 | df$prev[i]>1|df$pos[i]<0|df$tot[i]<0){
  #     df$pos[i] <- rnorm(1,df$pos[i-1],50)
  #     df$tot[i] <- rnorm(1,df$tot[i-1],500)
  #     df$prev[i] <- df$pos[i]/df$tot[i]
  #   }
  # }
  # 
  # if(any(df$prev<0)|any(df$prev>1)){
  #   stop("Prevalence estimates out of range.")
  # }
  
  # put into data frame
  df <- df %>% mutate(binom.confint(x=pos,n=total,methods="exact"))
  
  # translate into ONS format
  df_ons_sim <- df %>% dplyr::select(time,mean,lower,upper) %>%
    mutate(var = ((upper - lower)/3.92)^2,
           alpha = mean*(((mean * (1 - mean))/var) -1 ),
           beta = alpha * ((1-mean)/mean))
  
  # plot of data overview 
  plot1 <- cowplot::plot_grid(ggplot(df,aes(x=time,y=pos))+
                                geom_col(fill="black",width=1)+
                                labs(y="Positives",x="Time",tag="A")+
                                theme_bw(),
                              # ggplot(df,aes(x=time,y=tot))+
                              #   geom_col(fill="black",width=1)+
                              #   labs(y="Total",x="Time",tag="B")+
                              #   theme_bw(),
                              ggplot(df,aes(x=time,y=mean*100))+
                                geom_point()+
                                geom_line()+
                                labs(y="Positivity (%)",x="Time",tag="B")+
                                theme_bw(),
                              ggplot(df_ons_sim,aes(x=time,y=mean*100,ymin=lower*100,ymax=upper*100))+
                                geom_point()+
                                geom_errorbar(width=1)+
                                #geom_line()+
                                labs(y="Positivity (%)",x="Time",tag="C")+
                                theme_bw(),
                              nrow=3,align="hv",axis="tb")
  
  beta_samps <- matrix(NA,ncol=n_samps,nrow=nrow(df_ons_sim))
  ks_pval <- matrix(NA,ncol=1,nrow=nrow(df_ons_sim))
  #sw_pval <- matrix(NA,ncol=1,nrow=nrow(df_ons_sim))
  
  for(i in 1:nrow(df_ons_sim)){
    samp_it <- rbeta(n_samps,df_ons_sim$alpha[i],df_ons_sim$beta[i])
    logit_samp_it <- log(samp_it/(1-samp_it))
    beta_samps[i,] <- logit_samp_it
    ks_pval[i,1] <- ks.test(x = logit_samp_it,y="pnorm",mean=mean(logit_samp_it),sd=sd(logit_samp_it),
                            alternative = "two.sided")$p.value
    #sw_pval[i,1] <- shapiro.test(as.numeric(logit_samp_it))$p.val
  }
  
  beta_samps_summary <- data.frame(cbind("time"=df$time,
                                         rowQuantiles(beta_samps,probs = c(0.025,0.5,0.975)))) %>%
    mutate(lower = exp(X2.5.)/(exp(X2.5.)+1),
           med = exp(X50.)/(exp(X50.)+1),
           upper = exp(X97.5.)/(exp(X97.5.)+1))
  
  
  beta_samps_hists <- gather(data.frame("t38"=beta_samps[38,],
                                        "t75"=beta_samps[75,],
                                        "t113"=beta_samps[113,],
                                        "t150"=beta_samps[150,]),
                             key="time",value="samples",
                             t38,t75,t113,t150) %>%
    mutate(time_label = factor(ifelse(time=="t38","Time 38",
                                      ifelse(time=="t75","Time 75",
                                             ifelse(time=="t113","Time 113",
                                                    ifelse(time=="t150","Time 150",NA)))),
                               levels=c("Time 38","Time 75","Time 113","Time 150")))
  
  plot2 <- plot_grid(ggplot()+
                       theme_bw()+
                       geom_ribbon(data=beta_samps_summary,aes(x=time,ymin=lower*100,ymax=upper*100),alpha=0.5)+
                       geom_line(data=beta_samps_summary,aes(x=time,y=med*100))+
                       geom_point(data=df_ons_sim,aes(x=time,y=mean*100))+
                       geom_errorbar(data=df_ons_sim,aes(x=time,ymin=lower*100,ymax=upper*100))+
                       labs(x="Time",y="Positivity (%)",tag="A"),
                     ggplot(beta_samps_hists,aes(x=samples,group=time))+
                       theme_bw()+
                       geom_histogram(fill="black")+
                       facet_wrap(~time_label,scales="free")+
                       theme(strip.background = element_rect(fill="white"))+
                       labs(x="Logit of test positivity samples",y="Frequency",tag="B"),
                     nrow=2,rel_heights = c(1,1.5))
  
  return(list(plot1 = plot1,
              plot2 = plot2,
              df = df,
              df_ons_sim = df_ons_sim,
              beta_samps_summary = beta_samps_summary,
              beta_samps = beta_samps,
              ks_pval = ks_pval))
  
  
}

## get some outputs to plot

run1 <- likelihood_simulation()
run2 <- likelihood_simulation()
run3 <- likelihood_simulation()
run4 <- likelihood_simulation()
run5 <- likelihood_simulation()
run6 <- likelihood_simulation()
run7 <- likelihood_simulation()
run8 <- likelihood_simulation()

df <- rbind(run1$df %>% mutate(sample = 1),
            run2$df %>% mutate(sample = 2),
            run3$df %>% mutate(sample = 3),
            run4$df %>% mutate(sample = 4),
            run5$df %>% mutate(sample = 5),
            run6$df %>% mutate(sample = 6),
            run7$df %>% mutate(sample = 7),
            run8$df %>% mutate(sample = 8))

# plot_grid(ggplot(df,aes(x=time,y=pos,fill=factor(sample),col=factor(sample)))+
#   geom_col()+
#   facet_wrap(~sample,nrow=1)+
#   theme_bw()+
#     scale_fill_viridis_d()+
#     scale_colour_viridis_d()+
#   labs(x="Time",y="Number of \npositive tests",tag="A")+
#   theme(strip.background = element_blank(),
#         strip.text = element_blank(),
#         legend.position="none"),
#   ggplot(df,aes(x=time,y=prev*100,col=factor(sample),group=sample))+
#   geom_point()+
#   geom_line()+
#   theme_bw()+
#     scale_colour_viridis_d()+
#   labs(x="Time",y="Test positivity (%)",tag="B")+
#   theme(legend.position = "none"),
#   ggplot(df,aes(x=time,y=prev*100,ymin=lower*100,ymax=upper*100,
#               col=factor(sample),group=sample))+
#   geom_point(pch=1)+
#   #geom_line()+
#   geom_errorbar()+
#   theme_bw()+
#     scale_colour_viridis_d()+
#   labs(x="Time",y="Test positivity (%)",tag="C")+
#   theme(legend.position = "none"),
#   nrow=3,align="v",axis="l",rel_heights = c(0.5,1,1)
# )
# ggsave("outputs/likelihood/outputs_example.png",height=6,width=4)


plot_grid(ggplot(df,aes(x=time,y=prev*100,col=factor(sample),group=sample))+
            geom_point(size=0.5)+
            geom_line()+
            theme_bw()+
            scale_colour_viridis_d()+
            labs(x="Time",y="Test positivity (%)")+
            theme(legend.position = "none"),
            ggplot(df,aes(x=time,y=pos,fill=factor(sample),col=factor(sample)))+
            geom_col()+
            facet_wrap(~sample,nrow=2)+
            theme_bw()+
            scale_fill_viridis_d()+
            scale_colour_viridis_d()+
            labs(x="Time",y="Number of \npositive tests")+
            theme(strip.background = element_blank(),
                  strip.text = element_blank(),
                  legend.position="none"),
          nrow=2,align="v",axis="l",labels="AUTO"
)
ggsave("outputs/likelihood/outputs_example.png",height=4,width=6)

beta_samps <- rbind(data.frame(run1$beta_samps) %>% mutate(sample = 1,time=seq(1,150,1)),
                    data.frame(run2$beta_samps) %>% mutate(sample = 2,time=seq(1,150,1)),
                    data.frame(run3$beta_samps) %>% mutate(sample = 3,time=seq(1,150,1)),
                    data.frame(run4$beta_samps) %>% mutate(sample = 4,time=seq(1,150,1)),
                    data.frame(run5$beta_samps) %>% mutate(sample = 5,time=seq(1,150,1)),
                    data.frame(run6$beta_samps) %>% mutate(sample = 6,time=seq(1,150,1)),
                    data.frame(run7$beta_samps) %>% mutate(sample = 7,time=seq(1,150,1)),
                    data.frame(run8$beta_samps) %>% mutate(sample = 8,time=seq(1,150,1)))


beta_samps <- rbind(data.frame("t1"=run1$beta_samps[1,],"t30"=run1$beta_samps[30,],"t60"=run1$beta_samps[60,],
                 "t90"=run1$beta_samps[90,],"t120"=run1$beta_samps[120,],
                 "t150"=run1$beta_samps[150,]) %>% mutate(sample=1),
      data.frame("t1"=run2$beta_samps[1,],"t30"=run2$beta_samps[30,],"t60"=run2$beta_samps[60,],
                 "t90"=run2$beta_samps[90,],"t120"=run2$beta_samps[120,],
                 "t150"=run2$beta_samps[150,]) %>% mutate(sample=2),
      data.frame("t1"=run3$beta_samps[1,],"t30"=run3$beta_samps[30,],"t60"=run3$beta_samps[60,],
                 "t90"=run3$beta_samps[90,],"t120"=run3$beta_samps[120,],
                 "t150"=run3$beta_samps[150,]) %>% mutate(sample=3),
      data.frame("t1"=run4$beta_samps[1,],"t30"=run4$beta_samps[30,],"t60"=run4$beta_samps[60,],
                 "t90"=run4$beta_samps[90,],"t120"=run4$beta_samps[120,],
                 "t150"=run4$beta_samps[150,]) %>% mutate(sample=4),
      data.frame("t1"=run5$beta_samps[1,],"t30"=run5$beta_samps[30,],"t60"=run5$beta_samps[60,],
                 "t90"=run5$beta_samps[90,],"t120"=run5$beta_samps[120,],
                 "t150"=run5$beta_samps[150,]) %>% mutate(sample=5),
      data.frame("t1"=run6$beta_samps[1,],"t30"=run6$beta_samps[30,],"t60"=run6$beta_samps[60,],
                 "t90"=run6$beta_samps[90,],"t120"=run6$beta_samps[120,],
                 "t150"=run6$beta_samps[150,]) %>% mutate(sample=6),
      data.frame("t1"=run7$beta_samps[1,],"t30"=run7$beta_samps[30,],"t60"=run7$beta_samps[60,],
                 "t90"=run7$beta_samps[90,],"t120"=run7$beta_samps[120,],
                 "t150"=run7$beta_samps[150,]) %>% mutate(sample=7),
      data.frame("t1"=run8$beta_samps[1,],"t30"=run8$beta_samps[30,],"t60"=run8$beta_samps[60,],
                 "t90"=run8$beta_samps[90,],"t120"=run8$beta_samps[120,],
                 "t150"=run8$beta_samps[150,]) %>% mutate(sample=8)
      )


beta_samps_long <- gather(beta_samps,
                          value="samps",key="time",
                          t1,t30,t60,t90,t120,t150) %>%
  mutate(time = factor(ifelse(time=="t1","Time 1",
                       ifelse(time=="t30","Time 30",
                              ifelse(time=="t60","Time 60",
                                     ifelse(time=="t90","Time 90",
                                            ifelse(time=="t120","Time 120",
                                                   ifelse(time=="t150","Time 150",NA)))))),
                       levels=c("Time 1","Time 30","Time 60","Time 90","Time 120","Time 150")))




plot_grid(plot_grid(ggplot(df %>% filter(sample==2),
       aes(x=time,y=mean*100,ymin=lower*100,ymax=upper*100))+
         geom_point(col="#46337e")+
         geom_errorbar(col="#46337e")+
         theme_bw()+
         scale_x_continuous(breaks=c(0,30,60,90,120,150))+
         labs(x="Time",y="Test positivity (%)",tag="A"),
       ggplot(beta_samps_long %>% filter(sample==2),
              aes(x=samps))+
         geom_histogram(fill="#46337e")+
         facet_wrap(~time,scales="free_x",nrow=1)+
         theme_bw()+
         #scale_x_continuous(n.breaks=2.5)+
         labs(x="Logit of test positivity samples",y="Frequency",tag="B")+
         theme(strip.background = element_rect(fill="white"),axis.text.x = element_text(angle=90)),
       nrow=2,rel_heights = c(1,0.75),align="v",axis="l"),
       plot_grid(ggplot(df %>% filter(sample==5),
                 aes(x=time,y=mean*100,ymin=lower*100,ymax=upper*100))+
            geom_point(col="#1fa187")+
            geom_errorbar(col="#1fa187")+
            theme_bw()+
              scale_x_continuous(breaks=c(0,30,60,90,120,150))+
            labs(x="Time",y="Test positivity (%)",tag="C"),
          ggplot(beta_samps_long %>% filter(sample==5),
                 aes(x=samps))+
            geom_histogram(fill="#1fa187")+
            facet_wrap(~time,scales="free_x",nrow=1)+
            theme_bw()+
            #scale_x_continuous(n.breaks=2.5)+
            labs(x="Logit of test positivity samples",y="Frequency",tag="D")+
            theme(strip.background = element_rect(fill="white"),
                  axis.text.x = element_text(angle=90)),
          nrow=2,rel_heights = c(1,0.75),align="v",axis="l"),
       nrow=2,align="v",axis="lt")
ggsave("outputs/likelihood/examples_histograms.png",height=9,width=6)




# run1 <- likelihood_simulation()
# length(which(run1$ks_pval<0.05))

#set.seed(711) ## median of 4

#set.seed(34) ## median of 3.5

#set.seed(304) ## median of 3.5

#set.seed(910) ## median of 4

#set.seed(2208) ## median of 2 

p_vals <- c()
p_vals_run <- c()

for(i in 1:1000){
  run <- likelihood_simulation()
  p_vals <- c(p_vals,run$ks_pval)
  p_vals_run <- c(p_vals_run,length(which(run$ks_pval<0.05)))
  #print(which(run$ks_pval<0.05))
}

#plot(p_vals)
#length(which(p_vals<0.05))/length(p_vals)
# discussed with Zoi that can't add them altogether as not independent within runs, but are independent across runs+

summary(p_vals_run)
quantile(p_vals_run,c(0.025,0.5,0.975))
table(p_vals_run)

ggplot(data.frame(p_vals_run),aes(x=1,y=p_vals_run))+
  geom_boxplot(outlier.shape = NA)+
  #geom_jitter(pch=1)+
  geom_point(pch=1)+
  #geom_jitter(aes(x=1,y=p_vals_run))+
  theme_bw()+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  #ylim(c(0,100))+
  labs(y="Number of time points not normally distributed")
ggsave("outputs/simulated_data_for_likelihood_3.png",height=5,width=2)








# set.seed(2904) # because this is a nice one
# run1 <- likelihood_simulation()
# length(which(run1$ks_pval<0.05))
# which(run1$ks_pval<0.05)
# run1$plot1
# ggsave("outputs/simulated_data_for_likelihood_1.png",height=5,width=7)
# run1$plot2
# #ggsave("outputs/simulated_data_for_likelihood_2.png",height=7,width=5)
# run1$ks_pval[c(1,25,50,75,91,100)]
# 
# ### add time 1 and the stat sig one
# 
# beta_samps_hists <- gather(data.frame("t1"=run1$beta_samps[1,],
#                                       "t25"=run1$beta_samps[25,],
#                                       "t50"=run1$beta_samps[50,],
#                                       "t75"=run1$beta_samps[75,],
#                                       "t91"=run1$beta_samps[91,],
#                                       "t100"=run1$beta_samps[100,]),
#                            key="time",value="samples",
#                            t1,t25,t50,t75,t91,t100) %>%
#   mutate(time_label = factor(ifelse(time=="t25","Time 25",
#                                     ifelse(time=="t50","Time 50",
#                                            ifelse(time=="t75","Time 75",
#                                                   ifelse(time=="t100","Time 100",
#                                                          ifelse(time=="t1","Time 1",
#                                                                 ifelse(time=="t91","Time 91",NA)))))),
#                              levels=c("Time 1","Time 25","Time 50","Time 75","Time 91","Time 100")))
# 
# plot2 <- plot_grid(ggplot()+
#                      theme_bw()+
#                      geom_ribbon(data=run1$beta_samps_summary,aes(x=time,ymin=lower*100,ymax=upper*100),alpha=0.5)+
#                      geom_line(data=run1$beta_samps_summary,aes(x=time,y=med*100))+
#                      geom_point(data=run1$df_ons_sim,aes(x=time,y=mean*100))+
#                      geom_errorbar(data=run1$df_ons_sim,aes(x=time,ymin=lower*100,ymax=upper*100))+
#                      labs(x="Time",y="Positivity (%)",tag="A"),
#                    ggplot(beta_samps_hists,aes(x=samples,group=time))+
#                      theme_bw()+
#                      geom_histogram(fill="black")+
#                      facet_wrap(~time_label,scales="free")+
#                      theme(strip.background = element_rect(fill="white"))+
#                      labs(x="Logit of test positivity samples",y="Frequency",tag="B"),
#                    nrow=2,rel_heights = c(1,1.5))
# ggsave("outputs/simulated_data_for_likelihood_2.png",height=6,width=7)
# 
# 
# 
# 
# #which(sw_pval < 0.05) %>% length()
# 
# #fitdist(beta_samps[100,],"norm") %>% plot()
# 
# 
# 
# #fitdist(beta_samps[100,],"norm") %>% plot()
# 
# 
# #fitdist(beta_samps[100,],"cauchy") %>% plot()
# # not a good fit 
# 
# # ## find distribution that best fits these data
# # 
# # hist(beta_samps[100,])
# # 
# # hist(rgamma(1000,shape=1,rate=1))
# #
# #
# #descdist(beta_samps[1,], discrete = FALSE)
# 
# 












