## check if models converged or not 

library(tidyverse)
library(rstan)

# areaName <- "England"
# 
# fit20 <- readRDS("outputs/objects/fit_20_England.RDS")
# 
# summary(fit20)$summary
# 
# 
# 
# data.frame("parameter"=rownames(summary(fit20)$summary),summary(fit20)$summary)


areas <- c("England","Scotland","Wales","Northern Ireland")
knot_prop <- seq(20,40,10)

outputs <- c()

for(area in areas){
  for(k in knot_prop){
    
    #fit_it <- readRDS(paste0("outputs/objects/fit_",k,"_",area,".RDS"))
    # outputs_it <- data.frame("parameter"=rownames(summary(fit_it)$summary),summary(fit_it)$summary) %>%
    #   mutate(areaName = area,
    #          knot_prop = k)
    # rownames(outputs_it) <- c()
    outputs_it <- readRDS(paste0("outputs/objects/fit_convergence_",k,"_",area,".RDS")) %>%
      mutate(areaName = area,
             knot_prop = k)
    
    outputs <- rbind(outputs,outputs_it)
    
  }
}


summary(outputs)



ggplot(outputs,aes(x=knot_prop,y=Rhat,col=parameter))+
  geom_point()+
  theme(legend.position="none")



ggplot(outputs %>% filter(parameter=="tau"),aes(x=knot_prop,y=Rhat,col=areaName))+
  geom_point()+
  geom_line()


ggplot(outputs,aes(x=Rhat,y=n_eff))+
  geom_point()#+
  #geom_smooth(method="lm")


outputs %>% arrange(n_eff)





