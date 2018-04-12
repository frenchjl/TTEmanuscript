setwd('~/Documents/My Papers/2016-TTE Diagnostics/TTEmanuscript')
library(readr)
library(ggplot2)
library(dplyr)
library(survival)
library(purrr)

rm(list=ls())

odat <- read_csv('final_mod/run6/TTEdat_v11.csv') %>%
  filter(STIME==0)

ofit <- coxph(Surv(TIME,DV)~as.factor(dose)+x, data=odat)

load('TTE_simulated_data_model_runs1_8.RData')


sim_coef <- lapply(list(sim1,sim3,sim5), function(simdat){
  simdat %>% 
  split(.$simNumber) %>%
  map( ~ coef(coxph(Surv(TIME,DV)~as.factor(DOSE)+GENDER,data=.)))
})


sim_coef_tall <- lapply(sim_coef, function(ests) {
  data.frame(dose=rep(c(1,3,10),each=length(ests)),
             obs_log_hr = rep(coef(ofit)[1:3],each=length(ests)),
             sim_loghr = as.vector(t(bind_rows(ests))[,1:3]))
}) %>%
  bind_rows(.id='model') %>%
  ungroup()
sim_coef_tall = sim_coef_tall %>%
  mutate(model = factor(model,labels=c('Model 1','Model 3','Model 5')),
         group = factor(paste(sim_coef_tall$dose,'mg vs placebo'),
                        levels=paste(sort(unique(dose)),'mg vs placebo')))

log_hr_plot = sim_coef_tall %>%
  ggplot(aes(x=sim_loghr)) +
  geom_histogram(fill='white',col='black',bins=30) +
  #geom_vline(aes(xintercept=obs_log_hr),col='red') +
  facet_grid(model~group, labeller=label_value) +
  theme_bw() +
  labs(x='log hazard ratio') 

log_hr_plot

ggsave(log_hr_plot,height=8,width=6,filename = 'logHRplot_models135.png')
