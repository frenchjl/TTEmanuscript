setwd('~/Documents/My Papers/2016-TTE Diagnostics/TTEmanuscript')
library(readr)
library(ggplot2)
library(dplyr)
library(survival)
library(purrr)

odat <- read_csv('final_mod/run6/TTEdat_v11.csv') %>%
  filter(STIME==0)

ofit <- coxph(Surv(TIME,DV)~as.factor(dose)+x, data=odat)



simdat <- read_csv('final_mod/simulatedTimesForVPC.csv')

sim_coef <- simdat %>% 
  split(.$irep) %>%
  map( ~ coef(coxph(Surv(simTIME,simDV)~as.factor(DOSE)+X,data=.)))

sim_coef_tall <- data.frame(dose=rep(c(1,3,10),each=length(sim_coef)),
                            obs_log_hr = rep(coef(ofit)[1:3],each=length(sim_coef)),
                            sim_loghr = as.vector(t(bind_rows(sim_coef))[,1:3]))

sim_coef_tall %>%
  ggplot(aes(x=sim_loghr)) +
  geom_histogram(fill='white',col='black',bins=30) +
  geom_vline(aes(xintercept=obs_log_hr),col='red') +
  facet_grid(dose~., labeller='label_both') +
  theme_bw() +
  labs(x='log hazard ratio')
