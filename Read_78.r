#
# R script for working with Project Datasphere data
# Dataset: Colorec_AstraZe_2006_78
# Link: https://www.projectdatasphere.org/projectdatasphere/html/content/78
#
# Author: Jonathan L. French
# Modified: 5 January 2018

# Load libraries
library(haven)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(survminer)
library(muhaz)
library(heplots)

# Set project directory and read data
proj = '~/Documents/Projects/Project Datasphere/'
d = read_sas(file.path(proj,'AllProvidedFiles_78/Project Data Sphere - Cediranib Horizon III/Datasets/rdpsubj.sas7bdat'))

d %>% select(RANDCODE,OSCEN,OSTIM,OSEVENT,TIMETP,PFSEVENT)

# Filter to full analysis set
d = filter(d,FULL_SET==1)

# We will consider the following baseline covariates:
# - WHO Performance Status (0 vs 1+)
# - Baseline albumin (continuous)
# - Baseline age group (categorized)
# - baseline LDH (<> 1.5ULN)

# Create new variables
d = d %>%
  mutate(high_ldh = factor(LDH1_5, levels=0:1,
                          labels=c('Baseline LDH <= 1.5  ULN',
                                   'Baseline LDH > 1.5  ULN')),
         age_grp = factor(AGEGRP, levels=2:4,
                         labels=c('>=18 - <65 yrs',
                                  '>=65 - <75 yrs',
                                  '>=75 yrs')),
         alb_grp = cut(CRFALBCL, breaks=quantile(CRFALBCL),include.lowest=TRUE),
         alp_grp = cut(CRFALPCL, breaks=quantile(CRFALPCL),include.lowest=TRUE),
         vegf_grp = cut(BL_VEGFN, breaks=quantile(BL_VEGFN,na.rm=TRUE),
                        include.lowest=TRUE),
         who_status = factor(WHO_STR,levels=0:1,
                            labels=c('0','1+')),
         OSTIMyrs = OSTIM/365
  )


## Generate K-M Plots for EDA
         
# Overall
fit.os = survfit(Surv(OSTIMyrs,1-OSCEN)~1,data=d)
ggsurvplot(fit.os,break.time.by = 0.5,risk.table = TRUE)


fit.pfs = survfit(Surv(TIMETP,as.numeric(PFSEVENT))~1,data=d)
ggsurvplot(fit.pfs)


# By WHO Performance Status
d %>% count(who_status)
fit.os.who = survfit(Surv(OSTIMyrs,1-OSCEN)~who_status,data=d)
ggsurvplot(fit.os.who, break.time.by = 0.5,risk.table = TRUE)

# By Baseline Albumin (this is categorical, but continuous is available)
d %>% count(d$alb_grp)
fit.os.alb = survfit(Surv(OSTIMyrs,1-OSCEN)~alb_grp,data=d)
ggsurvplot(fit.os.alb, break.time.by = 0.5, risk.table = TRUE)

# By Baseline Alk Phos (this is categorical, but continuous is available)
d %>% count(d$alp_grp)
fit.os.alp = survfit(Surv(OSTIMyrs,1-OSCEN)~alp_grp,data=d)
ggsurvplot(fit.os.alp, break.time.by = 0.5, risk.table = TRUE)

# By age group
d %>% count(age_grp)
fit.os.age = survfit(Surv(OSTIMyrs,1-OSCEN)~age_grp,data=d)
ggsurvplot(fit.os.age, break.time.by = 0.5, risk.table = TRUE)

# By region
d %>% count(REGION)
fit.os.region = survfit(Surv(OSTIMyrs,1-OSCEN)~REGION,data=d)
ggsurvplot(fit.os.region)

# By LDH categorical
d %>% count(high_ldh)
fit.os.ldh = survfit(Surv(OSTIMyrs,1-OSCEN)~high_ldh,data=d)
ggsurvplot(fit.os.ldh, break.time.by = 0.5, risk.table = TRUE)


# By VEGF categorical
d %>% count(vegf_grp)
fit.os.vegf = survfit(Surv(OSTIMyrs,1-OSCEN)~vegf_grp,data=d)
ggsurvplot(fit.os.vegf, break.time.by = 0.5, risk.table = TRUE)

# By Prior Adjuvant therapy
d %>% count(PADJ_STR)
fit.os.padj = survfit(Surv(OSTIMyrs,1-OSCEN)~PADJ_STR,data=d)
ggsurvplot(fit.os.padj, break.time.by = 0.5, risk.table = TRUE)



# Fit models
dfit = filter(d,!is.na(high_ldh) &!is.na(BL_VEGFN))

dists = c('weibull','llogis','lnorm')
base_fits = lapply(dists, function(dist) {
  fargs = list(formula=as.formula('Surv(OSTIMyrs,1-OSCEN)~1'),
               data=dfit,
               dist=dist)
  do.call(flexsurvreg,fargs)
#  flexsurvreg(Surv(OSTIMyrs,1-OSCEN)~1,data=dfit, dist=dist)
})
names(base_fits) = dists
sapply(base_fits,AIC)


# Plot estimated hazard functions and observed
obs_haz = muhaz(times = dfit$OSTIMyrs, delta = 1-dfit$OSCEN)
# 
# est_haz = base_fits %>%
#   purrr::map_df(~as.data.frame(
#     summary(.x, 
#             type='hazard',
#             t=obs_haz$est.grid, 
#             ci=FALSE)
#   ),
#   .id='model')
# 
# est_haz = bind_rows(est_haz, 
#                     data.frame(model='observed',
#                                time=obs_haz$est.grid,
#                                est=obs_haz$haz.est))
# 
# 
# qplot(data=est_haz, x=time,y=est,group=model,geom='line',col=model)


# Create hazard-based VPC
haz_vpc_list = lapply(base_fits, function(model,nsim=500) {
  
  foo = parallel::mclapply(1:nsim, function(i) {
    sims = matrix(simFn(model, dat=mutate(dfit,censored=OSCEN, time=OSTIMyrs)),
                  ncol=2)
    data.frame(irep=i, time=sims[,1],event=sims[,2])
  }) %>% bind_rows()
  
  sim_haz = foo %>% 
    split(.$irep) %>%
    purrr::map(~muhaz(times = .$time, 
                      delta = .$event, 
                      max.time=max(obs_haz$est.grid)))
  
  quants_haz = sim_haz %>%
    purrr::map_df(~data.frame(time=.$est.grid, haz=.$haz.est),.id='irep') %>%
    group_by(time) %>%
    summarize(median = median(haz),
              qlo = quantile(haz, p=0.05),
              qhi = quantile(haz, p=0.95))
  
  p = ggplot(data=quants_haz, aes(x=time,y=median)) +
    geom_ribbon(aes(ymin=qlo,ymax=qhi),fill='red',alpha=0.2) +
    geom_line(col='red') +
    geom_line(data=filter(est_haz,model=='observed'),
              aes(x=time,y=est),col='black') +
    labs(x='Time (years)', y='Hazard') + 
    theme_bw()
  
  return(p)
})  


# Fit base and full models
fit.base = flexsurvreg(Surv(OSTIMyrs,1-OSCEN)~1,data=dfit, dist='llogis')
fit.full = update(fit.base, .~. + who_status+log(ALPCALC)+log(BL_VEGFN)+age_grp+high_ldh)

# Simulate data
nsim=500
foo = parallel::mclapply(1:nsim, function(i) {
  sims = matrix(simFn(fit.full, dat=mutate(dfit,censored=OSCEN, time=OSTIMyrs)),
                ncol=2)
  bind_cols(data.frame(irep=i, time=sims[,1],event=sims[,2]),
            select(dfit,who_status,ALPCALC,BL_VEGFN,age_grp,high_ldh,
                   alp_grp,vegf_grp,alb_grp))
}) %>% bind_rows()


# Make hazard-based VPCs

vpc_hazard <- function(formula,sim_data,obs_data,kern='e',...) {
  
  # Housekeeping
  split_variable = as.character(formula)[2]
  split_levels = levels(obs_data[[split_variable]])
  n_levels = length(split_levels)
  n_sims = length(unique(sim_data$irep))
  
  # Estimate hazard in observed data
  obs_haz = obs_data %>% 
    split(.[[split_variable]]) %>%
    purrr::map(~ muhaz(times = .[['OSTIMyrs']], delta = 1-.[['OSCEN']],
                       bw.method='local',kern=kern))
  
  est_haz = lapply(seq_along(obs_haz), function(i) {
    data.frame(time = obs_haz[[i]]$est.grid,
               haz = obs_haz[[i]]$haz.est,
               group = as.character(split_levels[i]))
  }) %>%
    bind_rows() %>%
    mutate(group = factor(group,levels=split_levels))
  
  max_time = max(purrr::map_dbl(obs_haz, ~max(.$est.grid)))
  
  # Estimate hazard functions in each simulated study
  sim_haz = sim_data %>% 
    split(list(.[['irep']],.[[split_variable]])) %>%
    purrr::map(~muhaz(times = .$time, 
                      delta = .$event, 
                      bw.method = 'local',
                      kern=kern,
                      max.time=max_time))
  
  # Aggregate across simulations
  quants_haz = lapply(1:n_levels, function(i) {
    sim_haz[((i-1)*n_sims+1):(i*n_sims)] %>%
      purrr::map_df(~data.frame(time=.$est.grid, haz=.$haz.est),.id='irep') %>%
      group_by(time) %>%
      summarize(median = median(haz),
                qlo = quantile(haz, p=0.05),
                qhi = quantile(haz, p=0.95)) %>%
      mutate(group = as.character(split_levels[i]))
  }) %>% bind_rows() %>%
    mutate(group = factor(group,levels=split_levels))

  # Plot  
  p = ggplot(data=quants_haz, aes(x=time,y=median)) +
    geom_ribbon(aes(ymin=qlo,ymax=qhi),fill='red',alpha=0.2) +
    geom_line(col='red') +
    geom_line(data=est_haz, aes(x=time,y=haz),col='black') +
    labs(x='Time (years)', y='Hazard') + 
    facet_wrap(~group) +
    theme_bw()
  
  return(p)
  
}

vpc_hazard(~alp_grp,foo,dfit)
vpc_hazard(~who_status,foo,dfit)
vpc_hazard(~high_ldh, foo, dfit)
vpc_hazard(~alb_grp, foo, dfit)



linpred = cbind(rep(1,nrow(dfit)),model.matrix(fit.base)) %*% as.numeric(coef(fit.base)[-1])

cumhaz = Hweibull(dfit$OSTIMyrs,shape = coef(fit.base)['shape'],scale = exp(linpred))
mres = 1-dfit$OSCEN - cumhaz

qplot(x=log(dfit$ALPCALC), y=mres) + geom_smooth()
