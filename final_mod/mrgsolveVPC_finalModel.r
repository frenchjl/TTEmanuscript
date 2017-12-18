#########################################
#
# R code for generating simulated
# TTE data using numerical integration
# of hazard in R.
#
#########################################

library(mrgsolve)
library(dplyr)
library(multidplyr)
library(readr)
library(XML)
#library(flexsurv)
library(survival)
library(survminer)

#'
#'  First, we read the dataset used for fitting model in NONMEM, keeping only the
#'  event records.
#'  
dd = read_csv('./run6/TTEdat_v11.csv',na='.',guess_max = 100000) %>%
  filter(EVID == 0 & STIME == 0) %>%
  rename(DOSE=dose, X=x, AUC=AUCss) 

#'
#' Let's look at a snippet of the data and plot the Kaplan-Meier estimates
#' of the survival curves.
#' 

head(dd)

fit_original = survfit(Surv(TIME,DV)~DOSE, data=dd)
survminer::ggsurvplot(fit_original)

#'
#' Define the model using mrgsolve notation.  We will point to the NONMEM
#' output in the $NMXML block.  For more information about defining models
#' using mrgsolve see the mrgsolve user's guide at
#' http://mrgsolve.github.io/user_guide/.
#' 

#path = './run6'

modelCode <- '
$NMXML
file = "./run6/run6.xml"
theta=TRUE, omega=TRUE,sigma=FALSE


$INIT 
CUMHAZ = 0 

$PARAM 
X=0, AUC=0

$MAIN

double EMAX = THETA1;
double EC50 = THETA2;
double COEFF = THETA3;
double DELTA = THETA4;
double LAMBDA = THETA5 * exp(ETA(1));

double COV = COEFF*X;                       // binary covariate effect  
double DEFF = EMAX * AUC / (EC50 + AUC);    // exposure effect

double DEL = 1E-6  ; // to keep from taking 0**power

//
// ODE block
//
$ODE
double BASE = DELTA*pow(LAMBDA,DELTA)*pow((SOLVERTIME+DEL),(DELTA-1));  // TV baseline weibell hazard 

dxdt_CUMHAZ = BASE * exp(DEFF + COV);  // Hazard model + covariate

//
// Table
//
$TABLE
double SURV = exp(-CUMHAZ); 

$CAPTURE CUMHAZ SURV 
'

#'
#' Next, we compile the mrgsolve model and interrogate the parameter values.
#' 
mrgsolveModel <- mcode(model="tte", code=modelCode,project='./run6/',preclean=TRUE)
param(mrgsolveModel)

#'
#'  Next, we integrate hazard function for each subject in the trial.
#'  We will get the cumulative hazard at a grid of times that are
#'  spaced every 0.5 hour.  Given the time scale of the original data
#'  this grid is sufficiently dense.
survivalCurves = mrgsolveModel %>% 
  idata_set(dd) %>% 
  mrgsim(delta=(1/7)/2, carry.out='X,DOSE') 

head(survivalCurves)

#'
#' Next, we define function for simulating nsims survival outcomes (T)
#' for a single survival curve.
simTimes = function(dd,nsims) {
  uvec = runif(nsims)
  simTTE = sapply(uvec, function(u) approx(x=dd$SURV,y=dd$time,xout=u)$y)
  return(data.frame(ID=dd$ID[1],irep=1:nsims,simTTE=simTTE))
}

#'
#' We will parallelize the simulations using the multidplyr package.
#' To do this, we will set-up a cluster, split the data across clusters, 
#' install dplyr library and copy simTimes function
# to each node

mycluster = create_cluster(cores = parallel::detectCores())

survivalCurves = survivalCurves %>%
  partition(ID,cluster= mycluster) %>%
  cluster_library('dplyr') %>%
  cluster_assign_value('simTimes',simTimes)
  
#'
#' Finally, we will simulate 1000 outcomes for each patient,
#' then recombine and rearrange by simulation number (irep) and ID.
#'
simulatedTimes =  survivalCurves %>%
  do(simTimes(.,nsims=1000)) %>%
  collect() %>%
  arrange(irep,ID)

#'
#' Now, we need to simulate the censoring times.  TO do this, we will
#' use a cubic spline estimate of the cumulative hazard (and survival
#' function) for the censoring process.  This is effectively a parametric
#' model which is a close approximation to a Kaplan-Meier estimate
#' of the survival function.  
#' 
#' For this example, we choose the number of internal knots for the spline 
#' through trial-and-error.  More formal methods for choosing the number
#' and location of the spline knots could have be used.  For this example,
#' four internal knots provides a good fit to the censoring distribution.
#' 

censorModel = flexsurv::flexsurvspline(Surv(TIME,1-DV)~1, data=dd, k=4)
plot(censorModel, xlab='Time', ylab='Survival function for censoring process')

#'
#' We will define an administrative censoring time to be the last observed
#' censoring time.  Subjects whose simulated censoring time is beyond the
#' last observed censoring time will be set to this administrative censoring
#' time.
#'  
adminCensorTime = max(dd$TIME[dd$DV==0])

#'
#' Next, we extract the estimated censoring survival values at a grid of 
#' times.  We will interpolate between these values and use the inverse CDF
#' method to simulate censoring times.
#' 
survEst = summary(censorModel, 
                  ci=FALSE, 
                  t=seq(0,max(dd$TIME),length=200)) %>%
  as.data.frame() 
ucens = runif(nrow(simulatedTimes))

simulatedTimes <- simulatedTimes %>%
  ungroup() %>%
  mutate(censorTime = approxfun(x=survEst$est, y=survEst$time)(ucens),
         censorTime = ifelse(is.na(censorTime), adminCensorTime, censorTime),
         simDV = as.numeric(simTTE < censorTime),
         simTIME = ifelse(simDV==1, simTTE, censorTime)) %>%
  left_join(dd %>% select(ID,DOSE,X))

head(simulatedTimes)

#'
#' As an example, we plot the simulated survival times from the first
#' simulated study.
#' 
fit1 = survfit(Surv(simTIME,simDV) ~ DOSE, data=filter(simulatedTimes,irep==1))

survminer::ggsurvplot(fit1)

save(simulatedTimes, file='simulatedTimesForVPC.Rdata')
write_csv(simulatedTimes,path = './simulatedTimesForVPC.csv',na='.',col_names = TRUE)
