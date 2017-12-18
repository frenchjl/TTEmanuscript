;Sim_start : add to simulation model   
;$SIZES NO=500 LIM6=500
;Sim_end

$PROBLEM    TTE model, Weibull + Emax + covariate (true model)

$INPUT      ID TIME DOSE CL X DV EVID STIME TCEN AUC

;-----------------------------------------------------data description 
; ID, subject identifier 
; TIME, in arbitrary time units (days, months, years, etc )
; DOSE=dose, actual dose 
; CL, actual clearance
; X=x, bianry covariate (0,1)
; DV, DV = 0 (no event observed = right censored (TRUE), DV = 1, an event occured at time = TIME
; EVID, EVID=3 reset the system at time zero/each new ID; EVID=0 indicates an observation
; STIME, flag which indicates if time was observed (STIME=0) or time is simulated (STIME=1)
; TCEN, censored event time (one time per ID)
; AUC, AUCss = dose[i]/CL[i]
;-----------------------------------------------------

$DATA    /pmx_bip/PMx_Playground/gBrown/PostDoc_project_ISoP/Manuscript/Response/NONMEM/DATA/TTEdat_v11.csv IGNORE=@

;Sim_start : remove from simulation model   
IGNORE=(STIME.EQ.2) ; simulated time, ignored for estimation
;IGNORE=(STIME.EQ.0) ; obsTTE time, ignored for simuation 
;Sim_end

$SUBROUTINE ADVAN=6 TOL=9
$MODEL      COMP=(HAZARD)


;--------------------------------------------------------------------------------------
;----------------------Initial estimates
;--------------------------------------------------------------------------------------

$THETA
(0, 2)   ; theta1, Emax (true = 1.663553) 
(0, 2)   ; theta2, EC50, gamma (true = 2)
(0, 0.6) ; theta3, coffiecent of x, (true = 0.6931472) 
(0, 3)   ; theta4, shape parameter delta (true = 3)
(0, 0.1) ; theta5, scale parameter lambda (true = 0.0884997)

$OMEGA 
0  FIX   ; place holder 

;--------------------------------------------------------------------------------------
;----------------------PK parameters
;--------------------------------------------------------------------------------------

$PK
  IF(NEWIND.NE.2) TP=0 

  EMAX = THETA(1)        
  EC50 = THETA(2)          
  COEFF = THETA(3)         
  DELTA = THETA(4)                
  LAMBDA = THETA(5) * EXP(ETA(1))              

  COV = COEFF*X         ; binary covariate effect  

;--------------------------------------------------------------------------------------
;----------------------Differential equation
;----------------------Typical Value weibull h0(t) = detla*(lambda^delta)* t^(delta-1) 
;--------------------------------------------------------------------------------------

$DES
DEL = 1E-6  ; to keep from taking 0**power
IMMED = 1 

BASE = DELTA*(LAMBDA**DELTA)*((T-TP)+DEL)**(DELTA-1)     ; TV baseline weibell hazard 

DEFF = (EMAX * AUC * IMMED) /(EC50 + (AUC * IMMED))      ; effect = Weibull + Emax + immediate  

DADT(1) = BASE * EXP(DEFF + COV)                         ; Hazard model + covariate

;--------------------------------------------------------------------------------------
;----------------------Predicted TTE  
;--------------------------------------------------------------------------------------

$ERROR

 IF(NEWIND.NE.2) OLDCHZ=0   ;reset the cumulative hazard
  CHZ = A(1)-OLDCHZ         ; cumulative hazard from previous time point in data set
                          
  OLDCHZ = A(1)             ;rename old cumulative hazard
  SUR = EXP(-CHZ)           ;survival probability

   DELX = 1E-6
   IMMEDX = 1 
   BASEX = DELTA*(LAMBDA**DELTA)*((TIME-TP)+DELX)**(DELTA-1)      ; TV baseline weibell hazard 
   DEFFX = (EMAX * AUC * IMMEDX) /(EC50 + (AUC * IMMEDX))         ; effect = Weibull + Emax  
 
   HAZNOW= BASEX * EXP(DEFFX + COV)                               ;  Hazard model + covariate
 
IF(DV.EQ.0) Y=SUR	                   ;right censored (prob of survival) = S(t)
IF(DV.EQ.1) Y=SUR*HAZNOW		       ;prob density function of event = h(t)*S(t) 


;------------------------------------------------------------------------
; Residuals where events DV = 1 and censoring DV = 0  
;------------------------------------------------------------------------


; Martingale residual: rM = (1-CENSOR) + log(SURV)  
MARTRES = (DV) - CHZ  

; deviance residual = sign(rM) * SQRT(-2*(rM + (1-CENS)*log(-log(SURV))))
 SIGNRM = 1
 IF (MARTRES < 0) SIGNRM = -1

IF (MDV.EQ.1) THEN
  DEVRES = 0
ELSE
  DEVRES = SIGNRM * SQRT(-2 * (MARTRES + (DV)*LOG(CHZ)))
ENDIF


;Simulation for model evaluation
IF (ICALL.EQ.4) THEN
CALL RANDOM (2,R)
   DV=0
   RTTE = 0
IF(TIME.EQ.24.0) RTTE=1 
IF(R.GE.SUR) THEN 
   DV=1
   RTTE = 1
ENDIF
ENDIF

;Sim_start : add/remove for simulation
$COVARIANCE PRINT=E
$ESTIMATION MAXEVAL=9999 METHOD=COND LAPLACE LIKE PRINT=1 SIGL=9  NSIG=3 MSFO=msfb_6
;$SIMULATION (5988566) (39978 UNIFORM) ONLYSIM NOPREDICTION SUB=100
;Sim_end   

$TABLE FILE=6.tab  ONEHEADER
ID TIME DOSE CL X DV AUC EVID MDV PRED CHZ SUR HAZNOW MARTRES DEVRES NOPRINT NOAPPEND 
  
$TABLE ID TIME SUR DOSE EVID NOPRINT ONEHEADER FILE=sdtab6
$TABLE ID CL AUC NOPRINT ONEHEADER FILE=patab6
$TABLE ID X ONEHEADER FILE=catab6
