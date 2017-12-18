
$NMXML
run = "run6"
project=path
theta=TRUE, omega=TRUE,sigma=FALSE


$INIT 
CUMHAZ = 0 

$PARAM 
DOSE=0, CLI=10, X1=0, AUC=0

$MAIN


double EMAX = THETA1;
double EC50 = THETA2;
double COEFF = THETA3;
double DELTA = THETA4;
double LAMBDA = THETA5 * exp(ETA1);

double DEFF = EMAX * AUC / (EC50 + AUC);

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
 
