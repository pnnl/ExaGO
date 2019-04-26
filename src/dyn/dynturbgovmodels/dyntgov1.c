#include "dyntgov1.h"
#include <private/dynturbgovmodelsimpl.h>

/*
  DYNEventMonitor_TGOV1 - Event monitoring routine for this turbine governor model
*/
PetscErrorCode DYNEventMonitor_TGOV1(DYNTurbgovModel dynturbgov, PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscScalar *fval)
{
  PetscErrorCode ierr;
  DYNTGOV1      tgov1;
  PetscInt       loc;
  PetscScalar    x2,dw;
  DYNGenModel    dyngen;
  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNTurbgovModelGetModelData(dynturbgov,(void**)&tgov1);CHKERRQ(ierr);

  /* Get the location for the first turbine governor variable */
  ierr = DYNTurbgovModelGetFirstVariableLocation(dynturbgov,&loc);CHKERRQ(ierr);

  x2  = x[loc+1];

  ierr = DYNTurbgovModelGetDynGen(dynturbgov,&dyngen);CHKERRQ(ierr);
  ierr = DYNGenModelGetSpeedDeviation(dyngen,t,x,&dw,NULL);CHKERRQ(ierr);

  if(!tgov1->X2atmax) fval[0] = tgov1->VMAX - x2;
  else fval[0] = ((tgov1->Pref - dw)/tgov1->R - x2)/tgov1->T1;

  if(!tgov1->X2atmin) fval[1] = tgov1->VMIN - x2;
  else fval[1] = ((tgov1->Pref - dw)/tgov1->R - x2)/tgov1->T1;

  PetscFunctionReturn(0);
}

/*
  DYNEventPostFunction_TGOV1 - Post event routine for this turbine governor model
*/
PetscErrorCode DYNEventPostFunction_TGOV1(DYNTurbgovModel dynturbgov, PetscInt nevents, PetscInt ev_list[], PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscBool forward_solve, PetscBool *solve_algebraic)
{
  PetscErrorCode ierr;
  DYNTGOV1      tgov1;
  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNTurbgovModelGetModelData(dynturbgov,(void**)&tgov1);CHKERRQ(ierr);

  if(ev_list[0] == 0) { /* Max X2 */
    if(!tgov1->X2atmax) {
      tgov1->X2atmax = 1;
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%4.3f s] Bus %d: Turbine-Governor Model TGOV1 for generator %s hit its max. V limit = %4.3f\n",t,tgov1->bus_i,tgov1->id,tgov1->VMAX);CHKERRQ(ierr);
    } else {
      tgov1->X2atmax = 0;
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%4.3f s] Bus %d: Turbine-Governor Model TGOV1 for generator %s released from its max. V limit = %4.3f\n",t,tgov1->bus_i,tgov1->id,tgov1->VMAX);CHKERRQ(ierr);
    }
  } else { /* Min X2 */
    if(!tgov1->X2atmin) {
      tgov1->X2atmin = 1;
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%4.3f s] Bus %d: Turbine-Governor Model TGOV1 for generator %s hit its min. V limit = %4.3f\n",t,tgov1->bus_i,tgov1->id,tgov1->VMIN);CHKERRQ(ierr);
    } else {
      tgov1->X2atmin = 0;
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%4.3f s] Bus %d: Turbine-Governor Model TGOV1 for generator %s released from its min. V limit = %4.3f\n",t,tgov1->bus_i,tgov1->id,tgov1->VMIN);CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(0);
}

/*
  DYNEventPostDGamma_TGOV1 - Compute derivative of gamma
*/
PetscErrorCode DYNEventPostDGamma_TGOV1(DYNTurbgovModel dynturbgov, PetscInt nevents, PetscInt ev_list[], PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscScalar *dgdx, PetscScalar *dgdp,PetscInt Valoc, PetscInt Vmloc, PetscInt Pgloc, PetscInt Qgloc)
{
  PetscErrorCode ierr;
  DYNTGOV1       tgov1;
  PetscInt       loc,gloc;
  DYNGenModel   dyngen;
  PetscInt      dwloc;
  PetscScalar   dPref_dPg;

  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNTurbgovModelGetModelData(dynturbgov,(void**)&tgov1);CHKERRQ(ierr);

  ierr = DYNTurbgovModelGetFirstVariableLocation(dynturbgov,&loc);CHKERRQ(ierr);

  ierr = DYNTurbgovModelGetDynGen(dynturbgov,&dyngen);CHKERRQ(ierr);
  ierr = DYNGenModelGetFirstVariableLocation(dyngen,&gloc);CHKERRQ(ierr);
  ierr = DYNGenModelGetSpeedDeviationLocation(dyngen,&dwloc);CHKERRQ(ierr); /* Note this is local */

  dPref_dPg = tgov1->R;
  if(ev_list[0] == 0) { /* VMAX */
    if(!tgov1->X2atmax) dgdx[loc+1] = -1;
    else {
      dgdx[loc+1] = -1/tgov1->T1;
      dgdx[gloc+dwloc] = -1/(tgov1->R*tgov1->T1);

      dgdp[Pgloc] = dPref_dPg/(tgov1->R*tgov1->T1);
    }
  } else{ /* VMIN */
    if(!tgov1->X2atmin) dgdx[loc+1] = -1;
    else {
      dgdx[loc+1] = -1/tgov1->T1;
      dgdx[gloc+dwloc] = -1/(tgov1->R*tgov1->T1);

      dgdp[Pgloc] = dPref_dPg/(tgov1->R*tgov1->T1);
    }
  }
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelSetEvent_TGOV1 - Sets the event info for the TGOV1 turbine governor model

  Input Parameters
. dynturbgov - the dynamic turbine governor model

  Output Parameters
+ nmons - the number of event monitors
. direction - the event directions
. terminate - flags for termination when the event is located
. eventfcn  - the function describing the event condition
- posteventfcn - An optional function that gets called when the event is located
*/
PetscErrorCode DYNTurbgovModelSetEvent_TGOV1(DYNTurbgovModel dynturbgov,PetscInt *nmons,PetscInt *direction,PetscBool *terminate,PetscErrorCode (**eventfcn)(DYNTurbgovModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*), PetscErrorCode (**posteventfcn)(DYNTurbgovModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscBool,PetscBool*),PetscErrorCode (**posteventdgamma)(DYNTurbgovModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscInt))
{

  PetscFunctionBegin;
  *nmons = 2;
  direction[0] = -1; direction[1] = 1;
  terminate[0] = terminate[1] = PETSC_FALSE;
  *eventfcn = DYNEventMonitor_TGOV1;
  *posteventfcn = DYNEventPostFunction_TGOV1;
  *posteventdgamma = DYNEventPostDGamma_TGOV1;
  PetscFunctionReturn(0);
}


/*
  DYNTurbgovModelReadData_TGOV1 - Parses the line containining the TGOV1 data and
  populates it in the TGOV1 data structure

  Input Parameters:
+ dynturbgov - the dynamic turbine governor struct
- line - The line in the dyr file

*/
PetscErrorCode DYNTurbgovModelReadData_TGOV1(DYNTurbgovModel dynturbgov,char* line,PetscScalar mbase,PetscScalar sbase)
{
  PetscErrorCode ierr;
  DYNTGOV1 tgov1;
  PetscScalar cratio=mbase/sbase;

  PetscFunctionBegin;
  ierr = DYNTurbgovModelGetModelData(dynturbgov,(void**)&tgov1);CHKERRQ(ierr);
  /* Read the TGOV1 data */
  /* IBUS, ’TGOV1’, I, R, T1, VMAX, VMIN, T2, T3, Dt/ */

  sscanf(line,"%d, 'TGOV1',%[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf", \
	 &tgov1->bus_i,tgov1->id,&tgov1->R,&tgov1->T1, &tgov1->VMAX,&tgov1->VMIN,&tgov1->T2, \
	 &tgov1->T3,&tgov1->DT);
	 
  tgov1->R    *= cratio;
  tgov1->VMAX *= cratio;
  tgov1->VMIN *= cratio;
  tgov1->DT *= cratio;
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelGetBusnumID_TGOV1 - Returns the bus number and ID associated with this turbgov model

  Input Parameters:
. dynturbgov - the dynamic turbgov model object

  Output Parameters:
+ busnum - the bus number 
- turbgovid  - the turbine governor ID
*/
PetscErrorCode DYNTurbgovModelGetBusnumID_TGOV1(DYNTurbgovModel dynturbgov,PetscInt *busnum,char **turbgovid)
{
  PetscErrorCode ierr;
  DYNTGOV1      tgov1;

  PetscFunctionBegin;
  ierr = DYNTurbgovModelGetModelData(dynturbgov,(void**)&tgov1);CHKERRQ(ierr);
  *busnum = tgov1->bus_i;
  *turbgovid  = tgov1->id;
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelDestroy_TGOV1 - Destroys the TGOV1 object data

  Input Parameters
. DYNTurbgovModel - the DYNTurbgovModel object

  Notes:
  Called when DYNTurbgovModelDestroy() is called
*/
PetscErrorCode DYNTurbgovModelDestroy_TGOV1(DYNTurbgovModel dynturbgov)
{
  PetscErrorCode ierr;
  DYNTGOV1      tgov1;

  PetscFunctionBegin;
  ierr = DYNTurbgovModelGetModelData(dynturbgov,(void**)&tgov1);CHKERRQ(ierr);
  ierr = PetscFree(tgov1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelGetNvar_TGOV1 - Returns the number of variables for the TGOV1 model

  Input parameters
. dynturbgov - the dynamic turbine governor object

  Output Parameters
. nvar - number of variables for this model
*/
PetscErrorCode DYNTurbgovModelGetNvar_TGOV1(DYNTurbgovModel dynturbgov,PetscInt *nvar)
{
  PetscFunctionBegin;
  *nvar = DYNTGOV1_nvar;
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelGetSizeof_TGOV1 - Returns the size of TGOV1 struct

  Input Parameters:
. dynturbgov - the dynamic turbgov model
 
  Output Parameters:
. size - size of the TGOV1 object obtained from sizeof()
*/
PetscErrorCode DYNTurbgovModelGetSizeof_TGOV1(DYNTurbgovModel dynturbgov,PetscInt *size)
{
  PetscFunctionBegin;
  *size = sizeof(struct _p_DYNTGOV1);
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelSetInitialConditions_TGOV1 - Sets the initial conditions (x(t0)) for the TGOV1 model

  Input Parameters:
+ dynturbgov - the DYNTurbgovModel object
. VD     - real component of the complex bus voltage
. VQ     - imaginary component of the complex bus voltage

  Output Parameters:
. x - the array for the variables for the bus on which the generator is incident.

  Notes:
   The location for the first variable of this turbine governor model can be obtained by DYNTurbgovModelGetFirstVariableLocation()
   and then the initial conditions populated in the x array using this location information
   x[loc] = 1st variable
   x[loc+1] = 2nd variable
*/
PetscErrorCode DYNTurbgovModelSetInitialConditions_TGOV1(DYNTurbgovModel dynturbgov,PetscScalar VD, PetscScalar VQ, PetscScalar *x)
{
  PetscErrorCode ierr;
  DYNTGOV1 tgov1;
  DYNGenModel dyngen;
  PetscScalar Pmech0,X1,X2;
  PetscInt    loc;

  PetscFunctionBegin;
  ierr = DYNTurbgovModelGetModelData(dynturbgov,(void**)&tgov1);CHKERRQ(ierr);
  ierr = DYNTurbgovModelGetFirstVariableLocation(dynturbgov,&loc);CHKERRQ(ierr);

  ierr = DYNTurbgovModelGetDynGen(dynturbgov,&dyngen);CHKERRQ(ierr);

  /* Get the initial mechanical power Pm0 from the generator */
  ierr = DYNGenModelGetInitialMechanicalPower(dyngen,&Pmech0);CHKERRQ(ierr);

  X2 = Pmech0;
  X1 = (1 - tgov1->T2/tgov1->T3)*X2;

  x[loc]    = X1;
  x[loc+1]  = X2;

  /* Set VMAX and VMIN flags appropriate to X2 values */
  //  if(PetscAbsScalar(tgov1->VMAX - X2) < PETSC_SMALL) tgov1->X2atmax = 1;
  //  if(PetscAbsScalar(tgov1->VMIN - X2) < PETSC_SMALL) tgov1->X2atmin = 1;

  //  if((tgov1->VMAX - X2) <= PETSC_SMALL) tgov1->X2atmax = 1;
  //  if((tgov1->VMIN - X2) >= PETSC_SMALL) tgov1->X2atmin = 1;

  /* Constants -> Pref */
  tgov1->Pref = tgov1->R*Pmech0;

  PetscFunctionReturn(0);
}

PetscErrorCode DYNTurbgovModelSetInitialConditionsP_TGOV1(DYNTurbgovModel dynturbgov,PetscScalar PG,PetscScalar QG,PetscScalar VA,PetscScalar VM,PetscInt PGloc,PetscInt QGloc,PetscInt VAloc,PetscInt VMloc,Mat ICp,PetscInt dynlocglob)
{
  DYNTGOV1       tgov1;
  PetscScalar    x_PG[2];
  PetscInt       rows[2];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DYNTurbgovModelGetModelData(dynturbgov,(void**)&tgov1);CHKERRQ(ierr);

  x_PG[0] = (1. - tgov1->T2/tgov1->T3);
  x_PG[1] = 1.;
  rows[0] = dynlocglob;
  rows[1] = dynlocglob+1;
  ierr = SetMatrixValues(ICp,2,rows,1,&PGloc,x_PG);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelGetEquationTypes_TGOV1 - Gets the number and indices of differential and algebraic equations
                                                   
  Input Parameters:
. dynturbgov - the DYNTurbgovModel object

  Output Parameters:
+ ndiff - number of differential equations
. nalg  - number of algebraic equations
- eqtype - an array of size nvar that has
            eqtype[i] = DIFF_EQ if equation i is a differential
            vartype[i] = ALG_VAR if equation i is an algebraic

  NOTES: The user does not need to create the array vartype
*/
PetscErrorCode DYNTurbgovModelGetEquationTypes_TGOV1(DYNTurbgovModel dynturbgov,PetscInt *ndiff,PetscInt *nalg, PetscInt *eqtype)
{
  PetscErrorCode ierr;
  DYNTGOV1 tgov1;

  PetscFunctionBegin;
  ierr = DYNTurbgovModelGetModelData(dynturbgov,(void**)&tgov1);CHKERRQ(ierr);
  *ndiff = 2;
  *nalg  = 0;
  eqtype[0] = eqtype[1] = DIFF_EQ;

  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelDAERHSJacobianP_TGOV1 - Computes the RHS Jacobian of the turbgov equations w.r.t. parameters

  Input Parameters:
+ dynturbgov     - the dynamic turbine governor model
. t          - the current time
. xdynturbgov    - turbine governor variables
. dynlocglob - starting location of turbine governor variables in the global vector 
- PGloc      - location of PG (generator real power) in the parameter vector

  Output Parameters:
. jacP       - matrix of partial derivatives of DAE equations w.r.t parameter PG

*/ 
PetscErrorCode DYNTurbgovModelDAERHSJacobianP_TGOV1(DYNTurbgovModel dynturbgov,PetscReal t,const PetscScalar *xdynturbgov,Mat jacP,PetscInt dynlocglob,PetscInt PGloc)
{
  PetscErrorCode ierr;
  DYNTGOV1      tgov1;
  PetscInt       loc;
  PetscInt       row,col;
  PetscScalar    value;
  PetscScalar    dPref_dPg;

  PetscFunctionBegin;

  ierr = DYNTurbgovModelGetModelData(dynturbgov,(void**)&tgov1);CHKERRQ(ierr);
  ierr = DYNTurbgovModelGetFirstVariableLocation(dynturbgov,&loc);CHKERRQ(ierr);

  if(!(tgov1->X2atmax || tgov1->X2atmin)) {  
    row = dynlocglob+1;
    col = PGloc;
    dPref_dPg = tgov1->R;
    value = dPref_dPg*(1/(tgov1->R*tgov1->T1));
    ierr = SetMatrixValues(jacP,1,&row,1,&col,&value);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelDAERHSJacobian_TGOV1 - Computes the RHS Jacobian of the TGOV1 DAE equations

  Input Parameters:
+ dynturbgov     - the dynamic turbine governor model
. t          - the current time
. VD         - real-part of complex bus voltage
. VQ         - imaginary-part of complex bus voltage
. x          - variables for the bus on which the generator is incident
. dynlocglob - starting location of turbine governor variables in the global vector 
- V_loc      - global location of VD and VQ V_loc[0] = VD_loc V_loc[1] = VQ_loc


  Output Parameters:
. J          - the Jacobian matrix

*/ 
PetscErrorCode DYNTurbgovModelDAERHSJacobian_TGOV1(DYNTurbgovModel dynturbgov,Mat J,PetscReal t,PetscScalar VD,PetscScalar VQ,PetscScalar *x,PetscInt dynlocglob,PetscInt V_loc[])
{
  PetscErrorCode ierr;
  DYNGenModel    dyngen;
  DYNTGOV1       tgov1;
  PetscInt       loc;
  PetscInt       row[2],col[4];
  PetscScalar    val[4];
  PetscInt       dwlocglob;

  PetscFunctionBegin;
  ierr = DYNTurbgovModelGetModelData(dynturbgov,(void**)&tgov1);CHKERRQ(ierr);
  ierr = DYNTurbgovModelGetFirstVariableLocation(dynturbgov,&loc);CHKERRQ(ierr);

  ierr = DYNTurbgovModelGetDynGen(dynturbgov,&dyngen);CHKERRQ(ierr);

  ierr = DYNGenModelGetSpeedDeviationGlobalLocation(dyngen,&dwlocglob);CHKERRQ(ierr);

  val[0] = val[1] = val[2] = val[3] = 0.0;

  row[0] = dynlocglob;
  col[0] = dynlocglob; col[1] = dynlocglob+1;
  val[0] = -1/tgov1->T3;
  val[1] = (1 - tgov1->T2/tgov1->T3)/tgov1->T3;

  ierr = SetMatrixValues(J,1,row,2,col,val);CHKERRQ(ierr);

  row[0] = dynlocglob+1;
  col[0] = dynlocglob+1; col[1] = dwlocglob;
  val[0] = -1/tgov1->T1;
  val[1] = -1/(tgov1->T1*tgov1->R);
  ierr = SetMatrixValues(J,1,row,2,col,val);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelDAERHSFunction_TGOV1 - Computes the rhs of the DAE function for the TGOV1 model 
                                and also returns the turbine governor currents in network refernce frame

  Input Parameters:
+ dynturbgov - the dynamic turbine governor model
. t      - the current time
. x      - array of the variables for this bus
. VD     - real-part of the bus voltage
- VQ     - imaginary part of the bus voltage

  Output Parameters:
. f      - array of rhs of DAE equations for the tgov1 model 

*/
PetscErrorCode DYNTurbgovModelDAERHSFunction_TGOV1(DYNTurbgovModel dynturbgov,PetscReal t,PetscScalar VD, PetscScalar VQ, PetscScalar *x,PetscScalar *f)
{
  PetscErrorCode ierr;
  DYNTGOV1      tgov1;
  DYNGenModel   dyngen;
  PetscInt       loc;
  PetscScalar    x1,x2,dw;


  PetscFunctionBegin;
  /* Get the model data */
  ierr = DYNTurbgovModelGetModelData(dynturbgov,(void**)&tgov1);CHKERRQ(ierr);

  /* Get the location for the first turbine governor variable */
  ierr = DYNTurbgovModelGetFirstVariableLocation(dynturbgov,&loc);CHKERRQ(ierr);

  /* Get Machine speed deviation */
  ierr = DYNTurbgovModelGetDynGen(dynturbgov,&dyngen);
  ierr = DYNGenModelGetSpeedDeviation(dyngen,t,x,&dw,NULL);CHKERRQ(ierr);

  x1  = x[loc];
  x2  = x[loc+1];

  f[loc] = (-x1 + (1 - tgov1->T2/tgov1->T3)*x2)/tgov1->T3;
  f[loc+1]   = ((tgov1->Pref - dw)/tgov1->R - x2)/tgov1->T1;
  if(tgov1->X2atmax) f[loc+1] = tgov1->VMAX - x2;
  else if(tgov1->X2atmin) f[loc+1] = tgov1->VMIN - x2;

  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelGetMechanicalPower - Returns the mechanical power input for the machine and its partial derivative w.r.t to
                                      the machine speed deviation

  Input Parameters:
+ dynturbgov - the dynamic turbine governor model
. t      - the current time
. xdyn   - array of the variables for this bus

  Output Parameters:
+ Pmech      - Mechanical power input for the machine
. dPmechddw  - Partial derivative of Pmech w.r.t. machine speed deviation
. dPmechdXtgov_num - Number of partial derivatives computed.
. dPmechdXtgov - An array of partial derivatives of the mechanical power w.r.t. the number
                 of turbine governor states given by dPmechdXtgov_num
- dPmechdXtgov_loc - locations for the partial derivatives with reference to the bus variable. Need to
                     get the first variable location using DYNTurbgovModelGetFirstVariableLocation() and
		     then add the displacement.

  Notes - dPmechdXtgov_num, dPmechdXtgov and dPmechdXtgov_loc are optionally returned if the passed-in arrays are not NULL.
          
*/
PetscErrorCode DYNTurbgovModelGetMechanicalPower_TGOV1(DYNTurbgovModel dynturbgov,PetscReal t, PetscScalar *xdyn,PetscScalar *Pmech,PetscScalar *dPmechddw, PetscInt *dPmechdXtgov_num,PetscScalar dPmechdXtgov[],PetscInt dPmechdXtgov_loc[])
{
  PetscErrorCode ierr;
  DYNTGOV1       tgov1;
  DYNGenModel    dyngen;
  PetscInt       loc;
  PetscScalar    x1,x2,dw;


  PetscFunctionBegin;
  /* Get the model data */
  ierr = DYNTurbgovModelGetModelData(dynturbgov,(void**)&tgov1);CHKERRQ(ierr);

  /* Get the location for the first turbine governor variable */
  ierr = DYNTurbgovModelGetFirstVariableLocation(dynturbgov,&loc);CHKERRQ(ierr);

  /* Get Machine speed deviation */
  ierr = DYNTurbgovModelGetDynGen(dynturbgov,&dyngen);
  ierr = DYNGenModelGetSpeedDeviation(dyngen,0,xdyn,&dw,NULL);CHKERRQ(ierr);

  x1  = xdyn[loc];
  x2 = xdyn[loc+1];

  *Pmech = x1 + tgov1->T2/tgov1->T3*x2 - tgov1->DT*dw;
  if(dPmechddw) *dPmechddw = -tgov1->DT;
  if(dPmechdXtgov_num) *dPmechdXtgov_num = 2;
  if(dPmechdXtgov) {
    dPmechdXtgov[0] = 1.0;
    dPmechdXtgov[1] = tgov1->T2/tgov1->T3;
  }
  if(dPmechdXtgov_loc) {
    dPmechdXtgov_loc[0] = loc;
    dPmechdXtgov_loc[1] = loc+1;
  }

  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelCreate_TGOV1 - Class constructor for TGOV1 turbine governor model
*/
PetscErrorCode DYNTurbgovModelCreate_TGOV1(DYNTurbgovModel dynturbgov)
{
  DYNTGOV1 tgov1;
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = PetscCalloc1(1,&tgov1);CHKERRQ(ierr);

  tgov1->X2atmax = tgov1->X2atmin = 0;
 
  dynturbgov->data = (void*)tgov1;

  /* Inherit the ops */
  dynturbgov->ops.readdata             = DYNTurbgovModelReadData_TGOV1;
  dynturbgov->ops.destroy              = DYNTurbgovModelDestroy_TGOV1;
  dynturbgov->ops.getnvar              = DYNTurbgovModelGetNvar_TGOV1;
  dynturbgov->ops.getsizeof            = DYNTurbgovModelGetSizeof_TGOV1;
  dynturbgov->ops.setinitialconditions = DYNTurbgovModelSetInitialConditions_TGOV1;
  dynturbgov->ops.setinitialconditionsp = DYNTurbgovModelSetInitialConditionsP_TGOV1;
  dynturbgov->ops.daerhsfunction       = DYNTurbgovModelDAERHSFunction_TGOV1;
  dynturbgov->ops.daerhsjacobian       = DYNTurbgovModelDAERHSJacobian_TGOV1;
  dynturbgov->ops.daerhsjacobianp      = DYNTurbgovModelDAERHSJacobianP_TGOV1;
  dynturbgov->ops.getequationtypes     = DYNTurbgovModelGetEquationTypes_TGOV1;
  dynturbgov->ops.getmechanicalpower   = DYNTurbgovModelGetMechanicalPower_TGOV1;
  dynturbgov->ops.getbusnumid          = DYNTurbgovModelGetBusnumID_TGOV1;
  dynturbgov->ops.setevent             = DYNTurbgovModelSetEvent_TGOV1;
  PetscFunctionReturn(0);
}


