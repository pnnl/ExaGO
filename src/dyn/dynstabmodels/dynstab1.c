#include "dynstab1.h"
#include <private/dynstabmodelsimpl.h>

/*
  DYNEventMonitor_STAB1 - Event monitoring routine for this stabiter model
*/
PetscErrorCode DYNEventMonitor_STAB1(DYNStabModel dynstab, PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscScalar *fval)
{
  PetscErrorCode ierr;
  DYNSTAB1      stab1;
  DYNGenModel    dyngen;
  PetscInt       loc;
  PetscScalar   x3;
  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNStabModelGetModelData(dynstab,(void**)&stab1);CHKERRQ(ierr);

  ierr = DYNStabModelGetDynGen(dynstab,&dyngen);CHKERRQ(ierr);
  ierr = DYNStabModelGetFirstVariableLocation(dynstab,&loc);CHKERRQ(ierr);

  x3 = x[loc+2];

  fval[0] = stab1->Hlim-x3;  
  fval[1] = -stab1->Hlim-x3; 
  PetscFunctionReturn(0);
}

/*
  DYNEventPostFunction_STAB1 - Post event routine for this stabiter model
*/
PetscErrorCode DYNEventPostFunction_STAB1(DYNStabModel dynstab, PetscInt nevents, PetscInt ev_list[], PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscBool forward_solve, PetscBool *solve_algebraic)
{
  PetscErrorCode ierr;
  DYNSTAB1      stab1;
  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNStabModelGetModelData(dynstab,(void**)&stab1);CHKERRQ(ierr);
  if(ev_list[0] == 0) { // Max x3 
    if(!stab1->VOTHSGatmax) stab1->VOTHSGatmax = 1;
    else stab1->VOTHSGatmax = 0;
  } else { // Min x3 
    if(!stab1->VOTHSGatmin) stab1->VOTHSGatmin = 1;
    else stab1->VOTHSGatmin = 0;
  }
  PetscFunctionReturn(0);
}

/*
  DYNEventPostDGamma_STAB1 - Compute derivative of gamma
*/
PetscErrorCode DYNEventPostDGamma_STAB1(DYNStabModel dynstab, PetscInt nevents, PetscInt ev_list[], PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscScalar *dgdx, PetscScalar *dgdp,PetscInt Valoc, PetscInt Vmloc, PetscInt Pgloc, PetscInt Qgloc)
{
  PetscErrorCode ierr;
  DYNSTAB1       stab1;
  PetscInt       loc;
  DYNGenModel   dyngen;

  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNStabModelGetModelData(dynstab,(void**)&stab1);CHKERRQ(ierr);

  ierr = DYNStabModelGetFirstVariableLocation(dynstab,&loc);CHKERRQ(ierr);

  ierr = DYNStabModelGetDynGen(dynstab,&dyngen);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNStabModelSetEvent_STAB1 - Sets the event info for the STAB1 stabiter model

  Input Parameters
. dynstab - the dynamic stabiter model

  Output Parameters
+ nmons - the number of event monitors
. direction - the event directions
. terminate - flags for termination when the event is located
. eventfcn  - the function describing the event condition
- posteventfcn - An optional function that gets called when the event is located
*/
PetscErrorCode DYNStabModelSetEvent_STAB1(DYNStabModel dynstab,PetscInt *nmons,PetscInt *direction,PetscBool *terminate,PetscErrorCode (**eventfcn)(DYNStabModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*), PetscErrorCode (**posteventfcn)(DYNStabModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscBool,PetscBool*),PetscErrorCode (**posteventdgamma)(DYNStabModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscInt))
{

  PetscFunctionBegin;
  *nmons = 2;
  //  direction[0] = 0; direction[1] = 0;
  direction[0] = -1; direction[1] = 1;
  terminate[0] = terminate[1] = PETSC_FALSE;
  *eventfcn = DYNEventMonitor_STAB1;
  *posteventfcn = DYNEventPostFunction_STAB1;
  *posteventdgamma = DYNEventPostDGamma_STAB1;
  PetscFunctionReturn(0);
}


/*
  DYNStabModelReadData_STAB1 - Parses the line containining the STAB1 data and
  populates it in the STAB1 data structure

  Input Parameters:
+ dynstab - the dynamic stabiter struct
- line - The line in the dyr file

*/
PetscErrorCode DYNStabModelReadData_STAB1(DYNStabModel dynstab,char* line)
{
  PetscErrorCode ierr;
  DYNSTAB1 stab1;

  PetscFunctionBegin;
  /*   IBUS, ’STAB1’, I, K/T, T, T1/T3, T3, T2/T4, T4,HLIM/ */
  ierr = DYNStabModelGetModelData(dynstab,(void**)&stab1);CHKERRQ(ierr);
	 
  /* Read the STAB1 data */
  sscanf(line,"%d, 'STAB1',%[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf", \
	 &stab1->bus_i,stab1->id,&stab1->K_T,&stab1->T,&stab1->T1_T3,&stab1->T3,\
	 &stab1->T2_T4,&stab1->T4,&stab1->Hlim);

  PetscFunctionReturn(0);
}

/*
  DYNStabModelGetBusnumID_STAB1 - Returns the bus number and ID associated with this stab model

  Input Parameters:
. dynstab - the dynamic stab model object

  Output Parameters:
+ busnum - the bus number 
- stabid  - the stabiter ID
*/
PetscErrorCode DYNStabModelGetBusnumID_STAB1(DYNStabModel dynstab,PetscInt *busnum,char **stabid)
{
  PetscErrorCode ierr;
  DYNSTAB1      stab1;

  PetscFunctionBegin;
  ierr = DYNStabModelGetModelData(dynstab,(void**)&stab1);CHKERRQ(ierr);
  *busnum = stab1->bus_i;
  *stabid  = stab1->id;
  PetscFunctionReturn(0);
}

/*
  DYNStabModelDestroy_STAB1 - Destroys the STAB1 object data

  Input Parameters
. DYNStabModel - the DYNStabModel object

  Notes:
  Called when DYNStabModelDestroy() is called
*/
PetscErrorCode DYNStabModelDestroy_STAB1(DYNStabModel dynstab)
{
  PetscErrorCode ierr;
  DYNSTAB1      stab1;

  PetscFunctionBegin;
  ierr = DYNStabModelGetModelData(dynstab,(void**)&stab1);CHKERRQ(ierr);
  ierr = PetscFree(stab1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNStabModelGetNvar_STAB1 - Returns the number of variables for the STAB1 model

  Input parameters
. dynstab - the dynamic stabiter object

  Output Parameters
. nvar - number of variables for this model
*/
PetscErrorCode DYNStabModelGetNvar_STAB1(DYNStabModel dynstab,PetscInt *nvar)
{
  PetscFunctionBegin;
  *nvar = DYNSTAB1_nvar;
  PetscFunctionReturn(0);
}

/*
  DYNStabModelGetSizeof_STAB1 - Returns the size of STAB1 struct

  Input Parameters:
. dynstab - the dynamic stab model
 
  Output Parameters:
. size - size of the STAB1 object obtained from sizeof()
*/
PetscErrorCode DYNStabModelGetSizeof_STAB1(DYNStabModel dynstab,PetscInt *size)
{
  PetscFunctionBegin;
  *size = sizeof(struct _p_DYNSTAB1);
  PetscFunctionReturn(0);
}

/*
  DYNStabModelSetInitialConditions_STAB1 - Sets the initial conditions (x(t0)) for the STAB1 model

  Input Parameters:
+ dynstab - the DYNStabModel object
. VD     - real component of the complex bus voltage
. VQ     - imaginary component of the complex bus voltage

  Output Parameters:
. x - the array for the variables for the bus on which the generator is incident.

  Notes:
   The location for the first variable of this stabiter model can be obtained by DYNStabModelGetFirstVariableLocation()
   and then the initial conditions populated in the x array using this location information
   x[loc] = 1st variable
   x[loc+1] = 2nd variable
*/
PetscErrorCode DYNStabModelSetInitialConditions_STAB1(DYNStabModel dynstab,PetscScalar VD, PetscScalar VQ, PetscScalar *x)
{
  PetscErrorCode ierr;
  DYNSTAB1 stab1;
  PetscInt loc;

  PetscFunctionBegin;
  ierr = DYNStabModelGetModelData(dynstab,(void**)&stab1);CHKERRQ(ierr);
  ierr = DYNStabModelGetFirstVariableLocation(dynstab,&loc);CHKERRQ(ierr);

  x[loc] = x[loc+1] = x[loc+2] = 0.0;

  PetscFunctionReturn(0);
}

/*
  DYNStabModelGetEquationTypes_STAB1 - Gets the number and indices of differential and algebraic equations
                                                   
  Input Parameters:
. dynstab - the DYNStabModel object

  Output Parameters:
+ ndiff - number of differential equations
. nalg  - number of algebraic equations
- eqtype - an array of size nvar that has
            eqtype[i] = DIFF_EQ if equation i is a differential
            vartype[i] = ALG_VAR if equation i is an algebraic

  NOTES: The user does not need to create the array vartype
*/
PetscErrorCode DYNStabModelGetEquationTypes_STAB1(DYNStabModel dynstab,PetscInt *ndiff,PetscInt *nalg, PetscInt *eqtype)
{
  PetscErrorCode ierr;
  DYNSTAB1 stab1;

  PetscFunctionBegin;
  ierr = DYNStabModelGetModelData(dynstab,(void**)&stab1);CHKERRQ(ierr);
  *ndiff = 3;
  *nalg  = 0;
  eqtype[0] = eqtype[1] = eqtype[2] = DIFF_EQ;

  PetscFunctionReturn(0);
}

/*
  DYNStabModelDAERHSJacobianP_STAB1 - Computes the RHS Jacobian of the stab equations w.r.t. parameters

  Input Parameters:
+ dynstab     - the dynamic stabiter model
. t          - the current time
. xdynstab    - stabiter variables
. dynlocglob - starting location of stabiter variables in the global vector 
- PGloc      - location of PG (generator real power) in the parameter vector

  Output Parameters:
. jacP       - matrix of partial derivatives of DAE equations w.r.t parameter PG

*/ 
PetscErrorCode DYNStabModelDAERHSJacobianP_STAB1(DYNStabModel dynstab,PetscReal t,const PetscScalar *xdynstab,Mat jacP,PetscInt dynlocglob,PetscInt PGloc)
{
  PetscErrorCode ierr;
  DYNSTAB1      stab1;

  PetscFunctionBegin;

  ierr = DYNStabModelGetModelData(dynstab,(void**)&stab1);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNStabModelDAERHSJacobian_STAB1 - Computes the RHS Jacobian of the STAB1 DAE equations

  Input Parameters:
+ dynstab     - the dynamic stabiter model
. t          - the current time
. VD         - real-part of complex bus voltage
. VQ         - imaginary-part of complex bus voltage
. x          - variables for the bus on which the generator is incident
. dynlocglob - starting location of stabiter variables in the global vector 
- V_loc      - global location of VD and VQ V_loc[0] = VD_loc V_loc[1] = VQ_loc


  Output Parameters:
. J          - the Jacobian matrix

*/ 
PetscErrorCode DYNStabModelDAERHSJacobian_STAB1(DYNStabModel dynstab,Mat J,PetscReal t,PetscScalar VD,PetscScalar VQ,PetscScalar *x,PetscInt dynlocglob,PetscInt V_loc[])
{
  PetscErrorCode ierr;
  DYNSTAB1       stab1;
  PetscInt       row[2],col[4];
  PetscScalar    val[4],dw;
  DYNGenModel    dyngen;
  PetscInt       dwlocglob;

  PetscFunctionBegin;

  ierr = DYNStabModelGetModelData(dynstab,(void**)&stab1);CHKERRQ(ierr);
  ierr = DYNStabModelGetDynGen(dynstab,&dyngen);CHKERRQ(ierr);

  ierr = DYNGenModelGetSpeedDeviation(dyngen,t,x,&dw,NULL);CHKERRQ(ierr);
  ierr = DYNGenModelGetSpeedDeviationGlobalLocation(dyngen,&dwlocglob);CHKERRQ(ierr);

  val[0] = val[1] = val[2] = val[3] = 0.0;

  row[0] = dynlocglob;
  col[0] = dynlocglob; col[1] = dwlocglob;
  val[0] = -1.0/stab1->T;
  val[1] = -stab1->K_T/stab1->T;

  ierr = SetMatrixValues(J,1,row,2,col,val);CHKERRQ(ierr);

  row[0] = dynlocglob+1;
  col[0] = dynlocglob; col[1] = dynlocglob + 1; col[2] = dwlocglob;
  val[0] = (1/stab1->T3)*(1 - stab1->T1_T3);
  val[1] = -1/stab1->T3;
  val[2] = (1/stab1->T3)*(1 - stab1->T1_T3)*stab1->K_T;


  ierr = SetMatrixValues(J,1,row,3,col,val);CHKERRQ(ierr);


  row[0] = dynlocglob+2;
  col[0] = dynlocglob; col[1] = dynlocglob + 1; col[2] = dynlocglob + 2; col[3] = dwlocglob;

  val[0] = (1/stab1->T4)*(1 - stab1->T2_T4)*stab1->T1_T3;
  val[1] = (1/stab1->T4)*(1 - stab1->T2_T4);
  val[2] = -1/stab1->T4;
  val[3] = (1/stab1->T4)*(1 - stab1->T2_T4)*stab1->T1_T3*stab1->K_T;

  ierr = SetMatrixValues(J,1,row,4,col,val);CHKERRQ(ierr);


  PetscFunctionReturn(0);
}

/*
  DYNStabModelDAERHSFunction_STAB1 - Computes the rhs of the DAE function for the STAB1 model 
                                and also returns the stabiter currents in network refernce frame

  Input Parameters:
+ dynstab - the dynamic stabiter model
. t      - the current time
. x      - array of the variables for this bus
. VD     - real-part of the bus voltage
- VQ     - imaginary part of the bus voltage

  Output Parameters:
. f      - array of rhs of DAE equations for the stab1 model 

*/
PetscErrorCode DYNStabModelDAERHSFunction_STAB1(DYNStabModel dynstab,PetscReal t,PetscScalar VD, PetscScalar VQ, PetscScalar *x,PetscScalar *f)
{
  PetscErrorCode ierr;
  DYNSTAB1      stab1;
  DYNGenModel   dyngen;
  PetscScalar   dw;
  PetscInt      loc;
  PetscScalar   x1,x2,x3;


  PetscFunctionBegin;
  /* Get the model data */
  ierr = DYNStabModelGetModelData(dynstab,(void**)&stab1);CHKERRQ(ierr);

  /* Get the location for the first stabilizer variable */
  ierr = DYNStabModelGetFirstVariableLocation(dynstab,&loc);CHKERRQ(ierr);

  x1 = x[loc];
  x2 = x[loc+1];
  x3 = x[loc+2];

  /* Get Machine speed deviation */
  ierr = DYNStabModelGetDynGen(dynstab,&dyngen);
  ierr = DYNGenModelGetSpeedDeviation(dyngen,t,x,&dw,NULL);CHKERRQ(ierr);

  f[loc]   = -x1/stab1->T -stab1->K_T*dw/stab1->T;
  f[loc+1] = -x2/stab1->T3 + (1/stab1->T3)*(1 - stab1->T1_T3)*(x1 + stab1->K_T*dw);
  f[loc+2] = -x3/stab1->T4 + (1/stab1->T4)*(1 - stab1->T2_T4)*(x2 + stab1->T1_T3*(x1 + stab1->K_T*dw));

  PetscFunctionReturn(0);
}

/*
  DYNStabModelGetVOTHSG_STAB1 - Returns the stabilizer signal for the exciter and optionally returns partial derivative w.r.t to
                                      the stabilizer variables and speed deviation

  Input Parameters:
+ dynstab - the dynamic stabilizer model
. t      - the current time
. xdyn   - array of the variables for this bus

  Output Parameters:
+ VOTHSG           - Stabilizer signal output
. dVOTHSGdXdyn_num - Number of partial derivatives of VOTHSG w.r.t. the dynamic variables at the bus 
. dVOTHSGdXdyn     - An array of partial derivatives of the stabilizer output signal VOTHSG w.r.t. the dynamic
                     variables at the bus
- dVOTHSGdXdyn_loc - Location of partial derivatives of the stabilizer output signal VOTHSG w.r.t. the dynamic
                     variables at the bus.  Need to get the first variable location using DYNStabModelGetFirstVariableLocation()
		     and then add the displacement.

  Notes - dVOTHSGdXdyn_num, dVOTHSGdXdyn and dVOTHSGdXdyn_loc are optionally returned if the passed-in arrays are not NULL.
          
*/
PetscErrorCode DYNStabModelGetVOTHSG_STAB1(DYNStabModel dynstab,PetscReal t, PetscScalar *xdyn,PetscScalar *VOTHSG,PetscInt *dVOTHSGdXdyn_num,PetscScalar dVOTHSGdXdyn[],PetscInt dVOTHSGdXdyn_loc[])
{
  PetscErrorCode ierr;
  DYNSTAB1       stab1;
  DYNGenModel    dyngen;
  PetscInt       loc,dwloc;
  PetscScalar    x1,x2,x3,dw;


  PetscFunctionBegin;
  /* Get the model data */
  ierr = DYNStabModelGetModelData(dynstab,(void**)&stab1);CHKERRQ(ierr);

  /* Get the location for the first stabilizer model variable */
  ierr = DYNStabModelGetFirstVariableLocation(dynstab,&loc);CHKERRQ(ierr);

  /* Get Machine speed deviation */
  ierr = DYNStabModelGetDynGen(dynstab,&dyngen);
  ierr = DYNGenModelGetSpeedDeviation(dyngen,0,xdyn,&dw,NULL);CHKERRQ(ierr);
  ierr = DYNGenModelGetSpeedDeviationLocation(dyngen,&dwloc);CHKERRQ(ierr);

  x1  = xdyn[loc];
  x2  = xdyn[loc+1];
  x3  = xdyn[loc+2];

  if(stab1->VOTHSGatmax) *VOTHSG = stab1->Hlim;
  else if(stab1->VOTHSGatmin) *VOTHSG = -stab1->Hlim;
  else *VOTHSG = x3 + stab1->T2_T4*(x2 + stab1->T1_T3*(x1 + stab1->K_T*dw));

  if(dVOTHSGdXdyn_num) {
    if(stab1->VOTHSGatmax || stab1->VOTHSGatmin) {
      *dVOTHSGdXdyn_num = 0;
      PetscFunctionReturn(0);
    } else {
      *dVOTHSGdXdyn_num = 4;
    }
  }
    
  if(dVOTHSGdXdyn_loc) {
    dVOTHSGdXdyn_loc[0] = dwloc; dVOTHSGdXdyn_loc[1] = loc;  dVOTHSGdXdyn_loc[2] = loc+1;  dVOTHSGdXdyn_loc[3] = loc+2;
  }

  if(dVOTHSGdXdyn) {
    dVOTHSGdXdyn[0] = stab1->T2_T4*stab1->T1_T3*stab1->K_T; /* dVOTHSG_ddw */
    dVOTHSGdXdyn[1] = stab1->T2_T4*stab1->T1_T3; /* dVOTHSG_dx1 */
    dVOTHSGdXdyn[2] = stab1->T2_T4; /* dVOTHSG_dx2 */
    dVOTHSGdXdyn[3] = 1.0;
  }

  PetscFunctionReturn(0);
}

/*
  DYNStabModelCreate_STAB1 - Class constructor for STAB1 stabiter model
*/
PetscErrorCode DYNStabModelCreate_STAB1(DYNStabModel dynstab)
{
  DYNSTAB1 stab1;
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = PetscCalloc1(1,&stab1);CHKERRQ(ierr);

  dynstab->data = (void*)stab1;

  stab1->VOTHSGatmax = stab1->VOTHSGatmin = 0;

  /* Inherit the ops */
  dynstab->ops.readdata             = DYNStabModelReadData_STAB1;
  dynstab->ops.destroy              = DYNStabModelDestroy_STAB1;
  dynstab->ops.getnvar              = DYNStabModelGetNvar_STAB1;
  dynstab->ops.getsizeof            = DYNStabModelGetSizeof_STAB1;
  dynstab->ops.setinitialconditions = DYNStabModelSetInitialConditions_STAB1;
  dynstab->ops.daerhsfunction       = DYNStabModelDAERHSFunction_STAB1;
  dynstab->ops.daerhsjacobian       = DYNStabModelDAERHSJacobian_STAB1;
  dynstab->ops.daerhsjacobianp      = DYNStabModelDAERHSJacobianP_STAB1;
  dynstab->ops.getequationtypes     = DYNStabModelGetEquationTypes_STAB1;
  dynstab->ops.getvothsg            = DYNStabModelGetVOTHSG_STAB1;
  dynstab->ops.getbusnumid          = DYNStabModelGetBusnumID_STAB1;
  dynstab->ops.setevent             = DYNStabModelSetEvent_STAB1;
  PetscFunctionReturn(0);
}


