#include "dyncv.h"
#include <private/dyngenmodelsimpl.h>

/*
  DYNEventMonitor_Cv - Event monitoring routine for this generator model
*/
PetscErrorCode DYNEventMonitor_Cv(DYNGenModel dyngen, PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscScalar *fval)
{
  PetscErrorCode ierr;
  DYNCv          cv;

  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNEventPostDGamma_Cv - Compute derivative of gamma
*/
PetscErrorCode DYNEventPostDGamma_Cv(DYNGenModel dyngen, PetscInt nevents, PetscInt ev_list[], PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscScalar *dgdx, PetscScalar *dgdp,PetscInt Valoc, PetscInt Vmloc, PetscInt Pgloc, PetscInt Qgloc)
{
  PetscErrorCode ierr;
  DYNCv      cv;

  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetFieldCurrent_Cv - Returns the field current Ifd

  Input Parameters:
. dyngen - the DYNGenModel object

  Output Parameters:
. Ifd - Field current Ifd
*/
PetscErrorCode DYNGenModelGetFieldCurrent_Cv(DYNGenModel dyngen,PetscScalar *x,PetscScalar *Ifd)
{
  PetscErrorCode ierr;
  DYNCv      cv;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNEventPostFunction_Cv - Post event routine for this generator model
*/
PetscErrorCode DYNEventPostFunction_Cv(DYNGenModel dyngen, PetscInt nevents, PetscInt ev_list[], PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscBool forward_solve, PetscBool *solve_algebraic)
{
  PetscErrorCode ierr;
  DYNCv      cv;

  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);

  *solve_algebraic = PETSC_FALSE;

  PetscFunctionReturn(0);
}

/*
  DYNGenModelSetEvent_Cv - Sets the event info for the CV generator model

  Input Parameters
. dyngen - the dynamic generator model

  Output Parameters
+ nmons - the number of event monitors
. direction - the event directions
. terminate - flags for termination when the event is located
. eventfcn  - the function describing the event condition
- posteventfcn - An optional function that gets called when the event is located
*/
PetscErrorCode DYNGenModelSetEvent_Cv(DYNGenModel dyngen,PetscInt *nmons,PetscInt *direction,PetscBool *terminate,PetscErrorCode (**eventfcn)(DYNGenModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*), PetscErrorCode (**posteventfcn)(DYNGenModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscBool,PetscBool*), PetscErrorCode (**posteventdgamma)(DYNGenModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscInt))
{

  PetscFunctionBegin;
  *nmons = 0;
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetInitialFieldVoltage_Cv - Returns the DC field voltage at t=t0.

  Input Parameters:
. DYNGenModel - dyngen

  Output Parameters:
. Efd0 - the initial field voltage

  Notes:
   The field voltage is calculated by the generator model when the initial conditions are computed
*/
PetscErrorCode DYNGenModelGetInitialFieldVoltage_Cv(DYNGenModel dyngen,PetscScalar *Efd0)
{
  PetscErrorCode ierr;
  DYNCv      dyncv;
  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&dyncv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetInitialMechanicalPower_Cv - Returns the mechanical power at t=t0.

  Input Parameters:
. DYNGenModel - dyngen

  Output Parameters:
. Pmech0 - the initial mechanical power

  Notes:
   The mechanical power is calculated by the generator model when the initial conditions are computed
*/
PetscErrorCode DYNGenModelGetInitialMechanicalPower_Cv(DYNGenModel dyngen,PetscScalar *Pmech0)
{
  PetscErrorCode ierr;
  DYNCv      dyncv;
  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&dyncv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelReadData_Cv - Parses the line containining the CV data and
  populates it in the Cv data structure

  Input Parameters:
+ dyngen - The generator base model
. line   - the file from the line
. mbase  - the machine base
- sbase  - the system base

  Notes: 
    The type of generator model should be set prior to calling this routine via DYNGenModelSetType()

    The conversion from machine base to system base,if applicable, should be done in this routine
*/
PetscErrorCode DYNGenModelReadData_Cv(DYNGenModel dyngen,char* line,PetscScalar mbase,PetscScalar sbase)
{
  PetscErrorCode ierr;
  DYNCv cv;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetBusnumID_Cv - Returns the bus number and ID associated with this generator model

  Input Parameters:
. dyngen - the dynamic generator model object

  Output Parameters:
+ busnum - the bus number 
- genid  - the generator ID
*/
PetscErrorCode DYNGenModelGetBusnumID_Cv(DYNGenModel dyngen,PetscInt *busnum,char **genid)
{
  PetscErrorCode ierr;
  DYNCv      cv;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);
  *busnum = cv->bus_i;
  *genid  = cv->id;
  PetscFunctionReturn(0);
}

/*
  DYNGenModelDestroy_Cv - Destroys the Cv object data

  Input Parameters
. DYNGenModel - the DYNGenModel object

  Notes:
  Called when DYNGenModelDestroy() is called
*/
PetscErrorCode DYNGenModelDestroy_Cv(DYNGenModel dyngen)
{
  PetscErrorCode ierr;
  DYNCv      cv;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);
  ierr = PetscFree(cv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetNvar_Cv - Returns the number of variables for the Cv model

  Input parameters
. dyngen - the dynamic generator object

  Output Parameters
. nvar - number of variables for this model
*/
PetscErrorCode DYNGenModelGetNvar_Cv(DYNGenModel dyngen,PetscInt *nvar)
{
  PetscFunctionBegin;
  *nvar = DYNCv_nvar;
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetSizeof_Cv - Returns the size of Cv struct

  Input Parameters:
. dyngen - the dynamic generator model
 
  Output Parameters:
. size - size of the Cv object obtained from sizeof()
*/
PetscErrorCode DYNGenModelGetSizeof_Cv(DYNGenModel dyngen,PetscInt *size)
{
  PetscFunctionBegin;
  *size = sizeof(struct _p_DYNCv);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetdFreqdState_Cv - Returns the partial derivative of the generator frequency w.r.t the machine state that governs it

  Input Parameters:
+ dyngen - the DYNGenModel object
. t      - the current time
. xdyn   - the state variables for this generator at time t

  Output Parameters:
+ dfreqdstate   - dfreq_dstate
- stateloc      - location w.r.t. to the bus variables
*/
PetscErrorCode DYNGenModelGetdFreqdState_Cv(DYNGenModel dyngen,PetscReal t, const PetscScalar* xdyn,PetscScalar* dfreqdstate,PetscInt *stateloc)
{
  PetscErrorCode ierr;
  DYNCv cv;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetFrequency_Cv - Returns the frequency of the generator at the current time.

  Input Parameters:
+ dyngen - the DYNGenModel object
. t      - the current time
. xdyn   - the state variables for this generator at time t

  Output Parameters:
. frequency   - the generator frequency in Hz

  Notes the locations to insert the values should be obtained by DYNGenModelGetFirstVariableLocation()
*/
PetscErrorCode DYNGenModelGetFrequency_Cv(DYNGenModel dyngen,PetscReal t,const PetscScalar *x,PetscScalar *frequency)
{
  PetscErrorCode ierr;
  DYNCv cv;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelSetFrequency_Cv - Set the frequency component when initializing the sensitivity variable lambda.

  Input Parameters:
+ dyngen - the DYNGenModel object
. xdyn   - the state variables for this generator at time t
. value  - the value to be assigned to frequency component

  Output Parameters:

  Notes the locations to insert the values should be obtained by DYNGenModelGetFirstVariableLocation()
*/
PetscErrorCode DYNGenModelSetFrequency_Cv(DYNGenModel dyngen,PetscScalar *x,PetscScalar value)
{
  PetscErrorCode ierr;
  DYNCv      cv;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetSpeedDeviation_Cv - Returns the speed deviation and its time derivative (optional) 
                                        of the generator at the current time.

  Input Parameters:
+ dyngen - the DYNGenModel object
. t      - the current time
. xdyn   - the state variables for this generator at time t

  Output Parameters:
+ dw   - the generator speed deviation
- dw_dt - time-derivative of speed deviation

*/
PetscErrorCode DYNGenModelGetSpeedDeviation_Cv(DYNGenModel dyngen,PetscReal t,const PetscScalar *x,PetscScalar *dw,PetscScalar *dw_dt)
{
  PetscErrorCode ierr;
  DYNCv cv;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetSpeedDeviationLocation_Cv - Returns the location of the speed deviation variable for this generator.

  Input Parameters:
. dyngen - the DYNGenModel object

  Output Parameters:
. dwloc   - the generator speed deviation variable location relative to the bus

  Notes the locations to insert the values should be obtained by DYNGenModelGetFirstVariableLocation()
*/
PetscErrorCode DYNGenModelGetSpeedDeviationLocation_Cv(DYNGenModel dyngen,PetscInt *dwloc)
{
  PetscErrorCode ierr;
  DYNCv cv;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DYNGenModelSetInitialFieldVoltageDiff_Cv(DYNGenModel dyngen,PetscScalar PG,PetscScalar QG,PetscScalar VA,PetscScalar VM,PetscScalar *Efd_PG,PetscScalar *Efd_QG,PetscScalar *Efd_VA,PetscScalar *Efd_VM)
{
  DYNCv      cv;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelSetInitialConditions_Cv - Sets the initial conditions (x(t0)) for the CV model

  Input Parameters:
+ dyngen - the DYNGenModel object
. PG     - generator real power output
. QG     - generator reactive power output
. VD     - real component of the complex bus voltage
. VQ     - imaginary component of the complex bus voltage

  Output Parameters:
. x - the initial conditions for the CV model

  Notes the locations to insert the values should be obtained by DYNGenModelGetFirstVariableLocation()
*/
PetscErrorCode DYNGenModelSetInitialConditions_Cv(DYNGenModel dyngen,PetscScalar PG, PetscScalar QG, PetscScalar VD, PetscScalar VQ, PetscScalar *x)
{
  PetscErrorCode ierr;
  DYNCv cv;
  PetscScalar Vm2;
  PetscScalar IGD,IGQ; /* Real and imaginary components of generator current */
  PetscInt    loc;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);
  ierr = DYNGenModelGetFirstVariableLocation(dyngen,&loc);CHKERRQ(ierr);

  Vm2 = VD*VD + VQ*VQ;

  IGD = (VD*PG + VQ*QG)/Vm2;
  IGQ = (VQ*PG - VD*QG)/Vm2;

  /* Variables */
  x[loc] = IGD;
  x[loc+1] = IGQ;

  /* Constants used in equations */
  cv->VDs = VD; cv->VQs = VQ;

  PetscFunctionReturn(0);
}

/*
  DYNGenModelSetInitialConditionsPCv - Sets the differentiation of initial conditions (x(t0)) for the CV model

  Input Parameters:
+ dyngen - the DYNGenModel object
. PG     - generator real power output
. QG     - generator reactive power output
. VA     - angle of the bus voltage
. VM     - magnitude of the complex bus voltage

  Output Parameters:
. x_PG - the derivative of initial conditions to PG
. x_QG - the derivative of initial conditions to QG
. x_VA - the derivative of initial conditions to VA
. x_VM - the derivative of initial conditions to VM

  Notes the locations to insert the values should be obtained by DYNGenModelGetFirstVariableLocation()
*/
PetscErrorCode DYNGenModelSetInitialConditionsP_Cv(DYNGenModel dyngen,PetscScalar PG,PetscScalar QG,PetscScalar VA,PetscScalar VM,PetscInt PGloc,PetscInt QGloc,PetscInt VAloc,PetscInt VMloc,Mat ICp,PetscInt dynlocglob)
{
  DYNCv      cv;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetEquationTypes_Cv - Gets the number and indices of differential and algebraic equations
                                                   
  Input Parameters:
. dyngen - the DYNGenModel object

  Output Parameters:
+ ndiff - number of differential equations
. nalg  - number of algebraic equations
- eqtype - an array of size nvar that has
            eqtype[i] = DIFF_EQ if equation i is a differential
            vartype[i] = ALG_VAR if equation i is an algebraic

  NOTES: The user does not need to create the array vartype
*/
PetscErrorCode DYNGenModelGetEquationTypes_Cv(DYNGenModel dyngen,PetscInt *ndiff,PetscInt *nalg, PetscInt *eqtype)
{
  DYNCv cv;
  PetscInt  nd=0,na=2;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);

  eqtype[0] = eqtype[1] = ALG_EQ;

  *ndiff = nd;
  *nalg  = na;

  PetscFunctionReturn(0);
}

/*
  DYNGenModelDAERHSJacobian_Cv - Computes the RHS Jacobian of the Cv DAE equations

  Input Parameters:
+ dyngen     - the dynamic generator model
. t          - the current time
. VD         - real-part of complex bus voltage
. VQ         - imaginary-part of complex bus voltage
. x          - all variables at the bus where this generator is incident
. dynlocglob - starting location of generator variables in the global vector 
. V_loc      - global location of VD and VQ V_loc[0] = VD_loc V_loc[1] = VQ_loc
- I_loc      - global location of ID and IQ I_loc[0] = ID_loc I_loc[1] = IQ_loc

  Output Parameters:
. J          - the Jacobian matrix

  Notes:
   For each bus, the network equations are ordered as [ID_loc;IQ_loc]

   The network equations are written as I_gen - I_net - I_load = 0. So the
   generator model should set dIG_dxdyn in the Jacobian

   The location for the first variable for the cv model should be obtained by
   DYNGenModelGetFirstVariableLocation(dyngen,&loc). x[loc] = 1st cv variable,
   x[loc+n-1] = nth cv variable
*/ 
PetscErrorCode DYNGenModelDAERHSJacobian_Cv(DYNGenModel dyngen,Mat J,PetscReal t,PetscScalar VD,PetscScalar VQ,PetscScalar *x,PetscInt dynlocglob,PetscInt V_loc[],PetscInt I_loc[])
{
  PetscErrorCode ierr;
  DYNCv      cv;
  PetscInt   row,col;
  PetscScalar val;
  PetscInt    loc;
  
  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);
  ierr = DYNGenModelGetFirstVariableLocation(dyngen,&loc);CHKERRQ(ierr);

  row = dynlocglob;
  col = V_loc[0];
  val = -1.0;
  ierr = SetMatrixValues(J,1,&row,1,&col,&val);CHKERRQ(ierr);
  row = dynlocglob+1;
  col = V_loc[1];
  ierr = SetMatrixValues(J,1,&row,1,&col,&val);CHKERRQ(ierr);

 /* dIGD_dxdyn */
  row = I_loc[0];
  col = dynlocglob;
  val = 1.0;
  ierr = SetMatrixValues(J,1,&row,1,&col,&val);CHKERRQ(ierr);
  row = I_loc[1];
  col = dynlocglob+1;
  ierr = SetMatrixValues(J,1,&row,1,&col,&val);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNGenModelDAERHSJacobianP_Cv - Sets the Jacobian of the Cv equations w.r.t. parameters

  Input Parameters:
+ dyngen - dyngenmodel object
. t      - current time
. x      - state variables for this bus
. dynlocglob - starting location of the variables for this dyngen model
. Valoc  - location of parameter voltage angle (Va)
. Vmloc  - location of parameter voltage magnitude (Vm)
. Pgloc  - location of parameter generator real power output (Pg)
- Qgloc  - location of parameter generator reactive power output (Qg)

  Output Parameters:
. jacP - Matrix of partial derivatives of DYN DAE equations w.r.t. parameters
*/
PetscErrorCode DYNGenModelDAERHSJacobianP_Cv(DYNGenModel dyngen, PetscReal t,const PetscScalar *x,Mat jacP,PetscInt dynlocglob,PetscInt Valoc,PetscInt Vmloc,PetscInt Pgloc,PetscInt Qgloc)
{
  PetscErrorCode ierr;
  DYNCv      cv;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelDAEFWDRHSJacobianP_Cv - Sets the Jacobian of the Cv equations w.r.t. parameters

  Input Parameters:
+ dyngen - dyngenmodel object
. t      - current time
. x      - state variables for this bus
. dynlocglob - starting location of the variables for this dyngen model
. Valoc  - location of parameter voltage angle (Va)
. Vmloc  - location of parameter voltage magnitude (Vm)
. Pgloc  - location of parameter generator real power output (Pg)
- Qgloc  - location of parameter generator reactive power output (Qg)

  Output Parameters:
. jacP - Matrix of partial derivatives of DYN DAE equations w.r.t. parameters
*/
PetscErrorCode DYNGenModelDAEFWDRHSJacobianP_Cv(DYNGenModel dyngen, PetscReal t,const PetscScalar *x,Vec *jacP,PetscInt dynlocglob,PetscInt Valoc,PetscInt Vmloc,PetscInt Pgloc,PetscInt Qgloc)
{
  PetscErrorCode ierr;
  DYNCv      cv;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelDAERHSFunction_Cv - Computes the rhs of the DAE function for the CV model 
                                and also returns the generator currents in network refernce frame

  Input Parameters:
+ dyngen - the dynamic generator model
. t      - the current time
. x      - array of all the variables for the bus on which this generator is incident
. VD     - real-part of the bus voltage
- VQ     - imaginary part of the bus voltage

  Output Parameters:
+ f      - array of rhs of DAE equations for the cv model 
. IGD    - real-part of generator current in network reference frame
. IGQ    - imaginary-part of generator current in network reference frame

  Notes:
   The locations for the variables (and the corresponding locations to insert entries in f)
   for this generator model should be obtained by DYNGenModelGetFirstVariableLocation()
*/
PetscErrorCode DYNGenModelDAERHSFunction_Cv(DYNGenModel dyngen,PetscReal t,PetscScalar VD, PetscScalar VQ, PetscScalar *x,PetscScalar *f,PetscScalar *IGD,PetscScalar *IGQ)
{
  PetscErrorCode ierr;
  DYNCv   cv;
  PetscInt    loc;
  
  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);
  ierr = DYNGenModelGetFirstVariableLocation(dyngen,&loc);CHKERRQ(ierr);

  
  *IGD   = x[loc];
  *IGQ  = x[loc+1];

  /* CV equations */
  f[loc]   = cv->VDs - VD;
  f[loc+1] = cv->VQs - VQ;

  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetCurrent_Cv - Returns the current injection

  Input Parameters:
. dyngen - the DYNGenModel object

  Output Parameters:
+ IGD - Real part of current injection
- IGQ - Reactive part of current injection
*/
PetscErrorCode DYNGenModelGetCurrent_Cv(DYNGenModel dyngen,PetscScalar VD, PetscScalar VQ, PetscScalar *xdyn, PetscScalar *IGD, PetscScalar *IGQ)
{
  DYNCv   cv;
  PetscInt    loc;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);
  ierr = DYNGenModelGetFirstVariableLocation(dyngen,&loc);CHKERRQ(ierr);

  *IGD  = xdyn[loc];
  *IGQ  = xdyn[loc+1];

  PetscFunctionReturn(0);
}

/*
  DYNGenModelSetVoltage_Cv - Sets the voltage for this generator model

  Input Parameters:
. dyngen - the DYNGenModel object

  Output Parameters:
+ VD - Real part of voltage
- VQ - Imaginary part of voltage

Notes: The generator voltage is the same as the bus voltage. The bus voltage is solved for
   at each time-step using the bus mismatch equations. The VD and VQ voltages set here are
   for controlling the internal parameters of some generator models, for e.g. constant voltage
   source model.
*/
PetscErrorCode DYNGenModelSetVoltage_Cv(DYNGenModel dyngen,PetscScalar VD, PetscScalar VQ)
{
  DYNCv   cv;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&cv);CHKERRQ(ierr);

  cv->VDs = VD;
  cv->VQs = VQ;

  PetscFunctionReturn(0);
}



/*
  DYNGenModelCreate_Cv - Class constructor for CV model
*/
PetscErrorCode DYNGenModelCreate_Cv(DYNGenModel dyngen)
{
  DYNCv cv;
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = PetscCalloc1(1,&cv);CHKERRQ(ierr);
  dyngen->data = (void*)cv;

  /* Inherit the ops */
  dyngen->ops.readdata             = DYNGenModelReadData_Cv;
  dyngen->ops.destroy              = DYNGenModelDestroy_Cv;
  dyngen->ops.getnvar              = DYNGenModelGetNvar_Cv;
  dyngen->ops.getsizeof            = DYNGenModelGetSizeof_Cv;
  dyngen->ops.setinitialconditions = DYNGenModelSetInitialConditions_Cv;
  dyngen->ops.setinitialconditionsp = DYNGenModelSetInitialConditionsP_Cv;
  dyngen->ops.setinitialfieldvoltagediff = DYNGenModelSetInitialFieldVoltageDiff_Cv;
  dyngen->ops.daerhsfunction       = DYNGenModelDAERHSFunction_Cv;
  dyngen->ops.daerhsjacobian       = DYNGenModelDAERHSJacobian_Cv;
  dyngen->ops.daerhsjacobianp      = DYNGenModelDAERHSJacobianP_Cv;
  dyngen->ops.daefwdrhsjacobianp   = DYNGenModelDAEFWDRHSJacobianP_Cv;
  dyngen->ops.getequationtypes     = DYNGenModelGetEquationTypes_Cv;
  dyngen->ops.getbusnumid          = DYNGenModelGetBusnumID_Cv;
  dyngen->ops.getfieldcurrent      = DYNGenModelGetFieldCurrent_Cv;
  dyngen->ops.getinitialfieldvoltage = DYNGenModelGetInitialFieldVoltage_Cv;
  dyngen->ops.getinitialmechanicalpower = DYNGenModelGetInitialMechanicalPower_Cv;
  dyngen->ops.getfrequency         = DYNGenModelGetFrequency_Cv;
  dyngen->ops.setfrequency         = DYNGenModelSetFrequency_Cv;
  dyngen->ops.getdfreqdstate       = DYNGenModelGetdFreqdState_Cv;
  dyngen->ops.getspeeddeviation    = DYNGenModelGetSpeedDeviation_Cv;
  dyngen->ops.getspeeddeviationlocation = DYNGenModelGetSpeedDeviationLocation_Cv;
  dyngen->ops.setevent             = DYNGenModelSetEvent_Cv;
  dyngen->ops.getcurrent           = DYNGenModelGetCurrent_Cv;
  dyngen->ops.setvoltage           = DYNGenModelSetVoltage_Cv;

  PetscFunctionReturn(0);
}

