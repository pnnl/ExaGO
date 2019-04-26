#include "dynpvd1.h"
#include <private/dyngenmodelsimpl.h>

/*
  DYNEventMonitor_Pvd1 - Event monitoring routine for this generator model
*/
PetscErrorCode DYNEventMonitor_Pvd1(DYNGenModel dyngen, PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscScalar *fval)
{
  PetscErrorCode ierr;
  DYNPvd1        pvd1;

  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNEventPostDGamma_Pvd1 - Compute derivative of gamma
*/
PetscErrorCode DYNEventPostDGamma_Pvd1(DYNGenModel dyngen, PetscInt nevents, PetscInt ev_list[], PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscScalar *dgdx, PetscScalar *dgdp,PetscInt Valoc, PetscInt Vmloc, PetscInt Pgloc, PetscInt Qgloc)
{
  PetscErrorCode ierr;
  DYNPvd1      pvd1;

  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetFieldCurrent_Pvd1 - Returns the field current Ifd

  Input Parameters:
. dyngen - the DYNGenModel object

  Output Parameters:
. Ifd - Field current Ifd
*/
PetscErrorCode DYNGenModelGetFieldCurrent_Pvd1(DYNGenModel dyngen,PetscScalar *x,PetscScalar *Ifd)
{
  PetscErrorCode ierr;
  DYNPvd1        pvd1;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNEventPostFunction_Pvd1 - Post event routine for this generator model
*/
PetscErrorCode DYNEventPostFunction_Pvd1(DYNGenModel dyngen, PetscInt nevents, PetscInt ev_list[], PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscBool forward_solve, PetscBool *solve_algebraic)
{
  PetscErrorCode ierr;
  DYNPvd1      pvd1;

  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNGenModelSetEvent_Pvd1 - Sets the event info for the PVD1 generator model

  Input Parameters
. dyngen - the dynamic generator model

  Output Parameters
+ nmons - the number of event monitors
. direction - the event directions
. terminate - flags for termination when the event is located
. eventfcn  - the function describing the event condition
- posteventfcn - An optional function that gets called when the event is located
*/
PetscErrorCode DYNGenModelSetEvent_Pvd1(DYNGenModel dyngen,PetscInt *nmons,PetscInt *direction,PetscBool *terminate,PetscErrorCode (**eventfcn)(DYNGenModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*), PetscErrorCode (**posteventfcn)(DYNGenModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscBool,PetscBool*), PetscErrorCode (**posteventdgamma)(DYNGenModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscInt))
{

  PetscFunctionBegin;
  *nmons = 2;
  direction[0] = -1; direction[1] = 1;
  terminate[0] = terminate[1] = PETSC_FALSE;
  *eventfcn = DYNEventMonitor_Pvd1;
  *posteventfcn = DYNEventPostFunction_Pvd1;
  *posteventdgamma = DYNEventPostDGamma_Pvd1;
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetInitialFieldVoltage_Pvd1 - Returns the DC field voltage at t=t0.

  Input Parameters:
. DYNGenModel - dyngen

  Output Parameters:
. Efd0 - the initial field voltage

  Notes:
   The field voltage is calculated by the generator model when the initial conditions are computed
*/
PetscErrorCode DYNGenModelGetInitialFieldVoltage_Pvd1(DYNGenModel dyngen,PetscScalar *Efd0)
{
  PetscErrorCode ierr;
  DYNPvd1      dynpvd1;
  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&dynpvd1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetInitialMechanicalPower_Pvd1 - Returns the mechanical power at t=t0.

  Input Parameters:
. DYNGenModel - dyngen

  Output Parameters:
. Pmech0 - the initial mechanical power

  Notes:
   The mechanical power is calculated by the generator model when the initial conditions are computed
*/
PetscErrorCode DYNGenModelGetInitialMechanicalPower_Pvd1(DYNGenModel dyngen,PetscScalar *Pmech0)
{
  PetscErrorCode ierr;
  DYNPvd1      dynpvd1;
  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&dynpvd1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelReadData_Pvd1 - Parses the line containining the PVD1 data and
  populates it in the Pvd1 data structure

  Input Parameters:
+ dyngen - The generator base model
. line   - the file from the line
. mbase  - the machine base
- sbase  - the system base

  Notes: 
    The type of generator model should be set prior to calling this routine via DYNGenModelSetType()

    The conversion from machine base to system base,if applicable, should be done in this routine
*/
PetscErrorCode DYNGenModelReadData_Pvd1(DYNGenModel dyngen,char* line,PetscScalar mbase,PetscScalar sbase)
{
  PetscErrorCode ierr;
  DYNPvd1 pvd1;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);
  /* Read the PVD1 data */
  /* Data format */
  /* IBUS, ’PVD1’, ID, PQFLAG, XC, QMX, QMN, V0, V1, DQDV, FDBD, DDN, IMAX, VT0, VT1, VT2, VT3, VRFLAG, FT0, FT1, FT2, FT3, FRFLAG, TG, TG, VTMAX, IVPNT1, IVPNT0, QMIN, KHV 
   */
  sscanf(line,"%d,'PVD1',%[^,],%d,%lf,%lf,%lf, \
               %lf,%lf,%lf,%lf,%lf,%lf, \
               %lf,%lf,%lf,%lf,%d, \
               %lf,%lf,%lf,%lf,%d, \
               %lf,%lf,%lf,%lf,%lf,%lf,%lf", 
	 &pvd1->bus_i,pvd1->id,&pvd1->pqflag,&pvd1->xc,&pvd1->qmx,&pvd1->qmn,
	 &pvd1->v0,&pvd1->v1,&pvd1->dqdv,&pvd1->fdbd,&pvd1->ddn,&pvd1->imax,
	 &pvd1->vt0,&pvd1->vt1,&pvd1->vt2,&pvd1->vt3,&pvd1->vrflag,
	 &pvd1->ft0,&pvd1->ft1,&pvd1->ft2,&pvd1->ft3,&pvd1->frflag,
	 &pvd1->tg,&pvd1->tf,&pvd1->vtmax,&pvd1->Ivpnt1,&pvd1->Ivpnt0,&pvd1->qmin,&pvd1->Khv);
   
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetBusnumID_Pvd1 - Returns the bus number and ID associated with this generator model

  Input Parameters:
. dyngen - the dynamic generator model object

  Output Parameters:
+ busnum - the bus number 
- genid  - the generator ID
*/
PetscErrorCode DYNGenModelGetBusnumID_Pvd1(DYNGenModel dyngen,PetscInt *busnum,char **genid)
{
  PetscErrorCode ierr;
  DYNPvd1      pvd1;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);
  *busnum = pvd1->bus_i;
  *genid  = pvd1->id;
  PetscFunctionReturn(0);
}

/*
  DYNGenModelDestroy_Pvd1 - Destroys the Pvd1 object data

  Input Parameters
. DYNGenModel - the DYNGenModel object

  Notes:
  Called when DYNGenModelDestroy() is called
*/
PetscErrorCode DYNGenModelDestroy_Pvd1(DYNGenModel dyngen)
{
  PetscErrorCode ierr;
  DYNPvd1      pvd1;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);
  ierr = PetscFree(pvd1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetNvar_Pvd1 - Returns the number of variables for the Pvd1 model

  Input parameters
. dyngen - the dynamic generator object

  Output Parameters
. nvar - number of variables for this model
*/
PetscErrorCode DYNGenModelGetNvar_Pvd1(DYNGenModel dyngen,PetscInt *nvar)
{
  PetscErrorCode ierr;
  DYNPvd1        pvd1;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);
  if(pvd1->tg < PETSC_SMALL) *nvar = 0;
  else *nvar = 2;
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetSizeof_Pvd1 - Returns the size of Pvd1 struct

  Input Parameters:
. dyngen - the dynamic generator model
 
  Output Parameters:
. size - size of the Pvd1 object obtained from sizeof()
*/
PetscErrorCode DYNGenModelGetSizeof_Pvd1(DYNGenModel dyngen,PetscInt *size)
{
  PetscFunctionBegin;
  *size = sizeof(struct _p_DYNPvd1);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetdFreqdState_Pvd1 - Returns the partial derivative of the generator frequency w.r.t the machine state that governs it

  Input Parameters:
+ dyngen - the DYNGenModel object
. t      - the current time
. xdyn   - the state variables for this generator at time t

  Output Parameters:
+ dfreqdstate   - dfreq_dstate
- stateloc      - location w.r.t. to the bus variables
*/
PetscErrorCode DYNGenModelGetdFreqdState_Pvd1(DYNGenModel dyngen,PetscReal t, const PetscScalar* xdyn,PetscScalar* dfreqdstate,PetscInt *stateloc)
{
  PetscErrorCode ierr;
  DYNPvd1 pvd1;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetFrequency_Pvd1 - Returns the frequency of the generator at the current time.

  Input Parameters:
+ dyngen - the DYNGenModel object
. t      - the current time
. xdyn   - the state variables for this generator at time t

  Output Parameters:
. frequency   - the generator frequency in Hz

  Notes the locations to insert the values should be obtained by DYNGenModelGetFirstVariableLocation()
*/
PetscErrorCode DYNGenModelGetFrequency_Pvd1(DYNGenModel dyngen,PetscReal t,const PetscScalar *x,PetscScalar *frequency)
{
  PetscErrorCode ierr;
  DYNPvd1 pvd1;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNGenModelSetFrequency_Pvd1 - Set the frequency component when initializing the sensitivity variable lambda.

  Input Parameters:
+ dyngen - the DYNGenModel object
. xdyn   - the state variables for this generator at time t
. value  - the value to be assigned to frequency component

  Output Parameters:

  Notes the locations to insert the values should be obtained by DYNGenModelGetFirstVariableLocation()
*/
PetscErrorCode DYNGenModelSetFrequency_Pvd1(DYNGenModel dyngen,PetscScalar *x,PetscScalar value)
{
  PetscErrorCode ierr;
  DYNPvd1      pvd1;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetSpeedDeviation_Pvd1 - Returns the speed deviation and its time derivative (optional) 
                                        of the generator at the current time.

  Input Parameters:
+ dyngen - the DYNGenModel object
. t      - the current time
. xdyn   - the state variables for this generator at time t

  Output Parameters:
+ dw   - the generator speed deviation
- dw_dt - time-derivative of speed deviation

*/
PetscErrorCode DYNGenModelGetSpeedDeviation_Pvd1(DYNGenModel dyngen,PetscReal t,const PetscScalar *x,PetscScalar *dw,PetscScalar *dw_dt)
{
  PetscErrorCode ierr;
  DYNPvd1 pvd1;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetSpeedDeviationLocation_Pvd1 - Returns the location of the speed deviation variable for this generator.

  Input Parameters:
. dyngen - the DYNGenModel object

  Output Parameters:
. dwloc   - the generator speed deviation variable location relative to the bus

  Notes the locations to insert the values should be obtained by DYNGenModelGetFirstVariableLocation()
*/
PetscErrorCode DYNGenModelGetSpeedDeviationLocation_Pvd1(DYNGenModel dyngen,PetscInt *dwloc)
{
  PetscErrorCode ierr;
  DYNPvd1 pvd1;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DYNGenModelSetInitialFieldVoltageDiff_Pvd1(DYNGenModel dyngen,PetscScalar PG,PetscScalar QG,PetscScalar VA,PetscScalar VM,PetscScalar *Efd_PG,PetscScalar *Efd_QG,PetscScalar *Efd_VA,PetscScalar *Efd_VM)
{
  DYNPvd1      pvd1;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNGenModelSetInitialConditions_Pvd1 - Sets the initial conditions (x(t0)) for the PVD1 model

  Input Parameters:
+ dyngen - the DYNGenModel object
. PG     - generator real power output
. QG     - generator reactive power output
. VD     - real component of the complex bus voltage
. VQ     - imaginary component of the complex bus voltage

  Output Parameters:
. x - the initial conditions for the PVD1 model

  Notes the locations to insert the values should be obtained by DYNGenModelGetFirstVariableLocation()
*/
PetscErrorCode DYNGenModelSetInitialConditions_Pvd1(DYNGenModel dyngen,PetscScalar PG, PetscScalar QG, PetscScalar VD, PetscScalar VQ, PetscScalar *x)
{
  PetscErrorCode ierr;
  DYNPvd1        pvd1;
  PetscScalar    It,theta,Vd,Ipcmd,Iqcmd;
  PetscInt       loc;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);
  ierr = DYNGenModelGetFirstVariableLocation(dyngen,&loc);CHKERRQ(ierr);

  theta = atan2(VQ,VD);

  /* Transform from network reference frame to reference frame at an angle -theta 
     so that Vd = Vm and Vq = 0 
  */
  Vd = VD*PetscCosScalar(-theta) - VQ*PetscSinScalar(-theta);

  /* Active and reactive currents */
  Ipcmd = PG/Vd; Iqcmd = -QG/Vd;
  
#if 0
  /* For debugging -- check if IGD and IGQ is correct */
  PetscScalar IGp, IGq,Vm,Vm2;
  PetscScalar IGD1,IGQ1,IGD,IGQ;

  IGp = Ipcmd;
  IGq = Iqcmd;

  /* Transform back to network reference frame */
  IGD = IGp*PetscCosScalar(-theta) + IGq*PetscSinScalar(-theta);
  IGQ = IGp*-PetscSinScalar(-theta) + IGq*PetscCosScalar(-theta);

  Vm  = PetscPowScalar(VD*VD + VQ*VQ,0.5);
  Vm2 = Vm*Vm;
  IGD1 = (VD*PG + VQ*QG)/Vm2;
  IGQ1 = (VQ*PG - VD*QG)/Vm2;
#endif

  It = PetscPowScalar(Ipcmd*Ipcmd+Iqcmd*Iqcmd,0.5);
  pvd1->imax *= It;
  x[loc]   =  Ipcmd;
  x[loc+1] =  Iqcmd;

  pvd1->Pref = PG;
  pvd1->Qref = QG;
  
  PetscFunctionReturn(0);
}

/*
  DYNGenModelSetInitialConditionsPPvd1 - Sets the differentiation of initial conditions (x(t0)) for the PVD1 model

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
PetscErrorCode DYNGenModelSetInitialConditionsP_Pvd1(DYNGenModel dyngen,PetscScalar PG,PetscScalar QG,PetscScalar VA,PetscScalar VM,PetscInt PGloc,PetscInt QGloc,PetscInt VAloc,PetscInt VMloc,Mat ICp,PetscInt dynlocglob)
{
  DYNPvd1      pvd1;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetEquationTypes_Pvd1 - Gets the number and indices of differential and algebraic equations
                                                   
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
PetscErrorCode DYNGenModelGetEquationTypes_Pvd1(DYNGenModel dyngen,PetscInt *ndiff,PetscInt *nalg, PetscInt *eqtype)
{
  DYNPvd1 pvd1;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);
  *ndiff = 2; *nalg = 0;
  eqtype[0] = eqtype[1] = DIFF_EQ;
  if(PetscAbsScalar(pvd1->tg) < PETSC_SMALL) {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Bus %d PVD1 model: Cannot have zero tg time constant");
  }

  PetscFunctionReturn(0);
}

/*
  DYNGenModelDAERHSJacobian_Pvd1 - Computes the RHS Jacobian of the Pvd1 DAE equations

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

   The location for the first variable for the pvd1 model should be obtained by
   DYNGenModelGetFirstVariableLocation(dyngen,&loc). x[loc] = 1st pvd1 variable,
   x[loc+n-1] = nth pvd1 variable
*/ 
PetscErrorCode DYNGenModelDAERHSJacobian_Pvd1(DYNGenModel dyngen,Mat J,PetscReal t,PetscScalar VD,PetscScalar VQ,PetscScalar *x,PetscInt dynlocglob,PetscInt V_loc[],PetscInt I_loc[])
{
  PetscErrorCode ierr;
  DYNPvd1      pvd1;
  
  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNGenModelDAERHSJacobianP_Pvd1 - Sets the Jacobian of the Pvd1 equations w.r.t. parameters

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
PetscErrorCode DYNGenModelDAERHSJacobianP_Pvd1(DYNGenModel dyngen, PetscReal t,const PetscScalar *x,Mat jacP,PetscInt dynlocglob,PetscInt Valoc,PetscInt Vmloc,PetscInt Pgloc,PetscInt Qgloc)
{
  PetscErrorCode ierr;
  DYNPvd1      pvd1;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelDAEFWDRHSJacobianP_Pvd1 - Sets the Jacobian of the Pvd1 equations w.r.t. parameters

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
PetscErrorCode DYNGenModelDAEFWDRHSJacobianP_Pvd1(DYNGenModel dyngen, PetscReal t,const PetscScalar *x,Vec *jacP,PetscInt dynlocglob,PetscInt Valoc,PetscInt Vmloc,PetscInt Pgloc,PetscInt Qgloc)
{
  PetscErrorCode ierr;
  DYNPvd1      pvd1;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNGenModelDAERHSFunction_Pvd1 - Computes the rhs of the DAE function for the PVD1 model 
                                and also returns the generator currents in network refernce frame

  Input Parameters:
+ dyngen - the dynamic generator model
. t      - the current time
. x      - array of all the variables for the bus on which this generator is incident
. VD     - real-part of the bus voltage
- VQ     - imaginary part of the bus voltage

  Output Parameters:
+ f      - array of rhs of DAE equations for the pvd1 model 
. IGD    - real-part of generator current in network reference frame
. IGQ    - imaginary-part of generator current in network reference frame

  Notes:
   The locations for the variables (and the corresponding locations to insert entries in f)
   for this generator model should be obtained by DYNGenModelGetFirstVariableLocation()
*/
PetscErrorCode DYNGenModelDAERHSFunction_Pvd1(DYNGenModel dyngen,PetscReal t,PetscScalar VD, PetscScalar VQ, PetscScalar *x,PetscScalar *f,PetscScalar *IGD,PetscScalar *IGQ)
{
  PetscErrorCode ierr;
  DYNPvd1   pvd1;
  PetscScalar Ip,Iq;
  PetscScalar Ipmax,Iqmax,Iqmin,Ipcmd,Iqcmd,Vd,theta,Vt;
  PetscInt    loc;
  
  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&pvd1);CHKERRQ(ierr);
  ierr = DYNGenModelGetFirstVariableLocation(dyngen,&loc);CHKERRQ(ierr);

  theta = atan2(VQ,VD);
  /* Transform from network reference frame to reference frame at an angle -theta 
     so that Vd = Vm and Vq = 0 
  */
  Vd = VD*PetscCosScalar(-theta) - VQ*PetscSinScalar(-theta);

  Vt = Vd;

  /* Active and reactive currents */
  Ipcmd = pvd1->Pref/PetscMax(0.01,Vt); 
  Iqcmd = pvd1->Qref/PetscMax(0.01,Vt);

  Ip = x[loc];
  Iq = x[loc+1];

  if(pvd1->pqflag == 0) { /* Q priority */
    Iqmax = pvd1->imax;
    Iqmin = -Iqmax;
    Iqcmd = PetscMax(Iqmin,PetscMin(Iqmax,Iqcmd));
    Ipmax = PetscPowScalar(pvd1->imax*pvd1->imax-Iqcmd*Iqcmd,0.5);
    Ipcmd = PetscMax(0,PetscMin(Ipmax,Ipcmd));
  } else { /* P priority */
    Ipmax = pvd1->imax;
    Ipcmd = PetscMax(-Ipmax,PetscMin(Ipmax,Ipcmd));
    Iqmax = PetscPowScalar(pvd1->imax*pvd1->imax-Ipcmd*Ipcmd,0.5);
    Iqmin = -Iqmax;
    Iqcmd = PetscMax(Iqmin,PetscMin(Iqmax,Iqcmd));
  }

  pvd1->fvl = pvd1->fvh = 1.0;
  /* Voltage-based tripping logic */
  if(Vt < pvd1->vt1) {
    if( Vt < pvd1->vt0) pvd1->fvl = 0; 
    else pvd1->fvl = (Vt - pvd1->vt0)/(pvd1->vt1 - pvd1->vt0);
  }

  if(Vt > pvd1->vt2) {
    if(Vt > pvd1->vt3) pvd1->fvh = 0; 
    else pvd1->fvh = (Vt - pvd1->vt3)/(pvd1->vt2 - pvd1->vt3);
  }
  pvd1->fvl = pvd1->fvh = 1.0;
 
  Ipcmd = pvd1->fvl*pvd1->fvh*Ipcmd;
  Iqcmd = pvd1->fvl*pvd1->fvh*Iqcmd;

  f[loc]   = (Ipcmd - Ip)/pvd1->tg;
  f[loc+1] = -(Iqcmd + Iq)/pvd1->tg; 

  /* Transform back to network reference frame */
  *IGD = Ip*PetscCosScalar(-theta)  + Iq*PetscSinScalar(-theta);
  *IGQ = Ip*-PetscSinScalar(-theta) + Iq*PetscCosScalar(-theta);
    
  PetscFunctionReturn(0);
}

/*
  DYNGenModelCreate_Pvd1 - Class constructor for PVD1 model
*/
PetscErrorCode DYNGenModelCreate_Pvd1(DYNGenModel dyngen)
{
  DYNPvd1 pvd1;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&pvd1);CHKERRQ(ierr);
  dyngen->data = (void*)pvd1;

  /* Inherit the ops */
  dyngen->ops.readdata             = DYNGenModelReadData_Pvd1;
  dyngen->ops.destroy              = DYNGenModelDestroy_Pvd1;
  dyngen->ops.getnvar              = DYNGenModelGetNvar_Pvd1;
  dyngen->ops.getsizeof            = DYNGenModelGetSizeof_Pvd1;
  dyngen->ops.setinitialconditions = DYNGenModelSetInitialConditions_Pvd1;
  dyngen->ops.setinitialconditionsp = DYNGenModelSetInitialConditionsP_Pvd1;
  dyngen->ops.setinitialfieldvoltagediff = DYNGenModelSetInitialFieldVoltageDiff_Pvd1;
  dyngen->ops.daerhsfunction       = DYNGenModelDAERHSFunction_Pvd1;
  dyngen->ops.daerhsjacobian       = DYNGenModelDAERHSJacobian_Pvd1;
  dyngen->ops.daerhsjacobianp      = DYNGenModelDAERHSJacobianP_Pvd1;
  dyngen->ops.daefwdrhsjacobianp   = DYNGenModelDAEFWDRHSJacobianP_Pvd1;
  dyngen->ops.getequationtypes     = DYNGenModelGetEquationTypes_Pvd1;
  dyngen->ops.getbusnumid          = DYNGenModelGetBusnumID_Pvd1;
  dyngen->ops.getfieldcurrent      = DYNGenModelGetFieldCurrent_Pvd1;
  dyngen->ops.getinitialfieldvoltage = DYNGenModelGetInitialFieldVoltage_Pvd1;
  dyngen->ops.getinitialmechanicalpower = DYNGenModelGetInitialMechanicalPower_Pvd1;
  dyngen->ops.getfrequency         = DYNGenModelGetFrequency_Pvd1;
  dyngen->ops.setfrequency         = DYNGenModelSetFrequency_Pvd1;
  dyngen->ops.getdfreqdstate       = DYNGenModelGetdFreqdState_Pvd1;
  dyngen->ops.getspeeddeviation    = DYNGenModelGetSpeedDeviation_Pvd1;
  dyngen->ops.getspeeddeviationlocation = DYNGenModelGetSpeedDeviationLocation_Pvd1;
  //  dyngen->ops.setevent             = DYNGenModelSetEvent_Pvd1;

  PetscFunctionReturn(0);
}

