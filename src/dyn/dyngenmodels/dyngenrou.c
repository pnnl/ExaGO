#include "dyngenrou.h"
#include <private/dyngenmodelsimpl.h>

/*
  DYNEventMonitor_Genrou - Event monitoring routine for this generator model
*/
PetscErrorCode DYNEventMonitor_Genrou(DYNGenModel dyngen, PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscScalar *fval)
{
  PetscErrorCode ierr;
  DYNGenrou      genrou;
  PetscScalar    f,fmin,fmax;

  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNGenModelGetModelData(dyngen,(void**)&genrou);CHKERRQ(ierr);

  ierr = DYNGenModelGetFrequency(dyngen,t,x,&f);CHKERRQ(ierr);
  ierr = DYNGenModelGetFrequencyLimits(dyngen,&fmax,&fmin);CHKERRQ(ierr);

  if(!genrou->freqatmax) {
    fval[0] = fmax - f;
  } else {
    fval[0] = 1; /* deactivated frequency limit check */
  }

  if(!genrou->freqatmin) fval[1] = fmin - f;
  else fval[1] = -1; /* deactivated frequency limit check */

  PetscFunctionReturn(0);
}

/*
  DYNEventPostDGamma_Genrou - Compute derivative of gamma
*/
PetscErrorCode DYNEventPostDGamma_Genrou(DYNGenModel dyngen, PetscInt nevents, PetscInt ev_list[], PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscScalar *dgdx, PetscScalar *dgdp,PetscInt Valoc, PetscInt Vmloc, PetscInt Pgloc, PetscInt Qgloc)
{
  PetscErrorCode ierr;
  DYNGenrou      genrou;
  PetscInt       loc;

  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNGenModelGetModelData(dyngen,(void**)&genrou);CHKERRQ(ierr);

  /* Get the location for the first exciter variable */
  ierr = DYNGenModelGetFirstVariableLocation(dyngen,&loc);CHKERRQ(ierr);

  /* at one time, only one event is triggered */
  if(ev_list[0] == 0) { /* Max freq. */
  }else{
  }
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetFieldCurrent_Genrou - Returns the field current Ifd

  Input Parameters:
. dyngen - the DYNGenModel object

  Output Parameters:
. Ifd - Field current Ifd
*/
PetscErrorCode DYNGenModelGetFieldCurrent_Genrou(DYNGenModel dyngen,PetscScalar *x,PetscScalar *Ifd)
{
  PetscErrorCode ierr;
  DYNGenrou      genrou;
  PetscInt       loc;
  PetscScalar    Eqp,psi1d,Id;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&genrou);CHKERRQ(ierr);

  ierr = DYNGenModelGetFirstVariableLocation(dyngen,&loc);CHKERRQ(ierr);

  Eqp   = x[loc];
  psi1d = x[loc+2];
  Id    = x[loc+6];

  *Ifd   = (Eqp + (genrou->Xd - genrou->Xdp)*(Id - (genrou->Xdp - genrou->Xddp)/PetscPowScalar((genrou->Xdp - genrou->Xl),2)*(psi1d + (genrou->Xdp - genrou->Xl)*Id - Eqp)))/(genrou->Xd - genrou->Xl);

  PetscFunctionReturn(0);
}

/*
  DYNEventPostFunction_Genrou - Post event routine for this generator model
*/
PetscErrorCode DYNEventPostFunction_Genrou(DYNGenModel dyngen, PetscInt nevents, PetscInt ev_list[], PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscBool forward_solve, PetscBool *solve_algebraic)
{
  PetscErrorCode ierr;
  DYNGenrou      genrou;
  PetscScalar    f,fmin,fmax;

  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNGenModelGetModelData(dyngen,(void**)&genrou);CHKERRQ(ierr);

  *solve_algebraic = PETSC_TRUE;
  ierr = DYNGenModelSetStatus(dyngen,0);CHKERRQ(ierr);

  ierr = DYNGenModelGetFrequencyLimits(dyngen,&fmax,&fmin);CHKERRQ(ierr);
  ierr = DYNGenModelGetFrequency(dyngen,t,x,&f);CHKERRQ(ierr);

  if(ev_list[0] == 0) { /* Max. freq */
    genrou->freqatmax = 1;
    ierr = PetscPrintf(PETSC_COMM_SELF,"[%4.3f s] Bus %d: Generator %s tripping by over-frequency protection fmax = %4.3f Hz\n",t,genrou->bus_i,genrou->id,fmax);CHKERRQ(ierr);
  } else { /* Min. freq. */
    genrou->freqatmin = 1;
    ierr = PetscPrintf(PETSC_COMM_SELF,"[%4.3f s] Bus %d: Generator %s tripping by under-frequency protection fmin = %4.3f Hz\n",t,genrou->bus_i,genrou->id,fmin);CHKERRQ(ierr);

  }

  PetscFunctionReturn(0);
}

/*
  DYNGenModelSetEvent_Genrou - Sets the event info for the GENROU generator model

  Input Parameters
. dyngen - the dynamic generator model

  Output Parameters
+ nmons - the number of event monitors
. direction - the event directions
. terminate - flags for termination when the event is located
. eventfcn  - the function describing the event condition
- posteventfcn - An optional function that gets called when the event is located
*/
PetscErrorCode DYNGenModelSetEvent_Genrou(DYNGenModel dyngen,PetscInt *nmons,PetscInt *direction,PetscBool *terminate,PetscErrorCode (**eventfcn)(DYNGenModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*), PetscErrorCode (**posteventfcn)(DYNGenModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscBool,PetscBool*), PetscErrorCode (**posteventdgamma)(DYNGenModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscInt))
{

  PetscFunctionBegin;
  *nmons = 2;
  direction[0] = -1; direction[1] = 1;
  terminate[0] = terminate[1] = PETSC_FALSE;
  *eventfcn = DYNEventMonitor_Genrou;
  *posteventfcn = DYNEventPostFunction_Genrou;
  *posteventdgamma = DYNEventPostDGamma_Genrou;
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetInitialFieldVoltage_Genrou - Returns the DC field voltage at t=t0.

  Input Parameters:
. DYNGenModel - dyngen

  Output Parameters:
. Efd0 - the initial field voltage

  Notes:
   The field voltage is calculated by the generator model when the initial conditions are computed
*/
PetscErrorCode DYNGenModelGetInitialFieldVoltage_Genrou(DYNGenModel dyngen,PetscScalar *Efd0)
{
  PetscErrorCode ierr;
  DYNGenrou      dyngenrou;
  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&dyngenrou);CHKERRQ(ierr);
  *Efd0 = dyngenrou->Efd;
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetInitialMechanicalPower_Genrou - Returns the mechanical power at t=t0.

  Input Parameters:
. DYNGenModel - dyngen

  Output Parameters:
. Pmech0 - the initial mechanical power

  Notes:
   The mechanical power is calculated by the generator model when the initial conditions are computed
*/
PetscErrorCode DYNGenModelGetInitialMechanicalPower_Genrou(DYNGenModel dyngen,PetscScalar *Pmech0)
{
  PetscErrorCode ierr;
  DYNGenrou      dyngenrou;
  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&dyngenrou);CHKERRQ(ierr);
  *Pmech0 = dyngenrou->Pm;
  PetscFunctionReturn(0);
}

/*
  DYNGenModelReadData_Genrou - Parses the line containining the GENROU data and
  populates it in the Genrou data structure

  Input Parameters:
+ dyngen - The generator base model
. line   - the file from the line
. mbase  - the machine base
- sbase  - the system base

  Notes: 
    The type of generator model should be set prior to calling this routine via DYNGenModelSetType()

    The conversion from machine base to system base,if applicable, should be done in this routine
*/
PetscErrorCode DYNGenModelReadData_Genrou(DYNGenModel dyngen,char* line,PetscScalar mbase,PetscScalar sbase)
{
  PetscErrorCode ierr;
  DYNGenrou genrou;
  PetscScalar cratio=mbase/sbase;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&genrou);CHKERRQ(ierr);
  /* Read the GENROU data */
  /* Data format */
  /* IBUS, ’GENROU’, I, T’do, T"do, T’qo, T"qo, H, D, Xd, Xq, X’d, X’q, X"d, Xl, S(1.0), S(1.2)/ */
  sscanf(line,"%d,'GENROU',%[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",&genrou->bus_i,genrou->id,&genrou->Td0p,&genrou->Td0dp,&genrou->Tq0p,&genrou->Tq0dp,&genrou->H,&genrou->D,&genrou->Xd,&genrou->Xq,&genrou->Xdp,&genrou->Xqp,&genrou->Xddp,&genrou->Xl,&genrou->S1,&genrou->S2);

  /* Convert from machine base to system base */
  genrou->H *= cratio;
  genrou->D *= cratio;

  genrou->Xd   /= cratio;
  genrou->Xq   /= cratio;
  genrou->Xdp  /= cratio;
  genrou->Xqp  /= cratio;
  genrou->Xddp /= cratio;
  genrou->Xl   /= cratio;

  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetBusnumID_Genrou - Returns the bus number and ID associated with this generator model

  Input Parameters:
. dyngen - the dynamic generator model object

  Output Parameters:
+ busnum - the bus number 
- genid  - the generator ID
*/
PetscErrorCode DYNGenModelGetBusnumID_Genrou(DYNGenModel dyngen,PetscInt *busnum,char **genid)
{
  PetscErrorCode ierr;
  DYNGenrou      genrou;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&genrou);CHKERRQ(ierr);
  *busnum = genrou->bus_i;
  *genid  = genrou->id;
  PetscFunctionReturn(0);
}

/*
  DYNGenModelDestroy_Genrou - Destroys the Genrou object data

  Input Parameters
. DYNGenModel - the DYNGenModel object

  Notes:
  Called when DYNGenModelDestroy() is called
*/
PetscErrorCode DYNGenModelDestroy_Genrou(DYNGenModel dyngen)
{
  PetscErrorCode ierr;
  DYNGenrou      genrou;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&genrou);CHKERRQ(ierr);
  ierr = PetscFree(genrou);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetNvar_Genrou - Returns the number of variables for the Genrou model

  Input parameters
. dyngen - the dynamic generator object

  Output Parameters
. nvar - number of variables for this model
*/
PetscErrorCode DYNGenModelGetNvar_Genrou(DYNGenModel dyngen,PetscInt *nvar)
{
  PetscFunctionBegin;
  *nvar = DYNGenrou_nvar;
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetSizeof_Genrou - Returns the size of Genrou struct

  Input Parameters:
. dyngen - the dynamic generator model
 
  Output Parameters:
. size - size of the Genrou object obtained from sizeof()
*/
PetscErrorCode DYNGenModelGetSizeof_Genrou(DYNGenModel dyngen,PetscInt *size)
{
  PetscFunctionBegin;
  *size = sizeof(struct _p_DYNGenrou);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetdFreqdState_Genrou - Returns the partial derivative of the generator frequency w.r.t the machine state that governs it

  Input Parameters:
+ dyngen - the DYNGenModel object
. t      - the current time
. xdyn   - the state variables for this generator at time t

  Output Parameters:
+ dfreqdstate   - dfreq_dstate
- stateloc      - location w.r.t. to the bus variables
*/
PetscErrorCode DYNGenModelGetdFreqdState_Genrou(DYNGenModel dyngen,PetscReal t, const PetscScalar* xdyn,PetscScalar* dfreqdstate,PetscInt *stateloc)
{
  PetscErrorCode ierr;
  DYNGenrou genrou;
  PetscInt    loc;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&genrou);CHKERRQ(ierr);

  ierr = DYNGenModelGetFirstVariableLocation(dyngen,&loc);CHKERRQ(ierr);

  *dfreqdstate = freq;
  if(stateloc) *stateloc = loc+5;

  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetFrequency_Genrou - Returns the frequency of the generator at the current time.

  Input Parameters:
+ dyngen - the DYNGenModel object
. t      - the current time
. xdyn   - the state variables for this generator at time t

  Output Parameters:
. frequency   - the generator frequency in Hz

  Notes the locations to insert the values should be obtained by DYNGenModelGetFirstVariableLocation()
*/
PetscErrorCode DYNGenModelGetFrequency_Genrou(DYNGenModel dyngen,PetscReal t,const PetscScalar *x,PetscScalar *frequency)
{
  PetscErrorCode ierr;
  DYNGenrou genrou;
  PetscScalar dw;
  PetscInt    loc;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&genrou);CHKERRQ(ierr);

  ierr = DYNGenModelGetFirstVariableLocation(dyngen,&loc);CHKERRQ(ierr);

  dw = x[loc+5];

  *frequency = freq*(1 + dw);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelSetFrequency_Genrou - Set the frequency component when initializing the sensitivity variable lambda.

  Input Parameters:
+ dyngen - the DYNGenModel object
. xdyn   - the state variables for this generator at time t
. value  - the value to be assigned to frequency component

  Output Parameters:

  Notes the locations to insert the values should be obtained by DYNGenModelGetFirstVariableLocation()
*/
PetscErrorCode DYNGenModelSetFrequency_Genrou(DYNGenModel dyngen,PetscScalar *x,PetscScalar value)
{
  PetscErrorCode ierr;
  DYNGenrou      genrou;
  PetscInt       loc;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&genrou);CHKERRQ(ierr);
  ierr = DYNGenModelGetFirstVariableLocation(dyngen,&loc);CHKERRQ(ierr);
  x[loc+5] = value;
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetSpeedDeviation_Genrou - Returns the speed deviation and its time derivative (optional) 
                                        of the generator at the current time.

  Input Parameters:
+ dyngen - the DYNGenModel object
. t      - the current time
. xdyn   - the state variables for this generator at time t

  Output Parameters:
+ dw   - the generator speed deviation
- dw_dt - time-derivative of speed deviation

*/
PetscErrorCode DYNGenModelGetSpeedDeviation_Genrou(DYNGenModel dyngen,PetscReal t,const PetscScalar *x,PetscScalar *dw,PetscScalar *dw_dt)
{
  PetscErrorCode ierr;
  DYNGenrou genrou;
  PetscInt    loc;
  DYNTurbgovModel dynturbgov;
  PetscScalar     Pmech,Eqp,Edp,Id,Iq,psi1d,psi2q;
  PetscScalar tempd1,tempd2,tempq1,tempq2;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&genrou);CHKERRQ(ierr);

  ierr = DYNGenModelGetFirstVariableLocation(dyngen,&loc);CHKERRQ(ierr);

  *dw = x[loc+5];

  if(dw_dt) {
    ierr = DYNGenModelGetDynTurbgov(dyngen,&dynturbgov);CHKERRQ(ierr);
    if(!dynturbgov) {
      Pmech = genrou->Pm;
    } else {
      ierr = DYNTurbgovModelGetMechanicalPower(dynturbgov,t,(PetscScalar*)x,&Pmech,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    }

    Eqp   = x[loc];
    Edp   = x[loc+1];
    psi1d = x[loc+2];
    psi2q = x[loc+3];
    Id    = x[loc+6];
    Iq    = x[loc+7];

    tempd1 = (genrou->Xddp - genrou->Xl)/(genrou->Xdp - genrou->Xl);
    tempd2 = (genrou->Xdp - genrou->Xddp)/(genrou->Xdp - genrou->Xl);
    tempq1 = (genrou->Xddp - genrou->Xl)/(genrou->Xqp - genrou->Xl);
    tempq2 = (genrou->Xqp - genrou->Xddp)/(genrou->Xqp - genrou->Xl);
    *dw_dt = ((Pmech-genrou->D*x[loc+5])/(1+x[loc+5]) - tempd1*Eqp*Iq - tempd2*psi1d*Iq - tempq1*Edp*Id + tempq2*psi2q*Id - (genrou->Xddp - genrou->Xddp)*Id*Iq)/(2*genrou->H);
  }
    

  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetSpeedDeviationLocation_Genrou - Returns the location of the speed deviation variable for this generator.

  Input Parameters:
. dyngen - the DYNGenModel object

  Output Parameters:
. dwloc   - the generator speed deviation variable location relative to the bus

  Notes the locations to insert the values should be obtained by DYNGenModelGetFirstVariableLocation()
*/
PetscErrorCode DYNGenModelGetSpeedDeviationLocation_Genrou(DYNGenModel dyngen,PetscInt *dwloc)
{
  PetscErrorCode ierr;
  DYNGenrou genrou;
  PetscInt    loc;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&genrou);CHKERRQ(ierr);

  ierr = DYNGenModelGetFirstVariableLocation(dyngen,&loc);CHKERRQ(ierr);

  *dwloc = loc+5;

  PetscFunctionReturn(0);
}

PetscErrorCode DYNGenModelSetInitialFieldVoltageDiff_Genrou(DYNGenModel dyngen,PetscScalar PG,PetscScalar QG,PetscScalar VA,PetscScalar VM,PetscScalar *Efd_PG,PetscScalar *Efd_QG,PetscScalar *Efd_VA,PetscScalar *Efd_VM)
{
  DYNGenrou      genrou;
  PetscScalar    VD,VQ,Vm2;
  PetscScalar    IGD,IGQ; /* Real and imaginary components of generator current */
  PetscScalar    delta,theta;
  PetscScalar    Vm2_VD,Vm2_VQ,IGD_PG,IGD_QG,IGD_VD,IGD_VQ,IGQ_PG,IGQ_QG,IGQ_VD,IGQ_VQ,delta_PG,delta_QG,delta_VD,delta_VQ,Id_PG,Id_QG,Id_VD,Id_VQ,Vq_PG,Vq_QG,Vq_VD,Vq_VQ,Eqp_PG,Eqp_QG,Eqp_VD,Eqp_VQ,delta_1,delta_2,VD_VA,VD_VM,VQ_VA,VQ_VM;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&genrou);CHKERRQ(ierr);

  VD = VM*PetscCosScalar(VA);
  VQ = VM*PetscSinScalar(VA);
  Vm2 = VD*VD + VQ*VQ;
  Vm2_VD = 2*VD;
  Vm2_VQ = 2*VQ;

  IGD = (VD*PG + VQ*QG)/Vm2;
  IGQ = (VQ*PG - VD*QG)/Vm2;
  IGD_PG = VD/Vm2;
  IGD_QG = VQ/Vm2;
  IGD_VD = PG/Vm2 - VD*PG*Vm2_VD/(Vm2*Vm2);
  IGD_VQ = QG/Vm2 - VQ*QG*Vm2_VQ/(Vm2*Vm2);
  IGQ_PG = VQ/Vm2;
  IGQ_QG = -VD/Vm2;
  IGQ_VD = -QG/Vm2 + VD*QG*Vm2_VD/(Vm2*Vm2);
  IGQ_VQ = PG/Vm2 - VQ*PG*Vm2_VQ/(Vm2*Vm2);

  delta = atan2(VQ + genrou->Xq*IGD,VD-genrou->Xq*IGQ);
  delta_1 = (VD - genrou->Xq*IGQ)/((VQ + genrou->Xq*IGD)*(VQ + genrou->Xq*IGD) + (VD-genrou->Xq*IGQ)*(VD-genrou->Xq*IGQ));
  delta_2 = -(VQ + genrou->Xq*IGD)/((VQ + genrou->Xq*IGD)*(VQ + genrou->Xq*IGD) + (VD-genrou->Xq*IGQ)*(VD-genrou->Xq*IGQ));
  delta_PG = delta_1*genrou->Xq*IGD_PG - delta_2*genrou->Xq*IGQ_PG;
  delta_QG = delta_1*genrou->Xq*IGD_QG - delta_2*genrou->Xq*IGQ_QG;
  delta_VD = delta_1*genrou->Xq*IGD_VD + delta_2*(1.-genrou->Xq*IGQ_VD);
  delta_VQ = delta_1*(1.+genrou->Xq*IGD_VQ) - delta_2*genrou->Xq*IGQ_VQ;

  theta = PETSC_PI/2.0 - delta;

  //Id = IGD*PetscCosScalar(theta) - IGQ*PetscSinScalar(theta);
   //Iq = IGD*PetscSinScalar(theta) + IGQ*PetscCosScalar(theta);
  Id_PG = IGD_PG*PetscCosScalar(theta) + IGD*PetscSinScalar(theta)*delta_PG - IGQ_PG*PetscSinScalar(theta) + IGQ*PetscCosScalar(theta)*delta_PG;
  Id_QG = IGD_QG*PetscCosScalar(theta) + IGD*PetscSinScalar(theta)*delta_QG - IGQ_QG*PetscSinScalar(theta) + IGQ*PetscCosScalar(theta)*delta_QG;
  Id_VD = IGD_VD*PetscCosScalar(theta) + IGD*PetscSinScalar(theta)*delta_VD - IGQ_VD*PetscSinScalar(theta) + IGQ*PetscCosScalar(theta)*delta_VD;
  Id_VQ = IGD_VQ*PetscCosScalar(theta) + IGD*PetscSinScalar(theta)*delta_VQ - IGQ_VQ*PetscSinScalar(theta) + IGQ*PetscCosScalar(theta)*delta_VQ;

  Vq_PG = - VD*PetscCosScalar(theta)*delta_PG + VQ*PetscSinScalar(theta)*delta_PG;
  Vq_QG = - VD*PetscCosScalar(theta)*delta_QG + VQ*PetscSinScalar(theta)*delta_QG;
  Vq_VD = PetscSinScalar(theta) - VD*PetscCosScalar(theta)*delta_VD + VQ*PetscSinScalar(theta)*delta_VD;
  Vq_VQ = - VD*PetscCosScalar(theta)*delta_VQ + PetscCosScalar(theta) + VQ*PetscSinScalar(theta)*delta_VQ;

  Eqp_PG = Vq_PG + genrou->Xdp*Id_PG;
  Eqp_QG = Vq_QG + genrou->Xdp*Id_QG;
  Eqp_VD = Vq_VD + genrou->Xdp*Id_VD;
  Eqp_VQ = Vq_VQ + genrou->Xdp*Id_VQ;

  //genrou->Efd = Eqp + (genrou->Xd - genrou->Xdp)*Id;
  *Efd_PG = Eqp_PG + (genrou->Xd-genrou->Xdp)*Id_PG;
  *Efd_QG = Eqp_QG + (genrou->Xd-genrou->Xdp)*Id_QG;
  VD_VA = -VM*PetscSinScalar(VA);
  VD_VM = PetscCosScalar(VA);
  VQ_VA = VM*PetscCosScalar(VA);
  VQ_VM = PetscSinScalar(VA);
  *Efd_VA = Eqp_VD*VD_VA + Eqp_VQ*VQ_VA + (genrou->Xd-genrou->Xdp)*(Id_VD*VD_VA + Id_VQ*VQ_VA);
  *Efd_VM = Eqp_VD*VD_VM + Eqp_VQ*VQ_VM + (genrou->Xd-genrou->Xdp)*(Id_VD*VD_VM + Id_VQ*VQ_VM);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelSetInitialConditions_Genrou - Sets the initial conditions (x(t0)) for the GENROU model

  Input Parameters:
+ dyngen - the DYNGenModel object
. PG     - generator real power output
. QG     - generator reactive power output
. VD     - real component of the complex bus voltage
. VQ     - imaginary component of the complex bus voltage

  Output Parameters:
. x - the initial conditions for the GENROU model

  Notes the locations to insert the values should be obtained by DYNGenModelGetFirstVariableLocation()
*/
PetscErrorCode DYNGenModelSetInitialConditions_Genrou(DYNGenModel dyngen,PetscScalar PG, PetscScalar QG, PetscScalar VD, PetscScalar VQ, PetscScalar *x)
{
  PetscErrorCode ierr;
  DYNGenrou genrou;
  PetscScalar Vm2;
  PetscScalar IGD,IGQ; /* Real and imaginary components of generator current */
  PetscScalar delta,theta;
  PetscScalar Id, Iq; /* genertor currents in synchronously rotating dq axis reference frame */
  PetscScalar Vd, Vq; /* terminal voltage in dq axis reference frame */
  PetscScalar Edp,Eqp; /* d and q axis transient emfs */
  PetscScalar psi1d,psi2q; /* d and q axis fluxes */
  PetscInt    loc;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&genrou);CHKERRQ(ierr);
  ierr = DYNGenModelGetFirstVariableLocation(dyngen,&loc);CHKERRQ(ierr);

  Vm2 = VD*VD + VQ*VQ;

  IGD = (VD*PG + VQ*QG)/Vm2;
  IGQ = (VQ*PG - VD*QG)/Vm2;

  delta = atan2(VQ + genrou->Xq*IGD,VD-genrou->Xq*IGQ);
  theta = PETSC_PI/2.0 - delta;

  Id = IGD*PetscCosScalar(theta) - IGQ*PetscSinScalar(theta);
  Iq = IGD*PetscSinScalar(theta) + IGQ*PetscCosScalar(theta);

  Vd = VD*PetscCosScalar(theta) - VQ*PetscSinScalar(theta);
  Vq = VD*PetscSinScalar(theta) + VQ*PetscCosScalar(theta);

  Edp = Vd - genrou->Xqp*Iq;
  Eqp = Vq + genrou->Xdp*Id;

  psi1d = Eqp  - (genrou->Xdp - genrou->Xl)*Id;
  psi2q = -Edp - (genrou->Xqp - genrou->Xl)*Iq;

  /* Variables */
  x[loc] = Eqp;
  x[loc+1] = Edp;
  x[loc+2] = psi1d;
  x[loc+3] = psi2q;
  x[loc+4] = delta;
  x[loc+5] = 0.0;  /* machine speed deviation \delta\omega */
  x[loc+6] = Id;
  x[loc+7] = Iq;

  /* Constants used in derivative evaluation */
  genrou->Efd = Eqp + (genrou->Xd - genrou->Xdp)*Id;
  genrou->Pm  = PG;
  PetscFunctionReturn(0);
}

/*
  DYNGenModelSetInitialConditionsPGenrou - Sets the differentiation of initial conditions (x(t0)) for the GENROU model

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
PetscErrorCode DYNGenModelSetInitialConditionsP_Genrou(DYNGenModel dyngen,PetscScalar PG,PetscScalar QG,PetscScalar VA,PetscScalar VM,PetscInt PGloc,PetscInt QGloc,PetscInt VAloc,PetscInt VMloc,Mat ICp,PetscInt dynlocglob)
{
  DYNGenrou      genrou;
  PetscScalar    VD,VQ,Vm2;
  PetscScalar    IGD,IGQ; /* Real and imaginary components of generator current */
  PetscScalar    delta,theta;
  PetscScalar    x_PG[8],x_QG[8],x_VA[8],x_VM[8];
  PetscScalar    Vm2_VD,Vm2_VQ,IGD_PG,IGD_QG,IGD_VD,IGD_VQ,IGQ_PG,IGQ_QG,IGQ_VD,IGQ_VQ,delta_PG,delta_QG,delta_VD,delta_VQ,Id_PG,Id_QG,Id_VD,Id_VQ,Iq_PG,Iq_QG,Iq_VD,Iq_VQ,Vd_PG,Vd_QG,Vd_VD,Vd_VQ,Vq_PG,Vq_QG,Vq_VD,Vq_VQ,Edp_PG,Edp_QG,Edp_VD,Edp_VQ,Eqp_PG,Eqp_QG,Eqp_VD,Eqp_VQ,delta_1,delta_2,VD_VA,VD_VM,VQ_VA,VQ_VM;
  PetscInt       i,rows[8];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&genrou);CHKERRQ(ierr);

  VD = VM*PetscCosScalar(VA);
  VQ = VM*PetscSinScalar(VA);
  Vm2 = VD*VD + VQ*VQ;
  Vm2_VD = 2.*VD;
  Vm2_VQ = 2.*VQ;

  IGD = (VD*PG + VQ*QG)/Vm2;
  IGQ = (VQ*PG - VD*QG)/Vm2;
  IGD_PG = VD/Vm2;
  IGD_QG = VQ/Vm2;
  IGD_VD = PG/Vm2 - VD*PG*Vm2_VD/(Vm2*Vm2);
  IGD_VQ = QG/Vm2 - VQ*QG*Vm2_VQ/(Vm2*Vm2);
  IGQ_PG = VQ/Vm2;
  IGQ_QG = -VD/Vm2;
  IGQ_VD = -QG/Vm2 + VD*QG*Vm2_VD/(Vm2*Vm2);
  IGQ_VQ = PG/Vm2 - VQ*PG*Vm2_VQ/(Vm2*Vm2);

  delta = atan2(VQ + genrou->Xq*IGD,VD-genrou->Xq*IGQ);
  delta_1 = (VD - genrou->Xq*IGQ)/((VQ + genrou->Xq*IGD)*(VQ + genrou->Xq*IGD) + (VD-genrou->Xq*IGQ)*(VD-genrou->Xq*IGQ));
  delta_2 = -(VQ + genrou->Xq*IGD)/((VQ + genrou->Xq*IGD)*(VQ + genrou->Xq*IGD) + (VD-genrou->Xq*IGQ)*(VD-genrou->Xq*IGQ));
  delta_PG = delta_1*genrou->Xq*IGD_PG - delta_2*genrou->Xq*IGQ_PG;
  delta_QG = delta_1*genrou->Xq*IGD_QG - delta_2*genrou->Xq*IGQ_QG;
  delta_VD = delta_1*genrou->Xq*IGD_VD + delta_2*(1.-genrou->Xq*IGQ_VD);
  delta_VQ = delta_1*(1.+genrou->Xq*IGD_VQ) - delta_2*genrou->Xq*IGQ_VQ;

  theta = PETSC_PI/2.0 - delta;

  //Id = IGD*PetscCosScalar(theta) - IGQ*PetscSinScalar(theta);
  //Iq = IGD*PetscSinScalar(theta) + IGQ*PetscCosScalar(theta);
  Id_PG = IGD_PG*PetscCosScalar(theta) + IGD*PetscSinScalar(theta)*delta_PG - IGQ_PG*PetscSinScalar(theta) + IGQ*PetscCosScalar(theta)*delta_PG;
  Id_QG = IGD_QG*PetscCosScalar(theta) + IGD*PetscSinScalar(theta)*delta_QG - IGQ_QG*PetscSinScalar(theta) + IGQ*PetscCosScalar(theta)*delta_QG;
  Id_VD = IGD_VD*PetscCosScalar(theta) + IGD*PetscSinScalar(theta)*delta_VD - IGQ_VD*PetscSinScalar(theta) + IGQ*PetscCosScalar(theta)*delta_VD;
  Id_VQ = IGD_VQ*PetscCosScalar(theta) + IGD*PetscSinScalar(theta)*delta_VQ - IGQ_VQ*PetscSinScalar(theta) + IGQ*PetscCosScalar(theta)*delta_VQ;
  Iq_PG = IGD_PG*PetscSinScalar(theta) - IGD*PetscCosScalar(theta)*delta_PG + IGQ_PG*PetscCosScalar(theta) + IGQ*PetscSinScalar(theta)*delta_PG;
  Iq_QG = IGD_QG*PetscSinScalar(theta) - IGD*PetscCosScalar(theta)*delta_QG + IGQ_QG*PetscCosScalar(theta) + IGQ*PetscSinScalar(theta)*delta_QG;
  Iq_VD = IGD_VD*PetscSinScalar(theta) - IGD*PetscCosScalar(theta)*delta_VD + IGQ_VD*PetscCosScalar(theta) + IGQ*PetscSinScalar(theta)*delta_VD;
  Iq_VQ = IGD_VQ*PetscSinScalar(theta) - IGD*PetscCosScalar(theta)*delta_VQ + IGQ_VQ*PetscCosScalar(theta) + IGQ*PetscSinScalar(theta)*delta_VQ;

  //Vd = VD*PetscCosScalar(theta) - VQ*PetscSinScalar(theta);
  //Vq = VD*PetscSinScalar(theta) + VQ*PetscCosScalar(theta);
  Vd_PG = VD*PetscSinScalar(theta)*delta_PG + VQ*PetscCosScalar(theta)*delta_PG;
  Vd_QG = VD*PetscSinScalar(theta)*delta_QG + VQ*PetscCosScalar(theta)*delta_QG;
  Vd_VD = PetscCosScalar(theta) + VD*PetscSinScalar(theta)*delta_VD + VQ*PetscCosScalar(theta)*delta_VD;
  Vd_VQ = VD*PetscSinScalar(theta)*delta_VQ - PetscSinScalar(theta) + VQ*PetscCosScalar(theta)*delta_VQ;
  Vq_PG = - VD*PetscCosScalar(theta)*delta_PG + VQ*PetscSinScalar(theta)*delta_PG;
  Vq_QG = - VD*PetscCosScalar(theta)*delta_QG + VQ*PetscSinScalar(theta)*delta_QG;
  Vq_VD = PetscSinScalar(theta) - VD*PetscCosScalar(theta)*delta_VD + VQ*PetscSinScalar(theta)*delta_VD;
  Vq_VQ = - VD*PetscCosScalar(theta)*delta_VQ + PetscCosScalar(theta) + VQ*PetscSinScalar(theta)*delta_VQ;

  //Edp = Vd - genrou->Xqp*Iq;
  //Eqp = Vq + genrou->Xdp*Id;
  Edp_PG = Vd_PG - genrou->Xqp*Iq_PG;
  Edp_QG = Vd_QG - genrou->Xqp*Iq_QG;
  Edp_VD = Vd_VD - genrou->Xqp*Iq_VD;
  Edp_VQ = Vd_VQ - genrou->Xqp*Iq_VQ;
  Eqp_PG = Vq_PG + genrou->Xdp*Id_PG;
  Eqp_QG = Vq_QG + genrou->Xdp*Id_QG;
  Eqp_VD = Vq_VD + genrou->Xdp*Id_VD;
  Eqp_VQ = Vq_VQ + genrou->Xdp*Id_VQ;

  x_PG[0] = Eqp_PG;
  x_PG[1] = Edp_PG;
  x_PG[2] = Eqp_PG - (genrou->Xdp - genrou->Xl)*Id_PG;
  x_PG[3] = -Edp_PG - (genrou->Xdp - genrou->Xl)*Iq_PG;
  x_PG[4] = delta_PG;
  x_PG[5] = 0.0;
  x_PG[6] = Id_PG;
  x_PG[7] = Iq_PG;

  x_QG[0] = Eqp_QG;
  x_QG[1] = Edp_QG;
  x_QG[2] = Eqp_QG - (genrou->Xdp - genrou->Xl)*Id_QG;
  x_QG[3] = -Edp_QG - (genrou->Xdp - genrou->Xl)*Iq_QG;
  x_QG[4] = delta_QG;
  x_QG[5] = 0.0;
  x_QG[6] = Id_QG;
  x_QG[7] = Iq_QG;

  //x_VD[0] = Eqp_VD;
  //x_VD[1] = Edp_VD;
  //x_VD[2] = Eqp_VD - (genrou->Xdp - genrou->Xl)*Id_VD;
  //x_VD[3] = -Edp_VD - (genrou-Xdp - genrou->Xl)*Iq_VD;
  //x_VD[4] = delta_VD;
  //x_VD[5] = 0.0;
  //x_VD[6] = Id_VD;
  //x_VD[7] = Iq_VD;

  //x_VQ[0] = Eqp_VQ;
  //x_VQ[1] = Edp_VQ;
  //x_VQ[2] = Eqp_VQ - (genrou->Xdp - genrou->Xl)*Id_VQ;
  //x_VQ[3] = -Edp_VQ - (genrou-Xdp - genrou->Xl)*Iq_VQ;
  //x_VQ[4] = delta_VQ;
  //x_VQ[5] = 0.0;
  //x_VQ[6] = Id_VQ;
  //x_VQ[7] = Iq_VQ;

  VD_VA = -VM*PetscSinScalar(VA);
  VD_VM = PetscCosScalar(VA);
  VQ_VA = VM*PetscCosScalar(VA);
  VQ_VM = PetscSinScalar(VA);
  x_VA[0]  = Eqp_VD*VD_VA + Eqp_VQ*VQ_VA;
  x_VA[1]  = Edp_VD*VD_VA + Edp_VQ*VQ_VA;
  x_VA[2]  = (Eqp_VD-(genrou->Xdp - genrou->Xl)*Id_VD)*VD_VA + (Eqp_VQ - (genrou->Xdp - genrou->Xl)*Id_VQ)*VQ_VA;
  x_VA[3]  = (-Edp_VD - (genrou->Xdp - genrou->Xl)*Iq_VD)*VD_VA + (-Edp_VQ - (genrou->Xdp - genrou->Xl)*Iq_VQ)*VQ_VA;
  x_VA[4]  = delta_VD*VD_VA + delta_VQ*VQ_VA;
  x_VA[5]  = 0.0;
  x_VA[6]  = Id_VD*VD_VA + Id_VQ*VQ_VA;
  x_VA[7]  = Iq_VD*VD_VA + Iq_VQ*VQ_VA;

  x_VM[0]  = Eqp_VD*VD_VM + Eqp_VQ*VQ_VM;
  x_VM[1]  = Edp_VD*VD_VM + Edp_VQ*VQ_VM;
  x_VM[2]  = (Eqp_VD-(genrou->Xdp - genrou->Xl)*Id_VD)*VD_VM + (Eqp_VQ - (genrou->Xdp - genrou->Xl)*Id_VQ)*VQ_VM;
  x_VM[3]  = (-Edp_VD - (genrou->Xdp - genrou->Xl)*Iq_VD)*VD_VM + (-Edp_VQ - (genrou->Xdp - genrou->Xl)*Iq_VQ)*VQ_VM;
  x_VM[4]  = delta_VD*VD_VM + delta_VQ*VQ_VM;
  x_VM[5]  = 0.0;
  x_VM[6]  = Id_VD*VD_VM + Id_VQ*VQ_VM;
  x_VM[7]  = Iq_VD*VD_VM + Iq_VQ*VQ_VM;

  for (i=0;i<8;i++) rows[i] = dynlocglob+i;
  ierr = SetMatrixValues(ICp,8,rows,1,&VAloc,x_VA);CHKERRQ(ierr);
  ierr = SetMatrixValues(ICp,8,rows,1,&VMloc,x_VM);CHKERRQ(ierr);
  ierr = SetMatrixValues(ICp,8,rows,1,&PGloc,x_PG);CHKERRQ(ierr);
  ierr = SetMatrixValues(ICp,8,rows,1,&QGloc,x_QG);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetEquationTypes_Genrou - Gets the number and indices of differential and algebraic equations
                                                   
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
PetscErrorCode DYNGenModelGetEquationTypes_Genrou(DYNGenModel dyngen,PetscInt *ndiff,PetscInt *nalg, PetscInt *eqtype)
{
  DYNGenrou genrou;
  PetscInt  nd=6,na=2;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&genrou);CHKERRQ(ierr);

  eqtype[0] = eqtype[1] = eqtype[4] = eqtype[5] = DIFF_EQ;
  eqtype[6] = eqtype[7] = ALG_EQ;
  if(PetscAbsScalar(genrou->Td0dp) < PETSC_SMALL) {
    nd--; na++;
    eqtype[2] = ALG_EQ;
    genrou->Td0dp = 1.0;
  } else {
    eqtype[2] = DIFF_EQ;
  }

  if(PetscAbsScalar(genrou->Tq0dp) < PETSC_SMALL) {
    nd--; na++;
    eqtype[3] = ALG_EQ;
    genrou->Tq0dp = 1.0;
  } else {
    eqtype[3] = DIFF_EQ;
  }

  *ndiff = nd;
  *nalg  = na;

  PetscFunctionReturn(0);
}

/*
  DYNGenModelDAERHSJacobian_Genrou - Computes the RHS Jacobian of the Genrou DAE equations

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

   The location for the first variable for the genrou model should be obtained by
   DYNGenModelGetFirstVariableLocation(dyngen,&loc). x[loc] = 1st genrou variable,
   x[loc+n-1] = nth genrou variable
*/ 
PetscErrorCode DYNGenModelDAERHSJacobian_Genrou(DYNGenModel dyngen,Mat J,PetscReal t,PetscScalar VD,PetscScalar VQ,PetscScalar *x,PetscInt dynlocglob,PetscInt V_loc[],PetscInt I_loc[])
{
  PetscErrorCode ierr;
  DYNGenrou      genrou;
  PetscScalar Eqp,Edp,psi1d,psi2q,delta,dw,Id,Iq,Efd,Pmech,dPmechddw=0.0;
  PetscInt    row[2],col[7];
  PetscScalar val[7];
  PetscInt    loc;
  DYNExcModel dynexc;
  DYNTurbgovModel dynturbgov;
  PetscInt    dPmechdxtgov_num=0,dPmechdxtgov_locglob[4];
  PetscScalar dPmechdxtgov[4];
  PetscInt    i,dEfddXexc_num,dEfddXexc_locglob[2];
  PetscScalar dEfddXexc[2];
  
  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&genrou);CHKERRQ(ierr);
  ierr = DYNGenModelGetFirstVariableLocation(dyngen,&loc);CHKERRQ(ierr);

  Eqp   = x[loc];
  Edp   = x[loc+1];
  psi1d = x[loc+2];
  psi2q = x[loc+3];
  delta = x[loc+4];
  dw    = x[loc+5];
  Id    = x[loc+6];
  Iq    = x[loc+7];

  ierr = DYNGenModelGetDynExc(dyngen,&dynexc);CHKERRQ(ierr);
  if(dynexc) {
    ierr = DYNExcModelGetFieldVoltage(dynexc,t,x,&Efd,&dEfddXexc_num,dEfddXexc,dEfddXexc_locglob);CHKERRQ(ierr);
  } else Efd = genrou->Efd;

  ierr = DYNGenModelGetDynTurbgov(dyngen,&dynturbgov);CHKERRQ(ierr);
  if(dynturbgov) {
    ierr = DYNTurbgovModelGetMechanicalPower(dynturbgov,t,x,&Pmech,&dPmechddw,&dPmechdxtgov_num,dPmechdxtgov,dPmechdxtgov_locglob);CHKERRQ(ierr);
  } else Pmech  = genrou->Pm;

  /*df[0]_dxdyn */
  row[0] = dynlocglob;
  col[0] = dynlocglob;    /* df0_dEqp */
  col[1] = dynlocglob+2;  /* df0_dpsi1d */
  col[2] = dynlocglob+6;  /* df0_dId */

  PetscScalar temp1,temp2,temp3;
  temp1 = genrou->Xd - genrou->Xdp;
  temp2 = (genrou->Xdp - genrou->Xddp)/PetscPowScalar((genrou->Xdp - genrou->Xl),2);
  temp3 = genrou->Xdp - genrou->Xl;
  val[0] = (-1 - temp1*-temp2*-1)/genrou->Td0p;
  val[1] = -temp1*-temp2/genrou->Td0p;
  val[2] = -temp1*(1 - temp2*temp3)/genrou->Td0p;

  ierr = SetMatrixValues(J,1,row,3,col,val);CHKERRQ(ierr);
  if(dynexc) {
    for(i=0; i < dEfddXexc_num; i++) dEfddXexc[i] *= 1/genrou->Td0p;
    ierr = SetMatrixValues(J,1,row,dEfddXexc_num,dEfddXexc_locglob,dEfddXexc);CHKERRQ(ierr);
  }

  /* df[1]_dxdyn */
  row[0] = dynlocglob+1;
  col[0] = dynlocglob+1; col[1] = dynlocglob+3; col[2] = dynlocglob+7;
  temp1 = genrou->Xq - genrou->Xqp;
  temp2 = (genrou->Xqp - genrou->Xddp)/PetscPowScalar((genrou->Xqp - genrou->Xl),2);
  temp3 = genrou->Xqp - genrou->Xl;
  val[0] = (-1 + temp1*-temp2)/genrou->Tq0p;
  val[1] = temp1*-temp2/genrou->Tq0p;
  val[2] = temp1*(1 - temp2*temp3)/genrou->Tq0p;

  ierr = SetMatrixValues(J,1,row,3,col,val);CHKERRQ(ierr);

  /* df[2]_dxdyn */
  row[0] = dynlocglob+2;
  col[0] = dynlocglob; col[1] = dynlocglob+2; col[2] = dynlocglob+6;
  val[0] = 1/genrou->Td0dp;
  val[1] = -1/genrou->Td0dp;
  val[2] = -(genrou->Xdp - genrou->Xl)/genrou->Td0dp;

  ierr = SetMatrixValues(J,1,row,3,col,val);CHKERRQ(ierr);

  /* df[3]_dxdyn */
  row[0] = dynlocglob+3;
  col[0] = dynlocglob+1; col[1] = dynlocglob+3; col[2] = dynlocglob+7;
  val[0] = -1/genrou->Tq0dp;
  val[1] = -1/genrou->Tq0dp;
  val[2] = -(genrou->Xqp - genrou->Xl)/genrou->Tq0dp;

  ierr = SetMatrixValues(J,1,row,3,col,val);CHKERRQ(ierr);

  /*df[4]_dxdyn */
  row[0] = dynlocglob+4;
  col[0] = dynlocglob+5;
  val[0] = (2*PETSC_PI*freq);
  ierr = SetMatrixValues(J,1,row,1,col,val);CHKERRQ(ierr);

  /*df[5]_dxdyn */
  PetscScalar tempd1,tempd2,tempq1,tempq2,const1;
  tempd1 = (genrou->Xddp - genrou->Xl)/(genrou->Xdp - genrou->Xl);
  tempd2 = (genrou->Xdp - genrou->Xddp)/(genrou->Xdp - genrou->Xl);
  tempq1 = (genrou->Xddp - genrou->Xl)/(genrou->Xqp - genrou->Xl);
  tempq2 = (genrou->Xqp - genrou->Xddp)/(genrou->Xqp - genrou->Xl);
  const1 = 1/(2*genrou->H);

  row[0] = dynlocglob + 5;
  col[0]  = dynlocglob; col[1] = dynlocglob+1; col[2] = dynlocglob+2; col[3] = dynlocglob+3; col[4] = dynlocglob+5; col[5] = dynlocglob+6; col[6] = dynlocglob+7;
  val[0] = -tempd1*Iq*const1;
  val[1] = -tempq1*Id*const1;
  val[2] = -tempd2*Iq*const1;
  val[3] =  tempq2*Id*const1;
  val[4] = (-(Pmech-genrou->D*dw)/PetscPowScalar((1+dw),2) + dPmechddw/(1+dw) - genrou->D/(1+dw))*const1;
  val[5] = (-tempq1*Edp + tempq2*psi2q)*const1;
  val[6] = (-tempd1*Eqp - tempd2*psi1d)*const1;

  ierr = SetMatrixValues(J,1,row,7,col,val);CHKERRQ(ierr);

  /* Partial derivatives w.r.t. turbine governor variables */
  if(dynturbgov) {
    for(i=0; i < dPmechdxtgov_num; i++) dPmechdxtgov[i] *= const1/(1+dw);
    ierr = SetMatrixValues(J,1,row,dPmechdxtgov_num,dPmechdxtgov_locglob,dPmechdxtgov);CHKERRQ(ierr);
  }

  /* df[6]_dxdyn and df[6]_dVDQ */
  row[0] = dynlocglob+6;
  col[0] = dynlocglob;  col[1] = dynlocglob+2; col[2] = dynlocglob+4; col[3] = dynlocglob+6; col[4] = V_loc[0]; col[5] = V_loc[1];
  val[0] = -tempd1;
  val[1] = -tempd2;
  val[2] = -VD*PetscSinScalar(delta) + VQ*PetscCosScalar(delta);
  val[3] = genrou->Xddp;
  val[4] = PetscCosScalar(delta);
  val[5] = PetscSinScalar(delta);

  ierr = SetMatrixValues(J,1,row,6,col,val);CHKERRQ(ierr);

  /* df[7]_dxdyn and df[6]_dVDQ */
  row[0] = dynlocglob+7;
  col[0] = dynlocglob+1;  col[1] = dynlocglob+3; col[2] = dynlocglob+4; col[3] = dynlocglob+7; col[4] = V_loc[0]; col[5] = V_loc[1];
  val[0] = -tempq1;
  val[1] = tempq2;
  val[2] = VD*PetscCosScalar(delta) + VQ*PetscSinScalar(delta);
  val[3] = -genrou->Xddp;
  val[4] = PetscSinScalar(delta);
  val[5] = -PetscCosScalar(delta);

  ierr = SetMatrixValues(J,1,row,6,col,val);CHKERRQ(ierr);

 /* dIGD_dxdyn */
  row[0] = I_loc[0];
  col[0] = dynlocglob+4; col[1] = dynlocglob+6; col[2] = dynlocglob+7;
  val[0] = (Id*PetscCosScalar(delta) - Iq*PetscSinScalar(delta));
  val[1] = PetscSinScalar(delta);
  val[2] = PetscCosScalar(delta);
 
  ierr = SetMatrixValues(J,1,row,3,col,val);CHKERRQ(ierr);

  /*dIGQ_dxdyn */
  row[0] = I_loc[1];
  col[0] = dynlocglob+4; col[1] = dynlocglob+6; col[2] = dynlocglob+7;
  val[0] = (Id*PetscSinScalar(delta) + Iq*PetscCosScalar(delta));
  val[1] = -PetscCosScalar(delta);
  val[2] =  PetscSinScalar(delta);

  ierr = SetMatrixValues(J,1,row,3,col,val);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNGenModelDAERHSJacobianP_Genrou - Sets the Jacobian of the Genrou equations w.r.t. parameters

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
PetscErrorCode DYNGenModelDAERHSJacobianP_Genrou(DYNGenModel dyngen, PetscReal t,const PetscScalar *x,Mat jacP,PetscInt dynlocglob,PetscInt Valoc,PetscInt Vmloc,PetscInt Pgloc,PetscInt Qgloc)
{
  PetscErrorCode ierr;
  DYNGenrou      genrou;
  PetscScalar    dw;
  PetscInt       row,col,loc;
  PetscScalar    value;
  DYNTurbgovModel dynturbgov;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&genrou);CHKERRQ(ierr);
  ierr = DYNGenModelGetDynTurbgov(dyngen,&dynturbgov);CHKERRQ(ierr);

  if(dynturbgov) {
    /* The turbine governor model will calculate this */ 
    PetscFunctionReturn(0);
  }

  ierr = DYNGenModelGetFirstVariableLocation(dyngen,&loc);CHKERRQ(ierr);

  dw = x[loc+5];

  row = dynlocglob + 5;
  col = Pgloc;
  value = (-1/((1+dw)*(1+dw)))/(2*genrou->H);

  ierr = SetMatrixValues(jacP,1,&row,1,&col,&value);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

/*
  DYNGenModelDAEFWDRHSJacobianP_Genrou - Sets the Jacobian of the Genrou equations w.r.t. parameters

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
PetscErrorCode DYNGenModelDAEFWDRHSJacobianP_Genrou(DYNGenModel dyngen, PetscReal t,const PetscScalar *x,Vec *jacP,PetscInt dynlocglob,PetscInt Valoc,PetscInt Vmloc,PetscInt Pgloc,PetscInt Qgloc)
{
  PetscErrorCode ierr;
  DYNGenrou      genrou;
  PetscScalar    dw;
  PetscInt       row,loc;
  //PetscInt       col;
  PetscScalar    value;

  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&genrou);CHKERRQ(ierr);
  ierr = DYNGenModelGetFirstVariableLocation(dyngen,&loc);CHKERRQ(ierr);

  dw = x[loc+5];

  row = dynlocglob + 5;
  //col = Pgloc;
  value = (-1/((1+dw)*(1+dw)))/(2*genrou->H);

  ierr = VecSetValue(jacP[Pgloc],row,value,INSERT_VALUES);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNGenModelDAERHSFunction_Genrou - Computes the rhs of the DAE function for the GENROU model 
                                and also returns the generator currents in network refernce frame

  Input Parameters:
+ dyngen - the dynamic generator model
. t      - the current time
. x      - array of all the variables for the bus on which this generator is incident
. VD     - real-part of the bus voltage
- VQ     - imaginary part of the bus voltage

  Output Parameters:
+ f      - array of rhs of DAE equations for the genrou model 
. IGD    - real-part of generator current in network reference frame
. IGQ    - imaginary-part of generator current in network reference frame

  Notes:
   The locations for the variables (and the corresponding locations to insert entries in f)
   for this generator model should be obtained by DYNGenModelGetFirstVariableLocation()
*/
PetscErrorCode DYNGenModelDAERHSFunction_Genrou(DYNGenModel dyngen,PetscReal t,PetscScalar VD, PetscScalar VQ, PetscScalar *x,PetscScalar *f,PetscScalar *IGD,PetscScalar *IGQ)
{
  PetscErrorCode ierr;
  DYNGenrou   genrou;
  PetscScalar Eqp,Edp,psi1d,psi2q,delta,dw,Id,Iq,Efd,Pmech;
  PetscScalar Vd,Vq;
  PetscInt    loc;
  DYNExcModel dynexc;
  DYNTurbgovModel dynturbgov;
  
  PetscFunctionBegin;
  ierr = DYNGenModelGetModelData(dyngen,(void**)&genrou);CHKERRQ(ierr);

  ierr = DYNGenModelGetFirstVariableLocation(dyngen,&loc);CHKERRQ(ierr);

  Eqp   = x[loc];
  Edp   = x[loc+1];
  psi1d = x[loc+2];
  psi2q = x[loc+3];
  delta = x[loc+4];
  dw    = x[loc+5];
  Id    = x[loc+6];
  Iq    = x[loc+7];

  ierr = DYNGenModelGetDynExc(dyngen,&dynexc);CHKERRQ(ierr);
  if(!dynexc) {
    Efd = genrou->Efd;
  } else {
    ierr = DYNExcModelGetFieldVoltage(dynexc,t,x,&Efd,NULL,NULL,NULL);CHKERRQ(ierr);
  }

  ierr = DYNGenModelGetDynTurbgov(dyngen,&dynturbgov);CHKERRQ(ierr);
  if(!dynturbgov) {
    Pmech = genrou->Pm;
  } else {
    ierr = DYNTurbgovModelGetMechanicalPower(dynturbgov,t,x,&Pmech,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  }

  /* GENROU differential equations */
  f[loc]   = (-Eqp - (genrou->Xd - genrou->Xdp)*(Id - (genrou->Xdp - genrou->Xddp)/PetscPowScalar((genrou->Xdp - genrou->Xl),2)*(psi1d + (genrou->Xdp - genrou->Xl)*Id - Eqp)) + Efd)/genrou->Td0p;
  f[loc+1] = (-Edp + (genrou->Xq - genrou->Xqp)*(Iq - (genrou->Xqp - genrou->Xddp)/PetscPowScalar((genrou->Xqp - genrou->Xl),2)*(psi2q + (genrou->Xqp - genrou->Xl)*Iq + Edp)))/genrou->Tq0p;
  f[loc+2] = (-psi1d + Eqp - (genrou->Xdp - genrou->Xl)*Id)/genrou->Td0dp;
  f[loc+3] = (-psi2q - Edp - (genrou->Xqp - genrou->Xl)*Iq)/genrou->Tq0dp;
  f[loc+4] = dw*(2*PETSC_PI*freq);

  PetscScalar tempd1,tempd2,tempq1,tempq2;
  tempd1 = (genrou->Xddp - genrou->Xl)/(genrou->Xdp - genrou->Xl);
  tempd2 = (genrou->Xdp - genrou->Xddp)/(genrou->Xdp - genrou->Xl);
  tempq1 = (genrou->Xddp - genrou->Xl)/(genrou->Xqp - genrou->Xl);
  tempq2 = (genrou->Xqp - genrou->Xddp)/(genrou->Xqp - genrou->Xl);
  f[loc+5] = ((Pmech-genrou->D*dw)/(1+dw) - tempd1*Eqp*Iq - tempd2*psi1d*Iq - tempq1*Edp*Id + tempq2*psi2q*Id - (genrou->Xddp - genrou->Xddp)*Id*Iq)/(2*genrou->H);

  Vd = VD*PetscSinScalar(delta) - VQ*PetscCosScalar(delta);
  Vq = VD*PetscCosScalar(delta) + VQ*PetscSinScalar(delta);

  f[loc+6] = genrou->Xddp*Id - tempd1*Eqp - tempd2*psi1d + Vq;
  f[loc+7] = -genrou->Xddp*Iq - tempq1*Edp + tempq2*psi2q + Vd;

  *IGD =  Id*PetscSinScalar(delta) + Iq*PetscCosScalar(delta);
  *IGQ = -Id*PetscCosScalar(delta) + Iq*PetscSinScalar(delta);

  PetscFunctionReturn(0);
}

/*
  DYNGenModelCreate_Genrou - Class constructor for GENROU model
*/
PetscErrorCode DYNGenModelCreate_Genrou(DYNGenModel dyngen)
{
  DYNGenrou genrou;
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = PetscCalloc1(1,&genrou);CHKERRQ(ierr);
  dyngen->data = (void*)genrou;

  /* Inherit the ops */
  dyngen->ops.readdata             = DYNGenModelReadData_Genrou;
  dyngen->ops.destroy              = DYNGenModelDestroy_Genrou;
  dyngen->ops.getnvar              = DYNGenModelGetNvar_Genrou;
  dyngen->ops.getsizeof            = DYNGenModelGetSizeof_Genrou;
  dyngen->ops.setinitialconditions = DYNGenModelSetInitialConditions_Genrou;
  dyngen->ops.setinitialconditionsp = DYNGenModelSetInitialConditionsP_Genrou;
  dyngen->ops.setinitialfieldvoltagediff = DYNGenModelSetInitialFieldVoltageDiff_Genrou;
  dyngen->ops.daerhsfunction       = DYNGenModelDAERHSFunction_Genrou;
  dyngen->ops.daerhsjacobian       = DYNGenModelDAERHSJacobian_Genrou;
  dyngen->ops.daerhsjacobianp      = DYNGenModelDAERHSJacobianP_Genrou;
  dyngen->ops.daefwdrhsjacobianp   = DYNGenModelDAEFWDRHSJacobianP_Genrou;
  dyngen->ops.getequationtypes     = DYNGenModelGetEquationTypes_Genrou;
  dyngen->ops.getbusnumid          = DYNGenModelGetBusnumID_Genrou;
  dyngen->ops.getfieldcurrent      = DYNGenModelGetFieldCurrent_Genrou;
  dyngen->ops.getinitialfieldvoltage = DYNGenModelGetInitialFieldVoltage_Genrou;
  dyngen->ops.getinitialmechanicalpower = DYNGenModelGetInitialMechanicalPower_Genrou;
  dyngen->ops.getfrequency         = DYNGenModelGetFrequency_Genrou;
  dyngen->ops.setfrequency         = DYNGenModelSetFrequency_Genrou;
  dyngen->ops.getdfreqdstate       = DYNGenModelGetdFreqdState_Genrou;
  dyngen->ops.getspeeddeviation    = DYNGenModelGetSpeedDeviation_Genrou;
  dyngen->ops.getspeeddeviationlocation = DYNGenModelGetSpeedDeviationLocation_Genrou;
  dyngen->ops.setevent             = DYNGenModelSetEvent_Genrou;

  genrou->freqatmax = genrou->freqatmin = 0;

  PetscFunctionReturn(0);
}

