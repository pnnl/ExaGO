#include "dyncompload.h"
#include <private/dynloadmodelsimpl.h>

/*
  DYNEventMonitor_Compload - Event monitoring routine for this load model
*/
PetscErrorCode DYNEventMonitor_Compload(DYNLoadModel dynload, PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscScalar *fval)
{
  PetscErrorCode ierr;
  DYNCOMPLOAD    compload;
  PetscScalar    Vm;

  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNLoadModelGetModelData(dynload,(void**)&compload);CHKERRQ(ierr);

  fval[0] = fval[1] = fval[2] = 999;

  if(compload->motstatus && compload->mot_on) {
    Vm = PetscSqrtScalar(VD*VD+VQ*VQ);
    if(compload->low_voltage_flag) {
      fval[0] = 999; /* Disable event */
      fval[1] = Vm - compload->Vi;
      fval[2] = (t - compload->t_timer_start) - compload->Ti*1/freq; /* Ti is in cycles, freq is defined in constants.h */
    } else {
      fval[0] = Vm - compload->Vi;
      fval[1] = 999;
      fval[2] = 999;
    }
  }
  PetscFunctionReturn(0);
}

/*
  DYNEventPostDGamma_Compload - Compute derivative of gamma
*/
PetscErrorCode DYNEventPostDGamma_Compload(DYNLoadModel dynload, PetscInt nevents, PetscInt ev_list[], PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscScalar *dgdx, PetscScalar *dgdp,PetscInt Valoc, PetscInt Vmloc, PetscInt Pgloc, PetscInt Qgloc)
{
  PetscErrorCode ierr;
  DYNCOMPLOAD      compload;
  PetscInt       loc;

  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNLoadModelGetModelData(dynload,(void**)&compload);CHKERRQ(ierr);

  /* Get the location for the first exciter variable */
  ierr = DYNLoadModelGetFirstVariableLocation(dynload,&loc);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNEventPostFunction_Compload - Post event routine for this load model
*/
PetscErrorCode DYNEventPostFunction_Compload(DYNLoadModel dynload, PetscInt nevents, PetscInt ev_list[], PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscBool forward_solve, PetscBool *solve_algebraic)
{
  PetscErrorCode ierr;
  DYNCOMPLOAD      compload;
  PetscScalar      Vm;

  PetscFunctionBegin;

  Vm = PetscSqrtScalar(VD*VD + VQ*VQ);

  /* Get the model data */
  ierr = DYNLoadModelGetModelData(dynload,(void**)&compload);CHKERRQ(ierr);

  if(compload->motstatus) {
    if(ev_list[0] == 0 && Vm < compload->Vi) { /* Vm < Vi */
      compload->t_timer_start = t;
      compload->low_voltage_flag = 1;
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%4.3f s] Bus %d: Low voltage triggering motor voltage relay, Vm %4.3f < Vthresh %4.3f\n",t,compload->bus_i,Vm,compload->Vi);CHKERRQ(ierr);
    } else if(ev_list[0] == 1 && Vm > compload->Vi) { /* Vm > Vi */
      compload->t_timer_start = 0;
      compload->low_voltage_flag = 0;
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%4.3f s] Bus %d: Voltage recovered resetting motor voltage relay, Vm %4.3f > Vthresh %4.3f\n",t,compload->bus_i,Vm,compload->Vi);CHKERRQ(ierr);
    } else if(ev_list[0] == 2 && Vm < compload->Vi) { /* Timer expired, disconnect motor */
      compload->motstatus = 0;
      compload->mot_on = 0;
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%4.3f s] Bus %d: Tripping induction motor load, Bus Vm = %4.3f < Threshold = %4.3f\n",t,compload->bus_i,Vm,compload->Vi);CHKERRQ(ierr);
      *solve_algebraic = PETSC_TRUE;
    }
  }
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelSetEvent_Compload - Sets the event info for the COMPLOAD load model

  Input Parameters
. dynload - the dynamic load model

  Output Parameters
+ nmons - the number of event monitors
. direction - the event directions
. terminate - flags for termination when the event is located
. eventfcn  - the function describing the event condition
- posteventfcn - An optional function that gets called when the event is located
*/
PetscErrorCode DYNLoadModelSetEvent_Compload(DYNLoadModel dynload,PetscInt *nmons,PetscInt *direction,PetscBool *terminate,PetscErrorCode (**eventfcn)(DYNLoadModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*), PetscErrorCode (**posteventfcn)(DYNLoadModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscBool,PetscBool*), PetscErrorCode (**posteventdgamma)(DYNLoadModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscInt))
{
  PetscErrorCode ierr;
  DYNCOMPLOAD compload;

  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&compload);CHKERRQ(ierr);

  if(compload->motstatus) {
    *nmons = 3;
    direction[0] = -1; direction[1] = 1; direction[2] = 0;
    terminate[0] = terminate[1] = terminate[2] = PETSC_FALSE;
    *eventfcn = DYNEventMonitor_Compload;
    *posteventfcn = DYNEventPostFunction_Compload;
    *posteventdgamma = DYNEventPostDGamma_Compload;
  } else *nmons = 0;
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelReadData_Compload - Parses the line containining the COMPLOAD data and
  populates it in the Compload data structure

  Input Parameters:
+ dynload - The load base model
. line   - the file from the line
. mbase  - the machine base (NULL if not applicable)
- sbase  - the system base (NULL if not applicable)

  Notes: 
    The type of load model should be set prior to calling this routine via DYNLoadModelSetType()

    The conversion from machine base to system base,if applicable, should be done in this routine
*/
PetscErrorCode DYNLoadModelReadData_Compload(DYNLoadModel dynload,char* line,PetscScalar Pd,PetscScalar sbase)
{
  PetscErrorCode ierr;
  DYNCOMPLOAD compload;
  PetscScalar cratio;

  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&compload);CHKERRQ(ierr);
  /*
    BUSI "COMPLOAD" ID PL QL IP IQ YP YQ IT RA XA XM R1 X1 R2 X2 E1 SE1 E2 SE2 MBASE PMULT \
    H VI TI TB D Tnom
  */
  sscanf(line,"%d,'COMPLOAD',%[^,],%lf,%lf,%lf,%lf,%lf,%lf,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",&compload->bus_i,compload->id,&compload->pl,&compload->ql,&compload->ip,&compload->iq,&compload->yp,&compload->yq,\
	 &compload->IT,&compload->Ra,&compload->Xa,&compload->Xm,&compload->R1,&compload->X1,&compload->R2,&compload->X2,\
	 &compload->E1,&compload->SE1,&compload->E2,&compload->SE2,&compload->Mbase,&compload->Pmult,\
	 &compload->H,&compload->Vi,&compload->Ti,&compload->Tb,&compload->D,&compload->Tnom);

  compload->Vm_thresh = 0.8;
  
  if(PetscAbsScalar(compload->R2-999) > 1e-6 || PetscAbsScalar(compload->X2-999) > 1e-6) {
    compload->cages = 2;
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Bus %d: No support for double cage induction motor model ",compload->bus_i);
  }

  if(PetscAbsScalar(compload->Mbase) < 1e-6) {
    /* Use PMULT*PL as the Mbase */
    cratio = compload->Pmult*Pd/sbase;
  } else {
    cratio = compload->Mbase/sbase;
  }
  compload->H *= cratio;
  compload->D *= cratio;

  compload->Ra /= cratio;
  compload->Xa /= cratio;
  compload->Xm /= cratio;
  
  compload->X0 = compload->Xa + compload->Xm;
  compload->Xp = compload->Xa + compload->X1*compload->Xm/(compload->X1 + compload->Xm);
  compload->Tp0 = (compload->X1 + compload->Xm)/(w_s*compload->R1);

  /* Convert from percentages to fraction */
  compload->pl /= 100.0;
  compload->ql /= 100.0;
  compload->ip /= 100.0;
  compload->iq /= 100.0;
  compload->yp /= 100.0;
  compload->yq /= 100.0;

  /* Compute motor parameters */
  compload->pmot = 1 - (compload->pl + compload->ip + compload->yp);
  if(PetscAbsScalar(compload->pmot) < 1e-6) compload->motstatus = PETSC_FALSE;
  else {
    compload->motstatus = PETSC_TRUE;
    compload->mot_on = 1;
  }

  compload->s0 = 0.05;

  PetscFunctionReturn(0);
}

/*
  DYNLoadModelGetBusnumID_Compload - Returns the bus number and ID associated with this load model

  Input Parameters:
. dynload - the dynamic load model object

  Output Parameters:
+ busnum - the bus number 
- loadid  - the load ID
*/
PetscErrorCode DYNLoadModelGetBusnumID_Compload(DYNLoadModel dynload,PetscInt *busnum,char **loadid)
{
  PetscErrorCode ierr;
  DYNCOMPLOAD         compload;

  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&compload);CHKERRQ(ierr);
  *busnum = compload->bus_i;
  *loadid  = compload->id;
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelDestroy_Compload - Destroys the Compload object data

  Input Parameters
. DYNLoadModel - the DYNLoadModel object

  Notes:
  Called when DYNLoadModelDestroy() is called
*/
PetscErrorCode DYNLoadModelDestroy_Compload(DYNLoadModel dynload)
{
  PetscErrorCode ierr;
  DYNCOMPLOAD         compload;

  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&compload);CHKERRQ(ierr);
  ierr = PetscFree(compload);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelGetNvar_Compload - Returns the number of variables for the Compload model

  Input parameters
. dynload - the dynamic load object

  Output Parameters
. nvar - number of variables for this model
*/
PetscErrorCode DYNLoadModelGetNvar_Compload(DYNLoadModel dynload,PetscInt *nvar)
{
  DYNCOMPLOAD         compload;
  PetscErrorCode      ierr;

  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&compload);CHKERRQ(ierr);


  if(!compload->motstatus) *nvar = 0;
  else {
    if(compload->cages == 1) *nvar = 3;
    else *nvar = 5;
  }
  
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelGetSizeof_Compload - Returns the size of Compload struct

  Input Parameters:
. dynload - the dynamic load model
 
  Output Parameters:
. size - size of the Compload object obtained from sizeof()
*/
PetscErrorCode DYNLoadModelGetSizeof_Compload(DYNLoadModel dynload,PetscInt *size)
{
  PetscFunctionBegin;
  *size = sizeof(struct _p_DYNCOMPLOAD);
  PetscFunctionReturn(0);
}

PetscErrorCode InitialSlipFunction(SNES snes, Vec X, Vec F,void *ctx)
{
  PetscErrorCode ierr;
  DYNCOMPLOAD    compload=(DYNCOMPLOAD)ctx;
  PetscScalar    *f;
  const PetscScalar *x;
  PetscScalar zr1,zi1,zr,zi,s;
  PetscScalar Xm=compload->Xm,Ra=compload->Ra,Xa=compload->Xa,R1=compload->R1,X1=compload->X1;
  PetscScalar Vm=compload->Vm0;
  PetscScalar Pmot;


  PetscFunctionBegin;

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(F,&f);CHKERRQ(ierr);

  s = x[0];

  if(compload->cages == 1) {
    /* Single-cage */
    zr1 = (Xm*R1/s*(Xm+X1) - X1*Xm*R1/s)/(R1*R1/(s*s) + (Xm+X1)*(Xm+X1));
    zi1 = (Xm*R1/s*R1/s + X1*Xm*(Xm+X1))/(R1*R1/(s*s) + (Xm+X1)*(Xm+X1));
    
    zr = Ra + zr1;
    zi = Xa + zi1;
    
    Pmot = Vm*Vm*zr/(zr*zr + zi*zi);

    f[0] = Pmot - compload->pmot;
  }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(F,&f);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}
  
/* Computes the initial slip for the motor */
static PetscErrorCode GetMotorInitialSlip(DYNCOMPLOAD compload)
{
  PetscErrorCode ierr;
  SNES           snes;
  Vec            X,F;
  Mat            J;
  SNESConvergedReason reason;
  PetscInt       row=0;

  PetscFunctionBegin;

  ierr = SNESCreate(PETSC_COMM_SELF,&snes);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_SELF,&X);CHKERRQ(ierr);
  ierr = VecSetSizes(X,1,1);CHKERRQ(ierr);
  ierr = VecSetFromOptions(X);CHKERRQ(ierr);
  ierr = VecDuplicate(X,&F);CHKERRQ(ierr);

  ierr = VecSetValue(X,0,0.001,INSERT_VALUES);CHKERRQ(ierr);

  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,1,1,1,PETSC_NULL,&J);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,J,J,SNESComputeJacobianDefault,PETSC_NULL);CHKERRQ(ierr);

  ierr = SNESSetFunction(snes,F,InitialSlipFunction,(void*)compload);CHKERRQ(ierr);

  ierr = SNESSetFromOptions(snes);

  ierr = SNESSolve(snes,NULL,X);CHKERRQ(ierr);
  ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);
  if(reason < 1) {
    SETERRQ1(PETSC_COMM_SELF,0,"Motor load initial slip calculation at bus %d did not converge",compload->bus_i);CHKERRQ(ierr);
  }
  ierr = VecGetValues(X,1,&row,&compload->s0);CHKERRQ(ierr);

  ierr = VecDestroy(&X);CHKERRQ(ierr);
  ierr = VecDestroy(&F);CHKERRQ(ierr);
  ierr = MatDestroy(&J);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
/*
  DYNLoadModelSetInitialConditions_Compload - Sets the initial conditions (x(t0)) for the COMPLOAD model

  Input Parameters:
+ dynload - the DYNLoadModel object
. Pl     - load real power output
. Ql     - load reactive power output
. VD     - real component of the complex bus voltage
. VQ     - imaginary component of the complex bus voltage

  Output Parameters:
. x - the initial conditions for the COMPLOAD model

  Notes the locations to insert the values should be obtained by DYNLoadModelGetFirstVariableLocation()
*/
PetscErrorCode DYNLoadModelSetInitialConditions_Compload(DYNLoadModel dynload,PetscScalar Pl, PetscScalar Ql, PetscScalar VD, PetscScalar VQ, PetscScalar *x)
{
  PetscErrorCode ierr;
  DYNCOMPLOAD    compload;
  PetscScalar Vm;
  PetscInt    loc;

  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&compload);CHKERRQ(ierr);

  Vm = PetscSqrtScalar(VD*VD + VQ*VQ);
  compload->Vm0 = Vm;
  ierr = DYNLoadModelGetFirstVariableLocation(dynload,&loc);CHKERRQ(ierr);

  if(compload->motstatus) {
    compload->pmot *= Pl;
    ierr = GetMotorInitialSlip(compload);CHKERRQ(ierr);

    PetscScalar zr1,zi1,zr,zi;
    PetscScalar Xm=compload->Xm,Ra=compload->Ra,Xa=compload->Xa,R1=compload->R1,X1=compload->X1;
    PetscScalar Vm=compload->Vm0,s=compload->s0;
    PetscScalar Pmot,Qmot,Id,Iq,Qcomp;

    if(compload->cages == 1) {
      /* Single-cage */
      zr1 = (Xm*R1/s*(Xm+X1) - X1*Xm*R1/s)/(R1*R1/(s*s) + (Xm+X1)*(Xm+X1));
      zi1 = (Xm*R1/s*R1/s + X1*Xm*(Xm+X1))/(R1*R1/(s*s) + (Xm+X1)*(Xm+X1));
    
      zr = Ra + zr1;
      zi = Xa + zi1;
    
      Pmot = Vm*Vm*zr/(zr*zr + zi*zi);
      Qmot = Vm*Vm*zi/(zr*zr + zi*zi);

      Id = (Pmot*VD + Qmot*VQ)/(Vm*Vm);
      Iq = (-Qmot*VD + Pmot*VQ)/(Vm*Vm);

      x[loc] = compload->s0; /* s */
      x[loc+1] = VD - Ra*Id + compload->Xp*Iq; /* Edp */
      x[loc+2] = VQ - Ra*Iq - compload->Xp*Id; /* Eqp */

      compload->Tnom = (x[loc+1]*Id + x[loc+2]*Iq)/PetscPowScalar(1-s,compload->D);
      
      /* Calculate compensating shunt */
      Qcomp = Ql*(1 - (compload->ql + compload->iq + compload->yq)) - Qmot;

      compload->Bcomp = -Qcomp/(Vm*Vm);
    }
  }

  /* Set up for the static part */
  compload->pl *= Pl; 
  compload->ql *= Ql;
  compload->ip *= Pl/Vm;
  compload->iq *= Ql/Vm;
  compload->yp *= Pl/(Vm*Vm);
  compload->yq *= Ql/(Vm*Vm);


  PetscFunctionReturn(0);
}

/*
  DYNLoadModelGetEquationTypes_Compload - Gets the number and indices of differential and algebraic equations
                                                   
  Input Parameters:
. dynload - the DYNLoadModel object

  Output Parameters:
+ ndiff - number of differential equations
. nalg  - number of algebraic equations
- eqtype - an array of size nvar that has
            eqtype[i] = DIFF_EQ if equation i is a differential
            vartype[i] = ALG_VAR if equation i is an algebraic

  NOTES: The user does not need to create the array vartype
*/
PetscErrorCode DYNLoadModelGetEquationTypes_Compload(DYNLoadModel dynload,PetscInt *ndiff,PetscInt *nalg, PetscInt *eqtype)
{
  DYNCOMPLOAD compload;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&compload);CHKERRQ(ierr);
  *nalg = 0;
  if(compload->motstatus == 0) *ndiff = 0;
  else {
    if(compload->cages == 1) {
      *ndiff = 3;
      eqtype[0] = eqtype[1] = eqtype[2] = DIFF_EQ;
    }
    else *ndiff = 5;
  }
      
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelDAERHSJacobian_Compload - Computes the RHS Jacobian of the Compload DAE equations

  Input Parameters:
+ dynload     - the dynamic load model
. t          - the current time
. VD         - real-part of comploadlex bus voltage
. VQ         - imaginary-part of comploadlex bus voltage
. x          - all variables at the bus where this load is incident
. dynlocglob - starting location of load variables in the global vector 
. V_loc      - global location of VD and VQ V_loc[0] = VD_loc V_loc[1] = VQ_loc
- I_loc      - global location of ID and IQ I_loc[0] = ID_loc I_loc[1] = IQ_loc

  Output Parameters:
. J          - the Jacobian matrix

  Notes:
   For each bus, the network equations are ordered as [ID_loc;IQ_loc]

   The network equations are written as I_gen - I_net - I_load = 0. So the
   load model should set dIG_dxdyn in the Jacobian

   The location for the first variable for the compload model should be obtained by
   DYNLoadModelGetFirstVariableLocation(dynload,&loc). x[loc] = 1st compload variable,
   x[loc+n-1] = nth compload variable
*/ 
PetscErrorCode DYNLoadModelDAERHSJacobian_Compload(DYNLoadModel dynload,Mat J,PetscReal t,PetscScalar VD,PetscScalar VQ,PetscScalar *x,PetscInt dynlocglob,PetscInt V_loc[],PetscInt I_loc[])
{
  PetscErrorCode ierr;
  DYNCOMPLOAD    compload;
  PetscScalar    Iloaddp,Iloadqp,Iloaddi,Iloadqi;
  PetscScalar    dIloaddp_dVD=0.0,dIloaddp_dVQ=0.0,dIloadqp_dVD=0.0,dIloadqp_dVQ=0.0;
  PetscScalar    dIloaddi_dVD=0.0,dIloaddi_dVQ=0.0,dIloadqi_dVD=0.0,dIloadqi_dVQ=0.0;
  PetscScalar    dIloaddz_dVD=0.0,dIloaddz_dVQ=0.0,dIloadqz_dVD=0.0,dIloadqz_dVQ=0.0;
  PetscScalar    dIloaddmot_dVD=0.0,dIloaddmot_dVQ=0.0,dIloadqmot_dVD=0.0,dIloadqmot_dVQ=0.0;
  PetscScalar    Vm,Vm2;
  PetscInt       row,col[6];
  PetscScalar    val[6];
  PetscScalar    Ra,Xp,X0,Tp0;
  PetscInt       loc;
  PetscScalar    s,Edp,Eqp,Id,Iq;

  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&compload);CHKERRQ(ierr);

  Vm = PetscSqrtScalar(VD*VD + VQ*VQ);
  Vm2 = Vm*Vm;
  
  if(Vm > compload->Vm_thresh) {
    Iloaddp = (compload->pl*VD + compload->ql*VQ)/Vm2;
    Iloadqp = (-compload->ql*VD + compload->pl*VQ)/Vm2;

    dIloaddp_dVD = compload->pl/Vm2  - 2*Iloaddp*VD/Vm2;
    dIloaddp_dVQ = compload->ql/Vm2  - 2*Iloaddp*VQ/Vm2;
    dIloadqp_dVD = -compload->ql/Vm2 - 2*Iloadqp*VD/Vm2;
    dIloadqp_dVQ = compload->pl/Vm2  - 2*Iloadqp*VQ/Vm2;

  } else {
    Iloaddp = (compload->pl*VD + compload->ql*VQ)/(compload->Vm_thresh*compload->Vm_thresh);
    Iloadqp = (-compload->ql*VD + compload->pl*VQ)/(compload->Vm_thresh*compload->Vm_thresh);

    dIloaddp_dVD = compload->pl/(compload->Vm_thresh*compload->Vm_thresh);
    dIloaddp_dVQ = compload->ql/(compload->Vm_thresh*compload->Vm_thresh);
    dIloadqp_dVD = -compload->ql/(compload->Vm_thresh*compload->Vm_thresh);
    dIloadqp_dVQ = compload->pl/(compload->Vm_thresh*compload->Vm_thresh);
  }
  
  if(Vm > 0.5) {
    Iloaddi = (compload->ip*VD + compload->iq*VQ)/Vm;
    Iloadqi = (-compload->iq*VD + compload->ip*VQ)/Vm;

    dIloaddi_dVD = compload->ip/Vm  - Iloaddi*VD/Vm2;
    dIloaddi_dVQ = compload->ip/Vm  - Iloaddi*VQ/Vm2;
    dIloadqi_dVD = -compload->ip/Vm - Iloadqi*VD/Vm2;
    dIloadqi_dVQ = compload->ip/Vm  - Iloadqi*VQ/Vm2;
  } else {
    Iloaddi = (compload->ip*VD + compload->iq*VQ)/0.5;
    Iloadqi = (-compload->iq*VD + compload->ip*VQ)/0.5;

    dIloaddi_dVD = compload->ip/0.5;
    dIloaddi_dVQ = compload->ip/0.5;
    dIloadqi_dVD = -compload->ip/0.5;
    dIloadqi_dVQ = compload->ip/0.5;
  }

  dIloaddz_dVD = compload->yp;
  dIloaddz_dVQ = compload->yq;
  dIloadqz_dVD = -compload->yq;
  dIloadqz_dVQ = compload->yp;

  if(compload->motstatus) {
    PetscScalar dId_dVD, dId_dVQ, dId_dEdp, dId_dEqp;
    PetscScalar dIq_dVD, dIq_dVQ, dIq_dEdp, dIq_dEqp;
    PetscScalar c1,c2;
    PetscScalar dTe_dEdp,dTe_dEqp,dTe_dVD,dTe_dVQ;

    Ra=compload->Ra;
    Xp=compload->Xp;
    X0=compload->X0;
    Tp0 = compload->Tp0;
    c1 = Ra/(Ra*Ra + Xp*Xp);
    c2 = Xp/(Ra*Ra + Xp*Xp);

    ierr = DYNLoadModelGetFirstVariableLocation(dynload,&loc);CHKERRQ(ierr);

    s = x[loc];
    Edp = x[loc+1];
    Eqp = x[loc+2];

    Id = c1*(VD - Edp) + c2*(VQ - Eqp);
    Iq = -c2*(VD - Edp) + c1*(VQ - Eqp);

    dId_dVD = c1;   dId_dVQ = c2;
    dIq_dVD = -c2;  dIq_dVQ = c1;

    dId_dEdp = -c1; dId_dEqp = -c2;
    dIq_dEdp = c2;  dIq_dEqp = -c1;

    dTe_dEdp = Edp*dId_dEdp + Id + Eqp*dIq_dEdp;
    dTe_dEqp = Edp*dId_dEqp + Eqp*dIq_dEqp + Iq;

    dTe_dVD = Edp*dId_dVD + Eqp*dIq_dVD;
    dTe_dVQ = Edp*dId_dVQ + Eqp*dIq_dVQ;

    dIloaddmot_dVD = dId_dVD;
    dIloaddmot_dVQ = dId_dVQ - compload->Bcomp;
    dIloadqmot_dVD = dIq_dVD + compload->Bcomp;
    dIloadqmot_dVQ = dIq_dVQ;

    row = dynlocglob;
    col[0] = dynlocglob; col[1] = dynlocglob+1; col[2] = dynlocglob+2; 
    col[3] = V_loc[0]; col[4] = V_loc[1];

    val[0] = -compload->D*compload->Tnom*PetscPowScalar(1-s,compload->D-1)/(2*compload->H);
    val[1] = -dTe_dEdp/(2*compload->H);
    val[2] = -dTe_dEqp/(2*compload->H);
    val[3] = -dTe_dVD/(2*compload->H);
    val[4] = -dTe_dVQ/(2*compload->H);

    ierr = SetMatrixValues(J,1,&row,5,col,val);CHKERRQ(ierr);

    row = dynlocglob + 1;
    val[0] =  w_s*Eqp;
    val[1] = -1*(1 + (X0 - Xp)*dIq_dEdp)/Tp0;
    val[2] =  w_s*s -1*(X0 - Xp)*dIq_dEqp/Tp0;
    val[3] = -1*(X0 - Xp)*dIq_dVD/Tp0;
    val[4] = -1*(X0 - Xp)*dIq_dVQ/Tp0;

    ierr = SetMatrixValues(J,1,&row,5,col,val);CHKERRQ(ierr);

    row = dynlocglob + 2;
    val[0] = -w_s*Edp;
    val[1] = -w_s*s + (X0 - Xp)*dId_dEdp/Tp0;
    val[2] = -1*(1 - (X0 - Xp)*dId_dEqp)/Tp0;
    val[3] = (X0 - Xp)*dId_dVD/Tp0;
    val[4] = (X0 - Xp)*dId_dVQ/Tp0;

    ierr = SetMatrixValues(J,1,&row,5,col,val);CHKERRQ(ierr);

    if(compload->mot_on) {
      /* dIloaddmot_dEdp, dIloaddmot_dEqp */
      row = I_loc[0];
      col[0] = dynlocglob+1; col[1] = dynlocglob+2;
      val[0] = -dId_dEdp; val[1] = -dId_dEqp;
      
      ierr = SetMatrixValues(J,1,&row,2,col,val);CHKERRQ(ierr);
      
      /* dIloadqmot_dEdp, dIloadqmot_dEqp */
      row = I_loc[1];
      col[0] = dynlocglob+1; col[1] = dynlocglob+2;
      val[0] = -dIq_dEdp; val[1] = -dIq_dEqp;
      
      ierr = SetMatrixValues(J,1,&row,2,col,val);CHKERRQ(ierr);
    }
  }

  row  = I_loc[0];
  col[0] = V_loc[0]; col[1] = V_loc[1];
  val[0] = -(dIloaddp_dVD + dIloaddi_dVD + dIloaddz_dVD + compload->mot_on*dIloaddmot_dVD);
  val[1] = -(dIloaddp_dVQ + dIloaddi_dVQ + dIloaddz_dVQ +  compload->mot_on*dIloaddmot_dVQ);

  ierr = SetMatrixValues(J,1,&row,2,col,val);CHKERRQ(ierr);

  row = I_loc[1];
  col[0] = V_loc[0]; col[1] = V_loc[1];
  val[0] = -(dIloadqp_dVD + dIloadqi_dVD + dIloadqz_dVD +  compload->mot_on*dIloadqmot_dVD);
  val[1] = -(dIloadqp_dVQ + dIloadqi_dVQ + dIloadqz_dVQ +  compload->mot_on*dIloadqmot_dVQ);
  
  ierr = SetMatrixValues(J,1,&row,2,col,val);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNLoadModelDAERHSJacobianP_Compload - Sets the Jacobian of the Compload equations w.r.t. parameters

  Input Parameters:
+ dynload - dynloadmodel object
. t      - current time
. x      - state variables for this bus
. dynlocglob - starting location of the variables for this dynload model
. Valoc  - location of parameter voltage angle (Va)
. Vmloc  - location of parameter voltage magnitude (Vm)
. Pgloc  - location of parameter load real power output (Pg)
- Qgloc  - location of parameter load reactive power output (Qg)

  Output Parameters:
. jacP - Matrix of partial derivatives of DYN DAE equations w.r.t. parameters
*/
PetscErrorCode DYNLoadModelDAERHSJacobianP_Compload(DYNLoadModel dynload, PetscReal t,const PetscScalar *x,Mat jacP,PetscInt dynlocglob,PetscInt valoc,PetscInt vmloc)
{
  DYNCOMPLOAD         compload;
  PetscInt       row[2],col[1];
  PetscScalar    VD,VQ,VM,ip_vm,iq_vm,yp_vm,yq_vm,Iloaddp_vm,Iloadqp_vm,Iloaddi_vm,Iloadqi_vm,Iloaddz_vm,Iloadqz_vm,val[2];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&compload);CHKERRQ(ierr);
  VD = x[0]; VQ = x[1];
  VM = PetscSqrtScalar(VD*VD + VQ*VQ);

  /* assume bus->vm == compload->Vm0 */
  ip_vm = -compload->ip/compload->Vm0;
  iq_vm = -compload->iq/compload->Vm0;
  yp_vm = -2.*compload->yp/compload->Vm0;
  yq_vm = -2.*compload->yq/compload->Vm0;

  if (VM > compload->Vm_thresh) {
    Iloaddp_vm = 0.;
    Iloadqp_vm = 0.;
  } else {
    Iloaddp_vm = -2.*(compload->pl*VD + compload->ql*VQ)/(compload->Vm0*compload->Vm0*compload->Vm0);
    Iloadqp_vm = -2.*(-compload->ql*VD + compload->pl*VQ)/(compload->Vm0*compload->Vm0*compload->Vm0);
  }

  if (VM > 0.5) {
    Iloaddi_vm = (ip_vm*VD + iq_vm*VQ)/VM;
    Iloadqi_vm = (-iq_vm*VD + ip_vm*VQ)/VM;
  } else {
    Iloaddi_vm = (ip_vm*VD + iq_vm*VQ)/compload->Vm0  - (compload->ip*VD + compload->iq*VQ)/(compload->Vm0*compload->Vm0);
    Iloadqi_vm = (-iq_vm*VD + ip_vm*VQ)/compload->Vm0 -(-compload->iq*VD + compload->ip*VQ)/(compload->Vm0*compload->Vm0);
  }

  Iloaddz_vm = yp_vm*VD + yq_vm*VQ;
  Iloadqz_vm = -yq_vm*VD + yp_vm*VQ;

  row[0] = dynlocglob; row[1] = dynlocglob+1;
  col[0] = vmloc;
  val[0] = -Iloadqp_vm - Iloadqi_vm - Iloadqz_vm;
  val[1] = -Iloaddp_vm - Iloaddi_vm - Iloaddz_vm;

  ierr = SetMatrixValues(jacP,2,row,1,col,val);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelDAEFWDRHSJacobianP_Compload - Sets the Jacobian of the Compload equations w.r.t. parameters

  Input Parameters:
+ dynload - dynloadmodel object
. t      - current time
. x      - state variables for this bus
. dynlocglob - starting location of the variables for this dynload model
. Valoc  - location of parameter voltage angle (Va)
. Vmloc  - location of parameter voltage magnitude (Vm)
. Pgloc  - location of parameter load real power output (Pg)
- Qgloc  - location of parameter load reactive power output (Qg)

  Output Parameters:
. jacP - Matrix of partial derivatives of DYN DAE equations w.r.t. parameters
*/
PetscErrorCode DYNLoadModelDAEFWDRHSJacobianP_Compload(DYNLoadModel dynload, PetscReal t,const PetscScalar *x,Vec *jacP,PetscInt dynlocglob,PetscInt Valoc,PetscInt Vmloc)
{
  PetscErrorCode ierr;
  DYNCOMPLOAD      compload;
  PetscInt    loc;

  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&compload);CHKERRQ(ierr);
  ierr = DYNLoadModelGetFirstVariableLocation(dynload,&loc);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNLoadModelDAERHSFunction_Compload - Computes the rhs of the DAE function for the COMPLOAD model 
                                and also returns the load currents in network refernce frame

  Input Parameters:
+ dynload - the dynamic load model
. t      - the current time
. x      - array of all the variables for the bus on which this load is incident
. VD     - real-part of the bus voltage
- VQ     - imaginary part of the bus voltage

  Output Parameters:
+ f      - array of rhs of DAE equations for the compload model 
. ILD    - real-part of load current in network reference frame
. ILQ    - imaginary-part of load current in network reference frame

  Notes:
   The locations for the variables (and the corresponding locations to insert entries in f)
   for this load model should be obtained by DYNLoadModelGetFirstVariableLocation()
*/
PetscErrorCode DYNLoadModelDAERHSFunction_Compload(DYNLoadModel dynload,PetscReal t,PetscScalar VD, PetscScalar VQ, PetscScalar *x,PetscScalar *f,PetscScalar *ILD,PetscScalar *ILQ)
{
  PetscErrorCode ierr;
  DYNCOMPLOAD    compload;
  PetscScalar    Iloaddp,Iloadqp,Iloaddi,Iloadqi,Iloaddz,Iloadqz,Iloaddmot=0.0,Iloadqmot=0.0;
  PetscScalar    Vm,Vm2,s,Edp,Eqp,Tm,Te,Id,Iq,c1,c2;
  PetscInt       loc;
  PetscScalar    Ra,Xp,X0,Tp0;
  
  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&compload);CHKERRQ(ierr);

  Vm = PetscSqrtScalar(VD*VD + VQ*VQ);
  Vm2 = Vm*Vm;

  if(compload->motstatus) {
    Ra=compload->Ra;
    Xp=compload->Xp;
    X0=compload->X0;
    Tp0 = compload->Tp0;

    ierr = DYNLoadModelGetFirstVariableLocation(dynload,&loc);CHKERRQ(ierr);

    s = x[loc];
    Edp = x[loc+1];
    Eqp = x[loc+2];

    Tm = compload->Tnom*PetscPowScalar(1-s,compload->D);
    
    c1 = Ra/(Ra*Ra + Xp*Xp);
    c2 = Xp/(Ra*Ra + Xp*Xp);
    
    Id = c1*(VD - Edp) + c2*(VQ - Eqp);
    Iq = -c2*(VD - Edp) + c1*(VQ - Eqp);

    s   = x[loc];
    Edp = x[loc+1];
    Eqp = x[loc+2];

    //    ierr = PetscPrintf(PETSC_COMM_SELF,"Motor %d: Id = %4.3f, Iq = %4.3f\n",compload->bus_i,Id,Iq);CHKERRQ(ierr);
    Te = Edp*Id + Eqp*Iq;
    
    f[loc] = (Tm - Te)/(2*compload->H);
    f[loc+1] =  w_s*s*Eqp - (Edp + (X0 - Xp)*Iq)/Tp0;
    f[loc+2] = -w_s*s*Edp - (Eqp - (X0 - Xp)*Id)/Tp0;

    Iloaddmot = compload->mot_on*(Id -  compload->Bcomp*VQ);
    Iloadqmot = compload->mot_on*(Iq +  compload->Bcomp*VD);

  }

  if(Vm > compload->Vm_thresh) {
    Iloaddp = (compload->pl*VD + compload->ql*VQ)/Vm2;
    Iloadqp = (-compload->ql*VD + compload->pl*VQ)/Vm2;
  } else {
    Iloaddp = (compload->pl*VD + compload->ql*VQ)/(compload->Vm_thresh*compload->Vm_thresh);
    Iloadqp = (-compload->ql*VD + compload->pl*VQ)/(compload->Vm_thresh*compload->Vm_thresh);
  }

  if(Vm > 0.5) {
    Iloaddi = (compload->ip*VD + compload->iq*VQ)/Vm;
    Iloadqi = (-compload->iq*VD + compload->ip*VQ)/Vm;
  } else {
    Iloaddi = (compload->ip*VD + compload->iq*VQ)/0.5;
    Iloadqi = (-compload->iq*VD + compload->ip*VQ)/0.5;
  }
    
  Iloaddz = compload->yp*VD + compload->yq*VQ;
  Iloadqz = -compload->yq*VD + compload->yp*VQ;
  
  *ILD = Iloaddp + Iloaddi + Iloaddz + Iloaddmot;
  *ILQ = Iloadqp + Iloadqi + Iloadqz + Iloadqmot;

  PetscFunctionReturn(0);
}

/*
  DYNLoadModelCreate_Compload - Class constructor for COMPLOAD model
*/
PetscErrorCode DYNLoadModelCreate_COMPLOAD(DYNLoadModel dynload)
{
  DYNCOMPLOAD compload;
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = PetscCalloc1(1,&compload);CHKERRQ(ierr);
  compload->cages = 1;
  compload->t_timer_start = 0;
  compload->low_voltage_flag = 0;

  /* Inherit the ops */
  dynload->ops.readdata             = DYNLoadModelReadData_Compload;
  dynload->ops.destroy              = DYNLoadModelDestroy_Compload;
  dynload->ops.getnvar              = DYNLoadModelGetNvar_Compload;
  dynload->ops.getsizeof            = DYNLoadModelGetSizeof_Compload;
  dynload->ops.setinitialconditions = DYNLoadModelSetInitialConditions_Compload;
  dynload->ops.daerhsfunction       = DYNLoadModelDAERHSFunction_Compload;
  dynload->ops.daerhsjacobian       = DYNLoadModelDAERHSJacobian_Compload;
  dynload->ops.daerhsjacobianp      = DYNLoadModelDAERHSJacobianP_Compload;
  dynload->ops.daefwdrhsjacobianp   = DYNLoadModelDAEFWDRHSJacobianP_Compload;
  dynload->ops.getequationtypes     = DYNLoadModelGetEquationTypes_Compload;
  dynload->ops.getbusnumid          = DYNLoadModelGetBusnumID_Compload;
  dynload->ops.setevent             = DYNLoadModelSetEvent_Compload;

  dynload->data = (void*)compload;  

  
  PetscFunctionReturn(0);
}

