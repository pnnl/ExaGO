#include "dynieeet1.h"
#include <private/dynexcmodelsimpl.h>

/*
  DYNEventMonitor_IEEET1 - Event monitoring routine for this exciter model
*/
PetscErrorCode DYNEventMonitor_IEEET1(DYNExcModel dynexc, PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscScalar *fval)
{
  PetscErrorCode ierr;
  DYNIEEET1      ieeet1;
  PetscInt       loc;
  PetscScalar    VT,VR,Efd,RF;
  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNExcModelGetModelData(dynexc,(void**)&ieeet1);CHKERRQ(ierr);

  /* Get the location for the first exciter variable */
  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);

  VT  = x[loc];
  VR  = x[loc+1];
  Efd = x[loc+2];
  RF  = x[loc+3];

  if(!ieeet1->VRatmax) fval[0] = ieeet1->VRmax - VR;
  else   fval[0] = (-VR + ieeet1->KA*RF - ieeet1->KA*ieeet1->KF/ieeet1->TF*Efd + ieeet1->KA*(ieeet1->Vref - VT))/ieeet1->TA;

  if(!ieeet1->VRatmin) fval[1] = ieeet1->VRmin - VR;
  else   fval[1] = (-VR + ieeet1->KA*RF - ieeet1->KA*ieeet1->KF/ieeet1->TF*Efd + ieeet1->KA*(ieeet1->Vref - VT))/ieeet1->TA;
  PetscFunctionReturn(0);
}

/*
  DYNEventPostDGamma_IEEET1 - Compute derivative of gamma
*/
PetscErrorCode DYNEventPostDGamma_IEEET1(DYNExcModel dynexc, PetscInt nevents, PetscInt ev_list[], PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscScalar *dgdx, PetscScalar *dgdp,PetscInt Valoc, PetscInt Vmloc, PetscInt Pgloc, PetscInt Qgloc)
{
  PetscErrorCode ierr;
  DYNIEEET1      ieeet1;
  PetscInt       loc;
  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNExcModelGetModelData(dynexc,(void**)&ieeet1);CHKERRQ(ierr);

  /* Get the location for the first exciter variable */
  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);

  /* at one time, only one event is triggered */
  if(ev_list[0] == 0) { /* Max VR */
    if(!ieeet1->VRatmax) dgdx[loc+1] = -1;
    else {
      dgdx[loc]   = -ieeet1->KA/ieeet1->TA;
      dgdx[loc+1] = -1./ieeet1->TA;
      dgdx[loc+2] = -ieeet1->KA*ieeet1->KF/ieeet1->TF/ieeet1->TA;
      dgdx[loc+3] = ieeet1->KA/ieeet1->TA;

      dgdp[Vmloc] = ieeet1->KA/ieeet1->TA;
    }
  }else{
   if(!ieeet1->VRatmin) dgdx[loc+1] = -1;
    else {
      dgdx[loc]   = -ieeet1->KA/ieeet1->TA;
      dgdx[loc+1] = -1./ieeet1->TA;
      dgdx[loc+2] = -ieeet1->KA*ieeet1->KF/ieeet1->TF/ieeet1->TA;
      dgdx[loc+3] = ieeet1->KA/ieeet1->TA;

      dgdp[Vmloc] = ieeet1->KA/ieeet1->TA;
    }
  }
  PetscFunctionReturn(0);
}

/*
  DYNEventPostFunction_IEEET1 - Post event routine for this exciter model
*/
PetscErrorCode DYNEventPostFunction_IEEET1(DYNExcModel dynexc, PetscInt nevents, PetscInt ev_list[], PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscBool forward_solve, PetscBool *solve_algebraic)
{
  PetscErrorCode ierr;
  DYNIEEET1      ieeet1;
  PetscFunctionBegin;

  //ierr = DYNEventPostDGamma_IEEET1(dynexc,nevents,ev_list,t,VD,VQ,x,dgdx,dgdp);CHKERRQ(ierr);

  /* Get the model data */
  ierr = DYNExcModelGetModelData(dynexc,(void**)&ieeet1);CHKERRQ(ierr);

  if(ev_list[0] == 0) { /* Max VR */
    if(!ieeet1->VRatmax) ieeet1->VRatmax = 1;
    else ieeet1->VRatmax = 0;
  } else { /* Min VR */
    if(!ieeet1->VRatmin) ieeet1->VRatmin = 1;
    else ieeet1->VRatmin = 0;
  }

  PetscFunctionReturn(0);
}


/*
  DYNExcModelSetEvent_IEEET1 - Sets the event info for the IEEET1 exciter model

  Input Parameters
. dynexc - the dynamic exciter model

  Output Parameters
+ nmons - the number of event monitors
. direction - the event directions
. terminate - flags for termination when the event is located
. eventfcn  - the function describing the event condition
- posteventfcn - An optional function that gets called when the event is located
*/
PetscErrorCode DYNExcModelSetEvent_IEEET1(DYNExcModel dynexc,PetscInt *nmons,PetscInt *direction,PetscBool *terminate,PetscErrorCode (**eventfcn)(DYNExcModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*), PetscErrorCode (**posteventfcn)(DYNExcModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscBool,PetscBool*), PetscErrorCode (**posteventdgamma)(DYNExcModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscInt))
{

  PetscFunctionBegin;
  *nmons = 2;
  direction[0] = -1; direction[1] = 1;
  terminate[0] = terminate[1] = PETSC_FALSE;
  *eventfcn = DYNEventMonitor_IEEET1;
  *posteventfcn = DYNEventPostFunction_IEEET1;
  *posteventdgamma = DYNEventPostDGamma_IEEET1;
  PetscFunctionReturn(0);
}


/*
  DYNExcModelReadData_IEEET1 - Parses the line containining the IEEET1 data and
  populates it in the IEEET1 data structure

  Input Parameters:
+ dynexc - the dynamic exciter struct
- line - The line in the dyr file

*/
PetscErrorCode DYNExcModelReadData_IEEET1(DYNExcModel dynexc,char* line)
{
  PetscErrorCode ierr;
  DYNIEEET1 ieeet1;
  PetscInt  unused;

  PetscFunctionBegin;
  ierr = DYNExcModelGetModelData(dynexc,(void**)&ieeet1);CHKERRQ(ierr);
  /* Read the IEEET1 data */
  sscanf(line,"%d, 'IEEET1',%[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%lf,%lf,%lf,%lf",\
	 &ieeet1->bus_i,ieeet1->id,&ieeet1->TR,&ieeet1->KA,&ieeet1->TA,&ieeet1->VRmax,\
	 &ieeet1->VRmin,&ieeet1->KE,&ieeet1->TE,&ieeet1->KF,&ieeet1->TF,&unused,\
	 &ieeet1->E1,&ieeet1->SE1,&ieeet1->E2,&ieeet1->SE2);

  PetscFunctionReturn(0);
}

/*
  DYNExcModelGetBusnumID_IEEET1 - Returns the bus number and ID associated with this exciter model

  Input Parameters:
. dynexc - the dynamic exciter model object

  Output Parameters:
+ busnum - the bus number 
- excid  - the exciter ID
*/
PetscErrorCode DYNExcModelGetBusnumID_IEEET1(DYNExcModel dynexc,PetscInt *busnum,char **excid)
{
  PetscErrorCode ierr;
  DYNIEEET1      ieeet1;

  PetscFunctionBegin;
  ierr = DYNExcModelGetModelData(dynexc,(void**)&ieeet1);CHKERRQ(ierr);
  *busnum = ieeet1->bus_i;
  *excid  = ieeet1->id;
  PetscFunctionReturn(0);
}

/*
  DYNExcModelDestroy_IEEET1 - Destroys the IEEET1 object data

  Input Parameters
. DYNExcModel - the DYNExcModel object

  Notes:
  Called when DYNExcModelDestroy() is called
*/
PetscErrorCode DYNExcModelDestroy_IEEET1(DYNExcModel dynexc)
{
  PetscErrorCode ierr;
  DYNIEEET1      ieeet1;

  PetscFunctionBegin;
  ierr = DYNExcModelGetModelData(dynexc,(void**)&ieeet1);CHKERRQ(ierr);
  ierr = PetscFree(ieeet1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNExcModelGetNvar_IEEET1 - Returns the number of variables for the IEEET1 model

  Input parameters
. dynexc - the dynamic exciter object

  Output Parameters
. nvar - number of variables for this model
*/
PetscErrorCode DYNExcModelGetNvar_IEEET1(DYNExcModel dynexc,PetscInt *nvar)
{
  PetscFunctionBegin;
  *nvar = DYNIEEET1_nvar;
  PetscFunctionReturn(0);
}

/*
  DYNExcModelGetSizeof_IEEET1 - Returns the size of IEEET1 struct

  Input Parameters:
. dynexc - the dynamic exciter model
 
  Output Parameters:
. size - size of the IEEET1 object obtained from sizeof()
*/
PetscErrorCode DYNExcModelGetSizeof_IEEET1(DYNExcModel dynexc,PetscInt *size)
{
  PetscFunctionBegin;
  *size = sizeof(struct _p_DYNIEEET1);
  PetscFunctionReturn(0);
}

/*
  DYNExcModelGetFieldVoltage - Returns the exciter field voltage and its partial derivative w.r.t to
                                  exciter variables (optional)

  Input Parameters:
+ dynexc - the exciter model
. t      - the current time
. xdyn   - array of the variables for this bus

  Output Parameters:
+ Efd           - Field voltage for the exciter
. dEfddXexc_num - Number of partial derivatives of Efd w.r.t. exciter variables.
. dEfddXexc     - An array of partial derivatives of the field voltage w.r.t. the number
                 of exciter states given by dEfddXexc_num
- dEfddXexc_loc - Global locations for the partial derivatives

  Notes - dEfddXexc_num, dEfddXexc and dEfddXexc_loc are optionally returned if the passed-in arrays are not NULL.
          
*/
PetscErrorCode DYNExcModelGetFieldVoltage_IEEET1(DYNExcModel dynexc,PetscReal t, PetscScalar *xdyn,PetscScalar *Efd,PetscInt *dEfddXexc_num,PetscScalar dEfddXexc[],PetscInt dEfddXexc_loc[])
{
  PetscErrorCode ierr;
  DYNIEEET1      ieeet1;
  PetscInt       loc;

  PetscFunctionBegin;
  /* Get Model data */
  ierr = DYNExcModelGetModelData(dynexc,(void**)&ieeet1);CHKERRQ(ierr);

  /* Get the location for the first variable for this exciter model */
  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);

  *Efd = xdyn[loc+2];
  if(dEfddXexc_num) *dEfddXexc_num   = 1;
  if(dEfddXexc)     dEfddXexc[0]     = 1.0;
  if(dEfddXexc_loc)  dEfddXexc_loc[0] = loc + 2;

  PetscFunctionReturn(0);
}

/*
  DYNExcModelSetInitialConditions_IEEET1 - Sets the initial conditions (x(t0)) for the IEEET1 model

  Input Parameters:
+ dynexc - the DYNExcModel object
. VD     - real component of the complex bus voltage
. VQ     - imaginary component of the complex bus voltage

  Output Parameters:
. x - the array for the variables for the bus on which the generator is incident.

  Notes:
   The location for the first variable of this exciter model can be obtained by DYNExcModelGetFirstVariableLocation()
   and then the initial conditions populated in the x array using this location information
   x[loc] = 1st variable
   x[loc+1] = 2nd variable
*/
PetscErrorCode DYNExcModelSetInitialConditions_IEEET1(DYNExcModel dynexc,PetscScalar VD, PetscScalar VQ, PetscScalar *x)
{
  PetscErrorCode ierr;
  DYNIEEET1 ieeet1;
  DYNGenModel dyngen;
  PetscScalar Vm;
  PetscScalar Efd,RF,VR,VT;
  PetscInt    loc;
  PetscScalar alpha,SE;

  PetscFunctionBegin;
  ierr = DYNExcModelGetModelData(dynexc,(void**)&ieeet1);CHKERRQ(ierr);
  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);

  ierr = DYNExcModelGetDynGen(dynexc,&dyngen);CHKERRQ(ierr);

  /* Voltage magnitude */
  Vm = PetscSqrtScalar(VD*VD + VQ*VQ);

  /* Get the field voltage Efd0 from the generator */
  ierr = DYNGenModelGetInitialFieldVoltage(dyngen,&Efd);CHKERRQ(ierr);

  RF = ieeet1->KF/ieeet1->TF*Efd;

  alpha = PetscSqrtScalar(ieeet1->SE1*ieeet1->E1/(ieeet1->SE2*ieeet1->E2));
  ieeet1->satA = (alpha*ieeet1->E2 - ieeet1->E1)/(alpha - 1);
  ieeet1->satB = ieeet1->SE1*ieeet1->E1/PetscPowScalar((ieeet1->E1-ieeet1->satA),2);
  ieeet1->Efdthresh = ieeet1->satA;

  SE  = 0.0;
  if (Efd > ieeet1->Efdthresh) SE = ieeet1->satB*PetscPowScalar((Efd - ieeet1->satA),2)/Efd;

  VR  =  (ieeet1->KE + SE)*Efd;
  VT = Vm;

  x[loc]    = VT;
  x[loc+1]  = VR;
  x[loc+2]  = Efd;
  x[loc+3]  = RF;

  /* Constants -> Reference voltage */
  ieeet1->Vref = Vm + VR/ieeet1->KA;

  /* Set the initial flags */
  ieeet1->VRatmax = ieeet1->VRatmin = 0;

  PetscFunctionReturn(0);
}

/*
  DYNExcModelSetInitialConditionsP_IEEET1 - Sets the differentiation of initial conditions (x(t0)) for the IEEET1 model

  Input Parameters:
+ dynexc - the DYNExcModel object
. VA     - angle of the bus voltage
. VM     - magnitude of the complex bus voltage

  Output Parameters:
. x_VA - the derivative of initial conditions to VA
. x_VM - the derivative of initial conditions to VMi

  Notes:
   The location for the first variable of this exciter model can be obtained by DYNExcModelGetFirstVariableLocation()
   and then the initial conditions populated in the x array using this location information
   x[loc] = 1st variable
   x[loc+1] = 2nd variable
*/
PetscErrorCode DYNExcModelSetInitialConditionsP_IEEET1(DYNExcModel dynexc,PetscScalar PG,PetscScalar QG,PetscScalar VA,PetscScalar VM,PetscInt PGloc,PetscInt QGloc,PetscInt VAloc,PetscInt VMloc,Mat ICp,PetscInt dynlocglob)
{
  DYNIEEET1   ieeet1;
  DYNGenModel dyngen;
  PetscScalar Efd,VT_VA,VT_VM,SE_Efd,VR_PG,VR_QG,VR_VA,VR_Efd,RF_PG,RF_QG,RF_VA,RF_VM,RF_Efd,Efd_PG,Efd_QG,Efd_VA,Efd_VM;
  PetscInt    i,rows[4];
  PetscScalar alpha,SE,x_PG[4],x_QG[4],x_VA[4],x_VM[4];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DYNExcModelGetModelData(dynexc,(void**)&ieeet1);CHKERRQ(ierr);
  ierr = DYNExcModelGetDynGen(dynexc,&dyngen);CHKERRQ(ierr);

  /* Get the field voltage Efd0 from the generator */
  ierr = DYNGenModelGetInitialFieldVoltage(dyngen,&Efd);CHKERRQ(ierr);
  ierr = DYNGenModelSetInitialFieldVoltageDiff(dyngen,PG,QG,VA,VM,&Efd_PG,&Efd_QG,&Efd_VA,&Efd_VM);CHKERRQ(ierr);

  alpha = PetscSqrtScalar(ieeet1->SE1*ieeet1->E1/(ieeet1->SE2*ieeet1->E2));
  ieeet1->satA = (alpha*ieeet1->E2 - ieeet1->E1)/(alpha - 1.);
  ieeet1->satB = ieeet1->SE1*ieeet1->E1/PetscPowScalar((ieeet1->E1-ieeet1->satA),2);
  ieeet1->Efdthresh = ieeet1->satA;

  //VT = PetscSqrtScalar(VD*VD + VQ*VQ);
  VT_VA = 0.;
  VT_VM = 1.;

  SE  = 0.;
  if (Efd > ieeet1->Efdthresh) SE = ieeet1->satB*PetscPowScalar((Efd - ieeet1->satA),2)/Efd;
  SE_Efd = 0.;
  if (Efd > ieeet1->Efdthresh) SE_Efd = ieeet1->satB*(1.-ieeet1->satA*ieeet1->satA/(Efd*Efd));
  //VR  =  (ieeet1->KE + SE)*Efd;
  VR_Efd = ieeet1->KE + SE + Efd*SE_Efd;
  VR_PG = VR_Efd*Efd_PG;
  VR_QG = VR_Efd*Efd_QG;
  VR_VA = VR_Efd*Efd_VA;

  //RF = ieeet1->KF/ieeet1->TF*Efd;
  RF_Efd = ieeet1->KF/ieeet1->TF;
  RF_PG = RF_Efd*Efd_PG;
  RF_QG = RF_Efd*Efd_QG;
  RF_VA = RF_Efd*Efd_VA;
  RF_VM = RF_Efd*Efd_VM;

  x_PG[0] = 0.;
  x_PG[1] = VR_PG;
  x_PG[2] = Efd_PG;
  x_PG[3] = RF_PG;

  x_QG[0] = 0.;
  x_QG[1] = VR_QG;
  x_QG[2] = Efd_QG;
  x_QG[3] = RF_QG;

  x_VA[0] = VT_VA;
  x_VA[1] = VR_VA;
  x_VA[2] = Efd_VA;
  x_VA[3] = RF_VA;

  x_VM[0] = VT_VM;
  x_VM[1] = VR_VA;
  x_VM[2] = Efd_VM;
  x_VM[3] = RF_VM;

  for (i=0;i<4;i++) rows[i] = dynlocglob+i;
  ierr = SetMatrixValues(ICp,4,rows,1,&VAloc,x_VA);CHKERRQ(ierr);
  ierr = SetMatrixValues(ICp,4,rows,1,&VMloc,x_VM);CHKERRQ(ierr);
  ierr = SetMatrixValues(ICp,4,rows,1,&PGloc,x_PG);CHKERRQ(ierr);
  ierr = SetMatrixValues(ICp,4,rows,1,&QGloc,x_QG);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNExcModelGetEquationTypes_IEEET1 - Gets the number and indices of differential and algebraic equations
                                                   
  Input Parameters:
. dynexc - the DYNExcModel object

  Output Parameters:
+ ndiff - number of differential equations
. nalg  - number of algebraic equations
- eqtype - an array of size nvar that has
            eqtype[i] = DIFF_EQ if equation i is a differential
            vartype[i] = ALG_VAR if equation i is an algebraic

  NOTES: The user does not need to create the array vartype
*/
PetscErrorCode DYNExcModelGetEquationTypes_IEEET1(DYNExcModel dynexc,PetscInt *ndiff,PetscInt *nalg, PetscInt *eqtype)
{
  PetscErrorCode ierr;
  DYNIEEET1 ieeet1;

  PetscFunctionBegin;
  ierr = DYNExcModelGetModelData(dynexc,(void**)&ieeet1);CHKERRQ(ierr);
  if(ieeet1->TR == 0) {
    *ndiff = 3;
    *nalg  = 1;
    eqtype[1] = eqtype[2] = eqtype[3] = DIFF_EQ;
    eqtype[0] = ALG_EQ;

    ieeet1->TR = 1.0;
  } else {
  *ndiff = 4;
  *nalg  = 0;
  eqtype[0] = eqtype[1] = eqtype[2] = eqtype[3] = DIFF_EQ;
  }

  PetscFunctionReturn(0);
}

/*
  DYNExcModelDAERHSJacobianP_IEEET1 - Computes the RHS Jacobian of the exciter equations w.r.t. parameters

  Input Parameters:
+ dynexc     - the dynamic exciter model
. t          - the current time
. xdynexc    - exciter variables
. dynlocglob - starting location of exciter variables in the global vector 
. Valoc      - location of Va (bus voltage angle) in the parameter vector
- Vmloc      - location of Vm (bus voltage magnitude) in the parameter vector

  Output Parameters:
. jacP       - matrix of partial derivatives of DAE equations w.r.t parameter Va,Vm

*/ 
PetscErrorCode DYNExcModelDAERHSJacobianP_IEEET1(DYNExcModel dynexc,PetscReal t,const PetscScalar *xdynexc,Mat jacP,PetscInt dynlocglob,PetscInt Valoc,PetscInt Vmloc)
{
  PetscErrorCode ierr;
  DYNIEEET1      ieeet1;
  PetscInt       loc;
  PetscInt       row,col;
  PetscScalar    value;

  PetscFunctionBegin;

  ierr = DYNExcModelGetModelData(dynexc,(void**)&ieeet1);CHKERRQ(ierr);
  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);

  row = dynlocglob+1;
  col = Vmloc;
  if(!(ieeet1->VRatmax || ieeet1->VRatmin)) {
    value = ieeet1->KA/ieeet1->TA;
    ierr = SetMatrixValues(jacP,1,&row,1,&col,&value);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  DYNExcModelDAEFWDRHSJacobianP_IEEET1 - Computes the RHS Jacobian of the exciter equations w.r.t. parameters

  Input Parameters:
+ dynexc     - the dynamic exciter model
. t          - the current time
. xdynexc    - exciter variables
. dynlocglob - starting location of exciter variables in the global vector
. Valoc      - location of Va (bus voltage angle) in the parameter vector
- Vmloc      - location of Vm (bus voltage magnitude) in the parameter vector

  Output Parameters:
. jacP       - matrix of partial derivatives of DAE equations w.r.t parameter Va,Vm

*/
PetscErrorCode DYNExcModelDAEFWDRHSJacobianP_IEEET1(DYNExcModel dynexc,PetscReal t,const PetscScalar *xdynexc,Vec *jacP,PetscInt dynlocglob,PetscInt Valoc,PetscInt Vmloc)
{
  PetscErrorCode ierr;
  DYNIEEET1      ieeet1;
  PetscInt       loc;
  PetscInt       row;
  PetscScalar    value;

  PetscFunctionBegin;

  ierr = DYNExcModelGetModelData(dynexc,(void**)&ieeet1);CHKERRQ(ierr);
  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);

  row = dynlocglob+1;
  if(!(ieeet1->VRatmax || ieeet1->VRatmin)) {
    value = ieeet1->KA/ieeet1->TA;
    ierr = VecSetValue(jacP[Vmloc],row,value,INSERT_VALUES);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  DYNExcModelDAERHSJacobian_IEEET1 - Computes the RHS Jacobian of the IEEET1 DAE equations

  Input Parameters:
+ dynexc     - the dynamic exciter model
. t          - the current time
. VD         - real-part of complex bus voltage
. VQ         - imaginary-part of complex bus voltage
. x          - variables for the bus on which the generator is incident
. dynlocglob - starting location of exciter variables in the global vector 
- V_loc      - global location of VD and VQ V_loc[0] = VD_loc V_loc[1] = VQ_loc


  Output Parameters:
. J          - the Jacobian matrix

  Notes:
   For each bus, the network equations are ordered as [ID_loc;IQ_loc]

*/
PetscErrorCode DYNExcModelDAERHSJacobian_IEEET1(DYNExcModel dynexc,Mat J,PetscReal t,PetscScalar VD,PetscScalar VQ,PetscScalar *x,PetscInt dynlocglob,PetscInt V_loc[])
{
  PetscErrorCode ierr;
  DYNIEEET1      ieeet1;
  PetscInt       loc,i;
  PetscScalar    Efd,Vm,SE,dSE_dEfd,dVm_dVD,dVm_dVQ;
  PetscInt       row[2],col[10];
  PetscScalar    val[10];
  DYNStabModel   dynstab;

  PetscFunctionBegin;
  ierr = DYNExcModelGetModelData(dynexc,(void**)&ieeet1);CHKERRQ(ierr);
  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);

  Vm = PetscSqrtScalar(VD*VD + VQ*VQ);

  Efd = x[loc+2];

  SE = 0.0;
  dSE_dEfd = 0.0;
  if (Efd > ieeet1->Efdthresh) {
    SE = ieeet1->satB*PetscPowScalar((Efd - ieeet1->satA),2)/Efd;
    dSE_dEfd = ieeet1->satB*(2*(Efd - ieeet1->satA)/Efd - (PetscPowScalar((Efd - ieeet1->satA),2)/(Efd*Efd)));
  }

  dVm_dVD = VD/Vm;
  dVm_dVQ = VQ/Vm;

  row[0] = dynlocglob;
  col[0] = dynlocglob; col[1] = V_loc[0]; col[2] = V_loc[1];
  val[0] = -1/ieeet1->TR;
  val[1] = dVm_dVD/ieeet1->TR;
  val[2] = dVm_dVQ/ieeet1->TR;

  ierr = SetMatrixValues(J,1,row,3,col,val);CHKERRQ(ierr);

  row[0] = dynlocglob+1;
  if(ieeet1->VRatmax || ieeet1->VRatmin) {
    val[0] = 1.0;
    col[0] = dynlocglob+1;
    ierr = SetMatrixValues(J,1,row,1,col,val);CHKERRQ(ierr);
  } else {
    col[0] = dynlocglob; col[1] = dynlocglob+1; col[2] = dynlocglob+2; col[3] = dynlocglob+3;
    val[0] = ieeet1->KA*-1/ieeet1->TA;
    val[1] = -1/ieeet1->TA;
    val[2] = -ieeet1->KA*ieeet1->KF/(ieeet1->TF*ieeet1->TA);
    val[3] = ieeet1->KA/ieeet1->TA;

    ierr = SetMatrixValues(J,1,row,4,col,val);CHKERRQ(ierr);

    ierr = DYNExcModelGetDynStab(dynexc,&dynstab);CHKERRQ(ierr);
    if(dynstab) {
      PetscInt dVOTHSGdXdyn_num,dVOTHSGdXdyn_loc[6];
      PetscScalar dVOTHSGdXdyn[6],VOTHSG;
      ierr = DYNStabModelGetVOTHSG(dynstab,t,x,&VOTHSG,&dVOTHSGdXdyn_num,dVOTHSGdXdyn,dVOTHSGdXdyn_loc);CHKERRQ(ierr);

      if(dVOTHSGdXdyn_num) {
	for(i=0; i < dVOTHSGdXdyn_num; i++) dVOTHSGdXdyn[i] *= ieeet1->KA/ieeet1->TA;
      
	ierr = SetMatrixValues(J,1,row,dVOTHSGdXdyn_num,dVOTHSGdXdyn_loc,dVOTHSGdXdyn);CHKERRQ(ierr);
      }
    } 
  }

  row[0] = dynlocglob+2;
  col[0] = dynlocglob+1; col[1] = dynlocglob+2;
  val[0] = 1/ieeet1->TE;
  val[1] = -(ieeet1->KE + SE + dSE_dEfd*Efd)/ieeet1->TE;

  ierr = SetMatrixValues(J,1,row,2,col,val);CHKERRQ(ierr);

  row[0] = dynlocglob+3;
  col[0] = dynlocglob+2; col[1] = dynlocglob+3;
  val[0] = ieeet1->KF/(ieeet1->TF*ieeet1->TF);
  val[1] = -1/ieeet1->TF;

  ierr = SetMatrixValues(J,1,row,2,col,val);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNExcModelDAERHSFunction_IEEET1 - Computes the rhs of the DAE function for the IEEET1 model 
                                and also returns the exciter currents in network refernce frame

  Input Parameters:
+ dynexc - the dynamic exciter model
. t      - the current time
. x      - array of the variables for this bus
. VD     - real-part of the bus voltage
- VQ     - imaginary part of the bus voltage

  Output Parameters:
+ f      - array of rhs of DAE equations for the ieeet1 model 
. IGD    - real-part of exciter current in network reference frame
. IGQ    - imaginary-part of exciter current in network reference frame
*/
PetscErrorCode DYNExcModelDAERHSFunction_IEEET1(DYNExcModel dynexc,PetscReal t,PetscScalar VD, PetscScalar VQ, PetscScalar *x,PetscScalar *f)
{
  PetscErrorCode ierr;
  DYNIEEET1      ieeet1;
  PetscInt       loc;
  PetscScalar    Vm;
  PetscScalar    VT,VR,Efd,RF;
  PetscScalar    SE;
  DYNStabModel   dynstab;
  PetscScalar    VOTHSG=0.0;

  PetscFunctionBegin;
  /* Get the model data */
  ierr = DYNExcModelGetModelData(dynexc,(void**)&ieeet1);CHKERRQ(ierr);

  /* Get the location for the first exciter variable */
  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);

  VT  = x[loc];
  VR  = x[loc+1];
  Efd = x[loc+2];
  RF  = x[loc+3];

  ierr = DYNExcModelGetDynStab(dynexc,&dynstab);CHKERRQ(ierr);
  if(dynstab) {
    ierr = DYNStabModelGetVOTHSG(dynstab,t,x,&VOTHSG,NULL,NULL,NULL);CHKERRQ(ierr);
  }

  Vm = PetscSqrtScalar(VD*VD + VQ*VQ);
  SE = 0;
  if (Efd > ieeet1->Efdthresh) {
    SE = ieeet1->satB*PetscPowScalar((Efd - ieeet1->satA),2)/Efd;
  }
  f[loc]   = (-VT + Vm)/ieeet1->TR;
  if(ieeet1->VRatmax) f[loc+1] = VR - ieeet1->VRmax;
  else if(ieeet1->VRatmin) f[loc+1] = VR - ieeet1->VRmin;
  else f[loc+1] = (-VR + ieeet1->KA*RF - ieeet1->KA*ieeet1->KF/ieeet1->TF*Efd + ieeet1->KA*(ieeet1->Vref + VOTHSG - VT))/ieeet1->TA;
  f[loc+2] = (-(ieeet1->KE + SE)*Efd + VR)/ieeet1->TE;
  f[loc+3] = (-RF + ieeet1->KF/ieeet1->TF*Efd)/ieeet1->TF;

  PetscFunctionReturn(0);
}

/*
  DYNExcModelCreate_IEEET1 - Class constructor for IEEET1 exciter model
*/
PetscErrorCode DYNExcModelCreate_IEEET1(DYNExcModel dynexc)
{
  DYNIEEET1 ieeet1;
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = PetscCalloc1(1,&ieeet1);CHKERRQ(ierr);

  ieeet1->VRatmax = ieeet1->VRatmin = 0;

  dynexc->data = (void*)ieeet1;

  /* Inherit the ops */
  dynexc->ops.readdata             = DYNExcModelReadData_IEEET1;
  dynexc->ops.destroy              = DYNExcModelDestroy_IEEET1;
  dynexc->ops.getnvar              = DYNExcModelGetNvar_IEEET1;
  dynexc->ops.getsizeof            = DYNExcModelGetSizeof_IEEET1;
  dynexc->ops.setinitialconditions = DYNExcModelSetInitialConditions_IEEET1;
  dynexc->ops.setinitialconditionsp = DYNExcModelSetInitialConditionsP_IEEET1;
  dynexc->ops.daerhsfunction       = DYNExcModelDAERHSFunction_IEEET1;
  dynexc->ops.daerhsjacobian       = DYNExcModelDAERHSJacobian_IEEET1;
  dynexc->ops.daerhsjacobianp      = DYNExcModelDAERHSJacobianP_IEEET1;
  dynexc->ops.daefwdrhsjacobianp   = DYNExcModelDAEFWDRHSJacobianP_IEEET1;
  dynexc->ops.getequationtypes     = DYNExcModelGetEquationTypes_IEEET1;
  dynexc->ops.getbusnumid          = DYNExcModelGetBusnumID_IEEET1;
  dynexc->ops.getfieldvoltage      = DYNExcModelGetFieldVoltage_IEEET1;
  dynexc->ops.setevent             = DYNExcModelSetEvent_IEEET1;
  PetscFunctionReturn(0);
}


