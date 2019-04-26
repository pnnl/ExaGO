#include "dynsexs.h"
#include <private/dynexcmodelsimpl.h>

/*
  DYNEventMonitor_SEXS - Event monitoring routine for this exciter model. Value of fval is checked outside the function to see if an event is triggered
*/
PetscErrorCode DYNEventMonitor_SEXS(DYNExcModel dynexc, PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscScalar *fval)
{
  PetscErrorCode ierr;
  PetscInt loc;
  DYNSEXS      sexs;
  PetscScalar x1,y1,Efd,Ec;

  DYNStabModel   dynstab;
  PetscScalar    VOTHSG=0.0;

  PetscFunctionBegin;

  /* Get the model data  */
  ierr = DYNExcModelGetModelData(dynexc,(void**)&sexs);CHKERRQ(ierr);

  /* Get the location for the first exciter variable */
  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);

  x1  =x[loc];
  Efd  = x[loc+1];

  ierr = DYNExcModelGetDynStab(dynexc,&dynstab);CHKERRQ(ierr);
  if(dynstab) {
    ierr = DYNStabModelGetVOTHSG(dynstab,t,x,&VOTHSG,NULL,NULL,NULL);CHKERRQ(ierr);
  }

  Ec = PetscSqrtScalar(VD*VD + VQ*VQ); /* Terminal voltage */
  
  y1 = x1 + sexs->TA_TB*(sexs->Vref + VOTHSG - Ec); /* Output of the lead lag block */

  if(!sexs->Eatmax) fval[0] = sexs->Emax - Efd;
  else fval[0] = (-Efd + sexs->K*y1)/sexs->TE;

  if(!sexs->Eatmin) fval[1] = sexs->Emin - Efd;
  else fval[1] =  (-Efd + sexs->K*y1)/sexs->TE;

  PetscFunctionReturn(0);
}

/*
  DYNEventPostDGamma_SEXS - Compute derivative of gamma
*/
PetscErrorCode DYNEventPostDGamma_SEXS(DYNExcModel dynexc, PetscInt nevents, PetscInt ev_list[], PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscScalar *dgdx, PetscScalar *dgdp,PetscInt Valoc, PetscInt Vmloc, PetscInt Pgloc, PetscInt Qgloc)
{
  PetscErrorCode ierr;
  DYNSEXS      sexs;
  PetscInt     loc;
  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNExcModelGetModelData(dynexc,(void**)&sexs);CHKERRQ(ierr);

  /* Get the location for the first exciter variable  */
  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);

  SETERRQ(PETSC_COMM_SELF,0,"Post Event Gamma function for exciter model SEXS not implemented yet");

  PetscFunctionReturn(0);
}

/*
  DYNEventPostFunction_SEXS - Post event routine for this exciter model
*/
PetscErrorCode DYNEventPostFunction_SEXS(DYNExcModel dynexc, PetscInt nevents, PetscInt ev_list[], PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscBool forward_solve, PetscBool *solve_algebraic)
{
  PetscErrorCode ierr;
  DYNSEXS      sexs;
  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNExcModelGetModelData(dynexc,(void**)&sexs);CHKERRQ(ierr);

  if(ev_list[0] == 0) { /* Max E */
    if(!sexs->Eatmax) {
      sexs->Eatmax = 1;
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%4.3f s] Bus %d: Exciter Model SEXS for generator %s hit its max. E limit = %4.3f\n",t,sexs->bus_i,sexs->id,sexs->Emax);CHKERRQ(ierr);
    } else {
      sexs->Eatmax = 0;
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%4.3f s] Bus %d: Exciter Model SEXS for generator %s released from max. E limit = %4.3f\n",t,sexs->bus_i,sexs->id,sexs->Emax);CHKERRQ(ierr);
    }
  } else { /* Min E */
    if(!sexs->Emin) {
      sexs->Eatmin = 1;
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%4.3f s] Bus %d: Exciter Model SEXS for generator %s hit its min. E limit = %4.3f\n",t,sexs->bus_i,sexs->id,sexs->Emin);CHKERRQ(ierr);
    } else {
      sexs->Eatmin = 0;
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%4.3f s] Bus %d: Exciter Model SEXS for generator %s released from min. E limit = %4.3f\n",t,sexs->bus_i,sexs->id,sexs->Emin);CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(0);
}


/*
  DYNExcModelSetEvent_SEXS - Sets the event info for the SEXS exciter model

  Input Parameters
. dynexc - the dynamic exciter model

  Output Parameters
+ nmons - the number of event monitors
. direction - the event directions
. terminate - flags for termination when the event is located
. eventfcn  - the function describing the event condition
- posteventfcn - An optional function that gets called when the event is located
*/
PetscErrorCode DYNExcModelSetEvent_SEXS(DYNExcModel dynexc,PetscInt *nmons,PetscInt *direction,PetscBool *terminate,PetscErrorCode (**eventfcn)(DYNExcModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*), PetscErrorCode (**posteventfcn)(DYNExcModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscBool,PetscBool*), PetscErrorCode (**posteventdgamma)(DYNExcModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscInt))
{

  PetscFunctionBegin;
  *nmons = 2;
  direction[0] = -1; direction[1] = 1;
  terminate[0] = terminate[1] = PETSC_FALSE;
  *eventfcn = DYNEventMonitor_SEXS;
  *posteventfcn = DYNEventPostFunction_SEXS;
  *posteventdgamma = DYNEventPostDGamma_SEXS;
  PetscFunctionReturn(0);
}


/*
  DYNExcModelReadData_SEXS - Parses the line containining the SEXS data and
  populates it in the SEXS data structure

  Input Parameters:
+ dynexc - the dynamic exciter struct
- line - The line in the dyr file

*/
PetscErrorCode DYNExcModelReadData_SEXS(DYNExcModel dynexc,char* line)
{
  PetscErrorCode ierr;
  DYNSEXS sexs;

  PetscFunctionBegin;
  ierr = DYNExcModelGetModelData(dynexc,(void**)&sexs);CHKERRQ(ierr);

  /* Read the SEXS data */ 
  //Data format: IBUS, ’SEXS’, I, TA/TB, TB, K, TE, EMIN, EMAX

  sscanf(line,"%d, 'SEXS',%[^,], %lf, %lf, %lf, %lf, %lf, %lf",\
	 &sexs->bus_i,sexs->id,&sexs->TA_TB,&sexs->TB,&sexs->K,&sexs->TE, \
	 &sexs->Emin,&sexs->Emax);

  PetscFunctionReturn(0);
}

/*
  DYNExcModelGetBusnumID_SEXS - Returns the bus number and ID associated with this exciter model

  Input Parameters:
. dynexc - the dynamic exciter model object

  Output Parameters:
+ busnum - the bus number 
- excid  - the exciter ID
*/
PetscErrorCode DYNExcModelGetBusnumID_SEXS(DYNExcModel dynexc,PetscInt *busnum,char **excid)
{
  PetscErrorCode ierr;
  DYNSEXS      sexs;

  PetscFunctionBegin;
  ierr = DYNExcModelGetModelData(dynexc,(void**)&sexs);CHKERRQ(ierr);
  *busnum = sexs->bus_i;
  *excid  = sexs->id;
  PetscFunctionReturn(0);
}

/*
  DYNExcModelDestroy_SEXS - Destroys the SEXS object data

  Input Parameters
. DYNExcModel - the DYNExcModel object

  Notes:
  Called when DYNExcModelDestroy() is called
*/
PetscErrorCode DYNExcModelDestroy_SEXS(DYNExcModel dynexc)
{
  PetscErrorCode ierr;
  DYNSEXS      sexs;

  PetscFunctionBegin;
  ierr = DYNExcModelGetModelData(dynexc,(void**)&sexs);CHKERRQ(ierr);
  ierr = PetscFree(sexs);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNExcModelGetNvar_SEXS - Returns the number of variables for the SEXS model

  Input parameters
. dynexc - the dynamic exciter object

  Output Parameters
. nvar - number of variables for this model
*/
PetscErrorCode DYNExcModelGetNvar_SEXS(DYNExcModel dynexc,PetscInt *nvar)
{
  PetscFunctionBegin;
  *nvar = DYNSEXS_nvar;
  PetscFunctionReturn(0);
}

/*
  DYNExcModelGetSizeof_SEXS - Returns the size of SEXS struct

  Input Parameters:
. dynexc - the dynamic exciter model
 
  Output Parameters:
. size - size of the SEXS object obtained from sizeof()
*/
PetscErrorCode DYNExcModelGetSizeof_SEXS(DYNExcModel dynexc,PetscInt *size)
{
  PetscFunctionBegin;
  *size = sizeof(struct _p_DYNSEXS);
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
PetscErrorCode DYNExcModelGetFieldVoltage_SEXS(DYNExcModel dynexc,PetscReal t, PetscScalar *xdyn,PetscScalar *Efd,PetscInt *dEfddXexc_num,PetscScalar dEfddXexc[],PetscInt dEfddXexc_loc[])
{
  PetscErrorCode ierr;
  DYNSEXS      sexs;
  PetscInt       loc;

  PetscFunctionBegin;
  /* Get Model data */
  ierr = DYNExcModelGetModelData(dynexc,(void**)&sexs);CHKERRQ(ierr);

  /* Get the location for the first variable for this exciter model */
  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);
  *Efd = xdyn[loc+1];

  if(dEfddXexc_num) *dEfddXexc_num   = 1;
  if(dEfddXexc) {
    dEfddXexc[0]     = 1.0;
    if(sexs->Eatmax || sexs->Eatmin) dEfddXexc[0] = 0.0;
  }
  if(dEfddXexc_loc)  dEfddXexc_loc[0] = loc+1;

  PetscFunctionReturn(0);
}

/*
  DYNExcModelSetInitialConditions_SEXS - Sets the initial conditions (x(t0)) for the SEXS model

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
PetscErrorCode DYNExcModelSetInitialConditions_SEXS(DYNExcModel dynexc,PetscScalar VD, PetscScalar VQ, PetscScalar *x)
{
  PetscErrorCode ierr;
  DYNSEXS sexs;
  DYNGenModel dyngen;
  PetscScalar Vm;
  PetscScalar  x1,Efd;
  PetscInt    loc;

  PetscFunctionBegin;
  ierr = DYNExcModelGetModelData(dynexc,(void**)&sexs);CHKERRQ(ierr);

  ierr = DYNExcModelGetDynGen(dynexc,&dyngen);CHKERRQ(ierr);

  /* Voltage magnitude */
  Vm = PetscSqrtScalar(VD*VD + VQ*VQ);

  /* Get the field voltage Efd0 from the generator */
  ierr = DYNGenModelGetInitialFieldVoltage(dyngen,&Efd);CHKERRQ(ierr);

  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);

  /* Set the reference voltage */
  sexs->Vref = Efd/sexs->K + Vm;
  x1 = (1 - sexs->TA_TB)*(sexs->Vref - Vm);

  x[loc]    = x1;
  x[loc+1]  = Efd;

  /* Set the initial flags */
  sexs->Eatmax= sexs->Eatmin = 0;

  PetscFunctionReturn(0);
}

/*
  DYNExcModelGetEquationTypes_SEXS - Gets the number and indices of differential and algebraic equations
                                                   
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
PetscErrorCode DYNExcModelGetEquationTypes_SEXS(DYNExcModel dynexc,PetscInt *ndiff,PetscInt *nalg, PetscInt *eqtype)
{
  PetscErrorCode ierr;
  DYNSEXS sexs;

  PetscFunctionBegin;
  ierr = DYNExcModelGetModelData(dynexc,(void**)&sexs);CHKERRQ(ierr);

  *ndiff = 2;
  *nalg  = 0;
  eqtype[0] = eqtype[1] = DIFF_EQ;

  PetscFunctionReturn(0);
}

/*
  DYNExcModelDAERHSJacobianP_SEXS - Computes the RHS Jacobian of the exciter equations w.r.t. parameters

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
PetscErrorCode DYNExcModelDAERHSJacobianP_SEXS(DYNExcModel dynexc,PetscReal t,const PetscScalar *xdynexc,Mat jacP,PetscInt dynlocglob,PetscInt Valoc,PetscInt Vmloc)
{
  PetscErrorCode ierr;
  DYNSEXS      sexs;
  PetscInt     loc;

  PetscFunctionBegin;

  ierr = DYNExcModelGetModelData(dynexc,(void**)&sexs);CHKERRQ(ierr);

  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);

  SETERRQ(PETSC_COMM_SELF,0,"No RHSJacobianP implemented for exciter model SEXS");
  PetscFunctionReturn(0);
}

/*
  DYNExcModelDAEFWDRHSJacobianP_SEXS - Computes the RHS Jacobian of the exciter equations w.r.t. parameters

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
PetscErrorCode DYNExcModelDAEFWDRHSJacobianP_SEXS(DYNExcModel dynexc,PetscReal t,const PetscScalar *xdynexc,Vec *jacP,PetscInt dynlocglob,PetscInt Valoc,PetscInt Vmloc)
{
  PetscErrorCode ierr;
  DYNSEXS      sexs;
  PetscInt     loc;

  PetscFunctionBegin;

  ierr = DYNExcModelGetModelData(dynexc,(void**)&sexs);CHKERRQ(ierr);

  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);

  SETERRQ(PETSC_COMM_SELF,0,"No DAEFWDRHSJacobianP implemented for exciter model SEXS");

  PetscFunctionReturn(0);
}

/*
  DYNExcModelDAERHSJacobian_SEXS - Computes the RHS Jacobian of the SEXS DAE equations

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

*/
PetscErrorCode DYNExcModelDAERHSJacobian_SEXS(DYNExcModel dynexc,Mat J,PetscReal t,PetscScalar VD,PetscScalar VQ,PetscScalar *x,PetscInt dynlocglob,PetscInt V_loc[])
{
  PetscErrorCode ierr;
  DYNSEXS      sexs;
  PetscInt       loc,i;
  PetscScalar    Vm,dVm_dVD,dVm_dVQ;
  PetscInt       row[2],col[10];
  PetscScalar    val[10];
  DYNStabModel   dynstab;

  PetscFunctionBegin;
  ierr = DYNExcModelGetModelData(dynexc,(void**)&sexs);CHKERRQ(ierr);

  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);

  Vm = PetscSqrtScalar(VD*VD + VQ*VQ);

  dVm_dVD = VD/Vm;
  dVm_dVQ = VQ/Vm;

  row[0] = dynlocglob;
  col[0] = dynlocglob; col[1] = V_loc[0]; col[2] = V_loc[1];
  val[0] = -1/sexs->TB;
  val[1] = -1*(1-sexs->TA_TB)*dVm_dVD/sexs->TB; 
  val[2] = -1*(1-sexs->TA_TB)*dVm_dVQ/sexs->TB;
  ierr = SetMatrixValues(J,1,row,3,col,val);CHKERRQ(ierr);
  ierr = DYNExcModelGetDynStab(dynexc,&dynstab);CHKERRQ(ierr);
  if(dynstab) {
    PetscInt dVOTHSGdXdyn_num,dVOTHSGdXdyn_loc[6];
    PetscScalar dVOTHSGdXdyn[6],VOTHSG;
    ierr = DYNStabModelGetVOTHSG(dynstab,t,x,&VOTHSG,&dVOTHSGdXdyn_num,dVOTHSGdXdyn,dVOTHSGdXdyn_loc);CHKERRQ(ierr);
    if(dVOTHSGdXdyn_num) {
      for(i=0; i < dVOTHSGdXdyn_num; i++) dVOTHSGdXdyn[i] *= (1-sexs->TA_TB)/sexs->TB;
      ierr = SetMatrixValues(J,1,row,dVOTHSGdXdyn_num,dVOTHSGdXdyn_loc,dVOTHSGdXdyn);CHKERRQ(ierr);
    }
  }

  row[0] = dynlocglob+1;
  if(sexs->Eatmax || sexs->Eatmin) {
    val[0] = 1.0;
    col[0] = dynlocglob+1;
    ierr = SetMatrixValues(J,1,row,1,col,val);CHKERRQ(ierr);
  } else {
    col[0] = dynlocglob;col[1]=dynlocglob+1;col[2]=V_loc[0];col[3]=V_loc[1];
    val[0] = sexs->K/sexs->TE;
    val[1] = -1/sexs->TE;
    val[2] = -dVm_dVD*sexs->K/sexs->TE*sexs->TA_TB;
    val[3] = -dVm_dVQ*sexs->K/sexs->TE*sexs->TA_TB;
    ierr = SetMatrixValues(J,1,row,4,col,val);CHKERRQ(ierr);
    if(dynstab) {
      PetscInt dVOTHSGdXdyn_num,dVOTHSGdXdyn_loc[6];
      PetscScalar dVOTHSGdXdyn[6],VOTHSG;
      ierr = DYNStabModelGetVOTHSG(dynstab,t,x,&VOTHSG,&dVOTHSGdXdyn_num,dVOTHSGdXdyn,dVOTHSGdXdyn_loc);CHKERRQ(ierr);
      if(dVOTHSGdXdyn_num) {
	for(i=0; i < dVOTHSGdXdyn_num; i++) dVOTHSGdXdyn[i] *= sexs->K/sexs->TE*sexs->TA_TB;
	ierr = SetMatrixValues(J,1,row,dVOTHSGdXdyn_num,dVOTHSGdXdyn_loc,dVOTHSGdXdyn);CHKERRQ(ierr);
      }
    }
  }

  PetscFunctionReturn(0);
}

/*
  DYNExcModelDAERHSFunction_SEXS - Computes the rhs of the DAE function for the SEXS model 
  
  Input Parameters:
+ dynexc - the dynamic exciter model
. t      - the current time
. x      - array of the variables for this bus
. VD     - real-part of the bus voltage
- VQ     - imaginary part of the bus voltage

  Output Parameters:
+ f      - array of rhs of DAE equations for the sexs model 
*/
PetscErrorCode DYNExcModelDAERHSFunction_SEXS(DYNExcModel dynexc,PetscReal t,PetscScalar VD, PetscScalar VQ, PetscScalar *x,PetscScalar *f)
{
  PetscErrorCode ierr;
  DYNSEXS      sexs;
  PetscInt       loc;
  PetscScalar    Vm;
  PetscScalar    x1,Efd,y1;

  DYNStabModel   dynstab;
  PetscScalar    VOTHSG=0.0;
  PetscFunctionBegin;
  /* Get the model data */
  ierr = DYNExcModelGetModelData(dynexc,(void**)&sexs);CHKERRQ(ierr);

  /* Get the location for the first exciter variable */
  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);

  x1  = x[loc];
  Efd  = x[loc+1];

  ierr = DYNExcModelGetDynStab(dynexc,&dynstab);CHKERRQ(ierr);
  if(dynstab) {
    ierr = DYNStabModelGetVOTHSG(dynstab,t,x,&VOTHSG,NULL,NULL,NULL);CHKERRQ(ierr);
  }

  Vm = PetscSqrtScalar(VD*VD + VQ*VQ);
  y1 = x1 + sexs->TA_TB*(sexs->Vref + VOTHSG - Vm);

  f[loc] = (-x1 + (1-sexs->TA_TB)*(sexs->Vref + VOTHSG - Vm))/sexs->TB;
  if(sexs->Eatmax) f[loc+1] = Efd - sexs->Emax;
  else if(sexs->Eatmin) f[loc+1] = Efd - sexs->Emin;
  else f[loc+1] = (-Efd + sexs->K*y1)/sexs->TE;

  PetscFunctionReturn(0);
}

/*
  DYNExcModelCreate_SEXS - Class constructor for SEXS exciter model
*/
PetscErrorCode DYNExcModelCreate_SEXS(DYNExcModel dynexc)
{
  DYNSEXS sexs;
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = PetscCalloc1(1,&sexs);CHKERRQ(ierr);
  sexs->Eatmax = sexs->Eatmin = 0;

  dynexc->data = (void*)sexs;

  /* Inherit the ops */
  dynexc->ops.readdata             = DYNExcModelReadData_SEXS;
  dynexc->ops.destroy              = DYNExcModelDestroy_SEXS;
  dynexc->ops.getnvar              = DYNExcModelGetNvar_SEXS;
  dynexc->ops.getsizeof            = DYNExcModelGetSizeof_SEXS;
  dynexc->ops.setinitialconditions = DYNExcModelSetInitialConditions_SEXS;
  dynexc->ops.daerhsfunction       = DYNExcModelDAERHSFunction_SEXS;
  dynexc->ops.daerhsjacobian       = DYNExcModelDAERHSJacobian_SEXS;
  dynexc->ops.daerhsjacobianp      = DYNExcModelDAERHSJacobianP_SEXS;
  dynexc->ops.daefwdrhsjacobianp   = DYNExcModelDAEFWDRHSJacobianP_SEXS;
  dynexc->ops.getequationtypes     = DYNExcModelGetEquationTypes_SEXS;
  dynexc->ops.getbusnumid          = DYNExcModelGetBusnumID_SEXS;
  dynexc->ops.getfieldvoltage      = DYNExcModelGetFieldVoltage_SEXS;
  dynexc->ops.setevent             = DYNExcModelSetEvent_SEXS;
  PetscFunctionReturn(0);
}


