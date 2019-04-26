#include "dynexst1.h"
#include <private/dynexcmodelsimpl.h>

/*
  DYNEventMonitor_EXST1 - Event monitoring routine for this exciter model. Value of fval is checked outside the function to see if an event is triggered
*/
PetscErrorCode DYNEventMonitor_EXST1(DYNExcModel dynexc, PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscScalar *fval)
{
  PetscErrorCode ierr;
  PetscInt loc;
  DYNEXST1      exst1;
  DYNGenModel dyngen;
  PetscScalar Ifd,Efd,VT,VF;

  DYNStabModel   dynstab;
  PetscScalar    VOTHSG=0.0;

  PetscFunctionBegin;

  // Get the model data 
  ierr = DYNExcModelGetModelData(dynexc,(void**)&exst1);CHKERRQ(ierr);

  // Name of elements of exst1: Vimin; Vimax; VRmax; VRmin; VIatmax;VIatmin;Efdatmax;Efdatmin;Efdmin;Efdmax;
  // Get the location for the first exciter variable
  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);

  Efd  =x[loc];
  VT  = x[loc+1];
  VF  = x[loc+3];

  ierr = DYNExcModelGetDynStab(dynexc,&dynstab);CHKERRQ(ierr);
  if(dynstab) {
    ierr = DYNStabModelGetVOTHSG(dynstab,t,x,&VOTHSG,NULL,NULL,NULL);CHKERRQ(ierr);
  }

  fval[0] = exst1->Vimax - (VOTHSG+exst1->Vref-VT-VF); // The signal Vi is the output of the adder block given by Vi=Vref-VT+VS-VF
  fval[1] = exst1->Vimin - (VOTHSG+exst1->Vref-VT-VF);

  //Get Ifd 
  ierr = DYNExcModelGetDynGen(dynexc,&dyngen);CHKERRQ(ierr);
  ierr = DYNGenModelGetFieldCurrent(dyngen,x,&Ifd);CHKERRQ(ierr); 
  // x was passed as a pointer to DYNEventMonitor_EXST1, hence there is no &.
  exst1->Efdmax=VT*exst1->VRmax-exst1->KC*Ifd;
  exst1->Efdmin=VT*exst1->VRmin-exst1->KC*Ifd;
  fval[2]=exst1->Efdmax-Efd;  
  fval[3]=exst1->Efdmin-Efd; 

  PetscFunctionReturn(0);
}

/*
  DYNEventPostDGamma_EXST1 - Compute derivative of gamma
*/
PetscErrorCode DYNEventPostDGamma_EXST1(DYNExcModel dynexc, PetscInt nevents, PetscInt ev_list[], PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscScalar *dgdx, PetscScalar *dgdp,PetscInt Valoc, PetscInt Vmloc, PetscInt Pgloc, PetscInt Qgloc)
{
  PetscErrorCode ierr;
  DYNEXST1      exst1;
  //  DYNStabModel   dynstab;
  //PetscScalar    Vm,dVm_dVD,dVm_dVQ;
  //  PetscInt       loc;
  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNExcModelGetModelData(dynexc,(void**)&exst1);CHKERRQ(ierr);
  /*
  // Get the location for the first exciter variable 
  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);

  for(i=0; i < nevents; i++) {
  if(ev_list[i] == 0) { // Max VI 
    if(exst1->VIatmax) {
      dgdx[loc+1]  = 1;
      dgdx[loc+3]  = 1;
      ierr = DYNExcModelGetDynStab(dynexc,&dynstab);CHKERRQ(ierr);
      if(dynstab) {
         PetscInt dVOTHSGdXdyn_num,dVOTHSGdXdyn_loc[6];
         PetscScalar dVOTHSGdXdyn[6],VOTHSG;
         ierr = DYNStabModelGetVOTHSG(dynstab,t,x,&VOTHSG,&dVOTHSGdXdyn_num,dVOTHSGdXdyn,dVOTHSGdXdyn_loc);CHKERRQ(ierr);
         if(dVOTHSGdXdyn_num) {
           ierr = VecSetValues(dgdx,dVOTHSGdXdyn_num,dVOTHSGdXdyn_loc,-1,INSERT_VALUES);CHKERRQ(ierr);
          }
       }  
      dgdp[Vmloc] = -1;
    }
  }else if(ev_list[i] == 1) {
   if(exst1->VIatmin){
      dgdx[loc+1]  = 1;
      dgdx[loc+3]  = 1;
      ierr = DYNExcModelGetDynStab(dynexc,&dynstab);CHKERRQ(ierr);
         if(dynstab) {
         PetscInt dVOTHSGdXdyn_num,dVOTHSGdXdyn_loc[6];
         PetscScalar dVOTHSGdXdyn[6],VOTHSG;
         ierr = DYNStabModelGetVOTHSG(dynstab,t,x,&VOTHSG,&dVOTHSGdXdyn_num,dVOTHSGdXdyn,dVOTHSGdXdyn_loc);CHKERRQ(ierr);
         if(dVOTHSGdXdyn_num) {
            ierr = VecSetValues(dgdx,dVOTHSGdXdyn_num,dVOTHSGdXdyn_loc,-1,INSERT_VALUES);CHKERRQ(ierr);
         }
      }  
      dgdp[Vmloc] = -1;
    }
  }
else if(ev_list[i] == 2) {
   if(exst1->Efdatmax){
      dgdx[loc]  = -1;
      dgdx[loc+1]  = exst1->VRmax; 
// dgdx[Ifdxloc] = -exst1->KC*dIFddx;
    }
  }else if(ev_list[i] == 3) {
   if(exst1->Efdatmin){
      dgdx[loc]  = -1;
      dgdx[loc+1]  = exst1->VRmin; 
// dgdx[Ifdxloc] = -exst1->KC*dIFddx;
    }
  }
}
*/
  PetscFunctionReturn(0);
}

/*
  DYNEventPostFunction_EXST1 - Post event routine for this exciter model
*/
PetscErrorCode DYNEventPostFunction_EXST1(DYNExcModel dynexc, PetscInt nevents, PetscInt ev_list[], PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscBool forward_solve, PetscBool *solve_algebraic)
{
  PetscErrorCode ierr;
  DYNEXST1      exst1;
  PetscInt      i;
  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNExcModelGetModelData(dynexc,(void**)&exst1);CHKERRQ(ierr);

  for(i=0; i < nevents; i++) {
    if(ev_list[i] == 0) { // Max VI 
      if(!exst1->VIatmax) {
	exst1->VIatmax = 1;
      } else {
	exst1->VIatmax = 0;
      }
    } else if(ev_list[i] == 1) { // Min VI
      if(!exst1->VIatmin) {
	exst1->VIatmin = 1;
      } else { 
	exst1->VIatmin = 0;
      }
    }

    if(ev_list[i] == 2) { // Max Efd 
      if(!exst1->Efdatmax) {
	exst1->Efdatmax = 1;
      } else {
	exst1->Efdatmax = 0;
      }
    } else if(ev_list[i] == 3) { // Min Efd
      if(!exst1->Efdatmin) {
	exst1->Efdatmin = 1;
      } else {
	exst1->Efdatmin = 0;
      }
    }
  }
  PetscFunctionReturn(0);
}


/*
  DYNExcModelSetEvent_EXST1 - Sets the event info for the EXST1 exciter model

  Input Parameters
. dynexc - the dynamic exciter model

  Output Parameters
+ nmons - the number of event monitors
. direction - the event directions
. terminate - flags for termination when the event is located
. eventfcn  - the function describing the event condition
- posteventfcn - An optional function that gets called when the event is located
*/
PetscErrorCode DYNExcModelSetEvent_EXST1(DYNExcModel dynexc,PetscInt *nmons,PetscInt *direction,PetscBool *terminate,PetscErrorCode (**eventfcn)(DYNExcModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*), PetscErrorCode (**posteventfcn)(DYNExcModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscBool,PetscBool*), PetscErrorCode (**posteventdgamma)(DYNExcModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscInt))
{

  PetscFunctionBegin;
  *nmons = 4;
  direction[0] = direction[1] = direction[2] = direction[3] = 0;
  terminate[0] = terminate[1] = terminate[2] = terminate[3] = PETSC_FALSE;
  *eventfcn = DYNEventMonitor_EXST1;
  *posteventfcn = DYNEventPostFunction_EXST1;
  *posteventdgamma = DYNEventPostDGamma_EXST1;
  PetscFunctionReturn(0);
}


/*
  DYNExcModelReadData_EXST1 - Parses the line containining the EXST1 data and
  populates it in the EXST1 data structure

  Input Parameters:
+ dynexc - the dynamic exciter struct
- line - The line in the dyr file

*/
PetscErrorCode DYNExcModelReadData_EXST1(DYNExcModel dynexc,char* line)
{
  PetscErrorCode ierr;
  DYNEXST1 exst1;

  PetscFunctionBegin;
  ierr = DYNExcModelGetModelData(dynexc,(void**)&exst1);CHKERRQ(ierr);

  /* Read the EXST1 data */ 
  //Data format: IBUS, ’EXST1’, I, TR, VIMAX, VIMIN, TC, TB, KA, TA, VRMAX, VRMIN, KC, KF, TF/
  //1, ’EXST1’,1 , 0.02, 13, 0, 1, 123, 0.1, 300, 11.8, 0, 1, 0.1, 123,/

  sscanf(line,"%d, 'EXST1',%[^,], %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf",\
	 &exst1->bus_i,exst1->id,&exst1->TR,&exst1->Vimax,&exst1->Vimin,&exst1->TC, \
	 &exst1->TB,&exst1->KA,&exst1->TA,&exst1->VRmax,&exst1->VRmin,&exst1->KC,\
	 &exst1->KF,&exst1->TF);

  PetscFunctionReturn(0);
}

/*
  DYNExcModelGetBusnumID_EXST1 - Returns the bus number and ID associated with this exciter model

  Input Parameters:
. dynexc - the dynamic exciter model object

  Output Parameters:
+ busnum - the bus number 
- excid  - the exciter ID
*/
PetscErrorCode DYNExcModelGetBusnumID_EXST1(DYNExcModel dynexc,PetscInt *busnum,char **excid)
{
  PetscErrorCode ierr;
  DYNEXST1      exst1;

  PetscFunctionBegin;
  ierr = DYNExcModelGetModelData(dynexc,(void**)&exst1);CHKERRQ(ierr);
  *busnum = exst1->bus_i;
  *excid  = exst1->id;
  PetscFunctionReturn(0);
}

/*
  DYNExcModelDestroy_EXST1 - Destroys the EXST1 object data

  Input Parameters
. DYNExcModel - the DYNExcModel object

  Notes:
  Called when DYNExcModelDestroy() is called
*/
PetscErrorCode DYNExcModelDestroy_EXST1(DYNExcModel dynexc)
{
  PetscErrorCode ierr;
  DYNEXST1      exst1;

  PetscFunctionBegin;
  ierr = DYNExcModelGetModelData(dynexc,(void**)&exst1);CHKERRQ(ierr);
  ierr = PetscFree(exst1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNExcModelGetNvar_EXST1 - Returns the number of variables for the EXST1 model

  Input parameters
. dynexc - the dynamic exciter object

  Output Parameters
. nvar - number of variables for this model
*/
PetscErrorCode DYNExcModelGetNvar_EXST1(DYNExcModel dynexc,PetscInt *nvar)
{
  PetscFunctionBegin;
  *nvar = DYNEXST1_nvar;
  PetscFunctionReturn(0);
}

/*
  DYNExcModelGetSizeof_EXST1 - Returns the size of EXST1 struct

  Input Parameters:
. dynexc - the dynamic exciter model
 
  Output Parameters:
. size - size of the EXST1 object obtained from sizeof()
*/
PetscErrorCode DYNExcModelGetSizeof_EXST1(DYNExcModel dynexc,PetscInt *size)
{
  PetscFunctionBegin;
  *size = sizeof(struct _p_DYNEXST1);
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
PetscErrorCode DYNExcModelGetFieldVoltage_EXST1(DYNExcModel dynexc,PetscReal t, PetscScalar *xdyn,PetscScalar *Efd,PetscInt *dEfddXexc_num,PetscScalar dEfddXexc[],PetscInt dEfddXexc_loc[])
{
  PetscErrorCode ierr;
  DYNEXST1      exst1;
  PetscInt       loc;

  PetscFunctionBegin;
  /* Get Model data */
  ierr = DYNExcModelGetModelData(dynexc,(void**)&exst1);CHKERRQ(ierr);

  /* Get the location for the first variable for this exciter model */
  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);
  *Efd = xdyn[loc];

  if(exst1->Efdatmax) *Efd = exst1->Efdmax; // Update Efd if limit is being hit
  if(exst1->Efdatmin) *Efd = exst1->Efdmin; // Update Efd if limit is being hit


  if(dEfddXexc_num) *dEfddXexc_num   = 1;
  if(dEfddXexc) {
    dEfddXexc[0]     = 1.0;
    if(exst1->Efdatmax || exst1->Efdatmin) dEfddXexc[0] = 0.0;
  }
  if(dEfddXexc_loc)  dEfddXexc_loc[0] = loc;

  PetscFunctionReturn(0);
}

/*
  DYNExcModelSetInitialConditions_EXST1 - Sets the initial conditions (x(t0)) for the EXST1 model

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
PetscErrorCode DYNExcModelSetInitialConditions_EXST1(DYNExcModel dynexc,PetscScalar VD, PetscScalar VQ, PetscScalar *x)
{
  PetscErrorCode ierr;
  DYNEXST1 exst1;
  DYNGenModel dyngen;
  PetscScalar Vm;
  PetscScalar  Ifd,Efd,VT,VF,X; //Vll will be needed if using transfer-function implementation of lead-lag block
  PetscInt    loc;

  PetscFunctionBegin;
  ierr = DYNExcModelGetModelData(dynexc,(void**)&exst1);CHKERRQ(ierr);

  ierr = DYNExcModelGetDynGen(dynexc,&dyngen);CHKERRQ(ierr);

  /* Voltage magnitude */
  Vm = PetscSqrtScalar(VD*VD + VQ*VQ);

  /* Get the field voltage Efd0 from the generator */
  ierr = DYNGenModelGetInitialFieldVoltage(dyngen,&Efd);CHKERRQ(ierr);

  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);
  VT = Vm;
  X=(Efd/exst1->KA)*(1-(exst1->TC/exst1->TB));// Vll=Efd/exst1->KA; //for Transfer function implementation of lead-lag
  VF=0;

  x[loc]    = Efd;
  x[loc+1]  = VT;
  x[loc+2] = X;   //x[loc+2]  = Vll; //for  Transfer function implementation of lead-lag
  x[loc+3]  = VF;

  //Set the reference voltage:
  exst1->Vref=Efd/exst1->KA+Vm;  //AT SOME POINT, SUBTRACT VS HERE TO GET Vref.
  /* Set the initial flags */
  exst1->VIatmax= exst1->VIatmin = 0;
  exst1->Efdatmax= exst1->Efdatmin = 0;

  ierr = DYNExcModelGetDynGen(dynexc,&dyngen);CHKERRQ(ierr);
  ierr = DYNGenModelGetFieldCurrent(dyngen,x,&Ifd);CHKERRQ(ierr); 
  exst1->Efdmax=VT*exst1->VRmax-exst1->KC*Ifd;
  exst1->Efdmin=VT*exst1->VRmin-exst1->KC*Ifd;

  PetscFunctionReturn(0);
}

/*
  DYNExcModelGetEquationTypes_EXST1 - Gets the number and indices of differential and algebraic equations
                                                   
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
PetscErrorCode DYNExcModelGetEquationTypes_EXST1(DYNExcModel dynexc,PetscInt *ndiff,PetscInt *nalg, PetscInt *eqtype)
{
  PetscErrorCode ierr;
  DYNEXST1 exst1;

  PetscFunctionBegin;
  ierr = DYNExcModelGetModelData(dynexc,(void**)&exst1);CHKERRQ(ierr);

  if(exst1->TR == 0) {
    *ndiff = 3;
    *nalg  = 1;
    eqtype[0] = DIFF_EQ; 
    eqtype[1] = ALG_EQ; 
    eqtype[2] = eqtype[3] = DIFF_EQ; 

    exst1->TR = 1.0;
  } else {
  *ndiff = 4;
  *nalg  = 0;
  eqtype[0] = eqtype[1] = eqtype[2] = eqtype[3] = DIFF_EQ;
  }

  PetscFunctionReturn(0);
}

/*
  DYNExcModelDAERHSJacobianP_EXST1 - Computes the RHS Jacobian of the exciter equations w.r.t. parameters

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
PetscErrorCode DYNExcModelDAERHSJacobianP_EXST1(DYNExcModel dynexc,PetscReal t,const PetscScalar *xdynexc,Mat jacP,PetscInt dynlocglob,PetscInt Valoc,PetscInt Vmloc)
{
  PetscErrorCode ierr;
  DYNEXST1      exst1;

  PetscFunctionBegin;

  ierr = DYNExcModelGetModelData(dynexc,(void**)&exst1);CHKERRQ(ierr);
  /*
  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);
  row[0] = dynlocglob;row[1]=dynlocglob+2;row[2]=dynlocglob+3;
  col = Vmloc;
  if(!(exst1->VIatmax || exst1->VIatmin)) {
    value[0] =exst1->KA*exst1->TC/(exst1->TB*exst1->TA);
    value[1]=(1-(exst1->TC/exst1->TB))/exst1->TB;
    value[2]=(exst1->KF*exst1->KA*exst1->TC)/(exst1->TF*exst1->TB*exst1->TA);
    ierr = SetMatrixValues(jacP,3,row,1,col,value);CHKERRQ(ierr);
  }
  */
  PetscFunctionReturn(0);
}

/*
  DYNExcModelDAEFWDRHSJacobianP_EXST1 - Computes the RHS Jacobian of the exciter equations w.r.t. parameters

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
PetscErrorCode DYNExcModelDAEFWDRHSJacobianP_EXST1(DYNExcModel dynexc,PetscReal t,const PetscScalar *xdynexc,Vec *jacP,PetscInt dynlocglob,PetscInt Valoc,PetscInt Vmloc)
{
  PetscErrorCode ierr;
  DYNEXST1      exst1;

  PetscFunctionBegin;

  ierr = DYNExcModelGetModelData(dynexc,(void**)&exst1);CHKERRQ(ierr);

  /*
  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);
  row[0] = dynlocglob;row[1]=dynlocglob+2;row[2]=dynlocglob+3;
  col = Vmloc;
  if(!(exst1->VIatmax || exst1->VIatmin)) {
    value[0] =exst1->KA*exst1->TC/(exst1->TB*exst1->TA);
    value[1]=(1-(exst1->TC/exst1->TB))/exst1->TB;
    value[2]=(exst1->KF*exst1->KA*exst1->TC)/(exst1->TF*exst1->TB*exst1->TA);
//    ierr = SetMatrixValues(jacP,3,row,1,col,value);CHKERRQ(ierr);
    ierr = VecSetValues(jacP[Vmloc],3,row,value,INSERT_VALUES);CHKERRQ(ierr);
  }
  */
  PetscFunctionReturn(0);
}

/*
  DYNExcModelDAERHSJacobian_EXST1 - Computes the RHS Jacobian of the EXST1 DAE equations

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
PetscErrorCode DYNExcModelDAERHSJacobian_EXST1(DYNExcModel dynexc,Mat J,PetscReal t,PetscScalar VD,PetscScalar VQ,PetscScalar *x,PetscInt dynlocglob,PetscInt V_loc[])
{
  PetscErrorCode ierr;
  DYNEXST1      exst1;
  PetscInt       loc,i;
  PetscScalar    Vm,dVm_dVD,dVm_dVQ;
  PetscInt       row[2],col[10];
  PetscScalar    val[10];
  DYNStabModel   dynstab;

  PetscFunctionBegin;
  ierr = DYNExcModelGetModelData(dynexc,(void**)&exst1);CHKERRQ(ierr);

  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);

  Vm = PetscSqrtScalar(VD*VD + VQ*VQ);

  dVm_dVD = VD/Vm;
  dVm_dVQ = VQ/Vm;

  row[0] = dynlocglob;
  if(exst1->VIatmax || exst1->VIatmin) {
    val[0] = -1/exst1->TA;val[1] = exst1->KA/exst1->TA;
    col[0] = dynlocglob;col[1]=dynlocglob+2;
    ierr = SetMatrixValues(J,1,row,2,col,val);CHKERRQ(ierr);
  } else {
  col[0] = dynlocglob;col[1]=dynlocglob+1;col[2]=dynlocglob+2;col[3]=dynlocglob+3;
  val[0] = -1/exst1->TA;
  val[1] = -exst1->KA*exst1->TC/(exst1->TB*exst1->TA);
  val[2] = exst1->KA/exst1->TA;
  val[3] = -exst1->KA*exst1->TC/(exst1->TB*exst1->TA);
  ierr = SetMatrixValues(J,1,row,4,col,val);CHKERRQ(ierr);
  ierr = DYNExcModelGetDynStab(dynexc,&dynstab);CHKERRQ(ierr);
  if(dynstab) {
    PetscInt dVOTHSGdXdyn_num,dVOTHSGdXdyn_loc[6];
    PetscScalar dVOTHSGdXdyn[6],VOTHSG;
    ierr = DYNStabModelGetVOTHSG(dynstab,t,x,&VOTHSG,&dVOTHSGdXdyn_num,dVOTHSGdXdyn,dVOTHSGdXdyn_loc);CHKERRQ(ierr);
    if(dVOTHSGdXdyn_num) {
	for(i=0; i < dVOTHSGdXdyn_num; i++) dVOTHSGdXdyn[i] *= exst1->KA*exst1->TC/(exst1->TB*exst1->TA);

	ierr = SetMatrixValues(J,1,row,dVOTHSGdXdyn_num,dVOTHSGdXdyn_loc,dVOTHSGdXdyn);CHKERRQ(ierr);
      }
    }
  }

  row[0] = dynlocglob+1;
  col[0] = dynlocglob+1;col[1]= V_loc[0]; col[2] = V_loc[1];
  val[0] = -1/exst1->TR;
  val[1] = dVm_dVD/exst1->TR; val[2] = dVm_dVQ/exst1->TR;
  ierr = SetMatrixValues(J,1,row,3,col,val);CHKERRQ(ierr);

  row[0] = dynlocglob+2;
  if(exst1->VIatmax || exst1->VIatmin) {
    val[0] = -1/exst1->TB;
    col[0] = dynlocglob+2;
    ierr = SetMatrixValues(J,1,row,1,col,val);CHKERRQ(ierr);
  } else {
  col[0] = dynlocglob+1; col[1] = dynlocglob+2; col[2] = dynlocglob+3;
  val[0] = -(1-(exst1->TC/exst1->TB))/exst1->TB;
  val[1] = -1/exst1->TB;
  val[2] = -(1-(exst1->TC/exst1->TB))/exst1->TB;
  ierr = SetMatrixValues(J,1,row,3,col,val);CHKERRQ(ierr);

  ierr = DYNExcModelGetDynStab(dynexc,&dynstab);CHKERRQ(ierr);
  if(dynstab) {
    PetscInt dVOTHSGdXdyn_num,dVOTHSGdXdyn_loc[6];
    PetscScalar dVOTHSGdXdyn[6],VOTHSG;
    ierr = DYNStabModelGetVOTHSG(dynstab,t,x,&VOTHSG,&dVOTHSGdXdyn_num,dVOTHSGdXdyn,dVOTHSGdXdyn_loc);CHKERRQ(ierr);
    if(dVOTHSGdXdyn_num) {
	for(i=0; i < dVOTHSGdXdyn_num; i++) dVOTHSGdXdyn[i] *= (1-(exst1->TC/exst1->TB))/exst1->TB;
     
	ierr = SetMatrixValues(J,1,row,dVOTHSGdXdyn_num,dVOTHSGdXdyn_loc,dVOTHSGdXdyn);CHKERRQ(ierr);
      }
    }  
 }

  row[0] = dynlocglob+3;
   if(exst1->VIatmax || exst1->VIatmin) {
    val[0] = -exst1->KF/(exst1->TF*exst1->TA);
    val[1] = exst1->KA*exst1->KF/(exst1->TF*exst1->TA);
    val[2] = -1/exst1->TF;
    col[0] = dynlocglob;col[1]=dynlocglob+2;col[2]=dynlocglob+3;
    ierr = SetMatrixValues(J,1,row,3,col,val);CHKERRQ(ierr);
  } else {
     col[0] = dynlocglob; col[1] = dynlocglob+1; col[2] = dynlocglob+2; col[3] = dynlocglob+3;
     val[0] = -exst1->KF/(exst1->TF*exst1->TA);
     val[1] = -(exst1->KF*exst1->KA*exst1->TC)/(exst1->TF*exst1->TB*exst1->TA);
     val[2] = exst1->KA*exst1->KF/(exst1->TF*exst1->TA);
     val[3] = -(exst1->KF*exst1->KA*exst1->TC/(exst1->TF*exst1->TB*exst1->TA)) - 1/(exst1->TF);
     ierr = SetMatrixValues(J,1,row,4,col,val);CHKERRQ(ierr);

     ierr = DYNExcModelGetDynStab(dynexc,&dynstab);CHKERRQ(ierr);
     if(dynstab) {
       PetscInt dVOTHSGdXdyn_num,dVOTHSGdXdyn_loc[6];
       PetscScalar dVOTHSGdXdyn[6],VOTHSG;
       ierr = DYNStabModelGetVOTHSG(dynstab,t,x,&VOTHSG,&dVOTHSGdXdyn_num,dVOTHSGdXdyn,dVOTHSGdXdyn_loc);CHKERRQ(ierr);
       if(dVOTHSGdXdyn_num) {
	 for(i=0; i < dVOTHSGdXdyn_num; i++) dVOTHSGdXdyn[i] *= (exst1->KF*exst1->KA*exst1->TC)/(exst1->TF*exst1->TB*exst1->TA);
	 
	 ierr = SetMatrixValues(J,1,row,dVOTHSGdXdyn_num,dVOTHSGdXdyn_loc,dVOTHSGdXdyn);CHKERRQ(ierr);
      }
    }  
 }

  PetscFunctionReturn(0);
}

/*
  DYNExcModelDAERHSFunction_EXST1 - Computes the rhs of the DAE function for the EXST1 model 

  Input Parameters:
+ dynexc - the dynamic exciter model
. t      - the current time
. x      - array of the variables for this bus
. VD     - real-part of the bus voltage
- VQ     - imaginary part of the bus voltage

  Output Parameters:
+ f      - array of rhs of DAE equations for the exst1 model 
*/
PetscErrorCode DYNExcModelDAERHSFunction_EXST1(DYNExcModel dynexc,PetscReal t,PetscScalar VD, PetscScalar VQ, PetscScalar *x,PetscScalar *f)
{
  PetscErrorCode ierr;
  DYNEXST1      exst1;
  PetscInt       loc;
  PetscScalar    Vm;
  PetscScalar    Efd,VT,X,VF;

  DYNStabModel   dynstab;
  PetscScalar    VOTHSG=0.0;
  PetscFunctionBegin;
  /* Get the model data */
  ierr = DYNExcModelGetModelData(dynexc,(void**)&exst1);CHKERRQ(ierr);

  /* Get the location for the first exciter variable */
  ierr = DYNExcModelGetFirstVariableLocation(dynexc,&loc);CHKERRQ(ierr);

  Efd  = x[loc];
  VT  = x[loc+1];
  X   = x[loc+2];
  VF  = x[loc+3];

  ierr = DYNExcModelGetDynStab(dynexc,&dynstab);CHKERRQ(ierr);
  if(dynstab) {
    ierr = DYNStabModelGetVOTHSG(dynstab,t,x,&VOTHSG,NULL,NULL,NULL);CHKERRQ(ierr);
  }

  Vm = PetscSqrtScalar(VD*VD + VQ*VQ);

  if(exst1->VIatmax)  f[loc]   = ( (exst1->KA*exst1->TC/exst1->TB)*exst1->Vimax + exst1->KA*X - Efd )/exst1->TA;
    else if(exst1->VIatmin)  f[loc] =  ( (exst1->KA*exst1->TC/exst1->TB)*exst1->Vimin + exst1->KA*X - Efd )/exst1->TA;
    else {   f[loc]   = ( (exst1->KA*exst1->TC/exst1->TB)*(VOTHSG+exst1->Vref-VT-VF) + exst1->KA*X - Efd )/exst1->TA;
    }
 
  f[loc+1] = (-VT + Vm)/exst1->TR;

    if(exst1->VIatmax) f[loc+2] = ( exst1->Vimax*(1-(exst1->TC/exst1->TB))-X )/exst1->TB;
    else if(exst1->VIatmin) f[loc+2] = ( exst1->Vimin*(1-(exst1->TC/exst1->TB))-X )/exst1->TB;
    else {  f[loc+2] = ( (1-(exst1->TC/exst1->TB))*(VOTHSG+exst1->Vref-VT-VF) - X )/exst1->TB;
    }


    if(exst1->VIatmax) {
      f[loc+3] = ( (exst1->KF/exst1->TA)*( (exst1->KA*exst1->TC/exst1->TB)*exst1->Vimax + exst1->KA*X - Efd ) - VF )/exst1->TF;
    }
    else if(exst1->VIatmin) {
      f[loc+3] =  ( (exst1->KF/exst1->TA)*( (exst1->KA*exst1->TC/exst1->TB)*exst1->Vimin + exst1->KA*X - Efd ) - VF )/exst1->TF;
    }
    else {  
      f[loc+3] =( (exst1->KF/exst1->TA)*( (exst1->KA*exst1->TC/exst1->TB)*(VOTHSG+exst1->Vref-VT-VF) + exst1->KA*X - Efd ) - VF )/exst1->TF;
    }

  PetscFunctionReturn(0);
}

/*
  DYNExcModelCreate_EXST1 - Class constructor for EXST1 exciter model
*/
PetscErrorCode DYNExcModelCreate_EXST1(DYNExcModel dynexc)
{
  DYNEXST1 exst1;
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = PetscCalloc1(1,&exst1);CHKERRQ(ierr);
  exst1->VIatmax = exst1->VIatmin = 0;
  exst1->Efdatmax = exst1->Efdatmin = 0;

  dynexc->data = (void*)exst1;

  /* Inherit the ops */
  dynexc->ops.readdata             = DYNExcModelReadData_EXST1;
  dynexc->ops.destroy              = DYNExcModelDestroy_EXST1;
  dynexc->ops.getnvar              = DYNExcModelGetNvar_EXST1;
  dynexc->ops.getsizeof            = DYNExcModelGetSizeof_EXST1;
  dynexc->ops.setinitialconditions = DYNExcModelSetInitialConditions_EXST1;
  dynexc->ops.daerhsfunction       = DYNExcModelDAERHSFunction_EXST1;
  dynexc->ops.daerhsjacobian       = DYNExcModelDAERHSJacobian_EXST1;
  dynexc->ops.daerhsjacobianp      = DYNExcModelDAERHSJacobianP_EXST1;
  dynexc->ops.daefwdrhsjacobianp   = DYNExcModelDAEFWDRHSJacobianP_EXST1;
  dynexc->ops.getequationtypes     = DYNExcModelGetEquationTypes_EXST1;
  dynexc->ops.getbusnumid          = DYNExcModelGetBusnumID_EXST1;
  dynexc->ops.getfieldvoltage      = DYNExcModelGetFieldVoltage_EXST1;
  dynexc->ops.setevent             = DYNExcModelSetEvent_EXST1;
  PetscFunctionReturn(0);
}


