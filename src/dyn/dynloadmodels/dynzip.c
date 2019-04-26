#include "dynzip.h"
#include <private/dynloadmodelsimpl.h>

/*
  DYNEventMonitor_Zip - Event monitoring routine for this load model
*/
PetscErrorCode DYNEventMonitor_Zip(DYNLoadModel dynload, PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscScalar *fval)
{
  PetscErrorCode ierr;
  DYNZIP         zip;

  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNLoadModelGetModelData(dynload,(void**)&zip);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNEventPostDGamma_Zip - Compute derivative of gamma
*/
PetscErrorCode DYNEventPostDGamma_Zip(DYNLoadModel dynload, PetscInt nevents, PetscInt ev_list[], PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscScalar *dgdx, PetscScalar *dgdp,PetscInt Valoc, PetscInt Vmloc, PetscInt Pgloc, PetscInt Qgloc)
{
  PetscErrorCode ierr;
  DYNZIP      zip;
  PetscInt       loc;

  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNLoadModelGetModelData(dynload,(void**)&zip);CHKERRQ(ierr);

  /* Get the location for the first exciter variable */
  ierr = DYNLoadModelGetFirstVariableLocation(dynload,&loc);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNEventPostFunction_Zip - Post event routine for this load model
*/
PetscErrorCode DYNEventPostFunction_Zip(DYNLoadModel dynload, PetscInt nevents, PetscInt ev_list[], PetscReal t, PetscScalar VD, PetscScalar VQ, PetscScalar *x, PetscBool forward_solve, PetscBool *solve_algebraic)
{
  PetscErrorCode ierr;
  DYNZIP      zip;

  PetscFunctionBegin;

  /* Get the model data */
  ierr = DYNLoadModelGetModelData(dynload,(void**)&zip);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNLoadModelSetEvent_Zip - Sets the event info for the ZIP load model

  Input Parameters
. dynload - the dynamic load model

  Output Parameters
+ nmons - the number of event monitors
. direction - the event directions
. terminate - flags for termination when the event is located
. eventfcn  - the function describing the event condition
- posteventfcn - An optional function that gets called when the event is located
*/
PetscErrorCode DYNLoadModelSetEvent_Zip(DYNLoadModel dynload,PetscInt *nmons,PetscInt *direction,PetscBool *terminate,PetscErrorCode (**eventfcn)(DYNLoadModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*), PetscErrorCode (**posteventfcn)(DYNLoadModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscBool,PetscBool*), PetscErrorCode (**posteventdgamma)(DYNLoadModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscInt))
{

  PetscFunctionBegin;
  *nmons = 2;
  direction[0] = -1; direction[1] = 1;
  terminate[0] = terminate[1] = PETSC_FALSE;
  *eventfcn = DYNEventMonitor_Zip;
  *posteventfcn = DYNEventPostFunction_Zip;
  *posteventdgamma = DYNEventPostDGamma_Zip;
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelReadData_Zip - Parses the line containining the ZIP data and
  populates it in the Zip data structure

  Input Parameters:
+ dynload - The load base model
. line   - the file from the line
. mbase  - the machine base (NULL if not applicable)
- sbase  - the system base (NULL if not applicable)

  Notes: 
    The type of load model should be set prior to calling this routine via DYNLoadModelSetType()

    The conversion from machine base to system base,if applicable, should be done in this routine
*/
PetscErrorCode DYNLoadModelReadData_Zip(DYNLoadModel dynload,char* line,PetscScalar mbase,PetscScalar sbase)
{
  PetscErrorCode ierr;
  DYNZIP zip;

  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&zip);CHKERRQ(ierr);

  sscanf(line,"%d,'ZIP',%[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf",&zip->bus_i,zip->id,&zip->pl,&zip->ql,&zip->ip,&zip->iq,&zip->yp,&zip->yq,&zip->Vm_thresh);

  if(PetscAbsScalar(zip->pl + zip->ip + zip->yp - 100.0) > PETSC_SMALL) SETERRQ4(PETSC_COMM_SELF,0,"Bus %d -- Inconsistent zip load model composition: cp = %f + ip = %f + yp = %f should equal 100.0",zip->bus_i,zip->pl,zip->ip,zip->yp);
  if(PetscAbsScalar(zip->ql + zip->iq + zip->yq - 100.0) > PETSC_SMALL) SETERRQ4(PETSC_COMM_SELF,0,"Bus %d -- Inconsistent zip load model composition: cq = %f + iq = %f + yq = %f should equal 100.0",zip->bus_i,zip->ql,zip->iq,zip->yq);
  

  PetscFunctionReturn(0);
}

/*
  DYNLoadModelGetBusnumID_Zip - Returns the bus number and ID associated with this load model

  Input Parameters:
. dynload - the dynamic load model object

  Output Parameters:
+ busnum - the bus number 
- loadid  - the load ID
*/
PetscErrorCode DYNLoadModelGetBusnumID_Zip(DYNLoadModel dynload,PetscInt *busnum,char **loadid)
{
  PetscErrorCode ierr;
  DYNZIP         zip;

  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&zip);CHKERRQ(ierr);
  *busnum = zip->bus_i;
  *loadid  = zip->id;
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelDestroy_Zip - Destroys the Zip object data

  Input Parameters
. DYNLoadModel - the DYNLoadModel object

  Notes:
  Called when DYNLoadModelDestroy() is called
*/
PetscErrorCode DYNLoadModelDestroy_Zip(DYNLoadModel dynload)
{
  PetscErrorCode ierr;
  DYNZIP         zip;

  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&zip);CHKERRQ(ierr);
  ierr = PetscFree(zip);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelGetNvar_Zip - Returns the number of variables for the Zip model

  Input parameters
. dynload - the dynamic load object

  Output Parameters
. nvar - number of variables for this model
*/
PetscErrorCode DYNLoadModelGetNvar_Zip(DYNLoadModel dynload,PetscInt *nvar)
{
  PetscFunctionBegin;
  *nvar = DYNZIP_nvar;
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelGetSizeof_Zip - Returns the size of Zip struct

  Input Parameters:
. dynload - the dynamic load model
 
  Output Parameters:
. size - size of the Zip object obtained from sizeof()
*/
PetscErrorCode DYNLoadModelGetSizeof_Zip(DYNLoadModel dynload,PetscInt *size)
{
  PetscFunctionBegin;
  *size = sizeof(struct _p_DYNZIP);
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelSetInitialConditions_Zip - Sets the initial conditions (x(t0)) for the ZIP model

  Input Parameters:
+ dynload - the DYNLoadModel object
. PD     - load real power output
. QD     - load reactive power output
. VD     - real component of the complex bus voltage
. VQ     - imaginary component of the complex bus voltage

  Output Parameters:
. x - the initial conditions for the ZIP model

  Notes the locations to insert the values should be obtained by DYNLoadModelGetFirstVariableLocation()
*/
PetscErrorCode DYNLoadModelSetInitialConditions_Zip(DYNLoadModel dynload,PetscScalar PD, PetscScalar QD, PetscScalar VD, PetscScalar VQ, PetscScalar *x)
{
  PetscErrorCode ierr;
  DYNZIP        zip;
  PetscScalar Vm;
  PetscInt    loc;

  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&zip);CHKERRQ(ierr);

  ierr = DYNLoadModelGetFirstVariableLocation(dynload,&loc);CHKERRQ(ierr);

  Vm = PetscSqrtScalar(VD*VD + VQ*VQ);

  zip->pl = zip->pl*PD/100.0;
  zip->ql = zip->ql*QD/100.0;
  zip->ip = zip->ip*PD/Vm/100.0;
  zip->iq = zip->iq*QD/Vm/100.0;
  zip->yp = zip->yp*PD/(Vm*Vm)/100.0;
  zip->yq = zip->yq*QD/(Vm*Vm)/100.0;

  zip->Vm0 = Vm;
  
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelGetEquationTypes_Zip - Gets the number and indices of differential and algebraic equations
                                                   
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
PetscErrorCode DYNLoadModelGetEquationTypes_Zip(DYNLoadModel dynload,PetscInt *ndiff,PetscInt *nalg, PetscInt *eqtype)
{
  DYNZIP zip;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&zip);CHKERRQ(ierr);

  *ndiff = 0;
  *nalg  = 0;

  PetscFunctionReturn(0);
}

/*
  DYNLoadModelDAERHSJacobian_Zip - Computes the RHS Jacobian of the Zip DAE equations

  Input Parameters:
+ dynload     - the dynamic load model
. t          - the current time
. VD         - real-part of complex bus voltage
. VQ         - imaginary-part of complex bus voltage
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

   The location for the first variable for the zip model should be obtained by
   DYNLoadModelGetFirstVariableLocation(dynload,&loc). x[loc] = 1st zip variable,
   x[loc+n-1] = nth zip variable
*/ 
PetscErrorCode DYNLoadModelDAERHSJacobian_Zip(DYNLoadModel dynload,Mat J,PetscReal t,PetscScalar VD,PetscScalar VQ,PetscScalar *x,PetscInt dynlocglob,PetscInt V_loc[],PetscInt I_loc[])
{
  PetscErrorCode ierr;
  DYNZIP      zip;
  PetscScalar    Iloaddp,Iloadqp,Iloaddi,Iloadqi;
  PetscScalar    dIloaddp_dVD,dIloaddp_dVQ,dIloadqp_dVD,dIloadqp_dVQ;
  PetscScalar    dIloaddi_dVD,dIloaddi_dVQ,dIloadqi_dVD,dIloadqi_dVQ;
  PetscScalar    dIloaddz_dVD,dIloaddz_dVQ,dIloadqz_dVD,dIloadqz_dVQ;
  PetscScalar    Vm,Vm2;
  PetscInt       row,col[2];
  PetscScalar    val[2];

  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&zip);CHKERRQ(ierr);

  Vm = PetscSqrtScalar(VD*VD + VQ*VQ);
  Vm2 = Vm*Vm;
  
  if(Vm > zip->Vm_thresh) {
    Iloaddp = (zip->pl*VD + zip->ql*VQ)/Vm2;
    Iloadqp = (-zip->ql*VD + zip->pl*VQ)/Vm2;

    dIloaddp_dVD = zip->pl/Vm2  - 2*Iloaddp*VD/Vm2;
    dIloaddp_dVQ = zip->ql/Vm2  - 2*Iloaddp*VQ/Vm2;
    dIloadqp_dVD = -zip->ql/Vm2 - 2*Iloadqp*VD/Vm2;
    dIloadqp_dVQ = zip->pl/Vm2  - 2*Iloadqp*VQ/Vm2;

  } else {
    Iloaddp = (zip->pl*VD + zip->ql*VQ)/(zip->Vm_thresh*zip->Vm_thresh);
    Iloadqp = (-zip->ql*VD + zip->pl*VQ)/(zip->Vm_thresh*zip->Vm_thresh);

    dIloaddp_dVD = zip->pl/(zip->Vm_thresh*zip->Vm_thresh);
    dIloaddp_dVQ = zip->ql/(zip->Vm_thresh*zip->Vm_thresh);
    dIloadqp_dVD = -zip->ql/(zip->Vm_thresh*zip->Vm_thresh);
    dIloadqp_dVQ = zip->pl/(zip->Vm_thresh*zip->Vm_thresh);
  }
  
  if(Vm > 0.5) {
    Iloaddi = (zip->ip*VD + zip->iq*VQ)/Vm;
    Iloadqi = (-zip->iq*VD + zip->ip*VQ)/Vm;

    dIloaddi_dVD = zip->ip/Vm  - Iloaddi*VD/Vm2;
    dIloaddi_dVQ = zip->ip/Vm  - Iloaddi*VQ/Vm2;
    dIloadqi_dVD = -zip->ip/Vm - Iloadqi*VD/Vm2;
    dIloadqi_dVQ = zip->ip/Vm  - Iloadqi*VQ/Vm2;
  } else {
    Iloaddi = (zip->ip*VD + zip->iq*VQ)/0.5;
    Iloadqi = (-zip->iq*VD + zip->ip*VQ)/0.5;

    dIloaddi_dVD = zip->ip/0.5;
    dIloaddi_dVQ = zip->ip/0.5;
    dIloadqi_dVD = -zip->ip/0.5;
    dIloadqi_dVQ = zip->ip/0.5;
  }

  dIloaddz_dVD = zip->yp;
  dIloaddz_dVQ = zip->yq;
  dIloadqz_dVD = -zip->yq;
  dIloadqz_dVQ = zip->yp;

  row  = I_loc[0];
  col[0] = V_loc[0]; col[1] = V_loc[1];
  val[0] = -(dIloaddp_dVD + dIloaddi_dVD + dIloaddz_dVD);
  val[1] = -(dIloaddp_dVQ + dIloaddi_dVQ + dIloaddz_dVQ);

  ierr = SetMatrixValues(J,1,&row,2,col,val);CHKERRQ(ierr);

  row = I_loc[1];
  col[0] = V_loc[0]; col[1] = V_loc[1];
  val[0] = -(dIloadqp_dVD + dIloadqi_dVD + dIloadqz_dVD);
  val[1] = -(dIloadqp_dVQ + dIloadqi_dVQ + dIloadqz_dVQ);
  
  ierr = SetMatrixValues(J,1,&row,2,col,val);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNLoadModelDAERHSJacobianP_Zip - Sets the Jacobian of the Zip equations w.r.t. parameters

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
PetscErrorCode DYNLoadModelDAERHSJacobianP_Zip(DYNLoadModel dynload, PetscReal t,const PetscScalar *x,Mat jacP,PetscInt dynlocglob,PetscInt valoc,PetscInt vmloc)
{
  DYNZIP         zip;
  PetscInt       row[2],col[1];
  PetscScalar    VD,VQ,VM,ip_vm,iq_vm,yp_vm,yq_vm,Iloaddp_vm,Iloadqp_vm,Iloaddi_vm,Iloadqi_vm,Iloaddz_vm,Iloadqz_vm,val[2];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&zip);CHKERRQ(ierr);
  VD = x[0]; VQ = x[1];
  VM = PetscSqrtScalar(VD*VD + VQ*VQ);

  /* assume bus->vm == zip->Vm0 */
  ip_vm = -zip->ip/zip->Vm0;
  iq_vm = -zip->iq/zip->Vm0;
  yp_vm = -2.*zip->yp/zip->Vm0;
  yq_vm = -2.*zip->yq/zip->Vm0;

  if (VM > zip->Vm_thresh) {
    Iloaddp_vm = 0.;
    Iloadqp_vm = 0.;
  } else {
    Iloaddp_vm = -2.*(zip->pl*VD + zip->ql*VQ)/(zip->Vm0*zip->Vm0*zip->Vm0);
    Iloadqp_vm = -2.*(-zip->ql*VD + zip->pl*VQ)/(zip->Vm0*zip->Vm0*zip->Vm0);
  }

  if (VM > 0.5) {
    Iloaddi_vm = (ip_vm*VD + iq_vm*VQ)/VM;
    Iloadqi_vm = (-iq_vm*VD + ip_vm*VQ)/VM;
  } else {
    Iloaddi_vm = (ip_vm*VD + iq_vm*VQ)/zip->Vm0  - (zip->ip*VD + zip->iq*VQ)/(zip->Vm0*zip->Vm0);
    Iloadqi_vm = (-iq_vm*VD + ip_vm*VQ)/zip->Vm0 -(-zip->iq*VD + zip->ip*VQ)/(zip->Vm0*zip->Vm0);
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
  DYNLoadModelDAEFWDRHSJacobianP_Zip - Sets the Jacobian of the Zip equations w.r.t. parameters

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
PetscErrorCode DYNLoadModelDAEFWDRHSJacobianP_Zip(DYNLoadModel dynload, PetscReal t,const PetscScalar *x,Vec *jacP,PetscInt dynlocglob,PetscInt Valoc,PetscInt Vmloc)
{
  PetscErrorCode ierr;
  DYNZIP      zip;
  PetscInt    loc;

  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&zip);CHKERRQ(ierr);
  ierr = DYNLoadModelGetFirstVariableLocation(dynload,&loc);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNLoadModelDAERHSFunction_Zip - Computes the rhs of the DAE function for the ZIP model 
                                and also returns the load currents in network refernce frame

  Input Parameters:
+ dynload - the dynamic load model
. t      - the current time
. x      - array of all the variables for the bus on which this load is incident
. VD     - real-part of the bus voltage
- VQ     - imaginary part of the bus voltage

  Output Parameters:
+ f      - array of rhs of DAE equations for the zip model 
. ILD    - real-part of load current in network reference frame
. ILQ    - imaginary-part of load current in network reference frame

  Notes:
   The locations for the variables (and the corresponding locations to insert entries in f)
   for this load model should be obtained by DYNLoadModelGetFirstVariableLocation()
*/
PetscErrorCode DYNLoadModelDAERHSFunction_Zip(DYNLoadModel dynload,PetscReal t,PetscScalar VD, PetscScalar VQ, PetscScalar *x,PetscScalar *f,PetscScalar *ILD,PetscScalar *ILQ)
{
  PetscErrorCode ierr;
  DYNZIP         zip;
  PetscScalar    Iloaddp,Iloadqp,Iloaddi,Iloadqi,Iloaddz,Iloadqz;
  PetscScalar    Vm,Vm2;
  
  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&zip);CHKERRQ(ierr);

  Vm = PetscSqrtScalar(VD*VD + VQ*VQ);
  Vm2 = Vm*Vm;

  if(Vm > zip->Vm_thresh) {
    Iloaddp = (zip->pl*VD + zip->ql*VQ)/Vm2;
    Iloadqp = (-zip->ql*VD + zip->pl*VQ)/Vm2;
  } else {
    Iloaddp = (zip->pl*VD + zip->ql*VQ)/(zip->Vm_thresh*zip->Vm_thresh);
    Iloadqp = (-zip->ql*VD + zip->pl*VQ)/(zip->Vm_thresh*zip->Vm_thresh);
  }

  if(Vm > 0.5) {
    Iloaddi = (zip->ip*VD + zip->iq*VQ)/Vm;
    Iloadqi = (-zip->iq*VD + zip->ip*VQ)/Vm;
  } else {
    Iloaddi = (zip->ip*VD + zip->iq*VQ)/0.5;
    Iloadqi = (-zip->iq*VD + zip->ip*VQ)/0.5;
  }
    
  Iloaddz = zip->yp*VD + zip->yq*VQ;
  Iloadqz = -zip->yq*VD + zip->yp*VQ;
  
  *ILD = Iloaddp + Iloaddi + Iloaddz;
  *ILQ = Iloadqp + Iloadqi + Iloadqz;

  PetscFunctionReturn(0);
}

PetscErrorCode DYNLoadModelSetConstantPowerLoad_Zip(DYNLoadModel dynload,PetscScalar Pl, PetscScalar Ql)
{
  PetscErrorCode ierr;
  DYNZIP         zip;
  
  PetscFunctionBegin;
  ierr = DYNLoadModelGetModelData(dynload,(void**)&zip);CHKERRQ(ierr);

  //  zip->pl = Pl;
  //  zip->ql = Ql;
  zip->yp = Pl;
  zip->yq = Ql;

  PetscFunctionReturn(0);
}


/*
  DYNLoadModelCreate_Zip - Class constructor for ZIP model
*/
PetscErrorCode DYNLoadModelCreate_ZIP(DYNLoadModel dynload)
{
  DYNZIP zip;
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = PetscCalloc1(1,&zip);CHKERRQ(ierr);
  dynload->data = (void*)zip;

  /* Inherit the ops */
  dynload->ops.readdata             = DYNLoadModelReadData_Zip;
  dynload->ops.destroy              = DYNLoadModelDestroy_Zip;
  dynload->ops.getnvar              = DYNLoadModelGetNvar_Zip;
  dynload->ops.getsizeof            = DYNLoadModelGetSizeof_Zip;
  dynload->ops.setinitialconditions = DYNLoadModelSetInitialConditions_Zip;
  dynload->ops.daerhsfunction       = DYNLoadModelDAERHSFunction_Zip;
  dynload->ops.daerhsjacobian       = DYNLoadModelDAERHSJacobian_Zip;
  dynload->ops.daerhsjacobianp      = DYNLoadModelDAERHSJacobianP_Zip;
  dynload->ops.daefwdrhsjacobianp   = DYNLoadModelDAEFWDRHSJacobianP_Zip;
  dynload->ops.getequationtypes     = DYNLoadModelGetEquationTypes_Zip;
  dynload->ops.getbusnumid          = DYNLoadModelGetBusnumID_Zip;
  //  dynload->ops.setevent             = DYNLoadModelSetEvent_Zip;
  dynload->ops.setconstantpowerload = DYNLoadModelSetConstantPowerLoad_Zip;

  PetscFunctionReturn(0);
}

