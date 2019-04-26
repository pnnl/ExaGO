#include <private/dynloadmodelsimpl.h>
#include <private/psimpl.h>

struct _p_DYNLoadModelList DYNLoadModelList[DYNLOADMODELSMAX];

PetscInt nloadmodelsregistered=0;

PetscBool  DYNLoadModelsRegisterAllCalled = PETSC_FALSE;

/*
  DYNLoadModelEventMonitor - Sets the event info for the load model

  Input Parameters
. dynload - the dynamic load model

*/
PetscErrorCode DYNLoadModelSetEventMonitor(DYNLoadModel dynload)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if(dynload->ops.setevent) {
    ierr = (*dynload->ops.setevent)(dynload,&dynload->numMonitors,dynload->direction,dynload->terminate,&dynload->eventfcn,&dynload->posteventfcn,&dynload->posteventdgamma);CHKERRQ(ierr);
  } else { 
    dynload->numMonitors = 0;
  }
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelSetStatus - Sets the load status (1 = ON, 0 = OFF)

  Inputs:
+ dynload - the DYNLoadModel object
- status - the load status
*/
PetscErrorCode DYNLoadModelSetStatus(DYNLoadModel dynload,PetscInt status)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PSLOADSetStatus(dynload->psload,status);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelDestroy - Destroys the dynload model object.

  Input Parameters
. dynload - the dynamic load model

  Notes: 
   This routine is only called for dynamic simulation. It calls the destroy method on
   the child class to free its data
*/
PetscErrorCode DYNLoadModelDestroy(DYNLoadModel dynload)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynload->ops.destroy)(dynload);CHKERRQ(ierr);
  ierr = PetscFree(dynload->eqtypes);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelRegister - Registers a dynamic load model

  Input Parameters:
+ sname     - model name (string)
. createfunction  - the class constructor
- sizeofstruct  - size of object (obtained with sizeof())
*/
PetscErrorCode DYNLoadModelRegister(const char sname[],PetscErrorCode (*createfunction)(DYNLoadModel),PetscInt sizeofstruct)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscFunctionBegin;
  for(i=0; i < nloadmodelsregistered;i++) {
    PetscBool match;
    ierr = PetscStrcmp(DYNLoadModelList[i].name,sname,&match);
    if(match) PetscFunctionReturn(0);
  }
  i = nloadmodelsregistered;
  ierr = PetscStrcpy(DYNLoadModelList[i].name,sname);CHKERRQ(ierr);
  DYNLoadModelList[i].create = createfunction;
  DYNLoadModelList[i].sizeofstruct = sizeofstruct;
  nloadmodelsregistered++;
  PetscFunctionReturn(0);
}

#include "dynzip.h"
#include "dyncompload.h"
extern PetscErrorCode DYNLoadModelCreate_ZIP(DYNLoadModel);
extern PetscErrorCode DYNLoadModelCreate_COMPLOAD(DYNLoadModel);

/*
  DYNLoadModelRegisterAll - Registers all built-in dynamic load models

*/
PetscErrorCode DYNLoadModelRegisterAll(void)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(!DYNLoadModelsRegisterAllCalled) DYNLoadModelsRegisterAllCalled = PETSC_TRUE;

  ierr = DYNLoadModelRegister("ZIP",DYNLoadModelCreate_ZIP,sizeof(struct _p_DYNZIP));CHKERRQ(ierr);
  ierr = DYNLoadModelRegister("COMPLOAD",DYNLoadModelCreate_COMPLOAD,sizeof(struct _p_DYNCOMPLOAD));CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNLoadModelGetModelType - Parses a line in dyr file and reads the model type (string identifier). If the line does not have
                            the registered load model then the type is set to NONE.

  Input Parameters:
. line - the line to parse

  Output Parameters:
. modeltype - type of the model (string name of the model)
*/
PetscErrorCode DYNLoadModelReadModelType(char *line, char loadtype[])
{
  PetscInt i;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for(i=0; i < nloadmodelsregistered; i++) {
    if(strstr(line,DYNLoadModelList[i].name) != NULL) {
      ierr = PetscStrcpy(loadtype,DYNLoadModelList[i].name);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
  }

  ierr = PetscStrcpy(loadtype,DYNLOADNONE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelGetSpeedDeviationLocation - Gets the location of the speed deviation variable for the load model relative to the bus variable starting location

  Input Parameters:
. dynload - the load model

  Output Parameters:
. loc - the location of the speed deviation

  Notes:
  The location of the speed deviation is relative to the bus variables. So, one should first obtain
  the first variable location for the load model implementation using DYNLoadModelGetFirstVariableLocation() 
  and then add the location of the dw variable to return loc. For example, if dw is the 3rd variable for this
  load model implementation then the following code should be used in the implementation
  ierr = DYNLoadModelGetFirstVariable(dynload,&firstvarloc);
  *loc = firstvarloc + 2;
*/
PetscErrorCode DYNLoadModelGetSpeedDeviationLocation(DYNLoadModel dynload, PetscInt *loc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynload->ops.getspeeddeviationlocation)(dynload,loc);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNLoadModelGetSpeedDeviationGlobalLocation - Gets the global location (location in the global vector) of the field voltage for the exciter model

  Input Parameters:
. dynload - the load model

  Output Parameters:
. locglob - the global location of the speed deviation variable
*/
PetscErrorCode DYNLoadModelGetSpeedDeviationGlobalLocation(DYNLoadModel dynload, PetscInt *locglob)
{
  PetscErrorCode ierr;
  PetscInt       dwloc;

  PetscFunctionBegin;

  /* Get the location of dw location relative to the bus */
  ierr = DYNLoadModelGetSpeedDeviationLocation(dynload,&dwloc);CHKERRQ(ierr);
  *locglob = dynload->bus->startlocglob + dwloc;

  PetscFunctionReturn(0);
}

/*
  DYNLoadModelGetSpeedDeviation - Returns the speed deviation and its time derivative
                                 (optional) of the load at the current time.

  Input Parameters:
+ dynload - the DYNLoadModel object
. t      - the current time
. xdyn   - the state variables for this load at time t

  Output Parameters:
+ spd     - the load speed deviation
- dspd_dt - time-derivative of load speed deviation (set to NULL if not required)
*/
PetscErrorCode DYNLoadModelGetSpeedDeviation(DYNLoadModel dynload,PetscReal t, const PetscScalar* xdyn,PetscScalar* spd, PetscScalar* dspd_dt)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynload->ops.getspeeddeviation)(dynload,t,xdyn,spd,dspd_dt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelGetFrequency - Returns the frequency of the load at the current time.

  Input Parameters:
+ dynload - the DYNLoadModel object
. t      - the current time
. xdyn   - the state variables for this load at time t

  Output Parameters:
. frequency   - the load frequency in Hz
*/
PetscErrorCode DYNLoadModelGetFrequency(DYNLoadModel dynload,PetscReal t, const PetscScalar* xdyn,PetscScalar* frequency)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynload->ops.getfrequency)(dynload,t,xdyn,frequency);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelSetFrequency - Set the frequency of the load at the current time.

  Input Parameters:
+ dynload - the DYNLoadModel object
. t      - the current time
. xdyn   - the state variables for this load at time t
. frequency = the value to be set

  Output Parameters:
*/
PetscErrorCode DYNLoadModelSetFrequency(DYNLoadModel dynload,PetscScalar* xdyn,PetscScalar frequency)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynload->ops.setfrequency)(dynload,xdyn,frequency);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelGetFrequencyLimits - Returns the min. and max. limits for the frequency of the load

  Input Parameters:
. dynload - the DYNLoadModel object


  Output Parameters:
+ fmax   - the max. load frequency allowed in Hz.
. fmin   - the min. load frequency allowed in Hz.
*/
PetscErrorCode DYNLoadModelGetFrequencyLimits(DYNLoadModel dynload,PetscScalar* fmax,PetscScalar* fmin)
{

  PetscFunctionBegin;
  *fmax = dynload->freq_max;
  *fmin = dynload->freq_min;
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelGetdFreqdState - Returns the partial derivative of the load frequency w.r.t the machine state that governs it

  Input Parameters:
+ dynload - the DYNLoadModel object
. t      - the current time
. xdyn   - the state variables for this load at time t

  Output Parameters:
+ dfreqdstate   - dfreq_dstate
- stateloc      - location w.r.t. to the bus variables
*/
PetscErrorCode DYNLoadModelGetdFreqdState(DYNLoadModel dynload,PetscReal t, const PetscScalar* xdyn,PetscScalar* dfreqdstate,PetscInt *stateloc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynload->ops.getdfreqdstate)(dynload,t,xdyn,dfreqdstate,stateloc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*
  DYNLoadModelGetBusnumID - Gets the bus number and load ID for this dynamic load model
                            from the data file. Calls the underlying implementation
  Input Parameters:
. dynload - the dynamic load model object

  Output Parameters:
+ busnum - bus number
. loadid  - the load ID
*/
PetscErrorCode DYNLoadModelGetBusnumID(DYNLoadModel dynload, PetscInt *busnum,char **loadid)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = (*dynload->ops.getbusnumid)(dynload,busnum,loadid);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelGetSizeof - Gets the size of the implementation type (obtained using sizeof())

  Input Parameters:
. dynload - the DYNLoadModel object

  Output Parameters:
. size - the size of the derived class object
*/
PetscErrorCode DYNLoadModelGetSizeof(DYNLoadModel dynload,PetscInt *size)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynload->ops.getsizeof)(dynload,size);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelGetNvar - Returns the number of variables for this dynamic load model

  Input Parameters:
. dynload - the DYNLoadModel object

  Output Parameters:
. Nvar - number of variables (dofs) for this load model
*/
PetscErrorCode DYNLoadModelGetNvar(DYNLoadModel dynload,PetscInt *Nvar)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynload->ops.getnvar)(dynload,Nvar);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelGetFirstVariableLocation - Gets the location for the first variable for the
                        dynload model from the bus variables array

  Input Parameters:
. dynload - the DYNLoadModel object

  Output Parameters
. loc - the location for the first variable 
*/
PetscErrorCode DYNLoadModelGetFirstVariableLocation(DYNLoadModel dynload,PetscInt *loc)
{
  PetscFunctionBegin;
  *loc = dynload->startloc;
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelGetInitialMechanicalPower - Returns the mechanical power at t=t0

  Input Parameters:
. dynload - the DYNLoadModel object

  Output Parameters:
. Pmech0 - Initial mechanical power at t=t0
*/
PetscErrorCode DYNLoadModelGetInitialMechanicalPower(DYNLoadModel dynload,PetscScalar *Pmech0)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynload->ops.getinitialmechanicalpower)(dynload,Pmech0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelSetType - Sets the type of the dynamic load model

  Input Parameters:
+ dynload - the dynamic load model object
- modeltype - type of the model
*/
PetscErrorCode DYNLoadModelSetType(DYNLoadModel dynload,DYNLoadModelType modeltype)
{
  PetscErrorCode ierr,(*r)(DYNLoadModel)=NULL;
  PetscInt       i;
  PetscFunctionBegin;
  PetscBool match;
  for(i=0;i < nloadmodelsregistered;i++) {
    ierr = PetscStrcmp(DYNLoadModelList[i].name,modeltype,&match);CHKERRQ(ierr);
    if(match) {
      r = DYNLoadModelList[i].create;
      break;
    }
  }

  if(!r) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown type for DYNLoadModel %s",modeltype);
  /* Null the function pointers */
  dynload->ops.readdata = 0;
  dynload->ops.getnvar = 0;
  dynload->ops.destroy = 0;
  dynload->ops.getsizeof = 0;
  dynload->ops.setinitialconditions = 0;
  dynload->ops.daerhsfunction=0;
  dynload->ops.daerhsjacobian=0;
  dynload->ops.getequationtypes=0;
  dynload->ops.getbusnumid=0;
  dynload->ops.getfrequency=0;
  dynload->ops.getspeeddeviation=0;
  dynload->ops.getspeeddeviationlocation=0;
  dynload->ops.getinitialmechanicalpower=0;
  dynload->ops.setevent=0;
  dynload->ops.setconstantpowerload=0;

  /* Copy the type name */
  ierr = PetscStrcpy(dynload->type,modeltype);CHKERRQ(ierr);

  /* Call the underlying implementation constructor */
  ierr = (*r)(dynload);CHKERRQ(ierr);
  /* Set the numnber of variables for this dynamic load model */
  //  ierr = DYNLoadModelGetNvar(dynload,&dynload->nvar);CHKERRQ(ierr);

  dynload->startloc = -1;
  dynload->nvar     =  0;
  dynload->ndiff    =  0;
  dynload->nalg     =  0;
  dynload->eqtypes  =  0;
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelCopy - Copies the dynload model

  Input Parameters:
. dynloadin - the dynload model to copy from

  Output Parameters:
. dynloadout - the dynload model to copy to

  Notes:
   If a new method is created for the dynloadmodel then that should be copied too.
   The underlying implementation (dynload->data) data pointer is copied.
*/
PetscErrorCode DYNLoadModelCopy(DYNLoadModel dynloadout,DYNLoadModel dynloadin)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscMemcpy(dynloadout,dynloadin,sizeof(struct _p_DYNLoadModel));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelGetModelData - Returns the model data associated with the particular implementation type

  Input Parameters:
. dynload - the DYNLoadModel object

  Output Parameters:
. modeldata - the model data struct
*/
PetscErrorCode DYNLoadModelGetModelData(DYNLoadModel dynload,void **modeldata)
{
  PetscFunctionBegin;
  *modeldata = (void*)dynload->data;
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelReadData - Parses the given line and reads the data for the dynamic load
  model

  Input Parameters:
+ dynload - The load base model
. line   - the file from the line
. mbase  - the machine base
- sbase  - the system base

  Notes: 
    The type of load model should be set prior to calling this routine via DYNLoadModelSetType()

    The conversion from machine base to system base,if applicable, should be done in this rouine
*/
PetscErrorCode DYNLoadModelReadData(DYNLoadModel dynload,char *line,PetscScalar mbase,PetscScalar sbase)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynload->ops.readdata)(dynload,line,mbase,sbase);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelSetInitialConditions - Sets the initial conditions for this load model

  Input Parameters:
+ dynload - the dynamic load model object
. Pg     - the real load power output
. Qg     - the reactive load power output
. VD     - real component of the bus voltage
- VQ     - imaginary component of the bus voltage

  Output Parameters:
. x   - the array of the variables for the bus on which this load is incident.

  Notes: 
   The initial conditions are populated in the array x

   Other constants, parameters can be also set during this function call in
   its implementation struct (for e.g. setting mechanical power Pm etc.)

   The locations to insert the variables for this load should be obtained from
   DYNLoadModelGetFirstVariableLocation()
   
   x contains all the variables for the load bus
*/
PetscErrorCode DYNLoadModelSetInitialConditions(DYNLoadModel dynload,PetscScalar Pg,PetscScalar Qg,PetscScalar VD,PetscScalar VQ,PetscScalar *x)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynload->ops.setinitialconditions)(dynload,Pg,Qg,VD,VQ,x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelDAERHSJacobian - Computes the RHS Jacobian of the DAE equations

  Input Parameters:
+ dynload     - the dynamic load model
. t          - the current time
. VD         - real-part of complex bus voltage
. VQ         - imaginary-part of complex bus voltage
. xdynload    - load variables
. dynlocglob - starting location of load variables in the global vector 
. V_loc      - global location of VD and VQ V_loc[0] = VD_loc V_loc[1] = VQ_loc
- I_loc      - global location of ID and IQ I_loc[0] = ID_loc I_loc[1] = IQ_loc

  Output Parameters:
. J          - the Jacobian matrix

  Notes:
   For each bus, the network equations are ordered as [IQ_loc;ID_loc]
*/ 
PetscErrorCode DYNLoadModelDAERHSJacobian(DYNLoadModel dynload,Mat J,PetscReal t,PetscScalar VD,PetscScalar VQ,PetscScalar *xdynload,PetscInt dynlocglob,PetscInt V_loc[],PetscInt I_loc[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynload->ops.daerhsjacobian)(dynload,J,t,VD,VQ,xdynload,dynlocglob,V_loc,I_loc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelDAERHSJacobianP - Sets the Jacobian of the DAE w.r.t. parameters

  Input Parameters:
+ dynload - dynloadmodel object
. t      - current time
. x      - state variables for this bus
. dynlocglob - starting location of the variables for this dynload model
. Valoc  - location of parameter voltage angle (Va)
- Vmloc  - location of parameter voltage magnitude (Vm)

  Output Parameters:
. jacP - Matrix of partial derivatives of DYN DAE equations w.r.t. parameters
*/
PetscErrorCode DYNLoadModelDAERHSJacobianP(DYNLoadModel dynload, PetscReal t,const PetscScalar *x,Mat jacP,PetscInt dynlocglob,PetscInt Valoc,PetscInt Vmloc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynload->ops.daerhsjacobianp)(dynload,t,x,jacP,dynlocglob,Valoc,Vmloc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelDAEFWDRHSJacobianP - Sets the Jacobian of the DAE w.r.t. parameters

  Input Parameters:
+ dynload - dynloadmodel object
. t      - current time
. x      - state variables for this bus
. dynlocglob - starting location of the variables for this dynload model
. Valoc  - location of parameter voltage angle (Va)
- Vmloc  - location of parameter voltage magnitude (Vm)

  Output Parameters:
. jacP - Matrix of partial derivatives of DYN DAE equations w.r.t. parameters
*/
PetscErrorCode DYNLoadModelDAEFWDRHSJacobianP(DYNLoadModel dynload, PetscReal t,const PetscScalar *x,Vec *jacP,PetscInt dynlocglob,PetscInt Valoc,PetscInt Vmloc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynload->ops.daefwdrhsjacobianp)(dynload,t,x,jacP,dynlocglob,Valoc,Vmloc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNLoadModelDAERHSFunction - Computes the RHS of the DAE function for the load model [f(x,y);g(x,y)] 
                              and also returns the load currents in network refernce frame

  Input Parameters:
+ dynload - the dynamic load model
. t      - the current time
. x      - array of the variables for the bus on which this load is incident
. VD     - real-part of bus voltage
. VQ     - imaginary part of bus voltage

  Output Parameters:
+ f      - array of rhs of DAE equations for this load
. IGD    - real-part of load current in network reference frame
. IGQ    - imaginary-part of load current in network reference frame

  Notes:
   The locations for the variables (and the corresponding locations to insert entries in f)
   for this load model should be obtained by DYNLoadModelGetFirstVariableLocation()

*/
PetscErrorCode DYNLoadModelDAERHSFunction(DYNLoadModel dynload,PetscReal t,PetscScalar VD, PetscScalar VQ, PetscScalar *x,PetscScalar *f,PetscScalar *IGD,PetscScalar *IGQ)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynload->ops.daerhsfunction)(dynload,t,VD,VQ,x,f,IGD,IGQ);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
			 
/*
  DYNLoadModelGetEquationTypes - Gets the number and indices of differential and algebraic equations for the load model
                                                   
  Input Parameters:
. dynload - the DYNLoadModel object

  Output Parameters:
+ ndiff - number of differential equations
. nalg  - number of algebraic equations
- eqtype - an array of size nvar that has
            eqtype[i] = DIFF_EQ if equation i is a differential
            eqtype[i] = ALG_EQ if equation  i is an algebraic

  NOTES: The user does not need to create the array vartype
*/
PetscErrorCode DYNLoadModelGetEquationTypes(DYNLoadModel dynload,PetscInt *ndiff,PetscInt *nalg, PetscInt *eqtype)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynload->ops.getequationtypes)(dynload,ndiff,nalg,eqtype);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DYNLoadModelSetConstantPowerLoad(DYNLoadModel dynload,PetscScalar Pl, PetscScalar Ql)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynload->ops.setconstantpowerload)(dynload,Pl,Ql);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
