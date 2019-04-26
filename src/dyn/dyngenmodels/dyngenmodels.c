#include <private/dyngenmodelsimpl.h>
#include <private/psimpl.h>

struct _p_DYNGenModelList DYNGenModelList[DYNGENMODELSMAX];

PetscInt ngenmodelsregistered=0;

PetscBool  DYNGenModelsRegisterAllCalled = PETSC_FALSE;

/*
  DYNGenModelEventMonitor - Sets the event info for the generator model

  Input Parameters
. dyngen - the dynamic generator model

*/
PetscErrorCode DYNGenModelSetEventMonitor(DYNGenModel dyngen)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if(dyngen->ops.setevent) {
    ierr = (*dyngen->ops.setevent)(dyngen,&dyngen->numMonitors,dyngen->direction,dyngen->terminate,&dyngen->eventfcn,&dyngen->posteventfcn,&dyngen->posteventdgamma);CHKERRQ(ierr);
  } else { 
    dyngen->numMonitors = 0;
  }
  PetscFunctionReturn(0);
}

/*
  DYNGenModelSetStatus - Sets the generator status (1 = ON, 0 = OFF)

  Inputs:
+ dyngen - the DYNGenModel object
- status - the generator status
*/
PetscErrorCode DYNGenModelSetStatus(DYNGenModel dyngen,PetscInt status)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PSGENSetStatus(dyngen->psgen,status);CHKERRQ(ierr);
  if(status) dyngen->bus->ngenON++;
  else dyngen->bus->ngenON--;
  PetscFunctionReturn(0);
}

/*
  DYNGenModelDestroy - Destroys the dyngen model object.

  Input Parameters
. dyngen - the dynamic generator model

  Notes: 
   This routine is only called for dynamic simulation. It calls the destroy method on
   the child class to free its data
*/
PetscErrorCode DYNGenModelDestroy(DYNGenModel dyngen)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dyngen->ops.destroy)(dyngen);CHKERRQ(ierr);
  ierr = PetscFree(dyngen->eqtypes);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelRegister - Registers a dynamic generator model

  Input Parameters:
+ sname     - model name (string)
. createfunction  - the class constructor
- sizeofstruct  - size of object (obtained with sizeof())
*/
PetscErrorCode DYNGenModelRegister(const char sname[],PetscErrorCode (*createfunction)(DYNGenModel),PetscInt sizeofstruct)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscFunctionBegin;
  for(i=0; i < ngenmodelsregistered;i++) {
    PetscBool match;
    ierr = PetscStrcmp(DYNGenModelList[i].name,sname,&match);
    if(match) PetscFunctionReturn(0);
  }
  i = ngenmodelsregistered;
  ierr = PetscStrcpy(DYNGenModelList[i].name,sname);CHKERRQ(ierr);
  DYNGenModelList[i].create = createfunction;
  DYNGenModelList[i].sizeofstruct = sizeofstruct;
  ngenmodelsregistered++;
  PetscFunctionReturn(0);
}

#include "dyngenrou.h"
#include "dynpvd1.h"
#include "dyncv.h"
extern PetscErrorCode DYNGenModelCreate_Genrou(DYNGenModel);
extern PetscErrorCode DYNGenModelCreate_Pvd1(DYNGenModel);
extern PetscErrorCode DYNGenModelCreate_Cv(DYNGenModel);

/*
  DYNGenModelRegisterAll - Registers all built-in dynamic generator models

*/
PetscErrorCode DYNGenModelRegisterAll(void)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(!DYNGenModelsRegisterAllCalled) DYNGenModelsRegisterAllCalled = PETSC_TRUE;

  ierr = DYNGenModelRegister("GENROU",DYNGenModelCreate_Genrou,sizeof(struct _p_DYNGenrou));CHKERRQ(ierr);
  ierr = DYNGenModelRegister("PVD1",DYNGenModelCreate_Pvd1,sizeof(struct _p_DYNPvd1));CHKERRQ(ierr);
  ierr = DYNGenModelRegister("CV",DYNGenModelCreate_Cv,sizeof(struct _p_DYNCv));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetDynExc - Returns the exciter associated with this generator

  Input Parameters:
. dyngen - the DYNGenModel object

  Output Parameters:
. dynexc - the dynexc object
*/
PetscErrorCode DYNGenModelGetDynExc(DYNGenModel dyngen,DYNExcModel *dynexc)
{
  PetscFunctionBegin;
  *dynexc = dyngen->dynexc;
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetDynTurbgov - Returns the turbine governor associated with this generator

  Input Parameters:
. dyngen - the DYNGenModel object

  Output Parameters:
. dynturbgov - the dynturbgov object
*/
PetscErrorCode DYNGenModelGetDynTurbgov(DYNGenModel dyngen,DYNTurbgovModel *dynturbgov)
{
  PetscFunctionBegin;
  *dynturbgov = dyngen->dynturbgov;
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetModelType - Parses a line in dyr file and reads the model type (string identifier). If the line does not have
                            the registered generator model then the type is set to NONE.

  Input Parameters:
. line - the line to parse

  Output Parameters:
. modeltype - type of the model (string name of the model)
*/
PetscErrorCode DYNGenModelReadModelType(char *line, char gentype[])
{
  PetscInt i;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for(i=0; i < ngenmodelsregistered; i++) {
    if(strstr(line,DYNGenModelList[i].name) != NULL) {
      ierr = PetscStrcpy(gentype,DYNGenModelList[i].name);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
  }

  ierr = PetscStrcpy(gentype,DYNGENNONE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetSpeedDeviationLocation - Gets the location of the speed deviation variable for the generator model relative to the bus variable starting location

  Input Parameters:
. dyngen - the generator model

  Output Parameters:
. loc - the location of the speed deviation

  Notes:
  The location of the speed deviation is relative to the bus variables. So, one should first obtain
  the first variable location for the generator model implementation using DYNGenModelGetFirstVariableLocation() 
  and then add the location of the dw variable to return loc. For example, if dw is the 3rd variable for this
  generator model implementation then the following code should be used in the implementation
  ierr = DYNGenModelGetFirstVariable(dyngen,&firstvarloc);
  *loc = firstvarloc + 2;
*/
PetscErrorCode DYNGenModelGetSpeedDeviationLocation(DYNGenModel dyngen, PetscInt *loc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dyngen->ops.getspeeddeviationlocation)(dyngen,loc);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetSpeedDeviationGlobalLocation - Gets the global location (location in the global vector) of the field voltage for the exciter model

  Input Parameters:
. dyngen - the generator model

  Output Parameters:
. locglob - the global location of the speed deviation variable
*/
PetscErrorCode DYNGenModelGetSpeedDeviationGlobalLocation(DYNGenModel dyngen, PetscInt *locglob)
{
  PetscErrorCode ierr;
  PetscInt       dwloc;

  PetscFunctionBegin;

  /* Get the location of dw location relative to the bus */
  ierr = DYNGenModelGetSpeedDeviationLocation(dyngen,&dwloc);CHKERRQ(ierr);
  *locglob = dyngen->bus->startlocglob + dwloc;

  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetSpeedDeviation - Returns the speed deviation and its time derivative
                                 (optional) of the generator at the current time.

  Input Parameters:
+ dyngen - the DYNGenModel object
. t      - the current time
. xdyn   - the state variables for this generator at time t

  Output Parameters:
+ spd     - the generator speed deviation
- dspd_dt - time-derivative of generator speed deviation (set to NULL if not required)
*/
PetscErrorCode DYNGenModelGetSpeedDeviation(DYNGenModel dyngen,PetscReal t, const PetscScalar* xdyn,PetscScalar* spd, PetscScalar* dspd_dt)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dyngen->ops.getspeeddeviation)(dyngen,t,xdyn,spd,dspd_dt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetFrequency - Returns the frequency of the generator at the current time.

  Input Parameters:
+ dyngen - the DYNGenModel object
. t      - the current time
. xdyn   - the state variables for this generator at time t

  Output Parameters:
. frequency   - the generator frequency in Hz
*/
PetscErrorCode DYNGenModelGetFrequency(DYNGenModel dyngen,PetscReal t, const PetscScalar* xdyn,PetscScalar* frequency)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dyngen->ops.getfrequency)(dyngen,t,xdyn,frequency);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelSetFrequency - Set the frequency of the generator at the current time.

  Input Parameters:
+ dyngen - the DYNGenModel object
. t      - the current time
. xdyn   - the state variables for this generator at time t
. frequency = the value to be set

  Output Parameters:
*/
PetscErrorCode DYNGenModelSetFrequency(DYNGenModel dyngen,PetscScalar* xdyn,PetscScalar frequency)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dyngen->ops.setfrequency)(dyngen,xdyn,frequency);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetFrequencyLimits - Returns the min. and max. limits for the frequency of the generator

  Input Parameters:
. dyngen - the DYNGenModel object


  Output Parameters:
+ fmax   - the max. generator frequency allowed in Hz.
. fmin   - the min. generator frequency allowed in Hz.
*/
PetscErrorCode DYNGenModelGetFrequencyLimits(DYNGenModel dyngen,PetscScalar* fmax,PetscScalar* fmin)
{

  PetscFunctionBegin;
  *fmax = dyngen->freq_max;
  *fmin = dyngen->freq_min;
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetdFreqdState - Returns the partial derivative of the generator frequency w.r.t the machine state that governs it

  Input Parameters:
+ dyngen - the DYNGenModel object
. t      - the current time
. xdyn   - the state variables for this generator at time t

  Output Parameters:
+ dfreqdstate   - dfreq_dstate
- stateloc      - location w.r.t. to the bus variables
*/
PetscErrorCode DYNGenModelGetdFreqdState(DYNGenModel dyngen,PetscReal t, const PetscScalar* xdyn,PetscScalar* dfreqdstate,PetscInt *stateloc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dyngen->ops.getdfreqdstate)(dyngen,t,xdyn,dfreqdstate,stateloc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*
  DYNGenModelGetBusnumID - Gets the bus number and generator ID for this dynamic generator model
                            from the data file. Calls the underlying implementation
  Input Parameters:
. dyngen - the dynamic generator model object

  Output Parameters:
+ busnum - bus number
. genid  - the generator ID
*/
PetscErrorCode DYNGenModelGetBusnumID(DYNGenModel dyngen, PetscInt *busnum,char **genid)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = (*dyngen->ops.getbusnumid)(dyngen,busnum,genid);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetSizeofType - Gets the size of the implementation type (obtained using sizeof())

  Input Parameters:
. dyngen - the DYNGenModel object

  Output Parameters:
. size - the size of the derived class object
*/
PetscErrorCode DYNGenModelGetSizeof(DYNGenModel dyngen,PetscInt *size)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dyngen->ops.getsizeof)(dyngen,size);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetNvar - Returns the number of variables for this dynamic generator model

  Input Parameters:
. dyngen - the DYNGenModel object

  Output Parameters:
. Nvar - number of variables (dofs) for this generator model
*/
PetscErrorCode DYNGenModelGetNvar(DYNGenModel dyngen,PetscInt *Nvar)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dyngen->ops.getnvar)(dyngen,Nvar);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetFirstVariableLocation - Gets the location for the first variable for the
                        dyngen model from the bus variables array

  Input Parameters:
. dyngen - the DYNGenModel object

  Output Parameters
. loc - the location for the first variable 
*/
PetscErrorCode DYNGenModelGetFirstVariableLocation(DYNGenModel dyngen,PetscInt *loc)
{
  PetscFunctionBegin;
  *loc = dyngen->startloc;
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetFieldCurrent - Returns the field current Ifd

  Input Parameters:
. dyngen - the DYNGenModel object

  Output Parameters:
. Ifd - Field current Ifd
*/
PetscErrorCode DYNGenModelGetFieldCurrent(DYNGenModel dyngen,PetscScalar *x,PetscScalar *Ifd)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dyngen->ops.getfieldcurrent)(dyngen,x,Ifd);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
           
/*
  DYNGenModelGetInitialFieldVoltage - Returns the field voltage at t=t0

  Input Parameters:
. dyngen - the DYNGenModel object

  Output Parameters:
. Efd0 - Initial field voltage at t=t0
*/
PetscErrorCode DYNGenModelGetInitialFieldVoltage(DYNGenModel dyngen,PetscScalar *Efd0)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dyngen->ops.getinitialfieldvoltage)(dyngen,Efd0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetInitialMechanicalPower - Returns the mechanical power at t=t0

  Input Parameters:
. dyngen - the DYNGenModel object

  Output Parameters:
. Pmech0 - Initial mechanical power at t=t0
*/
PetscErrorCode DYNGenModelGetInitialMechanicalPower(DYNGenModel dyngen,PetscScalar *Pmech0)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dyngen->ops.getinitialmechanicalpower)(dyngen,Pmech0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetCurrent - Returns the current injection

  Input Parameters:
. dyngen - the DYNGenModel object

  Output Parameters:
+ IGD - Real part of current injection
- IGQ - Reactive part of current injection
*/
PetscErrorCode DYNGenModelGetCurrent(DYNGenModel dyngen,PetscScalar VD, PetscScalar VQ, PetscScalar *xdyn, PetscScalar *IGD, PetscScalar *IGQ)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dyngen->ops.getcurrent)(dyngen,VD,VQ,xdyn,IGD,IGQ);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelSetVoltage - Sets the voltage for this generator model

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
PetscErrorCode DYNGenModelSetVoltage(DYNGenModel dyngen,PetscScalar VD, PetscScalar VQ)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dyngen->ops.setvoltage)(dyngen,VD,VQ);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*
  DYNGenModelSetType - Sets the type of the dynamic generator model

  Input Parameters:
+ dyngen - the dynamic generator model object
- modeltype - type of the model
*/
PetscErrorCode DYNGenModelSetType(DYNGenModel dyngen,DYNGenModelType modeltype)
{
  PetscErrorCode ierr,(*r)(DYNGenModel)=NULL;
  PetscInt       i;
  PetscFunctionBegin;
  PetscBool match;
  for(i=0;i < ngenmodelsregistered;i++) {
    ierr = PetscStrcmp(DYNGenModelList[i].name,modeltype,&match);CHKERRQ(ierr);
    if(match) {
      r = DYNGenModelList[i].create;
      break;
    }
  }

  dyngen->dynexc = 0;
  dyngen->dynturbgov = 0;
  if(!r) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown type for DYNGenModel %s",modeltype);
  /* Null the function pointers */
  dyngen->ops.readdata = 0;
  dyngen->ops.getnvar = 0;
  dyngen->ops.destroy = 0;
  dyngen->ops.getsizeof = 0;
  dyngen->ops.setinitialconditions = 0;
  dyngen->ops.setinitialconditionsp = 0;
  dyngen->ops.daerhsfunction=0;
  dyngen->ops.daerhsjacobian=0;
  dyngen->ops.getequationtypes=0;
  dyngen->ops.getbusnumid=0;
  dyngen->ops.getfrequency=0;
  dyngen->ops.getspeeddeviation=0;
  dyngen->ops.getspeeddeviationlocation=0;
  dyngen->ops.getinitialfieldvoltage=0;
  dyngen->ops.getinitialmechanicalpower=0;
  dyngen->ops.setevent=0;
  dyngen->ops.getcurrent=0;
  dyngen->ops.setvoltage=0;

  /* Copy the type name */
  ierr = PetscStrcpy(dyngen->type,modeltype);CHKERRQ(ierr);

  /* Call the underlying implementation constructor */
  ierr = (*r)(dyngen);CHKERRQ(ierr);
  /* Set the numnber of variables for this dynamic generator model */
  ierr = DYNGenModelGetNvar(dyngen,&dyngen->nvar);CHKERRQ(ierr);

  dyngen->startloc = -1;
  dyngen->ndiff    =  0;
  dyngen->nalg     =  0;
  dyngen->eqtypes  =  0;
  PetscFunctionReturn(0);
}

/*
  DYNGenModelCopy - Copies the dyngen model

  Input Parameters:
. dyngenin - the dyngen model to copy from

  Output Parameters:
. dyngenout - the dyngen model to copy to

  Notes:
   If a new method is created for the dyngenmodel then that should be copied too.
   The underlying implementation (dyngen->data) data pointer is copied.
*/
PetscErrorCode DYNGenModelCopy(DYNGenModel dyngenout,DYNGenModel dyngenin)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscMemcpy(dyngenout,dyngenin,sizeof(struct _p_DYNGenModel));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelGetModelData - Returns the model data associated with the particular implementation type

  Input Parameters:
. dyngen - the DYNGenModel object

  Output Parameters:
. modeldata - the model data struct
*/
PetscErrorCode DYNGenModelGetModelData(DYNGenModel dyngen,void **modeldata)
{
  PetscFunctionBegin;
  *modeldata = (void*)dyngen->data;
  PetscFunctionReturn(0);
}

/*
  DYNGenModelReadData - Parses the given line and reads the data for the dynamic generator
  model

  Input Parameters:
+ dyngen - The generator base model
. line   - the file from the line
. mbase  - the machine base
- sbase  - the system base

  Notes: 
    The type of generator model should be set prior to calling this routine via DYNGenModelSetType()

    The conversion from machine base to system base,if applicable, should be done in this rouine
*/
PetscErrorCode DYNGenModelReadData(DYNGenModel dyngen,char *line,PetscScalar mbase,PetscScalar sbase)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dyngen->ops.readdata)(dyngen,line,mbase,sbase);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelSetInitialConditions - Sets the initial conditions for this generator model

  Input Parameters:
+ dyngen - the dynamic generator model object
. Pg     - the real generator power output
. Qg     - the reactive generator power output
. VD     - real component of the bus voltage
- VQ     - imaginary component of the bus voltage

  Output Parameters:
. x   - the array of the variables for the bus on which this generator is incident.

  Notes: 
   The initial conditions are populated in the array x

   Other constants, parameters can be also set during this function call in
   its implementation struct (for e.g. setting mechanical power Pm etc.)

   The locations to insert the variables for this generator should be obtained from
   DYNGenModelGetFirstVariableLocation()
   
   x contains all the variables for the generator bus
*/
PetscErrorCode DYNGenModelSetInitialConditions(DYNGenModel dyngen,PetscScalar Pg,PetscScalar Qg,PetscScalar VD,PetscScalar VQ,PetscScalar *x)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dyngen->ops.setinitialconditions)(dyngen,Pg,Qg,VD,VQ,x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DYNGenModelSetInitialConditionsP(DYNGenModel dyngen,PetscScalar PG,PetscScalar QG,PetscScalar VA,PetscScalar VM,PetscInt PGloc,PetscInt QGloc,PetscInt VAloc,PetscInt VMloc,Mat icP,PetscInt dynlocglob)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr= (*dyngen->ops.setinitialconditionsp)(dyngen,PG,QG,VA,VM,PGloc,QGloc,VAloc,VMloc,icP,dynlocglob);
  PetscFunctionReturn(0);
}

PetscErrorCode DYNGenModelSetInitialFieldVoltageDiff(DYNGenModel dyngen,PetscScalar PG,PetscScalar QG,PetscScalar VA,PetscScalar VM,PetscScalar *Efd_PG,PetscScalar *Efd_QG,PetscScalar *Efd_VA,PetscScalar *Efd_VM)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr= (*dyngen->ops.setinitialfieldvoltagediff)(dyngen,PG,QG,VA,VM,Efd_PG,Efd_QG,Efd_VA,Efd_VM);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelDAERHSJacobian - Computes the RHS Jacobian of the DAE equations

  Input Parameters:
+ dyngen     - the dynamic generator model
. t          - the current time
. VD         - real-part of complex bus voltage
. VQ         - imaginary-part of complex bus voltage
. xdyngen    - generator variables
. dynlocglob - starting location of generator variables in the global vector 
. V_loc      - global location of VD and VQ V_loc[0] = VD_loc V_loc[1] = VQ_loc
- I_loc      - global location of ID and IQ I_loc[0] = ID_loc I_loc[1] = IQ_loc

  Output Parameters:
. J          - the Jacobian matrix

  Notes:
   For each bus, the network equations are ordered as [IQ_loc;ID_loc]
*/ 
PetscErrorCode DYNGenModelDAERHSJacobian(DYNGenModel dyngen,Mat J,PetscReal t,PetscScalar VD,PetscScalar VQ,PetscScalar *xdyngen,PetscInt dynlocglob,PetscInt V_loc[],PetscInt I_loc[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dyngen->ops.daerhsjacobian)(dyngen,J,t,VD,VQ,xdyngen,dynlocglob,V_loc,I_loc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelDAERHSJacobianP - Sets the Jacobian of the DAE w.r.t. parameters

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
PetscErrorCode DYNGenModelDAERHSJacobianP(DYNGenModel dyngen, PetscReal t,const PetscScalar *x,Mat jacP,PetscInt dynlocglob,PetscInt Valoc,PetscInt Vmloc,PetscInt Pgloc,PetscInt Qgloc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dyngen->ops.daerhsjacobianp)(dyngen,t,x,jacP,dynlocglob,Valoc,Vmloc,Pgloc,Qgloc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelDAEFWDRHSJacobianP - Sets the Jacobian of the DAE w.r.t. parameters

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
PetscErrorCode DYNGenModelDAEFWDRHSJacobianP(DYNGenModel dyngen, PetscReal t,const PetscScalar *x,Vec *jacP,PetscInt dynlocglob,PetscInt Valoc,PetscInt Vmloc,PetscInt Pgloc,PetscInt Qgloc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dyngen->ops.daefwdrhsjacobianp)(dyngen,t,x,jacP,dynlocglob,Valoc,Vmloc,Pgloc,Qgloc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGenModelDAERHSFunction - Computes the RHS of the DAE function for the generator model [f(x,y);g(x,y)] 
                              and also returns the generator currents in network refernce frame

  Input Parameters:
+ dyngen - the dynamic generator model
. t      - the current time
. x      - array of the variables for the bus on which this generator is incident
. VD     - real-part of bus voltage
. VQ     - imaginary part of bus voltage

  Output Parameters:
+ f      - array of rhs of DAE equations for this generator
. IGD    - real-part of generator current in network reference frame
. IGQ    - imaginary-part of generator current in network reference frame

  Notes:
   The locations for the variables (and the corresponding locations to insert entries in f)
   for this generator model should be obtained by DYNGenModelGetFirstVariableLocation()

*/
PetscErrorCode DYNGenModelDAERHSFunction(DYNGenModel dyngen,PetscReal t,PetscScalar VD, PetscScalar VQ, PetscScalar *x,PetscScalar *f,PetscScalar *IGD,PetscScalar *IGQ)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dyngen->ops.daerhsfunction)(dyngen,t,VD,VQ,x,f,IGD,IGQ);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
			 
/*
  DYNGenModelGetEquationTypes - Gets the number and indices of differential and algebraic equations for the generator model
                                                   
  Input Parameters:
. dyngen - the DYNGenModel object

  Output Parameters:
+ ndiff - number of differential equations
. nalg  - number of algebraic equations
- eqtype - an array of size nvar that has
            eqtype[i] = DIFF_EQ if equation i is a differential
            eqtype[i] = ALG_EQ if equation  i is an algebraic

  NOTES: The user does not need to create the array vartype
*/
PetscErrorCode DYNGenModelGetEquationTypes(DYNGenModel dyngen,PetscInt *ndiff,PetscInt *nalg, PetscInt *eqtype)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dyngen->ops.getequationtypes)(dyngen,ndiff,nalg,eqtype);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
