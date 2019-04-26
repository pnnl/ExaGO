#include <private/dynturbgovmodelsimpl.h>
#include <private/psimpl.h>

struct _p_DYNTurbgovModelList DYNTurbgovModelList[DYNTURBGOVMODELSMAX];

PetscInt nturbgovmodelsregistered=0;
PetscBool DYNTurbgovModelsRegisterAllCalled = PETSC_FALSE;

/*
  DYNTurbgovModelEventMonitor - Sets the event info for the turbgov model

  Input Parameters
. dynturbgov - the dynamic turbgov model

*/
PetscErrorCode DYNTurbgovModelSetEventMonitor(DYNTurbgovModel dynturbgov)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if(dynturbgov->ops.setevent) {
    ierr = (*dynturbgov->ops.setevent)(dynturbgov,&dynturbgov->numMonitors,dynturbgov->direction,dynturbgov->terminate,&dynturbgov->eventfcn,&dynturbgov->posteventfcn,&dynturbgov->posteventdgamma);CHKERRQ(ierr);
  } else { 
    dynturbgov->numMonitors = 0;
  }
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelDestroy - Destroys the dynturbgov model object.

  Input Parameters
. dynturbgov - the dynamic turbgov model

  Notes: 
   This routine is only called for dynamic simulation. It calls the destroy method on
   the child class to free its data
*/
PetscErrorCode DYNTurbgovModelDestroy(DYNTurbgovModel dynturbgov)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynturbgov->ops.destroy)(dynturbgov);CHKERRQ(ierr);
  ierr = PetscFree(dynturbgov->eqtypes);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelRegister - Registers a dynamic turbgov model

  Input Parameters:
+ sname     - model name (string)
. createfunction  - the class constructor
- sizeofstruct  - size of object (obtained with sizeof())
*/
PetscErrorCode DYNTurbgovModelRegister(const char sname[],PetscErrorCode (*createfunction)(DYNTurbgovModel),PetscInt sizeofstruct)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscFunctionBegin;
  for(i=0; i < nturbgovmodelsregistered;i++) {
    PetscBool match;
    ierr = PetscStrcmp(DYNTurbgovModelList[i].name,sname,&match);
    if(match) PetscFunctionReturn(0);
  }
  i = nturbgovmodelsregistered;
  ierr = PetscStrcpy(DYNTurbgovModelList[i].name,sname);CHKERRQ(ierr);
  DYNTurbgovModelList[i].create = createfunction;
  DYNTurbgovModelList[i].sizeofstruct = sizeofstruct;
  nturbgovmodelsregistered++;
  PetscFunctionReturn(0);
}

#include "dyntgov1.h"
extern PetscErrorCode DYNTurbgovModelCreate_TGOV1(DYNTurbgovModel);

/*
  DYNTurbgovModelRegisterAll - Registers all built-in turbgov models

*/
PetscErrorCode DYNTurbgovModelRegisterAll(void)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(!DYNTurbgovModelsRegisterAllCalled) DYNTurbgovModelsRegisterAllCalled = PETSC_TRUE;

  ierr = DYNTurbgovModelRegister("TGOV1",DYNTurbgovModelCreate_TGOV1,sizeof(struct _p_DYNTGOV1));CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelGetFirstVariableLocation - Gets the location for the first variable for the
                        dynturbgov model from the bus variables array

  Input Parameters:
. dynturbgov - the DYNTurbgovModel object

  Output Parameters
. loc - the location for the first variable 
*/
PetscErrorCode DYNTurbgovModelGetFirstVariableLocation(DYNTurbgovModel dynturbgov,PetscInt *loc)
{
  PetscFunctionBegin;
  *loc = dynturbgov->startloc;
  PetscFunctionReturn(0);
}


/*
  DYNTurbgovModelReadModelType - Parses a line in dyr file and reads the model type (string identifier). If the line does not have
                            the registered turbgov model then the type is set to NONE.

  Input Parameters:
. line - the line to parse

  Output Parameters:
. modeltype - type of the model (string name of the model)
*/
PetscErrorCode DYNTurbgovModelReadModelType(char *line, char gentype[])
{
  PetscInt i;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for(i=0; i < nturbgovmodelsregistered; i++) {
    if(strstr(line,DYNTurbgovModelList[i].name) != NULL) {
      ierr = PetscStrcpy(gentype,DYNTurbgovModelList[i].name);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
  }

  ierr = PetscStrcpy(gentype,DYNTURBGOVNONE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelGetBusnumID - Gets the bus number and turbgov ID for this dynamic turbgov model
                            from the data file. Calls the underlying implementation
  Input Parameters:
. dynturbgov - the dynamic turbgov model object

  Output Parameters:
+ busnum - bus number
. genid  - the turbgov ID
*/
PetscErrorCode DYNTurbgovModelGetBusnumID(DYNTurbgovModel dynturbgov, PetscInt *busnum,char **genid)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = (*dynturbgov->ops.getbusnumid)(dynturbgov,busnum,genid);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelGetSizeofType - Gets the size of the implementation type (obtained using sizeof())

  Input Parameters:
. dynturbgov - the DYNTurbgovModel object

  Output Parameters:
. size - the size of the derived class object
*/
PetscErrorCode DYNTurbgovModelGetSizeof(DYNTurbgovModel dynturbgov,PetscInt *size)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynturbgov->ops.getsizeof)(dynturbgov,size);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelGetNvar - Returns the number of variables for this dynamic turbgov model

  Input Parameters:
. dynturbgov - the DYNTurbgovModel object

  Output Parameters:
. Nvar - number of variables (dofs) for this turbgov model
*/
PetscErrorCode DYNTurbgovModelGetNvar(DYNTurbgovModel dynturbgov,PetscInt *Nvar)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynturbgov->ops.getnvar)(dynturbgov,Nvar);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelSetType - Sets the type of the dynamic turbgov model

  Input Parameters:
+ dynturbgov - the dynamic turbgov model object
- modeltype - type of the model
*/
PetscErrorCode DYNTurbgovModelSetType(DYNTurbgovModel dynturbgov,DYNTurbgovModelType modeltype)
{
  PetscErrorCode ierr,(*r)(DYNTurbgovModel)=NULL;
  PetscInt       i;
  PetscFunctionBegin;
  PetscBool match;
  for(i=0;i < nturbgovmodelsregistered;i++) {
    ierr = PetscStrcmp(DYNTurbgovModelList[i].name,modeltype,&match);CHKERRQ(ierr);
    if(match) {
      r = DYNTurbgovModelList[i].create;
      break;
    }
  }

  dynturbgov->dyngen = 0;
  if(!r) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown type for DYNTurbgovModel %s",modeltype);
  /* Null the function pointers */
  dynturbgov->ops.readdata = 0;
  dynturbgov->ops.getnvar = 0;
  dynturbgov->ops.destroy = 0;
  dynturbgov->ops.getsizeof = 0;
  dynturbgov->ops.setinitialconditions = 0;
  dynturbgov->ops.getmechanicalpower = 0;
  dynturbgov->ops.daerhsfunction=0;
  dynturbgov->ops.daerhsjacobian=0;
  dynturbgov->ops.getequationtypes=0;
  dynturbgov->ops.getbusnumid=0;
  dynturbgov->ops.setevent=0;

  /* Copy the type name */
  ierr = PetscStrcpy(dynturbgov->type,modeltype);CHKERRQ(ierr);

  /* Call the underlying implementation constructor */
  ierr = (*r)(dynturbgov);CHKERRQ(ierr);
  /* Set the numnber of variables for this dynamic turbgov model */
  ierr = DYNTurbgovModelGetNvar(dynturbgov,&dynturbgov->nvar);CHKERRQ(ierr);

  dynturbgov->startloc    = -1;
  dynturbgov->ndiff       =  0;
  dynturbgov->nalg        =  0;
  dynturbgov->eqtypes     =  0;
  dynturbgov->numMonitors = 0;
  dynturbgov->eventfcn    = NULL;
  dynturbgov->posteventfcn = NULL;
  dynturbgov->posteventdgamma = NULL;
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelCopy - Copies the dynturbgov model

  Input Parameters:
. dynturbgovin - the dynturbgov model to copy from

  Output Parameters:
. dynturbgovout - the dynturbgov model to copy to

  Notes:
   If a new method is created for the dynturbgovmodel then that should be copied too.
   The underlying implementation (dynturbgov->data) data pointer is copied.
*/
PetscErrorCode DYNTurbgovModelCopy(DYNTurbgovModel dynturbgovout,DYNTurbgovModel dynturbgovin)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscMemcpy(dynturbgovout,dynturbgovin,sizeof(struct _p_DYNTurbgovModel));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelGetModelData - Returns the model data associated with the particular implementation type

  Input Parameters:
. dynturbgov - the DYNTurbgovModel object

  Output Parameters:
. modeldata - the model data struct
*/
PetscErrorCode DYNTurbgovModelGetModelData(DYNTurbgovModel dynturbgov,void **modeldata)
{
  PetscFunctionBegin;
  *modeldata = (void*)dynturbgov->data;
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelGetDynGen - Returns the generator associated with this turbgov

  Input Parameters:
. dynturbgov - the DYNTurbgovModel object

  Output Parameters:
. dyngen - the dyngen object
*/
PetscErrorCode DYNTurbgovModelGetDynGen(DYNTurbgovModel dynturbgov,DYNGenModel *dyngen)
{
  PetscFunctionBegin;
  *dyngen = dynturbgov->dyngen;
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelReadData - Parses the given line and reads the data for the dynamic turbgov
  model

  Input Parameters:
+ dynturbgov - The turbgov base model
- line   - the file from the line
. mbase  - the machine base
- sbase  - the system base

  Notes: The type of turbgov model should be set prior to this call via DYNTurbgovModelSetType()
*/
PetscErrorCode DYNTurbgovModelReadData(DYNTurbgovModel dynturbgov,char *line,PetscScalar mbase,PetscScalar sbase)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynturbgov->ops.readdata)(dynturbgov,line,mbase,sbase);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelSetInitialConditions - Sets the initial conditions for this turbgov model

  Input Parameters:
+ dynturbgov - the dynamic turbgov model object
. VD     - real component of the bus voltage
- VQ     - imaginary component of the bus voltage

  Output Parameters:
. xdyn   - the initial conditions for this turbgov model

  Notes: 
   The initial conditions are populated in the xdyn array

   Other constants, parameters can be also set during this function call in
   its implementation struct (for e.g. setting reference voltage Vref and others)
*/
PetscErrorCode DYNTurbgovModelSetInitialConditions(DYNTurbgovModel dynturbgov,PetscScalar VD,PetscScalar VQ,PetscScalar *xdyn)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynturbgov->ops.setinitialconditions)(dynturbgov,VD,VQ,xdyn);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DYNTurbgovModelSetInitialConditionsP(DYNTurbgovModel dynturbgov,PetscScalar PG,PetscScalar QG,PetscScalar VA,PetscScalar VM,PetscInt PGloc,PetscInt QGloc,PetscInt VAloc,PetscInt VMloc,Mat ICp,PetscInt dynlocglob)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynturbgov->ops.setinitialconditionsp)(dynturbgov,PG,QG,VA,VM,PGloc,QGloc,VAloc,VMloc,ICp,dynlocglob);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelGetMechanicalPower - Returns the mechanical power input for the machine and its partial derivative w.r.t to
                                      the machine speed deviation

  Input Parameters:
+ dynturbgov - the dynamic turbgoviter model
. t      - the current time
. xdyn   - array of the variables for this bus

  Output Parameters:
+ Pmech      - Mechanical power input for the machine
. dPmechddw  - Partial derivative of Pmech w.r.t. machine speed deviation
. dPmechdXtgov_num - Number of partial derivatives computed.
. dPmechdXtgov - An array of partial derivatives of the mechanical power w.r.t. the number
                 of turbine governor states given by dPmechdXtgov_num
- dPmechdXtgov_loc - Global locations for the partial derivatives

  Notes - dPmechdXtgov_num, dPmechdXtgov and dPmechdXtgov_loc are optionally returned if the passed-in arrays are not NULL.
          
*/
PetscErrorCode DYNTurbgovModelGetMechanicalPower(DYNTurbgovModel dynturbgov,PetscReal t, PetscScalar *xdyn,PetscScalar *Pmech,PetscScalar *dPmechddw, PetscInt *dPmechdXtgov_num,PetscScalar dPmechdXtgov[],PetscInt dPmechdXtgov_loc[])
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  ierr = (*dynturbgov->ops.getmechanicalpower)(dynturbgov,t,xdyn,Pmech,dPmechddw,dPmechdXtgov_num,dPmechdXtgov,dPmechdXtgov_loc);CHKERRQ(ierr);
  if(dPmechdXtgov_num) {
    /* Conversion to global locations */
    for(i=0; i < *dPmechdXtgov_num; i++) {
      dPmechdXtgov_loc[i] = dynturbgov->dyngen->bus->startlocglob + dPmechdXtgov_loc[i];
    }
  }
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelDAERHSJacobian - Computes the RHS Jacobian of the DAE equations

  Input Parameters:
+ dynturbgov     - the dynamic turbgov model
. t          - the current time
. VD         - real-part of complex bus voltage
. VQ         - imaginary-part of complex bus voltage
. xdynturbgov    - turbgov variables
. dynlocglob - starting location of turbgov variables in the global vector 
- V_loc      - global location of VD and VQ V_loc[0] = VD_loc V_loc[1] = VQ_loc


  Output Parameters:
. J          - the Jacobian matrix

  Notes:
   For each bus, the network equations are ordered as [IQ_loc;ID_loc]
*/ 
PetscErrorCode DYNTurbgovModelDAERHSJacobian(DYNTurbgovModel dynturbgov,Mat J,PetscReal t,PetscScalar VD,PetscScalar VQ,PetscScalar *xdynturbgov,PetscInt dynlocglob,PetscInt V_loc[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynturbgov->ops.daerhsjacobian)(dynturbgov,J,t,VD,VQ,xdynturbgov,dynlocglob,V_loc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelDAERHSJacobianP - Computes the RHS Jacobian of the DAE equations w.r.t. parameters

  Input Parameters:
+ dynturbgov     - the dynamic turbgov model
. t          - the current time
. xdynturbgov    - turbgov variables
. dynlocglob - starting location of turbgov variables in the global vector 
- PGloc      - location of PG (generator real power) in the parameter vector

  Output Parameters:
. jacP       - matrix of partial derivatives of DAE equations w.r.t parameter PG

*/ 
PetscErrorCode DYNTurbgovModelDAERHSJacobianP(DYNTurbgovModel dynturbgov,PetscReal t,const PetscScalar *xdynturbgov,Mat jacP,PetscInt dynlocglob,PetscInt PGloc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynturbgov->ops.daerhsjacobianp)(dynturbgov,t,xdynturbgov,jacP,dynlocglob,PGloc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNTurbgovModelDAERHSFunction - Computes the RHS of the DAE function for the turbgov model [f(x,y);g(x,y)] 
                              and also returns the turbgov currents in network refernce frame

  Input Parameters:
+ dynturbgov - the dynamic turbgov model
. t      - the current time
. x      - array of variables for the bus on which the generator is incident
. VD     - real-part of bus voltage
. VQ     - imaginary part of bus voltage

  Output Parameters:
+ f      - the residual array in which the turbgov equation residuals are inserted. The location
           of the turbgov variables can be obtained using DYNTurbgovModelGetFirstVariableLocation(). Using
           the first location, the residuals can be inserted as
           f[loc] = 1st turbgov equation residual
           f[loc+1] = 2nd ..
         ... and so on
*/
PetscErrorCode DYNTurbgovModelDAERHSFunction(DYNTurbgovModel dynturbgov,PetscReal t,PetscScalar VD, PetscScalar VQ, PetscScalar *x,PetscScalar *f)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynturbgov->ops.daerhsfunction)(dynturbgov,t,VD,VQ,x,f);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
			 
/*
  DYNTurbgovModelGetEquationTypes - Gets the number and indices of differential and algebraic equations for the turbgov model
                                                   
  Input Parameters:
. dynturbgov - the DYNTurbgovModel object

  Output Parameters:
+ ndiff - number of differential equations
. nalg  - number of algebraic equations
- eqtype - an array of size nvar that has
            eqtype[i] = DIFF_EQ if equation i is a differential
            eqtype[i] = ALG_EQ if equation  i is an algebraic

  NOTES: The user does not need to create the array vartype
*/
PetscErrorCode DYNTurbgovModelGetEquationTypes(DYNTurbgovModel dynturbgov,PetscInt *ndiff,PetscInt *nalg, PetscInt *eqtype)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynturbgov->ops.getequationtypes)(dynturbgov,ndiff,nalg,eqtype);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
