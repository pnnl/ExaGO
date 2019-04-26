#include <private/dynstabmodelsimpl.h>
#include <private/psimpl.h>

struct _p_DYNStabModelList DYNStabModelList[DYNSTABMODELSMAX];

PetscInt nstabmodelsregistered=0;
PetscBool DYNStabModelsRegisterAllCalled = PETSC_FALSE;

/*
  DYNStabModelEventMonitor - Sets the event info for the stab model

  Input Parameters
. dynstab - the dynamic stab model

*/
PetscErrorCode DYNStabModelSetEventMonitor(DYNStabModel dynstab)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if(dynstab->ops.setevent) {
    ierr = (*dynstab->ops.setevent)(dynstab,&dynstab->numMonitors,dynstab->direction,dynstab->terminate,&dynstab->eventfcn,&dynstab->posteventfcn,&dynstab->posteventdgamma);CHKERRQ(ierr);
  } else { 
    dynstab->numMonitors = 0;
  }
  PetscFunctionReturn(0);
}

/*
  DYNStabModelDestroy - Destroys the dynstab model object.

  Input Parameters
. dynstab - the dynamic stab model

  Notes: 
   This routine is only called for dynamic simulation. It calls the destroy method on
   the child class to free its data
*/
PetscErrorCode DYNStabModelDestroy(DYNStabModel dynstab)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynstab->ops.destroy)(dynstab);CHKERRQ(ierr);
  ierr = PetscFree(dynstab->eqtypes);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNStabModelRegister - Registers a dynamic stab model

  Input Parameters:
+ sname     - model name (string)
. createfunction  - the class constructor
- sizeofstruct  - size of object (obtained with sizeof())
*/
PetscErrorCode DYNStabModelRegister(const char sname[],PetscErrorCode (*createfunction)(DYNStabModel),PetscInt sizeofstruct)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscFunctionBegin;
  for(i=0; i < nstabmodelsregistered;i++) {
    PetscBool match;
    ierr = PetscStrcmp(DYNStabModelList[i].name,sname,&match);
    if(match) PetscFunctionReturn(0);
  }
  i = nstabmodelsregistered;
  ierr = PetscStrcpy(DYNStabModelList[i].name,sname);CHKERRQ(ierr);
  DYNStabModelList[i].create = createfunction;
  DYNStabModelList[i].sizeofstruct = sizeofstruct;
  nstabmodelsregistered++;
  PetscFunctionReturn(0);
}

#include "dynstab1.h"
extern PetscErrorCode DYNStabModelCreate_STAB1(DYNStabModel);

/*
  DYNStabModelRegisterAll - Registers all built-in stab models

*/
PetscErrorCode DYNStabModelRegisterAll(void)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(!DYNStabModelsRegisterAllCalled) DYNStabModelsRegisterAllCalled = PETSC_TRUE;

  ierr = DYNStabModelRegister("STAB1",DYNStabModelCreate_STAB1,sizeof(struct _p_DYNSTAB1));CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNStabModelGetVOTHSG - Returns the stabilizer signal for the exciter and optionally returns partial derivative w.r.t to
                                      the stabilizer variables and speed deviation

  Input Parameters:
+ dynstab - the dynamic stabilizer model
. t      - the current time
. xdyn   - array of the variables for this bus

  Output Parameters:
+ VOTHSG           - Stabilizer signal output
. dVOTHSGdXdyn_num - Number of partial derivatives of VOTHSG w.r.t. the dynamic variables at the bus 
. dVOTHSGdXdyn     - An array of partial derivatives of the stabilizer output signal VOTHSG w.r.t. the dynamic
                     variables at the bus
- dVOTHSGdXdyn_loc - Location of partial derivatives of the stabilizer output signal VOTHSG w.r.t. the dynamic
                     variables at the bus.  Need to get the first variable location using DYNStabModelGetFirstVariableLocation()
		     and then add the displacement.

  Notes - dVOTHSGdXdyn_num, dVOTHSGdXdyn and dVOTHSGdXdyn_loc are optionally returned if the passed-in arrays are not NULL.
          
*/
PetscErrorCode DYNStabModelGetVOTHSG(DYNStabModel dynstab,PetscReal t, PetscScalar *xdyn,PetscScalar *VOTHSG,PetscInt *dVOTHSGdXdyn_num,PetscScalar dVOTHSGdXdyn[],PetscInt dVOTHSGdXdyn_loc[])
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  ierr = (*dynstab->ops.getvothsg)(dynstab,t,xdyn,VOTHSG,dVOTHSGdXdyn_num,dVOTHSGdXdyn,dVOTHSGdXdyn_loc);CHKERRQ(ierr);
  if(dVOTHSGdXdyn_num) {
    /* Conversion to global locations */
    for(i=0; i < *dVOTHSGdXdyn_num; i++) {
      dVOTHSGdXdyn_loc[i] = dynstab->dyngen->bus->startlocglob + dVOTHSGdXdyn_loc[i];
    }
  }
  PetscFunctionReturn(0);
}

/*
  DYNStabModelGetFirstVariableLocation - Gets the location for the first variable for the
                        dynstab model from the bus variables array

  Input Parameters:
. dynstab - the DYNStabModel object

  Output Parameters
. loc - the location for the first variable 
*/
PetscErrorCode DYNStabModelGetFirstVariableLocation(DYNStabModel dynstab,PetscInt *loc)
{
  PetscFunctionBegin;
  *loc = dynstab->startloc;
  PetscFunctionReturn(0);
}


/*
  DYNStabModelReadModelType - Parses a line in dyr file and reads the model type (string identifier). If the line does not have
                            the registered stab model then the type is set to NONE.

  Input Parameters:
. line - the line to parse

  Output Parameters:
. modeltype - type of the model (string name of the model)
*/
PetscErrorCode DYNStabModelReadModelType(char *line, char gentype[])
{
  PetscInt i;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for(i=0; i < nstabmodelsregistered; i++) {
    if(strstr(line,DYNStabModelList[i].name) != NULL) {
      ierr = PetscStrcpy(gentype,DYNStabModelList[i].name);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
  }

  ierr = PetscStrcpy(gentype,DYNSTABNONE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNStabModelGetBusnumID - Gets the bus number and stab ID for this dynamic stab model
                            from the data file. Calls the underlying implementation
  Input Parameters:
. dynstab - the dynamic stab model object

  Output Parameters:
+ busnum - bus number
. genid  - the stab ID
*/
PetscErrorCode DYNStabModelGetBusnumID(DYNStabModel dynstab, PetscInt *busnum,char **genid)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = (*dynstab->ops.getbusnumid)(dynstab,busnum,genid);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNStabModelGetSizeofType - Gets the size of the implementation type (obtained using sizeof())

  Input Parameters:
. dynstab - the DYNStabModel object

  Output Parameters:
. size - the size of the derived class object
*/
PetscErrorCode DYNStabModelGetSizeof(DYNStabModel dynstab,PetscInt *size)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynstab->ops.getsizeof)(dynstab,size);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNStabModelGetNvar - Returns the number of variables for this dynamic stab model

  Input Parameters:
. dynstab - the DYNStabModel object

  Output Parameters:
. Nvar - number of variables (dofs) for this stab model
*/
PetscErrorCode DYNStabModelGetNvar(DYNStabModel dynstab,PetscInt *Nvar)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynstab->ops.getnvar)(dynstab,Nvar);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNStabModelSetType - Sets the type of the dynamic stab model

  Input Parameters:
+ dynstab - the dynamic stab model object
- modeltype - type of the model
*/
PetscErrorCode DYNStabModelSetType(DYNStabModel dynstab,DYNStabModelType modeltype)
{
  PetscErrorCode ierr,(*r)(DYNStabModel)=NULL;
  PetscInt       i;
  PetscFunctionBegin;
  PetscBool match;
  for(i=0;i < nstabmodelsregistered;i++) {
    ierr = PetscStrcmp(DYNStabModelList[i].name,modeltype,&match);CHKERRQ(ierr);
    if(match) {
      r = DYNStabModelList[i].create;
      break;
    }
  }

  dynstab->dyngen = 0;
  if(!r) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown type for DYNStabModel %s",modeltype);
  /* Null the function pointers */
  dynstab->ops.readdata = 0;
  dynstab->ops.getnvar = 0;
  dynstab->ops.destroy = 0;
  dynstab->ops.getsizeof = 0;
  dynstab->ops.setinitialconditions = 0;
  dynstab->ops.getvothsg = 0;
  dynstab->ops.daerhsfunction=0;
  dynstab->ops.daerhsjacobian=0;
  dynstab->ops.getequationtypes=0;
  dynstab->ops.getbusnumid=0;
  dynstab->ops.setevent=0;

  /* Copy the type name */
  ierr = PetscStrcpy(dynstab->type,modeltype);CHKERRQ(ierr);

  /* Call the underlying implementation constructor */
  ierr = (*r)(dynstab);CHKERRQ(ierr);
  /* Set the numnber of variables for this dynamic stab model */
  ierr = DYNStabModelGetNvar(dynstab,&dynstab->nvar);CHKERRQ(ierr);

  dynstab->startloc    = -1;
  dynstab->ndiff       =  0;
  dynstab->nalg        =  0;
  dynstab->eqtypes     =  0;
  dynstab->numMonitors = 0;
  dynstab->eventfcn    = NULL;
  dynstab->posteventfcn = NULL;
  dynstab->posteventdgamma = NULL;
  PetscFunctionReturn(0);
}

/*
  DYNStabModelCopy - Copies the dynstab model

  Input Parameters:
. dynstabin - the dynstab model to copy from

  Output Parameters:
. dynstabout - the dynstab model to copy to

  Notes:
   If a new method is created for the dynstabmodel then that should be copied too.
   The underlying implementation (dynstab->data) data pointer is copied.
*/
PetscErrorCode DYNStabModelCopy(DYNStabModel dynstabout,DYNStabModel dynstabin)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscMemcpy(dynstabout,dynstabin,sizeof(struct _p_DYNStabModel));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNStabModelGetModelData - Returns the model data associated with the particular implementation type

  Input Parameters:
. dynstab - the DYNStabModel object

  Output Parameters:
. modeldata - the model data struct
*/
PetscErrorCode DYNStabModelGetModelData(DYNStabModel dynstab,void **modeldata)
{
  PetscFunctionBegin;
  *modeldata = (void*)dynstab->data;
  PetscFunctionReturn(0);
}

/*
  DYNStabModelGetDynGen - Returns the generator associated with this stab

  Input Parameters:
. dynstab - the DYNStabModel object

  Output Parameters:
. dyngen - the dyngen object
*/
PetscErrorCode DYNStabModelGetDynGen(DYNStabModel dynstab,DYNGenModel *dyngen)
{
  PetscFunctionBegin;
  *dyngen = dynstab->dyngen;
  PetscFunctionReturn(0);
}

/*
  DYNStabModelReadData - Parses the given line and reads the data for the dynamic stab
  model

  Input Parameters:
+ dynstab - The stab base model
- line   - the file from the line

  Notes: The type of stab model should be set prior to this call via DYNStabModelSetType()
*/
PetscErrorCode DYNStabModelReadData(DYNStabModel dynstab,char *line)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynstab->ops.readdata)(dynstab,line);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNStabModelSetInitialConditions - Sets the initial conditions for this stab model

  Input Parameters:
+ dynstab - the dynamic stab model object
. VD     - real component of the bus voltage
- VQ     - imaginary component of the bus voltage

  Output Parameters:
. xdyn   - the initial conditions for this stab model

  Notes: 
   The initial conditions are populated in the xdyn array

   Other constants, parameters can be also set during this function call in
   its implementation struct (for e.g. setting reference voltage Vref and others)
*/
PetscErrorCode DYNStabModelSetInitialConditions(DYNStabModel dynstab,PetscScalar VD,PetscScalar VQ,PetscScalar *xdyn)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynstab->ops.setinitialconditions)(dynstab,VD,VQ,xdyn);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNStabModelDAERHSJacobian - Computes the RHS Jacobian of the DAE equations

  Input Parameters:
+ dynstab     - the dynamic stab model
. t          - the current time
. VD         - real-part of complex bus voltage
. VQ         - imaginary-part of complex bus voltage
. xdynstab    - stab variables
. dynlocglob - starting location of stab variables in the global vector 
- V_loc      - global location of VD and VQ V_loc[0] = VD_loc V_loc[1] = VQ_loc


  Output Parameters:
. J          - the Jacobian matrix

  Notes:
   For each bus, the network equations are ordered as [IQ_loc;ID_loc]
*/ 
PetscErrorCode DYNStabModelDAERHSJacobian(DYNStabModel dynstab,Mat J,PetscReal t,PetscScalar VD,PetscScalar VQ,PetscScalar *xdynstab,PetscInt dynlocglob,PetscInt V_loc[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynstab->ops.daerhsjacobian)(dynstab,J,t,VD,VQ,xdynstab,dynlocglob,V_loc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNStabModelDAERHSJacobianP - Computes the RHS Jacobian of the DAE equations w.r.t. parameters

  Input Parameters:
+ dynstab     - the dynamic stab model
. t          - the current time
. xdynstab    - stab variables
. dynlocglob - starting location of stab variables in the global vector 
- PGloc      - location of PG (generator real power) in the parameter vector

  Output Parameters:
. jacP       - matrix of partial derivatives of DAE equations w.r.t parameter PG

*/ 
PetscErrorCode DYNStabModelDAERHSJacobianP(DYNStabModel dynstab,PetscReal t,const PetscScalar *xdynstab,Mat jacP,PetscInt dynlocglob,PetscInt PGloc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynstab->ops.daerhsjacobianp)(dynstab,t,xdynstab,jacP,dynlocglob,PGloc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNStabModelDAERHSFunction - Computes the RHS of the DAE function for the stab model [f(x,y);g(x,y)] 
                              and also returns the stab currents in network refernce frame

  Input Parameters:
+ dynstab - the dynamic stab model
. t      - the current time
. x      - array of variables for the bus on which the generator is incident
. VD     - real-part of bus voltage
. VQ     - imaginary part of bus voltage

  Output Parameters:
+ f      - the residual array in which the stab equation residuals are inserted. The location
           of the stab variables can be obtained using DYNStabModelGetFirstVariableLocation(). Using
           the first location, the residuals can be inserted as
           f[loc] = 1st stab equation residual
           f[loc+1] = 2nd ..
         ... and so on
*/
PetscErrorCode DYNStabModelDAERHSFunction(DYNStabModel dynstab,PetscReal t,PetscScalar VD, PetscScalar VQ, PetscScalar *x,PetscScalar *f)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynstab->ops.daerhsfunction)(dynstab,t,VD,VQ,x,f);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
			 
/*
  DYNStabModelGetEquationTypes - Gets the number and indices of differential and algebraic equations for the stab model
                                                   
  Input Parameters:
. dynstab - the DYNStabModel object

  Output Parameters:
+ ndiff - number of differential equations
. nalg  - number of algebraic equations
- eqtype - an array of size nvar that has
            eqtype[i] = DIFF_EQ if equation i is a differential
            eqtype[i] = ALG_EQ if equation  i is an algebraic

  NOTES: The user does not need to create the array vartype
*/
PetscErrorCode DYNStabModelGetEquationTypes(DYNStabModel dynstab,PetscInt *ndiff,PetscInt *nalg, PetscInt *eqtype)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynstab->ops.getequationtypes)(dynstab,ndiff,nalg,eqtype);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
