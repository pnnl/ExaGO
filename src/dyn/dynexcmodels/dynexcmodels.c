#include <private/dynexcmodelsimpl.h>
#include <private/psimpl.h>

struct _p_DYNExcModelList DYNExcModelList[DYNEXCMODELSMAX];

PetscInt nexcmodelsregistered=0;
PetscBool DYNExcModelsRegisterAllCalled = PETSC_FALSE;

/*
  DYNExcModelEventMonitor - Sets the event info for the exciter model

  Input Parameters
. dynexc - the dynamic exciter model

*/
PetscErrorCode DYNExcModelSetEventMonitor(DYNExcModel dynexc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if(dynexc->ops.setevent) {
    ierr = (*dynexc->ops.setevent)(dynexc,&dynexc->numMonitors,dynexc->direction,dynexc->terminate,&dynexc->eventfcn,&dynexc->posteventfcn,&dynexc->posteventdgamma);CHKERRQ(ierr);
  } else { 
    dynexc->numMonitors = 0;
  }
  PetscFunctionReturn(0);
}

/*
  DYNExcModelDestroy - Destroys the dynexc model object.

  Input Parameters
. dynexc - the dynamic exciter model

  Notes: 
   This routine is only called for dynamic simulation. It calls the destroy method on
   the child class to free its data
*/
PetscErrorCode DYNExcModelDestroy(DYNExcModel dynexc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynexc->ops.destroy)(dynexc);CHKERRQ(ierr);
  ierr = PetscFree(dynexc->eqtypes);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNExcModelRegister - Registers a dynamic exciter model

  Input Parameters:
+ sname     - model name (string)
. createfunction  - the class constructor
- sizeofstruct  - size of object (obtained with sizeof())
*/
PetscErrorCode DYNExcModelRegister(const char sname[],PetscErrorCode (*createfunction)(DYNExcModel),PetscInt sizeofstruct)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscFunctionBegin;
  for(i=0; i < nexcmodelsregistered;i++) {
    PetscBool match;
    ierr = PetscStrcmp(DYNExcModelList[i].name,sname,&match);
    if(match) PetscFunctionReturn(0);
  }
  i = nexcmodelsregistered;
  ierr = PetscStrcpy(DYNExcModelList[i].name,sname);CHKERRQ(ierr);
  DYNExcModelList[i].create = createfunction;
  DYNExcModelList[i].sizeofstruct = sizeofstruct;
  nexcmodelsregistered++;
  PetscFunctionReturn(0);
}

#include "dynieeet1.h"
#include "dynexst1.h"
#include "dynsexs.h"
extern PetscErrorCode DYNExcModelCreate_IEEET1(DYNExcModel);
extern PetscErrorCode DYNExcModelCreate_EXST1(DYNExcModel);
extern PetscErrorCode DYNExcModelCreate_SEXS(DYNExcModel);

/*
  DYNExcModelRegisterAll - Registers all built-in exciter models

*/
PetscErrorCode DYNExcModelRegisterAll(void)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(!DYNExcModelsRegisterAllCalled) DYNExcModelsRegisterAllCalled = PETSC_TRUE;

  ierr = DYNExcModelRegister("IEEET1",DYNExcModelCreate_IEEET1,sizeof(struct _p_DYNIEEET1));CHKERRQ(ierr);
  ierr = DYNExcModelRegister("EXST1",DYNExcModelCreate_EXST1,sizeof(struct _p_DYNEXST1));CHKERRQ(ierr);
  ierr = DYNExcModelRegister("SEXS",DYNExcModelCreate_SEXS,sizeof(struct _p_DYNSEXS));CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNExcModelGetFirstVariableLocation - Gets the location for the first variable for the
                        dynexc model from the bus variables array

  Input Parameters:
. dynexc - the DYNExcModel object

  Output Parameters
. loc - the location for the first variable 
*/
PetscErrorCode DYNExcModelGetFirstVariableLocation(DYNExcModel dynexc,PetscInt *loc)
{
  PetscFunctionBegin;
  *loc = dynexc->startloc;
  PetscFunctionReturn(0);
}


/*
  DYNExcModelReadModelType - Parses a line in dyr file and reads the model type (string identifier). If the line does not have
                            the registered exciter model then the type is set to NONE.

  Input Parameters:
. line - the line to parse

  Output Parameters:
. modeltype - type of the model (string name of the model)
*/
PetscErrorCode DYNExcModelReadModelType(char *line, char gentype[])
{
  PetscInt i;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for(i=0; i < nexcmodelsregistered; i++) {
    if(strstr(line,DYNExcModelList[i].name) != NULL) {
      ierr = PetscStrcpy(gentype,DYNExcModelList[i].name);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
  }

  ierr = PetscStrcpy(gentype,DYNEXCNONE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNExcModelGetBusnumID - Gets the bus number and exciter ID for this dynamic exciter model
                            from the data file. Calls the underlying implementation
  Input Parameters:
. dynexc - the dynamic exciter model object

  Output Parameters:
+ busnum - bus number
. genid  - the exciter ID
*/
PetscErrorCode DYNExcModelGetBusnumID(DYNExcModel dynexc, PetscInt *busnum,char **genid)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = (*dynexc->ops.getbusnumid)(dynexc,busnum,genid);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNExcModelGetSizeofType - Gets the size of the implementation type (obtained using sizeof())

  Input Parameters:
. dynexc - the DYNExcModel object

  Output Parameters:
. size - the size of the derived class object
*/
PetscErrorCode DYNExcModelGetSizeof(DYNExcModel dynexc,PetscInt *size)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynexc->ops.getsizeof)(dynexc,size);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNExcModelGetNvar - Returns the number of variables for this dynamic exciter model

  Input Parameters:
. dynexc - the DYNExcModel object

  Output Parameters:
. Nvar - number of variables (dofs) for this exciter model
*/
PetscErrorCode DYNExcModelGetNvar(DYNExcModel dynexc,PetscInt *Nvar)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynexc->ops.getnvar)(dynexc,Nvar);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNExcModelSetType - Sets the type of the dynamic exciter model

  Input Parameters:
+ dynexc - the dynamic exciter model object
- modeltype - type of the model
*/
PetscErrorCode DYNExcModelSetType(DYNExcModel dynexc,DYNExcModelType modeltype)
{
  PetscErrorCode ierr,(*r)(DYNExcModel)=NULL;
  PetscInt       i;
  PetscFunctionBegin;
  PetscBool match;
  for(i=0;i < nexcmodelsregistered;i++) {
    ierr = PetscStrcmp(DYNExcModelList[i].name,modeltype,&match);CHKERRQ(ierr);
    if(match) {
      r = DYNExcModelList[i].create;
      break;
    }
  }

  dynexc->dyngen = 0;
  dynexc->dynstab = 0;

  if(!r) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown type for DYNExcModel %s",modeltype);
  /* Null the function pointers */
  dynexc->ops.readdata = 0;
  dynexc->ops.getnvar = 0;
  dynexc->ops.destroy = 0;
  dynexc->ops.getsizeof = 0;
  dynexc->ops.setinitialconditions = 0;
  dynexc->ops.daerhsfunction=0;
  dynexc->ops.daerhsjacobian=0;
  dynexc->ops.getequationtypes=0;
  dynexc->ops.getbusnumid=0;
  dynexc->ops.getfieldvoltage=0;
  dynexc->ops.setevent=0;

  /* Copy the type name */
  ierr = PetscStrcpy(dynexc->type,modeltype);CHKERRQ(ierr);

  /* Call the underlying implementation constructor */
  ierr = (*r)(dynexc);CHKERRQ(ierr);
  /* Set the numnber of variables for this dynamic exciter model */
  ierr = DYNExcModelGetNvar(dynexc,&dynexc->nvar);CHKERRQ(ierr);

  dynexc->startloc    = -1;
  dynexc->ndiff       =  0;
  dynexc->nalg        =  0;
  dynexc->eqtypes     =  0;
  dynexc->numMonitors = 0;
  dynexc->eventfcn    = NULL;
  dynexc->posteventfcn = NULL;
  dynexc->posteventdgamma = NULL;
  PetscFunctionReturn(0);
}

/*
  DYNExcModelCopy - Copies the dynexc model

  Input Parameters:
. dynexcin - the dynexc model to copy from

  Output Parameters:
. dynexcout - the dynexc model to copy to

  Notes:
   If a new method is created for the dynexcmodel then that should be copied too.
   The underlying implementation (dynexc->data) data pointer is copied.
*/
PetscErrorCode DYNExcModelCopy(DYNExcModel dynexcout,DYNExcModel dynexcin)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscMemcpy(dynexcout,dynexcin,sizeof(struct _p_DYNExcModel));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNExcModelGetModelData - Returns the model data associated with the particular implementation type

  Input Parameters:
. dynexc - the DYNExcModel object

  Output Parameters:
. modeldata - the model data struct
*/
PetscErrorCode DYNExcModelGetModelData(DYNExcModel dynexc,void **modeldata)
{
  PetscFunctionBegin;
  *modeldata = (void*)dynexc->data;
  PetscFunctionReturn(0);
}

/*
  DYNExcModelGetDynGen - Returns the generator associated with this exciter

  Input Parameters:
. dynexc - the DYNExcModel object

  Output Parameters:
. dyngen - the dyngen object
*/
PetscErrorCode DYNExcModelGetDynGen(DYNExcModel dynexc,DYNGenModel *dyngen)
{
  PetscFunctionBegin;
  *dyngen = dynexc->dyngen;
  PetscFunctionReturn(0);
}

/*
  DYNExcModelGetDynGen - Returns the stabilizer associated with this exciter

  Input Parameters:
. dynexc - the DYNExcModel object

  Output Parameters:
. dynstab - the dynstab object
*/
PetscErrorCode DYNExcModelGetDynStab(DYNExcModel dynexc,DYNStabModel *dynstab)
{
  PetscFunctionBegin;
  *dynstab = dynexc->dynstab;
  PetscFunctionReturn(0);
}

/*
  DYNExcModelReadData - Parses the given line and reads the data for the dynamic exciter
  model

  Input Parameters:
+ dynexc - The exciter base model
- line   - the file from the line

  Notes: The type of exciter model should be set prior to this call via DYNExcModelSetType()
*/
PetscErrorCode DYNExcModelReadData(DYNExcModel dynexc,char *line)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynexc->ops.readdata)(dynexc,line);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNExcModelSetInitialConditions - Sets the initial conditions for this exciter model

  Input Parameters:
+ dynexc - the dynamic exciter model object
. VD     - real component of the bus voltage
- VQ     - imaginary component of the bus voltage

  Output Parameters:
. xdyn   - the initial conditions for this exciter model

  Notes: 
   The initial conditions are populated in the xdyn array

   Other constants, parameters can be also set during this function call in
   its implementation struct (for e.g. setting reference voltage Vref and others)
*/
PetscErrorCode DYNExcModelSetInitialConditions(DYNExcModel dynexc,PetscScalar VD,PetscScalar VQ,PetscScalar *xdyn)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynexc->ops.setinitialconditions)(dynexc,VD,VQ,xdyn);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DYNExcModelSetInitialConditionsP(DYNExcModel dynexc,PetscScalar PG,PetscScalar QG,PetscScalar VA,PetscScalar VM,PetscInt PGloc,PetscInt QGloc,PetscInt VAloc,PetscInt VMloc,Mat ICp,PetscInt dynlocglob)
{
  //  PetscErrorCode ierr;

  PetscFunctionBegin;
  //ierr = (*dynexc->ops.setinitialconditionsp)(dynexc,PG,QG,VA,VM,PGloc,QGloc,VAloc,VMloc,ICp,dynlocglob);CHKERRQ(ierr);
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
PetscErrorCode DYNExcModelGetFieldVoltage(DYNExcModel dynexc,PetscReal t, PetscScalar *xdyn,PetscScalar *Efd,PetscInt *dEfddXexc_num,PetscScalar dEfddXexc[],PetscInt dEfddXexc_loc[])
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  ierr = (*dynexc->ops.getfieldvoltage)(dynexc,t,xdyn,Efd,dEfddXexc_num,dEfddXexc,dEfddXexc_loc);CHKERRQ(ierr);
  if(dEfddXexc_num) {
    /* Conversion to global locations */
    for(i=0; i < *dEfddXexc_num; i++) {
      dEfddXexc_loc[i] = dynexc->dyngen->bus->startlocglob + dEfddXexc_loc[i];
    }
  }
  PetscFunctionReturn(0);
}

/*
  DYNExcModelDAERHSJacobian - Computes the RHS Jacobian of the DAE equations

  Input Parameters:
+ dynexc     - the dynamic exciter model
. t          - the current time
. VD         - real-part of complex bus voltage
. VQ         - imaginary-part of complex bus voltage
. xdynexc    - exciter variables
. dynlocglob - starting location of exciter variables in the global vector 
- V_loc      - global location of VD and VQ V_loc[0] = VD_loc V_loc[1] = VQ_loc


  Output Parameters:
. J          - the Jacobian matrix

  Notes:
   For each bus, the network equations are ordered as [IQ_loc;ID_loc]
*/ 
PetscErrorCode DYNExcModelDAERHSJacobian(DYNExcModel dynexc,Mat J,PetscReal t,PetscScalar VD,PetscScalar VQ,PetscScalar *xdynexc,PetscInt dynlocglob,PetscInt V_loc[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynexc->ops.daerhsjacobian)(dynexc,J,t,VD,VQ,xdynexc,dynlocglob,V_loc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNExcModelDAERHSJacobianP - Computes the RHS Jacobian of the DAE equations w.r.t. parameters

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
PetscErrorCode DYNExcModelDAERHSJacobianP(DYNExcModel dynexc,PetscReal t,const PetscScalar *xdynexc,Mat jacP,PetscInt dynlocglob,PetscInt Valoc,PetscInt Vmloc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynexc->ops.daerhsjacobianp)(dynexc,t,xdynexc,jacP,dynlocglob,Valoc,Vmloc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNExcModelDAEFWDRHSJacobianP - Computes the RHS Jacobian of the DAE equations w.r.t. parameters

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
PetscErrorCode DYNExcModelDAEFWDRHSJacobianP(DYNExcModel dynexc,PetscReal t,const PetscScalar *xdynexc,Vec *jacP,PetscInt dynlocglob,PetscInt Valoc,PetscInt Vmloc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynexc->ops.daefwdrhsjacobianp)(dynexc,t,xdynexc,jacP,dynlocglob,Valoc,Vmloc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNExcModelDAERHSFunction - Computes the RHS of the DAE function for the exciter model [f(x,y);g(x,y)] 
                              and also returns the exciter currents in network refernce frame

  Input Parameters:
+ dynexc - the dynamic exciter model
. t      - the current time
. x      - array of variables for the bus on which the generator is incident
. VD     - real-part of bus voltage
. VQ     - imaginary part of bus voltage

  Output Parameters:
+ f      - the residual array in which the exciter equation residuals are inserted. The location
           of the exciter variables can be obtained using DYNExcModelGetFirstVariableLocation(). Using
           the first location, the residuals can be inserted as
           f[loc] = 1st exciter equation residual
           f[loc+1] = 2nd ..
         ... and so on
*/
PetscErrorCode DYNExcModelDAERHSFunction(DYNExcModel dynexc,PetscReal t,PetscScalar VD, PetscScalar VQ, PetscScalar *x,PetscScalar *f)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynexc->ops.daerhsfunction)(dynexc,t,VD,VQ,x,f);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
			 
/*
  DYNExcModelGetEquationTypes - Gets the number and indices of differential and algebraic equations for the exciter model
                                                   
  Input Parameters:
. dynexc - the DYNExcModel object

  Output Parameters:
+ ndiff - number of differential equations
. nalg  - number of algebraic equations
- eqtype - an array of size nvar that has
            eqtype[i] = DIFF_EQ if equation i is a differential
            eqtype[i] = ALG_EQ if equation  i is an algebraic

  NOTES: The user does not need to create the array vartype
*/
PetscErrorCode DYNExcModelGetEquationTypes(DYNExcModel dynexc,PetscInt *ndiff,PetscInt *nalg, PetscInt *eqtype)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynexc->ops.getequationtypes)(dynexc,ndiff,nalg,eqtype);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
