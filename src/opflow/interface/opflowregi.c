#include <private/opflowimpl.h>

struct _p_OPFLOWFormulationList OPFLOWFormulationList[8];

static PetscInt nformulationsregistered=0;

struct _p_OPFLOWSolverList OPFLOWSolverList[8];

static PetscInt nsolversregistered=0;

static PetscBool OPFLOWFormulationRegisterAllCalled = PETSC_FALSE;

static PetscBool OPFLOWSolverRegisterAllCalled = PETSC_FALSE;

/*
  OPFLOWFormulationRegister - Registers an OPFLOW formulation

  Input Parameters:
+ sname     - formuation name (string)
- createfunction  - the class constructor
*/
PetscErrorCode OPFLOWFormulationRegister(const char sname[],PetscErrorCode (*createfunction)(OPFLOW))
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscFunctionBegin;
  for(i=0; i < nformulationsregistered;i++) {
    PetscBool match;
    ierr = PetscStrcmp(OPFLOWFormulationList[i].name,sname,&match);
    if(match) PetscFunctionReturn(0);
  }
  i = nformulationsregistered;
  ierr = PetscStrcpy(OPFLOWFormulationList[i].name,sname);CHKERRQ(ierr);
  OPFLOWFormulationList[i].create = createfunction;
  nformulationsregistered++;
  PetscFunctionReturn(0);
}

/*
  OpflowsolverRegister - Registers an OPFLOW formulation

  Input Parameters:
+ sname     - formuation name (string)
- createfunction  - the class constructor
*/
PetscErrorCode OPFLOWSolverRegister(const char sname[],PetscErrorCode (*createfunction)(OPFLOW))
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscFunctionBegin;
  for(i=0; i < nsolversregistered;i++) {
    PetscBool match;
    ierr = PetscStrcmp(OPFLOWSolverList[i].name,sname,&match);
    if(match) PetscFunctionReturn(0);
  }
  i = nsolversregistered;
  ierr = PetscStrcpy(OPFLOWSolverList[i].name,sname);CHKERRQ(ierr);
  OPFLOWSolverList[i].create = createfunction;
  nsolversregistered++;
  PetscFunctionReturn(0);
}

extern PetscErrorCode OPFLOWFormulationCreate_PBPOL(OPFLOW);

/*
  OPFLOWFormulationRegisterAll - Registers all built OPFLOW formulations
*/
PetscErrorCode OPFLOWFormulationRegisterAll(void)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(OPFLOWFormulationRegisterAllCalled) PetscFunctionReturn(0);

  ierr = OPFLOWFormulationRegister("POWER_BALANCE_POLAR",OPFLOWFormulationCreate_PBPOL);CHKERRQ(ierr);
  OPFLOWFormulationRegisterAllCalled = PETSC_TRUE;

  PetscFunctionReturn(0);
}

/*
  OPFLOWSolverRegisterAll - Registers all OPFLOW solvers
*/
PetscErrorCode OPFLOWSolverRegisterAll(void)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(OPFLOWSolverRegisterAllCalled) PetscFunctionReturn(0);

  //  ierr = OPFLOWSolverRegister("ipopt",OPFLOWSolverCreate_IPOPT);CHKERRQ(ierr);
  OPFLOWSolverRegisterAllCalled = PETSC_TRUE;

  PetscFunctionReturn(0);
}

/*
  OPFLOWSetSolver - Sets the solver for OPFLOW

  Input Parameters:
+ opflow - opflow application object
- solvername - name of the solver
*/
PetscErrorCode OPFLOWSetSolver(OPFLOW opflow,const char* solvername)
{
  PetscErrorCode ierr,(*r)(OPFLOW)=NULL;
  PetscInt       i;
  PetscFunctionBegin;
  PetscBool match;
  for(i=0;i < nsolversregistered;i++) {
    ierr = PetscStrcmp(OPFLOWSolverList[i].name,solvername,&match);CHKERRQ(ierr);
    if(match) {
      r = OPFLOWSolverList[i].create;
      break;
    }
  }

  if(!r) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown type for OPFLOW Solver %s",solvername);

  /* Null the function pointers */
  opflow->solverops.destroy = 0;

  opflow->formops.destroy   = 0;

  /* Call the underlying implementation constructor */
  ierr = (*r)(opflow);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  OPFLOWSetFormulation - Sets the formulation for OPFLOW

  Input Parameters:
+ opflow - opflow application object
- solvername - name of the formulation
*/
PetscErrorCode OPFLOWSetFormulation(OPFLOW opflow,const char* formulationname)
{
  PetscErrorCode ierr,(*r)(OPFLOW)=NULL;
  PetscInt       i;
  PetscFunctionBegin;
  PetscBool match;
  for(i=0;i < nformulationsregistered;i++) {
    ierr = PetscStrcmp(OPFLOWFormulationList[i].name,formulationname,&match);CHKERRQ(ierr);
    if(match) {
      r = OPFLOWFormulationList[i].create;
      break;
    }
  }

  if(!r) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown type for OPFLOW Formulation %s",formulationname);

  /* Null the function pointers */

  /* Call the underlying implementation constructor */
  ierr = (*r)(opflow);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
