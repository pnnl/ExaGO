#include <private/opflowimpl.h>

/*
  OPFLOWFormulationRegister - Registers an OPFLOW formulation

  Input Parameters:
+ sname     - formuation name (string)
- createfunction  - the class constructor
*/
PetscErrorCode OPFLOWFormulationRegister(OPFLOW opflow,const char sname[],PetscErrorCode (*createfunction)(OPFLOW))
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscFunctionBegin;
  for(i=0; i < opflow->nformulationsregistered;i++) {
    PetscBool match;
    ierr = PetscStrcmp(opflow->OPFLOWFormulationList[i].name,sname,&match);
    if(match) PetscFunctionReturn(0);
  }
  i = opflow->nformulationsregistered;
  ierr = PetscStrcpy(opflow->OPFLOWFormulationList[i].name,sname);CHKERRQ(ierr);
  opflow->OPFLOWFormulationList[i].create = createfunction;
  opflow->nformulationsregistered++;
  PetscFunctionReturn(0);
}

/*
  OPFLOWSolverRegister - Registers an OPFLOW formulation

  Input Parameters:
+ sname     - formuation name (string)
- createfunction  - the class constructor
*/
PetscErrorCode OPFLOWSolverRegister(OPFLOW opflow,const char sname[],PetscErrorCode (*createfunction)(OPFLOW))
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscFunctionBegin;
  for(i=0; i < opflow->nsolversregistered;i++) {
    PetscBool match;
    ierr = PetscStrcmp(opflow->OPFLOWSolverList[i].name,sname,&match);
    if(match) PetscFunctionReturn(0);
  }
  i = opflow->nsolversregistered;
  ierr = PetscStrcpy(opflow->OPFLOWSolverList[i].name,sname);CHKERRQ(ierr);
  opflow->OPFLOWSolverList[i].create = createfunction;
  opflow->nsolversregistered++;
  PetscFunctionReturn(0);
}

extern PetscErrorCode OPFLOWFormulationCreate_PBPOL(OPFLOW);

/*
  OPFLOWFormulationRegisterAll - Registers all built OPFLOW formulations
*/
PetscErrorCode OPFLOWFormulationRegisterAll(OPFLOW opflow)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(opflow->OPFLOWFormulationRegisterAllCalled) PetscFunctionReturn(0);

  ierr = OPFLOWFormulationRegister(opflow,OPFLOWFORMULATION_PBPOL,OPFLOWFormulationCreate_PBPOL);CHKERRQ(ierr);
  opflow->OPFLOWFormulationRegisterAllCalled = PETSC_TRUE;

  PetscFunctionReturn(0);
}

extern PetscErrorCode OPFLOWSolverCreate_IPOPT(OPFLOW);
/*
  OPFLOWSolverRegisterAll - Registers all OPFLOW solvers
*/
PetscErrorCode OPFLOWSolverRegisterAll(OPFLOW opflow)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(opflow->OPFLOWSolverRegisterAllCalled) PetscFunctionReturn(0);

  ierr = OPFLOWSolverRegister(opflow,OPFLOWSOLVER_IPOPT,OPFLOWSolverCreate_IPOPT);CHKERRQ(ierr);
  opflow->OPFLOWSolverRegisterAllCalled = PETSC_TRUE;

  PetscFunctionReturn(0);
}
