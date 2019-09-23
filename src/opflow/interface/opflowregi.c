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
  for(i=0;i < opflow->nsolversregistered;i++) {
    ierr = PetscStrcmp(opflow->OPFLOWSolverList[i].name,solvername,&match);CHKERRQ(ierr);
    if(match) {
      r = opflow->OPFLOWSolverList[i].create;
      break;
    }
  }

  if(!r) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown type for OPFLOW Solver %s",solvername);

  /* Initialize (Null) the function pointers */
  opflow->solverops.destroy = 0;
  opflow->solverops.solve   = 0;
  opflow->solverops.setup   = 0;

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
  for(i=0;i < opflow->nformulationsregistered;i++) {
    ierr = PetscStrcmp(opflow->OPFLOWFormulationList[i].name,formulationname,&match);CHKERRQ(ierr);
    if(match) {
      r = opflow->OPFLOWFormulationList[i].create;
      break;
    }
  }

  if(!r) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown type for OPFLOW Formulation %s",formulationname);

  /* Null the function pointers */
  opflow->formops.destroy                        = 0;
  opflow->formops.setvariablebounds              = 0;
  opflow->formops.setconstraintbounds            = 0;
  opflow->formops.setvariableandconstraintbounds = 0;
  opflow->formops.setinitialguess                = 0;
  opflow->formops.computeequalityconstraints     = 0;
  opflow->formops.computeinequalityconstraints   = 0;
  opflow->formops.computeconstraints             = 0;
  opflow->formops.computeobjandgradient          = 0;
  opflow->formops.computeobjective               = 0;
  opflow->formops.computegradient                = 0;
  opflow->formops.computejacobian                = 0;
  opflow->formops.computeequalityconstraintjacobian = 0;
  opflow->formops.computeinequalityconstraintjacobian = 0;
  opflow->formops.computehessian                 = 0;

  /* Call the underlying implementation constructor */
  ierr = (*r)(opflow);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
