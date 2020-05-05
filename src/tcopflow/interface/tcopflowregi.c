#include <scopflow_config.h>
#include <private/tcopflowimpl.h>

/*
  TCOPFLOWSolverRegister - Registers an TCOPFLOW formulation

  Input Parameters:
+ sname     - formuation name (string)
- createfunction  - the class constructor
*/
PetscErrorCode TCOPFLOWSolverRegister(TCOPFLOW tcopflow,const char sname[],PetscErrorCode (*createfunction)(TCOPFLOW))
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscFunctionBegin;
  for(i=0; i < tcopflow->nsolversregistered;i++) {
    PetscBool match;
    ierr = PetscStrcmp(tcopflow->TCOPFLOWSolverList[i].name,sname,&match);
    if(match) PetscFunctionReturn(0);
  }
  i = tcopflow->nsolversregistered;
  ierr = PetscStrcpy(tcopflow->TCOPFLOWSolverList[i].name,sname);CHKERRQ(ierr);
  tcopflow->TCOPFLOWSolverList[i].create = createfunction;
  tcopflow->nsolversregistered++;
  PetscFunctionReturn(0);
}

#if defined(SCOPFLOW_HAVE_IPOPT)
extern PetscErrorCode TCOPFLOWSolverCreate_IPOPT(TCOPFLOW);
#endif

/*
  TCOPFLOWSolverRegisterAll - Registers all TCOPFLOW solvers
*/
PetscErrorCode TCOPFLOWSolverRegisterAll(TCOPFLOW tcopflow)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(tcopflow->TCOPFLOWSolverRegisterAllCalled) PetscFunctionReturn(0);

#if defined(SCOPFLOW_HAVE_IPOPT)
  ierr = TCOPFLOWSolverRegister(tcopflow,TCOPFLOWSOLVER_IPOPT,TCOPFLOWSolverCreate_IPOPT);CHKERRQ(ierr);
#endif

  tcopflow->TCOPFLOWSolverRegisterAllCalled = PETSC_TRUE;

  PetscFunctionReturn(0);
}
