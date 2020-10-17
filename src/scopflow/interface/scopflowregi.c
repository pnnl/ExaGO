#include <exago_config.h>
#include <private/scopflowimpl.h>

/*
  SCOPFLOWSolverRegister - Registers an SCOPFLOW model

  Input Parameters:
+ sname           - model name (string)
- createfunction  - the class constructor
*/
PetscErrorCode SCOPFLOWSolverRegister(SCOPFLOW scopflow,const char sname[],PetscErrorCode (*createfunction)(SCOPFLOW))
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscFunctionBegin;
  for(i=0; i < scopflow->nsolversregistered;i++) {
    PetscBool match;
    ierr = PetscStrcmp(scopflow->SCOPFLOWSolverList[i].name,sname,&match);
    if(match) PetscFunctionReturn(0);
  }
  i = scopflow->nsolversregistered;
  ierr = PetscStrcpy(scopflow->SCOPFLOWSolverList[i].name,sname);CHKERRQ(ierr);
  scopflow->SCOPFLOWSolverList[i].create = createfunction;
  scopflow->nsolversregistered++;
  PetscFunctionReturn(0);
}

#if defined(EXAGO_HAVE_IPOPT)
extern PetscErrorCode SCOPFLOWSolverCreate_IPOPT(SCOPFLOW);
#endif

#if defined(EXAGO_HAVE_PIPS)
extern PetscErrorCode SCOPFLOWSolverCreate_PIPS(SCOPFLOW);
#endif

extern PetscErrorCode SCOPFLOWSolverCreate_EMPAR(SCOPFLOW);

/*
  SCOPFLOWSolverRegisterAll - Registers all SCOPFLOW solvers
*/
PetscErrorCode SCOPFLOWSolverRegisterAll(SCOPFLOW scopflow)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(scopflow->SCOPFLOWSolverRegisterAllCalled) PetscFunctionReturn(0);

#if defined(EXAGO_HAVE_IPOPT)
  ierr = SCOPFLOWSolverRegister(scopflow,SCOPFLOWSOLVER_IPOPT,SCOPFLOWSolverCreate_IPOPT);CHKERRQ(ierr);
#endif

#if defined(EXAGO_HAVE_PIPS)
  ierr = SCOPFLOWSolverRegister(scopflow,SCOPFLOWSOLVER_PIPS,SCOPFLOWSolverCreate_PIPS);CHKERRQ(ierr);
#endif
  ierr = SCOPFLOWSolverRegister(scopflow,SCOPFLOWSOLVER_EMPAR,SCOPFLOWSolverCreate_EMPAR);CHKERRQ(ierr);
  scopflow->SCOPFLOWSolverRegisterAllCalled = PETSC_TRUE;

  PetscFunctionReturn(0);
}
