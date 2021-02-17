#include <exago_config.h>
#include <private/scopflowimpl.h>

/*
  SCOPFLOWModelRegister - Registers a SCOPFLOW model

  Input Parameters:
+ sname     - model name (string)
- createfunction  - the class constructor
*/
PetscErrorCode SCOPFLOWModelRegister(SCOPFLOW scopflow,const char sname[],PetscErrorCode (*createfunction)(SCOPFLOW))
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscFunctionBegin;
  for(i=0; i < scopflow->nmodelsregistered;i++) {
    PetscBool match;
    ierr = PetscStrcmp(scopflow->SCOPFLOWModelList[i].name,sname,&match);
    if(match) PetscFunctionReturn(0);
  }
  if(scopflow->nmodelsregistered == SCOPFLOWMODELSMAX) {
    SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,"Cannot register %s OPFLOW model, maximum limit %d reached\n",sname,SCOPFLOWMODELSMAX);
  }
  i = scopflow->nmodelsregistered;
  ierr = PetscStrcpy(scopflow->SCOPFLOWModelList[i].name,sname);CHKERRQ(ierr);
  scopflow->SCOPFLOWModelList[i].create = createfunction;
  scopflow->nmodelsregistered++;
  PetscFunctionReturn(0);
}

extern PetscErrorCode SCOPFLOWModelCreate_GENRAMP(SCOPFLOW);
extern PetscErrorCode SCOPFLOWModelCreate_GENRAMPT(SCOPFLOW);

/*
  SCOPFLOWModelRegisterAll - Registers all built-in SCOPFLOW models
*/
PetscErrorCode SCOPFLOWModelRegisterAll(SCOPFLOW scopflow)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(scopflow->SCOPFLOWModelRegisterAllCalled) PetscFunctionReturn(0);

  ierr = SCOPFLOWModelRegister(scopflow,SCOPFLOWMODEL_GENRAMP,SCOPFLOWModelCreate_GENRAMP);CHKERRQ(ierr);
  ierr = SCOPFLOWModelRegister(scopflow,SCOPFLOWMODEL_GENRAMPT,SCOPFLOWModelCreate_GENRAMPT);CHKERRQ(ierr);

  scopflow->SCOPFLOWModelRegisterAllCalled = PETSC_TRUE;

  PetscFunctionReturn(0);
}

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

#if defined(EXAGO_ENABLE_IPOPT)
extern PetscErrorCode SCOPFLOWSolverCreate_IPOPT(SCOPFLOW);
extern PetscErrorCode SCOPFLOWSolverCreate_IPOPTNEW(SCOPFLOW);
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

#if defined(EXAGO_ENABLE_IPOPT)
  ierr = SCOPFLOWSolverRegister(scopflow,SCOPFLOWSOLVER_IPOPT,SCOPFLOWSolverCreate_IPOPT);CHKERRQ(ierr);
  ierr = SCOPFLOWSolverRegister(scopflow,SCOPFLOWSOLVER_IPOPTNEW,SCOPFLOWSolverCreate_IPOPTNEW);CHKERRQ(ierr);

#endif

  ierr = SCOPFLOWSolverRegister(scopflow,SCOPFLOWSOLVER_EMPAR,SCOPFLOWSolverCreate_EMPAR);CHKERRQ(ierr);
  scopflow->SCOPFLOWSolverRegisterAllCalled = PETSC_TRUE;

  PetscFunctionReturn(0);
}
