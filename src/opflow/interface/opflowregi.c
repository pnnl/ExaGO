
#include <private/opflowimpl.h>
#include <exago_config.h>

/*
  OPFLOWModelRegister - Registers an OPFLOW model

  Input Parameters:
+ sname     - model name (string)
- createfunction  - the class constructor
*/
PetscErrorCode OPFLOWModelRegister(OPFLOW opflow,const char sname[],PetscErrorCode (*createfunction)(OPFLOW))
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscFunctionBegin;
  for(i=0; i < opflow->nmodelsregistered;i++) {
    PetscBool match;
    ierr = PetscStrcmp(opflow->OPFLOWModelList[i].name,sname,&match);
    if(match) PetscFunctionReturn(0);
  }
  if(opflow->nmodelsregistered == OPFLOWMODELSMAX) {
    SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,"Cannot register %s OPFLOW model, maximum limit %d reached\n",sname,OPFLOWMODELSMAX);
  }
  i = opflow->nmodelsregistered;
  ierr = PetscStrcpy(opflow->OPFLOWModelList[i].name,sname);CHKERRQ(ierr);
  opflow->OPFLOWModelList[i].create = createfunction;
  opflow->nmodelsregistered++;
  PetscFunctionReturn(0);
}

/*
  OPFLOWSolverRegister - Registers an OPFLOW model

  Input Parameters:
+ sname     - solver name (string)
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
  if(opflow->nsolversregistered == OPFLOWSOLVERSMAX) {
    SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,"Cannot register %s OPFLOW solver, maximum limit %d reached\n",sname,OPFLOWSOLVERSMAX);
  }
  i = opflow->nsolversregistered;
  ierr = PetscStrcpy(opflow->OPFLOWSolverList[i].name,sname);CHKERRQ(ierr);
  opflow->OPFLOWSolverList[i].create = createfunction;
  opflow->nsolversregistered++;
  PetscFunctionReturn(0);
}

extern PetscErrorCode OPFLOWModelCreate_PBPOL(OPFLOW);
extern PetscErrorCode OPFLOWModelCreate_PBCAR(OPFLOW);
extern PetscErrorCode OPFLOWModelCreate_IBCAR(OPFLOW);
extern PetscErrorCode OPFLOWModelCreate_IBCAR2(OPFLOW);
extern PetscErrorCode OPFLOWModelCreate_PBPOL2(OPFLOW);

#if defined(EXAGO_HAVE_HIOP)
extern PetscErrorCode OPFLOWModelCreate_PBPOLHIOP(OPFLOW);
#endif

#if defined(EXAGO_HAVE_RAJA)
extern PetscErrorCode OPFLOWModelCreate_PBPOLRAJAHIOP(OPFLOW);
extern PetscErrorCode OPFLOWModelCreate_PBPOLRAJA(OPFLOW);
#endif

/*
  OPFLOWModelRegisterAll - Registers all built OPFLOW models
*/
PetscErrorCode OPFLOWModelRegisterAll(OPFLOW opflow)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(opflow->OPFLOWModelRegisterAllCalled) PetscFunctionReturn(0);

  ierr = OPFLOWModelRegister(opflow,OPFLOWMODEL_PBPOL,OPFLOWModelCreate_PBPOL);CHKERRQ(ierr);
  ierr = OPFLOWModelRegister(opflow,OPFLOWMODEL_PBCAR,OPFLOWModelCreate_PBCAR);CHKERRQ(ierr);
  ierr = OPFLOWModelRegister(opflow,OPFLOWMODEL_IBCAR,OPFLOWModelCreate_IBCAR);CHKERRQ(ierr);
  ierr = OPFLOWModelRegister(opflow,OPFLOWMODEL_IBCAR2,OPFLOWModelCreate_IBCAR2);CHKERRQ(ierr);
  ierr = OPFLOWModelRegister(opflow,OPFLOWMODEL_PBPOL2,OPFLOWModelCreate_PBPOL2);CHKERRQ(ierr);

#if defined(EXAGO_HAVE_HIOP)
  ierr = OPFLOWModelRegister(opflow,OPFLOWMODEL_PBPOLHIOP,OPFLOWModelCreate_PBPOLHIOP);CHKERRQ(ierr);
#endif

#if defined(EXAGO_HAVE_RAJA)
  ierr = OPFLOWModelRegister(opflow,OPFLOWMODEL_PBPOLRAJAHIOP,OPFLOWModelCreate_PBPOLRAJAHIOP);CHKERRQ(ierr);
  ierr = OPFLOWModelRegister(opflow,OPFLOWMODEL_PBPOLRAJA,OPFLOWModelCreate_PBPOLRAJA);CHKERRQ(ierr);
#endif

  opflow->OPFLOWModelRegisterAllCalled = PETSC_TRUE;

  PetscFunctionReturn(0);
}

#if defined(EXAGO_HAVE_IPOPT)
extern PetscErrorCode OPFLOWSolverCreate_IPOPT(OPFLOW);
#endif
extern PetscErrorCode OPFLOWSolverCreate_TAO(OPFLOW);
#if defined(EXAGO_HAVE_HIOP)
extern PetscErrorCode OPFLOWSolverCreate_HIOP(OPFLOW);
extern PetscErrorCode OPFLOWSolverCreate_HIOPNEW(OPFLOW);
#endif

/*
  OPFLOWSolverRegisterAll - Registers all OPFLOW solvers
*/
PetscErrorCode OPFLOWSolverRegisterAll(OPFLOW opflow)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(opflow->OPFLOWSolverRegisterAllCalled) PetscFunctionReturn(0);

#if defined(EXAGO_HAVE_IPOPT)
  ierr = OPFLOWSolverRegister(opflow,OPFLOWSOLVER_IPOPT,OPFLOWSolverCreate_IPOPT);CHKERRQ(ierr);
#endif
  ierr = OPFLOWSolverRegister(opflow,OPFLOWSOLVER_TAO,OPFLOWSolverCreate_TAO);CHKERRQ(ierr);
#if defined(EXAGO_HAVE_HIOP)
  ierr = OPFLOWSolverRegister(opflow,OPFLOWSOLVER_HIOP,OPFLOWSolverCreate_HIOP);CHKERRQ(ierr);
  ierr = OPFLOWSolverRegister(opflow,OPFLOWSOLVER_HIOPNEW,OPFLOWSolverCreate_HIOPNEW);CHKERRQ(ierr);
#endif
  opflow->OPFLOWSolverRegisterAllCalled = PETSC_TRUE;

  PetscFunctionReturn(0);
}
