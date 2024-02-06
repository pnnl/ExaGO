
#include <exago_config.h>
#include <private/opflowimpl.h>

/*
  OPFLOWModelRegister - Registers an OPFLOW model

  Input Parameters:
+ sname     - model name (string)
- createfunction  - the class constructor
*/
PetscErrorCode OPFLOWModelRegister(OPFLOW opflow, const char sname[],
                                   PetscErrorCode (*createfunction)(OPFLOW)) {
  PetscErrorCode ierr;
  PetscInt i;
  PetscFunctionBegin;
  for (i = 0; i < opflow->nmodelsregistered; i++) {
    PetscBool match;
    ierr = PetscStrcmp(opflow->OPFLOWModelList[i].name, sname, &match);
    if (match)
      PetscFunctionReturn(0);
  }
  if (opflow->nmodelsregistered == OPFLOWMODELSMAX) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
            "Cannot register %s OPFLOW model, maximum limit %d reached\n",
            sname, OPFLOWMODELSMAX);
  }
  i = opflow->nmodelsregistered;
  ierr = PetscStrcpy(opflow->OPFLOWModelList[i].name, sname);
  CHKERRQ(ierr);
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
PetscErrorCode OPFLOWSolverRegister(OPFLOW opflow, const char sname[],
                                    PetscErrorCode (*createfunction)(OPFLOW)) {
  PetscErrorCode ierr;
  PetscInt i;
  PetscFunctionBegin;
  for (i = 0; i < opflow->nsolversregistered; i++) {
    PetscBool match;
    ierr = PetscStrcmp(opflow->OPFLOWSolverList[i].name, sname, &match);
    if (match)
      PetscFunctionReturn(0);
  }
  if (opflow->nsolversregistered == OPFLOWSOLVERSMAX) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
            "Cannot register %s OPFLOW solver, maximum limit %d reached\n",
            sname, OPFLOWSOLVERSMAX);
  }
  i = opflow->nsolversregistered;
  ierr = PetscStrcpy(opflow->OPFLOWSolverList[i].name, sname);
  CHKERRQ(ierr);
  opflow->OPFLOWSolverList[i].create = createfunction;
  opflow->nsolversregistered++;
  PetscFunctionReturn(0);
}

extern PetscErrorCode OPFLOWModelCreate_PBPOL(OPFLOW);
extern PetscErrorCode OPFLOWModelCreate_PBCAR(OPFLOW);
extern PetscErrorCode OPFLOWModelCreate_IBCAR(OPFLOW);
extern PetscErrorCode OPFLOWModelCreate_IBCAR2(OPFLOW);
extern PetscErrorCode OPFLOWModelCreate_DCOPF(OPFLOW);
extern PetscErrorCode OPFLOWModelCreate_PBPOLUC(OPFLOW);

#if defined(EXAGO_ENABLE_HIOP)
extern PetscErrorCode OPFLOWModelCreate_PBPOLHIOP(OPFLOW);
#endif

#if defined(EXAGO_ENABLE_RAJA)
extern PetscErrorCode OPFLOWModelCreate_PBPOLRAJAHIOP(OPFLOW);
#if defined(EXAGO_ENABLE_HIOP_SPARSE)
extern PetscErrorCode OPFLOWModelCreate_PBPOLRAJAHIOPSPARSE(OPFLOW);
#endif
#endif

/*
  OPFLOWModelRegisterAll - Registers all built OPFLOW models
*/
PetscErrorCode OPFLOWModelRegisterAll(OPFLOW opflow) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if (opflow->OPFLOWModelRegisterAllCalled)
    PetscFunctionReturn(0);

#if defined(EXAGO_ENABLE_IPOPT) || defined(EXAGO_ENABLE_HIOP_SPARSE)
  ierr =
      OPFLOWModelRegister(opflow, OPFLOWMODEL_PBPOL, OPFLOWModelCreate_PBPOL);
  CHKERRQ(ierr);
#endif
  ierr =
      OPFLOWModelRegister(opflow, OPFLOWMODEL_PBCAR, OPFLOWModelCreate_PBCAR);
  CHKERRQ(ierr);
  ierr =
      OPFLOWModelRegister(opflow, OPFLOWMODEL_IBCAR, OPFLOWModelCreate_IBCAR);
  CHKERRQ(ierr);
  ierr =
      OPFLOWModelRegister(opflow, OPFLOWMODEL_IBCAR2, OPFLOWModelCreate_IBCAR2);
  CHKERRQ(ierr);
  ierr = OPFLOWModelRegister(opflow, OPFLOWMODEL_DCOPF, OPFLOWModelCreate_DCOPF);
  CHKERRQ(ierr);

  ierr = OPFLOWModelRegister(opflow, OPFLOWMODEL_PBPOLUC, OPFLOWModelCreate_PBPOLUC);
  CHKERRQ(ierr);


#if defined(EXAGO_ENABLE_HIOP)
  ierr = OPFLOWModelRegister(opflow, OPFLOWMODEL_PBPOLHIOP,
                             OPFLOWModelCreate_PBPOLHIOP);
  CHKERRQ(ierr);
#endif

#if defined(EXAGO_ENABLE_RAJA)
  ierr = OPFLOWModelRegister(opflow, OPFLOWMODEL_PBPOLRAJAHIOP,
                             OPFLOWModelCreate_PBPOLRAJAHIOP);
  CHKERRQ(ierr);
#if defined(EXAGO_ENABLE_HIOP_SPARSE)
  ierr = OPFLOWModelRegister(opflow, OPFLOWMODEL_PBPOLRAJAHIOPSPARSE,
                             OPFLOWModelCreate_PBPOLRAJAHIOPSPARSE);
  CHKERRQ(ierr);
#endif
#endif

  opflow->OPFLOWModelRegisterAllCalled = PETSC_TRUE;

  PetscFunctionReturn(0);
}

#if defined(EXAGO_ENABLE_IPOPT)
extern PetscErrorCode OPFLOWSolverCreate_IPOPT(OPFLOW);
#endif
#if defined(EXAGO_ENABLE_HIOP)
extern PetscErrorCode OPFLOWSolverCreate_HIOP(OPFLOW);
#if defined(EXAGO_ENABLE_HIOP_SPARSE)
extern PetscErrorCode OPFLOWSolverCreate_HIOPSPARSE(OPFLOW);
#if defined(EXAGO_ENABLE_RAJA)
extern PetscErrorCode OPFLOWSolverCreate_HIOPSPARSEGPU(OPFLOW);
#endif
#endif
#endif

/*
  OPFLOWSolverRegisterAll - Registers all OPFLOW solvers
*/
PetscErrorCode OPFLOWSolverRegisterAll(OPFLOW opflow) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if (opflow->OPFLOWSolverRegisterAllCalled)
    PetscFunctionReturn(0);

#if defined(EXAGO_ENABLE_IPOPT)
  ierr = OPFLOWSolverRegister(opflow, OPFLOWSOLVER_IPOPT,
                              OPFLOWSolverCreate_IPOPT);
  CHKERRQ(ierr);
#endif
#if defined(EXAGO_ENABLE_HIOP)
  ierr =
      OPFLOWSolverRegister(opflow, OPFLOWSOLVER_HIOP, OPFLOWSolverCreate_HIOP);
  CHKERRQ(ierr);
#if defined(EXAGO_ENABLE_HIOP_SPARSE)
  ierr = OPFLOWSolverRegister(opflow, OPFLOWSOLVER_HIOPSPARSE,
                              OPFLOWSolverCreate_HIOPSPARSE);
  CHKERRQ(ierr);
#if defined(EXAGO_ENABLE_RAJA)
  ierr = OPFLOWSolverRegister(opflow, OPFLOWSOLVER_HIOPSPARSEGPU,
                              OPFLOWSolverCreate_HIOPSPARSEGPU);
  CHKERRQ(ierr);
#endif
#endif
#endif
  opflow->OPFLOWSolverRegisterAllCalled = PETSC_TRUE;

  PetscFunctionReturn(0);
}
