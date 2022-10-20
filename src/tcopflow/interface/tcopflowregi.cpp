#include <exago_config.h>
#include <private/tcopflowimpl.h>

/*
  TCOPFLOWModelRegister - Registers a TCOPFLOW model

  Input Parameters:
+ sname     - model name (string)
- createfunction  - the class constructor
*/
PetscErrorCode
TCOPFLOWModelRegister(TCOPFLOW tcopflow, const char sname[],
                      PetscErrorCode (*createfunction)(TCOPFLOW)) {
  PetscErrorCode ierr;
  PetscInt i;
  PetscFunctionBegin;
  for (i = 0; i < tcopflow->nmodelsregistered; i++) {
    PetscBool match;
    ierr = PetscStrcmp(tcopflow->TCOPFLOWModelList[i].name, sname, &match);
    if (match)
      PetscFunctionReturn(0);
  }
  if (tcopflow->nmodelsregistered == TCOPFLOWMODELSMAX) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
             "Cannot register %s OPFLOW model, maximum limit %d reached\n",
             sname, TCOPFLOWMODELSMAX);
  }
  i = tcopflow->nmodelsregistered;
  ierr = PetscStrcpy(tcopflow->TCOPFLOWModelList[i].name, sname);
  CHKERRQ(ierr);
  tcopflow->TCOPFLOWModelList[i].create = createfunction;
  tcopflow->nmodelsregistered++;
  PetscFunctionReturn(0);
}

extern PetscErrorCode TCOPFLOWModelCreate_GENRAMP(TCOPFLOW);

/*
  TCOPFLOWModelRegisterAll - Registers all built-in TCOPFLOW models
*/
PetscErrorCode TCOPFLOWModelRegisterAll(TCOPFLOW tcopflow) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if (tcopflow->TCOPFLOWModelRegisterAllCalled)
    PetscFunctionReturn(0);

  ierr = TCOPFLOWModelRegister(tcopflow, TCOPFLOWMODEL_GENRAMP,
                               TCOPFLOWModelCreate_GENRAMP);
  CHKERRQ(ierr);

  tcopflow->TCOPFLOWModelRegisterAllCalled = PETSC_TRUE;

  PetscFunctionReturn(0);
}

/*
  TCOPFLOWSolverRegister - Registers an TCOPFLOW formulation

  Input Parameters:
+ sname     - formuation name (string)
- createfunction  - the class constructor
*/
PetscErrorCode
TCOPFLOWSolverRegister(TCOPFLOW tcopflow, const char sname[],
                       PetscErrorCode (*createfunction)(TCOPFLOW)) {
  PetscErrorCode ierr;
  PetscInt i;
  PetscFunctionBegin;
  for (i = 0; i < tcopflow->nsolversregistered; i++) {
    PetscBool match;
    ierr = PetscStrcmp(tcopflow->TCOPFLOWSolverList[i].name, sname, &match);
    if (match)
      PetscFunctionReturn(0);
  }
  i = tcopflow->nsolversregistered;
  ierr = PetscStrcpy(tcopflow->TCOPFLOWSolverList[i].name, sname);
  CHKERRQ(ierr);
  tcopflow->TCOPFLOWSolverList[i].create = createfunction;
  tcopflow->nsolversregistered++;
  PetscFunctionReturn(0);
}

#if defined(EXAGO_ENABLE_IPOPT)
extern PetscErrorCode TCOPFLOWSolverCreate_IPOPT(TCOPFLOW);
#endif

/*
  TCOPFLOWSolverRegisterAll - Registers all TCOPFLOW solvers
*/
PetscErrorCode TCOPFLOWSolverRegisterAll(TCOPFLOW tcopflow) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if (tcopflow->TCOPFLOWSolverRegisterAllCalled)
    PetscFunctionReturn(0);

#if defined(EXAGO_ENABLE_IPOPT)
  ierr = TCOPFLOWSolverRegister(tcopflow, TCOPFLOWSOLVER_IPOPT,
                                TCOPFLOWSolverCreate_IPOPT);
  CHKERRQ(ierr);
#endif

  tcopflow->TCOPFLOWSolverRegisterAllCalled = PETSC_TRUE;

  PetscFunctionReturn(0);
}
