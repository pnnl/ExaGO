#include <exago_config.h>
#include <private/sopflowimpl.h>

/*
  SOPFLOWModelRegister - Registers a SOPFLOW model

  Input Parameters:
+ sname     - model name (string)
- createfunction  - the class constructor
*/
PetscErrorCode SOPFLOWModelRegister(SOPFLOW sopflow, const char sname[],
                                    PetscErrorCode (*createfunction)(SOPFLOW)) {
  PetscErrorCode ierr;
  PetscInt i;
  PetscFunctionBegin;
  for (i = 0; i < sopflow->nmodelsregistered; i++) {
    PetscBool match;
    ierr = PetscStrcmp(sopflow->SOPFLOWModelList[i].name, sname, &match);
    if (match)
      PetscFunctionReturn(0);
  }
  if (sopflow->nmodelsregistered == SOPFLOWMODELSMAX) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
             "Cannot register %s OPFLOW model, maximum limit %d reached\n",
             sname, SOPFLOWMODELSMAX);
  }
  i = sopflow->nmodelsregistered;
  ierr = PetscStrcpy(sopflow->SOPFLOWModelList[i].name, sname);
  CHKERRQ(ierr);
  sopflow->SOPFLOWModelList[i].create = createfunction;
  sopflow->nmodelsregistered++;
  PetscFunctionReturn(0);
}

extern PetscErrorCode SOPFLOWModelCreate_GENRAMP(SOPFLOW);
extern PetscErrorCode SOPFLOWModelCreate_GENRAMPC(SOPFLOW);

/*
  SOPFLOWModelRegisterAll - Registers all built-in SOPFLOW models
*/
PetscErrorCode SOPFLOWModelRegisterAll(SOPFLOW sopflow) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if (sopflow->SOPFLOWModelRegisterAllCalled)
    PetscFunctionReturn(0);

  ierr = SOPFLOWModelRegister(sopflow, SOPFLOWMODEL_GENRAMP,
                              SOPFLOWModelCreate_GENRAMP);
  CHKERRQ(ierr);
  ierr = SOPFLOWModelRegister(sopflow, SOPFLOWMODEL_GENRAMPC,
                              SOPFLOWModelCreate_GENRAMPC);
  CHKERRQ(ierr);

  sopflow->SOPFLOWModelRegisterAllCalled = PETSC_TRUE;

  PetscFunctionReturn(0);
}

/*
  SOPFLOWSolverRegister - Registers an SOPFLOW model

  Input Parameters:
+ sname           - model name (string)
- createfunction  - the class constructor
*/
PetscErrorCode
SOPFLOWSolverRegister(SOPFLOW sopflow, const char sname[],
                      PetscErrorCode (*createfunction)(SOPFLOW)) {
  PetscErrorCode ierr;
  PetscInt i;
  PetscFunctionBegin;
  for (i = 0; i < sopflow->nsolversregistered; i++) {
    PetscBool match;
    ierr = PetscStrcmp(sopflow->SOPFLOWSolverList[i].name, sname, &match);
    if (match)
      PetscFunctionReturn(0);
  }
  i = sopflow->nsolversregistered;
  ierr = PetscStrcpy(sopflow->SOPFLOWSolverList[i].name, sname);
  CHKERRQ(ierr);
  sopflow->SOPFLOWSolverList[i].create = createfunction;
  sopflow->nsolversregistered++;
  PetscFunctionReturn(0);
}

#if defined(EXAGO_ENABLE_IPOPT)
extern PetscErrorCode SOPFLOWSolverCreate_IPOPT(SOPFLOW);
#endif

extern PetscErrorCode SOPFLOWSolverCreate_EMPAR(SOPFLOW);

#if defined(EXAGO_ENABLE_HIOP)
extern PetscErrorCode SOPFLOWSolverCreate_HIOP(SOPFLOW);
#endif

/*
  SOPFLOWSolverRegisterAll - Registers all SOPFLOW solvers
*/
PetscErrorCode SOPFLOWSolverRegisterAll(SOPFLOW sopflow) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if (sopflow->SOPFLOWSolverRegisterAllCalled)
    PetscFunctionReturn(0);

#if defined(EXAGO_ENABLE_IPOPT)
  ierr = SOPFLOWSolverRegister(sopflow, SOPFLOWSOLVER_IPOPT,
                               SOPFLOWSolverCreate_IPOPT);
  CHKERRQ(ierr);
#endif

  ierr = SOPFLOWSolverRegister(sopflow, SOPFLOWSOLVER_EMPAR,
                               SOPFLOWSolverCreate_EMPAR);
  CHKERRQ(ierr);

#if defined(EXAGO_ENABLE_HIOP)
  ierr = SOPFLOWSolverRegister(sopflow, SOPFLOWSOLVER_HIOP,
                               SOPFLOWSolverCreate_HIOP);
  CHKERRQ(ierr);
#endif

  sopflow->SOPFLOWSolverRegisterAllCalled = PETSC_TRUE;

  PetscFunctionReturn(0);
}
