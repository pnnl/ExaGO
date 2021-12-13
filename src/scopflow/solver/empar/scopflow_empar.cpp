#include <exago_config.h>

#include "scopflow_empar.h"
#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>
#include <private/tcopflowimpl.h>

PetscErrorCode SCOPFLOWSolverSolve_EMPAR(SCOPFLOW scopflow) {
  PetscErrorCode ierr;
  PetscInt c;
  TCOPFLOW tcopflow;
  OPFLOW opflow;

  PetscFunctionBegin;

  if (!scopflow->ismultiperiod) {
    for (c = 0; c < scopflow->nc; c++) {
      opflow = scopflow->opflows[c];
      ierr = OPFLOWSolve(opflow);
    }
  } else {
    for (c = 0; c < scopflow->nc; c++) {
      tcopflow = scopflow->tcopflows[c];
      ierr = TCOPFLOWSolve(tcopflow);
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverDestroy_EMPAR(SCOPFLOW scopflow) {
  PetscErrorCode ierr;
  SCOPFLOWSolver_EMPAR empar = (SCOPFLOWSolver_EMPAR)scopflow->solver;

  PetscFunctionBegin;

  ierr = PetscFree(empar);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetObjective_EMPAR(SCOPFLOW scopflow,
                                                PetscReal *obj) {
  PetscErrorCode ierr;
  PetscReal temp;
  PetscFunctionBegin;
  if (!scopflow->comm->rank) {
    if (!scopflow->ismultiperiod) {
      temp = scopflow->opflows[0]->obj;
    } else {
      temp = scopflow->tcopflows[0]->obj;
    }
  }
  // ierr = MPI_Bcast(&temp,1,MPI_REAL,0,scopflow->comm->type);CHKERRQ(ierr);
  ierr = MPI_Bcast(&temp, 1, MPI_DOUBLE, 0, scopflow->comm->type);
  CHKERRQ(ierr);
  scopflow->obj = temp;
  *obj = scopflow->obj;
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetSolution_EMPAR(SCOPFLOW scopflow,
                                               PetscInt cont_num, Vec *X) {
  TCOPFLOW tcopflow;
  OPFLOW opflow;

  PetscFunctionBegin;

  if (scopflow->cstart <= cont_num && cont_num < scopflow->cend) {
    if (!scopflow->ismultiperiod) {
      opflow = scopflow->opflows[cont_num - scopflow->cstart];
      *X = opflow->X;
    } else {
      tcopflow = scopflow->tcopflows[cont_num - scopflow->cstart];
      *X = tcopflow->X;
    }
  } else {
    *X = NULL;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConstraints_EMPAR(SCOPFLOW scopflow,
                                                  PetscInt cont_num, Vec *G) {
  TCOPFLOW tcopflow;
  OPFLOW opflow;

  PetscFunctionBegin;

  if (scopflow->cstart <= cont_num && cont_num < scopflow->cend) {
    if (!scopflow->ismultiperiod) {
      opflow = scopflow->opflows[cont_num - scopflow->cstart];
      *G = opflow->G;
    } else {
      tcopflow = scopflow->tcopflows[cont_num - scopflow->cstart];
      *G = tcopflow->G;
    }
  } else {
    *G = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConstraintMultipliers_EMPAR(SCOPFLOW scopflow,
                                                            PetscInt cont_num,
                                                            Vec *Lambda) {
  TCOPFLOW tcopflow;
  OPFLOW opflow;

  PetscFunctionBegin;
  if (scopflow->cstart <= cont_num && cont_num < scopflow->cend) {
    if (!scopflow->ismultiperiod) {
      opflow = scopflow->opflows[cont_num - scopflow->cstart];
      *Lambda = opflow->Lambda;
    } else {
      tcopflow = scopflow->tcopflows[cont_num - scopflow->cstart];
      *Lambda = tcopflow->Lambda;
    }
  } else {
    *Lambda = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConvergenceStatus_EMPAR(SCOPFLOW scopflow,
                                                        PetscBool *status) {
  PetscFunctionBegin;

  *status = PETSC_TRUE;

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverSetUp_EMPAR(SCOPFLOW scopflow) {
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverCreate_EMPAR(SCOPFLOW scopflow) {
  PetscErrorCode ierr;
  SCOPFLOWSolver_EMPAR empar;

  PetscFunctionBegin;

  ierr = PetscCalloc1(1, &empar);
  CHKERRQ(ierr);

  scopflow->solver = empar;

  scopflow->solverops.setup = SCOPFLOWSolverSetUp_EMPAR;
  scopflow->solverops.solve = SCOPFLOWSolverSolve_EMPAR;
  scopflow->solverops.destroy = SCOPFLOWSolverDestroy_EMPAR;
  scopflow->solverops.getobjective = SCOPFLOWSolverGetObjective_EMPAR;
  scopflow->solverops.getsolution = SCOPFLOWSolverGetSolution_EMPAR;
  scopflow->solverops.getconvergencestatus =
      SCOPFLOWSolverGetConvergenceStatus_EMPAR;
  scopflow->solverops.getconstraints = SCOPFLOWSolverGetConstraints_EMPAR;
  scopflow->solverops.getconstraintmultipliers =
      SCOPFLOWSolverGetConstraintMultipliers_EMPAR;

  PetscFunctionReturn(0);
}
