#include <exago_config.h>

#include "sopflow_empar.h"
#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>
#include <private/sopflowimpl.h>

PetscErrorCode SOPFLOWSolverSolve_EMPAR(SOPFLOW sopflow) {
  PetscErrorCode ierr;
  PetscInt s;
  SCOPFLOW scopflow;
  OPFLOW opflow;

  PetscFunctionBegin;

  if (!sopflow->ismulticontingency) {
    for (s = 0; s < sopflow->ns; s++) {
      opflow = sopflow->opflows[s];
      ierr = OPFLOWSolve(opflow);
    }
  } else {
    for (s = 0; s < sopflow->ns; s++) {
      scopflow = sopflow->scopflows[s];
      ierr = SCOPFLOWSolve(scopflow);
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverDestroy_EMPAR(SOPFLOW sopflow) {
  PetscErrorCode ierr;
  SOPFLOWSolver_EMPAR empar = (SOPFLOWSolver_EMPAR)sopflow->solver;

  PetscFunctionBegin;

  ierr = PetscFree(empar);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetBaseObjective_EMPAR(SOPFLOW sopflow,
                                                   PetscReal *obj) {
  PetscErrorCode ierr;
  PetscReal temp = 0.0;
  PetscFunctionBegin;
  if (!sopflow->comm->rank) {
    if (!sopflow->ismulticontingency) {
      temp = sopflow->opflows[0]->obj;
    } else {
      temp = sopflow->scopflows[0]->objbase;
    }
  }

  // ierr = MPI_Bcast(&temp,1,MPI_REAL,0,sopflow->comm->type);CHKERRQ(ierr);
  ierr = MPI_Bcast(&temp, 1, MPI_DOUBLE, 0, sopflow->comm->type);
  CHKERRQ(ierr);
  sopflow->objbase = temp;
  *obj = sopflow->objbase;
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetTotalObjective_EMPAR(SOPFLOW sopflow,
                                                    PetscReal *obj) {
  PetscErrorCode ierr;
  PetscReal temp = 0.0;
  PetscFunctionBegin;
  if (!sopflow->ismulticontingency) {
    for (int i = 0; i < sopflow->ns; i++) {
      temp += sopflow->opflows[i]->obj;
    }
  } else {
    temp = sopflow->scopflows[0]->objtot;
  }

  // ierr = MPI_Bcast(&temp,1,MPI_REAL,0,sopflow->comm->type);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&temp, &sopflow->objtot, 1, MPI_DOUBLE, MPI_SUM,
                       sopflow->comm->type);
  CHKERRQ(ierr);
  *obj = sopflow->objtot;
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetSolution_EMPAR(SOPFLOW sopflow,
                                              PetscInt scen_num, Vec *X) {
  SCOPFLOW scopflow;
  OPFLOW opflow;

  PetscFunctionBegin;

  if (sopflow->sstart <= scen_num && scen_num < sopflow->send) {
    if (!sopflow->ismulticontingency) {
      opflow = sopflow->opflows[scen_num - sopflow->sstart];
      *X = opflow->X;
    } else {
      scopflow = sopflow->scopflows[scen_num - sopflow->sstart];
      *X = scopflow->X;
    }
  } else {
    *X = NULL;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetConstraints_EMPAR(SOPFLOW sopflow,
                                                 PetscInt scen_num, Vec *G) {
  SCOPFLOW scopflow;
  OPFLOW opflow;

  PetscFunctionBegin;

  if (sopflow->sstart <= scen_num && scen_num < sopflow->send) {
    if (!sopflow->ismulticontingency) {
      opflow = sopflow->opflows[scen_num - sopflow->sstart];
      *G = opflow->G;
    } else {
      scopflow = sopflow->scopflows[scen_num - sopflow->sstart];
      *G = scopflow->G;
    }
  } else {
    *G = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetConstraintMultipliers_EMPAR(SOPFLOW sopflow,
                                                           PetscInt scen_num,
                                                           Vec *Lambda) {
  SCOPFLOW scopflow;
  OPFLOW opflow;

  PetscFunctionBegin;
  if (sopflow->sstart <= scen_num && scen_num < sopflow->send) {
    if (!sopflow->ismulticontingency) {
      opflow = sopflow->opflows[scen_num - sopflow->sstart];
      *Lambda = opflow->Lambda;
    } else {
      scopflow = sopflow->scopflows[scen_num - sopflow->sstart];
      *Lambda = scopflow->Lambda;
    }
  } else {
    *Lambda = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetConvergenceStatus_EMPAR(SOPFLOW sopflow,
                                                       PetscBool *status) {
  PetscFunctionBegin;

  *status = PETSC_TRUE;

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverSetUp_EMPAR(SOPFLOW sopflow) {
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverCreate_EMPAR(SOPFLOW sopflow) {
  PetscErrorCode ierr;
  SOPFLOWSolver_EMPAR empar;

  PetscFunctionBegin;

  ierr = PetscCalloc1(1, &empar);
  CHKERRQ(ierr);

  sopflow->solver = empar;

  sopflow->solverops.setup = SOPFLOWSolverSetUp_EMPAR;
  sopflow->solverops.solve = SOPFLOWSolverSolve_EMPAR;
  sopflow->solverops.destroy = SOPFLOWSolverDestroy_EMPAR;
  sopflow->solverops.getbaseobjective = SOPFLOWSolverGetBaseObjective_EMPAR;
  sopflow->solverops.gettotalobjective = SOPFLOWSolverGetTotalObjective_EMPAR;
  sopflow->solverops.getsolution = SOPFLOWSolverGetSolution_EMPAR;
  sopflow->solverops.getconvergencestatus =
      SOPFLOWSolverGetConvergenceStatus_EMPAR;
  sopflow->solverops.getconstraints = SOPFLOWSolverGetConstraints_EMPAR;
  sopflow->solverops.getconstraintmultipliers =
      SOPFLOWSolverGetConstraintMultipliers_EMPAR;

  PetscFunctionReturn(0);
}
