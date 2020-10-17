#include <exago_config.h>

#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>
#include "scopflow-empar.h"

PetscErrorCode SCOPFLOWSolverSolve_EMPAR(SCOPFLOW scopflow)
{
  PetscErrorCode      ierr;
  PetscInt            c;
  OPFLOW              opflow;

  PetscFunctionBegin;


  for(c=0; c < scopflow->nc; c++) {
    opflow = scopflow->opflows[c];
    ierr   = OPFLOWSolve(opflow);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverDestroy_EMPAR(SCOPFLOW scopflow)
{
  PetscErrorCode     ierr;
  SCOPFLOWSolver_EMPAR empar = (SCOPFLOWSolver_EMPAR)scopflow->solver;

  PetscFunctionBegin;

  ierr = PetscFree(empar);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetObjective_EMPAR(SCOPFLOW scopflow,PetscReal *obj)
{
  PetscErrorCode ierr;
  PetscReal      temp;
  PetscFunctionBegin;
  if(!scopflow->comm->rank) {
    temp = scopflow->opflows[0]->obj;
  }
  ierr = MPI_Bcast(&temp,1,MPI_REAL,0,scopflow->comm->type);CHKERRQ(ierr);
  scopflow->obj = temp;
  *obj = scopflow->obj;
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetSolution_EMPAR(SCOPFLOW scopflow,PetscInt cont_num,Vec *X)
{
  OPFLOW         opflow;

  PetscFunctionBegin;

  if(scopflow->cstart <= cont_num && cont_num < scopflow->cend) {
    opflow = scopflow->opflows[cont_num-scopflow->cstart];
    *X = opflow->X; 
  } else {
    *X = NULL;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConstraints_EMPAR(SCOPFLOW scopflow,PetscInt cont_num,Vec *G)
{
  OPFLOW              opflow;

  PetscFunctionBegin;

  if(scopflow->cstart <= cont_num && cont_num < scopflow->cend) {
    opflow = scopflow->opflows[cont_num-scopflow->cstart];
    *G = opflow->G; 
  } else {
    *G = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConstraintMultipliers_EMPAR(SCOPFLOW scopflow,PetscInt cont_num,Vec *Lambda)
{
  OPFLOW              opflow;

  PetscFunctionBegin;

  if(scopflow->cstart <= cont_num && cont_num < scopflow->cend) {
    opflow = scopflow->opflows[cont_num-scopflow->cstart];
    *Lambda = opflow->Lambda; 
  } else {
    *Lambda = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConvergenceStatus_EMPAR(SCOPFLOW scopflow,PetscBool *status)
{
  PetscFunctionBegin;

  *status = PETSC_TRUE;

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverSetUp_EMPAR(SCOPFLOW scopflow)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverCreate_EMPAR(SCOPFLOW scopflow)
{
  PetscErrorCode ierr;
  SCOPFLOWSolver_EMPAR empar;
  
  PetscFunctionBegin;

  ierr = PetscCalloc1(1,&empar);CHKERRQ(ierr);

  scopflow->solver = empar;

  scopflow->solverops.setup = SCOPFLOWSolverSetUp_EMPAR;
  scopflow->solverops.solve = SCOPFLOWSolverSolve_EMPAR;
  scopflow->solverops.destroy = SCOPFLOWSolverDestroy_EMPAR;
  scopflow->solverops.getobjective = SCOPFLOWSolverGetObjective_EMPAR;
  scopflow->solverops.getsolution  = SCOPFLOWSolverGetSolution_EMPAR;
  scopflow->solverops.getconvergencestatus = SCOPFLOWSolverGetConvergenceStatus_EMPAR;
  scopflow->solverops.getconstraints = SCOPFLOWSolverGetConstraints_EMPAR;
  scopflow->solverops.getconstraintmultipliers = SCOPFLOWSolverGetConstraintMultipliers_EMPAR;

  PetscFunctionReturn(0);
}
