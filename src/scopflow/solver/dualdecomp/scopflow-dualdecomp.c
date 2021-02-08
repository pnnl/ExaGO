#include <exago_config.h>

#include <private/opflowimpl.h>
#include <private/tcopflowimpl.h>
#include <private/scopflowimpl.h>
#include "scopflow-dualdecomp.h"

PetscErrorCode SCOPFLOWSolverSolve_DUALDECOMP(SCOPFLOW scopflow)
{
  PetscErrorCode      ierr;
  PetscInt            c;
  TCOPFLOW            tcopflow;
  OPFLOW              opflow;

  PetscFunctionBegin;

  if(!scopflow->ismultiperiod) {
    for(c=0; c < scopflow->nc; c++) {
      opflow = scopflow->opflows[c];
      ierr   = OPFLOWSolve(opflow);
    }
  } else {
    for(c=0; c < scopflow->nc; c++) {
      tcopflow = scopflow->tcopflows[c];
      ierr   = TCOPFLOWSolve(tcopflow);
    }
  }    

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverDestroy_DUALDECOMP(SCOPFLOW scopflow)
{
  PetscErrorCode     ierr;
  SCOPFLOWSolver_DUALDECOMP dualdecomp = (SCOPFLOWSolver_DUALDECOMP)scopflow->solver;

  PetscFunctionBegin;

  ierr = PetscFree(dualdecomp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetObjective_DUALDECOMP(SCOPFLOW scopflow,PetscReal *obj)
{
  PetscErrorCode ierr;
  PetscReal      temp;
  PetscFunctionBegin;
  if(!scopflow->comm->rank) {
    if(!scopflow->ismultiperiod) {
      temp = scopflow->opflows[0]->obj;
    } else {
      temp = scopflow->tcopflows[0]->obj;
    }
  }
  // ierr = MPI_Bcast(&temp,1,MPI_REAL,0,scopflow->comm->type);CHKERRQ(ierr);
  ierr = MPI_Bcast(&temp,1,MPI_DOUBLE,0,scopflow->comm->type);CHKERRQ(ierr);
  scopflow->obj = temp;
  *obj = scopflow->obj;
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetSolution_DUALDECOMP(SCOPFLOW scopflow,PetscInt cont_num,Vec *X)
{
  TCOPFLOW       tcopflow;
  OPFLOW         opflow;

  PetscFunctionBegin;

  if(scopflow->cstart <= cont_num && cont_num < scopflow->cend) {
    if(!scopflow->ismultiperiod) {
      opflow = scopflow->opflows[cont_num-scopflow->cstart];
      *X = opflow->X; 
    } else {
      tcopflow = scopflow->tcopflows[cont_num-scopflow->cstart];
      *X = tcopflow->X; 
    }
  } else {
    *X = NULL;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConstraints_DUALDECOMP(SCOPFLOW scopflow,PetscInt cont_num,Vec *G)
{
  TCOPFLOW            tcopflow;
  OPFLOW              opflow;

  PetscFunctionBegin;

  if(scopflow->cstart <= cont_num && cont_num < scopflow->cend) {
    if(!scopflow->ismultiperiod) {
      opflow = scopflow->opflows[cont_num-scopflow->cstart];
      *G = opflow->G; 
    } else {
      tcopflow = scopflow->tcopflows[cont_num-scopflow->cstart];
      *G = tcopflow->G; 
    }      
  } else {
    *G = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConstraintMultipliers_DUALDECOMP(SCOPFLOW scopflow,PetscInt cont_num,Vec *Lambda)
{
  TCOPFLOW            tcopflow;
  OPFLOW              opflow;

  PetscFunctionBegin;
  if(scopflow->cstart <= cont_num && cont_num < scopflow->cend) {
    if(!scopflow->ismultiperiod) {
      opflow = scopflow->opflows[cont_num-scopflow->cstart];
      *Lambda = opflow->Lambda; 
    } else {
      tcopflow = scopflow->tcopflows[cont_num-scopflow->cstart];
      *Lambda = tcopflow->Lambda; 
    }      
  } else {
    *Lambda = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConvergenceStatus_DUALDECOMP(SCOPFLOW scopflow,PetscBool *status)
{
  PetscFunctionBegin;

  *status = PETSC_TRUE;

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverSetUp_DUALDECOMP(SCOPFLOW scopflow)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverCreate_DUALDECOMP(SCOPFLOW scopflow)
{
  PetscErrorCode ierr;
  SCOPFLOWSolver_DUALDECOMP dualdecomp;
  
  PetscFunctionBegin;

  ierr = PetscCalloc1(1,&dualdecomp);CHKERRQ(ierr);

  scopflow->solver = dualdecomp;

  scopflow->solverops.setup = SCOPFLOWSolverSetUp_DUALDECOMP;
  scopflow->solverops.solve = SCOPFLOWSolverSolve_DUALDECOMP;
  scopflow->solverops.destroy = SCOPFLOWSolverDestroy_DUALDECOMP;
  scopflow->solverops.getobjective = SCOPFLOWSolverGetObjective_DUALDECOMP;
  scopflow->solverops.getsolution  = SCOPFLOWSolverGetSolution_DUALDECOMP;
  scopflow->solverops.getconvergencestatus = SCOPFLOWSolverGetConvergenceStatus_DUALDECOMP;
  scopflow->solverops.getconstraints = SCOPFLOWSolverGetConstraints_DUALDECOMP;
  scopflow->solverops.getconstraintmultipliers = SCOPFLOWSolverGetConstraintMultipliers_DUALDECOMP;

  PetscFunctionReturn(0);
}
