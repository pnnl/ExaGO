#include <exago_config.h>

#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>
#include <private/sopflowimpl.h>
#include "sopflow-empar.h"

PetscErrorCode SOPFLOWSolverSolve_EMPAR(SOPFLOW sopflow)
{
  PetscErrorCode      ierr;
  PetscInt            s;
  SCOPFLOW            scopflow;
  OPFLOW              opflow;

  PetscFunctionBegin;

  if(!sopflow->ismulticontingency) {
    for(s=0; s < sopflow->ns; s++) {
      opflow = sopflow->opflows[s];
      ierr   = OPFLOWSolve(opflow);
    }
  } else {
    for(s=0; s < sopflow->ns; s++) {
      scopflow = sopflow->scopflows[s];
      ierr   = SCOPFLOWSolve(scopflow);
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverDestroy_EMPAR(SOPFLOW sopflow)
{
  PetscErrorCode     ierr;
  SOPFLOWSolver_EMPAR empar = (SOPFLOWSolver_EMPAR)sopflow->solver;

  PetscFunctionBegin;

  ierr = PetscFree(empar);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetObjective_EMPAR(SOPFLOW sopflow,PetscReal *obj)
{
  PetscErrorCode ierr;
  PetscReal      temp;
  PetscFunctionBegin;
  if(!sopflow->comm->rank) {
    if(!sopflow->ismulticontingency) {
      temp = sopflow->opflows[0]->obj;
    } else {
      temp = sopflow->scopflows[0]->obj;
    }
  }

  // ierr = MPI_Bcast(&temp,1,MPI_REAL,0,sopflow->comm->type);CHKERRQ(ierr);
  ierr = MPI_Bcast(&temp,1,MPI_DOUBLE,0,sopflow->comm->type);CHKERRQ(ierr);
  sopflow->obj = temp;
  *obj = sopflow->obj;
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetSolution_EMPAR(SOPFLOW sopflow,PetscInt scen_num,Vec *X)
{
  SCOPFLOW       scopflow;
  OPFLOW         opflow;

  PetscFunctionBegin;

  if(sopflow->sstart <= scen_num && scen_num < sopflow->send) {
    if(!sopflow->ismulticontingency) {
      opflow = sopflow->opflows[scen_num-sopflow->sstart];
      *X = opflow->X; 
    } else {
      scopflow = sopflow->scopflows[scen_num-sopflow->sstart];
      *X = scopflow->X; 
    }
  } else {
    *X = NULL;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetConstraints_EMPAR(SOPFLOW sopflow,PetscInt scen_num,Vec *G)
{
  SCOPFLOW            scopflow;
  OPFLOW              opflow;

  PetscFunctionBegin;

  if(sopflow->sstart <= scen_num && scen_num < sopflow->send) {
    if(!sopflow->ismulticontingency) {
      opflow = sopflow->opflows[scen_num-sopflow->sstart];
      *G = opflow->G; 
    } else {
      scopflow = sopflow->scopflows[scen_num-sopflow->sstart];
      *G = scopflow->G; 
    }
  } else {
    *G = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetConstraintMultipliers_EMPAR(SOPFLOW sopflow,PetscInt scen_num,Vec *Lambda)
{
  SCOPFLOW            scopflow;
  OPFLOW              opflow;

  PetscFunctionBegin;
  if(sopflow->sstart <= scen_num && scen_num < sopflow->send) {
    if(!sopflow->ismulticontingency) {
      opflow = sopflow->opflows[scen_num-sopflow->sstart];
      *Lambda = opflow->Lambda; 
    } else {
      scopflow = sopflow->scopflows[scen_num-sopflow->sstart];
      *Lambda = scopflow->Lambda; 
    }
  } else {
    *Lambda = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetConvergenceStatus_EMPAR(SOPFLOW sopflow,PetscBool *status)
{
  PetscFunctionBegin;

  *status = PETSC_TRUE;

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverSetUp_EMPAR(SOPFLOW sopflow)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverCreate_EMPAR(SOPFLOW sopflow)
{
  PetscErrorCode ierr;
  SOPFLOWSolver_EMPAR empar;
  
  PetscFunctionBegin;

  ierr = PetscCalloc1(1,&empar);CHKERRQ(ierr);

  sopflow->solver = empar;

  sopflow->solverops.setup = SOPFLOWSolverSetUp_EMPAR;
  sopflow->solverops.solve = SOPFLOWSolverSolve_EMPAR;
  sopflow->solverops.destroy = SOPFLOWSolverDestroy_EMPAR;
  sopflow->solverops.getobjective = SOPFLOWSolverGetObjective_EMPAR;
  sopflow->solverops.getsolution  = SOPFLOWSolverGetSolution_EMPAR;
  sopflow->solverops.getconvergencestatus = SOPFLOWSolverGetConvergenceStatus_EMPAR;
  sopflow->solverops.getconstraints = SOPFLOWSolverGetConstraints_EMPAR;
  sopflow->solverops.getconstraintmultipliers = SOPFLOWSolverGetConstraintMultipliers_EMPAR;

  PetscFunctionReturn(0);
}
