#include <private/opflowimpl.h>
#include "pbpol.h"

PetscErrorCode OPFLOWFormulationDestroy_PBPOL(OPFLOW opflow)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscFree(opflow->formulation);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetVariableBounds_PBPOL(OPFLOW opflow,Vec Xl,Vec Xu)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetConstraintBounds_PBPOL(OPFLOW opflow,Vec Gl,Vec Gu)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetInitialGuess_PBPOL(OPFLOW opflow,Vec X)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeEqualityConstraints_PBPOL(OPFLOW opflow,Vec X,Vec Ge)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeInequalityConstraints_PBPOL(OPFLOW opflow,Vec X,Vec Gi)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeConstraints_PBPOL(OPFLOW opflow,Vec X,Vec G)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeObjandGradient_PBPOL(OPFLOW opflow,Vec X,PetscScalar *obj,Vec Grad)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeObjective_PBPOL(OPFLOW opflow,Vec X,PetscScalar *obj)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeGradient_PBPOL(OPFLOW opflow,Vec X,Vec Grad)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWFormulationSetNumVariables_PBPOL(OPFLOW opflow,PetscInt *busnvar,PetscInt *branchnvar)
{
  PetscInt i,ngen;
  PS       ps=opflow->ps;
  PSBUS    bus;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  /* No variables for the branches */
  for(i=0; i < ps->nbranch; i++) {
    branchnvar[i] = 0;
  }

  /* Variables for the buses */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetNGen(bus,&ngen);CHKERRQ(ierr);
    /* Number of variables = 2 + 2*ngen (2 variables for voltage + Pg, Qg for each gen) */
    busnvar[i] = 2 + 2*ngen;
  }

  PetscFunctionReturn(0);
}


PetscErrorCode OPFLOWFormulationCreate_PBPOL(OPFLOW opflow)
{
  PBPOL pbpol;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  
  ierr = PetscCalloc1(1,&pbpol);CHKERRQ(ierr);

  opflow->formulation = pbpol;

  /* Inherit Ops */
  opflow->formops.destroy = OPFLOWFormulationDestroy_PBPOL;
  opflow->formops.setnumvariables = OPFLOWFormulationSetNumVariables_PBPOL;
  opflow->formops.setvariablebounds = OPFLOWSetVariableBounds_PBPOL;
  opflow->formops.setconstraintbounds = OPFLOWSetConstraintBounds_PBPOL;
  opflow->formops.setinitialguess = OPFLOWSetInitialGuess_PBPOL;
  opflow->formops.computeequalityconstraints = OPFLOWComputeEqualityConstraints_PBPOL;
  opflow->formops.computeinequalityconstraints = OPFLOWComputeInequalityConstraints_PBPOL;
  opflow->formops.computeobjandgradient = OPFLOWComputeObjandGradient_PBPOL;
  opflow->formops.computeobjective = OPFLOWComputeObjective_PBPOL;
  opflow->formops.computegradient  = OPFLOWComputeGradient_PBPOL;
  
  PetscFunctionReturn(0);
}
