#include <private/opflowimpl.h>
#include "pbpol.h"

PetscErrorCode OPFLOWDestroy_PBPOL(OPFLOW opflow)
{
  PetscFunctionBegin;

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


PetscErrorCode OPFLOWFormulationCreate_PBPOL(OPFLOW opflow)
{
  PetscFunctionBegin;
  
  /* Inherit Ops */
  opflow->formops.destroy = OPFLOWDestroy_PBPOL;
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
