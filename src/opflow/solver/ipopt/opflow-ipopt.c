#if defined(SCOPFLOW_HAVE_IPOPT)

#include <private/opflowimpl.h>
#include "opflow-ipopt.h"

PetscErrorCode OPFLOWSolverDestroy_IPOPT(OPFLOW opflow)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscFree(opflow->solver);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverSetUp_IPOPT(OPFLOW opflow)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  PetscFunctionReturn(0);
}


PetscErrorCode OPFLOWSolverCreate_IPOPT(OPFLOW opflow)
{
  PetscErrorCode ierr;
  OPFLOWSolver_IPOPT ipopt;
  
  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&ipopt);CHKERRQ(ierr);

  opflow->solver = ipopt;

  opflow->solverops.setup = OPFLOWSolverSetUp_IPOPT;
  opflow->solverops.destroy = OPFLOWSolverDestroy_IPOPT;


  PetscFunctionReturn(0);
}

#endif
