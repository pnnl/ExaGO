#include <private/opflowimpl.h>
#include "opflow-tao.h"

PetscErrorCode OPFLOWObjectiveandGradientFunction_TAO(Tao nlp,Vec X,PetscScalar *obj,Vec grad,void *ctx)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=(OPFLOW)ctx;

  PetscFunctionBegin;
  *obj = 0.0;
  ierr = (*opflow->formops.computeobjandgradient)(opflow,opflow->X,obj,grad);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode OPFLOWEqualityConstraintsFunction_TAO(Tao nlp,Vec X,Vec Ge,void *ctx)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=(OPFLOW)ctx;

  PetscFunctionBegin;
  ierr = (*opflow->formops.computeequalityconstraints)(opflow,X,Ge);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWInequalityConstraintsFunction_TAO(Tao nlp,Vec X,Vec Gi,void *ctx)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=(OPFLOW)ctx;

  PetscFunctionBegin;
  ierr = (*opflow->formops.computeinequalityconstraints)(opflow,X,Gi);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWEqualityConstraintsJacobian_TAO(Tao nlp, Vec X, Mat Je, Mat Je_pre, void *ctx)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=(OPFLOW)ctx;

  PetscFunctionBegin;
  ierr = (*opflow->formops.computeequalityconstraintjacobian)(opflow,X,Je);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWInequalityConstraintsJacobian_TAO(Tao nlp, Vec X, Mat Ji, Mat Ji_pre, void *ctx)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=(OPFLOW)ctx;

  PetscFunctionBegin;
  ierr = (*opflow->formops.computeinequalityconstraintjacobian)(opflow,X,Ji);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode OPFLOWHessian_TAO(Tao nlp,Vec X, Mat H, Mat H_pre, void *ctx)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=(OPFLOW)ctx;

  PetscFunctionBegin;
  ierr = (*opflow->formops.computehessian)(opflow,X,opflow->Lambda,H);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverSolve_TAO(OPFLOW opflow)
{
  PetscErrorCode     ierr;
  OPFLOWSolver_TAO tao=opflow->solver;

  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverDestroy_TAO(OPFLOW opflow)
{
  PetscErrorCode     ierr;
  OPFLOWSolver_TAO tao=opflow->solver;

  PetscFunctionBegin;

  if(tao->nlp) {
    TaoDestroy(&tao->nlp);
  }

  ierr = PetscFree(tao);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverSetUp_TAO(OPFLOW opflow)
{
  PetscErrorCode   ierr;
  OPFLOWSolver_TAO tao=(OPFLOWSolver_TAO)opflow->solver;
  PetscFunctionBegin;

  /* Create Tao solver and set type */
  ierr = TaoCreate(opflow->comm->type,&tao->nlp);CHKERRQ(ierr);
  ierr = TaoSetType(tao->nlp,TAOPDIPM);CHKERRQ(ierr);
  ierr = TaoSetOptionsPrefix(tao->nlp,"opflow_");CHKERRQ(ierr);

  /* Set Callback routines */

  /* Objective and gradient */
  ierr = TaoSetObjectiveAndGradientRoutine(tao->nlp,OPFLOWObjectiveandGradientFunction_TAO,(void*)opflow);CHKERRQ(ierr);

  /* Equality Constraints */
  ierr = TaoSetEqualityConstraintsRoutine(tao->nlp,opflow->Ge,OPFLOWEqualityConstraintsFunction_TAO,(void*)opflow);CHKERRQ(ierr);

  /* Inequality Constraints */
  ierr = TaoSetInequalityConstraintsRoutine(opflow->nlp,opflow->Gi,OPFLOWInequalityConstraintsFunction_TAO,(void*)opflow);CHKERRQ(ierr);

  /* Equality Jacobian */
  ierr = TaoSetJacobianEqualityRoutine(opflow->nlp,opflow->Jac_Ge,opflow->Jac_Ge,OPFLOWEqualityConstraintsJacobian,(void*)opflow);CHKERRQ(ierr);

  /* Inequality Jacobian */
  ierr = TaoSetJacobianInequalityRoutine(opflow->nlp,opflow->Jac_Gi,opflow->Jac_Gi,OPFLOWInequalityConstraintsJacobian,(void*)opflow);CHKERRQ(ierr);
  ierr = TaoSetFromOptions(opflow->nlp);CHKERRQ(ierr);

  ierr = TaoSetHessianRoutine(opflow->nlp,NULL,NULL,OPFLOWHessian,NULL);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverCreate_TAO(OPFLOW opflow)
{
  PetscErrorCode ierr;
  OPFLOWSolver_TAO tao;
  
  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&tao);CHKERRQ(ierr);

  tao->nlp = NULL;
  opflow->solver = tao;

  opflow->solverops.setup = OPFLOWSolverSetUp_TAO;
  opflow->solverops.solve = OPFLOWSolverSolve_TAO;
  opflow->solverops.destroy = OPFLOWSolverDestroy_TAO;

  PetscFunctionReturn(0);
}
