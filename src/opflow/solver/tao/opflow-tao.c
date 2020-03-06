#include <private/opflowimpl.h>
#include "opflow-tao.h"

PetscErrorCode OPFLOWObjectiveandGradientFunction_TAO(Tao nlp,Vec X,PetscScalar *obj,Vec grad,void *ctx)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=(OPFLOW)ctx;

  PetscFunctionBegin;
  *obj = 0.0;
  ierr = (*opflow->formops.computeobjandgradient)(opflow,X,obj,grad);CHKERRQ(ierr);
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
  PetscInt       gloc=opflow->nconeq,i;
  PetscScalar    *gi,*gu;

  PetscFunctionBegin;
  ierr = (*opflow->formops.computeinequalityconstraints)(opflow,X,Gi);CHKERRQ(ierr);
  ierr = VecGetArray(Gi,&gi);CHKERRQ(ierr);
  ierr = VecGetArray(opflow->Gu,&gu);CHKERRQ(ierr);
  for(i=0; i < opflow->nconineq; i++) {
    gi[i] = gu[gloc+i] - gi[i];
  }
  ierr = VecRestoreArray(Gi,&gi);CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->Gu,&gu);CHKERRQ(ierr);

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
  ierr = MatScale(Ji,-1.0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode OPFLOWHessian_TAO(Tao nlp,Vec X, Mat H, Mat H_pre, void *ctx)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=(OPFLOW)ctx;
  Vec            DE,DI;

  PetscFunctionBegin;
  ierr = TaoGetDualVariables(nlp,&DE,&DI);CHKERRQ(ierr);
  //  ierr = VecScale(opflow->Lambdai,-1.0);CHKERRQ(ierr);
  ierr = (*opflow->formops.computehessian)(opflow,X,DE,DI,H);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverSolve_TAO(OPFLOW opflow)
{
  PetscErrorCode     ierr;
  OPFLOWSolver_TAO   tao=(OPFLOWSolver_TAO)opflow->solver;

  PetscFunctionBegin;

  ierr = TaoSetInitialVector(tao->nlp,opflow->X);CHKERRQ(ierr);

  ierr = TaoSetVariableBounds(tao->nlp,opflow->Xl,opflow->Xu);CHKERRQ(ierr);

  ierr = TaoSolve(tao->nlp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverDestroy_TAO(OPFLOW opflow)
{
  PetscErrorCode     ierr;
  OPFLOWSolver_TAO tao = (OPFLOWSolver_TAO)opflow->solver;

  PetscFunctionBegin;

  ierr = VecDestroy(&tao->Glineq);CHKERRQ(ierr);
  ierr = VecDestroy(&tao->Guineq);CHKERRQ(ierr);

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
  /*  ierr = TaoSetType(tao->nlp,TAOPDIPM);CHKERRQ(ierr); */
  ierr = TaoSetOptionsPrefix(tao->nlp,"opflow_");CHKERRQ(ierr);

  /* Set Callback routines */
  /* Objective and gradient */
  ierr = TaoSetObjectiveAndGradientRoutine(tao->nlp,OPFLOWObjectiveandGradientFunction_TAO,(void*)opflow);CHKERRQ(ierr);

  /* Equality Constraints */
  ierr = TaoSetEqualityConstraintsRoutine(tao->nlp,opflow->Ge,OPFLOWEqualityConstraintsFunction_TAO,(void*)opflow);CHKERRQ(ierr);

  /* Inequality Constraints */
  ierr = TaoSetInequalityConstraintsRoutine(tao->nlp,opflow->Gi,OPFLOWInequalityConstraintsFunction_TAO,(void*)opflow);CHKERRQ(ierr);

  /* Equality Jacobian */
  ierr = TaoSetJacobianEqualityRoutine(tao->nlp,opflow->Jac_Ge,opflow->Jac_Ge,OPFLOWEqualityConstraintsJacobian_TAO,(void*)opflow);CHKERRQ(ierr);

  /* Inequality Jacobian */
  ierr = TaoSetJacobianInequalityRoutine(tao->nlp,opflow->Jac_Gi,opflow->Jac_Gi,OPFLOWInequalityConstraintsJacobian_TAO,(void*)opflow);CHKERRQ(ierr);
  ierr = TaoSetFromOptions(tao->nlp);CHKERRQ(ierr);

  /* Set Hessian routine */
  ierr = TaoSetHessianRoutine(tao->nlp,opflow->Hes,opflow->Hes,OPFLOWHessian_TAO,(void*)opflow);CHKERRQ(ierr);

  /* Create vectors for inequality constraint bounds */
  ierr = VecDuplicate(opflow->Gi,&tao->Glineq);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->Gu,&tao->Guineq);CHKERRQ(ierr);

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
