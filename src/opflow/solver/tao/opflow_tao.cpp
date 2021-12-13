#include "opflow_tao.h"
#include <private/opflowimpl.h>

PetscErrorCode OPFLOWObjectiveandGradientFunction_TAO(Tao nlp, Vec X,
                                                      PetscScalar *obj,
                                                      Vec grad, void *ctx) {
  PetscErrorCode ierr;
  OPFLOW opflow = (OPFLOW)ctx;

  PetscFunctionBegin;
  *obj = 0.0;
  ierr = (*opflow->modelops.computeobjandgradient)(opflow, X, obj, grad);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWEqualityConstraintsFunction_TAO(Tao nlp, Vec X, Vec Ge,
                                                     void *ctx) {
  PetscErrorCode ierr;
  OPFLOW opflow = (OPFLOW)ctx;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->eqconslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeequalityconstraints)(opflow, X, Ge);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->eqconslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWInequalityConstraintsFunction_TAO(Tao nlp, Vec X, Vec Gi,
                                                       void *ctx) {
  PetscErrorCode ierr;
  OPFLOW opflow = (OPFLOW)ctx;
  PetscInt gloc = opflow->nconeq, i;
  PetscScalar *gi, *gu;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->ineqconslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeinequalityconstraints)(opflow, X, Gi);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->ineqconslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = VecGetArray(Gi, &gi);
  CHKERRQ(ierr);
  ierr = VecGetArray(opflow->Gu, &gu);
  CHKERRQ(ierr);
  for (i = 0; i < opflow->nconineq; i++) {
    gi[i] = gu[gloc + i] - gi[i];
  }
  ierr = VecRestoreArray(Gi, &gi);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->Gu, &gu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWEqualityConstraintsJacobian_TAO(Tao nlp, Vec X, Mat Je,
                                                     Mat Je_pre, void *ctx) {
  PetscErrorCode ierr;
  OPFLOW opflow = (OPFLOW)ctx;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->eqconsjaclogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeequalityconstraintjacobian)(opflow, X, Je);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->eqconsjaclogger, 0, 0, 0, 0);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWInequalityConstraintsJacobian_TAO(Tao nlp, Vec X, Mat Ji,
                                                       Mat Ji_pre, void *ctx) {
  PetscErrorCode ierr;
  OPFLOW opflow = (OPFLOW)ctx;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->ineqconsjaclogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeinequalityconstraintjacobian)(opflow, X, Ji);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->ineqconsjaclogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = MatScale(Ji, -1.0);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWHessian_TAO(Tao nlp, Vec X, Mat H, Mat H_pre, void *ctx) {
  PetscErrorCode ierr;
  OPFLOW opflow = (OPFLOW)ctx;
  Vec DE, DI;

  PetscFunctionBegin;
  ierr = TaoGetDualVariables(nlp, &DE, &DI);
  CHKERRQ(ierr);
  ierr = PetscLogEventBegin(opflow->hesslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  //  ierr = VecScale(opflow->Lambdai,-1.0);CHKERRQ(ierr);
  ierr = (*opflow->modelops.computehessian)(opflow, X, DE, DI, H);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->hesslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverSolve_TAO(OPFLOW opflow) {
  PetscErrorCode ierr;
  OPFLOWSolver_TAO tao = (OPFLOWSolver_TAO)opflow->solver;

  PetscFunctionBegin;

  ierr = TaoSetInitialVector(tao->nlp, opflow->X);
  CHKERRQ(ierr);

  ierr = TaoSetVariableBounds(tao->nlp, opflow->Xl, opflow->Xu);
  CHKERRQ(ierr);

  ierr = TaoSolve(tao->nlp);
  CHKERRQ(ierr);

  /* Get the iteration count */
  ierr = TaoGetTotalIterationNumber(tao->nlp, &opflow->numits);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverDestroy_TAO(OPFLOW opflow) {
  PetscErrorCode ierr;
  OPFLOWSolver_TAO tao = (OPFLOWSolver_TAO)opflow->solver;

  PetscFunctionBegin;

  if (opflow->Nconineq) {
    ierr = VecDestroy(&tao->Glineq);
    CHKERRQ(ierr);
    ierr = VecDestroy(&tao->Guineq);
    CHKERRQ(ierr);
  }

  if (tao->nlp) {
    TaoDestroy(&tao->nlp);
  }

  ierr = PetscFree(tao);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverSetUp_TAO(OPFLOW opflow) {
  PetscErrorCode ierr;
  OPFLOWSolver_TAO tao = (OPFLOWSolver_TAO)opflow->solver;
  PetscFunctionBegin;

  /* Create Tao solver and set type */
  ierr = TaoCreate(opflow->comm->type, &tao->nlp);
  CHKERRQ(ierr);
  /*  ierr = TaoSetType(tao->nlp,TAOPDIPM);CHKERRQ(ierr); */
  ierr = TaoSetOptionsPrefix(tao->nlp, "opflow_");
  CHKERRQ(ierr);

  /* Set Callback routines */
  /* Objective and gradient */
  ierr = TaoSetObjectiveAndGradientRoutine(
      tao->nlp, OPFLOWObjectiveandGradientFunction_TAO, (void *)opflow);
  CHKERRQ(ierr);

  /* Equality Constraints */
  ierr = TaoSetEqualityConstraintsRoutine(tao->nlp, opflow->Ge,
                                          OPFLOWEqualityConstraintsFunction_TAO,
                                          (void *)opflow);
  CHKERRQ(ierr);

  /* Inequality Constraints */
  ierr = TaoSetInequalityConstraintsRoutine(
      tao->nlp, opflow->Gi, OPFLOWInequalityConstraintsFunction_TAO,
      (void *)opflow);
  CHKERRQ(ierr);

  /* Equality Jacobian */
  ierr = TaoSetJacobianEqualityRoutine(tao->nlp, opflow->Jac_Ge, opflow->Jac_Ge,
                                       OPFLOWEqualityConstraintsJacobian_TAO,
                                       (void *)opflow);
  CHKERRQ(ierr);

  /* Inequality Jacobian */
  ierr = TaoSetJacobianInequalityRoutine(
      tao->nlp, opflow->Jac_Gi, opflow->Jac_Gi,
      OPFLOWInequalityConstraintsJacobian_TAO, (void *)opflow);
  CHKERRQ(ierr);
  ierr = TaoSetFromOptions(tao->nlp);
  CHKERRQ(ierr);

  /* Set Hessian routine */
  ierr = TaoSetHessianRoutine(tao->nlp, opflow->Hes, opflow->Hes,
                              OPFLOWHessian_TAO, (void *)opflow);
  CHKERRQ(ierr);

  /* Create vectors for inequality constraint bounds */
  if (opflow->Nconineq) {
    ierr = VecDuplicate(opflow->Gi, &tao->Glineq);
    CHKERRQ(ierr);
    ierr = VecDuplicate(opflow->Gu, &tao->Guineq);
    CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetObjective_TAO(OPFLOW opflow, PetscReal *obj) {
  PetscErrorCode ierr;
  OPFLOWSolver_TAO tao = (OPFLOWSolver_TAO)opflow->solver;

  PetscFunctionBegin;
  ierr = TaoGetObjective(tao->nlp, obj);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetSolution_TAO(OPFLOW opflow, Vec *X) {
  PetscErrorCode ierr;
  OPFLOWSolver_TAO tao = (OPFLOWSolver_TAO)opflow->solver;

  PetscFunctionBegin;
  ierr = TaoGetSolutionVector(tao->nlp, X);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConstraints_TAO(OPFLOW opflow, Vec *G) {
  PetscErrorCode ierr;
  PetscScalar *g, *ge, *gi, *gu;
  PetscInt gloc = opflow->nconeq, i;

  PetscFunctionBegin;
  ierr = VecGetArray(opflow->G, &g);
  CHKERRQ(ierr);
  ierr = VecGetArray(opflow->Ge, &ge);
  CHKERRQ(ierr);

  ierr = PetscMemcpy(g, ge, opflow->nconeq * sizeof(PetscScalar));
  if (opflow->Nconineq) {
    ierr = VecGetArray(opflow->Gi, &gi);
    CHKERRQ(ierr);
    ierr = VecGetArray(opflow->Gu, &gu);
    CHKERRQ(ierr);
    /* TAO needs the equations in the form g(x) > 0 as opposed
       to gl <= g(x) <= gu which is used by IPOPT and others.
       So, the OPFLOW TAO interface uses gu - g(x) as the
       inequality constraint (See OPFLOWComputeInequalityConstraints_TAO
       above). So, we need to do the manipulation below to retrieve
       g(x)
    */
    for (i = 0; i < opflow->nconineq; i++) {
      g[gloc + i] = -(gi[i] - gu[gloc + i]);
    }
    ierr = VecRestoreArray(opflow->Gi, &gi);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(opflow->Gu, &gu);
    CHKERRQ(ierr);
  }
  ierr = VecRestoreArray(opflow->G, &g);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->Ge, &ge);
  CHKERRQ(ierr);

  *G = opflow->G;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConstraintMultipliers_TAO(OPFLOW opflow,
                                                        Vec *Lambda) {
  PetscErrorCode ierr;
  OPFLOWSolver_TAO tao = (OPFLOWSolver_TAO)opflow->solver;
  Vec DE, DI;
  PetscScalar *lambda, *lambdae, *lambdai;

  PetscFunctionBegin;
  ierr = TaoGetDualVariables(tao->nlp, &DE, &DI);
  CHKERRQ(ierr);
  ierr = VecGetArray(opflow->Lambda, &lambda);
  CHKERRQ(ierr);
  ierr = VecGetArray(DE, &lambdae);
  CHKERRQ(ierr);

  ierr = PetscMemcpy(lambda, lambdae, opflow->nconeq * sizeof(PetscScalar));
  if (opflow->Nconineq) {
    ierr = VecGetArray(DI, &lambdai);
    CHKERRQ(ierr);
    /* TAO needs the equations in the form g(x) > 0 as opposed
       to gl <= g(x) <= gu which is used by IPOPT and others.
       So, the OPFLOW TAO interface uses gu - g(x) as the
       inequality constraint (See OPFLOWComputeInequalityConstraints_TAO
       above). So, the lagrange multipliers obtained may be different than
       IPOPT.
    */
    ierr = PetscMemcpy(lambda + opflow->nconeq, lambdai,
                       opflow->nconineq * sizeof(PetscScalar));

    ierr = VecRestoreArray(DI, &lambdai);
    CHKERRQ(ierr);
  }
  ierr = VecRestoreArray(opflow->Lambda, &lambda);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(DE, &lambdae);
  CHKERRQ(ierr);

  *Lambda = opflow->Lambda;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConvergenceStatus_TAO(OPFLOW opflow,
                                                    PetscBool *status) {
  PetscErrorCode ierr;
  OPFLOWSolver_TAO tao = (OPFLOWSolver_TAO)opflow->solver;
  TaoConvergedReason convergedreason;

  PetscFunctionBegin;
  ierr = TaoGetConvergedReason(tao->nlp, &convergedreason);
  CHKERRQ(ierr);
  if (convergedreason > 0)
    *status = PETSC_TRUE; /* See IpReturnCodes_inc.h in IPOPT. The first two
                             denote convergence */
  else
    *status = PETSC_FALSE;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverCreate_TAO(OPFLOW opflow) {
  PetscErrorCode ierr;
  OPFLOWSolver_TAO tao;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1, &tao);
  CHKERRQ(ierr);

  tao->nlp = NULL;
  opflow->solver = tao;

  opflow->solverops.setup = OPFLOWSolverSetUp_TAO;
  opflow->solverops.solve = OPFLOWSolverSolve_TAO;
  opflow->solverops.destroy = OPFLOWSolverDestroy_TAO;
  opflow->solverops.getobjective = OPFLOWSolverGetObjective_TAO;
  opflow->solverops.getsolution = OPFLOWSolverGetSolution_TAO;
  opflow->solverops.getconvergencestatus = OPFLOWSolverGetConvergenceStatus_TAO;
  opflow->solverops.getconstraints = OPFLOWSolverGetConstraints_TAO;
  opflow->solverops.getconstraintmultipliers =
      OPFLOWSolverGetConstraintMultipliers_TAO;

  PetscFunctionReturn(0);
}
