#include <exago_config.h>
#if defined(EXAGO_ENABLE_IPOPT)

#include "opflow_ipopt.h"
#include <private/opflowimpl.h>

/* IPOPT callback functions */
Bool eval_opflow_f(PetscInt n, PetscScalar *x, Bool new_x,
                   PetscScalar *obj_value, UserDataPtr user_data) {
  PetscErrorCode ierr;
  OPFLOW opflow = (OPFLOW)user_data;

  *obj_value = 0.0;
  ierr = VecPlaceArray(opflow->X, x);
  CHKERRQ(ierr);
  ierr = OPFLOWComputeObjective(opflow, opflow->X, obj_value);
  CHKERRQ(ierr);
  ierr = VecResetArray(opflow->X);
  CHKERRQ(ierr);

  if (opflow->modelops.computeauxobjective) {
    ierr = (*opflow->modelops.computeauxobjective)(opflow, (const double *)x,
                                                   obj_value, opflow->userctx);
    CHKERRQ(ierr);
  }
  return TRUE;
}

Bool eval_opflow_grad_f(PetscInt n, PetscScalar *x, Bool new_x,
                        PetscScalar *grad_f, UserDataPtr user_data) {
  PetscErrorCode ierr;
  OPFLOW opflow = (OPFLOW)user_data;

  ierr = VecPlaceArray(opflow->X, x);
  CHKERRQ(ierr);
  ierr = VecPlaceArray(opflow->gradobj, grad_f);
  CHKERRQ(ierr);
  ierr = VecSet(opflow->gradobj, 0.0);
  CHKERRQ(ierr);
  ierr = OPFLOWComputeGradient(opflow, opflow->X, opflow->gradobj);
  CHKERRQ(ierr);
  ierr = VecResetArray(opflow->X);
  CHKERRQ(ierr);
  ierr = VecResetArray(opflow->gradobj);
  CHKERRQ(ierr);

  if (opflow->modelops.computeauxgradient) {
    ierr = (*opflow->modelops.computeauxgradient)(opflow, (const double *)x,
                                                  grad_f, opflow->userctx);
    CHKERRQ(ierr);
  }

  return TRUE;
}

Bool eval_opflow_g(PetscInt n, PetscScalar *x, Bool new_x, PetscInt m,
                   PetscScalar *g, UserDataPtr user_data) {
  PetscErrorCode ierr;
  OPFLOW opflow = (OPFLOW)user_data;

  ierr = VecPlaceArray(opflow->X, x);
  CHKERRQ(ierr);

  /* Equality constraints */
  ierr = VecPlaceArray(opflow->Ge, g);
  CHKERRQ(ierr);
  ierr = PetscLogEventBegin(opflow->eqconslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = OPFLOWComputeEqualityConstraints(opflow, opflow->X, opflow->Ge);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->eqconslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = VecResetArray(opflow->Ge);
  CHKERRQ(ierr);

  if (opflow->Nconineq) {
    /* Inequality constraints */
    ierr = VecPlaceArray(opflow->Gi, g + opflow->nconeq);
    CHKERRQ(ierr);
    ierr = PetscLogEventBegin(opflow->ineqconslogger, 0, 0, 0, 0);
    CHKERRQ(ierr);
    ierr = OPFLOWComputeInequalityConstraints(opflow, opflow->X, opflow->Gi);
    CHKERRQ(ierr);
    ierr = PetscLogEventEnd(opflow->ineqconslogger, 0, 0, 0, 0);
    CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Gi);
    CHKERRQ(ierr);
  }

  ierr = VecResetArray(opflow->X);
  CHKERRQ(ierr);

  return TRUE;
}

Bool eval_opflow_jac_g(PetscInt n, PetscScalar *x, Bool new_x, PetscInt m,
                       PetscInt nele_jac, PetscInt *iRow, PetscInt *jCol,
                       PetscScalar *values, UserDataPtr user_data) {

  PetscErrorCode ierr;
  OPFLOW opflow = (OPFLOW)user_data;
  PetscInt *iRowstart = iRow, *jColstart = jCol;
  PetscInt roffset, coffset;
  PetscInt nrow, ncol;
  PetscInt nvals;
  const PetscInt *cols;
  const PetscScalar *vals;
  PetscInt i, j;

  if (values == NULL) {
    /* Set locations only */

    roffset = 0;
    coffset = 0;

    /* Equality constrained Jacobian */
    ierr = PetscLogEventBegin(opflow->eqconsjaclogger, 0, 0, 0, 0);
    CHKERRQ(ierr);
    ierr = (*opflow->modelops.computeequalityconstraintjacobian)(
        opflow, opflow->X, opflow->Jac_Ge);
    CHKERRQ(ierr);
    ierr = PetscLogEventEnd(opflow->eqconsjaclogger, 0, 0, 0, 0);
    CHKERRQ(ierr);
    ierr = MatGetSize(opflow->Jac_Ge, &nrow, &ncol);
    CHKERRQ(ierr);

    /* Copy over locations to triplet format */
    for (i = 0; i < nrow; i++) {
      ierr = MatGetRow(opflow->Jac_Ge, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
      for (j = 0; j < nvals; j++) {
        iRowstart[j] = roffset + i;
        jColstart[j] = coffset + cols[j];
      }
      /* Increment iRow,jCol pointers */
      iRowstart += nvals;
      jColstart += nvals;
      ierr = MatRestoreRow(opflow->Jac_Ge, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
    }

    if (opflow->Nconineq) {
      /* Inequality constrained Jacobian */
      roffset = opflow->nconeq;
      ierr = PetscLogEventBegin(opflow->ineqconsjaclogger, 0, 0, 0, 0);
      CHKERRQ(ierr);
      ierr = (*opflow->modelops.computeinequalityconstraintjacobian)(
          opflow, opflow->X, opflow->Jac_Gi);
      CHKERRQ(ierr);
      ierr = PetscLogEventEnd(opflow->ineqconsjaclogger, 0, 0, 0, 0);
      CHKERRQ(ierr);

      ierr = MatGetSize(opflow->Jac_Gi, &nrow, &ncol);
      CHKERRQ(ierr);
      /* Copy over locations to triplet format */
      for (i = 0; i < nrow; i++) {
        ierr = MatGetRow(opflow->Jac_Gi, i, &nvals, &cols, &vals);
        CHKERRQ(ierr);
        for (j = 0; j < nvals; j++) {
          iRowstart[j] = roffset + i;
          jColstart[j] = coffset + cols[j];
        }
        /* Increment iRow,jCol pointers */
        iRowstart += nvals;
        jColstart += nvals;
      }
    }
  } else {
    ierr = VecPlaceArray(opflow->X, x);
    CHKERRQ(ierr);
    /* Compute equality constraint jacobian */
    ierr = PetscLogEventBegin(opflow->eqconsjaclogger, 0, 0, 0, 0);
    CHKERRQ(ierr);
    ierr = (*opflow->modelops.computeequalityconstraintjacobian)(
        opflow, opflow->X, opflow->Jac_Ge);
    CHKERRQ(ierr);
    ierr = PetscLogEventEnd(opflow->eqconsjaclogger, 0, 0, 0, 0);
    CHKERRQ(ierr);

    ierr = MatGetSize(opflow->Jac_Ge, &nrow, &ncol);
    CHKERRQ(ierr);
    /* Copy over values */
    for (i = 0; i < nrow; i++) {
      ierr = MatGetRow(opflow->Jac_Ge, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
      for (j = 0; j < nvals; j++) {
        values[j] = vals[j];
      }
      values += nvals;
      ierr = MatRestoreRow(opflow->Jac_Ge, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
    }

    if (opflow->Nconineq) {
      ierr = PetscLogEventBegin(opflow->ineqconsjaclogger, 0, 0, 0, 0);
      CHKERRQ(ierr);
      /* Compute inequality constraint jacobian */
      ierr = (*opflow->modelops.computeinequalityconstraintjacobian)(
          opflow, opflow->X, opflow->Jac_Gi);
      CHKERRQ(ierr);
      ierr = PetscLogEventEnd(opflow->ineqconsjaclogger, 0, 0, 0, 0);
      CHKERRQ(ierr);

      ierr = MatGetSize(opflow->Jac_Gi, &nrow, &ncol);
      CHKERRQ(ierr);

      /* Copy over values */
      for (i = 0; i < nrow; i++) {
        ierr = MatGetRow(opflow->Jac_Gi, i, &nvals, &cols, &vals);
        CHKERRQ(ierr);
        for (j = 0; j < nvals; j++) {
          values[j] = vals[j];
        }
        values += nvals;
        ierr = MatRestoreRow(opflow->Jac_Gi, i, &nvals, &cols, &vals);
        CHKERRQ(ierr);
      }
    }
    ierr = VecResetArray(opflow->X);
    CHKERRQ(ierr);
  }

  return TRUE;
}

Bool eval_opflow_h(PetscInt n, PetscScalar *x, Bool new_x,
                   PetscScalar obj_factor, PetscInt m, PetscScalar *lambda,
                   Bool new_lambda, PetscInt nele_hess, PetscInt *iRow,
                   PetscInt *jCol, PetscScalar *values, UserDataPtr user_data) {
  PetscErrorCode ierr;
  PetscInt nrow;
  OPFLOW opflow = (OPFLOW)user_data;
  PetscInt nvals;
  const PetscInt *cols;
  const PetscScalar *vals;
  PetscInt i, j;
  PetscInt ctr = 0;

  opflow->obj_factor = obj_factor;

  if (values == NULL) {
    ierr = PetscLogEventBegin(opflow->hesslogger, 0, 0, 0, 0);
    CHKERRQ(ierr);
    ierr = (*opflow->modelops.computehessian)(
        opflow, opflow->X, opflow->Lambdae, opflow->Lambdai, opflow->Hes);
    CHKERRQ(ierr);
    ierr = PetscLogEventEnd(opflow->hesslogger, 0, 0, 0, 0);
    CHKERRQ(ierr);
    ierr = MatGetSize(opflow->Hes, &nrow, &nrow);
    CHKERRQ(ierr);

    /* Copy over locations to triplet format */
    /* Note that IPOPT
       requires a lower diagonal Hessian (see note
       https://coin-or.github.io/Ipopt/INTERFACES.html#INTERFACE_CODE) Hence, we
       only add lower diagonal locations
    */
    for (i = 0; i < nrow; i++) {
      ierr = MatGetRow(opflow->Hes, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
      ctr = 0;
      for (j = 0; j < nvals; j++) {
        if (cols[j] >= i) { /* upper triangle */
          /* save as lower triangle locations */
          iRow[ctr] = cols[j];
          jCol[ctr] = i;
          ctr++;
        }
      }
      iRow += ctr;
      jCol += ctr;
      ierr = MatRestoreRow(opflow->Hes, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
    }

  } else {
    ierr = VecPlaceArray(opflow->X, x);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Lambdae, lambda);
    CHKERRQ(ierr);
    if (opflow->Nconineq) {
      ierr = VecPlaceArray(opflow->Lambdai, lambda + opflow->nconeq);
      CHKERRQ(ierr);
    }

    ierr = PetscLogEventBegin(opflow->hesslogger, 0, 0, 0, 0);
    CHKERRQ(ierr);
    /* Compute Hessian */
    ierr = (*opflow->modelops.computehessian)(
        opflow, opflow->X, opflow->Lambdae, opflow->Lambdai, opflow->Hes);
    CHKERRQ(ierr);
    ierr = PetscLogEventEnd(opflow->hesslogger, 0, 0, 0, 0);
    CHKERRQ(ierr);

    /* Auxillary hessian */
    if (opflow->modelops.computeauxhessian) {
      ierr = (*opflow->modelops.computeauxhessian)(opflow, x, opflow->Hes,
                                                   opflow->userctx);
      CHKERRQ(ierr);
      ierr = MatAssemblyBegin(opflow->Hes, MAT_FINAL_ASSEMBLY);
      CHKERRQ(ierr);
      ierr = MatAssemblyEnd(opflow->Hes, MAT_FINAL_ASSEMBLY);
      CHKERRQ(ierr);
    }

    ierr = VecResetArray(opflow->X);
    CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Lambdae);
    CHKERRQ(ierr);
    if (opflow->Nconineq) {
      ierr = VecResetArray(opflow->Lambdai);
      CHKERRQ(ierr);
    }

    /* Copy over values */
    ierr = MatGetSize(opflow->Hes, &nrow, &nrow);
    CHKERRQ(ierr);
    for (i = 0; i < nrow; i++) {
      ierr = MatGetRow(opflow->Hes, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
      ctr = 0;
      for (j = 0; j < nvals; j++) {
        if (cols[j] >= i) { /* Upper triangle values (same as lower triangle) */
          values[ctr] = vals[j];
          ctr++;
        }
      }
      values += ctr;
      ierr = MatRestoreRow(opflow->Hes, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
    }
  }

  return 1;
}

Bool OPFLOWSolverMonitor_IPOPT(Index alg_mod, Index iter_count,
                               Number obj_value, Number inf_pr, Number inf_du,
                               Number mu, Number d_norm,
                               Number regularization_size, Number alpha_du,
                               Number alpha_pr, Index ls_trials,
                               UserDataPtr user_data) {
  OPFLOW opflow = (OPFLOW)user_data;
  opflow->numits = iter_count;
  return 1;
}

PetscErrorCode OPFLOWSolverSolve_IPOPT(OPFLOW opflow) {
  PetscErrorCode ierr;
  OPFLOWSolver_IPOPT ipopt = (OPFLOWSolver_IPOPT)opflow->solver;
  PetscScalar *x, *g, *lambda;

  PetscFunctionBegin;

  ierr = VecGetArray(opflow->X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArray(opflow->G, &g);
  CHKERRQ(ierr);
  ierr = VecGetArray(opflow->Lambda, &lambda);
  CHKERRQ(ierr);
  /* Solve */
  ierr = PetscLogEventBegin(opflow->solvelogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ipopt->solve_status =
      IpoptSolve(ipopt->nlp, x, g, &opflow->obj, lambda, NULL, NULL, opflow);
  ierr = PetscLogEventEnd(opflow->solvelogger, 0, 0, 0, 0);
  CHKERRQ(ierr);

  ierr = VecRestoreArray(opflow->X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->G, &g);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->Lambda, &lambda);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverDestroy_IPOPT(OPFLOW opflow) {
  PetscErrorCode ierr;
  OPFLOWSolver_IPOPT ipopt = (OPFLOWSolver_IPOPT)opflow->solver;

  PetscFunctionBegin;

  if (ipopt->nlp) {
    FreeIpoptProblem(ipopt->nlp);
  }

  ierr = PetscFree(ipopt);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetObjective_IPOPT(OPFLOW opflow, PetscReal *obj) {
  PetscFunctionBegin;
  *obj = opflow->obj;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetSolution_IPOPT(OPFLOW opflow, Vec *X) {
  PetscFunctionBegin;
  *X = opflow->X;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConstraints_IPOPT(OPFLOW opflow, Vec *G) {
  PetscFunctionBegin;
  *G = opflow->G;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConstraintMultipliers_IPOPT(OPFLOW opflow,
                                                          Vec *Lambda) {
  PetscFunctionBegin;
  *Lambda = opflow->Lambda;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConvergenceStatus_IPOPT(OPFLOW opflow,
                                                      PetscBool *status) {
  OPFLOWSolver_IPOPT ipopt = (OPFLOWSolver_IPOPT)opflow->solver;

  PetscFunctionBegin;
  if (ipopt->solve_status == 0 || ipopt->solve_status == 1)
    *status = PETSC_TRUE; /* See IpReturnCodes_inc.h in IPOPT. The first two
                             denote convergence */
  else
    *status = PETSC_FALSE;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverSetUp_IPOPT(OPFLOW opflow) {
  PetscErrorCode ierr;
  OPFLOWSolver_IPOPT ipopt = (OPFLOWSolver_IPOPT)opflow->solver;
  PetscScalar *xl, *xu, *gl, *gu;
  MatInfo info_eq, info_ineq, info_hes;

  PetscFunctionBegin;

  /* Compute nonzeros for the Jacobian */
  /* Equality constraint Jacobian */
  ierr = PetscLogEventBegin(opflow->eqconsjaclogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeequalityconstraintjacobian)(
      opflow, opflow->X, opflow->Jac_Ge);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->eqconsjaclogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = MatSetOption(opflow->Jac_Ge, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  CHKERRQ(ierr);

  ierr = MatGetInfo(opflow->Jac_Ge, MAT_LOCAL, &info_eq);
  CHKERRQ(ierr);
  ipopt->nnz_jac_ge = info_eq.nz_used;

  ipopt->nnz_jac_gi = 0;
  if (opflow->Nconineq) {
    ierr = PetscLogEventBegin(opflow->ineqconsjaclogger, 0, 0, 0, 0);
    CHKERRQ(ierr);
    ierr = (*opflow->modelops.computeinequalityconstraintjacobian)(
        opflow, opflow->X, opflow->Jac_Gi);
    CHKERRQ(ierr);
    ierr = PetscLogEventEnd(opflow->ineqconsjaclogger, 0, 0, 0, 0);
    CHKERRQ(ierr);
    ierr =
        MatSetOption(opflow->Jac_Gi, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
    CHKERRQ(ierr);

    ierr = MatGetInfo(opflow->Jac_Gi, MAT_LOCAL, &info_ineq);
    CHKERRQ(ierr);
    ipopt->nnz_jac_gi = info_ineq.nz_used;
  }
  ipopt->nnz_jac_g = ipopt->nnz_jac_ge + ipopt->nnz_jac_gi;

  /* Compute non-zeros for Hessian */
  ierr = PetscLogEventBegin(opflow->hesslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computehessian)(opflow, opflow->X, opflow->Lambdae,
                                            opflow->Lambdai, opflow->Hes);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->hesslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = MatSetOption(opflow->Hes, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  CHKERRQ(ierr);

  ierr = MatGetInfo(opflow->Hes, MAT_LOCAL, &info_hes);
  CHKERRQ(ierr);
  ipopt->nnz_hes = (info_hes.nz_used - opflow->nx) / 2 + opflow->nx;

  /* Create IPOPT solver instance */
  ierr = VecGetArray(opflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecGetArray(opflow->Xu, &xu);
  CHKERRQ(ierr);
  ierr = VecGetArray(opflow->Gl, &gl);
  CHKERRQ(ierr);
  ierr = VecGetArray(opflow->Gu, &gu);
  CHKERRQ(ierr);

  /* Create IPOPT problem */
  ipopt->nlp = CreateIpoptProblem(
      opflow->nx, xl, xu, opflow->ncon, gl, gu, ipopt->nnz_jac_g,
      ipopt->nnz_hes, 0, &eval_opflow_f, &eval_opflow_g, &eval_opflow_grad_f,
      &eval_opflow_jac_g, &eval_opflow_h);

  ierr = VecRestoreArray(opflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->Xu, &xu);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->Gl, &gl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->Gu, &gu);
  CHKERRQ(ierr);

  /* IPOPT solver options */
  {
    AddIpoptNumOption(ipopt->nlp, (char *)"tol", opflow->tolerance);
    AddIpoptIntOption(ipopt->nlp, (char *)"max_iter", 5000);
    AddIpoptStrOption(ipopt->nlp, (char *)"mu_strategy", (char *)"monotone");
    AddIpoptStrOption(ipopt->nlp, (char *)"fixed_variable_treatment",
                      (char *)"relax_bounds");
    AddIpoptStrOption(ipopt->nlp, (char *)"inf_pr_output", (char *)"internal");
    AddIpoptNumOption(ipopt->nlp, (char *)"constr_mult_init_max", 0.0);
    AddIpoptNumOption(ipopt->nlp, (char *)"residual_ratio_max", 1e3);
    AddIpoptNumOption(ipopt->nlp, (char *)"residual_ratio_singular", 1e4);
  }

  /* Add intermediate callback function to get the solver info
     Note this is called by IPOPT at each iteration
  */
  SetIntermediateCallback(ipopt->nlp, OPFLOWSolverMonitor_IPOPT);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverCreate_IPOPT(OPFLOW opflow) {
  PetscErrorCode ierr;
  OPFLOWSolver_IPOPT ipopt;

  PetscFunctionBegin;

  if (opflow->comm->size > 1)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,
            "IPOPT solver does not support execution in parallel\n");
  ierr = PetscCalloc1(1, &ipopt);
  CHKERRQ(ierr);

  ipopt->nlp = NULL;
  ipopt->nnz_jac_g = 0;
  ipopt->nnz_hes = 0;
  opflow->solver = ipopt;

  opflow->solverops.setup = OPFLOWSolverSetUp_IPOPT;
  opflow->solverops.solve = OPFLOWSolverSolve_IPOPT;
  opflow->solverops.destroy = OPFLOWSolverDestroy_IPOPT;
  opflow->solverops.getobjective = OPFLOWSolverGetObjective_IPOPT;
  opflow->solverops.getsolution = OPFLOWSolverGetSolution_IPOPT;
  opflow->solverops.getconvergencestatus =
      OPFLOWSolverGetConvergenceStatus_IPOPT;
  opflow->solverops.getconstraints = OPFLOWSolverGetConstraints_IPOPT;
  opflow->solverops.getconstraintmultipliers =
      OPFLOWSolverGetConstraintMultipliers_IPOPT;

  PetscFunctionReturn(0);
}

#endif
