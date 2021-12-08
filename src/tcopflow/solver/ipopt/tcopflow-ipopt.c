#include <exago_config.h>
#if defined(EXAGO_ENABLE_IPOPT)

#include "tcopflow-ipopt.h"
#include <private/opflowimpl.h>
#include <private/tcopflowimpl.h>

/* IPOPT callback functions */
Bool eval_tcopflow_f(PetscInt n, PetscScalar *x, Bool new_x,
                     PetscScalar *obj_value, UserDataPtr user_data) {
  PetscErrorCode ierr;
  TCOPFLOW tcopflow = (TCOPFLOW)user_data;

  *obj_value = 0.0;

  ierr = VecPlaceArray(tcopflow->X, x);
  CHKERRQ(ierr);
  ierr =
      (*tcopflow->modelops.computeobjective)(tcopflow, tcopflow->X, obj_value);
  CHKERRQ(ierr);
  ierr = VecResetArray(tcopflow->X);
  CHKERRQ(ierr);

  return TRUE;
}

Bool eval_tcopflow_grad_f(PetscInt n, PetscScalar *x, Bool new_x,
                          PetscScalar *grad_f, UserDataPtr user_data) {

  PetscErrorCode ierr;
  TCOPFLOW tcopflow = (TCOPFLOW)user_data;

  ierr = VecPlaceArray(tcopflow->X, x);
  CHKERRQ(ierr);
  ierr = VecPlaceArray(tcopflow->gradobj, grad_f);
  CHKERRQ(ierr);
  ierr = (*tcopflow->modelops.computegradient)(tcopflow, tcopflow->X,
                                               tcopflow->gradobj);
  CHKERRQ(ierr);
  ierr = VecResetArray(tcopflow->X);
  CHKERRQ(ierr);
  ierr = VecResetArray(tcopflow->gradobj);
  CHKERRQ(ierr);

  return TRUE;
}

Bool eval_tcopflow_g(PetscInt n, PetscScalar *x, Bool new_x, PetscInt m,
                     PetscScalar *g, UserDataPtr user_data) {
  PetscErrorCode ierr;
  TCOPFLOW tcopflow = (TCOPFLOW)user_data;

  ierr = VecPlaceArray(tcopflow->X, x);
  CHKERRQ(ierr);
  ierr = VecPlaceArray(tcopflow->G, g);
  CHKERRQ(ierr);
  ierr = (*tcopflow->modelops.computeconstraints)(tcopflow, tcopflow->X,
                                                  tcopflow->G);
  CHKERRQ(ierr);
  ierr = VecResetArray(tcopflow->X);
  CHKERRQ(ierr);
  ierr = VecResetArray(tcopflow->G);
  CHKERRQ(ierr);

  return TRUE;
}

Bool eval_tcopflow_jac_g(PetscInt n, PetscScalar *x, Bool new_x, PetscInt m,
                         PetscInt nele_jac, PetscInt *iRow, PetscInt *jCol,
                         PetscScalar *values, UserDataPtr user_data) {

  PetscErrorCode ierr;
  TCOPFLOW tcopflow = (TCOPFLOW)user_data;
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

    /* Compute Jacobian */
    ierr = (*tcopflow->modelops.computejacobian)(tcopflow, tcopflow->X,
                                                 tcopflow->Jac);
    CHKERRQ(ierr);

    ierr = MatGetSize(tcopflow->Jac, &nrow, &ncol);
    CHKERRQ(ierr);

    /* Copy over locations to triplet format */
    for (i = 0; i < nrow; i++) {
      ierr = MatGetRow(tcopflow->Jac, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
      for (j = 0; j < nvals; j++) {
        iRowstart[j] = roffset + i;
        jColstart[j] = coffset + cols[j];
      }
      /* Increment iRow,jCol pointers */
      iRowstart += nvals;
      jColstart += nvals;
      ierr = MatRestoreRow(tcopflow->Jac, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
    }
  } else {
    ierr = VecPlaceArray(tcopflow->X, x);
    CHKERRQ(ierr);
    /* Compute jacobian */
    ierr = (*tcopflow->modelops.computejacobian)(tcopflow, tcopflow->X,
                                                 tcopflow->Jac);
    CHKERRQ(ierr);

    ierr = MatGetSize(tcopflow->Jac, &nrow, &ncol);
    CHKERRQ(ierr);
    /* Copy over values */
    for (i = 0; i < nrow; i++) {
      ierr = MatGetRow(tcopflow->Jac, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
      for (j = 0; j < nvals; j++) {
        values[j] = vals[j];
      }
      values += nvals;
      ierr = MatRestoreRow(tcopflow->Jac, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
    }
    ierr = VecResetArray(tcopflow->X);
    CHKERRQ(ierr);
  }

  return TRUE;
}

Bool eval_tcopflow_h(PetscInt n, PetscScalar *x, Bool new_x,
                     PetscScalar obj_factor, PetscInt m, PetscScalar *lambda,
                     Bool new_lambda, PetscInt nele_hess, PetscInt *iRow,
                     PetscInt *jCol, PetscScalar *values,
                     UserDataPtr user_data) {
  PetscErrorCode ierr;
  PetscInt nrow;
  TCOPFLOW tcopflow = (TCOPFLOW)user_data;
  PetscInt j, k;
  const PetscInt *cols;
  const PetscScalar *vals;
  PetscInt ctr = 0;
  PetscInt nvals;

  tcopflow->obj_factor = obj_factor;

  ierr = MatGetSize(tcopflow->Hes, &nrow, NULL);
  CHKERRQ(ierr);

  if (values == NULL) {
    /* Set locations only */

    /* Compute Hessian */
    ierr = (*tcopflow->modelops.computehessian)(
        tcopflow, tcopflow->X, tcopflow->Lambda, tcopflow->Hes);
    CHKERRQ(ierr);

    /* Copy over locations to triplet format */
    /* Note that IPOPT
       requires a lower diagonal Hessian (see note
       https://coin-or.github.io/Ipopt/INTERFACES.html#INTERFACE_CODE) Hence, we
       only add lower diagonal locations
    */
    for (j = 0; j < nrow; j++) {
      ierr = MatGetRow(tcopflow->Hes, j, &nvals, &cols, &vals);
      CHKERRQ(ierr);
      ctr = 0;
      for (k = 0; k < nvals; k++) {
        if (cols[k] >= j) { /* upper triangle */
                            /* save as lower triangle locations */
          iRow[ctr] = cols[k];
          jCol[ctr] = j;
          ctr++;
        }
      }
      iRow += ctr;
      jCol += ctr;
      ierr = MatRestoreRow(tcopflow->Hes, j, &nvals, &cols, &vals);
      CHKERRQ(ierr);
    }
  } else {
    /* Copy values */
    ierr = VecPlaceArray(tcopflow->X, x);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(tcopflow->Lambda, lambda);
    CHKERRQ(ierr);

    /* Compute Hessian */
    ierr = (*tcopflow->modelops.computehessian)(
        tcopflow, tcopflow->X, tcopflow->Lambda, tcopflow->Hes);
    CHKERRQ(ierr);

    /* Copy over values to triplet format */
    /* Note that IPOPT
       requires a lower diagonal Hessian (see note
       https://coin-or.github.io/Ipopt/INTERFACES.html#INTERFACE_CODE) Hence, we
       only add lower diagonal locations
    */
    for (j = 0; j < nrow; j++) {
      ierr = MatGetRow(tcopflow->Hes, j, &nvals, &cols, &vals);
      CHKERRQ(ierr);
      ctr = 0;
      for (k = 0; k < nvals; k++) {
        if (cols[k] >= j) { /* upper triangle */
          /* save as lower triangle locations */
          values[ctr] = vals[k];
          ctr++;
        }
      }
      values += ctr;
      ierr = MatRestoreRow(tcopflow->Hes, j, &nvals, &cols, &vals);
      CHKERRQ(ierr);
    }

    ierr = VecResetArray(tcopflow->X);
    CHKERRQ(ierr);
    ierr = VecResetArray(tcopflow->Lambda);
    CHKERRQ(ierr);
  }

  return TRUE;
}

Bool TCOPFLOWSolverMonitor_IPOPT(Index alg_mod, Index iter_count,
                                 Number obj_value, Number inf_pr, Number inf_du,
                                 Number mu, Number d_norm,
                                 Number regularization_size, Number alpha_du,
                                 Number alpha_pr, Index ls_trials,
                                 UserDataPtr user_data) {
  TCOPFLOW tcopflow = (TCOPFLOW)user_data;
  tcopflow->numiter = iter_count;
  return 1;
}

PetscErrorCode TCOPFLOWSolverSolve_IPOPT(TCOPFLOW tcopflow) {
  PetscErrorCode ierr;
  TCOPFLOWSolver_IPOPT tcopflowipopt = (TCOPFLOWSolver_IPOPT)tcopflow->solver;
  MatInfo info_jac, info_hes;
  PetscScalar *x, *g, *xl, *xu, *gl, *gu, *lam;

  PetscFunctionBegin;

  tcopflowipopt->nnz_jac_g = tcopflowipopt->nnz_hes = 0;

  /* Compute nonzeros for the Jacobian */
  ierr = (*tcopflow->modelops.computejacobian)(tcopflow, tcopflow->X,
                                               tcopflow->Jac);
  CHKERRQ(ierr);
  ierr = MatSetOption(tcopflow->Jac, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  CHKERRQ(ierr);

  ierr = MatGetInfo(tcopflow->Jac, MAT_LOCAL, &info_jac);
  CHKERRQ(ierr);
  tcopflowipopt->nnz_jac_g = info_jac.nz_used;

  /* Compute non-zeros for Hessian */
  ierr = (*tcopflow->modelops.computehessian)(tcopflow, tcopflow->X,
                                              tcopflow->Lambda, tcopflow->Hes);
  CHKERRQ(ierr);
  ierr = MatSetOption(tcopflow->Hes, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  CHKERRQ(ierr);

  ierr = MatGetInfo(tcopflow->Hes, MAT_LOCAL, &info_hes);
  CHKERRQ(ierr);
  tcopflowipopt->nnz_hes = (info_hes.nz_used - tcopflow->Nx) / 2 + tcopflow->Nx;

  ierr = VecGetArray(tcopflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecGetArray(tcopflow->Xu, &xu);
  CHKERRQ(ierr);
  ierr = VecGetArray(tcopflow->Gl, &gl);
  CHKERRQ(ierr);
  ierr = VecGetArray(tcopflow->Gu, &gu);
  CHKERRQ(ierr);

  /* Create IPOPT problem */
  tcopflowipopt->nlp = CreateIpoptProblem(
      tcopflow->Nx, xl, xu, tcopflow->Ncon, gl, gu, tcopflowipopt->nnz_jac_g,
      tcopflowipopt->nnz_hes, 0, &eval_tcopflow_f, &eval_tcopflow_g,
      &eval_tcopflow_grad_f, &eval_tcopflow_jac_g, &eval_tcopflow_h);

  ierr = VecRestoreArray(tcopflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(tcopflow->Xu, &xu);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(tcopflow->Gl, &gl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(tcopflow->Gu, &gu);
  CHKERRQ(ierr);

  /* IPOPT solver options */
  AddIpoptNumOption(tcopflowipopt->nlp, (char *)"tol", tcopflow->tolerance);
  AddIpoptIntOption(tcopflowipopt->nlp, (char *)"max_iter", 5000);
  AddIpoptStrOption(tcopflowipopt->nlp, (char *)"mu_strategy",
                    (char *)"monotone");
  AddIpoptStrOption(tcopflowipopt->nlp, (char *)"fixed_variable_treatment",
                    (char *)"relax_bounds");
  AddIpoptStrOption(tcopflowipopt->nlp, (char *)"inf_pr_output",
                    (char *)"internal");
  AddIpoptNumOption(tcopflowipopt->nlp, (char *)"constr_mult_init_max", 0.0);
  AddIpoptNumOption(tcopflowipopt->nlp, (char *)"residual_ratio_max", 1e3);
  AddIpoptNumOption(tcopflowipopt->nlp, (char *)"residual_ratio_singular", 1e4);

  ierr = VecGetArray(tcopflow->X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArray(tcopflow->G, &g);
  CHKERRQ(ierr);
  ierr = VecGetArray(tcopflow->Lambda, &lam);
  CHKERRQ(ierr);

  /* Add intermediate callback to get solver info
     This is called by IPOPT each iteration
  */
  SetIntermediateCallback(tcopflowipopt->nlp, TCOPFLOWSolverMonitor_IPOPT);

  /* Solve */
  tcopflowipopt->solve_status = IpoptSolve(
      tcopflowipopt->nlp, x, g, &tcopflow->obj, lam, NULL, NULL, tcopflow);

  ierr = VecRestoreArray(tcopflow->X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(tcopflow->G, &g);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(tcopflow->Lambda, &lam);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode TCOPFLOWSolverDestroy_IPOPT(TCOPFLOW tcopflow) {
  PetscErrorCode ierr;
  TCOPFLOWSolver_IPOPT ipopt = (TCOPFLOWSolver_IPOPT)tcopflow->solver;

  PetscFunctionBegin;

  if (ipopt->nlp) {
    FreeIpoptProblem(ipopt->nlp);
    ipopt->nlp = NULL;
  }

  ierr = PetscFree(ipopt);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TCOPFLOWSolverSetUp_IPOPT(TCOPFLOW tcopflow) {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

PetscErrorCode TCOPFLOWSolverGetObjective_IPOPT(TCOPFLOW tcopflow,
                                                PetscReal *obj) {
  PetscFunctionBegin;
  *obj = tcopflow->obj;
  PetscFunctionReturn(0);
}

PetscErrorCode TCOPFLOWSolverGetSolution_IPOPT(TCOPFLOW tcopflow,
                                               PetscInt t_num, Vec *X) {
  PetscErrorCode ierr;
  OPFLOW opflow = tcopflow->opflows[t_num];
  Vec Xi = opflow->X;
  PetscInt nxi = opflow->nx;
  PetscScalar *xi, *x;
  PetscInt ix = tcopflow->xstarti[t_num];

  PetscFunctionBegin;
  ierr = VecGetArray(Xi, &xi);
  CHKERRQ(ierr);
  ierr = VecGetArray(tcopflow->X, &x);
  CHKERRQ(ierr);

  ierr = PetscArraycpy(xi, x + ix, nxi);
  CHKERRQ(ierr);

  ierr = VecRestoreArray(Xi, &xi);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(tcopflow->X, &x);
  CHKERRQ(ierr);

  *X = Xi;
  PetscFunctionReturn(0);
}

PetscErrorCode TCOPFLOWSolverGetConstraints_IPOPT(TCOPFLOW tcopflow,
                                                  PetscInt t_num, Vec *G) {
  PetscErrorCode ierr;
  OPFLOW opflow = tcopflow->opflows[t_num];
  Vec Gi = opflow->G;
  PetscInt ngi = opflow->ncon;
  PetscScalar *gi, *g;
  PetscInt ig = tcopflow->gstarti[t_num];

  PetscFunctionBegin;
  ierr = VecGetArray(Gi, &gi);
  CHKERRQ(ierr);
  ierr = VecGetArray(tcopflow->G, &g);
  CHKERRQ(ierr);

  ierr = PetscArraycpy(gi, g + ig, ngi);
  CHKERRQ(ierr);

  ierr = VecRestoreArray(Gi, &gi);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(tcopflow->G, &g);
  CHKERRQ(ierr);

  *G = Gi;

  PetscFunctionReturn(0);
}

PetscErrorCode TCOPFLOWSolverGetConstraintMultipliers_IPOPT(TCOPFLOW tcopflow,
                                                            PetscInt t_num,
                                                            Vec *Lambda) {
  PetscErrorCode ierr;
  OPFLOW opflow = tcopflow->opflows[t_num];
  Vec Lambdai = opflow->Lambda;
  PetscInt ngi = opflow->ncon;
  PetscScalar *lambdai, *lambda;
  PetscInt ig = tcopflow->gstarti[t_num];

  PetscFunctionBegin;
  ierr = VecGetArray(Lambdai, &lambdai);
  CHKERRQ(ierr);
  ierr = VecGetArray(tcopflow->Lambda, &lambda);
  CHKERRQ(ierr);

  ierr = PetscArraycpy(lambdai, lambda + ig, ngi);
  CHKERRQ(ierr);

  ierr = VecRestoreArray(Lambdai, &lambdai);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(tcopflow->Lambda, &lambda);
  CHKERRQ(ierr);

  *Lambda = Lambdai;

  PetscFunctionReturn(0);
}

PetscErrorCode TCOPFLOWSolverGetConvergenceStatus_IPOPT(TCOPFLOW tcopflow,
                                                        PetscBool *status) {
  TCOPFLOWSolver_IPOPT ipopt = (TCOPFLOWSolver_IPOPT)tcopflow->solver;

  PetscFunctionBegin;
  if (ipopt->solve_status < 2)
    *status = PETSC_TRUE; /* See IpReturnCodes_inc.h in IPOPT. The first two
                             denote convergence */
  else
    *status = PETSC_FALSE;

  PetscFunctionReturn(0);
}

PetscErrorCode TCOPFLOWSolverCreate_IPOPT(TCOPFLOW tcopflow) {
  PetscErrorCode ierr;
  TCOPFLOWSolver_IPOPT ipopt;

  PetscFunctionBegin;

  if (tcopflow->comm->size > 1)
    SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_SUP,
             "IPOPT solver does not support execution in parallel\n",
             tcopflow->comm->size);
  ierr = PetscCalloc1(1, &ipopt);
  CHKERRQ(ierr);

  ipopt->nlp = NULL;
  ipopt->nnz_jac_g = 0;
  ipopt->nnz_hes = 0;
  tcopflow->solver = ipopt;

  tcopflow->solverops.setup = TCOPFLOWSolverSetUp_IPOPT;
  tcopflow->solverops.solve = TCOPFLOWSolverSolve_IPOPT;
  tcopflow->solverops.destroy = TCOPFLOWSolverDestroy_IPOPT;
  tcopflow->solverops.getobjective = TCOPFLOWSolverGetObjective_IPOPT;
  tcopflow->solverops.getsolution = TCOPFLOWSolverGetSolution_IPOPT;
  tcopflow->solverops.getconvergencestatus =
      TCOPFLOWSolverGetConvergenceStatus_IPOPT;
  tcopflow->solverops.getconstraints = TCOPFLOWSolverGetConstraints_IPOPT;
  tcopflow->solverops.getconstraintmultipliers =
      TCOPFLOWSolverGetConstraintMultipliers_IPOPT;

  PetscFunctionReturn(0);
}

#endif
