#include <exago_config.h>
#if defined(EXAGO_ENABLE_HIOP)
#if defined(EXAGO_ENABLE_HIOP_SPARSE)

#include "opflow_hiopsparse.h"
#include <private/opflowimpl.h>

OPFLOWHIOPSPARSEInterface::OPFLOWHIOPSPARSEInterface(OPFLOW opflowin) {
  opflow = opflowin;
}

bool OPFLOWHIOPSPARSEInterface::get_prob_sizes(hiop::size_type &n,
                                               hiop::size_type &m) {
  n = opflow->nx;
  m = opflow->ncon;
  return true;
}

bool OPFLOWHIOPSPARSEInterface::get_vars_info(const hiop::size_type &n,
                                              double *xlow, double *xupp,
                                              NonlinearityType *type) {
  PetscErrorCode ierr;
  PetscInt i;
  const PetscScalar *xl, *xu;

  ierr = (*opflow->modelops.setvariablebounds)(opflow, opflow->Xl, opflow->Xu);
  CHKERRQ(ierr);
  if (opflow->modelops.updatevariablebounds) {
    ierr = (*opflow->modelops.updatevariablebounds)(
        opflow, opflow->Xl, opflow->Xu, opflow->userctx);
    CHKERRQ(ierr);
  }

  ierr = VecGetArrayRead(opflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecGetArrayRead(opflow->Xu, &xu);
  CHKERRQ(ierr);

  for (i = 0; i < n; i++) {
    xlow[i] = xl[i];
    xupp[i] = xu[i];
    type[i] = hiopNonlinear;
  }

  ierr = VecRestoreArrayRead(opflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(opflow->Xu, &xu);
  CHKERRQ(ierr);

  return true;
}

bool OPFLOWHIOPSPARSEInterface::get_cons_info(const hiop::size_type &m,
                                              double *clow, double *cupp,
                                              NonlinearityType *type) {
  PetscInt i;
  PetscErrorCode ierr;

  ierr = VecPlaceArray(opflow->Gl, clow);
  CHKERRQ(ierr);
  ierr = VecPlaceArray(opflow->Gu, cupp);
  CHKERRQ(ierr);

  ierr =
      (*opflow->modelops.setconstraintbounds)(opflow, opflow->Gl, opflow->Gu);
  CHKERRQ(ierr);

  ierr = VecResetArray(opflow->Gl);
  CHKERRQ(ierr);
  ierr = VecResetArray(opflow->Gu);
  CHKERRQ(ierr);

  for (i = 0; i < m; i++)
    type[i] = hiopNonlinear;

  return true;
}

bool OPFLOWHIOPSPARSEInterface::get_sparse_blocks_info(
    hiop::size_type &nx, hiop::size_type &nnz_sparse_Jaceq,
    hiop::size_type &nnz_sparse_Jacineq,
    hiop::size_type &nnz_sparse_Hess_Lagr) {
  PetscErrorCode ierr;
  PetscScalar *xl, *xu, *gl, *gu;
  MatInfo info_eq, info_ineq, info_hes;

  nx = opflow->nx;

  /* Compute nonzeros for the Jacobian */
  /* Equality constraint Jacobian */
  ierr = (*opflow->modelops.computeequalityconstraintjacobian)(
      opflow, opflow->X, opflow->Jac_Ge);
  CHKERRQ(ierr);
  ierr = MatSetOption(opflow->Jac_Ge, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  CHKERRQ(ierr);

  ierr = MatGetInfo(opflow->Jac_Ge, MAT_LOCAL, &info_eq);
  CHKERRQ(ierr);

  nnz_sparse_Jaceq = info_eq.nz_used;

  nnz_sparse_Jacineq = 0;
  if (opflow->Nconineq) {
    ierr = (*opflow->modelops.computeinequalityconstraintjacobian)(
        opflow, opflow->X, opflow->Jac_Gi);
    CHKERRQ(ierr);
    ierr =
        MatSetOption(opflow->Jac_Gi, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
    CHKERRQ(ierr);

    ierr = MatGetInfo(opflow->Jac_Gi, MAT_LOCAL, &info_ineq);
    CHKERRQ(ierr);

    nnz_sparse_Jacineq = info_ineq.nz_used;
  }

  /* Compute non-zeros for Hessian */
  ierr = (*opflow->modelops.computehessian)(opflow, opflow->X, opflow->Lambdae,
                                            opflow->Lambdai, opflow->Hes);
  CHKERRQ(ierr);
  ierr = MatSetOption(opflow->Hes, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  CHKERRQ(ierr);

  ierr = MatGetInfo(opflow->Hes, MAT_LOCAL, &info_hes);
  CHKERRQ(ierr);

  nnz_sparse_Hess_Lagr = (info_hes.nz_used - opflow->nx) / 2 + opflow->nx;

  return true;
}

bool OPFLOWHIOPSPARSEInterface::eval_f(const hiop::size_type &n,
                                       const double *x, bool new_x,
                                       double &obj_value) {
  PetscErrorCode ierr;
  PetscScalar *xarr;

  obj_value = 0.0;
  ierr = VecPlaceArray(opflow->X, x);
  CHKERRQ(ierr);

  /* Compute objective */
  ierr = (*opflow->modelops.computeobjective)(opflow, opflow->X, &obj_value);
  CHKERRQ(ierr);

  ierr = VecResetArray(opflow->X);

  if (opflow->modelops.computeauxobjective) {
    ierr = (*opflow->modelops.computeauxobjective)(opflow, (const double *)x,
                                                   &obj_value, opflow->userctx);
    CHKERRQ(ierr);
  }

  CHKERRQ(ierr);

  return true;
}

bool OPFLOWHIOPSPARSEInterface::eval_cons(const hiop::size_type &n,
                                          const hiop::size_type &m,
                                          const hiop::size_type &num_cons,
                                          const hiop::size_type *idx_cons,
                                          const double *x, bool new_x,
                                          double *cons) {
  return false;
}

bool OPFLOWHIOPSPARSEInterface::eval_cons(const hiop::size_type &n,
                                          const hiop::size_type &m,
                                          const double *x, bool new_x,
                                          double *cons) {
  PetscErrorCode ierr;

  ierr = VecPlaceArray(opflow->X, x);
  CHKERRQ(ierr);

  /* Equality constaints */
  ierr = VecPlaceArray(opflow->Ge, cons);
  CHKERRQ(ierr);
  ierr = VecSet(opflow->Ge, 0.0);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeequalityconstraints)(opflow, opflow->X,
                                                        opflow->Ge);
  CHKERRQ(ierr);
  ierr = VecResetArray(opflow->Ge);
  CHKERRQ(ierr);

  if (opflow->nconineq) {
    /* Inequality constraints */
    ierr = VecPlaceArray(opflow->Gi, cons + opflow->nconeq);
    CHKERRQ(ierr);
    ierr = VecSet(opflow->Gi, 0.0);
    CHKERRQ(ierr);
    ierr = (*opflow->modelops.computeinequalityconstraints)(opflow, opflow->X,
                                                            opflow->Gi);
    CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Gi);
    CHKERRQ(ierr);
  }

  ierr = VecResetArray(opflow->X);
  CHKERRQ(ierr);

  return true;
}

bool OPFLOWHIOPSPARSEInterface::eval_grad_f(const hiop::size_type &n,
                                            const double *x, bool new_x,
                                            double *gradf) {
  PetscErrorCode ierr;

  ierr = VecPlaceArray(opflow->X, x);
  CHKERRQ(ierr);
  ierr = VecPlaceArray(opflow->gradobj, gradf);
  CHKERRQ(ierr);

  ierr = VecSet(opflow->gradobj, 0.0);
  CHKERRQ(ierr);
  ierr =
      (*opflow->modelops.computegradient)(opflow, opflow->X, opflow->gradobj);
  CHKERRQ(ierr);

  ierr = VecResetArray(opflow->X);
  CHKERRQ(ierr);
  ierr = VecResetArray(opflow->gradobj);
  CHKERRQ(ierr);

  if (opflow->modelops.computeauxgradient) {
    ierr = (*opflow->modelops.computeauxgradient)(opflow, (const double *)x,
                                                  gradf, opflow->userctx);
    CHKERRQ(ierr);
  }

  return true;
}

bool OPFLOWHIOPSPARSEInterface::eval_Jac_cons(
    const hiop::size_type &n, const hiop::size_type &m,
    const hiop::size_type &num_cons, const hiop::size_type *idx_cons,
    const double *x, bool new_x, const hiop::size_type &nnzJacS,
    hiop::index_type *iJacS, hiop::index_type *jJacS, double *MJacS) {
  return false;
}

bool OPFLOWHIOPSPARSEInterface::eval_Jac_cons(
    const hiop::size_type &n, const hiop::size_type &m, const double *x,
    bool new_x, const hiop::size_type &nnzJacS, hiop::index_type *iRow,
    hiop::index_type *jCol, double *values) {
  PetscErrorCode ierr;
  PetscInt *iRowstart = iRow, *jColstart = jCol;
  PetscInt roffset, coffset;
  PetscInt nrow, ncol;
  PetscInt nvals;
  const PetscInt *cols;
  const PetscScalar *vals;
  PetscInt i, j;

  if (iRow != 0 && jCol != 0) {
    /* Set locations only */

    roffset = 0;
    coffset = 0;

    /* Equality constrained Jacobian */
    ierr = (*opflow->modelops.computeequalityconstraintjacobian)(
        opflow, opflow->X, opflow->Jac_Ge);
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

      ierr = (*opflow->modelops.computeinequalityconstraintjacobian)(
          opflow, opflow->X, opflow->Jac_Gi);
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
    ierr = (*opflow->modelops.computeequalityconstraintjacobian)(
        opflow, opflow->X, opflow->Jac_Ge);
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
      /* Compute inequality constraint jacobian */
      ierr = (*opflow->modelops.computeinequalityconstraintjacobian)(
          opflow, opflow->X, opflow->Jac_Gi);
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
  return true;
}

bool OPFLOWHIOPSPARSEInterface::eval_Hess_Lagr(
    const hiop::size_type &n, const hiop::size_type &m, const double *x,
    bool new_x, const double &obj_factor, const double *lambda, bool new_lambda,
    const hiop::size_type &nnzHSS, hiop::index_type *iRow,
    hiop::index_type *jCol, double *values) {
  PetscErrorCode ierr;
  PetscInt nrow;
  PetscInt nvals;
  const PetscInt *cols;
  const PetscScalar *vals;
  PetscInt i, j;
  PetscInt ctr = 0;

  opflow->obj_factor = obj_factor;

  if (iRow != NULL && jCol != NULL) {
    ierr = (*opflow->modelops.computehessian)(
        opflow, opflow->X, opflow->Lambdae, opflow->Lambdai, opflow->Hes);
    CHKERRQ(ierr);
    ierr = MatGetSize(opflow->Hes, &nrow, &nrow);
    CHKERRQ(ierr);

    /* Copy over locations to triplet format */
    /* Note that HIOP requires a upper triangular Hessian as oppposed
       to IPOPT which requires a lower triangular Hessian
    */
    for (i = 0; i < nrow; i++) {
      ierr = MatGetRow(opflow->Hes, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
      ctr = 0;
      for (j = 0; j < nvals; j++) {
        if (cols[j] >= i) { /* upper triangle */
          /* save as upper triangle locations */
          iRow[ctr] = i;
          jCol[ctr] = cols[j];
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

    /* Compute Hessian */
    ierr = (*opflow->modelops.computehessian)(
        opflow, opflow->X, opflow->Lambdae, opflow->Lambdai, opflow->Hes);
    CHKERRQ(ierr);

    if (opflow->modelops.computeauxhessian) {
      ierr = (*opflow->modelops.computeauxhessian)(opflow, x, opflow->Hes,
                                                   opflow->userctx);
      CHKERRQ(ierr);
      ierr = MatAssemblyBegin(opflow->Hes, MAT_FINAL_ASSEMBLY);
      CHKERRQ(ierr);
      ierr = MatAssemblyEnd(opflow->Hes, MAT_FINAL_ASSEMBLY);
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

    ierr = VecResetArray(opflow->X);
    CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Lambdae);
    CHKERRQ(ierr);
    if (opflow->Nconineq) {
      ierr = VecResetArray(opflow->Lambdai);
      CHKERRQ(ierr);
    }
  }

  return true;
}

bool OPFLOWHIOPSPARSEInterface::get_starting_point(
    const hiop::size_type &global_n, double *x0) {
  PetscErrorCode ierr;
  const PetscScalar *xarr;

  /* Set initial guess */
  ierr = (*opflow->modelops.setinitialguess)(opflow, opflow->X);
  CHKERRQ(ierr);

  ierr = VecGetArrayRead(opflow->X, &xarr);
  CHKERRQ(ierr);
  memcpy(x0, xarr, opflow->nx * sizeof(double));
  ierr = VecRestoreArrayRead(opflow->X, &xarr);
  CHKERRQ(ierr);

  return true;
}

void OPFLOWHIOPSPARSEInterface::solution_callback(
    hiop::hiopSolveStatus status, int n, const double *xsol, const double *z_L,
    const double *z_U, int m, const double *gsol, const double *lamsol,
    double obj_value) {
  PetscErrorCode ierr;
  OPFLOWSolver_HIOPSPARSE hiop = (OPFLOWSolver_HIOPSPARSE)opflow->solver;
  PetscScalar *x, *lam, *g;

  /* Copy over solution details */
  hiop->status = status;
  opflow->obj = obj_value;

  ierr = VecGetArray(opflow->X, &x);
  CHKERRV(ierr);
  memcpy(x, xsol, opflow->nx * sizeof(double));
  ierr = VecRestoreArray(opflow->X, &x);
  CHKERRV(ierr);

  if (lamsol) {
    /* HIOP returns a NULL for lamsol - probably lamsol needs to be added to
     HIOP. Need to remove this condition once it is fixed
    */
    ierr = VecGetArray(opflow->Lambda, &lam);
    CHKERRV(ierr);
    ierr =
        PetscMemcpy(lam, (double *)lamsol, opflow->ncon * sizeof(PetscScalar));
    CHKERRV(ierr);
    ierr = VecRestoreArray(opflow->Lambda, &lam);
    CHKERRV(ierr);
  } else {
    ierr = VecSet(opflow->Lambda, -9999.0);
    CHKERRV(ierr);
  }

  if (gsol) {
    /* Same situation as lamsol - gsol is NULL */
    ierr = VecGetArray(opflow->G, &g);
    CHKERRV(ierr);
    ierr = PetscMemcpy(g, (double *)gsol, opflow->ncon * sizeof(PetscScalar));
    CHKERRV(ierr);
    ierr = VecRestoreArray(opflow->G, &g);
    CHKERRV(ierr);
  }
}

bool OPFLOWHIOPSPARSEInterface::iterate_callback(
    int iter, double obj_value, double logbar_obj_value, int n, const double *x,
    const double *z_L, const double *z_U, int m_ineq, const double *s, int m,
    const double *g, const double *lambda, double inf_pr, double inf_du,
    double onenorm_pr_, double mu, double alpha_du, double alpha_pr,
    int ls_trials) {
  opflow->numits = iter;
  return true;
}

PetscErrorCode OPFLOWSolverSetUp_HIOPSPARSE(OPFLOW opflow) {
  PetscErrorCode ierr;
  OPFLOWSolver_HIOPSPARSE hiop = (OPFLOWSolver_HIOPSPARSE)opflow->solver;
  PetscBool flg1;
  int verbose_level = 3;

  PetscFunctionBegin;

  hiop->nlp = new OPFLOWHIOPSPARSEInterface(opflow);
  hiop->sp = new hiop::hiopNlpSparseIneq(*hiop->nlp);

  ierr = PetscOptionsBegin(opflow->comm->type, NULL, "HIOP options", NULL);
  CHKERRQ(ierr);
  ierr = PetscOptionsInt("-hiop_verbosity_level",
                         "HIOP verbosity level (Integer 0 to 12)", "",
                         verbose_level, &verbose_level, NULL);
  CHKERRQ(ierr);
#if defined(EXAGO_ENABLE_IPOPT)
  hiop->ipopt_debug = PETSC_FALSE;
  ierr = PetscOptionsBool("-hiop_ipopt_debug",
                          "Flag enabling debugging HIOP code with IPOPT", "",
                          hiop->ipopt_debug, &hiop->ipopt_debug, NULL);
  CHKERRQ(ierr);
#endif
  PetscOptionsEnd();

#if defined(EXAGO_ENABLE_IPOPT)
  // IPOPT Adapter
  if (hiop->ipopt_debug) {
    std::cout << "using IPOPT adapter...\n\n";
    hiop->ipoptTNLP = new hiop::hiopSparse2IpoptTNLP(hiop->nlp);
    hiop->ipoptApp = new Ipopt::IpoptApplication();

    // Using options included in HiOp's IpoptAdapter_driver.cpp
    hiop->ipoptApp->Options()->SetStringValue("recalc_y", "no");
    hiop->ipoptApp->Options()->SetStringValue("mu_strategy", "monotone");
    hiop->ipoptApp->Options()->SetNumericValue("bound_frac", 1e-8);
    hiop->ipoptApp->Options()->SetNumericValue("bound_push", 1e-8);
    hiop->ipoptApp->Options()->SetNumericValue("bound_relax_factor", 0.);
    hiop->ipoptApp->Options()->SetNumericValue("constr_mult_init_max", 0.001);
    hiop->ipoptApp->Options()->SetStringValue("derivative_test",
                                              "second-order");

    Ipopt::ApplicationReturnStatus status = hiop->ipoptApp->Initialize();

    if (status != Solve_Succeeded) {
      std::cout << std::endl
                << std::endl
                << "*** Error during initialization!" << std::endl;
      return (int)status;
    }
    PetscFunctionReturn(0);
  }
#endif

  // hiop->sp->options->SetNumericValue("fixed_var_perturb",1.0e-6);
  // hiop->sp->options->SetStringValue("mem_space","host");

  hiop->sp->options->SetStringValue("dualsUpdateType", "linear");
  //  hiop->sp->options->SetStringValue("dualsInitialization", "zero");
  hiop->sp->options->SetStringValue("fixed_var", "relax");
  hiop->sp->options->SetStringValue("Hessian", "analytical_exact");
  hiop->sp->options->SetStringValue("KKTLinsys", "xdycyd");
  hiop->sp->options->SetStringValue("compute_mode", "cpu");
  hiop->sp->options->SetIntegerValue("verbosity_level", verbose_level);
  hiop->sp->options->SetNumericValue("mu0", 1e-1);
  hiop->sp->options->SetNumericValue("tolerance", opflow->tolerance);
  hiop->sp->options->SetNumericValue("bound_relax_perturb", 1e-8);
  //  hiop->sp->options->SetStringValue("scaling_type", "none");

  hiop->solver = new hiop::hiopAlgFilterIPMNewton(hiop->sp);

  /* Error if model is not power balance hiop */
  ierr = PetscStrcmp(opflow->modelname, OPFLOWMODEL_PBPOL, &flg1);
  CHKERRQ(ierr);
  if (!flg1) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
            "Only power balance polar model allowed\n Run with -opflow_model "
            "POWER_BALANCE_POLAR\n");
    exit(1);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverSolve_HIOPSPARSE(OPFLOW opflow) {
  OPFLOWSolver_HIOPSPARSE hiop = (OPFLOWSolver_HIOPSPARSE)opflow->solver;

  PetscFunctionBegin;
#if defined(EXAGO_ENABLE_IPOPT)
  if (!hiop->ipopt_debug) {
    hiop->status = hiop->solver->run();
  } else { // Ipopt Adapter
    std::cout << "Solving with IPOPT adapter...\n\n";
    ApplicationReturnStatus status =
        hiop->ipoptApp->OptimizeTNLP(hiop->ipoptTNLP);

    if (status == Solve_Succeeded) {
      std::cout << std::endl
                << std::endl
                << "*** The problem solved!" << std::endl;
    } else {
      std::cout << std::endl
                << std::endl
                << "*** The problem FAILED!" << std::endl;
      PetscFunctionReturn(1);
    }
  }
#else
  hiop->status = hiop->solver->run();
#endif

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConvergenceStatus_HIOPSPARSE(OPFLOW opflow,
                                                           PetscBool *status) {
  OPFLOWSolver_HIOPSPARSE hiop = (OPFLOWSolver_HIOPSPARSE)opflow->solver;

  PetscFunctionBegin;
  if (hiop->status < 3)
    *status = PETSC_TRUE; /* See hiopInterface.hpp. The first three denote
                             convergence */
  else
    *status = PETSC_FALSE;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetObjective_HIOPSPARSE(OPFLOW opflow,
                                                   PetscReal *obj) {
  PetscFunctionBegin;
  *obj = opflow->obj;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetSolution_HIOPSPARSE(OPFLOW opflow, Vec *X) {
  PetscFunctionBegin;
  *X = opflow->X;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConstraints_HIOPSPARSE(OPFLOW opflow, Vec *G) {
  PetscFunctionBegin;
  *G = opflow->G;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConstraintMultipliers_HIOPSPARSE(OPFLOW opflow,
                                                               Vec *Lambda) {
  PetscFunctionBegin;
  *Lambda = opflow->Lambda;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverDestroy_HIOPSPARSE(OPFLOW opflow) {
  PetscErrorCode ierr;
  OPFLOWSolver_HIOPSPARSE hiop = (OPFLOWSolver_HIOPSPARSE)opflow->solver;

  PetscFunctionBegin;

  delete hiop->sp;
  delete hiop->nlp;

  ierr = PetscFree(hiop);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverCreate_HIOPSPARSE(OPFLOW opflow) {
  PetscErrorCode ierr;
  OPFLOWSolver_HIOPSPARSE hiop;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1, &hiop);
  CHKERRQ(ierr);

  opflow->solver = hiop;

  opflow->solverops.setup = OPFLOWSolverSetUp_HIOPSPARSE;
  opflow->solverops.solve = OPFLOWSolverSolve_HIOPSPARSE;
  opflow->solverops.destroy = OPFLOWSolverDestroy_HIOPSPARSE;
  opflow->solverops.getobjective = OPFLOWSolverGetObjective_HIOPSPARSE;
  opflow->solverops.getconvergencestatus =
      OPFLOWSolverGetConvergenceStatus_HIOPSPARSE;
  opflow->solverops.getsolution = OPFLOWSolverGetSolution_HIOPSPARSE;
  opflow->solverops.getconstraints = OPFLOWSolverGetConstraints_HIOPSPARSE;
  opflow->solverops.getconstraintmultipliers =
      OPFLOWSolverGetConstraintMultipliers_HIOPSPARSE;

  opflow->eqcons_scaling = 1e6;
  PetscFunctionReturn(0);
}

#endif
#endif
