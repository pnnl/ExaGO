#include <exago_config.h>
#if defined(EXAGO_ENABLE_HIOP)
#if defined(EXAGO_ENABLE_HIOP_SPARSE)

#include <private/opflowimpl.h>
#include "opflow_hiopsparsegpu.hpp"

OPFLOWHIOPSPARSEGPUInterface::OPFLOWHIOPSPARSEGPUInterface(OPFLOW opflowin) {
  opflow = opflowin;
}

bool OPFLOWHIOPSPARSEGPUInterface::get_prob_sizes(hiop::size_type &n, hiop::size_type &m) {
  n = opflow->nx;
  m = opflow->ncon;
  return true;
}

bool OPFLOWHIOPSPARSEGPUInterface::get_vars_info(const hiop::size_type &n,
                                                 double *xlow, double *xupp,
                                                 NonlinearityType *type) {
  PetscErrorCode ierr;
  PetscInt i;

  ierr = (*opflow->modelops.setvariableboundsarray)(opflow, xlow, xupp);
  CHKERRQ(ierr);

  for (i = 0; i < n; i++) {
    type[i] = hiopNonlinear;
  }

  return true;
}

bool OPFLOWHIOPSPARSEGPUInterface::get_cons_info(const hiop::size_type &m,
                                                 double *clow, double *cupp,
                                                 NonlinearityType *type) {
  PetscInt i;
  PetscErrorCode ierr;

  ierr = (*opflow->modelops.setconstraintboundsarray)(opflow, clow, cupp);
  CHKERRQ(ierr);

  for (i = 0; i < m; i++)
    type[i] = hiopNonlinear;

  return true;
}

bool OPFLOWHIOPSPARSEGPUInterface::get_sparse_blocks_info(
							  hiop::size_type &nx, hiop::size_type &nnz_sparse_Jaceq, hiop::size_type &nnz_sparse_Jacineq,
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

  nnz_sparse_Jaceq = opflow->nnz_eqjacsp = info_eq.nz_used;

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

    nnz_sparse_Jacineq = opflow->nnz_ineqjacsp = info_ineq.nz_used;
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

  opflow->nnz_hesssp = nnz_sparse_Hess_Lagr;

  return true;
}

bool OPFLOWHIOPSPARSEGPUInterface::eval_f(const hiop::size_type &n, const double *x,
                                          bool new_x, double &obj_value) {
  PetscErrorCode ierr;

  obj_value = 0.0;

  /* Compute objective */
  ierr = (*opflow->modelops.computeobjectivearray)(opflow, x, &obj_value);
  CHKERRQ(ierr);

  return true;
}

bool OPFLOWHIOPSPARSEGPUInterface::eval_cons(const hiop::size_type &n,
                                          const hiop::size_type &m,
                                          const hiop::size_type &num_cons,
                                          const hiop::size_type *idx_cons,
                                          const double *x, bool new_x,
                                          double *cons) {
  return false;
}

bool OPFLOWHIOPSPARSEGPUInterface::eval_cons(const hiop::size_type &n,
                                             const hiop::size_type &m,
                                             const double *x, bool new_x,
                                             double *cons) {
  PetscErrorCode ierr;

  /* Equality constaints */
  ierr = PetscLogEventBegin(opflow->eqconslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  
  ierr = (*opflow->modelops.computeequalityconstraintsarray)(opflow, x, cons);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->eqconslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);

  if (opflow->nconineq) {
    ierr = PetscLogEventBegin(opflow->ineqconslogger, 0, 0, 0, 0);
    CHKERRQ(ierr);
    /* Inequality constraints */
    ierr = (*opflow->modelops.computeinequalityconstraintsarray)(opflow, x,
								 cons + opflow->nconeq);
    CHKERRQ(ierr);
    ierr = PetscLogEventEnd(opflow->ineqconslogger, 0, 0, 0, 0);
    CHKERRQ(ierr);
  }

  return true;
}

bool OPFLOWHIOPSPARSEGPUInterface::eval_grad_f(const hiop::size_type &n,
                                               const double *x, bool new_x,
                                               double *gradf) {
  PetscErrorCode ierr;

  ierr = PetscLogEventBegin(opflow->gradlogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computegradientarray)(opflow, x, gradf);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->gradlogger, 0, 0, 0, 0);
  CHKERRQ(ierr);

  return true;
}

bool OPFLOWHIOPSPARSEGPUInterface::eval_Jac_cons(
    const hiop::size_type &n, const hiop::size_type &m, const hiop::size_type &num_cons,
    const hiop::size_type *idx_cons, const double *x, bool new_x, const hiop::size_type &nnzJacS,
    hiop::index_type *iJacS_dev, hiop::index_type *jJacS_dev, double *MJacS_dev) {
  return false;
}

bool OPFLOWHIOPSPARSEGPUInterface::eval_Jac_cons(const hiop::size_type &n,
                                                 const hiop::size_type &m,
                                                 const double *x, bool new_x,
                                                 const hiop::size_type &nnzJacS, hiop::index_type *iJacS_dev,
                                                 hiop::index_type *jJacS_dev, double *MJacS_dev) {
  PetscErrorCode ierr;

  /* Sparse Jacobian */
  ierr = (*opflow->modelops.computesparseequalityconstraintjacobianhiop)(
           opflow, x, iJacS_dev, jJacS_dev, MJacS_dev);
  CHKERRQ(ierr);

  if (opflow->nconineq) {
    /* Sparse Inequality constraint Jacobian */
    ierr = (*opflow->modelops.computesparseinequalityconstraintjacobianhiop)(opflow, x, iJacS_dev, jJacS_dev, MJacS_dev);
    CHKERRQ(ierr);
  }

  return true;
}

bool OPFLOWHIOPSPARSEGPUInterface::eval_Hess_Lagr(
    const hiop::size_type &n, const hiop::size_type &m, const double *x_dev, bool new_x,
    const double &obj_factor, const double *lambda_dev, bool new_lambda,
    const hiop::size_type &nnzHSS, hiop::index_type *iHSS_dev, hiop::index_type *jHSS_dev, double *MHSS_dev) {
  PetscErrorCode ierr;

  opflow->obj_factor = obj_factor;

  /* Compute sparse hessian */
  ierr = PetscLogEventBegin(opflow->sparsehesslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computesparsehessianhiop)(opflow, x_dev, lambda_dev, iHSS_dev,
                                                      jHSS_dev, MHSS_dev);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->sparsehesslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);

  return true;
}

bool OPFLOWHIOPSPARSEGPUInterface::get_starting_point(const hiop::size_type &global_n,
                                                      double *x0) {
  PetscErrorCode ierr;
  const PetscScalar *xarr;

  /* Set initial guess */
  ierr = (*opflow->modelops.setinitialguessarray)(opflow, x0);
  CHKERRQ(ierr);

  return true;
}

void OPFLOWHIOPSPARSEGPUInterface::solution_callback(
    hiop::hiopSolveStatus status, int n, const double *xsol, const double *z_L,
    const double *z_U, int m, const double *gsol, const double *lamsol,
    double obj_value) {
  PetscErrorCode ierr;
  OPFLOWSolver_HIOPSPARSEGPU hiop = (OPFLOWSolver_HIOPSPARSEGPU)opflow->solver;
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
        PetscMemcpy((double *)lamsol, lam, opflow->ncon * sizeof(PetscScalar));
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
    ierr = PetscMemcpy((double *)gsol, g, opflow->ncon * sizeof(PetscScalar));
    CHKERRV(ierr);
    ierr = VecRestoreArray(opflow->G, &g);
    CHKERRV(ierr);
  }
}

bool OPFLOWHIOPSPARSEGPUInterface::iterate_callback(
    int iter, double obj_value, double logbar_obj_value, int n, const double *x,
    const double *z_L, const double *z_U, int m_ineq, const double *s, int m,
    const double *g, const double *lambda, double inf_pr, double inf_du,
    double onenorm_pr_, double mu, double alpha_du, double alpha_pr,
    int ls_trials) {
  opflow->numits = iter;
  return true;
}

PetscErrorCode OPFLOWSolverSetUp_HIOPSPARSEGPU(OPFLOW opflow) {
  PetscErrorCode ierr;
  OPFLOWSolver_HIOPSPARSEGPU hiop = (OPFLOWSolver_HIOPSPARSEGPU)opflow->solver;
  PetscBool flg1;
  int verbose_level = 3;

  PetscFunctionBegin;

  hiop->nlp = new OPFLOWHIOPSPARSEGPUInterface(opflow);
  hiop->sp = new hiop::hiopNlpSparse(*hiop->nlp);

  hiop->ipopt_debug = PETSC_FALSE;

  PetscOptionsBegin(opflow->comm->type, NULL, "HIOP options", NULL);

  ierr = PetscOptionsInt("-hiop_verbosity_level",
                         "HIOP verbosity level (Integer 0 to 12)", "",
                         verbose_level, &verbose_level, NULL);
  CHKERRQ(ierr);
#if defined(EXAGO_ENABLE_IPOPT)
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

  hiop->sp->options->SetStringValue("mem_space","device");
  hiop->sp->options->SetStringValue("compute_mode", "gpu");

  hiop->sp->options->SetStringValue("dualsInitialization", "zero");
  hiop->sp->options->SetStringValue("duals_init", "zero");

  hiop->sp->options->SetStringValue("fact_acceptor", "inertia_free");
  hiop->sp->options->SetStringValue("linsol_mode", "speculative");
  hiop->sp->options->SetStringValue("fixed_var", "relax");
  hiop->sp->options->SetStringValue("Hessian", "analytical_exact");
  hiop->sp->options->SetStringValue("KKTLinsys", "xdycyd");

  hiop->sp->options->SetIntegerValue("verbosity_level", verbose_level);
  hiop->sp->options->SetNumericValue("mu0", 1e-1);
  hiop->sp->options->SetNumericValue("tolerance", opflow->tolerance);
  hiop->sp->options->SetNumericValue("bound_relax_perturb", 1e-4);
  hiop->sp->options->SetStringValue("scaling_type", "none");

  hiop->solver = new hiop::hiopAlgFilterIPMNewton(hiop->sp);

  /* Error if model is not power balance hiop */
  ierr = PetscStrcmp(opflow->modelname.c_str(), OPFLOWMODEL_PBPOLRAJAHIOPSPARSE, &flg1);
  CHKERRQ(ierr);
  if (!flg1) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
            "Only PBPOLRAJAHIOPSPARSE model allowed with solver HIOPSPARSEGPU \n Run with -opflow_model "
            "PBPOLRAJAHIOPSPARSE\n");
    exit(1);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverSolve_HIOPSPARSEGPU(OPFLOW opflow) {
  OPFLOWSolver_HIOPSPARSEGPU hiop = (OPFLOWSolver_HIOPSPARSEGPU)opflow->solver;

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

PetscErrorCode
OPFLOWSolverGetConvergenceStatus_HIOPSPARSEGPU(OPFLOW opflow,
                                               PetscBool *status) {
  OPFLOWSolver_HIOPSPARSEGPU hiop = (OPFLOWSolver_HIOPSPARSEGPU)opflow->solver;

  PetscFunctionBegin;
  if (hiop->status < 3)
    *status = PETSC_TRUE; /* See hiopInterface.hpp. The first three denote
                             convergence */
  else
    *status = PETSC_FALSE;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetObjective_HIOPSPARSEGPU(OPFLOW opflow,
                                                      PetscReal *obj) {
  PetscFunctionBegin;
  *obj = opflow->obj;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetSolution_HIOPSPARSEGPU(OPFLOW opflow, Vec *X) {
  PetscFunctionBegin;
  *X = opflow->X;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConstraints_HIOPSPARSEGPU(OPFLOW opflow, Vec *G) {
  PetscFunctionBegin;
  *G = opflow->G;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConstraintMultipliers_HIOPSPARSEGPU(OPFLOW opflow,
                                                                  Vec *Lambda) {
  PetscFunctionBegin;
  *Lambda = opflow->Lambda;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverDestroy_HIOPSPARSEGPU(OPFLOW opflow) {
  PetscErrorCode ierr;
  OPFLOWSolver_HIOPSPARSEGPU hiop = (OPFLOWSolver_HIOPSPARSEGPU)opflow->solver;

  PetscFunctionBegin;

  delete hiop->solver;
  delete hiop->sp;
  delete hiop->nlp;

  ierr = PetscFree(hiop);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverCreate_HIOPSPARSEGPU(OPFLOW opflow) {
  PetscErrorCode ierr;
  OPFLOWSolver_HIOPSPARSEGPU hiop;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1, &hiop);
  CHKERRQ(ierr);

  opflow->solver = hiop;

  opflow->solverops.setup = OPFLOWSolverSetUp_HIOPSPARSEGPU;
  opflow->solverops.solve = OPFLOWSolverSolve_HIOPSPARSEGPU;
  opflow->solverops.destroy = OPFLOWSolverDestroy_HIOPSPARSEGPU;
  opflow->solverops.getobjective = OPFLOWSolverGetObjective_HIOPSPARSEGPU;
  opflow->solverops.getconvergencestatus =
      OPFLOWSolverGetConvergenceStatus_HIOPSPARSEGPU;
  opflow->solverops.getsolution = OPFLOWSolverGetSolution_HIOPSPARSEGPU;
  opflow->solverops.getconstraints = OPFLOWSolverGetConstraints_HIOPSPARSEGPU;
  opflow->solverops.getconstraintmultipliers =
      OPFLOWSolverGetConstraintMultipliers_HIOPSPARSEGPU;

  PetscFunctionReturn(0);
}

#endif
#endif
