
#include <exago_config.h>
#if defined(EXAGO_ENABLE_HIOP)
#include "opflow_hiop.h"
#include <private/opflowimpl.h>

typedef enum { AUTO = 0, CPU = 1, HYBRID = 2, GPU = 3 } HIOPComputeMode;
const char *HIOPComputeModeChoices[] = {
    "auto", "cpu", "hybrid", "gpu", "HIOPComputeModeChoices", "", 0};

const char *HIOPMemSpaceChoices[] = {
    "default", "host", "um", "device", "HIOPMemSpaceChoices", "", 0};

OPFLOWHIOPInterface::OPFLOWHIOPInterface(OPFLOW opflowin) { opflow = opflowin; }

bool OPFLOWHIOPInterface::get_prob_sizes(hiop::size_type &n,
                                         hiop::size_type &m) {
  n = opflow->nx;
  m = opflow->ncon;
  return true;
}

bool OPFLOWHIOPInterface::get_vars_info(const hiop::size_type &n, double *xlow,
                                        double *xupp, NonlinearityType *type) {
  PetscInt i;
  PetscErrorCode ierr;
  //  PetscPrintf(MPI_COMM_SELF,"Enter get_vars_info\n");

  ierr = (*opflow->modelops.setvariableboundsarray)(opflow, xlow, xupp);
  CHKERRQ(ierr);

  for (i = 0; i < n; i++) {
    type[i] = hiopNonlinear;
  }
  //  PetscPrintf(MPI_COMM_SELF,"Exit get_vars_info\n");
  return true;
}

bool OPFLOWHIOPInterface::get_cons_info(const hiop::size_type &m, double *clow,
                                        double *cupp, NonlinearityType *type) {
  PetscInt i;
  PetscErrorCode ierr;
  //  PetscPrintf(MPI_COMM_SELF,"Enter get_cons_info \n");

  ierr = (*opflow->modelops.setconstraintboundsarray)(opflow, clow, cupp);
  CHKERRQ(ierr);

  for (i = 0; i < m; i++)
    type[i] = hiopNonlinear;

  //  PetscPrintf(MPI_COMM_SELF,"Exit get_cons_info \n");

  return true;
}

bool OPFLOWHIOPInterface::get_sparse_dense_blocks_info(
    int &nx_sparse, int &nx_dense, int &nnz_sparse_Jace, int &nnz_sparse_Jaci,
    int &nnz_sparse_Hess_Lagr_SS, int &nnz_sparse_Hess_Lagr_SD) {
  //  PetscPrintf(MPI_COMM_SELF,"Enter sparse_dense_blocks_info \n");

  nx_sparse = opflow->nxsparse;
  nx_dense = opflow->nxdense;

  nnz_sparse_Jace = opflow->nnz_eqjacsp;
  nnz_sparse_Jaci = opflow->nnz_ineqjacsp;
  nnz_sparse_Hess_Lagr_SS = opflow->nnz_hesssp;
  nnz_sparse_Hess_Lagr_SD = 0;

  //  PetscPrintf(MPI_COMM_SELF,"Enter sparse_dense_blocks_info \n");

  return true;
}

bool OPFLOWHIOPInterface::eval_f(const hiop::size_type &n, const double *x,
                                 bool new_x, double &obj_value) {
  PetscErrorCode ierr;

  //  PetscPrintf(MPI_COMM_SELF,"Enter eval_f \n");

  obj_value = 0.0;

  /* Compute objective */
  ierr = (*opflow->modelops.computeobjectivearray)(opflow, x, &obj_value);
  CHKERRQ(ierr);

  //  PetscPrintf(MPI_COMM_SELF,"Exit eval_f \n");

  return true;
}

bool OPFLOWHIOPInterface::eval_cons(const hiop::size_type &n,
                                    const hiop::size_type &m,
                                    const hiop::size_type &num_cons,
                                    const hiop::size_type *idx_cons,
                                    const double *x, bool new_x, double *cons) {
  PetscErrorCode ierr;
  OPFLOWSolver_HIOP hiop = (OPFLOWSolver_HIOP)opflow->solver;

  //  PetscPrintf(MPI_COMM_SELF,"Enter eval_cons \n");

  if (!num_cons) {
    hiop->cons_call++;
    return true;
  }

  if (hiop->cons_call % 2 == 0) {
    /* Equality constaints */
    ierr = PetscLogEventBegin(opflow->eqconslogger, 0, 0, 0, 0);
    CHKERRQ(ierr);

    ierr = (*opflow->modelops.computeequalityconstraintsarray)(opflow, x, cons);
    CHKERRQ(ierr);
    ierr = PetscLogEventEnd(opflow->eqconslogger, 0, 0, 0, 0);
    CHKERRQ(ierr);
    hiop->cons_call++;
  } else {
    if (opflow->nconineq) {
      ierr = PetscLogEventBegin(opflow->ineqconslogger, 0, 0, 0, 0);
      CHKERRQ(ierr);
      /* Inequality constraints */
      ierr = (*opflow->modelops.computeinequalityconstraintsarray)(opflow, x,
                                                                   cons);
      CHKERRQ(ierr);
      ierr = PetscLogEventEnd(opflow->ineqconslogger, 0, 0, 0, 0);
      CHKERRQ(ierr);
      hiop->cons_call++;
    }
  }

  //  PetscPrintf(MPI_COMM_SELF,"Exit eval_cons \n");
  return true;
}

bool OPFLOWHIOPInterface::eval_grad_f(const hiop::size_type &n, const double *x,
                                      bool new_x, double *gradf) {
  PetscErrorCode ierr;
  //  PetscPrintf(MPI_COMM_SELF,"Enter eval_grad_f \n");
  ierr = PetscLogEventBegin(opflow->gradlogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computegradientarray)(opflow, x, gradf);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->gradlogger, 0, 0, 0, 0);
  CHKERRQ(ierr);

  //  PetscPrintf(MPI_COMM_SELF,"Exit eval_grad_f \n");

  return true;
}

bool OPFLOWHIOPInterface::eval_Jac_cons(
    const hiop::size_type &n, const hiop::size_type &m,
    const hiop::size_type &num_cons, const hiop::size_type *idx_cons,
    const double *x, bool new_x, const hiop::size_type &nsparse,
    const hiop::size_type &ndense, const int &nnzJacS, int *iJacS, int *jJacS,
    double *MJacS, double *JacD) {
  PetscErrorCode ierr;
  OPFLOWSolver_HIOP hiop = (OPFLOWSolver_HIOP)opflow->solver;
  //  PetscPrintf(MPI_COMM_SELF,"Enter eval_Jac_cons \n");

  if (!num_cons) {
    hiop->cons_call++;
    return true;
  }

  if (hiop->cons_call % 2 == 0) {
    /* Equality constraints */
    //    PetscPrintf(MPI_COMM_SELF,"Came here eq. \n");

    /* Sparse Jacobian */
    ierr = (*opflow->modelops.computesparseequalityconstraintjacobianhiop)(
        opflow, x, iJacS, jJacS, MJacS);
    CHKERRQ(ierr);

    ierr = PetscLogEventBegin(opflow->denseeqconsjaclogger, 0, 0, 0, 0);
    CHKERRQ(ierr);

    /* Dense equality constraint Jacobian */
    ierr = (*opflow->modelops.computedenseequalityconstraintjacobianhiop)(
        opflow, x, JacD);
    CHKERRQ(ierr);

    ierr = PetscLogEventEnd(opflow->denseeqconsjaclogger, 0, 0, 0, 0);
    CHKERRQ(ierr);
    hiop->cons_call++;
  } else {
    /* Dense inequality constraint Jacobian */
    //    PetscPrintf(MPI_COMM_SELF,"Came here ineq. \n");

    if (opflow->nconineq) {
      /* Sparse Inequality constraint Jacobian */
      ierr = (*opflow->modelops.computesparseinequalityconstraintjacobianhiop)(
          opflow, x, iJacS, jJacS, MJacS);
      CHKERRQ(ierr);

      ierr = PetscLogEventBegin(opflow->denseineqconsjaclogger, 0, 0, 0, 0);
      CHKERRQ(ierr);
      ierr = (*opflow->modelops.computedenseinequalityconstraintjacobianhiop)(
          opflow, x, JacD);
      CHKERRQ(ierr);
      ierr = PetscLogEventEnd(opflow->denseineqconsjaclogger, 0, 0, 0, 0);
      CHKERRQ(ierr);
      hiop->cons_call++;
    }
  }
  //  PetscPrintf(MPI_COMM_SELF,"Exit eval_Jac_cons \n");

  return true;
}

bool OPFLOWHIOPInterface::eval_Hess_Lagr(
    const hiop::size_type &n, const hiop::size_type &m, const double *x,
    bool new_x, const double &obj_factor, const double *lambda, bool new_lambda,
    const hiop::size_type &nsparse, const hiop::size_type &ndense,
    const int &nnzHSS, int *iHSS, int *jHSS, double *MHSS, double *HDD,
    int &nnzHSD, int *iHSD, int *jHSD, double *MHSD) {
  PetscErrorCode ierr;
  //  PetscPrintf(MPI_COMM_SELF,"Enter eval_Hess_Lagr \n");

  opflow->obj_factor = obj_factor;

  /* Compute sparse hessian */
  ierr = PetscLogEventBegin(opflow->sparsehesslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computesparsehessianhiop)(opflow, x, lambda, iHSS,
                                                      jHSS, MHSS);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->sparsehesslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);

  ierr = PetscLogEventBegin(opflow->densehesslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  /* Compute dense hessian */
  ierr = (*opflow->modelops.computedensehessianhiop)(opflow, x, lambda, HDD);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->densehesslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  //  PetscPrintf(MPI_COMM_SELF,"Exit eval_Hess_Lagr \n");

  return true;
}

bool OPFLOWHIOPInterface::get_starting_point(const hiop::size_type &global_n,
                                             double *x0) {
  PetscErrorCode ierr;
  //  PetscPrintf(MPI_COMM_SELF,"Enter get_starting_point \n");

  ierr = (*opflow->modelops.setinitialguessarray)(opflow, x0);
  CHKERRQ(ierr);

  //  PetscPrintf(MPI_COMM_SELF,"Exit get_starting_point \n");

  return true;
}

bool OPFLOWHIOPInterface::iterate_callback(
    int iter, double obj_value, double logbar_obj_value, int n, const double *x,
    const double *z_L, const double *z_U, int m_ineq, const double *s, int m,
    const double *g, const double *lambda, double inf_pr, double inf_du,
    double mu, double onenorm_pr_, double alpha_du, double alpha_pr,
    int ls_trials) {
  opflow->numits = iter;
  return true;
}

void OPFLOWHIOPInterface::solution_callback(
    hiop::hiopSolveStatus status, int n, const double *xsol, const double *z_L,
    const double *z_U, int m, const double *gsol, const double *lamsol,
    double obj_value) {
  PetscErrorCode ierr;
  OPFLOWSolver_HIOP hiop = (OPFLOWSolver_HIOP)opflow->solver;

  /* Copy over solution details */
  hiop->status = status;
  opflow->obj = obj_value;

  if (opflow->modelops.solutioncallbackhiop) {
    ierr = (*opflow->modelops.solutioncallbackhiop)(opflow, xsol, z_L, z_U,
                                                    gsol, lamsol, obj_value);
    CHKERRV(ierr);
  }
}

PetscErrorCode OPFLOWSolverSetUp_HIOP(OPFLOW opflow) {
  PetscErrorCode ierr;
  OPFLOWSolver_HIOP hiop = (OPFLOWSolver_HIOP)opflow->solver;
  PetscBool ismodelpbpolhiop, ismodelpbpolrajahiop;
  HIOPComputeMode compute_mode = AUTO;
#ifndef EXAGO_ENABLE_GPU
  opflow->mem_space = DEFAULT;
#else
  opflow->mem_space = DEVICE;
#endif
  hiop->cons_call = 0;

  int verbose_level = OPFLOWOptions::hiop_verbosity_level.default_value;
  PetscBool mode_set = PETSC_FALSE;

  PetscFunctionBegin;

  hiop->nlp = new OPFLOWHIOPInterface(opflow);
  hiop->mds = new hiop::hiopNlpMDS(*hiop->nlp);

  PetscOptionsBegin(opflow->comm->type, NULL, "HIOP options", NULL);

  ierr = PetscOptionsEnum(OPFLOWOptions::hiop_compute_mode.opt.c_str(),
                          OPFLOWOptions::hiop_compute_mode.desc.c_str(), "",
                          HIOPComputeModeChoices, (PetscEnum)compute_mode,
                          (PetscEnum *)&compute_mode, &mode_set);
  CHKERRQ(ierr);

  if (mode_set == PETSC_FALSE) {
    hiop->mds->options->SetStringValue("compute_mode",
                                       opflow->_p_hiop_compute_mode.c_str());
  } else {
    hiop->mds->options->SetStringValue("compute_mode",
                                       HIOPComputeModeChoices[compute_mode]);
  }

  ierr = PetscOptionsInt(OPFLOWOptions::hiop_verbosity_level.opt.c_str(),
                         OPFLOWOptions::hiop_verbosity_level.desc.c_str(), "",
                         verbose_level, &opflow->_p_hiop_verbosity_level, NULL);
  CHKERRQ(ierr);

  ierr = PetscOptionsEnum(OPFLOWOptions::hiop_mem_space.opt.c_str(),
                          OPFLOWOptions::hiop_mem_space.desc.c_str(), "",
                          HIOPMemSpaceChoices, (PetscEnum)opflow->mem_space,
                          (PetscEnum *)&opflow->mem_space, &mode_set);
  CHKERRQ(ierr);

#if defined(EXAGO_ENABLE_IPOPT)
  hiop->ipopt_debug =
      static_cast<PetscBool>(OPFLOWOptions::hiop_ipopt_debug.default_value);
  ierr = PetscOptionsBool(OPFLOWOptions::hiop_ipopt_debug.opt.c_str(),
                          OPFLOWOptions::hiop_ipopt_debug.desc.c_str(), "",
                          hiop->ipopt_debug, &hiop->ipopt_debug, NULL);
  CHKERRQ(ierr);
#endif
  PetscOptionsEnd();

#if defined(EXAGO_ENABLE_IPOPT)
  // IPOPT Adapter
  if (hiop->ipopt_debug) {
    std::cout << "using IPOPT adapter...\n\n";
    hiop->ipoptTNLP = new hiop::hiopMDS2IpoptTNLP(hiop->nlp);
    hiop->ipoptApp = new Ipopt::IpoptApplication();

    // Using options included in HiOp's IpoptAdapter_driver.cpp
    hiop->ipoptApp->Options()->SetStringValue("recalc_y", "no");
    hiop->ipoptApp->Options()->SetStringValue("mu_strategy", "monotone");
    hiop->ipoptApp->Options()->SetNumericValue("bound_frac", 1e-8);
    hiop->ipoptApp->Options()->SetNumericValue("bound_push", 1e-8);
    hiop->ipoptApp->Options()->SetNumericValue("bound_relax_factor", 0.);
    hiop->ipoptApp->Options()->SetNumericValue("bound_relax_perturb", 1e-8);
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

  hiop->mds->options->SetStringValue("duals_update_type", "linear");
  hiop->mds->options->SetStringValue("duals_init", "zero");

  hiop->mds->options->SetStringValue("fixed_var", "relax");
  hiop->mds->options->SetStringValue("Hessian", "analytical_exact");
  hiop->mds->options->SetStringValue("KKTLinsys", "xdycyd");
  hiop->mds->options->SetIntegerValue("verbosity_level",
                                      opflow->_p_hiop_verbosity_level);
  hiop->mds->options->SetNumericValue("mu0", 1e-1);
  hiop->mds->options->SetNumericValue("tolerance", opflow->tolerance);
  hiop->mds->options->SetNumericValue("bound_relax_perturb", 1e-6);
  hiop->mds->options->SetStringValue("scaling_type", "none");

  /* Error if model is not power balance hiop or power balance raja hiop */
  ismodelpbpolhiop =
      static_cast<PetscBool>(opflow->modelname == OPFLOWMODEL_PBPOLHIOP);
#if defined(EXAGO_ENABLE_RAJA)
  ismodelpbpolrajahiop =
      static_cast<PetscBool>(opflow->modelname == OPFLOWMODEL_PBPOLRAJAHIOP);
#endif
  if (!ismodelpbpolhiop && !ismodelpbpolrajahiop) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
            "%s opflow model not supported with HIOP\n", opflow->modelname);
    PetscFunctionReturn(1);
    exit(0);
  }

  if (ismodelpbpolrajahiop) {
    hiop->mds->options->SetStringValue("mem_space",
                                       HIOPMemSpaceChoices[opflow->mem_space]);
  }

  //  ierr = PetscPrintf(MPI_COMM_SELF,"Came in OPFLOWSetUp\n");CHKERRQ(ierr);
  hiop->solver = new hiop::hiopAlgFilterIPMNewton(hiop->mds);

  //  ierr = PetscPrintf(MPI_COMM_SELF,"Exit OPFLOWSetUp\n");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverSolve_HIOP(OPFLOW opflow) {
  OPFLOWSolver_HIOP hiop = (OPFLOWSolver_HIOP)opflow->solver;

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
    }
  }
#else
  hiop->status = hiop->solver->run();
#endif

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConvergenceStatus_HIOP(OPFLOW opflow,
                                                     PetscBool *status) {
  OPFLOWSolver_HIOP hiop = (OPFLOWSolver_HIOP)opflow->solver;

  PetscFunctionBegin;
  if (hiop->status < 3)
    *status = PETSC_TRUE; /* See hiopInterface.hpp. The first three denote
                             convergence */
  else
    *status = PETSC_FALSE;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetObjective_HIOP(OPFLOW opflow, PetscReal *obj) {
  PetscFunctionBegin;
  *obj = opflow->obj;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetSolution_HIOP(OPFLOW opflow, Vec *X) {
  PetscFunctionBegin;
  *X = opflow->X;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConstraints_HIOP(OPFLOW opflow, Vec *G) {
  PetscFunctionBegin;
  *G = opflow->G;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConstraintMultipliers_HIOP(OPFLOW opflow,
                                                         Vec *Lambda) {
  PetscFunctionBegin;
  *Lambda = opflow->Lambda;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverDestroy_HIOP(OPFLOW opflow) {
  PetscErrorCode ierr;
  OPFLOWSolver_HIOP hiop = (OPFLOWSolver_HIOP)opflow->solver;

  PetscFunctionBegin;

  delete hiop->solver;
  delete hiop->mds;
  delete hiop->nlp;

  ierr = PetscFree(hiop);
  CHKERRQ(ierr);

#ifdef EXAGO_ENABLE_GPU
  magma_finalize();
#endif

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverCreate_HIOP(OPFLOW opflow) {
  PetscErrorCode ierr;
  OPFLOWSolver_HIOP hiop;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1, &hiop);
  CHKERRQ(ierr);

  opflow->solver = hiop;

  opflow->solverops.setup = OPFLOWSolverSetUp_HIOP;
  opflow->solverops.solve = OPFLOWSolverSolve_HIOP;
  opflow->solverops.destroy = OPFLOWSolverDestroy_HIOP;
  opflow->solverops.getobjective = OPFLOWSolverGetObjective_HIOP;
  opflow->solverops.getconvergencestatus =
      OPFLOWSolverGetConvergenceStatus_HIOP;
  opflow->solverops.getsolution = OPFLOWSolverGetSolution_HIOP;
  opflow->solverops.getconstraints = OPFLOWSolverGetConstraints_HIOP;
  opflow->solverops.getconstraintmultipliers =
      OPFLOWSolverGetConstraintMultipliers_HIOP;

#ifdef EXAGO_ENABLE_GPU
  magma_init();
#endif
  PetscFunctionReturn(0);
}

#endif
