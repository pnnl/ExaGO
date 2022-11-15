#include <iostream>
#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>
#include <private/tcopflowimpl.h>
#include <utils.h>

/**
 * @brief Creates a security constrained optimal power flow application object
 *
 * @param[in] mpicomm the MPI communicator
 * @param[in] scopflowout pointer to the security constrained optimal power flow
 * application object
 */
PetscErrorCode SCOPFLOWCreate(MPI_Comm mpicomm, SCOPFLOW *scopflowout) {
  PetscErrorCode ierr;
  SCOPFLOW scopflow;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1, &scopflow);
  CHKERRQ(ierr);

  ierr = COMMCreate(mpicomm, &scopflow->comm);
  CHKERRQ(ierr);

  scopflow->Nconeq = scopflow->nconeq = 0;
  scopflow->Nconineq = scopflow->nconineq = 0;
  scopflow->Ncon = scopflow->ncon = 0;
  scopflow->Nx = scopflow->nx = 0;
  scopflow->Nc = scopflow->nc = SCOPFLOWOptions::Nc.default_value;
  scopflow->Gi = NULL;
  scopflow->Lambdai = NULL;

  scopflow->obj_factor = 1.0;
  scopflow->objbase = scopflow->objtot = 0.0;

  scopflow->mode = SCOPFLOWOptions::mode.default_value;
  scopflow->tolerance = 1e-6;

  scopflow->scen = NULL;

  scopflow->ismultiperiod = SCOPFLOWOptions::enable_multiperiod.default_value;

  scopflow->opflow0 = NULL;
  scopflow->opflows = NULL;

  scopflow->solver = NULL;
  scopflow->model = NULL;
  scopflow->solverset = PETSC_FALSE;

  scopflow->nmodelsregistered = 0;
  scopflow->SCOPFLOWModelRegisterAllCalled = PETSC_FALSE;

  /* Register all models */
  ierr = SCOPFLOWModelRegisterAll(scopflow);

  scopflow->nsolversregistered = 0;
  scopflow->SCOPFLOWSolverRegisterAllCalled = PETSC_FALSE;

  scopflow->enable_powerimbalance_variables =
      SCOPFLOWOptions::enable_powerimbalance.default_value;
  scopflow->ignore_lineflow_constraints =
      SCOPFLOWOptions::ignore_lineflow_constraints.default_value;

  /* Register all solvers */
  ierr = SCOPFLOWSolverRegisterAll(scopflow);

  /* Run-time options */
  scopflow->iscoupling = SCOPFLOWOptions::iscoupling.default_value;

  scopflow->ctgclist = NULL;
  scopflow->ctgcfileset = PETSC_FALSE;

  /* Default subproblem model and solver */
  (void)std::strncpy(scopflow->subproblem_model,
                     SCOPFLOWOptions::subproblem_model.default_value.c_str(),
                     sizeof(scopflow->subproblem_model));
  (void)std::strncpy(scopflow->subproblem_solver,
                     SCOPFLOWOptions::subproblem_solver.default_value.c_str(),
                     sizeof(scopflow->subproblem_solver));
#ifdef EXAGO_ENABLE_HIOP
  (void)std::strncpy(scopflow->compute_mode,
                     SCOPFLOWOptions::compute_mode.default_value.c_str(),
                     sizeof(scopflow->compute_mode));
  scopflow->verbosity_level = SCOPFLOWOptions::verbosity_level.default_value;
#endif

  scopflow->setupcalled = PETSC_FALSE;
  *scopflowout = scopflow;

  ExaGOLog(EXAGO_LOG_INFO, "{}", "SCOPFLOW: Application created");
  PetscFunctionReturn(0);
}

/**
 * @brief Destroys the security constrained optimal power flow application
 * object
 * @param[in] scopflow the scopflow object to destroy
 */
PetscErrorCode SCOPFLOWDestroy(SCOPFLOW *scopflow) {
  PetscErrorCode ierr;
  PetscInt c;

  PetscFunctionBegin;

  /* Solution vector */
  ierr = VecDestroy(&(*scopflow)->X);
  CHKERRQ(ierr);
  ierr = VecDestroy(&(*scopflow)->gradobj);
  CHKERRQ(ierr);

  /* Lower and upper bounds on X */
  ierr = VecDestroy(&(*scopflow)->Xl);
  CHKERRQ(ierr);
  ierr = VecDestroy(&(*scopflow)->Xu);
  CHKERRQ(ierr);

  /* Constraints vector */
  ierr = VecDestroy(&(*scopflow)->G);
  CHKERRQ(ierr);

  ierr = VecDestroy(&(*scopflow)->Gl);
  CHKERRQ(ierr);
  ierr = VecDestroy(&(*scopflow)->Gu);
  CHKERRQ(ierr);
  ierr = VecDestroy(&(*scopflow)->Lambda);
  CHKERRQ(ierr);

  if ((*scopflow)->comm->size == 1) {
    /* Jacobian of constraints and Hessian */
    ierr = MatDestroy(&(*scopflow)->Jac);
    CHKERRQ(ierr);
    ierr = MatDestroy(&(*scopflow)->Hes);
    CHKERRQ(ierr);
  }

  if ((*scopflow)->solverops.destroy) {
    ierr = ((*scopflow)->solverops.destroy)(*scopflow);
  }

  if ((*scopflow)->modelops.destroy) {
    ierr = ((*scopflow)->modelops.destroy)(*scopflow);
  }

  /* Destroy TCOPFLOW or OPFLOW objects */
  if (!(*scopflow)->ismultiperiod) {
    if ((*scopflow)->opflow0 != NULL) {
      ierr = OPFLOWDestroy(&(*scopflow)->opflow0);
      CHKERRQ(ierr);
    }
    for (c = 0; c < (*scopflow)->nc; c++) {
      ierr = OPFLOWDestroy(&(*scopflow)->opflows[c]);
      CHKERRQ(ierr);
    }
    ierr = PetscFree((*scopflow)->opflows);
    CHKERRQ(ierr);
  } else {
    for (c = 0; c < (*scopflow)->nc; c++) {
      ierr = TCOPFLOWDestroy(&(*scopflow)->tcopflows[c]);
      CHKERRQ(ierr);
    }
    ierr = PetscFree((*scopflow)->tcopflows);
    CHKERRQ(ierr);
  }

  ierr = PetscFree((*scopflow)->xstarti);
  CHKERRQ(ierr);
  ierr = PetscFree((*scopflow)->gstarti);
  CHKERRQ(ierr);
  ierr = PetscFree((*scopflow)->nxi);
  CHKERRQ(ierr);
  ierr = PetscFree((*scopflow)->ngi);
  CHKERRQ(ierr);
  ierr = PetscFree((*scopflow)->nconeqcoup);
  CHKERRQ(ierr);
  ierr = PetscFree((*scopflow)->nconineqcoup);
  CHKERRQ(ierr);
  if ((*scopflow)->Nc > 1) {
    ierr = ContingencyListDestroy(&(*scopflow)->ctgclist);
    CHKERRQ(ierr);
  }

  ierr = COMMDestroy(&(*scopflow)->comm);
  CHKERRQ(ierr);

  ierr = PetscFree(*scopflow);
  CHKERRQ(ierr);
  //  *scopflow = 0;
  PetscFunctionReturn(0);
}

/**
 * @brief Enable/Disable multi-period SCOPLOW
 *
 * @param[in] scopflow the scopflow application object
 * @param[in] ismultperiod PETSC_FALSE for single-period, PETSC_TRUE otherwise
 */
PetscErrorCode SCOPFLOWEnableMultiPeriod(SCOPFLOW scopflow,
                                         PetscBool ismultiperiod) {
  PetscFunctionBegin;
  scopflow->ismultiperiod = ismultiperiod;
  PetscFunctionReturn(0);
}

/**
 * @brief Enable/Disable power imbalance variables
 *
 * @param[in] scopflow the scopflow application object
 * @param[in] use_power_imbalance PETSC_TRUE if using power imbalance variables,
 *                                false otherwise.
 */
PetscErrorCode
SCOPFLOWEnablePowerImbalanceVariables(SCOPFLOW scopflow,
                                      PetscBool use_power_imbalance) {
  PetscFunctionBegin;
  scopflow->enable_powerimbalance_variables = use_power_imbalance;
  PetscFunctionReturn(0);
}

/**
 * @brief Enable/Disable line flow constraints
 *
 * @param[in] scopflow the scopflow application object
 * @param[in] ignore_constraints PETSC_FALSE if constraints,
 *                                PETSC_TRUE otherwise.
 */
PetscErrorCode SCOPFLOWIgnoreLineflowConstraints(SCOPFLOW scopflow,
                                                 PetscBool ignore_constraints) {
  PetscFunctionBegin;
  scopflow->ignore_lineflow_constraints = ignore_constraints;
  PetscFunctionReturn(0);
}

/**
 * @brief Sets the model for SCOPFLOW
 *
 * @param[in] scopflow scopflow application object
 * @param[in] modelname name of the model
 */
PetscErrorCode SCOPFLOWSetModel(SCOPFLOW scopflow, const char *modelname) {
  PetscErrorCode ierr, (*r)(SCOPFLOW) = NULL;
  PetscInt i;
  PetscFunctionBegin;
  PetscBool match;
  for (i = 0; i < scopflow->nmodelsregistered; i++) {
    ierr = PetscStrcmp(scopflow->SCOPFLOWModelList[i].name, modelname, &match);
    CHKERRQ(ierr);
    if (match) {
      r = scopflow->SCOPFLOWModelList[i].create;
      break;
    }
  }

  if (!r)
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE,
            "Unknown type for SCOPFLOW Model %s", modelname);

  /* Null the function pointers */
  scopflow->modelops.destroy = 0;
  scopflow->modelops.setup = 0;
  scopflow->modelops.setnumvariablesandconstraints = 0;
  scopflow->modelops.setvariablebounds = 0;
  scopflow->modelops.setconstraintbounds = 0;
  scopflow->modelops.setvariableandconstraintbounds = 0;
  scopflow->modelops.setinitialguess = 0;
  scopflow->modelops.computeconstraints = 0;
  scopflow->modelops.computejacobian = 0;
  scopflow->modelops.computehessian = 0;
  scopflow->modelops.computeobjandgradient = 0;
  scopflow->modelops.computebaseobjective = 0;
  scopflow->modelops.computetotalobjective = 0;
  scopflow->modelops.computegradient = 0;

  ierr = PetscStrcpy(scopflow->modelname, modelname);
  CHKERRQ(ierr);
  /* Call the underlying implementation constructor */
  ierr = (*r)(scopflow);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/**
 * @brief Sets the solver for SCOPFLOW
 *
 * @param[in] scopflow the scopflow application object
 * @param[in] solvername name of the solver
 */
PetscErrorCode SCOPFLOWSetSolver(SCOPFLOW scopflow, const char *solvername) {
  PetscErrorCode ierr, (*r)(SCOPFLOW) = NULL;
  PetscInt i;
  PetscFunctionBegin;
  PetscBool match;
  for (i = 0; i < scopflow->nsolversregistered; i++) {
    ierr =
        PetscStrcmp(scopflow->SCOPFLOWSolverList[i].name, solvername, &match);
    CHKERRQ(ierr);
    if (match) {
      r = scopflow->SCOPFLOWSolverList[i].create;
      break;
    }
  }

  if (!r)
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE,
            "Unknown type for SCOPFLOW Solver %s", solvername);

  /* Initialize (Null) the function pointers */
  scopflow->solverops.destroy = 0;
  scopflow->solverops.solve = 0;
  scopflow->solverops.setup = 0;

  ierr = PetscStrcpy(scopflow->solvername, solvername);
  CHKERRQ(ierr);
  /* Call the underlying implementation constructor */
  ierr = (*r)(scopflow);
  CHKERRQ(ierr);
  scopflow->solverset = PETSC_TRUE;

  PetscFunctionReturn(0);
}

/**
 * @brief Sets and reads the network data
 *
 * @param[in] scopflow the scopflow object
 * @param[in] netfile the name of the network file
 *
 * @note The input data must be in MATPOWER data format
 */
PetscErrorCode SCOPFLOWSetNetworkData(SCOPFLOW scopflow, const char netfile[]) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscMemcpy(scopflow->netfile, netfile,
                     PETSC_MAX_PATH_LEN * sizeof(char));
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/**
 * @brief Sets the pload data
 *
 * @param[in] scopflow the scopflow object
 * @param[in] profile the name of the pload file
 *
 * @note The input data must be in MATPOWER data format
 */
PetscErrorCode SCOPFLOWSetPLoadData(SCOPFLOW scopflow, const char profile[]) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr =
      PetscMemcpy(scopflow->pload, profile, PETSC_MAX_PATH_LEN * sizeof(char));
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/**
 * @brief Sets the qload data
 *
 * @param[in] scopflow the scopflow object
 * @param[in] profile the name of the qload file
 *
 * @note The input data must be in MATPOWER data format
 */
PetscErrorCode SCOPFLOWSetQLoadData(SCOPFLOW scopflow, const char profile[]) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr =
      PetscMemcpy(scopflow->qload, profile, PETSC_MAX_PATH_LEN * sizeof(char));
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/**
 * @brief Sets the WindGen profile data
 *
 * @param[in] scopflow the scopflow object
 * @param[in] profile the name of the profile
 *
 * @note The input data must be in MATPOWER data format
 */
PetscErrorCode SCOPFLOWSetWindGenProfile(SCOPFLOW scopflow,
                                         const char profile[]) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscMemcpy(scopflow->windgen, profile,
                     PETSC_MAX_PATH_LEN * sizeof(char));
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* This is kind of a hack to update the variable bounds for OPFLOW based on the
 * mode SCOPFLOW uses
 */
PetscErrorCode SCOPFLOWUpdateOPFLOWVariableBounds(OPFLOW opflow, Vec Xl, Vec Xu,
                                                  void *ctx) {
  PetscErrorCode ierr;
  SCOPFLOW scopflow = (SCOPFLOW)ctx;

  PetscFunctionBegin;
  if (opflow->has_gensetpoint) {
    /* Modify the bounds on ramping variables */
    PetscInt j, k;
    PS ps = opflow->ps;
    PSBUS bus;
    PSGEN gen;
    PetscScalar *xl, *xu;

    ierr = VecGetArray(Xl, &xl);
    CHKERRQ(ierr);
    ierr = VecGetArray(Xu, &xu);
    CHKERRQ(ierr);
    for (j = 0; j < ps->nbus; j++) {
      bus = &ps->bus[j];
      for (k = 0; k < bus->ngen; k++) {
        ierr = PSBUSGetGen(bus, k, &gen);
        CHKERRQ(ierr);
        if (!gen->status || gen->isrenewable)
          continue;
        if (scopflow->mode == 0) {
          /* Only ref. bus responsible for make-up power for contingencies */
          if (bus->ide != REF_BUS) {
            xl[opflow->idxn2sd_map[gen->startxpdevloc]] =
                xu[opflow->idxn2sd_map[gen->startxpdevloc]] = 0.0;
          } else {
            // Reference bus can supply full output
            xl[opflow->idxn2sd_map[gen->startxpdevloc]] = gen->pb - gen->pt;
            xu[opflow->idxn2sd_map[gen->startxpdevloc]] = gen->pt - gen->pb;
          }
        } else {
          xl[opflow->idxn2sd_map[gen->startxpdevloc]] = -gen->ramp_rate_30min;
          xu[opflow->idxn2sd_map[gen->startxpdevloc]] = gen->ramp_rate_30min;
        }
      }
    }
    ierr = VecRestoreArray(Xl, &xl);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(Xu, &xu);
    CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
/*
  SCOPFLOWSetUp - Sets up an security constrained optimal power flow application
object

  Input Parameters:
. scopflow - the SCOPFLOW object

  Notes:
  This routine sets up the SCOPFLOW object and the underlying PS object. It
  also distributes the PS object when used in parallel.
*/
PetscErrorCode SCOPFLOWSetUp(SCOPFLOW scopflow) {
  PetscErrorCode ierr;
  PetscInt c, i, j;
  PS ps;
  OPFLOW opflow;

  char ploadprofile[PETSC_MAX_PATH_LEN];
  char qloadprofile[PETSC_MAX_PATH_LEN];
  char windgenprofile[PETSC_MAX_PATH_LEN];
  PetscBool flg1 = PETSC_FALSE, flg2 = PETSC_FALSE, flg3 = PETSC_FALSE;

  char scopflowsolvername[max_solver_name_len];
  char scopflowmodelname[max_model_name_len];
  PetscBool scopflowsolverset = PETSC_FALSE;

  (void)std::strncpy(scopflowmodelname,
                     SCOPFLOWOptions::model.default_value.c_str(),
                     sizeof(scopflowmodelname));
  (void)std::strncpy(scopflowsolvername,
                     SCOPFLOWOptions::solver.default_value.c_str(),
                     sizeof(scopflowsolvername));

  PetscFunctionBegin;

  PetscOptionsBegin(scopflow->comm->type, NULL, "SCOPFLOW options", NULL);

  ierr = PetscOptionsString(
      SCOPFLOWOptions::model.opt.c_str(), SCOPFLOWOptions::model.desc.c_str(),
      "", scopflowmodelname, scopflowmodelname, max_model_name_len, NULL);
  CHKERRQ(ierr);
  ierr = PetscOptionsString(SCOPFLOWOptions::solver.opt.c_str(),
                            SCOPFLOWOptions::solver.desc.c_str(), "",
                            scopflowsolvername, scopflowsolvername,
                            max_solver_name_len, &scopflowsolverset);
  CHKERRQ(ierr);
  ierr =
      PetscOptionsString(SCOPFLOWOptions::subproblem_model.opt.c_str(),
                         SCOPFLOWOptions::subproblem_model.desc.c_str(), "",
                         scopflow->subproblem_model, scopflow->subproblem_model,
                         max_model_name_len, NULL);
  CHKERRQ(ierr);
  ierr = PetscOptionsString(SCOPFLOWOptions::subproblem_solver.opt.c_str(),
                            SCOPFLOWOptions::subproblem_solver.desc.c_str(), "",
                            scopflow->subproblem_solver,
                            scopflow->subproblem_solver, max_solver_name_len,
                            NULL);
  CHKERRQ(ierr);
#ifdef EXAGO_ENABLE_HIOP
  ierr = PetscOptionsString(SCOPFLOWOptions::compute_mode.opt.c_str(),
                            SCOPFLOWOptions::compute_mode.desc.c_str(), "",
                            scopflow->compute_mode, scopflow->compute_mode,
                            max_model_name_len, NULL);
  CHKERRQ(ierr);
  ierr = PetscOptionsInt(SCOPFLOWOptions::verbosity_level.opt.c_str(),
                         SCOPFLOWOptions::verbosity_level.desc.c_str(), "",
                         scopflow->verbosity_level, &scopflow->verbosity_level,
                         NULL);
  CHKERRQ(ierr);
#endif

  ierr = PetscOptionsBool(SCOPFLOWOptions::iscoupling.opt.c_str(),
                          SCOPFLOWOptions::iscoupling.desc.c_str(), "",
                          scopflow->iscoupling, &scopflow->iscoupling, NULL);
  CHKERRQ(ierr);
  ierr = PetscOptionsInt(SCOPFLOWOptions::Nc.opt.c_str(),
                         SCOPFLOWOptions::Nc.desc.c_str(), "", scopflow->Nc,
                         &scopflow->Nc, NULL);
  CHKERRQ(ierr);
  ierr = PetscOptionsInt(SCOPFLOWOptions::mode.opt.c_str(),
                         SCOPFLOWOptions::mode.desc.c_str(), "", scopflow->mode,
                         &scopflow->mode, NULL);
  CHKERRQ(ierr);
  ierr =
      PetscOptionsBool(SCOPFLOWOptions::enable_multiperiod.opt.c_str(),
                       SCOPFLOWOptions::enable_multiperiod.desc.c_str(), "",
                       scopflow->ismultiperiod, &scopflow->ismultiperiod, NULL);
  CHKERRQ(ierr);
  if (scopflow->ismultiperiod) {
    /* Set loadp,loadq, and windgen files */
    ierr = PetscOptionsString(SCOPFLOWOptions::ploadprofile.opt.c_str(),
                              SCOPFLOWOptions::ploadprofile.desc.c_str(), "",
                              ploadprofile, ploadprofile, PETSC_MAX_PATH_LEN,
                              &flg1);
    CHKERRQ(ierr);
    ierr = PetscOptionsString(SCOPFLOWOptions::qloadprofile.opt.c_str(),
                              SCOPFLOWOptions::qloadprofile.desc.c_str(), "",
                              qloadprofile, qloadprofile, PETSC_MAX_PATH_LEN,
                              &flg2);
    CHKERRQ(ierr);
    ierr = PetscOptionsString(SCOPFLOWOptions::windgenprofile.opt.c_str(),
                              SCOPFLOWOptions::windgenprofile.desc.c_str(), "",
                              windgenprofile, windgenprofile,
                              PETSC_MAX_PATH_LEN, &flg3);
    CHKERRQ(ierr);
    ierr = PetscOptionsReal(SCOPFLOWOptions::dT.opt.c_str(),
                            SCOPFLOWOptions::dT.desc.c_str(), "", scopflow->dT,
                            &scopflow->dT, NULL);
    CHKERRQ(ierr);
    ierr = PetscOptionsReal(SCOPFLOWOptions::duration.opt.c_str(),
                            SCOPFLOWOptions::duration.desc.c_str(), "",
                            scopflow->duration, &scopflow->duration, NULL);
    CHKERRQ(ierr);
  }
  ierr = PetscOptionsReal(SCOPFLOWOptions::tolerance.opt.c_str(),
                          SCOPFLOWOptions::tolerance.desc.c_str(), "",
                          scopflow->tolerance, &scopflow->tolerance, NULL);
  CHKERRQ(ierr);
  PetscOptionsEnd();

  if (scopflow->ctgcfileset) {
    if (scopflow->Nc < 0)
      scopflow->Nc = MAX_CONTINGENCIES;
    else
      scopflow->Nc += 1;

    if (scopflow->Nc > 1) {
      /* Create contingency list object */
      ierr = ContingencyListCreate(scopflow->Nc, &scopflow->ctgclist);
      CHKERRQ(ierr);
      ierr = ContingencyListSetData(
          scopflow->ctgclist, scopflow->ctgcfileformat, scopflow->ctgcfile);
      CHKERRQ(ierr);
      ierr = ContingencyListReadData(scopflow->ctgclist, &scopflow->Nc);
      CHKERRQ(ierr);
      scopflow->Nc += 1;
    }
  } else {
    if (scopflow->Nc == -1)
      scopflow->Nc = 1;
  }

  int q = scopflow->Nc / scopflow->comm->size;
  int d = scopflow->Nc % scopflow->comm->size;
  if (d) {
    scopflow->nc = q + ((scopflow->comm->rank < d) ? 1 : 0);
  } else {
    scopflow->nc = q;
  }
  ierr = MPI_Scan(&scopflow->nc, &scopflow->cend, 1, MPIU_INT, MPI_SUM,
                  scopflow->comm->type);
  CHKERRQ(ierr);
  scopflow->cstart = scopflow->cend - scopflow->nc;

  ExaGOLog(EXAGO_LOG_INFO,
           "SCOPFLOW running with {:d} subproblems (base case + {:d} "
           "contingencies)",
           scopflow->nc, scopflow->nc - 1);
  //  ExaGOLogUseEveryRank(PETSC_TRUE);
  //  ExaGOLog(EXAGO_LOG_INFO,"Rank %d has %d contingencies, range [%d --
  //  %d]\n",scopflow->comm->rank,scopflow->nc,scopflow->cstart,scopflow->cend);
  //  ExaGOLogUseEveryRank(PETSC_FALSE);

  /* Set model */
  if (!scopflow->ismultiperiod) {
    ierr = SCOPFLOWSetModel(scopflow, scopflowmodelname);
    CHKERRQ(ierr);
  } else {
    ierr = SCOPFLOWSetModel(scopflow, "GENRAMPT");
    CHKERRQ(ierr);
  }

  /* Set solver */
  if (scopflowsolverset) {
    if (scopflow->solver)
      ierr = (*scopflow->solverops.destroy)(scopflow);
    ierr = SCOPFLOWSetSolver(scopflow, scopflowsolvername);
    CHKERRQ(ierr);
    ExaGOLog(EXAGO_LOG_INFO, "SCOPFLOW: Using {} solver", scopflowsolvername);
  } else {
    if (!scopflow->solver) {
      ierr = SCOPFLOWSetSolver(scopflow, SCOPFLOWSOLVER_IPOPT);
      CHKERRQ(ierr);
      ExaGOLog(EXAGO_LOG_INFO, "SCOPFLOW: Using {} solver",
               SCOPFLOWSOLVER_IPOPT);
    }
  }

  if (!scopflow->ismultiperiod) {
    /* Create OPFLOW objects */
    ierr = PetscCalloc1(scopflow->nc, &scopflow->opflows);
    CHKERRQ(ierr);

    /* Create base-case OPFLOW */
    ierr = OPFLOWCreate(PETSC_COMM_SELF, &scopflow->opflow0);
    CHKERRQ(ierr);
    ierr = OPFLOWHasBusPowerImbalance(
        scopflow->opflow0, scopflow->enable_powerimbalance_variables);
    CHKERRQ(ierr);
    ierr = OPFLOWIgnoreLineflowConstraints(
        scopflow->opflow0, scopflow->ignore_lineflow_constraints);
    CHKERRQ(ierr);
    /* Base-case problem model should be POWER_BALANCE_POLAR */
    ierr = OPFLOWSetModel(scopflow->opflow0, OPFLOWMODEL_PBPOL);
    CHKERRQ(ierr);
    /* Base-case problem solver should be IPOPT */
    ierr = OPFLOWSetSolver(scopflow->opflow0, OPFLOWSOLVER_IPOPT);
    CHKERRQ(ierr);
    ierr = OPFLOWReadMatPowerData(scopflow->opflow0, scopflow->netfile);
    CHKERRQ(ierr);
    ierr = OPFLOWSetInitializationType(scopflow->opflow0, scopflow->type);
    CHKERRQ(ierr);
    ierr = OPFLOWSetGenBusVoltageType(scopflow->opflow0,
                                      scopflow->genbusvoltagetype);
#ifdef EXAGO_ENABLE_HIOP
    ierr = OPFLOWSetHIOPComputeMode(scopflow->opflow0, scopflow->compute_mode);
    CHKERRQ(ierr);
    ierr = OPFLOWSetHIOPVerbosityLevel(scopflow->opflow0,
                                       scopflow->verbosity_level);
    CHKERRQ(ierr);
#endif
    CHKERRQ(ierr);

    ierr = PSSetUp(scopflow->opflow0->ps);
    CHKERRQ(ierr);
    if (scopflow->scen) {
      // Scenario set by SOPFLOW, apply it */
      ierr = PSApplyScenario(scopflow->opflow0->ps, *scopflow->scen);
      CHKERRQ(ierr);
    }
    ierr = OPFLOWSetUp(scopflow->opflow0);
    CHKERRQ(ierr);

    ierr = OPFLOWSetObjectiveType(scopflow->opflow0, MIN_GEN_COST);
    CHKERRQ(ierr);

    for (c = 0; c < scopflow->nc; c++) {
      ierr = OPFLOWCreate(PETSC_COMM_SELF, &scopflow->opflows[c]);
      CHKERRQ(ierr);
      ierr = OPFLOWHasBusPowerImbalance(
          scopflow->opflows[c], scopflow->enable_powerimbalance_variables);
      CHKERRQ(ierr);
      ierr = OPFLOWIgnoreLineflowConstraints(
          scopflow->opflows[c], scopflow->ignore_lineflow_constraints);
      CHKERRQ(ierr);
      ierr = OPFLOWSetModel(scopflow->opflows[c], scopflow->subproblem_model);
      CHKERRQ(ierr);
      ierr = OPFLOWSetSolver(scopflow->opflows[c], scopflow->subproblem_solver);
      CHKERRQ(ierr);
      ierr = OPFLOWSetInitializationType(scopflow->opflows[c], scopflow->type);
      CHKERRQ(ierr);
      ierr = OPFLOWSetGenBusVoltageType(scopflow->opflows[c],
                                        scopflow->genbusvoltagetype);
      CHKERRQ(ierr);
#ifdef EXAGO_ENABLE_HIOP
      ierr = OPFLOWSetHIOPComputeMode(scopflow->opflows[c],
                                      scopflow->compute_mode);
      CHKERRQ(ierr);
      ierr = OPFLOWSetHIOPVerbosityLevel(scopflow->opflows[c],
                                         scopflow->verbosity_level);
      CHKERRQ(ierr);
#endif

      ierr = OPFLOWReadMatPowerData(scopflow->opflows[c], scopflow->netfile);
      CHKERRQ(ierr);
      /* Set up the PS object for opflow */
      ps = scopflow->opflows[c]->ps;
      ierr = PSSetUp(ps);
      CHKERRQ(ierr);

      if (scopflow->scen) {
        // Scenario set by SOPFLOW, apply it */
        ierr = PSApplyScenario(ps, *scopflow->scen);
        CHKERRQ(ierr);
      }

      /* Set contingencies */
      if (scopflow->ctgcfileset && scopflow->Nc > 1) {
        Contingency ctgc = scopflow->ctgclist->cont[scopflow->cstart + c];
        for (j = 0; j < ctgc.noutages; j++) {
          if (ctgc.outagelist[j].type == GEN_OUTAGE) {
            PetscInt gbus = ctgc.outagelist[j].bus;
            char *gid = ctgc.outagelist[j].id;
            PetscInt status = ctgc.outagelist[j].status;
            ierr = PSSetGenStatus(ps, gbus, gid, status);
            CHKERRQ(ierr);
          }
          if (ctgc.outagelist[j].type == BR_OUTAGE) {
            PetscInt fbus = ctgc.outagelist[j].fbus;
            PetscInt tbus = ctgc.outagelist[j].tbus;
            char *brid = ctgc.outagelist[j].id;
            PetscInt status = ctgc.outagelist[j].status;
            ierr = PSSetLineStatus(ps, fbus, tbus, brid, status);
            CHKERRQ(ierr);
          }
        }
      }

      if (scopflow->cstart + c == 0) { /* First stage */
        ierr = OPFLOWSetModel(scopflow->opflows[c], OPFLOWMODEL_PBPOL);
        CHKERRQ(ierr);
        ierr = OPFLOWSetSolver(scopflow->opflows[c], OPFLOWSOLVER_IPOPT);
        CHKERRQ(ierr);
        ierr = OPFLOWSetObjectiveType(scopflow->opflows[c], MIN_GEN_COST);
        CHKERRQ(ierr);
      } else { /* Second stages */
        ierr = OPFLOWHasGenSetPoint(scopflow->opflows[c], PETSC_TRUE);
        CHKERRQ(ierr); /* Activates ramping variables */
        ierr = OPFLOWSetModel(scopflow->opflows[c], OPFLOWMODEL_PBPOL);
        CHKERRQ(ierr);
        ierr = OPFLOWSetSolver(scopflow->opflows[c], OPFLOWSOLVER_IPOPT);
        CHKERRQ(ierr);
        ierr = OPFLOWSetObjectiveType(scopflow->opflows[c], NO_OBJ);
        CHKERRQ(ierr);
        ierr = OPFLOWSetUpdateVariableBoundsFunction(
            scopflow->opflows[c], SCOPFLOWUpdateOPFLOWVariableBounds,
            (void *)scopflow);
      }

      ierr = OPFLOWSetUp(scopflow->opflows[c]);
      CHKERRQ(ierr);
    }
  } else {
    TCOPFLOW tcopflow;

    /* Create TCOPFLOW objects */
    ierr = PetscCalloc1(scopflow->nc, &scopflow->tcopflows);
    CHKERRQ(ierr);
    for (c = 0; c < scopflow->nc; c++) {
      ierr = TCOPFLOWCreate(PETSC_COMM_SELF, &scopflow->tcopflows[c]);
      CHKERRQ(ierr);
      tcopflow = scopflow->tcopflows[c];

      ierr = TCOPFLOWSetNetworkData(tcopflow, scopflow->netfile);
      CHKERRQ(ierr);
      // CLI option overrides any existing setting in scopflow
      if (flg1) {
        ierr = SCOPFLOWSetPLoadData(scopflow, ploadprofile);
        CHKERRQ(ierr);
      }
      if (flg2) {
        ierr = SCOPFLOWSetQLoadData(scopflow, qloadprofile);
        CHKERRQ(ierr);
      }
      if (flg3) {
        ierr = SCOPFLOWSetWindGenProfile(scopflow, windgenprofile);
        CHKERRQ(ierr);
      }

      ierr =
          TCOPFLOWSetLoadProfiles(tcopflow, scopflow->pload, scopflow->qload);
      CHKERRQ(ierr);
      ierr = TCOPFLOWSetWindGenProfiles(tcopflow, scopflow->windgen);
      CHKERRQ(ierr);
      ierr = TCOPFLOWSetTimeStepandDuration(tcopflow, scopflow->dT,
                                            scopflow->duration);
      CHKERRQ(ierr);

      if (scopflow->scen) {
        /* SOPFLOW has set scenario, pass it to TCOPFLOW */
        ierr = TCOPFLOWSetScenario(tcopflow, scopflow->scen);
        CHKERRQ(ierr);
      }

      /* Set contingencies */
      if (scopflow->ctgcfileset) {
        Contingency ctgc = scopflow->ctgclist->cont[scopflow->cstart + c];
        // Set this contingency with TCOPFLOW
        ierr = TCOPFLOWSetContingency(tcopflow, &ctgc);
        CHKERRQ(ierr);
      }

      ierr = TCOPFLOWSetUp(tcopflow);
      CHKERRQ(ierr);
    }
  }

  ierr = PetscCalloc1(scopflow->nc, &scopflow->nconeqcoup);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->nc, &scopflow->nconineqcoup);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->nc, &scopflow->nxi);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->nc, &scopflow->ngi);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->nc, &scopflow->xstarti);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->nc, &scopflow->gstarti);
  CHKERRQ(ierr);

  /* Set number of variables and constraints */
  ierr = (*scopflow->modelops.setnumvariablesandconstraints)(
      scopflow, scopflow->nxi, scopflow->ngi, scopflow->nconeqcoup,
      scopflow->nconineqcoup);

  if (!scopflow->ismultiperiod) {
    scopflow->nx = scopflow->nxi[0];
    scopflow->ncon = scopflow->ngi[0];
    scopflow->nconcoup = scopflow->nconeqcoup[0] + scopflow->nconineqcoup[0];
    opflow = scopflow->opflows[0];
    scopflow->nconeq = opflow->nconeq;
    scopflow->nconineq = opflow->nconineq;

    for (i = 1; i < scopflow->nc; i++) {
      scopflow->xstarti[i] = scopflow->xstarti[i - 1] + scopflow->nxi[i - 1];
      scopflow->gstarti[i] = scopflow->gstarti[i - 1] + scopflow->ngi[i - 1];
      scopflow->nx += scopflow->nxi[i];
      scopflow->ncon += scopflow->ngi[i];
      scopflow->nconcoup += scopflow->nconeqcoup[i] + scopflow->nconineqcoup[i];
      opflow = scopflow->opflows[i];
      scopflow->nconeq += opflow->nconeq;
      scopflow->nconineq += opflow->nconineq;
    }
  } else {
    TCOPFLOW tcopflow;
    scopflow->nx = scopflow->nxi[0];
    scopflow->ncon = scopflow->ngi[0];
    scopflow->nconcoup = scopflow->nconineqcoup[0];
    tcopflow = scopflow->tcopflows[0];
    scopflow->nconeq = tcopflow->Nconeq;
    scopflow->nconineq = tcopflow->Nconineq + tcopflow->Nconcoup;

    for (i = 1; i < scopflow->nc; i++) {
      scopflow->xstarti[i] = scopflow->xstarti[i - 1] + scopflow->nxi[i - 1];
      scopflow->gstarti[i] = scopflow->gstarti[i - 1] + scopflow->ngi[i - 1];
      scopflow->nx += scopflow->nxi[i];
      scopflow->ncon += scopflow->ngi[i];
      scopflow->nconcoup += scopflow->nconineqcoup[i];
      tcopflow = scopflow->tcopflows[i];
      scopflow->nconeq += tcopflow->Nconeq;
      scopflow->nconineq += tcopflow->Nconineq + tcopflow->Nconcoup;
    }
  }

  /* Create vector X */
  ierr = VecCreate(scopflow->comm->type, &scopflow->X);
  CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->X, scopflow->nx, PETSC_DECIDE);
  CHKERRQ(ierr);
  ierr = VecSetFromOptions(scopflow->X);
  CHKERRQ(ierr);
  ierr = VecGetSize(scopflow->X, &scopflow->Nx);
  CHKERRQ(ierr);

  ierr = VecDuplicate(scopflow->X, &scopflow->Xl);
  CHKERRQ(ierr);
  ierr = VecDuplicate(scopflow->X, &scopflow->Xu);
  CHKERRQ(ierr);
  ierr = VecDuplicate(scopflow->X, &scopflow->gradobj);
  CHKERRQ(ierr);

  /* Vector for constraints */
  ierr = VecCreate(scopflow->comm->type, &scopflow->G);
  CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->G, scopflow->ncon, PETSC_DECIDE);
  CHKERRQ(ierr);

  ierr = VecSetFromOptions(scopflow->G);
  CHKERRQ(ierr);

  ierr = VecGetSize(scopflow->G, &scopflow->Ncon);
  CHKERRQ(ierr);

  /* Constraint bounds vectors  */
  ierr = VecDuplicate(scopflow->G, &scopflow->Gl);
  CHKERRQ(ierr);
  ierr = VecDuplicate(scopflow->G, &scopflow->Gu);
  CHKERRQ(ierr);

  /* The matrices are not used in parallel, so we don't need to create them */
  if (scopflow->comm->size == 1) {
    /* Constraint Jacobian */
    ierr = MatCreate(scopflow->comm->type, &scopflow->Jac);
    CHKERRQ(ierr);
    ierr = MatSetSizes(scopflow->Jac, scopflow->ncon, scopflow->nx,
                       scopflow->Ncon, scopflow->Nx);
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(scopflow->Jac);
    CHKERRQ(ierr);

    /* Assume 10% sparsity */
    ierr = MatSeqAIJSetPreallocation(scopflow->Jac,
                                     (PetscInt)(0.1 * scopflow->nx), NULL);
    CHKERRQ(ierr);

    ierr = MatSetOption(scopflow->Jac, MAT_NEW_NONZERO_ALLOCATION_ERR,
                        PETSC_FALSE);
    CHKERRQ(ierr);

    /* Hessian */
    ierr = MatCreate(scopflow->comm->type, &scopflow->Hes);
    CHKERRQ(ierr);
    ierr = MatSetSizes(scopflow->Hes, scopflow->nx, scopflow->nx, scopflow->Nx,
                       scopflow->Nx);
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(scopflow->Hes);
    CHKERRQ(ierr);
    /* Assume 10% sparsity */
    ierr = MatSeqAIJSetPreallocation(scopflow->Hes,
                                     (PetscInt)(0.1 * scopflow->nx), NULL);
    CHKERRQ(ierr);
    ierr = MatSetOption(scopflow->Hes, MAT_NEW_NONZERO_ALLOCATION_ERR,
                        PETSC_FALSE);
    CHKERRQ(ierr);
  }

  /* Lagrangian multipliers */
  ierr = VecDuplicate(scopflow->G, &scopflow->Lambda);
  CHKERRQ(ierr);

  ierr = (*scopflow->solverops.setup)(scopflow);
  CHKERRQ(ierr);
  ExaGOLog(EXAGO_LOG_INFO, "{}", "SCOPFLOW: Setup completed");

  scopflow->setupcalled = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/**
 * @brief Solves the AC security constrained optimal power flow
 *
 * @param[in] scopflow the security constrained optimal power flow application
 * object
 */
PetscErrorCode SCOPFLOWSolve(SCOPFLOW scopflow) {
  PetscErrorCode ierr;
  PetscBool issolver_ipopt;

  PetscFunctionBegin;

  if (!scopflow->setupcalled) {
    ierr = SCOPFLOWSetUp(scopflow);
  }

  ierr = PetscStrcmp(scopflow->solvername, "IPOPT", &issolver_ipopt);
  CHKERRQ(ierr);

  if (issolver_ipopt) { /* Don't need to do this if solver is not ipopt */
    /* Set bounds on variables */
    if (scopflow->modelops.setvariablebounds) {
      ierr = (*scopflow->modelops.setvariablebounds)(scopflow, scopflow->Xl,
                                                     scopflow->Xu);
      CHKERRQ(ierr);
    }

    /* Set bounds on constraints */
    if (scopflow->modelops.setconstraintbounds) {
      ierr = (*scopflow->modelops.setconstraintbounds)(scopflow, scopflow->Gl,
                                                       scopflow->Gu);
      CHKERRQ(ierr);
    }

    /* Set initial guess */
    if (scopflow->modelops.setinitialguess) {
      ierr = (*scopflow->modelops.setinitialguess)(scopflow, scopflow->X);
      CHKERRQ(ierr);
    }

    ierr = VecSet(scopflow->Lambda, 1.0);
    CHKERRQ(ierr);
  }

  /* Solve */
  ierr = (*scopflow->solverops.solve)(scopflow);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/**
 * @brief Returns the objective function value for the base case
 *
 * @param[in] scopflow the scopflow object
 * @param[in] objbase the objective function value for the base case
 *
 * @note Should be called after the optimization finishes
 */
PetscErrorCode SCOPFLOWGetBaseObjective(SCOPFLOW scopflow, PetscReal *objbase) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*scopflow->solverops.getbaseobjective)(scopflow, objbase);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/**
 * @brief Returns the total objective function value
 *
 * @param[in] scopflow the scopflow object
 * @param[in] objtot the total objective function value
 *
 * @note Should be called after the optimization finishes
 */
PetscErrorCode SCOPFLOWGetTotalObjective(SCOPFLOW scopflow, PetscReal *objtot) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*scopflow->solverops.gettotalobjective)(scopflow, objtot);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/**
 * @brief Returns the SCOPFLOW solution for a given contingency
 *
 * @param[in] scopflow the scopflow object
 * @param[in] contnum Contingency number (0 for base/no-contingency)
 * @param[in] X the scopflow solution for the given contingency
 *
 * @note Should be called after the optimization finishes
 */
PetscErrorCode SCOPFLOWGetSolution(SCOPFLOW scopflow, PetscInt contnum,
                                   Vec *X) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*scopflow->solverops.getsolution)(scopflow, contnum, X);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/**
 * @brief Returns the SCOPFLOW constraints for a given contingency
 *
 * @param[in] scopflow the scopflow object
 * @param[in] contnum contingency number (0 for base/no-contingency)
 * @param[in] G the scopflow constraints
 *
 * @note Should be called after the optimization finishes.
 * Equality constraints first followed by inequality constraints
 */
PetscErrorCode SCOPFLOWGetConstraints(SCOPFLOW scopflow, PetscInt contnum,
                                      Vec *G) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*scopflow->solverops.getconstraints)(scopflow, contnum, G);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/**
 * @brief Returns the SCOPFLOW constraint multipliers for a given contingency
 *
 * @param[in] scopflow the scopflow object
 * @param[in] contnum  Contingency number (0 for base/no-contingency)
 * @param[in] Lambda pointer to vector object
 *
 * @note Should be called after the optimization finishes.
 * Equality constraint multipliers first followed by inequality constraint
 * multipliers
 */
PetscErrorCode SCOPFLOWGetConstraintMultipliers(SCOPFLOW scopflow,
                                                PetscInt contnum, Vec *Lambda) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*scopflow->solverops.getconstraintmultipliers)(scopflow, contnum,
                                                         Lambda);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/**
 * @brief Check SCOPFLOW convergence
 *
 * @param[in] scopflow the scopflow object
 * @param[in] status PETSC_TRUE if converged, PETSC_FALSE otherwise
 *
 * @note Should be called after the optimization finishes
 */
PetscErrorCode SCOPFLOWGetConvergenceStatus(SCOPFLOW scopflow,
                                            PetscBool *status) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*scopflow->solverops.getconvergencestatus)(scopflow, status);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/**
 * @brief Gets the number of variables and constraints
 *
 * @param[in] scopflow the scopflow object
 * @param[in] nx number of variables
 * @param[in] ncon number of constraints
 *
 * @note Must be called after SCOPFLOWSetUp
 */
PetscErrorCode SCOPFLOWGetNumVariablesandConstraints(SCOPFLOW scopflow,
                                                     PetscInt *Nx,
                                                     PetscInt *Ncon) {
  PetscFunctionBegin;
  *Nx = scopflow->nx;
  *Ncon = scopflow->ncon;
  PetscFunctionReturn(0);
}

/**
 * @brief Get the number of iterations for a given solver
 *
 * @param[in] scopflow application object
 * @param [in] iter pointer to solver iteration variable
 */
PetscErrorCode SCOPFLOWGetNumIterations(SCOPFLOW scopflow, PetscInt *iter) {
  PetscErrorCode ierr;
  *iter = scopflow->numiter;
  PetscFunctionReturn(0);
}

/**
 * @brief Set the time step for SCOPFLOW
 *
 * @param[in] scopflow application object
 * @param[in] dT solver timestep
 */
PetscErrorCode SCOPFLOWSetTimeStep(SCOPFLOW scopflow, PetscReal dT) {
  PetscFunctionBegin;
  scopflow->dT = dT;
  PetscFunctionReturn(0);
}

/**
 * @brief Set the problem duration for SCOPFLOW
 *
 * @param[in] scopflow application object
 * @param[in] duration application duration
 */
PetscErrorCode SCOPFLOWSetDuration(SCOPFLOW scopflow, PetscReal duration) {
  PetscFunctionBegin;
  scopflow->duration = duration;
  PetscFunctionReturn(0);
}

/**
 * @brief Set the time-step and problem duration for SCOPFLOW
 *
 * @param[in] scopflow application object
 * @param[in] dT time-step
 * @param[in] duration problem duration
 */
PetscErrorCode SCOPFLOWSetTimeStepandDuration(SCOPFLOW scopflow, PetscReal dT,
                                              PetscReal duration) {
  PetscFunctionBegin;
  scopflow->dT = dT;
  scopflow->duration = duration;
  PetscFunctionReturn(0);
}

/**
 * @brief Set the solver tolerance for SCOPFLOW
 *
 * @param[in] scopflow application object
 * @param[in] tol solver tolerance
 */
PetscErrorCode SCOPFLOWSetTolerance(SCOPFLOW scopflow, PetscReal tol) {
  PetscFunctionBegin;
  scopflow->tolerance = tol;
  PetscFunctionReturn(0);
}

/**
 * @brief Get the solver tolerance for SCOPFLOW
 *
 * @param[in] scopflow application object
 * @param[in] tol pointer to solver tolerance variable
 */
PetscErrorCode SCOPFLOWGetTolerance(SCOPFLOW scopflow, PetscReal *tol) {
  PetscFunctionBegin;
  *tol = scopflow->tolerance;
  PetscFunctionReturn(0);
}

/**
 * @brief Set SCOPFLOW initialization type
 *
 * @param[in] scopflow application object
 * @param[in] type initialization type for underlying OPFLOW structs
 */
PetscErrorCode SCOPFLOWSetInitilizationType(SCOPFLOW scopflow,
                                            OPFLOWInitializationType type) {
  PetscFunctionBegin;
  scopflow->type = type;
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSetScenario(SCOPFLOW scopflow, Scenario *scen) {
  PetscFunctionBegin;
  scopflow->scen = scen;
  PetscFunctionReturn(0);
}

/**
 * @brief Set SCOPFLOW number of contingencies
 *
 * @param[in] scopflow application object
 * @param[in] num_cont number of contingencies
 */
PetscErrorCode SCOPFLOWSetNumContingencies(SCOPFLOW scopflow,
                                           PetscInt num_cont) {
  PetscFunctionBegin;
  scopflow->Nc = num_cont;
  PetscFunctionReturn(0);
}

/**
 * @brief Set SCOPFLOW genbusvoltage type
 *
 * @param[in] scopflow application object
 * @param[in] type genbusvoltage type for underlying OPFLOW structs
 */
PetscErrorCode SCOPFLOWSetGenBusVoltageType(SCOPFLOW scopflow,
                                            OPFLOWGenBusVoltageType type) {
  PetscFunctionBegin;
  scopflow->genbusvoltagetype = type;
  PetscFunctionReturn(0);
}

/**
 * @brief Set SCOPFLOW mode
 *
 * @param[in] scopflow application object
 * @param[in] mode for scopflow. 0 = preventive, 1 = corrective
 */
PetscErrorCode SCOPFLOWSetMode(SCOPFLOW scopflow, PetscInt mode) {
  PetscFunctionBegin;
  scopflow->mode = mode;
  PetscFunctionReturn(0);
}

/*
  SCOPFLOWSetLoadProfiles - Sets the data files for time-varying load profiles

  Input Parameter
+  scopflow - The SCOPFLOW object
.  ploadproffile - The name of the file with real power load variationt
-  qloadproffile - The name of the file with reactive power load variationt
*/
PetscErrorCode SCOPFLOWSetLoadProfiles(SCOPFLOW scopflow,
                                       const char ploadprofile[],
                                       const char qloadprofile[]) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = SCOPFLOWSetPLoadData(scopflow, ploadprofile);
  CHKERRQ(ierr);
  ierr = SCOPFLOWSetQLoadData(scopflow, qloadprofile);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
 * @brief Set the contingency data file for SCOPFLOW
 *
 * @param[in] scopflow application object
 * @param[in] contingency file format
 * @param[in] name of the contingency list file
 */
PetscErrorCode
SCOPFLOWSetContingencyData(SCOPFLOW scopflow,
                           ContingencyFileInputFormat ctgcfileformat,
                           const char ctgcfile[]) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscMemcpy(scopflow->ctgcfile, ctgcfile,
                     PETSC_MAX_PATH_LEN * sizeof(char));
  CHKERRQ(ierr);
  scopflow->ctgcfileformat = ctgcfileformat;
  scopflow->ctgcfileset = PETSC_TRUE;

  PetscFunctionReturn(0);
}

/*
  SCOPFLOWSetSubproblemModel - Set subproblem model

  Inputs
+ scopflow - scopflow object
- model    - model name

  Options
  -scopflow_subproblem_model

  Notes: This is used with HIOP solver only
*/
PetscErrorCode SCOPFLOWSetSubproblemModel(SCOPFLOW scopflow,
                                          const char modelname[]) {
  PetscErrorCode ierr;
  ierr = PetscStrcpy(scopflow->subproblem_model, modelname);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SCOPFLOWSetSubproblemSolver - Set subproblem solver

  Inputs
+ scopflow - scopflow object
- solver   - solver name

  Option Name
  -scopflow_subproblem_solver
  Notes: This is used with HIOP solver only
*/
PetscErrorCode SCOPFLOWSetSubproblemSolver(SCOPFLOW scopflow,
                                           const char solvername[]) {
  PetscErrorCode ierr;
  ierr = PetscStrcpy(scopflow->subproblem_solver, solvername);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SCOPFLOWSetSubproblemComputeMode - Set subproblem solver compute mode

  Inputs
+ scopflow - scopflow object
- mode   - subproblem solver compute mode

  Option Name
  -hiop_compute_mode
  Notes: This is used with HIOP solver only
*/
PetscErrorCode SCOPFLOWSetComputeMode(SCOPFLOW scopflow, const char mode[]) {
  PetscErrorCode ierr;
  ierr = PetscStrcpy(scopflow->compute_mode, mode);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SOPFLOWSetSubproblemVerbosityLevel - Set subproblem solver verbosity level

  Inputs
+ sopflow - sopflow object
- level   - subproblem solver verbosity level

  Option Name
  -hiop_verbosity_level
  Notes: This is used with HIOP solver only
*/
PetscErrorCode SCOPFLOWSetVerbosityLevel(SCOPFLOW scopflow, int level) {
  PetscErrorCode ierr;
  scopflow->verbosity_level = level;
  PetscFunctionReturn(0);
}
