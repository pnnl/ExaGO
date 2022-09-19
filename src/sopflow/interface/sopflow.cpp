#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>
#include <private/sopflowimpl.h>
#include <utils.h>

/*
  SOPFLOWCreate - Creates a stochastic optimal power flow application object

  Input Parameters
. mpicomm - The MPI communicator

  Output Parameters
. sopflowout - The stochastic optimal power flow application object
*/
PetscErrorCode SOPFLOWCreate(MPI_Comm mpicomm, SOPFLOW *sopflowout) {
  PetscErrorCode ierr;
  SOPFLOW sopflow;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1, &sopflow);
  CHKERRQ(ierr);

  ierr = COMMCreate(mpicomm, &sopflow->comm);
  CHKERRQ(ierr);

  sopflow->Nconeq = sopflow->nconeq = 0;
  sopflow->Nconineq = sopflow->nconineq = 0;
  sopflow->Ncon = sopflow->ncon = 0;
  sopflow->Nx = sopflow->nx = 0;
  sopflow->Ns = sopflow->ns = SOPFLOWOptions::Ns.default_value;
  sopflow->Gi = NULL;
  sopflow->Lambdai = NULL;

  sopflow->obj_factor = 1.0;
  sopflow->objbase = sopflow->objtot = 0.0;
  sopflow->tolerance = SOPFLOWOptions::tolerance.default_value;

  sopflow->solver = NULL;
  sopflow->model = NULL;

  /* Default subproblemmodel and solver */
  (void)std::strncpy(sopflow->subproblem_model,
                     SOPFLOWOptions::subproblem_model.default_value.c_str(),
                     sizeof(sopflow->subproblem_model));
  (void)std::strncpy(sopflow->subproblem_solver,
                     SOPFLOWOptions::subproblem_model.default_value.c_str(),
                     sizeof(sopflow->subproblem_solver));

  sopflow->mode = SOPFLOWOptions::mode.default_value;

  sopflow->ismulticontingency =
      SOPFLOWOptions::enable_multicontingency.default_value;
  sopflow->Nc = SOPFLOWOptions::Nc.default_value;
  sopflow->ismultiperiod = PETSC_FALSE;

  sopflow->flatten_contingencies =
      SOPFLOWOptions::flatten_contingencies.default_value;
  sopflow->ctgclist = NULL;
  sopflow->ctgcfileset = PETSC_FALSE;

  sopflow->nmodelsregistered = 0;
  sopflow->SOPFLOWModelRegisterAllCalled = PETSC_FALSE;

  /* Register all models */
  ierr = SOPFLOWModelRegisterAll(sopflow);

  sopflow->nsolversregistered = 0;
  sopflow->SOPFLOWSolverRegisterAllCalled = PETSC_FALSE;

  /* Register all solvers */
  ierr = SOPFLOWSolverRegisterAll(sopflow);

  /* Run-time options */
  sopflow->iscoupling = SOPFLOWOptions::iscoupling.default_value;

#if defined(EXAGO_ENABLE_HIOP)
  /* HIOP options */
  strcpy(sopflow->compute_mode, "auto");
  sopflow->verbosity_level = 0;
#endif

  sopflow->scenfileset = PETSC_FALSE;
  sopflow->scenunctype = NONE;
  sopflow->setupcalled = PETSC_FALSE;
  *sopflowout = sopflow;

  ExaGOLog(EXAGO_LOG_INFO, "{}", "SOPFLOW: Application created");
  PetscFunctionReturn(0);
}

/*
  SOPFLOWDestroy - Destroys the stochastic optimal power flow application object

  Input Parameter
. sopflow - The SOPFLOW object to destroy
*/
PetscErrorCode SOPFLOWDestroy(SOPFLOW *sopflow) {
  PetscErrorCode ierr;
  PetscInt s, i;

  PetscFunctionBegin;

  /* Solution vector */
  ierr = VecDestroy(&(*sopflow)->X);
  CHKERRQ(ierr);
  ierr = VecDestroy(&(*sopflow)->gradobj);
  CHKERRQ(ierr);

  /* Lower and upper bounds on X */
  ierr = VecDestroy(&(*sopflow)->Xl);
  CHKERRQ(ierr);
  ierr = VecDestroy(&(*sopflow)->Xu);
  CHKERRQ(ierr);

  /* Constraints vector */
  ierr = VecDestroy(&(*sopflow)->G);
  CHKERRQ(ierr);

  ierr = VecDestroy(&(*sopflow)->Gl);
  CHKERRQ(ierr);
  ierr = VecDestroy(&(*sopflow)->Gu);
  CHKERRQ(ierr);
  ierr = VecDestroy(&(*sopflow)->Lambda);
  CHKERRQ(ierr);

  /* Jacobian of constraints and Hessian */
  if ((*sopflow)->comm->size == 1) {
    ierr = MatDestroy(&(*sopflow)->Jac);
    CHKERRQ(ierr);
    ierr = MatDestroy(&(*sopflow)->Hes);
    CHKERRQ(ierr);
  }

  if ((*sopflow)->solverops.destroy) {
    ierr = ((*sopflow)->solverops.destroy)(*sopflow);
  }

  if ((*sopflow)->modelops.destroy) {
    ierr = ((*sopflow)->modelops.destroy)(*sopflow);
  }

  /* Destroy SCOPFLOW or OPFLOW objects */
  if ((*sopflow)->ismulticontingency) {
    for (s = 0; s < (*sopflow)->ns; s++) {
      ierr = SCOPFLOWDestroy(&(*sopflow)->scopflows[s]);
      CHKERRQ(ierr);
    }
    ierr = PetscFree((*sopflow)->scopflows);
    CHKERRQ(ierr);
  } else {
    ierr = OPFLOWDestroy(&(*sopflow)->opflow0);
    CHKERRQ(ierr);
    for (s = 0; s < (*sopflow)->ns; s++) {
      ierr = OPFLOWDestroy(&(*sopflow)->opflows[s]);
      CHKERRQ(ierr);
    }
    ierr = PetscFree((*sopflow)->opflows);
    CHKERRQ(ierr);

    ierr = PetscFree((*sopflow)->scen_num);
    CHKERRQ(ierr);
    if ((*sopflow)->flatten_contingencies) {
      if ((*sopflow)->Nc > 1 && (*sopflow)->ctgclist) {
        ierr = ContingencyListDestroy(&(*sopflow)->ctgclist);
        CHKERRQ(ierr);
      }
      ierr = PetscFree((*sopflow)->cont_num);
      CHKERRQ(ierr);
    }
  }

  ierr = PetscFree((*sopflow)->xstarti);
  CHKERRQ(ierr);
  ierr = PetscFree((*sopflow)->gstarti);
  CHKERRQ(ierr);
  ierr = PetscFree((*sopflow)->nxi);
  CHKERRQ(ierr);
  ierr = PetscFree((*sopflow)->ngi);
  CHKERRQ(ierr);
  ierr = PetscFree((*sopflow)->nconeqcoup);
  CHKERRQ(ierr);
  ierr = PetscFree((*sopflow)->nconineqcoup);
  CHKERRQ(ierr);

  /* Destroy scenario list */
  for (s = 0; s < (*sopflow)->Ns; s++) {
    for (i = 0; i < (*sopflow)->scenlist.scen[s].nforecast; i++) {
      ierr = PetscFree((*sopflow)->scenlist.scen[s].forecastlist[i].buses);
      CHKERRQ(ierr);
      for (int j = 0; j < (*sopflow)->scenlist.scen[s].forecastlist[i].nele;
           j++) {
        ierr = PetscFree((*sopflow)->scenlist.scen[s].forecastlist[i].id[j]);
        CHKERRQ(ierr);
      }
      ierr = PetscFree((*sopflow)->scenlist.scen[s].forecastlist[i].id);
      CHKERRQ(ierr);
      ierr = PetscFree((*sopflow)->scenlist.scen[s].forecastlist[i].val);
      CHKERRQ(ierr);
    }
  }

  ierr = PetscFree((*sopflow)->scenlist.scen);
  CHKERRQ(ierr);

  MPI_Comm_free(&(*sopflow)->subcomm);

  ierr = COMMDestroy(&(*sopflow)->comm);
  CHKERRQ(ierr);

  ierr = PetscFree(*sopflow);
  CHKERRQ(ierr);
  //  *sopflow = 0;
  PetscFunctionReturn(0);
}

/*
  SOPFLOWEnableMultiContingency - Enable/Disable multi-contingency SOPLOW

  Input Parameters:
+ sopflow - sopflow application object
- ismulticontingency - PETSC_FALSE for no contingencies, PETSC_TRUE otherwise
*/
PetscErrorCode SOPFLOWEnableMultiContingency(SOPFLOW sopflow,
                                             PetscBool ismulticontingency) {
  PetscFunctionBegin;
  sopflow->ismulticontingency = ismulticontingency;
  PetscFunctionReturn(0);
}

/*
  SOPFLOWFlattenContingencies - Flatten or not flatten contingencies

  Input Parameters:
+ sopflow - sopflow application object
- flatten - PETSC_FALSE if not flattening contingencies, PETSC_TRUE otherwise
*/
PetscErrorCode SOPFLOWFlattenContingencies(SOPFLOW sopflow, PetscBool flatten) {
  PetscFunctionBegin;
  sopflow->flatten_contingencies = flatten;
  PetscFunctionReturn(0);
}

/*
  SOPFLOWSetModel - Sets the model for SOPFLOW

  Input Parameters:
+ sopflow - opflow application object
- modelname - name of the model
*/
PetscErrorCode SOPFLOWSetModel(SOPFLOW sopflow, const char *modelname) {
  PetscErrorCode ierr, (*r)(SOPFLOW) = NULL;
  PetscInt i;
  PetscFunctionBegin;
  PetscBool match;
  for (i = 0; i < sopflow->nmodelsregistered; i++) {
    ierr = PetscStrcmp(sopflow->SOPFLOWModelList[i].name, modelname, &match);
    CHKERRQ(ierr);
    if (match) {
      r = sopflow->SOPFLOWModelList[i].create;
      break;
    }
  }

  if (!r)
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE,
             "Unknown type for SOPFLOW Model %s", modelname);

  /* Null the function pointers */
  sopflow->modelops.destroy = 0;
  sopflow->modelops.setup = 0;
  sopflow->modelops.setnumvariablesandconstraints = 0;
  sopflow->modelops.setvariablebounds = 0;
  sopflow->modelops.setconstraintbounds = 0;
  sopflow->modelops.setvariableandconstraintbounds = 0;
  sopflow->modelops.setinitialguess = 0;
  sopflow->modelops.computeconstraints = 0;
  sopflow->modelops.computejacobian = 0;
  sopflow->modelops.computehessian = 0;
  sopflow->modelops.computeobjandgradient = 0;
  sopflow->modelops.computebaseobjective = 0;
  sopflow->modelops.computetotalobjective = 0;
  sopflow->modelops.computegradient = 0;

  ierr = PetscStrcpy(sopflow->modelname, modelname);
  CHKERRQ(ierr);
  /* Call the underlying implementation constructor */
  ierr = (*r)(sopflow);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
   SOPFLOWSetTolerance - Set the solver tolerance

  Input Parameters:
+ sopflow - sopflow application object
- tol    - solver tolerance
*/
PetscErrorCode SOPFLOWSetTolerance(SOPFLOW sopflow, PetscReal tol) {
  PetscFunctionBegin;
  sopflow->tolerance = tol;
  PetscFunctionReturn(0);
}

/**
 * SOPFLOWGetTolerance - Get the solver tolerance
 * Input Parameters:
 * sopflow - sopflow application object
 * tol    - pointer to solver tolerance variable that will be set
 */
PetscErrorCode SOPFLOWGetTolerance(SOPFLOW sopflow, PetscReal *tol) {
  PetscFunctionBegin;
  *tol = sopflow->tolerance;
  PetscFunctionReturn(0);
}

/*
   SOPFLOWSetInitializationType - Set initialization type for underlying OPFLOW

  Input Parameters:
+ sopflow - sopflow application object
- type    - OPFLOW initialization type
*/
PetscErrorCode SOPFLOWSetInitializationType(SOPFLOW sopflow,
                                            OPFLOWInitializationType type) {
  PetscFunctionBegin;
  sopflow->initialization_type = type;
  PetscFunctionReturn(0);
}

/*
   SOPFLOWSetIgnoreLineflowConstraints - Set ignore lineflow constraints flag
for underlying OPFLOW

  Input Parameters:
+ sopflow - sopflow application object
- flag    - OPFLOW ignore lineflow constraints flag
*/
PetscErrorCode SOPFLOWSetIgnoreLineflowConstraints(SOPFLOW sopflow,
                                                   PetscBool flag) {
  PetscFunctionBegin;
  sopflow->ignore_lineflow_constraints = flag;
  PetscFunctionReturn(0);
}

/*
   SOPFLOWSetGenBusVoltageType - Set gen bus voltage type for underlying OPFLOW

  Input Parameters:
+ sopflow - sopflow application object
- type    - OPFLOW gen bus voltage type
*/
PetscErrorCode SOPFLOWSetGenBusVoltageType(SOPFLOW sopflow,
                                           OPFLOWGenBusVoltageType type) {
  PetscFunctionBegin;
  sopflow->gen_bus_voltage_type = type;
  PetscFunctionReturn(0);
}

/*
  SOPFLOWSetSolver - Sets the solver for SOPFLOW

  Input Parameters:
+ sopflow - sopflow application object
- solvername - name of the solver
*/
PetscErrorCode SOPFLOWSetSolver(SOPFLOW sopflow, const char *solvername) {
  PetscErrorCode ierr, (*r)(SOPFLOW) = NULL;
  PetscInt i;
  PetscFunctionBegin;
  PetscBool match;
  for (i = 0; i < sopflow->nsolversregistered; i++) {
    ierr = PetscStrcmp(sopflow->SOPFLOWSolverList[i].name, solvername, &match);
    CHKERRQ(ierr);
    if (match) {
      r = sopflow->SOPFLOWSolverList[i].create;
      break;
    }
  }

  if (!r)
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE,
             "Unknown type for SOPFLOW Solver %s", solvername);

  /* Initialize (Null) the function pointers */
  sopflow->solverops.destroy = 0;
  sopflow->solverops.solve = 0;
  sopflow->solverops.setup = 0;

  ierr = PetscStrcpy(sopflow->solvername, solvername);
  CHKERRQ(ierr);
  /* Call the underlying implementation constructor */
  ierr = (*r)(sopflow);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/**
 * @brief Set SOPFLOW number of scenarios
 *
 * @param[in] sopflow application object
 * @param[in] num_scen number of scenarios
 */
PetscErrorCode SOPFLOWSetNumScenarios(SOPFLOW sopflow, PetscInt num_scen) {
  PetscFunctionBegin;
  sopflow->Ns = num_scen;
  PetscFunctionReturn(0);
}

/*
  SOPFLOWSetNetworkData - Sets and reads the network data

  Input Parameter
+  sopflow - The SOPFLOW object
-  netfile - The name of the network file

  Notes: The input data must be in MATPOWER data format
*/
PetscErrorCode SOPFLOWSetNetworkData(SOPFLOW sopflow, const char netfile[]) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr =
      PetscMemcpy(sopflow->netfile, netfile, PETSC_MAX_PATH_LEN * sizeof(char));
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  SOPFLOWSetContingencyData - Sets and reads the contingency data

  Input Parameter
+  sopflow - The SOPFLOW object
-  ctgfile - The name of the ctgfile file
*/
PetscErrorCode SOPFLOWSetContingencyData(SOPFLOW sopflow,
                                         const char ctgfile[]) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr =
      PetscMemcpy(sopflow->ctgfile, ctgfile, PETSC_MAX_PATH_LEN * sizeof(char));
  CHKERRQ(ierr);
  sopflow->ctgcfileset = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/*
  SOPFLOWSetWindGenProfile - Sets and reads the windgen data

  Input Parameter
+  sopflow - The SOPFLOW object
-  windgen - The name of the network file
*/
PetscErrorCode SOPFLOWSetWindGenProfile(SOPFLOW sopflow, const char windgen[]) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr =
      PetscMemcpy(sopflow->windgen, windgen, PETSC_MAX_PATH_LEN * sizeof(char));
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* This is kind of a hack to update the variable bounds for OPFLOW based on the
 * mode SOPFLOW uses
 */
PetscErrorCode SOPFLOWUpdateOPFLOWVariableBounds(OPFLOW opflow, Vec Xl, Vec Xu,
                                                 void *ctx) {
  PetscErrorCode ierr;
  SOPFLOW sopflow = (SOPFLOW)ctx;

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
        if (sopflow->mode == 0) {
          /* Only ref. bus responsible for make-up power for contingencies */
          if (bus->ide == REF_BUS) {
            // Ref. bus can supply full output
            xl[opflow->idxn2sd_map[gen->startxpdevloc]] = gen->pb - gen->pt;
            xu[opflow->idxn2sd_map[gen->startxpdevloc]] = gen->pt - gen->pb;
          } else {
            xl[opflow->idxn2sd_map[gen->startxpdevloc]] =
                xu[opflow->idxn2sd_map[gen->startxpdevloc]] = 0.0;
          }
        } else {
          xl[opflow->idxn2sd_map[gen->startxpdevloc]] =
              -gen->ramp_rate_30min; // gen->pb - gen->pt
          xu[opflow->idxn2sd_map[gen->startxpdevloc]] =
              gen->ramp_rate_30min; // gen->pt - gen->pb
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

/**
 * @brief Set the contingency data file for SOPFLOW
 *
 * @param[in] sopflow application object
 * @param[in] contingency file format
 * @param[in] name of the contingency list file
 */
PetscErrorCode
SOPFLOWSetContingencyData(SOPFLOW sopflow,
                          ContingencyFileInputFormat ctgcfileformat,
                          const char ctgcfile[]) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscMemcpy(sopflow->ctgcfile, ctgcfile,
                     PETSC_MAX_PATH_LEN * sizeof(char));
  CHKERRQ(ierr);
  sopflow->ctgcfileformat = ctgcfileformat;
  sopflow->ctgcfileset = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/*
  SOPFLOWSetUp - Sets up an stochastic optimal power flow application object

  Input Parameters:
. sopflow - the SOPFLOW object

  Notes:
  This routine sets up the SOPFLOW object and the underlying PS object. It
  also distributes the PS object when used in parallel.
*/
PetscErrorCode SOPFLOWSetUp(SOPFLOW sopflow) {
  PetscErrorCode ierr;
  PetscBool sopflowsolverset;
  PetscInt c, s, i, j;
  PS ps;
  OPFLOW opflow;
  char ctgcfile[PETSC_MAX_PATH_LEN];
  char windgen[PETSC_MAX_PATH_LEN];
  PetscBool flgctgc = PETSC_FALSE;
  PetscBool flgwindgen = PETSC_FALSE;
  PetscBool issopflowsolverhiop;

  char opflowmodelname[max_model_name_len];
  char sopflowsolvername[max_solver_name_len];
  char sopflowmodelname[max_model_name_len];

  (void)std::strncpy(sopflowmodelname,
                     SOPFLOWOptions::sopflow_model.default_value.c_str(),
                     sizeof(sopflowmodelname));
  (void)std::strncpy(sopflowsolvername,
                     SOPFLOWOptions::sopflow_solver.default_value.c_str(),
                     sizeof(sopflowsolvername));
  (void)std::strncpy(opflowmodelname,
                     SOPFLOWOptions::opflow_model.default_value.c_str(),
                     sizeof(opflowmodelname));

  PetscFunctionBegin;

  ierr = PetscOptionsBegin(sopflow->comm->type, NULL, "SOPFLOW options", NULL);
  CHKERRQ(ierr);
  ierr = PetscOptionsString(SOPFLOWOptions::sopflow_model.opt.c_str(),
                            SOPFLOWOptions::sopflow_model.desc.c_str(), "",
                            sopflowmodelname, sopflowmodelname,
                            max_model_name_len, NULL);
  CHKERRQ(ierr);
  ierr = PetscOptionsString(SOPFLOWOptions::sopflow_solver.opt.c_str(),
                            SOPFLOWOptions::sopflow_solver.desc.c_str(), "",
                            sopflowsolvername, sopflowsolvername,
                            max_solver_name_len, &sopflowsolverset);
  CHKERRQ(ierr);
  ierr =
      PetscOptionsString(SOPFLOWOptions::subproblem_model.opt.c_str(),
                         SOPFLOWOptions::subproblem_model.desc.c_str(), "",
                         sopflow->subproblem_model, sopflow->subproblem_model,
                         max_model_name_len, NULL);
  CHKERRQ(ierr);
  ierr =
      PetscOptionsString(SOPFLOWOptions::subproblem_solver.opt.c_str(),
                         SOPFLOWOptions::subproblem_solver.desc.c_str(), "",
                         sopflow->subproblem_solver, sopflow->subproblem_solver,
                         max_solver_name_len, NULL);
  CHKERRQ(ierr);

  ierr = PetscOptionsBool(SOPFLOWOptions::iscoupling.opt.c_str(),
                          SOPFLOWOptions::iscoupling.desc.c_str(), "",
                          sopflow->iscoupling, &sopflow->iscoupling, NULL);
  CHKERRQ(ierr);
  ierr = PetscOptionsInt(SOPFLOWOptions::Ns.opt.c_str(),
                         SOPFLOWOptions::Ns.desc.c_str(), "", sopflow->Ns,
                         &sopflow->Ns, NULL);
  CHKERRQ(ierr);
  ierr = PetscOptionsInt(SOPFLOWOptions::mode.opt.c_str(),
                         SOPFLOWOptions::mode.desc.c_str(), "", sopflow->mode,
                         &sopflow->mode, NULL);
  CHKERRQ(ierr);

  ierr = PetscOptionsBool(SOPFLOWOptions::enable_multicontingency.opt.c_str(),
                          SOPFLOWOptions::enable_multicontingency.desc.c_str(),
                          "", sopflow->ismulticontingency,
                          &sopflow->ismulticontingency, NULL);
  CHKERRQ(ierr);

  ierr = PetscOptionsBool(SOPFLOWOptions::flatten_contingencies.opt.c_str(),
                          SOPFLOWOptions::flatten_contingencies.desc.c_str(),
                          "", sopflow->flatten_contingencies,
                          &sopflow->flatten_contingencies, NULL);
  CHKERRQ(ierr);

  ierr = PetscOptionsString(SOPFLOWOptions::ctgcfile.opt.c_str(),
                            SOPFLOWOptions::ctgcfile.desc.c_str(), "", ctgcfile,
                            ctgcfile, PETSC_MAX_PATH_LEN, &flgctgc);
  CHKERRQ(ierr);
  ierr = PetscOptionsInt(SOPFLOWOptions::Nc.opt.c_str(),
                         SOPFLOWOptions::Nc.desc.c_str(), "", sopflow->Nc,
                         &sopflow->Nc, NULL);
  CHKERRQ(ierr);
  ierr = PetscOptionsString(SOPFLOWOptions::windgen.opt.c_str(),
                            SOPFLOWOptions::windgen.desc.c_str(), "", windgen,
                            windgen, PETSC_MAX_PATH_LEN, &flgwindgen);
  CHKERRQ(ierr);

  ierr = PetscOptionsReal(SOPFLOWOptions::tolerance.opt.c_str(),
                          SOPFLOWOptions::tolerance.desc.c_str(), "",
                          sopflow->tolerance, &sopflow->tolerance, NULL);
  CHKERRQ(ierr);
  PetscOptionsEnd();

  if (sopflow->Ns == 0)
    SETERRQ(PETSC_COMM_SELF, 0, "Number of scenarios should be greater than 0");

  if (sopflow->scenfileset) {
    if (sopflow->Ns == -1) {
      ierr = SOPFLOWGetNumScenarios(sopflow, sopflow->scenfileformat,
                                    sopflow->scenfile, &sopflow->Ns);
      CHKERRQ(ierr);
    }

    ierr = PetscCalloc1(sopflow->Ns, &sopflow->scenlist.scen);
    CHKERRQ(ierr);

    for (s = 0; s < sopflow->Ns; s++)
      sopflow->scenlist.scen->nforecast = 0;
    ierr = SOPFLOWReadScenarioData(sopflow, sopflow->scenfileformat,
                                   sopflow->scenfile);
    CHKERRQ(ierr);
    sopflow->Ns = sopflow->scenlist.Nscen;
  } else {
    if (sopflow->Ns == -1)
      sopflow->Ns = 1;
  }

  if (sopflow->ismulticontingency && sopflow->flatten_contingencies) {
    sopflow->ismulticontingency =
        PETSC_FALSE; // Collapse contingencies into scenarios

    if (flgctgc && !sopflow->ctgcfileset) {
      /* Need to remove hard coded native format later */
      ierr = SOPFLOWSetContingencyData(sopflow, NATIVE, ctgcfile);
      CHKERRQ(ierr);
    }

    if (sopflow->Nc < 0)
      sopflow->Nc = MAX_CONTINGENCIES;
    else
      sopflow->Nc += 1;

    if (sopflow->Nc > 1) {
      /* Create contingency list object */
      ierr = ContingencyListCreate(sopflow->Nc, &sopflow->ctgclist);
      CHKERRQ(ierr);
      ierr = ContingencyListSetData(sopflow->ctgclist, sopflow->ctgcfileformat,
                                    sopflow->ctgcfile);
      CHKERRQ(ierr);
      ierr = ContingencyListReadData(sopflow->ctgclist, &sopflow->Nc);
      CHKERRQ(ierr);
      sopflow->Nc += 1;
    } else {
      if (sopflow->Nc == -1)
        sopflow->Nc = 1;
    }
  } else {
    sopflow->Nc = 1;
  }

  if (sopflow->ismulticontingency) {
    /* Calculate the sizes of the subcommunicator */
    /* Here, we compute a group of ranks that will work on each
       scenario. This grouping of ranks is useful so that the
       underlying scopflow operates in parallel using the
       rank group for its scenario. As an example, assume
       there 4 ranks (sopflow->comm->size) with 2 scenaarios (sopflow->Ns)
       then we create 2 subcommunicators (one for each scenario) withh
       the first subcommunicator using ranks 0 and 1 and the second
       one using ranks 2 and 3. We use MPI_Comm_split*() to split
       the outer communicator (sopflow->comm->type) into subcommunicators
       (sopflow->subcomm)
    */
    int *color, *ns, *sstart, *send;
    color = (int *)malloc(sopflow->comm->size * sizeof(int));
    ns = (int *)malloc(sopflow->comm->size * sizeof(int));
    sstart = (int *)malloc(sopflow->comm->size * sizeof(int));
    send = (int *)malloc(sopflow->comm->size * sizeof(int));

    for (i = 0; i < sopflow->comm->size; i++) {
      if (sopflow->comm->size > sopflow->Ns) {
        ns[i] = 1;
        if (sopflow->comm->size % sopflow->Ns == 0)
          color[i] = i / (sopflow->comm->size / sopflow->Ns);
        else {
          color[i] =
              PetscMin(i / (sopflow->comm->size / sopflow->Ns),
                       sopflow->Ns); /* This is not an optimal distribution. It
                                        can be further optimized */
        }
      } else {
        ns[i] = sopflow->Ns / sopflow->comm->size;
        if (i >= (sopflow->comm->size - sopflow->Ns % sopflow->comm->size))
          ns[i] += 1;
        color[i] = i;
      }
    }
    sopflow->ns = ns[sopflow->comm->rank];

    sstart[0] = 0;
    send[0] = sstart[0] + ns[0];
    int color_i = color[0];
    i = 0;

    /*
        if(!sopflow->comm->rank) {
        ierr = PetscPrintf(PETSC_COMM_SELF,"Rank[%d]: color = %d, ns = %d,
       sstart= %d, send=%d\n",i,color[i],ns[i],sstart[i],send[i]);CHKERRQ(ierr);
        }
    */

    for (i = 1; i < sopflow->comm->size; i++) {
      if (color[i] == color_i) { /* Same group - copy start and end */
        sstart[i] = sstart[i - 1];
        send[i] = send[i - 1];
      } else { /* New group - update sstart and send */
        sstart[i] = send[i - 1];
        send[i] = sstart[i] + ns[i];
        color_i = color[i];
      }

      /*
        if(!sopflow->comm->rank) {
        ierr = PetscPrintf(PETSC_COMM_SELF,"Rank[%d]: color = %d, ns = %d,
        sstart= %d,
        send=%d\n",i,color[i],ns[i],sstart[i],send[i]);CHKERRQ(ierr);
        }
      */
    }

    sopflow->sstart = sstart[sopflow->comm->rank];
    sopflow->send = send[sopflow->comm->rank];

    MPI_Barrier(sopflow->comm->type);
    ExaGOLog(EXAGO_LOG_INFO, "SOPFLOW running with %d scenarios\n",
             sopflow->Ns);
    //  ExaGOLogUseEveryRank(PETSC_TRUE);
    //  ExaGOLog(EXAGO_LOG_INFO,"Rank %d scenario range [%d --
    //  %d]\n",sopflow->comm->rank,sopflow->sstart,sopflow->send);
    //  ExaGOLogUseEveryRank(PETSC_FALSE);

    /* Create subcommunicators to manage scopflows */
    MPI_Comm_split(sopflow->comm->type, color[sopflow->comm->rank],
                   sopflow->comm->rank, &sopflow->subcomm);

    free(color);
    free(sstart);
    free(send);
    free(ns);
  } else {
    int q = (sopflow->Ns * sopflow->Nc) / sopflow->comm->size;
    int d = (sopflow->Ns * sopflow->Nc) % sopflow->comm->size;

    if (d) {
      sopflow->ns = q + ((sopflow->comm->rank < d) ? 1 : 0);
    } else {
      sopflow->ns = q;
    }
    ierr = MPI_Scan(&sopflow->ns, &sopflow->send, 1, MPIU_INT, MPI_SUM,
                    sopflow->comm->type);
    CHKERRQ(ierr);
    sopflow->sstart = sopflow->send - sopflow->ns;

    ierr = PetscCalloc1(sopflow->ns, &sopflow->scen_num);
    CHKERRQ(ierr);
    if (sopflow->flatten_contingencies) {
      ierr = PetscCalloc1(sopflow->ns, &sopflow->cont_num);
      CHKERRQ(ierr);

      for (s = 0; s < sopflow->ns; s++) {
        sopflow->scen_num[s] = (sopflow->sstart + s) / sopflow->Nc;
        sopflow->cont_num[s] = (sopflow->sstart + s) % sopflow->Nc;
        ierr = PetscPrintf(PETSC_COMM_SELF,
                           "Rank[%d],s = %d, scen_num = %d, cont_num = %d\n",
                           sopflow->comm->rank, sopflow->sstart + s,
                           sopflow->scen_num[s], sopflow->cont_num[s]);
      }
    } else {
      for (s = 0; s < sopflow->ns; s++) {
        sopflow->scen_num[s] = sopflow->sstart + s;

        ierr = PetscPrintf(PETSC_COMM_SELF, "Rank[%d],s = %d, scen_num = %d\n",
                           sopflow->comm->rank, sopflow->sstart + s,
                           sopflow->scen_num[s]);
      }
    }
  }

  /* Set Model */
  if (!sopflow->ismulticontingency) {
    ierr = SOPFLOWSetModel(sopflow, sopflowmodelname);
    CHKERRQ(ierr);
  } else {
    ierr = SOPFLOWSetModel(sopflow, "GENRAMPC");
    CHKERRQ(ierr);
  }

  /* Set solver */
  if (sopflowsolverset) {
    if (sopflow->solver)
      ierr = (*sopflow->solverops.destroy)(sopflow);
    ierr = SOPFLOWSetSolver(sopflow, sopflowsolvername);
    CHKERRQ(ierr);
    ExaGOLog(EXAGO_LOG_INFO, "SOPFLOW: Using {} solver\n", sopflowsolvername);
  } else {
    if (!sopflow->solver) {
      ierr = SOPFLOWSetSolver(sopflow, SOPFLOWSOLVER_IPOPT);
      CHKERRQ(ierr);
      ExaGOLog(EXAGO_LOG_INFO, "SOPFLOW: Using {} solver\n",
               SOPFLOWSOLVER_IPOPT);
    }
  }

  ierr = PetscStrcmp(sopflow->solvername, "HIOP", &issopflowsolverhiop);
  CHKERRQ(ierr);

  if (!sopflow->ismulticontingency) {

    /* Create base-case OPFLOW, only used when solver is HIOP */
    ierr = OPFLOWCreate(PETSC_COMM_SELF, &sopflow->opflow0);
    CHKERRQ(ierr);
    ierr = OPFLOWSetModel(sopflow->opflow0, OPFLOWMODEL_PBPOL);
    CHKERRQ(ierr);
    ierr = OPFLOWSetInitializationType(sopflow->opflow0,
                                       sopflow->initialization_type);
    CHKERRQ(ierr);
    ierr = OPFLOWReadMatPowerData(sopflow->opflow0, sopflow->netfile);
    CHKERRQ(ierr);
    ierr = PSSetUp(sopflow->opflow0->ps);
    CHKERRQ(ierr);
    if (sopflow->scenfileset) {
      ierr = PSApplyScenario(sopflow->opflow0->ps, sopflow->scenlist.scen[0]);
      CHKERRQ(ierr);
    }
    ierr = OPFLOWHasGenSetPoint(sopflow->opflow0, PETSC_TRUE);
    CHKERRQ(ierr);

    ierr = OPFLOWSetUp(sopflow->opflow0);
    CHKERRQ(ierr);

    /* Create OPFLOW objects */
    ierr = PetscCalloc1(sopflow->ns, &sopflow->opflows);
    CHKERRQ(ierr);
    for (s = 0; s < sopflow->ns; s++) {
      ierr = OPFLOWCreate(PETSC_COMM_SELF, &sopflow->opflows[s]);
      CHKERRQ(ierr);
      ierr = OPFLOWSetModel(sopflow->opflows[s], OPFLOWMODEL_PBPOL);
      CHKERRQ(ierr);
      //    ierr =
      //    OPFLOWSetSolver(sopflow->opflows[c],opflowsolvername);CHKERRQ(ierr);

      ierr = OPFLOWSetInitializationType(sopflow->opflows[s],
                                         sopflow->initialization_type);
      CHKERRQ(ierr);
      ierr = OPFLOWSetGenBusVoltageType(sopflow->opflows[s],
                                        sopflow->gen_bus_voltage_type);
      CHKERRQ(ierr);
      ierr = OPFLOWIgnoreLineflowConstraints(
          sopflow->opflows[s], sopflow->ignore_lineflow_constraints);
      CHKERRQ(ierr);
      ierr = OPFLOWReadMatPowerData(sopflow->opflows[s], sopflow->netfile);
      CHKERRQ(ierr);
      /* Set up the PS object for opflow */
      ps = sopflow->opflows[s]->ps;
      ierr = PSSetUp(ps);
      CHKERRQ(ierr);

      if (sopflow->scenfileset) {
        /* Apply scenario */
        ierr =
            PSApplyScenario(ps, sopflow->scenlist.scen[sopflow->scen_num[s]]);
        CHKERRQ(ierr);
        /* Set scenario probability */
        ierr =
            OPFLOWSetWeight(sopflow->opflows[s],
                            sopflow->scenlist.scen[sopflow->scen_num[s]].prob);
      }

      if (flgctgc && sopflow->flatten_contingencies && sopflow->Nc > 1) {
        ierr = PSApplyContingency(
            ps, sopflow->ctgclist->cont[sopflow->cont_num[s]]);
        CHKERRQ(ierr);
      }

      ierr = OPFLOWHasGenSetPoint(sopflow->opflows[s], PETSC_TRUE);
      CHKERRQ(ierr);
      if (sopflow->sstart + s != 0) {
        ierr = OPFLOWSetUpdateVariableBoundsFunction(
            sopflow->opflows[s], SOPFLOWUpdateOPFLOWVariableBounds,
            (void *)sopflow);
      }

      if (flgctgc && sopflow->flatten_contingencies) {
        if (sopflow->cont_num[s] != 0) {
          //	  ierr =
          // OPFLOWSetObjectiveType(sopflow->opflows[s],NO_OBJ);CHKERRQ(ierr);
        }
      }

      /* Set subproblem parameters */
      if (issopflowsolverhiop) {
        ierr = OPFLOWSetModel(sopflow->opflows[s], sopflow->subproblem_model);
        CHKERRQ(ierr);
        ierr = OPFLOWSetSolver(sopflow->opflows[s], sopflow->subproblem_solver);
        CHKERRQ(ierr);
        ierr = OPFLOWSetHIOPComputeMode(sopflow->opflows[s],
                                        sopflow->compute_mode);
        CHKERRQ(ierr);
        ierr = OPFLOWSetHIOPVerbosityLevel(sopflow->opflows[s],
                                           sopflow->verbosity_level);
        CHKERRQ(ierr);
      } 

      ierr = OPFLOWSetUp(sopflow->opflows[s]);
      CHKERRQ(ierr);
      
    }

  } else {

    /* Create SCOPFLOW objects */
    ierr = PetscCalloc1(sopflow->ns, &sopflow->scopflows);
    CHKERRQ(ierr);

    /* Override windgen and ctgc files if CLI option set */
    if (flgctgc) {
      ierr = SOPFLOWSetContingencyData(sopflow, ctgcfile);
      CHKERRQ(ierr);
    }
    if (flgwindgen) {
      ierr = SOPFLOWSetWindGenProfile(sopflow, windgen);
      CHKERRQ(ierr);
    }

    for (s = 0; s < sopflow->ns; s++) {
      /* Create SCOPFLOW object */
      ierr = SCOPFLOWCreate(sopflow->subcomm, &sopflow->scopflows[s]);
      CHKERRQ(ierr);
      /* Set Network data */
      ierr = SCOPFLOWSetNetworkData(sopflow->scopflows[s], sopflow->netfile);
      CHKERRQ(ierr);
      /* Set contingency data */
      /* Should not set hard coded native format */

      ierr = SCOPFLOWSetContingencyData(sopflow->scopflows[s], NATIVE,
                                        sopflow->ctgfile);
      CHKERRQ(ierr);
      ierr = SCOPFLOWSetInitilizationType(sopflow->scopflows[s],
                                          sopflow->initialization_type);
      CHKERRQ(ierr);

      ierr = SCOPFLOWSetNumContingencies(sopflow->scopflows[s], sopflow->Nc);
      CHKERRQ(ierr);
      ierr = SCOPFLOWSetTimeStepandDuration(sopflow->scopflows[s], sopflow->dT,
                                            sopflow->duration);
      CHKERRQ(ierr);
      ierr = SCOPFLOWSetLoadProfiles(
          sopflow->scopflows[s], sopflow->ploadprofile, sopflow->qloadprofile);
      CHKERRQ(ierr);
      ierr = SCOPFLOWSetWindGenProfile(sopflow->scopflows[s], sopflow->windgen);
      if (sopflow->scenfileset) {
        ierr = SCOPFLOWSetScenario(sopflow->scopflows[s],
                                   &sopflow->scenlist.scen[s]);
        CHKERRQ(ierr);
      }

      /* Set up */
      ierr = SCOPFLOWSetUp(sopflow->scopflows[s]);
      CHKERRQ(ierr);
    }
  }

  ierr = PetscCalloc1(sopflow->ns, &sopflow->nconeqcoup);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(sopflow->ns, &sopflow->nconineqcoup);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(sopflow->ns, &sopflow->nxi);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(sopflow->ns, &sopflow->ngi);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(sopflow->ns, &sopflow->xstarti);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(sopflow->ns, &sopflow->gstarti);
  CHKERRQ(ierr);

  /* Set number of variables and constraints */
  ierr = (*sopflow->modelops.setnumvariablesandconstraints)(
      sopflow, sopflow->nxi, sopflow->ngi, sopflow->nconeqcoup,
      sopflow->nconineqcoup);

  if (!sopflow->ismulticontingency) {
    sopflow->nx = sopflow->nxi[0];
    sopflow->ncon = sopflow->ngi[0];
    sopflow->nconcoup = sopflow->nconeqcoup[0] + sopflow->nconineqcoup[0];
    opflow = sopflow->opflows[0];
    sopflow->nconeq = opflow->nconeq;
    sopflow->nconineq = opflow->nconineq;

    for (i = 1; i < sopflow->ns; i++) {
      sopflow->xstarti[i] = sopflow->xstarti[i - 1] + sopflow->nxi[i - 1];
      sopflow->gstarti[i] = sopflow->gstarti[i - 1] + sopflow->ngi[i - 1];
      sopflow->nx += sopflow->nxi[i];
      sopflow->ncon += sopflow->ngi[i];
      sopflow->nconcoup += sopflow->nconeqcoup[i] + sopflow->nconineqcoup[i];
      opflow = sopflow->opflows[i];
      sopflow->nconeq += opflow->nconeq;
      sopflow->nconineq += opflow->nconineq;
    }
  } else {
    SCOPFLOW scopflow;
    sopflow->nx = sopflow->nxi[0];
    sopflow->ncon = sopflow->ngi[0];
    sopflow->nconcoup = sopflow->nconineqcoup[0];
    scopflow = sopflow->scopflows[0];
    sopflow->nconeq = scopflow->nconeq;
    sopflow->nconineq = scopflow->nconineq + scopflow->nconcoup;

    for (i = 1; i < sopflow->ns; i++) {
      sopflow->xstarti[i] = sopflow->xstarti[i - 1] + sopflow->nxi[i - 1];
      sopflow->gstarti[i] = sopflow->gstarti[i - 1] + sopflow->ngi[i - 1];
      sopflow->nx += sopflow->nxi[i];
      sopflow->ncon += sopflow->ngi[i];
      sopflow->nconcoup += sopflow->nconineqcoup[i];
      scopflow = sopflow->scopflows[i];
      sopflow->nconeq += scopflow->nconeq;
      sopflow->nconineq += scopflow->nconineq + scopflow->nconcoup;
    }
  }

  /* Create vector X */
  ierr = VecCreate(sopflow->comm->type, &sopflow->X);
  CHKERRQ(ierr);
  ierr = VecSetSizes(sopflow->X, sopflow->nx, PETSC_DECIDE);
  CHKERRQ(ierr);
  ierr = VecSetFromOptions(sopflow->X);
  CHKERRQ(ierr);
  ierr = VecGetSize(sopflow->X, &sopflow->Nx);
  CHKERRQ(ierr);

  ierr = VecDuplicate(sopflow->X, &sopflow->Xl);
  CHKERRQ(ierr);
  ierr = VecDuplicate(sopflow->X, &sopflow->Xu);
  CHKERRQ(ierr);
  ierr = VecDuplicate(sopflow->X, &sopflow->gradobj);
  CHKERRQ(ierr);

  /* vector for constraints */
  ierr = VecCreate(sopflow->comm->type, &sopflow->G);
  CHKERRQ(ierr);
  ierr = VecSetSizes(sopflow->G, sopflow->ncon, PETSC_DECIDE);
  CHKERRQ(ierr);
  ierr = VecSetFromOptions(sopflow->G);
  CHKERRQ(ierr);
  ierr = VecGetSize(sopflow->G, &sopflow->Ncon);
  CHKERRQ(ierr);

  /* Constraint bounds vectors  */
  ierr = VecDuplicate(sopflow->G, &sopflow->Gl);
  CHKERRQ(ierr);
  ierr = VecDuplicate(sopflow->G, &sopflow->Gu);
  CHKERRQ(ierr);

  if (sopflow->comm->size == 1) {
    /* Constraint Jacobian */
    ierr = MatCreate(sopflow->comm->type, &sopflow->Jac);
    CHKERRQ(ierr);
    ierr = MatSetSizes(sopflow->Jac, sopflow->ncon, sopflow->nx, sopflow->Ncon,
                       sopflow->Nx);
    CHKERRQ(ierr);
    //    ierr = MatSetFromOptions(sopflow->Jac);
    //    CHKERRQ(ierr);
    /* Assume 10% sparsity */
    ierr = MatSeqAIJSetPreallocation(sopflow->Jac,
                                     (PetscInt)(0.1 * sopflow->Nx), NULL);
    CHKERRQ(ierr);
    ierr =
        MatSetOption(sopflow->Jac, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    CHKERRQ(ierr);

    /* Hessian */
    ierr = MatCreate(sopflow->comm->type, &sopflow->Hes);
    CHKERRQ(ierr);
    ierr = MatSetSizes(sopflow->Hes, sopflow->nx, sopflow->nx, sopflow->Nx,
                       sopflow->Nx);
    CHKERRQ(ierr);
    /*    ierr = MatSetFromOptions(sopflow->Hes);
    CHKERRQ(ierr);
    */
    /* Assume 10% sparsity */
    ierr = MatSeqAIJSetPreallocation(sopflow->Hes,
                                     (PetscInt)(0.1 * sopflow->Nx), NULL);
    CHKERRQ(ierr);
    ierr =
        MatSetOption(sopflow->Hes, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    CHKERRQ(ierr);
  }

  /* Lagrangian multipliers */
  ierr = VecDuplicate(sopflow->G, &sopflow->Lambda);
  CHKERRQ(ierr);

  ierr = (*sopflow->solverops.setup)(sopflow);
  CHKERRQ(ierr);
  ExaGOLog(EXAGO_LOG_INFO, "{}", "SOPFLOW: Setup completed");

  sopflow->setupcalled = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/*
  SOPFLOWSolve - Solves the AC stochastic optimal power flow

  Input Parameters:
. sopflow - the stochastic optimal power flow application object
*/
PetscErrorCode SOPFLOWSolve(SOPFLOW sopflow) {
  PetscErrorCode ierr;
  PetscBool issolver_ipopt;

  PetscFunctionBegin;

  if (!sopflow->setupcalled) {
    ierr = SOPFLOWSetUp(sopflow);
  }

  ierr = PetscStrcmp(sopflow->solvername, "IPOPT", &issolver_ipopt);
  CHKERRQ(ierr);

  if (issolver_ipopt) {
    /* Set bounds on variables */
    if (sopflow->modelops.setvariablebounds) {
      ierr = (*sopflow->modelops.setvariablebounds)(sopflow, sopflow->Xl,
                                                    sopflow->Xu);
      CHKERRQ(ierr);
    }

    /* Set bounds on constraints */
    if (sopflow->modelops.setconstraintbounds) {
      ierr = (*sopflow->modelops.setconstraintbounds)(sopflow, sopflow->Gl,
                                                      sopflow->Gu);
      CHKERRQ(ierr);
    }

    /* Set initial guess */
    if (sopflow->modelops.setinitialguess) {
      ierr = (*sopflow->modelops.setinitialguess)(sopflow, sopflow->X);
      CHKERRQ(ierr);
    }

    ierr = VecSet(sopflow->Lambda, 1.0);
    CHKERRQ(ierr);
  }

  /* Solve */
  ierr = (*sopflow->solverops.solve)(sopflow);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  SOPFLOWGetBaseObjective - Returns the objective function value for the
base-case

  Input Parameters:
+ SOPFLOW - the SOPFLOW object
- objbase    - the objective function value for the base case

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode SOPFLOWGetBaseObjective(SOPFLOW sopflow, PetscReal *objbase) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*sopflow->solverops.getbaseobjective)(sopflow, objbase);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SOPFLOWGetTotalObjective - Returns the total objective function value

  Input Parameters:
+ SOPFLOW - the SOPFLOW object
- objtot  - the total objective function value

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode SOPFLOWGetTotalObjective(SOPFLOW sopflow, PetscReal *objtot) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*sopflow->solverops.gettotalobjective)(sopflow, objtot);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SOPFLOWGetBaseSolution - Returns the SOPFLOW solution for a given scenario

  Input Parameters:
+ SOPFLOW - the SOPFLOW object
. scennum  - Scenario number (0 for base)
- X        - the sopflow solution for the given scenario

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode SOPFLOWGetSolution(SOPFLOW sopflow, PetscInt scennum, Vec *X) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*sopflow->solverops.getsolution)(sopflow, scennum, X);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SOPFLOWGetConstraints - Returns the SOPFLOW constraints for a given scenario

  Input Parameters:
+ SOPFLOW - the SOPFLOW object
. scennum  - Scenario number (0 for base)
- G    - the sopflow constraints

  Notes: Should be called after the optimization finishes.
         Equality constraints first followed by inequality constraints
*/
PetscErrorCode SOPFLOWGetConstraints(SOPFLOW sopflow, PetscInt scennum,
                                     Vec *G) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*sopflow->solverops.getconstraints)(sopflow, scennum, G);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SOPFLOWGetConstraintMultipliers - Returns the SOPFLOW constraint multipliers
for a given scenario

  Input Parameters:
+ SOPFLOW - the SOPFLOW object
. scennum  - Scenario number (0 for base)
- G    - the sopflow constraint lagrange multipliers

  Notes: Should be called after the optimization finishes.
    Equality constraint multipliers first followed by inequality constraint
multipliers
*/
PetscErrorCode SOPFLOWGetConstraintMultipliers(SOPFLOW sopflow,
                                               PetscInt scennum, Vec *Lambda) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr =
      (*sopflow->solverops.getconstraintmultipliers)(sopflow, scennum, Lambda);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SOPFLOWGetConvergenceStatus - Did SOPFLOW converge?

  Input Parameters:
+ SOPFLOW - the SOPFLOW object
- status - PETSC_TRUE if converged, PETSC_FALSE otherwise

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode SOPFLOWGetConvergenceStatus(SOPFLOW sopflow, PetscBool *status) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*sopflow->solverops.getconvergencestatus)(sopflow, status);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SOPFLOWGetNumIterations - Returns the number of iterations for given solver
*/
PetscErrorCode SOPFLOWGetNumIterations(SOPFLOW sopflow, PetscInt *iter) {
  PetscErrorCode ierr;
  *iter = sopflow->numiter;
  PetscFunctionReturn(0);
}

/*
  SOPFLOWSetNumContingencies - Sets the number of contingencies
*/
PetscErrorCode SOPFLOWSetNumContingencies(SOPFLOW sopflow, PetscInt numctgc) {
  PetscFunctionBegin;
  sopflow->Nc = numctgc;
  PetscFunctionReturn(0);
}

/*
  SOPFLOWSetTimeStepandDuration - Sets the time-step and optimization horizon
for multi-period multi-contingecy SOPFLOW

  Input Parameters:
+ sopflow - the SOPFLOW object
. dT       - time step in minutes
- duration - duration (horizon) in hours
*/
PetscErrorCode SOPFLOWSetTimeStepandDuration(SOPFLOW sopflow, PetscReal dT,
                                             PetscReal duration) {
  PetscFunctionBegin;
  sopflow->dT = dT;
  sopflow->duration = duration;
  PetscFunctionReturn(0);
}

/*
  SOPFLOWSetLoadProfiles - Sets the data files for time-varying load profiles

  Input Parameter
+  sopflow - The SOPFLOW object
.  ploadproffile - The name of the file with real power load variationt
-  qloadproffile - The name of the file with reactive power load variationt
*/
PetscErrorCode SOPFLOWSetLoadProfiles(SOPFLOW sopflow,
                                      const char ploadprofile[],
                                      const char qloadprofile[]) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if (ploadprofile) {
    ierr = PetscMemcpy(sopflow->ploadprofile, ploadprofile,
                       PETSC_MAX_PATH_LEN * sizeof(char));
    CHKERRQ(ierr);
  }

  if (qloadprofile) {
    ierr = PetscMemcpy(sopflow->qloadprofile, qloadprofile,
                       PETSC_MAX_PATH_LEN * sizeof(char));
    CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/*
  SOPFLOWSetSubproblemModel - Set subproblem model

  Inputs
+ sopflow - sopflow object
- model   - model name

  Options
  -sopflow_subproblem_model

  Notes: This is used with HIOP solver only
*/
PetscErrorCode SOPFLOWSetSubproblemModel(SOPFLOW sopflow,
                                         const char modelname[]) {
  PetscErrorCode ierr;
  ierr = PetscStrcpy(sopflow->subproblem_model, modelname);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SOPFLOWSetSubproblemSolver - Set subproblem solver

  Inputs
+ sopflow - sopflow object
- solver   - solver name

  Option Name
  -sopflow_subproblem_solver
  Notes: This is used with HIOP solver only
*/
PetscErrorCode SOPFLOWSetSubproblemSolver(SOPFLOW sopflow,
                                          const char solvername[]) {
  PetscErrorCode ierr;
  ierr = PetscStrcpy(sopflow->subproblem_solver, solvername);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SOPFLOWSetSubproblemComputeMode - Set subproblem solver compute mode

  Inputs
+ sopflow - sopflow object
- mode   - subproblem solver compute mode

  Option Name
  -hiop_compute_mode
  Notes: This is used with HIOP solver only
*/
PetscErrorCode SOPFLOWSetSubproblemComputeMode(SOPFLOW sopflow,
                                               const char mode[]) {
  PetscErrorCode ierr;
  ierr = PetscStrcpy(sopflow->compute_mode, mode);
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
PetscErrorCode SOPFLOWSetSubproblemVerbosityLevel(SOPFLOW sopflow, int level) {
  PetscErrorCode ierr;
  sopflow->verbosity_level = level;
  PetscFunctionReturn(0);
}
