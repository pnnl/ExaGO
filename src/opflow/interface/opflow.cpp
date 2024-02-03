#include <exago_config.h>
#include <petsc/private/dmnetworkimpl.h>
#include <private/opflowimpl.h>
#include <private/pflowimpl.h>

const char *const OPFLOWInitializationTypes[] = {
    "MIDPOINT",    "FROMFILE", "ACPF",
    "FLATSTART",   "DCOPF",    "OPFLOWInitializationType",
    "OPFLOWINIT_", NULL};

const char *const OPFLOWObjectiveTypes[] = {"MIN_GEN_COST",
                                            "MIN_GENSETPOINT_DEVIATION",
                                            "NO_OBJ",
                                            "OPFLOWObjectiveType",
                                            "",
                                            NULL};

const char *const OPFLOWGenBusVoltageTypes[] = {"VARIABLE_WITHIN_BOUNDS",
                                                "FIXED_WITHIN_QBOUNDS",
                                                "FIXED_AT_SETPOINT",
                                                "OPFLOWGenBusVoltageType",
                                                "",
                                                NULL};

const char *const OPFLOWOutputFormatTypes[] = {
    "MATPOWER", "CSV", "JSON", "MINIMAL", "OutputFormat", "", NULL};

void swap_dm(DM *dm1, DM *dm2) {
  DM temp = *dm1;
  *dm1 = *dm2;
  *dm2 = temp;
}

/* Converts an array xin in natural ordering to an array xout in sparse-dense
   ordering
   Used only with HIOP solver
*/
PetscErrorCode OPFLOWNaturalToSpDense(OPFLOW opflow, const double *xin,
                                      double *xout) {
  int i;

  PetscFunctionBegin;
  for (i = 0; i < opflow->nx; i++) {
    xout[opflow->idxn2sd_map[i]] = xin[i];
  }

  PetscFunctionReturn(0);
}

/* Converts an array xin in sparse dense ordering to an array xout in natural
   ordering
   Used only with HIOP solver
*/
PetscErrorCode OPFLOWSpDenseToNatural(OPFLOW opflow, const double *xin,
                                      double *xout) {
  int i;

  PetscFunctionBegin;
  for (i = 0; i < opflow->nx; i++) {
    xout[i] = xin[opflow->idxn2sd_map[i]];
  }
  PetscFunctionReturn(0);
}

/* Sets the list of lines monitored */
PetscErrorCode OPFLOWGetLinesMonitored(OPFLOW opflow) {
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  PSLINE line;
  PetscInt i, j;

  PetscFunctionBegin;
  if (opflow->nlinesmon) {
    /* Number of lines monitored already set through file. */
    PetscFunctionReturn(0);
  }
  if (opflow->nlinekvmon == -1) {
    opflow->nlinekvmon = opflow->ps->nkvlevels;
    ierr = PetscMemcpy(opflow->linekvmon, ps->kvlevels,
                       opflow->nlinekvmon * sizeof(PetscScalar));
  }
  ierr = PetscMalloc1(ps->nline, &opflow->linesmon);
  CHKERRQ(ierr);
  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    if (!line->status || line->rateA > 1e5)
      continue;
    for (j = 0; j < opflow->nlinekvmon; j++) {
      if (PetscAbsScalar(line->kvlevel - opflow->linekvmon[j]) < 1e-5) {
        opflow->linesmon[opflow->nlinesmon++] = i;
        break;
      }
    }
  }
  PetscFunctionReturn(0);
}

/*
   OPFLOWSetOutputFormat - Sets the format for solving the solution

   Input Parameters:
+  opflow - the opflow object
-  otype  - output format type

   Command-line option: -opflow_output_format
   Notes: Must be called before OPFLOWSetUp()
*/
PetscErrorCode OPFLOWSetOutputFormat(OPFLOW opflow, OutputFormat otype) {
  PetscFunctionBegin;
  opflow->outputformat = otype;
  PetscFunctionReturn(0);
}

/*
   OPFLOWSetGenBusVoltageType - Sets the voltage control mode for generator
buses

   Input Parameters:
+  opflow - the opflow object
-  vtype  - voltage control type

   Command-line option: -opflow_genbusvoltage
   Notes: Must be called before OPFLOWSetUp()
*/
PetscErrorCode OPFLOWSetGenBusVoltageType(OPFLOW opflow,
                                          OPFLOWGenBusVoltageType vtype) {
  PetscFunctionBegin;
  opflow->genbusvoltagetype = vtype;
  PetscFunctionReturn(0);
}

/*
   OPFLOWGetGenBusVoltageType - Sets the voltage control mode for generator
buses

   Input Parameters:
+  opflow - the opflow object
-  vtype  - voltage control type pointer

   Command-line option: -opflow_genbusvoltage
*/
PetscErrorCode OPFLOWGetGenBusVoltageType(OPFLOW opflow,
                                          OPFLOWGenBusVoltageType *vtype) {
  PetscFunctionBegin;
  *vtype = opflow->genbusvoltagetype;
  PetscFunctionReturn(0);
}

/*
  OPFLOWSkipOptions - Skip run-time options?

  Inputs:
+ opflow - OPFLOW object
- skip_options - Skip run-time options?

  Notes: Should be called before OPFLOWSetUp
*/
PetscErrorCode OPFLOWSkipOptions(OPFLOW opflow, PetscBool skip_options) {
  PetscFunctionBegin;
  opflow->skip_options = skip_options;
  PetscFunctionReturn(0);
}
/*
  OPFLOWHasLoadLoss - Use load loss in the OPFLOW formulation

  Input Parameters:
+ opflow - OPFLOW object
- hasloadloss - Use load loss?

  Command-line option: -opflow_include_loadloss_variables
  Notes: Should be called before OPFLOWSetUp
*/
PetscErrorCode OPFLOWHasLoadLoss(OPFLOW opflow, PetscBool hasloadloss) {
  PetscFunctionBegin;
  opflow->include_loadloss_variables = hasloadloss;
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetHasLoadLoss - Use load loss in the OPFLOW formulation

  Input Parameters:
+ opflow - OPFLOW object
- hasloadloss - Use load loss?

*/
PetscErrorCode OPFLOWGetHasLoadLoss(OPFLOW opflow, PetscBool *hasloadloss) {
  PetscFunctionBegin;
  *hasloadloss = opflow->include_loadloss_variables;
  PetscFunctionReturn(0);
}

/*
  OPFLOWHasBusPowerImbalance - Use bus power imbalance in the OPFLOW formulation

  Input Parameters:
+ opflow - OPFLOW object
- hasbuspowerimbalance - Use power imbalance?

  Command-line option: -opflow_include_powerimbalance_variables
  Notes: Should be called before OPFLOWSetUp

*/
PetscErrorCode OPFLOWHasBusPowerImbalance(OPFLOW opflow,
                                          PetscBool hasbuspowerimbalance) {
  PetscFunctionBegin;
  opflow->include_powerimbalance_variables = hasbuspowerimbalance;
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetHasBusPowerImbalance - Use bus power imbalance in the OPFLOW
formulation

  Input Parameters:
+ opflow - OPFLOW object
- hasbuspowerimbalance - Use power imbalance?
*/
PetscErrorCode OPFLOWGetHasBusPowerImbalance(OPFLOW opflow,
                                             PetscBool *hasbuspowerimbalance) {
  PetscFunctionBegin;
  *hasbuspowerimbalance = opflow->include_powerimbalance_variables;
  PetscFunctionReturn(0);
}

/*
  OPFLOWSetLoadLossPenalty - Set load loss penalty

  Input Parameters:
+ opflow - OPFLOW object
- penalty - penalty for load loss (for each load)

  Command-line option: -opflow_loadloss_penalty
  Notes: Should be called before OPFLOWSetUp

*/
PetscErrorCode OPFLOWSetLoadLossPenalty(OPFLOW opflow, PetscReal penalty) {
  PetscFunctionBegin;
  opflow->loadloss_penalty = penalty;
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetLoadLossPenalty - Get load loss penalty

  Input Parameters:
+ opflow - OPFLOW object
- penalty - penalty for load loss (for each load)

*/
PetscErrorCode OPFLOWGetLoadLossPenalty(OPFLOW opflow, PetscReal *penalty) {
  PetscFunctionBegin;
  *penalty = opflow->loadloss_penalty;
  PetscFunctionReturn(0);
}

/*
  OPFLOWSetBusPowerImbalancePenalty - Set bus power imbalance penalty

  Input Parameters:
+ opflow - OPFLOW object
- penalty - penalty for power imbalance (at each bus)

  Command-line option: -opflow_powerimbalance_penalty
  Notes: Should be called before OPFLOWSetUp

*/
PetscErrorCode OPFLOWSetBusPowerImbalancePenalty(OPFLOW opflow,
                                                 PetscReal penalty) {
  PetscFunctionBegin;
  opflow->powerimbalance_penalty = penalty;
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetBusPowerImbalancePenalty - Get bus power imbalance penalty

  Input Parameters:
+ opflow - OPFLOW object
- penalty - penalty for power imbalance (at each bus)
*/
PetscErrorCode OPFLOWGetBusPowerImbalancePenalty(OPFLOW opflow,
                                                 PetscReal *penalty) {
  PetscFunctionBegin;
  *penalty = opflow->powerimbalance_penalty;
  PetscFunctionReturn(0);
}

/*
  OPFLOWSetHIOPComputeMode - Set compute mode for HIOP solver
  Input parameters
+ opflow - OPFLOW object
- mode - mode for HIOP solver

  Command-line option: -hiop_compute_mode
*/
PetscErrorCode OPFLOWSetHIOPComputeMode(OPFLOW opflow, const std::string mode) {
  PetscFunctionBegin;
  opflow->_p_hiop_compute_mode = mode;
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetHIOPComputeMode - Get compute mode for HIOP solver
  Input parameters
+ opflow - OPFLOW object
- mode - mode for HIOP solver
*/
PetscErrorCode OPFLOWGetHIOPComputeMode(OPFLOW opflow, std::string *mode) {
  PetscFunctionBegin;
  *mode = opflow->_p_hiop_compute_mode;
  PetscFunctionReturn(0);
}

/*
  OPFLOWSetHIOPMemSpace - Set mem space for HIOP solver
  Input parameters
+ opflow - OPFLOW object
- mem_space - memory space for HIOP solver

  Command-line option: -hiop_mem_space
*/
PetscErrorCode OPFLOWSetHIOPMemSpace(OPFLOW opflow,
                                     const std::string mem_space) {
  PetscFunctionBegin;
  if (mem_space == "DEFAULT")
    opflow->mem_space = static_cast<HIOPMemSpace>(0);
  else if (mem_space == "HOST")
    opflow->mem_space = static_cast<HIOPMemSpace>(1);
  else if (mem_space == "UM")
    opflow->mem_space = static_cast<HIOPMemSpace>(2);
  else if (mem_space == "DEVICE")
    opflow->mem_space = static_cast<HIOPMemSpace>(3);
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetHIOPMemSpace - Get mem space for HIOP solver
  Input parameters
+ opflow - OPFLOW object
- mem_space - memory space for HIOP solver
*/
PetscErrorCode OPFLOWGetHIOPMemSpace(OPFLOW opflow, std::string *mem_space) {
  PetscFunctionBegin;
  if (opflow->mem_space == 0)
    *mem_space = "DEFAULT";
  else if (opflow->mem_space == 1)
    *mem_space = "HOST";
  else if (opflow->mem_space == 2)
    *mem_space = "UM";
  else if (opflow->mem_space == 3)
    *mem_space = "DEVICE";
  PetscFunctionReturn(0);
}

/*
  OPFLOWSetHIOPVerbosityLevel - Set verbosity level for HIOP solver
  Input parameters
+ opflow - OPFLOW object
- level - verbositiy level for HIOP solver

  Command-line option: -hiop_compute_mode
*/
PetscErrorCode OPFLOWSetHIOPVerbosityLevel(OPFLOW opflow, int level) {
  PetscFunctionBegin;
  opflow->_p_hiop_verbosity_level = level;
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetHIOPVerbosityLevel - Get verbosity level for HIOP solver
  Input parameters
+ opflow - OPFLOW object
- level - verbositiy level for HIOP solver
*/
PetscErrorCode OPFLOWGetHIOPVerbosityLevel(OPFLOW opflow, int *level) {
  PetscFunctionBegin;
  *level = opflow->_p_hiop_verbosity_level;
  PetscFunctionReturn(0);
}

/*
  OPFLOWHasGenSetPoint - Use gen. set point in the OPFLOW formulation

  Input Parameters:
+ opflow - OPFLOW object
- hassetpoint - Use set-point?

  Notes: Should be called before OPFLOWSetUp
*/
PetscErrorCode OPFLOWHasGenSetPoint(OPFLOW opflow, PetscBool hassetpoint) {
  PetscFunctionBegin;
  opflow->has_gensetpoint = hassetpoint;
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetHasGenSetPoint - Use gen. set point in the OPFLOW formulation

  Input Parameters:
+ opflow - OPFLOW object
- hassetpoint - Use set-point?

*/
PetscErrorCode OPFLOWGetHasGenSetPoint(OPFLOW opflow, PetscBool *hassetpoint) {
  PetscFunctionBegin;
  *hassetpoint = opflow->has_gensetpoint;
  PetscFunctionReturn(0);
}

/*
  OPFLOWUseAGC - Uses AGC to proportionally redispatch the generators

  Input Parameters:
+ opflow - OPFLOW object
- useagc - Use AGC?

  Notes: Should be called before OPFLOWSetUp
*/
PetscErrorCode OPFLOWUseAGC(OPFLOW opflow, PetscBool useagc) {
  PetscFunctionBegin;
  opflow->use_agc = useagc;
  if (useagc)
    opflow->has_gensetpoint = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetUseAGC - Uses AGC to proportionally redispatch the generators

  Input Parameters:
+ opflow - OPFLOW object
- useagc - Use AGC?
*/
PetscErrorCode OPFLOWGetUseAGC(OPFLOW opflow, PetscBool *useagc) {
  PetscFunctionBegin;
  *useagc = opflow->use_agc;
  PetscFunctionReturn(0);
}
/*
  OPFLOWGetSizes - Returns the size of solution vector, and constraints

  Input Parameters:
- opflow - the OPFLOW object

  Output Parameters:
+ nx - size of X vector
. nconeq - size of equality constraints vector
- nconineq - size of inequality constraints vector

  Notes:
  Should be called after OPFLOWSetUp()
*/
PetscErrorCode OPFLOWGetSizes(OPFLOW opflow, int *nx, int *nconeq,
                              int *nconineq) {
  PetscFunctionBegin;
  *nx = opflow->nx;
  *nconeq = opflow->nconeq;
  *nconineq = opflow->nconineq;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWGetVariableOrdering(OPFLOW opflow, int **ordering) {
  PetscFunctionBegin;
  *ordering = opflow->idxn2sd_map;
  PetscFunctionReturn(0);
}

/*
   OPFLOWInitializeDCOPF - DCOPF initialization of OPFLOW
*/
PetscErrorCode OPFLOWInitializeDCOPF(OPFLOW opflow, PetscBool *converged) {
  PetscErrorCode ierr;
  OPFLOW dcopflow; // DCOPF object

  PetscFunctionBegin;

  ierr = OPFLOWCreate(opflow->comm->type, &dcopflow);
  CHKERRQ(ierr);

  ierr = OPFLOWReadMatPowerData(dcopflow, opflow->ps->net_file_name);
  CHKERRQ(ierr);

  ierr = OPFLOWSetModel(dcopflow, OPFLOWMODEL_DCOPF);
  CHKERRQ(ierr);

#if defined(EXAGO_ENABLE_IPOPT)
  ierr = OPFLOWSetSolver(dcopflow, OPFLOWSOLVER_IPOPT);
  CHKERRQ(ierr);
#else
  ierr = OPFLOWSetSolver(dcopflow, OPFLOWSOLVER_HIOPSPARSE);
  CHKERRQ(ierr);
#endif

  /* Set options */
  ierr = OPFLOWHasLoadLoss(dcopflow, opflow->include_loadloss_variables);
  CHKERRQ(ierr);
  ierr = OPFLOWSetLoadLossPenalty(dcopflow, opflow->loadloss_penalty);
  CHKERRQ(ierr);

  ierr = OPFLOWHasBusPowerImbalance(dcopflow,
                                    opflow->include_powerimbalance_variables);
  CHKERRQ(ierr);
  ierr = OPFLOWSetBusPowerImbalancePenalty(dcopflow,
                                           opflow->powerimbalance_penalty);
  CHKERRQ(ierr);

  ierr = OPFLOWHasGenSetPoint(dcopflow, opflow->has_gensetpoint);
  CHKERRQ(ierr);

  ierr = OPFLOWUseAGC(dcopflow, opflow->use_agc);
  CHKERRQ(ierr);

  ierr = OPFLOWIgnoreLineflowConstraints(dcopflow,
                                         opflow->ignore_lineflow_constraints);
  CHKERRQ(ierr);
  ierr =
      OPFLOWAllowLineflowViolation(dcopflow, opflow->allow_lineflow_violation);
  CHKERRQ(ierr);

  ierr = OPFLOWSetTolerance(dcopflow, opflow->tolerance);
  CHKERRQ(ierr);

  ierr = OPFLOWSetObjectiveType(dcopflow, opflow->objectivetype);
  CHKERRQ(ierr);

  ierr = OPFLOWSkipOptions(dcopflow, PETSC_TRUE);
  CHKERRQ(ierr);

  ierr = OPFLOWSolve(dcopflow);
  CHKERRQ(ierr);

  ierr = OPFLOWSolutionToPS(dcopflow);
  CHKERRQ(ierr);
  ierr = OPFLOWGetConvergenceStatus(dcopflow, converged);
  CHKERRQ(ierr);

  ierr = PSCopy(dcopflow->ps, opflow->ps);
  CHKERRQ(ierr);

  ierr = OPFLOWDestroy(&dcopflow);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  OPFLOWSetUpInitPflow - Sets up the power flow solver object used in obtaining
  initial conditions.

  Note: This power flow solver object shares the underlying power system object
  and all the data. It is created at the end of OPFLOWSetUp. So, it already has
  the distributed ps object. The only difference is that it uses a different
  PetscSection for storing the degrees of freedom that gets associated with the
  dmnetwork (and plex)
*/
PetscErrorCode OPFLOWSetUpInitPflow(OPFLOW opflow) {
  PetscErrorCode ierr;
  PetscInt vStart, vEnd, eStart, eEnd;
  DM networkdm, plexdm;
  PetscInt i;

  PetscFunctionBegin;

  networkdm = opflow->ps->networkdm;

  ierr = PFLOWCreate(opflow->comm->type, &opflow->initpflow);
  CHKERRQ(ierr);

  /* PFLOW creates a new PS object. Destroy it so that we can associate the
     PS from OPFLOW with initpflow
  */
  ierr = PSDestroy(&opflow->initpflow->ps);
  CHKERRQ(ierr);

  opflow->initpflow->ps = opflow->ps;
  /* Increase the reference count for opflow->ps */
  ierr = PSIncreaseReferenceCount(opflow->ps);
  CHKERRQ(ierr);

  ierr = PetscSectionCreate(opflow->comm->type, &opflow->initpflowpsection);
  CHKERRQ(ierr);

  ierr = DMNetworkGetVertexRange(networkdm, &vStart, &vEnd);
  CHKERRQ(ierr);
  ierr = DMNetworkGetEdgeRange(networkdm, &eStart, &eEnd);
  CHKERRQ(ierr);

  ierr = PetscSectionSetChart(opflow->initpflowpsection, eStart, vEnd);
  CHKERRQ(ierr);

  for (i = vStart; i < vEnd; i++) {
    /* Two variables at each bus/vertex */
    ierr = PetscSectionSetDof(opflow->initpflowpsection, i, 2);
    CHKERRQ(ierr);
  }
  ierr = PetscSectionSetUp(opflow->initpflowpsection);
  CHKERRQ(ierr);

  /* Clone DM to be used with initpflow */
  ierr = DMClone(networkdm, &opflow->initpflowdm);
  CHKERRQ(ierr);

  /* Set initpflowdm in opflow->ps->networkdm, the previous networkdm get
   * stashed in opflow->initpflowdm */
  swap_dm(&opflow->ps->networkdm, &opflow->initpflowdm);
  networkdm = opflow->ps->networkdm;

  /* Get the plex dm */
  ierr = DMNetworkGetPlex(networkdm, &plexdm);
  CHKERRQ(ierr);

  /* Get default sections associated with this plex */
  ierr = DMGetSection(plexdm, &opflow->defaultsection);
  CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)opflow->defaultsection);
  CHKERRQ(ierr);

  ierr = DMGetGlobalSection(plexdm, &opflow->defaultglobalsection);
  CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)opflow->defaultglobalsection);
  CHKERRQ(ierr);

  /* Set the new section created for initial power flow */
  ierr = DMSetSection(plexdm, opflow->initpflowpsection);
  CHKERRQ(ierr);
  ierr = DMGetGlobalSection(plexdm, &opflow->initpflowpglobsection);
  CHKERRQ(ierr);

  /* Set up PFLOW object. Note pflow->ps will not be set up again as it has
     been already set up by opflow
  */
  ierr = PFLOWSetUp(opflow->initpflow);
  CHKERRQ(ierr);

  /*
     Update the ref. counts for init pflow sections so that they do not
     get destroyed when DMSetDefaultSection is called
  */
  // ierr =
  // PetscObjectReference((PetscObject)opflow->initpflowpsection);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)opflow->initpflowpglobsection);
  CHKERRQ(ierr);

  /* Reset the sections */
  ierr = DMSetSection(plexdm, opflow->defaultsection);
  CHKERRQ(ierr);
  ierr = DMSetGlobalSection(plexdm, opflow->defaultglobalsection);
  CHKERRQ(ierr);

  /* Reset dm */
  swap_dm(&opflow->ps->networkdm, &opflow->initpflowdm);

  PetscFunctionReturn(0);
}

/*
  OPFLOWComputePrePflow - Computes the steady-state power flow for initializing
optimal power flow

  Input Parameters:
. opflow - the OPFLOW object
*/
PetscErrorCode OPFLOWComputePrePflow(OPFLOW opflow, PetscBool *converged) {
  PetscErrorCode ierr;
  DM plexdm;
  PFLOW initpflow = opflow->initpflow;
  PS ps = opflow->ps;

  PetscFunctionBegin;
  /* Set initpflowdm for solving the power flow */
  swap_dm(&opflow->ps->networkdm, &opflow->initpflowdm);

  ierr = DMNetworkGetPlex(opflow->ps->networkdm, &plexdm);
  CHKERRQ(ierr);

  /* Increase the ref. counts for the default sections so that they do not get
     destroyed when DMSetDefaultXXX is called with the initpflowxxx sections
  */
  ierr = PetscObjectReference((PetscObject)opflow->defaultsection);
  CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)opflow->defaultglobalsection);
  CHKERRQ(ierr);

  /* Set the new section created for initial power flow. */
  ierr = DMSetSection(plexdm, opflow->initpflowpsection);
  CHKERRQ(ierr);
  ierr = DMSetGlobalSection(plexdm, opflow->initpflowpglobsection);
  CHKERRQ(ierr);

  ierr = DMCreateSectionSF(plexdm, opflow->initpflowpsection,
                           opflow->initpflowpglobsection);
  CHKERRQ(ierr);

  /* Reset the edge and bus starting locations of variables */
  ierr = PSSetEdgeandBusStartLoc(ps);
  CHKERRQ(ierr);

  /* Solve */
  ierr = PFLOWSolve(initpflow);
  CHKERRQ(ierr);
  ierr = PFLOWConverged(initpflow, converged);
  CHKERRQ(ierr);

  /* Update bus and gen structs in pflow->ps */
  ierr = PFLOWSolutionToPS(initpflow);
  CHKERRQ(ierr);

  /*
     Update the ref. counts for init pflow sections so that they do not
     get destroyed when DMSetSection is called
  */
  ierr = PetscObjectReference((PetscObject)opflow->initpflowpsection);
  CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)opflow->initpflowpglobsection);
  CHKERRQ(ierr);

  /* Reset the sections */
  ierr = DMSetSection(plexdm, opflow->defaultsection);
  CHKERRQ(ierr);
  ierr = DMSetGlobalSection(plexdm, opflow->defaultglobalsection);
  CHKERRQ(ierr);

  /* Reset SF */
  ierr = DMCreateSectionSF(plexdm, opflow->defaultsection,
                           opflow->defaultglobalsection);
  CHKERRQ(ierr);

  /* Reset the bus and edge starting locations for variables */
  ierr = PSSetEdgeandBusStartLoc(ps);
  CHKERRQ(ierr);

  /* Reset dm */
  swap_dm(&opflow->initpflowdm, &opflow->ps->networkdm);
  PetscFunctionReturn(0);
}

/*
  OPFLOWCreate - Creates an optimal power flow application object

  Input Parameters
. mpicomm - The MPI communicator

  Output Parameters
. opflowout - The optimal power flow application object

  If any of these default values are updated, the corresponding option in
  include/opflow.h must be updated as well.
*/
PetscErrorCode OPFLOWCreate(MPI_Comm mpicomm, OPFLOW *opflowout) {
  PetscErrorCode ierr;
  OPFLOW opflow;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1, &opflow);
  CHKERRQ(ierr);

  ierr = COMMCreate(mpicomm, &opflow->comm);
  CHKERRQ(ierr);

  ierr = PSCreate(mpicomm, &opflow->ps);
  CHKERRQ(ierr);

  /* Set the application with the PS object */
  ierr = PSSetApplication(opflow->ps, opflow, APP_ACOPF);
  CHKERRQ(ierr);

  opflow->modelname = OPFLOWOptions::model.opt;
  opflow->solvername = OPFLOWOptions::solver.opt;

  opflow->Nconeq = opflow->nconeq = -1;
  opflow->Nconineq = opflow->nconineq = -1;
  opflow->Ncon = opflow->ncon = -1;
  opflow->Nx = opflow->nx = -1;
  opflow->Gi = NULL;
  opflow->Lambdai = NULL;

  opflow->obj_factor = 1.0;
  opflow->obj = 0.0;
  opflow->obj_gencost =
      PETSC_TRUE; /* Generation cost minimization ON by default */
  opflow->solutiontops = PETSC_FALSE;
  opflow->tolerance = OPFLOWOptions::tolerance.default_value;

  opflow->has_gensetpoint = OPFLOWOptions::has_gensetpoint.default_value;
  opflow->use_agc = OPFLOWOptions::use_agc.default_value;

  opflow->solver = NULL;
  opflow->model = NULL;

  opflow->initializationtype = static_cast<OPFLOWInitializationType>(
      OPFLOWOptions::initialization.ToEnum(OPFLOWInitializationTypes, 4,
                                           /*prefix=*/"OPFLOWINIT_"));

  opflow->objectivetype = static_cast<OPFLOWObjectiveType>(
      OPFLOWOptions::objective.ToEnum(OPFLOWObjectiveTypes, 3));

  opflow->genbusvoltagetype = static_cast<OPFLOWGenBusVoltageType>(
      OPFLOWOptions::genbusvoltage.ToEnum(OPFLOWGenBusVoltageTypes, 3));

  opflow->outputformat = static_cast<OutputFormat>(
      OPFLOWOptions::outputformat.ToEnum(OPFLOWOutputFormatTypes, 4));

  opflow->nmodelsregistered = opflow->nsolversregistered = 0;
  opflow->OPFLOWModelRegisterAllCalled = opflow->OPFLOWSolverRegisterAllCalled =
      PETSC_FALSE;

  /* Register all models */
  ierr = OPFLOWModelRegisterAll(opflow);

  /* Register all solvers */
  ierr = OPFLOWSolverRegisterAll(opflow);

  /* Run-time options */
  opflow->ignore_lineflow_constraints =
      OPFLOWOptions::ignore_lineflow_constraints.default_value;
  opflow->allow_lineflow_violation =
      OPFLOWOptions::allow_lineflow_violation.default_value;
  opflow->lineflowviolation_penalty =
      OPFLOWOptions::lineflowviolation_penalty.default_value;
  opflow->include_loadloss_variables =
      OPFLOWOptions::include_loadloss_variables.default_value;
  opflow->include_powerimbalance_variables =
      OPFLOWOptions::include_powerimbalance_variables.default_value;

  opflow->loadloss_penalty = OPFLOWOptions::loadloss_penalty.default_value;
  opflow->powerimbalance_penalty =
      OPFLOWOptions::powerimbalance_penalty.default_value;

  opflow->spdnordering = PETSC_FALSE;

  opflow->nlinekvmon = -1;
  ierr = PetscMalloc1(MAX_KV_LEVELS, &opflow->linekvmon);
  CHKERRQ(ierr);
  opflow->nlinesmon = 0;
  opflow->linesmon = NULL;

  opflow->_p_hiop_compute_mode = "auto";
  opflow->_p_hiop_verbosity_level = 0;

  opflow->skip_options = PETSC_FALSE;

  opflow->weight = 1.0;
  opflow->setupcalled = PETSC_FALSE;

  *opflowout = opflow;

  //  ierr = PetscPrintf(opflow->comm->type,"OPFLOW: Application created\n");
  PetscFunctionReturn(0);
}

/*
  OPFLOWSetWeight - Set the weight (probability) for this optimal power flow

  Input Parameter
+ opflow - optimal power flow object
- weight - weighting factor for the system conditions in (0,1) (default = 1.0)

  Notes:
    This routine sets the weight for a given optimal power flow. The weight
    is used for weighing the objective function weight*f(x). The values needs to
be between 0 and 1 (both included)
*/
PetscErrorCode OPFLOWSetWeight(OPFLOW opflow, PetscScalar weight) {
  PetscFunctionBegin;
  opflow->weight = weight;
  PetscFunctionReturn(0);
}

/*
  OPFLOWDestroy - Destroys the optimal power flow application object

  Input Parameter
. opflow - The OPFLOW object to destroy
*/
PetscErrorCode OPFLOWDestroy(OPFLOW *opflow) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PSDestroy(&(*opflow)->ps);
  CHKERRQ(ierr);

  if (!(*opflow)->setupcalled) {
    PetscFunctionReturn(0);
  }

  if ((*opflow)->initializationtype == OPFLOWINIT_ACPF) {
    /* Destroy objects created for initial power flow */
    ierr = PFLOWDestroy(&(*opflow)->initpflow);
    CHKERRQ(ierr);

    ierr = DMDestroy(&(*opflow)->initpflowdm);
    CHKERRQ(ierr);

    PetscObjectDereference((PetscObject)((*opflow)->initpflowpsection));
    PetscObjectDereference((PetscObject)((*opflow)->initpflowpglobsection));
    ierr = PetscSectionDestroy(&(*opflow)->initpflowpsection);
    CHKERRQ(ierr);
    ierr = PetscSectionDestroy(&(*opflow)->initpflowpglobsection);
    CHKERRQ(ierr);
  }

  ierr = COMMDestroy(&(*opflow)->comm);
  CHKERRQ(ierr);

  /* Solution vector */
  ierr = VecDestroy(&(*opflow)->X);
  CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->localX);
  CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->gradobj);
  CHKERRQ(ierr);

  /* Lower and upper bounds on X */
  ierr = VecDestroy(&(*opflow)->Xl);
  CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Xu);
  CHKERRQ(ierr);

  /* Constraints vector */
  ierr = VecDestroy(&(*opflow)->G);
  CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Ge);
  CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Gelocal);
  CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Lambdae);
  CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Lambdaelocal);
  CHKERRQ(ierr);
  if ((*opflow)->nconineq) {
    ierr = VecDestroy(&(*opflow)->Gi);
    CHKERRQ(ierr);
    ierr = VecDestroy(&(*opflow)->Lambdai);
    CHKERRQ(ierr);
  }
  ierr = VecDestroy(&(*opflow)->Gl);
  CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Gu);
  CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Lambda);
  CHKERRQ(ierr);

  /* Index sets and vecscatter for equality constraints */
  ierr = ISDestroy(&(*opflow)->isconeqlocal);
  CHKERRQ(ierr);
  ierr = ISDestroy(&(*opflow)->isconeqglob);
  CHKERRQ(ierr);
  ierr = VecScatterDestroy(&(*opflow)->scattereqcon);
  CHKERRQ(ierr);

  /* Jacobian of constraints */
  ierr = MatDestroy(&(*opflow)->Jac);
  CHKERRQ(ierr);
  ierr = MatDestroy(&(*opflow)->Jac_Ge);
  CHKERRQ(ierr);
  ierr = MatDestroy(&(*opflow)->Jac_Gi);
  CHKERRQ(ierr);

  ierr = MatDestroy(&(*opflow)->Hes);
  CHKERRQ(ierr);

  if ((*opflow)->modelops.destroy) {
    ierr = ((*opflow)->modelops.destroy)(*opflow);
  }

  if ((*opflow)->solverops.destroy) {
    ierr = ((*opflow)->solverops.destroy)(*opflow);
  }

  ierr = PetscFree((*opflow)->idxn2sd_map);
  CHKERRQ(ierr);
  ierr = PetscFree((*opflow)->busnvararray);
  CHKERRQ(ierr);
  ierr = PetscFree((*opflow)->branchnvararray);
  CHKERRQ(ierr);
  ierr = PetscFree((*opflow)->eqconglobloc);
  CHKERRQ(ierr);
  ierr = PetscFree((*opflow)->linekvmon);
  CHKERRQ(ierr);
  ierr = PetscFree((*opflow)->linesmon);
  CHKERRQ(ierr);

  ierr = PetscFree(*opflow);
  CHKERRQ(ierr);
  *opflow = 0;
  PetscFunctionReturn(0);
}

/*
   OPFLOWSetTolerance - Set the solver tolerance

  Input Parameters:
+ opflow - opflow application object
- tol    - solver tolerance
*/
PetscErrorCode OPFLOWSetTolerance(OPFLOW opflow, PetscReal tol) {
  PetscFunctionBegin;
  opflow->tolerance = tol;
  PetscFunctionReturn(0);
}

/**
 * OPFLOWGetNumIterations - Get the solver iterations
 * Input Parameters:
 * opflow - opflow application object
 * tol    - pointer to solver iterations variable that will be set
 */
PetscErrorCode OPFLOWGetNumIterations(OPFLOW opflow, PetscInt *its) {
  PetscFunctionBegin;
  *its = opflow->numits;
  PetscFunctionReturn(0);
}

/**
 * OPFLOWGetTolerance - Get the solver tolerance
 * Input Parameters:
 * opflow - opflow application object
 * tol    - pointer to solver tolerance variable that will be set
 */
PetscErrorCode OPFLOWGetTolerance(OPFLOW opflow, PetscReal *tol) {
  PetscFunctionBegin;
  *tol = opflow->tolerance;
  PetscFunctionReturn(0);
}

/*
  OPFLOWSetSolver - Sets the solver for OPFLOW

  Input Parameters:
+ opflow - opflow application object
- solvername - name of the solver

  Note: must be called before OPFLOWSetUp
*/
PetscErrorCode OPFLOWSetSolver(OPFLOW opflow, const std::string solvername) {
  PetscErrorCode ierr, (*r)(OPFLOW) = NULL;
  PetscInt i;
  PetscFunctionBegin;
  PetscBool match;
  for (i = 0; i < opflow->nsolversregistered; i++) {
    ierr = PetscStrcmp(opflow->OPFLOWSolverList[i].name, solvername.c_str(),
                       &match);
    CHKERRQ(ierr);
    if (match) {
      r = opflow->OPFLOWSolverList[i].create;
      break;
    }
  }

  if (!r)
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE,
            "Unknown type for OPFLOW Solver %s. You may need to rebuild ExaGO "
            "to support this solver.",
            solvername.c_str());

  /* Initialize (Null) the function pointers */
  opflow->solverops.destroy = 0;
  opflow->solverops.solve = 0;
  opflow->solverops.setup = 0;

  opflow->solvername = solvername;
  /* Call the underlying implementation constructor */
  ierr = (*r)(opflow);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  OPFLOWGetSolver - Gets the solver for OPFLOW

  Input Parameters:
+ opflow - opflow application object
- solvername - name of the solver
*/
PetscErrorCode OPFLOWGetSolver(OPFLOW opflow, std::string *solvername) {
  PetscFunctionBegin;
  *solvername = opflow->solvername;
  PetscFunctionReturn(0);
}

/*
  OPFLOWSetModel - Sets the model for OPFLOW

  Input Parameters:
+ opflow - opflow application object
- modelname - name of the model

  Note: must be called before OPFLOWSetUp
*/
PetscErrorCode OPFLOWSetModel(OPFLOW opflow, const std::string modelname) {
  PetscErrorCode ierr, (*r)(OPFLOW) = NULL;
  PetscInt i;
  PetscFunctionBegin;
  PetscBool match;
  for (i = 0; i < opflow->nmodelsregistered; i++) {
    ierr =
        PetscStrcmp(opflow->OPFLOWModelList[i].name, modelname.c_str(), &match);
    CHKERRQ(ierr);
    if (match) {
      r = opflow->OPFLOWModelList[i].create;
      break;
    }
  }

  if (!r)
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE,
            "Unknown type for OPFLOW Model %s. You may need to rebuild ExaGO "
            "with the solver that supports this model.",
            modelname.c_str());

  /* Null the function pointers */
  opflow->modelops.destroy = 0;
  opflow->modelops.setup = 0;
  opflow->modelops.setnumvariables = 0;
  opflow->modelops.setnumconstraints = 0;
  opflow->modelops.updatevariablebounds = 0;
  opflow->modelops.setvariablebounds = 0;
  opflow->modelops.setvariableboundsarray = 0;
  opflow->modelops.setconstraintbounds = 0;
  opflow->modelops.setconstraintboundsarray = 0;
  opflow->modelops.setvariableandconstraintbounds = 0;
  opflow->modelops.setvariableandconstraintboundsarray = 0;
  opflow->modelops.setinitialguess = 0;
  opflow->modelops.setinitialguessarray = 0;
  opflow->modelops.computeequalityconstraints = 0;
  opflow->modelops.computeequalityconstraintsarray = 0;
  opflow->modelops.computeinequalityconstraints = 0;
  opflow->modelops.computeinequalityconstraintsarray = 0;
  opflow->modelops.computeconstraints = 0;
  opflow->modelops.computeconstraintsarray = 0;
  opflow->modelops.computeequalityconstraintjacobian = 0;
  opflow->modelops.computeinequalityconstraintjacobian = 0;
  opflow->modelops.computehessian = 0;
  opflow->modelops.computeobjandgradient = 0;
  opflow->modelops.computeobjective = 0;
  opflow->modelops.computeobjectivearray = 0;
  opflow->modelops.computegradient = 0;
  opflow->modelops.computegradientarray = 0;
  opflow->modelops.computejacobian = 0;
  opflow->modelops.solutiontops = 0;
  opflow->modelops.computesparseequalityconstraintjacobianhiop = 0;
  opflow->modelops.computesparseinequalityconstraintjacobianhiop = 0;
  opflow->modelops.computesparsehessianhiop = 0;
  opflow->modelops.computedenseequalityconstraintjacobianhiop = 0;
  opflow->modelops.computedenseinequalityconstraintjacobianhiop = 0;
  opflow->modelops.computedensehessianhiop = 0;
  opflow->modelops.solutioncallbackhiop = 0;
  opflow->modelops.computeauxobjective = 0;
  opflow->modelops.computeauxgradient = 0;
  opflow->modelops.computeauxhessian = 0;

  opflow->modelname = modelname;
  /* Call the underlying implementation constructor */
  ierr = (*r)(opflow);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  OPFLOWGetModel - gets the model for OPFLOW

  Input Parameters:
+ opflow - opflow application object
- modelname - name of the model
*/
PetscErrorCode OPFLOWGetModel(OPFLOW opflow, std::string *modelname) {
  PetscFunctionBegin;
  *modelname = opflow->modelname;
  PetscFunctionReturn(0);
}

/*
  OPFLOWReadMatPowerData - Reads the network data given in MATPOWER data format

  Input Parameter
+  opflow - The OPFLOW object
-  netfile - The name of the network file

*/
PetscErrorCode OPFLOWReadMatPowerData(OPFLOW opflow, const char netfile[]) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* Read MatPower data file and populate the PS data structure */
  ierr = PSReadMatPowerData(opflow->ps, netfile);
  CHKERRQ(ierr);
  //  ierr = PetscPrintf(opflow->comm->type,"OPFLOW: Finished reading network
  //  data file %s\n",netfile);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWSetGICData - Sets the GIC data file

  Input Parameter
+  opflow - The OPFLOW object
-  gicfile - The name of the GIC data file

  Notes: The GIC data file is only used for visualization. It contains the
substation geospatial coordinates and mapping of buses to substations. See the
Electric Grid Data Repository files for examples of GIC data files (given for
CTIVSg cases)
*/
PetscErrorCode OPFLOWSetGICData(OPFLOW opflow, const char gicfile[]) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* Set GIC data file */
  ierr = PSSetGICData(opflow->ps, gicfile);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   OPFLOWSetNumConstraints - Sets the number of constraints for the OPFLOW
problem

   Input Parameters:
.  OPFLOW - the opflow application object

   Output Parameters:
+  busnconeq - number of equality constraints at each bus
.  branchnconeq - number of equality constraints at each branch
.  nconeq   -  local number of equality constraints
-  nconineq -  local number of inequality constraints

*/
PetscErrorCode OPFLOWSetNumConstraints(OPFLOW opflow, PetscInt *branchnconeq,
                                       PetscInt *busnconeq, PetscInt *nconeq,
                                       PetscInt *nconineq) {
  PetscErrorCode ierr;
  PetscInt i, vStart, vEnd, eStart, eEnd;
  DM networkdm = opflow->ps->networkdm;
  PS ps = opflow->ps;
  DM plexdm;
  PetscSection buseqconsection, buseqconglobsection;
  PetscSection varsection, varglobsection;
  PetscInt nconeqloc = 0;

  PetscFunctionBegin;
  ierr = (*opflow->modelops.setnumconstraints)(opflow, branchnconeq, busnconeq,
                                               nconeq, nconineq);

  ierr = PetscSectionCreate(opflow->comm->type, &buseqconsection);
  CHKERRQ(ierr);

  ierr = DMNetworkGetVertexRange(networkdm, &vStart, &vEnd);
  CHKERRQ(ierr);
  ierr = DMNetworkGetEdgeRange(networkdm, &eStart, &eEnd);
  CHKERRQ(ierr);

  ierr = PetscSectionSetChart(buseqconsection, eStart, vEnd);
  CHKERRQ(ierr);

  for (i = 0; i < ps->nline; i++) {
    ierr = PetscSectionSetDof(buseqconsection, eStart + i, 0);
  }

  for (i = 0; i < ps->nbus; i++) {
    ierr = PetscSectionSetDof(buseqconsection, vStart + i, busnconeq[i]);
    CHKERRQ(ierr);
    nconeqloc += busnconeq[i];
  }
  ierr = PetscSectionSetUp(buseqconsection);
  CHKERRQ(ierr);

  /* What we want to do here is have the DM set up the global section
     for the equality constraints (busnconeqglobsection) to get the global
     indices for equality constraints. We will use it later when operating
     with equality constraints.
     To do so, we need to
     i) Get the local and global variables section (varsection and
     globvarsection) from the DM. Increase the reference count for the sections
     so that they do not get destroyed when the new section is set in ii ii) Set
     the section for equality constraints (busnconeqsection) on the DM iii) Have
     the DM generate the Global section and save it to busnconeqglobsection iv)
     Copy over the global indices from busnconeqglobsection. v) Reset the
     variable sections on the DM and decrement the reference count otherwise
         there will be a memory leak.
  */
  /* (i) */
  ierr = DMNetworkGetPlex(networkdm, &plexdm);
  CHKERRQ(ierr);
  ierr = DMGetSection(plexdm, &varsection);
  CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)varsection);
  CHKERRQ(ierr);
  ierr = DMGetGlobalSection(plexdm, &varglobsection);
  CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)varglobsection);
  CHKERRQ(ierr);

  /*  (ii) */
  ierr = DMSetSection(plexdm, buseqconsection);
  CHKERRQ(ierr);

  /* (iii) */
  ierr = DMGetGlobalSection(plexdm, &buseqconglobsection);
  CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)buseqconglobsection);
  CHKERRQ(ierr);

  /* (iv) Set the global indices */
  ierr = PetscCalloc1(ps->nbus, &opflow->eqconglobloc);
  CHKERRQ(ierr);
  for (i = 0; i < ps->nbus; i++) {
    ierr = DMNetworkGetGlobalVecOffset(networkdm, vStart + i, ALL_COMPONENTS,
                                       &opflow->eqconglobloc[i]);
    CHKERRQ(ierr);
  }

  /* (v) */
  ierr = DMSetSection(plexdm, varsection);
  CHKERRQ(ierr);
  ierr = DMSetGlobalSection(plexdm, varglobsection);
  CHKERRQ(ierr);
  ierr = PetscObjectDereference((PetscObject)varsection);
  CHKERRQ(ierr);
  ierr = PetscObjectDereference((PetscObject)varglobsection);
  CHKERRQ(ierr);

  ierr = PetscSectionDestroy(&buseqconsection);
  CHKERRQ(ierr);
  ierr = PetscSectionDestroy(&buseqconglobsection);
  CHKERRQ(ierr);

  /* Create VecScatter for scattering values from global eq. constraints vector
   * to local eq. constraints vector. */
  PetscInt *eqconloc, *eqconglobloc, ctr = 0, j;
  ierr = PetscCalloc1(nconeqloc, &eqconloc);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(nconeqloc, &eqconglobloc);
  CHKERRQ(ierr);
  for (i = 0; i < ps->nbus; i++) {
    for (j = 0; j < busnconeq[i]; j++) {
      eqconloc[ctr] = ctr;
      eqconglobloc[ctr] = opflow->eqconglobloc[i] + j;
      ctr++;
    }
  }

  ierr = ISCreateGeneral(opflow->comm->type, nconeqloc, eqconloc,
                         PETSC_OWN_POINTER, &opflow->isconeqlocal);
  CHKERRQ(ierr);
  ierr = ISCreateGeneral(opflow->comm->type, nconeqloc, eqconglobloc,
                         PETSC_OWN_POINTER, &opflow->isconeqglob);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
   OPFLOWSetNumVariables - Sets the number of variables for the OPFLOW problem

   Input Parameters:
. opflow - OPFLOW application object

   Output Parameters:
+ busnvararray    - array of number of variables at each bus
. branchnvararray - array of number of variables at each branch
- nx              - total number of variables
*/
PetscErrorCode OPFLOWSetNumVariables(OPFLOW opflow, PetscInt *busnvararray,
                                     PetscInt *branchnvararray, PetscInt *nx) {
  PetscErrorCode ierr;
  PetscSection varsection, oldvarsection;
  PetscInt i, vStart, vEnd, eStart, eEnd;
  DM networkdm = opflow->ps->networkdm;
  PS ps = opflow->ps;
  DM plexdm;
  PetscSection globalvarsection;

  PetscFunctionBegin;

  ierr = PetscSectionCreate(opflow->comm->type, &varsection);
  CHKERRQ(ierr);

  ierr = DMNetworkGetVertexRange(networkdm, &vStart, &vEnd);
  CHKERRQ(ierr);
  ierr = DMNetworkGetEdgeRange(networkdm, &eStart, &eEnd);
  CHKERRQ(ierr);

  ierr = PetscSectionSetChart(varsection, eStart, vEnd);
  CHKERRQ(ierr);

  ierr = (*opflow->modelops.setnumvariables)(opflow, busnvararray,
                                             branchnvararray, nx);

  for (i = 0; i < ps->nline; i++) {
    ierr = PetscSectionSetDof(varsection, eStart + i, branchnvararray[i]);
  }

  for (i = 0; i < ps->nbus; i++) {
    ierr = PetscSectionSetDof(varsection, vStart + i, busnvararray[i]);
    CHKERRQ(ierr);
  }
  ierr = PetscSectionSetUp(varsection);
  CHKERRQ(ierr);

  ierr = DMNetworkGetPlex(networkdm, &plexdm);
  CHKERRQ(ierr);
  ierr = DMGetSection(plexdm, &oldvarsection);
  CHKERRQ(ierr);
  ierr = PetscSectionDestroy(&oldvarsection);
  CHKERRQ(ierr);
  /* This is hack to null the networkdm->DofSection pointer otherwise
     DMDestroy_Network throws an error. The issue is that though the
     section is destroyed, there is a dangling pointer for
     networkdm->DofSection. This hack needs to be fixed correctly via changes to
     DMNetwork
  */
  DM_Network *networkdmdata = (DM_Network *)(networkdm->data);
  networkdmdata->DofSection = 0;

  /* Set the section with number of variables */
  ierr = DMSetSection(plexdm, varsection);
  CHKERRQ(ierr);
  ierr = DMGetGlobalSection(plexdm, &globalvarsection);
  CHKERRQ(ierr);
  networkdmdata->DofSection = varsection;

  /* Update starting locations for variables at each line */
  for (i = 0; i < ps->nline; i++) {
    ierr = DMNetworkGetLocalVecOffset(networkdm, eStart + i, ALL_COMPONENTS,
                                      &ps->line[i].startloc);
    CHKERRQ(ierr);
    ierr = DMNetworkGetGlobalVecOffset(networkdm, eStart + i, ALL_COMPONENTS,
                                       &ps->line[i].startlocglob);
    CHKERRQ(ierr);
  }

  /* Update starting locations for variables at each bus */
  for (i = 0; i < ps->nbus; i++) {
    ierr = DMNetworkGetLocalVecOffset(networkdm, vStart + i, ALL_COMPONENTS,
                                      &ps->bus[i].startloc);
    CHKERRQ(ierr);
    ierr = DMNetworkGetGlobalVecOffset(networkdm, vStart + i, ALL_COMPONENTS,
                                       &ps->bus[i].startlocglob);
    CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/*
  OPFLOWSetUp - Sets up an optimal power flow application object

  Input Parameters:
. opflow - the OPFLOW object

  Notes:
  This routine sets up the OPFLOW object and the underlying PS object. It
  also distributes the PS object when used in parallel.
*/
PetscErrorCode OPFLOWSetUp(OPFLOW opflow) {
  PetscErrorCode ierr;
  char modelname[32], solvername[32];
  PetscBool modelset = PETSC_FALSE, solverset = PETSC_FALSE;
  PS ps = opflow->ps;
  PetscInt *branchnconeq, *busnconeq;
  PetscInt sendbuf[3], recvbuf[3], i;
  PetscInt nlinekvlevels = MAX_KV_LEVELS;

  PetscFunctionBegin;

  /* Read run-time options */
  if (!opflow->skip_options) {
    PetscOptionsBegin(opflow->comm->type, NULL, "OPFLOW options", NULL);

    ierr = PetscOptionsString(OPFLOWOptions::model.opt.c_str(),
                              OPFLOWOptions::model.desc.c_str(), "", modelname,
                              modelname, 32, &modelset);
    CHKERRQ(ierr);

    ierr = PetscOptionsString(OPFLOWOptions::solver.opt.c_str(),
                              OPFLOWOptions::solver.desc.c_str(), "",
                              solvername, solvername, 32, &solverset);
    CHKERRQ(ierr);
    ierr = PetscOptionsEnum(OPFLOWOptions::initialization.opt.c_str(),
                            OPFLOWOptions::initialization.desc.c_str(), "",
                            OPFLOWInitializationTypes,
                            (PetscEnum)opflow->initializationtype,
                            (PetscEnum *)&opflow->initializationtype, NULL);
    CHKERRQ(ierr);
    ierr =
        PetscOptionsEnum(OPFLOWOptions::objective.opt.c_str(),
                         OPFLOWOptions::objective.desc.c_str(), "",
                         OPFLOWObjectiveTypes, (PetscEnum)opflow->objectivetype,
                         (PetscEnum *)&opflow->objectivetype, NULL);
    CHKERRQ(ierr);
    ierr = PetscOptionsEnum(OPFLOWOptions::genbusvoltage.opt.c_str(),
                            OPFLOWOptions::genbusvoltage.desc.c_str(), "",
                            OPFLOWGenBusVoltageTypes,
                            (PetscEnum)opflow->genbusvoltagetype,
                            (PetscEnum *)&opflow->genbusvoltagetype, NULL);
    CHKERRQ(ierr);

    ierr = PetscOptionsEnum(OPFLOWOptions::outputformat.opt.c_str(),
                            OPFLOWOptions::outputformat.desc.c_str(), "",
                            OPFLOWOutputFormatTypes,
                            (PetscEnum)opflow->outputformat,
                            (PetscEnum *)&opflow->outputformat, NULL);
    CHKERRQ(ierr);

    ierr = PetscOptionsBool(OPFLOWOptions::has_gensetpoint.opt.c_str(),
                            OPFLOWOptions::has_gensetpoint.desc.c_str(), "",
                            opflow->has_gensetpoint, &opflow->has_gensetpoint,
                            NULL);
    CHKERRQ(ierr);
    ierr = PetscOptionsBool(OPFLOWOptions::use_agc.opt.c_str(),
                            OPFLOWOptions::use_agc.desc.c_str(), "",
                            opflow->use_agc, &opflow->use_agc, NULL);
    CHKERRQ(ierr);
    if (opflow->use_agc)
      opflow->has_gensetpoint = PETSC_TRUE;
    ierr = PetscOptionsReal(OPFLOWOptions::tolerance.opt.c_str(),
                            OPFLOWOptions::tolerance.desc.c_str(), "",
                            opflow->tolerance, &opflow->tolerance, NULL);
    CHKERRQ(ierr);

    if (opflow->objectivetype == MIN_GEN_COST) {
      opflow->obj_gencost = PETSC_TRUE;
    } else if (opflow->objectivetype == MIN_GENSETPOINT_DEVIATION) {
      opflow->obj_gencost = PETSC_FALSE;
      opflow->has_gensetpoint = PETSC_TRUE;
    }

    /* Ignore line flow constraints? */
    ierr = PetscOptionsBool(
        OPFLOWOptions::ignore_lineflow_constraints.opt.c_str(),
        OPFLOWOptions::ignore_lineflow_constraints.desc.c_str(), "",
        opflow->ignore_lineflow_constraints,
        &opflow->ignore_lineflow_constraints, NULL);
    CHKERRQ(ierr);

    /* Allow line flow violations ? */
    ierr =
        PetscOptionsBool(OPFLOWOptions::allow_lineflow_violation.opt.c_str(),
                         OPFLOWOptions::allow_lineflow_violation.desc.c_str(),
                         "", opflow->allow_lineflow_violation,
                         &opflow->allow_lineflow_violation, NULL);
    CHKERRQ(ierr);

    /* Line flow penalty */
    ierr =
        PetscOptionsReal(OPFLOWOptions::lineflowviolation_penalty.opt.c_str(),
                         OPFLOWOptions::lineflowviolation_penalty.desc.c_str(),
                         "", opflow->lineflowviolation_penalty,
                         &opflow->lineflowviolation_penalty, NULL);

    ierr =
        PetscOptionsBool(OPFLOWOptions::include_loadloss_variables.opt.c_str(),
                         OPFLOWOptions::include_loadloss_variables.desc.c_str(),
                         "", opflow->include_loadloss_variables,
                         &opflow->include_loadloss_variables, NULL);
    CHKERRQ(ierr);
    ierr = PetscOptionsReal(OPFLOWOptions::loadloss_penalty.opt.c_str(),
                            OPFLOWOptions::loadloss_penalty.desc.c_str(), "",
                            opflow->loadloss_penalty, &opflow->loadloss_penalty,
                            NULL);
    CHKERRQ(ierr);
    ierr = PetscOptionsBool(
        OPFLOWOptions::include_powerimbalance_variables.opt.c_str(),
        OPFLOWOptions::include_powerimbalance_variables.desc.c_str(), "",
        opflow->include_powerimbalance_variables,
        &opflow->include_powerimbalance_variables, NULL);
    CHKERRQ(ierr);
    ierr = PetscOptionsReal(OPFLOWOptions::powerimbalance_penalty.opt.c_str(),
                            OPFLOWOptions::powerimbalance_penalty.desc.c_str(),
                            "", opflow->powerimbalance_penalty,
                            &opflow->powerimbalance_penalty, NULL);
    CHKERRQ(ierr);

    if (!opflow->ignore_lineflow_constraints) {
      ierr = PetscOptionsScalarArray("-opflow_monitor_line_kvlevels",
                                     "KV levels for lines to monitor", "",
                                     opflow->linekvmon, &nlinekvlevels, NULL);
      CHKERRQ(ierr);
      opflow->nlinekvmon = nlinekvlevels;
      if (!opflow->nlinekvmon)
        opflow->nlinekvmon = -1; /* Value not set by option so use -1 to
                                  indicate to select all kvlevels */
    }
    PetscOptionsEnd();
  }

  /* Set up underlying PS object */
  ierr = OPFLOWSetUpPS(opflow);
  CHKERRQ(ierr);

  /* Set model if CLI argument */
  if (modelset) {
    if (opflow->model)
      ierr = (*opflow->modelops.destroy)(opflow);
    ierr = OPFLOWSetModel(opflow, modelname);
    CHKERRQ(ierr);
    //    ierr = PetscPrintf(opflow->comm->type,"OPFLOW: Using %s
    //    model\n",modelname);CHKERRQ(ierr);
  } else {
    // Specify default arguments if model not set
    if (!opflow->model) {
      ierr = OPFLOWSetModel(opflow, OPFLOWOptions::model.default_value.c_str());
      CHKERRQ(ierr);
    }
  }

  /* Set solver if CLI argument */
  if (solverset) {
    if (opflow->solver)
      ierr = (*opflow->solverops.destroy)(opflow);
    ierr = OPFLOWSetSolver(opflow, solvername);
    CHKERRQ(ierr);
    //    ierr = PetscPrintf(opflow->comm->type,"OPFLOW: Using %s
    //    solver\n",solvername);CHKERRQ(ierr);
  } else {
    // Specify default arguments if solver not set
    if (!opflow->solver) {

      // Default values may have OPFLOWSOLVER_ prepended to the string ExaGO
      // uses internally, so we strip it before assigning (if found)
      auto default_solver =
          std::string(OPFLOWOptions::solver.default_value.c_str());
      auto underscore = default_solver.find("_");

      if (underscore != std::string::npos)
        default_solver.erase(0, ++underscore);

      ierr = OPFLOWSetSolver(opflow, default_solver.c_str());
      CHKERRQ(ierr);
    }
  }
  /* Once model and solver are set, we should check for compatibility */
  ierr = OPFLOWCheckModelSolverCompatibility(opflow);
  CHKERRQ(ierr);

  /* Get list of monitored lines. These will be included in
     the inequality constraints for line flows
  */
  if (!opflow->ignore_lineflow_constraints) {
    ierr = OPFLOWGetLinesMonitored(opflow);
    CHKERRQ(ierr);
  }
  if (opflow->include_loadloss_variables) {
    PS ps;
    PSBUS bus;
    PSLOAD load;

    ierr = OPFLOWGetPS(opflow, &ps);
    CHKERRQ(ierr);

    for (int i = 0; i < ps->nbus; i++) {
      bus = &(ps->bus[i]);
      for (int l = 0; l < bus->nload; l++) {
        ierr = PSBUSGetLoad(bus, l, &load);
        CHKERRQ(ierr);
        if (load->loss_cost == BOGUSLOSSCOST) {
          load->loss_cost = opflow->loadloss_penalty;
        }
      }
    }
  }

  ierr = PetscCalloc1(ps->nbus, &opflow->busnvararray);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->nline, &opflow->branchnvararray);
  CHKERRQ(ierr);

  /* Set up number of variables for branches and buses */
  ierr = OPFLOWSetNumVariables(opflow, opflow->busnvararray,
                               opflow->branchnvararray, &opflow->nx);
  CHKERRQ(ierr);

  ierr = PetscCalloc1(ps->nline, &branchnconeq);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->nbus, &busnconeq);
  CHKERRQ(ierr);
  /* Set up number of equality and inequality constraints and
     number of equality constraints at each bus */

  /* Set number of constraints */
  ierr = OPFLOWSetNumConstraints(opflow, branchnconeq, busnconeq,
                                 &opflow->nconeq, &opflow->nconineq);
  CHKERRQ(ierr);
  opflow->ncon = opflow->nconeq + opflow->nconineq;

  ierr = PetscFree(branchnconeq);
  CHKERRQ(ierr);
  ierr = PetscFree(busnconeq);
  CHKERRQ(ierr);

  sendbuf[0] = opflow->nx;
  sendbuf[1] = opflow->nconeq;
  sendbuf[2] = opflow->nconineq;
  ierr =
      MPI_Allreduce(sendbuf, recvbuf, 3, MPIU_INT, MPI_SUM, opflow->comm->type);
  CHKERRQ(ierr);
  opflow->Nx = recvbuf[0];
  opflow->Nconeq = recvbuf[1];
  opflow->Nconineq = recvbuf[2];
  opflow->Ncon = opflow->Nconeq + opflow->Nconineq;
  //  ierr = PetscPrintf(PETSC_COMM_SELF,"OPFLOW: Rank %d: nx = %d nconeq = %d,
  //  nconineq = %d, ncon =
  //  %d\n",opflow->comm->rank,opflow->nx,opflow->nconeq,opflow->nconineq,opflow->ncon);CHKERRQ(ierr);

  /* Set vertex local to global ordering */
  ierr = DMNetworkSetVertexLocalToGlobalOrdering(opflow->ps->networkdm);
  CHKERRQ(ierr);

  /* Create solution vector and upper/lower bounds */
  ierr = PSCreateGlobalVector(opflow->ps, &opflow->X);
  CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X, &opflow->gradobj);
  CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X, &opflow->Xl);
  CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X, &opflow->Xu);
  CHKERRQ(ierr);
  /* Default variable bounds are -inf <= x <= +inf */
  ierr = VecSet(opflow->Xl, PETSC_NINFINITY);
  CHKERRQ(ierr);
  ierr = VecSet(opflow->Xu, PETSC_INFINITY);
  CHKERRQ(ierr);

  /* Vectors for constraints */
  ierr = VecCreate(ps->comm->type, &opflow->G);
  CHKERRQ(ierr);
  ierr = VecSetSizes(opflow->G, opflow->ncon, opflow->Ncon);
  CHKERRQ(ierr);
  ierr = VecSetFromOptions(opflow->G);
  CHKERRQ(ierr);

  ierr = VecDuplicate(opflow->G, &opflow->Gl);
  CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->G, &opflow->Gu);
  CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->G, &opflow->Lambda);
  CHKERRQ(ierr);

  ierr = VecCreate(ps->comm->type, &opflow->Ge);
  CHKERRQ(ierr);
  ierr = VecSetSizes(opflow->Ge, opflow->nconeq, opflow->Nconeq);
  CHKERRQ(ierr);
  ierr = VecSetFromOptions(opflow->Ge);
  CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->Ge, &opflow->Lambdae);
  CHKERRQ(ierr);

  PetscInt nconeqlocal;

  ierr = ISGetSize(opflow->isconeqlocal, &nconeqlocal);
  CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_SELF, &opflow->Gelocal);
  CHKERRQ(ierr);
  ierr = VecSetSizes(opflow->Gelocal, nconeqlocal, nconeqlocal);
  CHKERRQ(ierr);
  ierr = VecSetFromOptions(opflow->Gelocal);
  CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->Gelocal, &opflow->Lambdaelocal);
  CHKERRQ(ierr);

  ierr = VecScatterCreate(opflow->Ge, opflow->isconeqglob, opflow->Gelocal,
                          opflow->isconeqlocal, &opflow->scattereqcon);
  CHKERRQ(ierr);

  if (opflow->Nconineq) {
    ierr = VecCreate(ps->comm->type, &opflow->Gi);
    CHKERRQ(ierr);
    ierr = VecSetSizes(opflow->Gi, opflow->nconineq, opflow->Nconineq);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(opflow->Gi);
    CHKERRQ(ierr);
    ierr = VecDuplicate(opflow->Gi, &opflow->Lambdai);
    CHKERRQ(ierr);
  }

  /* Create equality and inequality constraint Jacobian matrices */
  ierr = MatCreate(opflow->comm->type, &opflow->Jac_Ge);
  CHKERRQ(ierr);
  ierr = MatSetSizes(opflow->Jac_Ge, opflow->nconeq, opflow->nx, opflow->Nconeq,
                     opflow->Nx);
  CHKERRQ(ierr);
  ierr = MatSetFromOptions(opflow->Jac_Ge);
  CHKERRQ(ierr);
  /* Assume 10% sparsity */
  PetscInt nzrow;
  if (opflow->nx < 1000) {
    /* Small case, assume 10% sparsity */
    nzrow = (PetscInt)(0.1 * opflow->nx);
  } else {
    /* Bigger case, assume 2% sparsity */
    nzrow = (PetscInt)(0.02 * opflow->nx);
  }
  ierr = MatSeqAIJSetPreallocation(opflow->Jac_Ge, nzrow, NULL);
  CHKERRQ(ierr);
  ierr =
      MatSetOption(opflow->Jac_Ge, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  CHKERRQ(ierr);

  if (opflow->Nconineq) {
    ierr = MatCreate(opflow->comm->type, &opflow->Jac_Gi);
    CHKERRQ(ierr);
    ierr = MatSetSizes(opflow->Jac_Gi, opflow->nconineq, opflow->nx,
                       opflow->Nconineq, opflow->Nx);
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(opflow->Jac_Gi);
    CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(opflow->Jac_Gi, 6, NULL);
    CHKERRQ(ierr);
    ierr = MatSetOption(opflow->Jac_Gi, MAT_NEW_NONZERO_ALLOCATION_ERR,
                        PETSC_FALSE);
    CHKERRQ(ierr);
  }

  /* Create Hessian */
  ierr = PSCreateMatrix(opflow->ps, &opflow->Hes);
  CHKERRQ(ierr);

  /* Create natural to sparse dense ordering mapping (needed for some models)
     Here, we create the mapping array and fill up natural ordering. The model
     set up function should change the mapping to what it needs
  */
  ierr = PetscMalloc1(opflow->nx, &opflow->idxn2sd_map);
  CHKERRQ(ierr);
  for (i = 0; i < opflow->nx; i++)
    opflow->idxn2sd_map[i] = i;

  /* Model set up */
  if (opflow->modelops.setup) {
    ierr = (*opflow->modelops.setup)(opflow);
    CHKERRQ(ierr);
  }

  /* Set bounds on variables */
  ierr = OPFLOWComputeVariableBounds(opflow, opflow->Xl, opflow->Xu);
  CHKERRQ(ierr);

  /* Set bounds on constraints */
  ierr = OPFLOWComputeConstraintBounds(opflow, opflow->Gl, opflow->Gu);
  CHKERRQ(ierr);

  /* Set initial guess */
  ierr = OPFLOWSetInitialGuess(opflow, opflow->X, opflow->Lambda);
  CHKERRQ(ierr);

  /* Initial guess for multipliers */
  /*  ierr = VecSet(opflow->Lambda, 1.0);
  CHKERRQ(ierr);
  ierr = VecSet(opflow->Lambdae, 1.0);
  CHKERRQ(ierr);
  if (opflow->Nconineq) {
    ierr = VecSet(opflow->Lambdai, 1.0);
    CHKERRQ(ierr);
  }
  */

  /* Solver set up */
  if (opflow->solverops.setup) {
    ierr = (*opflow->solverops.setup)(opflow);
    CHKERRQ(ierr);
  }
  //  ierr = PetscPrintf(opflow->comm->type,"OPFLOW: Setup
  //  completed\n");CHKERRQ(ierr);

  /* Register events for logging */
  ierr = PetscLogEventRegister("OPFLOWObj", 0, &opflow->objlogger);
  CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OPFLOWGrad", 0, &opflow->gradlogger);
  CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OPFLOWEqCons", 0, &opflow->eqconslogger);
  CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OPFLOWIneqCons", 0, &opflow->ineqconslogger);
  CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OPFLOWEqConsJac", 0, &opflow->eqconsjaclogger);
  CHKERRQ(ierr);
  ierr =
      PetscLogEventRegister("OPFLOWIneqConsJac", 0, &opflow->ineqconsjaclogger);
  CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OPFLOWHess", 0, &opflow->hesslogger);
  CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OPFLOWSolve", 0, &opflow->solvelogger);
  CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OPFLOWDenseHess", 0, &opflow->densehesslogger);
  CHKERRQ(ierr);
  ierr =
      PetscLogEventRegister("OPFLOWSparseHess", 0, &opflow->sparsehesslogger);
  CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OPFLOWDenseIneqConsJac", 0,
                               &opflow->denseineqconsjaclogger);
  CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OPFLOWDenseEqConsJac", 0,
                               &opflow->denseeqconsjaclogger);
  CHKERRQ(ierr);

  ierr = PetscLogEventRegister("OPFLOWSaveSolution", 0, &opflow->outputlogger);
  CHKERRQ(ierr);

  opflow->solve_real_time = 0.0;
  opflow->solve_cpu_time = 0.0;

  /* Compute area participation factors */
  ierr = PSComputeParticipationFactors(ps);
  CHKERRQ(ierr);
  opflow->setupcalled = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/* OPFLOWSetInitialGuess - Sets the initial guess for OPFLOW

   Input Parameters:
.  opflow - the optimal power flow application object

   Output Parameters:
+  X - initial guess (primal)
-  Lambda - Lagrange multipliers for constraints

*/
PetscErrorCode OPFLOWSetInitialGuess(OPFLOW opflow, Vec X, Vec Lambda) {
  PetscErrorCode ierr;
  PetscBool converged = PETSC_FALSE; // Solver convergence status

  PetscFunctionBegin;

  /* Set initial guess */
  switch (opflow->initializationtype) {
  case OPFLOWINIT_MIDPOINT:
  case OPFLOWINIT_FROMFILE:
  case OPFLOWINIT_FLATSTART:
    if (opflow->modelops.setinitialguess) {
      ierr = (*opflow->modelops.setinitialguess)(opflow, X, Lambda);
      CHKERRQ(ierr);
    } else if (opflow->modelops.setinitialguessarray) {
      PetscScalar *x0;
      ierr = VecGetArray(X, &x0);
      CHKERRQ(ierr);
      ierr = (*opflow->modelops.setinitialguessarray)(opflow, x0);
      CHKERRQ(ierr);
      ierr = VecRestoreArray(X, &x0);
      CHKERRQ(ierr);
    }
    break;
  case OPFLOWINIT_ACPF:
    ierr = OPFLOWSetUpInitPflow(opflow);              // Set up power flow
    ierr = OPFLOWComputePrePflow(opflow, &converged); // Solve power flow
    CHKERRQ(ierr);
    if (!converged) {
      SETERRQ(PETSC_COMM_SELF, 0,
              "AC power flow initialization did not converge\n");
    }
    if (opflow->modelops.setinitialguess) {
      ierr = (*opflow->modelops.setinitialguess)(opflow, X, Lambda);
      CHKERRQ(ierr);
    } else if (opflow->modelops.setinitialguessarray) {
      PetscScalar *x0;
      ierr = VecGetArray(X, &x0);
      CHKERRQ(ierr);
      ierr = (*opflow->modelops.setinitialguessarray)(opflow, x0);
      CHKERRQ(ierr);
      ierr = VecRestoreArray(X, &x0);
      CHKERRQ(ierr);
    }
    break;
  case OPFLOWINIT_DCOPF:
    ierr = OPFLOWInitializeDCOPF(opflow, &converged); // DCOPF initialization
    CHKERRQ(ierr);
    if (!converged) {
      SETERRQ(PETSC_COMM_SELF, 0, "DCOPF initialization did not converge\n");
    }
    if (opflow->modelops.setinitialguess) {
      ierr = (*opflow->modelops.setinitialguess)(opflow, X, Lambda);
      CHKERRQ(ierr);
    } else if (opflow->modelops.setinitialguessarray) {
      PetscScalar *x0;
      ierr = VecGetArray(X, &x0);
      CHKERRQ(ierr);
      ierr = (*opflow->modelops.setinitialguessarray)(opflow, x0);
      CHKERRQ(ierr);
      ierr = VecRestoreArray(X, &x0);
      CHKERRQ(ierr);
    }
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
            "Unknown OPFLOW initialization type\n");
  }

  PetscFunctionReturn(0);
}

/*
  OPFLOWSolve - Solves the AC optimal power flow

  Input Parameters:
. opflow - the optimal power flow application object
*/
PetscErrorCode OPFLOWSolve(OPFLOW opflow) {
  PetscErrorCode ierr;
  PetscLogDouble real1 = 0.0, real2 = 0.0;
  PetscLogDouble cpu1 = 0.0, cpu2 = 0.0;
  PetscFunctionBegin;

  ierr = PetscTime(&real1);
  CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&cpu1);
  CHKERRQ(ierr);

  if (!opflow->setupcalled) {
    ierr = OPFLOWSetUp(opflow);
  }

  /* Solve */
  ierr = PetscLogEventBegin(opflow->solvelogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = (*opflow->solverops.solve)(opflow);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->solvelogger, 0, 0, 0, 0);
  CHKERRQ(ierr);

  //  ierr = VecView(opflow->X,0);CHKERRQ(ierr);

  ierr = PetscTime(&real2);
  CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&cpu2);
  CHKERRQ(ierr);

  opflow->solve_real_time = real2 - real1;
  opflow->solve_cpu_time = cpu2 - cpu1;

  PetscFunctionReturn(0);
}

/*
  OPFLOWGetObjective - Returns the objective function value

  Input Parameters:
- OPFLOW - the OPFLOW object

  Output Parameters:
- obj    - the objective function value

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode OPFLOWGetObjective(OPFLOW opflow, PetscReal *obj) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*opflow->solverops.getobjective)(opflow, obj);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetVariableBounds - Returns the variable bounds

  Input Parameters:
- OPFLOW - the OPFLOW object

  Output Parameters:
+ Xl     - lower bound on X
- Xu     - upper bound on X

  Notes: Should be called after the optimization finishes. This function merely
         returns the already computed lower and upper bound vectors.
*/
PetscErrorCode OPFLOWGetVariableBounds(OPFLOW opflow, Vec *Xl, Vec *Xu) {
  PetscFunctionBegin;
  *Xl = opflow->Xl;
  *Xu = opflow->Xu;
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeObjective - Computes the objective function value

  Input Parameters:
+ OPFLOW - the OPFLOW object
- X      - the solution vector

  Output Parameters:
- obj    - the objective function value

*/
PetscErrorCode OPFLOWComputeObjective(OPFLOW opflow, Vec X, PetscReal *obj) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->objlogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeobjective)(opflow, X, obj);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->objlogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeObjectiveArray - Computes the objective function value

  Input Parameters:
+ OPFLOW - the OPFLOW object
- x      - the solution array

  Output Parameters:
- obj    - the objective function value

*/
PetscErrorCode OPFLOWComputeObjectiveArray(OPFLOW opflow, const double *x,
                                           double *obj) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->objlogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeobjectivearray)(opflow, x, obj);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->objlogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeVariableBounds - Computes the bounds on X

  Input Parameters:
- OPFLOW - the OPFLOW object

  Output Parameters:
+ Xl     - lower bound on X
- Xu     - upper bound on X

  Notes: Should be called after the optimization is set up.
*/
PetscErrorCode OPFLOWComputeVariableBounds(OPFLOW opflow, Vec Xl, Vec Xu) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (opflow->modelops.setvariablebounds) {
    ierr = (*opflow->modelops.setvariablebounds)(opflow, Xl, Xu);
    CHKERRQ(ierr);
  } else if (opflow->modelops.setvariableboundsarray) {
    PetscScalar *xl, *xu;
    ierr = VecGetArray(Xl, &xl);
    CHKERRQ(ierr);
    ierr = VecGetArray(Xu, &xu);
    CHKERRQ(ierr);
    ierr = (*opflow->modelops.setvariableboundsarray)(opflow, xl, xu);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(Xl, &xl);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(Xu, &xu);
    CHKERRQ(ierr);
  }

  if (opflow->modelops.updatevariablebounds) {
    ierr = (*opflow->modelops.updatevariablebounds)(opflow, Xl, Xu,
                                                    opflow->userctx);
    CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeGradient - Computes the gradient of the objective function

  Input Parameters:
+ OPFLOW - the OPFLOW object
. X      - the solution vector
- grad    - the gradient of the objective function

  Notes: Should be called after the optimization is set up
*/
PetscErrorCode OPFLOWComputeGradient(OPFLOW opflow, Vec X, Vec grad) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->gradlogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computegradient)(opflow, X, grad);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->gradlogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeGradientArray - Computes the gradient of the objective function
(array version)

  Input Parameters:
+ OPFLOW - the OPFLOW object
. x      - the solution array
- grad    - the gradient of the objective function

  Notes: Should be called after the optimization is set up
*/
PetscErrorCode OPFLOWComputeGradientArray(OPFLOW opflow, const double *x,
                                          double *grad) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->gradlogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computegradientarray)(opflow, x, grad);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->gradlogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetSolution - Returns the OPFLOW solution

  Input Parameters:
+ OPFLOW - the OPFLOW object
- X      - the opflow solution

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode OPFLOWGetSolution(OPFLOW opflow, Vec *X) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*opflow->solverops.getsolution)(opflow, X);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetConstraints - Returns the OPFLOW constraints

  Input Parameters:
+ OPFLOW - the OPFLOW object

  Output Parameters:
- G    - the opflow constraints

  Notes: Should be called after the optimization finishes.
         Equality constraints first followed by inequality constraints
*/
PetscErrorCode OPFLOWGetConstraints(OPFLOW opflow, Vec *G) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*opflow->solverops.getconstraints)(opflow, G);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeEqualityConstraints - Computes OPFLOW equality constraints
  Input Parameters:
+ OPFLOW - the OPFLOW object
. X      - the solution vector

  Output Parameters:
- Ge      - OPFLOW equallity constraints

  Notes: Should be called after the optimization is set up.

*/
PetscErrorCode OPFLOWComputeEqualityConstraints(OPFLOW opflow, Vec X, Vec Ge) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->eqconslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeequalityconstraints)(opflow, X, Ge);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->eqconslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeInequalityConstraints - Computes OPFLOW inequality constraints
  Input Parameters:
+ OPFLOW - the OPFLOW object
. X      - the solution vector

  Output Parameters:
- Gi      - OPFLOW inequallity constraints

  Notes: Should be called after the optimization is set up.
*/
PetscErrorCode OPFLOWComputeInequalityConstraints(OPFLOW opflow, Vec X,
                                                  Vec Gi) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->ineqconslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeinequalityconstraints)(opflow, X, Gi);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->ineqconslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeEqualityConstraintsArray - Computes OPFLOW equality constraints
(array version)

  Input Parameters:
+ OPFLOW - the OPFLOW object
. x      - the solution array

  Output Parameters:
- ge      - OPFLOW equality constraints array

  Notes: Should be called after the optimization is set up.
*/
PetscErrorCode OPFLOWComputeEqualityConstraintsArray(OPFLOW opflow,
                                                     const double *x,
                                                     double *ge) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->eqconslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeequalityconstraintsarray)(opflow, x, ge);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->eqconslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeInequalityConstraintsArray - Computes OPFLOW inequality
constraints (array version)

  Input Parameters:
+ OPFLOW - the OPFLOW object
. x      - the solution array

  Output Parameters:
- gi      - OPFLOW inequality constraints array

  Notes: Should be called after the optimization is set up.
*/
PetscErrorCode OPFLOWComputeInequalityConstraintsArray(OPFLOW opflow,
                                                       const double *x,
                                                       double *gi) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->ineqconslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeinequalityconstraintsarray)(opflow, x, gi);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->ineqconslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeConstraints - Computes the OPFLOW constraints

  Input Parameters:
+ OPFLOW - the OPFLOW object
. X      - the solution vector

  Output Parameters:
- G    - the opflow constraints

  Notes: Should be called after the optimization is set up. G has
         equality constraints followed by inequality constraints
*/
PetscErrorCode OPFLOWComputeConstraints(OPFLOW opflow, Vec X, Vec G) {
  PetscErrorCode ierr;
  PetscScalar *g;

  PetscFunctionBegin;
  ierr = VecGetArray(G, &g);
  CHKERRQ(ierr);

  ierr = VecPlaceArray(opflow->Ge, g);
  CHKERRQ(ierr);
  ierr = OPFLOWComputeEqualityConstraints(opflow, X, opflow->Ge);
  CHKERRQ(ierr);
  ierr = VecResetArray(opflow->Ge);
  CHKERRQ(ierr);

  if (opflow->Nconineq) {
    ierr = VecPlaceArray(opflow->Gi, g + opflow->nconeq);
    CHKERRQ(ierr);
    ierr = OPFLOWComputeInequalityConstraints(opflow, X, opflow->Gi);
    CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Gi);
    CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(G, &g);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeHessian - Computes the OPFLOW Hessian

  Input Parameters:
+ OPFLOW - the OPFLOW object
. X      - the solution vector
. Lambda - Lagrange multiplier
- obj_factor - scaling factor for objective part of Hessian (used by IPOPT)

  Output Parameters:
- H    - the Hessian

  Notes: Should be called after the optimization is set up. Lambda has
multipliers - equality constraints followed by inequality constraints
*/
PetscErrorCode OPFLOWComputeHessian(OPFLOW opflow, Vec X, Vec Lambda,
                                    PetscScalar obj_factor, Mat H) {
  PetscErrorCode ierr;
  PetscScalar *lambda;

  PetscFunctionBegin;

  opflow->obj_factor = obj_factor;

  /* IPOPT passes the scaled Lagrangian multipliers for Hessian calculation. The
     scaling factor is the Hessian objective factor it uses in the Hessian
     calculation
  */
  ierr = VecScale(Lambda, obj_factor);
  CHKERRQ(ierr);

  ierr = VecGetArray(Lambda, &lambda);
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
  ierr = (*opflow->modelops.computehessian)(opflow, X, opflow->Lambdae,
                                            opflow->Lambdai, H);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->hesslogger, 0, 0, 0, 0);
  CHKERRQ(ierr);

  ierr = VecResetArray(opflow->Lambdae);
  CHKERRQ(ierr);

  if (opflow->Nconineq) {
    ierr = VecResetArray(opflow->Lambdai);
    CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(Lambda, &lambda);
  CHKERRQ(ierr);

  /* Rescale Lambda back to its original value */
  ierr = VecScale(Lambda, 1. / obj_factor);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeConstraintJacobian - Computes the OPFLOW constraint Jacobians

  Input Parameters:
+ OPFLOW - the OPFLOW object
. X      - the solution vector

  Output Parameters:
+ Jeq    - equality constraint Jacobian
- Jineq  - Inequality constraint Jacobian

  Notes: Should be called after the optimization is set up.
*/
PetscErrorCode OPFLOWComputeConstraintJacobian(OPFLOW opflow, Vec X, Mat Jeq,
                                               Mat Jineq) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MatZeroEntries(Jeq);
  CHKERRQ(ierr);
  ierr = PetscLogEventBegin(opflow->eqconsjaclogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeequalityconstraintjacobian)(opflow, X, Jeq);
  CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->eqconsjaclogger, 0, 0, 0, 0);
  CHKERRQ(ierr);

  if (opflow->Nconineq) {
    ierr = MatZeroEntries(Jineq);
    CHKERRQ(ierr);
    ierr = PetscLogEventBegin(opflow->ineqconsjaclogger, 0, 0, 0, 0);
    CHKERRQ(ierr);
    ierr = (*opflow->modelops.computeinequalityconstraintjacobian)(opflow, X,
                                                                   Jineq);
    CHKERRQ(ierr);
    ierr = PetscLogEventEnd(opflow->ineqconsjaclogger, 0, 0, 0, 0);
    CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/*
  OPFLOWGetConstraintBounds - Returns the OPFLOW constraint bounds

  Input Parameters:
+ OPFLOW - the OPFLOW object

  Output Parameters:
+ Gl    - lower bound on constraints
- Gu    - upper bound on constraints

  Notes: Should be called after the optimization finishes.
         Equality constraint bounds first followed by inequality constraints
*/
PetscErrorCode OPFLOWGetConstraintBounds(OPFLOW opflow, Vec *Gl, Vec *Gu) {

  PetscFunctionBegin;
  *Gl = opflow->Gl;
  *Gu = opflow->Gu;
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetConstraintJacobian - Returns the OPFLOW equality and inequality
constraint jacobian

  Input Parameters:
+ OPFLOW - the OPFLOW object

  Output Parameters:
+ Jeq    - equality constraint Jacobian
- Jineq  - inequality constraint Jacobian

  Notes: Should be called after the optimization finishes.
         Jineq is set to NULL if no inequality constraints exists
*/
PetscErrorCode OPFLOWGetConstraintJacobian(OPFLOW opflow, Mat *Jeq,
                                           Mat *Jineq) {

  PetscFunctionBegin;

  *Jeq = opflow->Jac_Ge;
  if (opflow->nconineq)
    *Jineq = opflow->Jac_Gi;
  else
    *Jineq = NULL;

  PetscFunctionReturn(0);
}

/*
  OPFLOWGetHessian - Returns the OPFLOW Hessian

  Input Parameters:
+ OPFLOW - the OPFLOW object

  Output Parameters:
+ Hess    - Hessian matrix
- obj_factor - factor used in scaling objective Hessian (used by IPOPT)
  Notes: Should be called after the optimization finishes.
*/
PetscErrorCode OPFLOWGetHessian(OPFLOW opflow, Mat *Hess,
                                PetscScalar *obj_factor) {

  PetscFunctionBegin;

  *Hess = opflow->Hes;
  *obj_factor = opflow->obj_factor;

  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeConstraintBounds - Computes the bounds on constraints

  Input Parameters:
- OPFLOW - the OPFLOW object

  Output Parameters:
+ Gl     - lower bound on constraints
- Gu     - upper bound on constraints

  Notes: Should be called after the optimization is set up.
*/
PetscErrorCode OPFLOWComputeConstraintBounds(OPFLOW opflow, Vec Gl, Vec Gu) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (opflow->modelops.setconstraintbounds) {
    ierr = (*opflow->modelops.setconstraintbounds)(opflow, Gl, Gu);
    CHKERRQ(ierr);
  } else if (opflow->modelops.setconstraintboundsarray) {
    PetscScalar *gl, *gu;
    ierr = VecGetArray(Gl, &gl);
    CHKERRQ(ierr);
    ierr = VecGetArray(Gu, &gu);
    CHKERRQ(ierr);
    ierr = (*opflow->modelops.setconstraintboundsarray)(opflow, gl, gu);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(Gl, &gl);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(Gu, &gu);
    CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetConstraintMultipliers - Returns the OPFLOW constraint multipliers

  Input Parameters:
+ OPFLOW - the OPFLOW object

  Output Parameters:
- G    - the opflow constraint lagrange multipliers

  Notes: Should be called after the optimization finishes.
    Equality constraint multipliers first followed by inequality constraint
multipliers
*/
PetscErrorCode OPFLOWGetConstraintMultipliers(OPFLOW opflow, Vec *Lambda) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*opflow->solverops.getconstraintmultipliers)(opflow, Lambda);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetConvergenceStatus - Did OPFLOW converge?

  Input Parameters:
+ OPFLOW - the OPFLOW object
- status - PETSC_TRUE if converged, PETSC_FALSE otherwise

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode OPFLOWGetConvergenceStatus(OPFLOW opflow, PetscBool *status) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*opflow->solverops.getconvergencestatus)(opflow, status);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWSolutionToPS - Updates the PS struct from OPFLOW solution

  Input Parameters:
. opflow - the OPFLOW object

  Notes: Updates the different fields in the PS struct from the OPFLOW solution
*/
PetscErrorCode OPFLOWSolutionToPS(OPFLOW opflow) {
  PetscBool conv_status;
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = (*opflow->modelops.solutiontops)(opflow);

  /* Save the objective to PS */
  opflow->ps->opflowobj = opflow->obj;

  /* Save convergence status to PS */
  ierr = OPFLOWGetConvergenceStatus(opflow, &conv_status);
  opflow->ps->opflow_converged = conv_status;

  ierr = OPFLOWSetSummaryStats(opflow);

  opflow->solutiontops = PETSC_TRUE;
  PetscFunctionReturn(ierr);
}

/*
  OPFLOWSetObjectiveType - Sets the objective function for OPFLOW

  Input Parameters
+ opflow - the opflow object
- objtype - type of objective function

  Command-line option: -opflow_objective <obj_type>

  Notes: Must be called before OPFLOWSetUp
*/
PetscErrorCode OPFLOWSetObjectiveType(OPFLOW opflow,
                                      OPFLOWObjectiveType objtype) {
  PetscFunctionBegin;
  opflow->objectivetype = objtype;
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetObjectiveType - Gets the objective function for OPFLOW

  Input Parameters
+ opflow - the opflow object
- objtype - type of objective function

  Notes: Must be called after OPFLOWSetUp
*/
PetscErrorCode OPFLOWGetObjectiveType(OPFLOW opflow,
                                      OPFLOWObjectiveType *objtype) {
  PetscFunctionBegin;
  *objtype = opflow->objectivetype;
  PetscFunctionReturn(0);
}

/* OPFLOWGetPS - Gets the underlying PS object

  Input Parameters
. opflow - the OPFLOW object

  Output Parameters
. ps - the ps object

  Notes: This function returns the PS object that holds the network data. Using
the PS object one can make changes to the network parameters. A typical case is
         changing some network parameters before solving opflow.

         OPFLOWSetUpPS() must be called before OPFLOWGetPS()
*/
PetscErrorCode OPFLOWGetPS(OPFLOW opflow, PS *ps) {
  PetscFunctionBegin;
  if (!opflow->ps->setupcalled)
    SETERRQ(PETSC_COMM_SELF, 0,
            "OPFLOWSetUpPS() must be called before calling OPFLOWGetPS()\n");
  if (ps)
    *ps = opflow->ps;
  PetscFunctionReturn(0);
}

/* OPFLOWSetUpPS - Sets the underlying PS network object to be used by OPFLOW

  Input Parameters
. opflow - the OPFLOW object

  Notes: This function is an intermediate function that can be called for
setting up the PS network object prior to solving OPFLOW. A typical use-case is
some network parameter needs changing before solving opflow. In such case, the
work flow would be

  1. OPFLOWCreate();
  2. OPFLOWReadMatPowerData();
  3. OPFLOWSetUpPS();
  4. OPFLOWGetPS();
  ... change the network data by the PS object retrieved via OPFLOWGetPS().
  5. OPFLOWSolve();

 Skip steps 3 and 4 if no network changes are needed to be done.
*/
PetscErrorCode OPFLOWSetUpPS(OPFLOW opflow) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PSSetUp(opflow->ps);
  CHKERRQ(ierr);
  /* set individual load costs if necessary */

  PetscFunctionReturn(0);
}

/*
  OPFLOWSetAuxillaryObjective - Set additional objective and gradient functions

  Input Parameters:
+ opflow - opflow object
. objfunc - objective function
. gradfunc - associated function
- userctx - user struct to pass data needed for objective and gradient
calculation

  Notes:
    objfunc and gradfunc are called at every iteration of the optimization

    objfunc and should have the following signatures:
        PetscErrorCode MyObjectiveFunction(OPFLOW opflow,const double *x, double
*obj, void* userctx) PetscErrorCode MyGradientFunction(OPFLOW opflow,const
double *x, double *grad, void* userctx)
*/

PetscErrorCode OPFLOWSetAuxillaryObjective(OPFLOW opflow,
                                           OPFLOWAuxObjectiveFunction objfunc,
                                           OPFLOWAuxGradientFunction gradfunc,
                                           OPFLOWAuxHessianFunction hessfunc,
                                           void *userctx) {
  PetscFunctionBegin;
  opflow->modelops.computeauxobjective = objfunc;
  opflow->modelops.computeauxgradient = gradfunc;
  opflow->modelops.computeauxhessian = hessfunc;
  opflow->userctx = userctx;
  PetscFunctionReturn(0);
}

/* Sets a user-defined function to update/change the vaariable bounds.
   Called in ComputeVariableBounds after OPFLOW computes its variable bounds
*/
PetscErrorCode OPFLOWSetUpdateVariableBoundsFunction(
    OPFLOW opflow, PetscErrorCode (*updatefunc)(OPFLOW, Vec, Vec, void *),
    void *ctx) {

  PetscFunctionBegin;
  opflow->modelops.updatevariablebounds = updatefunc;
  opflow->userctx = ctx;
  PetscFunctionReturn(0);
}

/**
 * @brief Sets intialization type for opflow
 *
 * @param[in] opflow the OPFLOW object.
 * @param[in] type the initialization type.
 *
 * @note Must be called before OPFLOWSetUp
 */
PetscErrorCode OPFLOWSetInitializationType(OPFLOW opflow,
                                           OPFLOWInitializationType type) {
  PetscFunctionBegin;
  opflow->initializationtype = type;
  PetscFunctionReturn(0);
}

/**
 * @brief Gets intialization type for opflow
 *
 * @param[in] opflow the OPFLOW object.
 * @param[in] type the initialization type.
 */
PetscErrorCode OPFLOWGetInitializationType(OPFLOW opflow,
                                           OPFLOWInitializationType *type) {
  PetscFunctionBegin;
  *type = opflow->initializationtype;
  PetscFunctionReturn(0);
}

/**
 * @brief Sets ignore_lineflow_constraints for opflow
 *
 * @param[in] opflow the OPFLOW object.
 * @param[in] set value for opflow->ignore_lineflow_constraints.
 *
 * @note Must be called before OPFLOWSetUp
 */
PetscErrorCode OPFLOWIgnoreLineflowConstraints(OPFLOW opflow, PetscBool set) {
  PetscFunctionBegin;
  opflow->ignore_lineflow_constraints = set;
  PetscFunctionReturn(0);
}

/**
 * @brief Gets ignore_lineflow_constraints
 *
 * @param[in]  opflow the OPFLOW object.
 * @param[out] value for opflow->ignore_lineflow_constraints.
 */
PetscErrorCode OPFLOWGetIgnoreLineflowConstraints(OPFLOW opflow,
                                                  PetscBool *set) {
  PetscFunctionBegin;
  *set = opflow->ignore_lineflow_constraints;
  PetscFunctionReturn(0);
}

/**
 * @brief Allow line flow limit violations for opflow
 *
 * @param[in] opflow the OPFLOW object.
 * @param[in] set value for opflow->allow_lineflow_violation.
 *
 * @note Must be called before OPFLOWSetUp
 */
PetscErrorCode OPFLOWAllowLineflowViolation(OPFLOW opflow, PetscBool set) {
  PetscFunctionBegin;
  opflow->allow_lineflow_violation = set;
  PetscFunctionReturn(0);
}

/**
 * @brief Gets the value set for allow_lineflow_violation
 *
 * @param[in]  opflow the OPFLOW object.
 * @param[out] value for opflow->ignore_lineflow_constraints.
 */
PetscErrorCode OPFLOWGetAllowLineFlowViolation(OPFLOW opflow, PetscBool *set) {
  PetscFunctionBegin;
  *set = opflow->allow_lineflow_violation;
  PetscFunctionReturn(0);
}

/**
 * @brief Set penalty for line flow limit violation
 *
 * @param [in] opflow the OPFLOW object
 * @param [in] line flow limit violation penalty
 *
 *
 * @note Must be called before OPFLOWSetUp

*/
PetscErrorCode OPFLOWSetLineFlowViolationPenalty(OPFLOW opflow,
                                                 PetscReal penalty) {
  PetscFunctionBegin;
  opflow->lineflowviolation_penalty = penalty;
  PetscFunctionReturn(0);
}

/**
 * @brief Returns the penalty set for line flow limit violation
 *
 * @param [in] opflow the OPFLOW object
 * @param [out] penalty set for line flow limit violation
 */
PetscErrorCode OPFLOWGetLineFlowViolationPenalty(OPFLOW opflow,
                                                 PetscReal *penalty) {
  PetscFunctionBegin;
  *penalty = opflow->lineflowviolation_penalty;
  PetscFunctionReturn(0);
}

/*
  OPFLOWMonitorLines - List of lines to monitor. The flows for these lines
                       are included as inequality constraints in OPFLOW

 Input Parameter:
+ opflow      - OPFLOW object
. nkvlevels   - Number of kvlevels to monitor (Use -1 to monitor all kvlevels)
. kvlevels    - line kvlevels to monitor
- monitorfile - File with list of lines to monitor.

  Notes:
    The lines to monitor are either specified through a file OR by
    kvlevels, but not both. Use NULL for monitorfile if file is not set.
    If monitorfile is given then the kvlevels are ignored.

*/
PetscErrorCode OPFLOWSetLinesMonitored(OPFLOW opflow, PetscInt nkvlevels,
                                       const PetscScalar *kvlevels,
                                       const char *monitorfile) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  if (monitorfile != NULL) {
    SETERRQ(opflow->comm->type, PETSC_ERR_SUP,
            "Providing line list via file not yet supported");
  }

  if (nkvlevels < 0) {
    opflow->nlinekvmon = opflow->ps->nkvlevels;
    ierr = PetscMemcpy(opflow->linekvmon, opflow->ps->kvlevels,
                       opflow->nlinekvmon * sizeof(PetscScalar));
    ExaGOCheckError(ierr);
  } else if (nkvlevels == 0) {
    opflow->ignore_lineflow_constraints = PETSC_TRUE;
    opflow->nlinekvmon = 0;
    opflow->linesmon = NULL;
  } else {
    opflow->nlinekvmon = nkvlevels;
    ierr = PetscMemcpy(opflow->linekvmon, kvlevels,
                       nkvlevels * sizeof(PetscScalar));
    ExaGOCheckError(ierr);
  }

  PetscFunctionReturn(ierr);
}

/*
  OPFLOWSetSummaryStats - Sets the summary stats for the OPFLOW run

  Inputs:
. opflow - the OPFLOW object
*/
PetscErrorCode OPFLOWSetSummaryStats(OPFLOW opflow) {
  PetscErrorCode ierr;
  PS ps = opflow->ps;

  PetscFunctionBegin;

  ierr = PSComputeSummaryStats(ps);
  CHKERRQ(ierr);

  ps->solve_real_time = opflow->solve_real_time;
  ps->solve_cpu_time = opflow->solve_cpu_time;

  PetscFunctionReturn(0);
}

/**
 * @brief Check used internally to verify model + solver compatibility
 *
 * @param[in] opflow model to verify
 *
 * @note Models and solvers can only be defined per what exago is built with.
 */
PetscErrorCode OPFLOWCheckModelSolverCompatibility(OPFLOW opflow) {
  PetscFunctionBegin;
#if defined(EXAGO_ENABLE_IPOPT)
  PetscBool ipopt, ipopt_pbpol, ipopt_pbcar, ipopt_ibcar, ipopt_ibcar2,
      ipopt_dcopf;
  ipopt = static_cast<PetscBool>(opflow->solvername == OPFLOWSOLVER_IPOPT);
  ipopt_pbpol = static_cast<PetscBool>(opflow->modelname == OPFLOWMODEL_PBPOL);
  ipopt_pbcar = static_cast<PetscBool>(opflow->modelname == OPFLOWMODEL_PBCAR);
  ipopt_ibcar = static_cast<PetscBool>(opflow->modelname == OPFLOWMODEL_IBCAR);
  ipopt_ibcar2 =
      static_cast<PetscBool>(opflow->modelname == OPFLOWMODEL_IBCAR2);
  ipopt_dcopf = static_cast<PetscBool>(opflow->modelname == "DCOPF");
  if (ipopt && !(ipopt_pbpol || ipopt_pbcar || ipopt_ibcar || ipopt_ibcar2 ||
                 ipopt_dcopf)) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
            "OPFLOW solver IPOPT incompatible with model %s",
            opflow->modelname.c_str());
  }
#endif
#if defined(EXAGO_ENABLE_HIOP)
  PetscBool hiop, hiop_pbpol;
  hiop = static_cast<PetscBool>(opflow->solvername == OPFLOWSOLVER_HIOP);
  hiop_pbpol =
      static_cast<PetscBool>(opflow->modelname == OPFLOWMODEL_PBPOLHIOP);
#if defined(EXAGO_ENABLE_RAJA)
  PetscBool rajahiop_pbpol;
  rajahiop_pbpol =
      static_cast<PetscBool>(opflow->modelname == OPFLOWMODEL_PBPOLRAJAHIOP);
#else
  PetscBool rajahiop_pbpol = PETSC_FALSE;
#endif // RAJA
  if (hiop && !(hiop_pbpol || rajahiop_pbpol)) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
            "OPFLOW solver HIOP incompatible with model %s",
            (opflow->modelname).c_str());
  }
#if defined(EXAGO_ENABLE_HIOP_SPARSE)
  PetscBool hiop_sparse;
  PetscBool hiop_sparse_pbpol;
  PetscBool hiop_sparse_dcopf;
  hiop_sparse =
      static_cast<PetscBool>(opflow->solvername == OPFLOWSOLVER_HIOPSPARSE);
  hiop_sparse_pbpol =
      static_cast<PetscBool>(opflow->modelname == OPFLOWMODEL_PBPOL);
  hiop_sparse_dcopf = static_cast<PetscBool>(opflow->modelname == "DCOPF");
  if (hiop_sparse && !(hiop_sparse_pbpol || hiop_sparse_dcopf)) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
            "OPFLOW solver HIOPSPARSE incompatible with model %s",
            (opflow->modelname).c_str());
  }
#if defined(EXAGO_ENABLE_RAJA)
  PetscBool hiop_sparsegpu;
  PetscBool pbpolrajahiopsparse;

  hiop_sparsegpu =
      static_cast<PetscBool>(opflow->solvername == OPFLOWSOLVER_HIOPSPARSEGPU);
  pbpolrajahiopsparse = static_cast<PetscBool>(opflow->modelname ==
                                               OPFLOWMODEL_PBPOLRAJAHIOPSPARSE);

  if (hiop_sparsegpu && !pbpolrajahiopsparse) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
            "OPFLOW solver HIOPSPARSE incompatible with model %s",
            (opflow->modelname).c_str());
  }
#endif // EXAGO_ENABLE_RAJA
#endif // HIOP_SPARSE
  /* FIXME: The PBPOLRAJAHIOP has trouble with individual load
     shedding, so avoid using when load shedding is enabled */
  if (opflow->include_loadloss_variables) {
    if (rajahiop_pbpol) {
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
              "OPFLOW load shedding incompatible with model %s",
              (opflow->modelname).c_str());
    }
  }
#endif // HIOP
  PetscFunctionReturn(0);
}
