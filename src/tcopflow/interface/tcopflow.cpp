#include <private/opflowimpl.h>
#include <private/tcopflowimpl.h>
#include <utils.h>

/*
  TCOPFLOWSetTimeStepandDuration - Sets the time-step and optimization horizon

  Input Parameters:
+ tcopflow - the TCOPFLOW object
. dT       - time step in minutes
- duration - duration (horizon) in hours
*/
PetscErrorCode TCOPFLOWSetTimeStepandDuration(TCOPFLOW tcopflow, PetscReal dT,
                                              PetscReal duration) {
  PetscFunctionBegin;
  tcopflow->dT = dT;
  tcopflow->duration = duration;
  PetscFunctionReturn(0);
}

/*
  TCOPFLOWSetLoadProfiles - Sets the data files for time-varying load profiles

  Input Parameter
+  tcopflow - The TCOPFLOW object
.  ploadproffile - The name of the file with real power load variationt
-  qloadproffile - The name of the file with reactive power load variationt
*/
PetscErrorCode TCOPFLOWSetLoadProfiles(TCOPFLOW tcopflow,
                                       const char ploadprofile[],
                                       const char qloadprofile[]) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if (ploadprofile) {
    ierr = PetscMemcpy(tcopflow->ploadprofile, ploadprofile,
                       PETSC_MAX_PATH_LEN * sizeof(char));
    CHKERRQ(ierr);
    tcopflow->ploadprofileset = PETSC_TRUE;
  }

  if (qloadprofile) {
    ierr = PetscMemcpy(tcopflow->qloadprofile, qloadprofile,
                       PETSC_MAX_PATH_LEN * sizeof(char));
    CHKERRQ(ierr);
    tcopflow->qloadprofileset = PETSC_TRUE;
  }

  PetscFunctionReturn(0);
}

/*
  TCOPFLOWWindGenProfiles - Sets the data file for wind generation profiles

  Input Parameter
+  tcopflow - The TCOPFLOW object
-  windgenproffile - The name of the file with wind generation profile
*/
PetscErrorCode TCOPFLOWSetWindGenProfiles(TCOPFLOW tcopflow,
                                          const char windgenprofile[]) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if (windgenprofile) {
    ierr = PetscMemcpy(tcopflow->windgenprofile, windgenprofile,
                       PETSC_MAX_PATH_LEN * sizeof(char));
    CHKERRQ(ierr);
    tcopflow->windgenprofileset = PETSC_TRUE;
  }

  PetscFunctionReturn(0);
}

/*
  TCOPFLOWCreate - Creates a multi-period optimal power flow application object

  Input Parameters
. mpicomm - The MPI communicator

  Output Parameters
. tcopflowout - The multi-period optimal power flow application object
*/
PetscErrorCode TCOPFLOWCreate(MPI_Comm mpicomm, TCOPFLOW *tcopflowout) {
  PetscErrorCode ierr;
  TCOPFLOW tcopflow;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1, &tcopflow);
  CHKERRQ(ierr);

  ierr = COMMCreate(mpicomm, &tcopflow->comm);
  CHKERRQ(ierr);

  tcopflow->Nconeq = tcopflow->nconeq = 0;
  tcopflow->Nconineq = tcopflow->nconineq = 0;
  tcopflow->Ncon = tcopflow->ncon = 0;
  tcopflow->Nx = tcopflow->nx = 0;

  tcopflow->obj_factor = 1.0;
  tcopflow->obj = 0.0;

  tcopflow->duration = 1.0;
  tcopflow->dT = 60.0;
  tcopflow->Nt = 1;
  tcopflow->ploadprofileset = PETSC_FALSE;
  tcopflow->qloadprofileset = PETSC_FALSE;
  tcopflow->windgenprofileset = PETSC_FALSE;

  tcopflow->solver = NULL;
  tcopflow->model = NULL;

  tcopflow->ctgc = NULL;
  tcopflow->scen = NULL;

  tcopflow->nmodelsregistered = 0;
  tcopflow->TCOPFLOWModelRegisterAllCalled = PETSC_FALSE;

  /* Register all models */
  ierr = TCOPFLOWModelRegisterAll(tcopflow);

  tcopflow->nsolversregistered = 0;
  tcopflow->TCOPFLOWSolverRegisterAllCalled = PETSC_FALSE;

  /* Register all solvers */
  ierr = TCOPFLOWSolverRegisterAll(tcopflow);

  /* Run-time optiont */
  tcopflow->iscoupling = PETSC_TRUE;

  tcopflow->setupcalled = PETSC_FALSE;
  *tcopflowout = tcopflow;

  ExaGOLog(EXAGO_LOG_INFO, "{}", "TCOPFLOW: Application created");
  PetscFunctionReturn(0);
}

/*
  TCOPFLOWDestroy - Destroys the multi-period optimal power flow application
object

  Input Parameter
. tcopflow - The TCOPFLOW object to destroy
*/
PetscErrorCode TCOPFLOWDestroy(TCOPFLOW *tcopflow) {
  PetscErrorCode ierr;
  PetscInt i;

  PetscFunctionBegin;
  ierr = COMMDestroy(&(*tcopflow)->comm);
  CHKERRQ(ierr);

  /* Solution vector */
  ierr = VecDestroy(&(*tcopflow)->X);
  CHKERRQ(ierr);

  /* Lower and upper bounds on X */
  ierr = VecDestroy(&(*tcopflow)->Xl);
  CHKERRQ(ierr);
  ierr = VecDestroy(&(*tcopflow)->Xu);
  CHKERRQ(ierr);

  /* Destroy gradient vector */
  ierr = VecDestroy(&(*tcopflow)->gradobj);
  CHKERRQ(ierr);

  /* Constraints vector */
  ierr = VecDestroy(&(*tcopflow)->G);
  CHKERRQ(ierr);

  ierr = VecDestroy(&(*tcopflow)->Gl);
  CHKERRQ(ierr);
  ierr = VecDestroy(&(*tcopflow)->Gu);
  CHKERRQ(ierr);

  ierr = VecDestroy(&(*tcopflow)->Lambda);
  CHKERRQ(ierr);

  ierr = MatDestroy(&(*tcopflow)->Jac);
  CHKERRQ(ierr);
  ierr = MatDestroy(&(*tcopflow)->Hes);
  CHKERRQ(ierr);

  if ((*tcopflow)->solverops.destroy) {
    ierr = ((*tcopflow)->solverops.destroy)(*tcopflow);
  }

  if ((*tcopflow)->modelops.destroy) {
    ierr = ((*tcopflow)->modelops.destroy)(*tcopflow);
  }

  /* Destroy OPFLOW objects */
  for (i = 0; i < (*tcopflow)->Nt; i++) {
    ierr = OPFLOWDestroy(&(*tcopflow)->opflows[i]);
    CHKERRQ(ierr);
  }

  ierr = PetscFree((*tcopflow)->xstarti);
  CHKERRQ(ierr);
  ierr = PetscFree((*tcopflow)->gstarti);
  CHKERRQ(ierr);
  ierr = PetscFree((*tcopflow)->nxi);
  CHKERRQ(ierr);
  ierr = PetscFree((*tcopflow)->ngi);
  CHKERRQ(ierr);
  ierr = PetscFree((*tcopflow)->nconineqcoup);
  CHKERRQ(ierr);
  ierr = PetscFree((*tcopflow)->opflows);
  CHKERRQ(ierr);
  ierr = PetscFree(*tcopflow);
  CHKERRQ(ierr);
  //  *tcopflow = 0;
  PetscFunctionReturn(0);
}

/*
  TCOPFLOWSetModel - Sets the model for TCOPFLOW

  Input Parameters:
+ tcopflow - opflow application object
- modelname - name of the model
*/
PetscErrorCode TCOPFLOWSetModel(TCOPFLOW tcopflow, const char *modelname) {
  PetscErrorCode ierr, (*r)(TCOPFLOW) = NULL;
  PetscInt i;
  PetscFunctionBegin;
  PetscBool match;
  for (i = 0; i < tcopflow->nmodelsregistered; i++) {
    ierr = PetscStrcmp(tcopflow->TCOPFLOWModelList[i].name, modelname, &match);
    CHKERRQ(ierr);
    if (match) {
      r = tcopflow->TCOPFLOWModelList[i].create;
      break;
    }
  }

  if (!r)
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE,
             "Unknown type for TCOPFLOW Model %s", modelname);

  /* Null the function pointers */
  tcopflow->modelops.destroy = 0;
  tcopflow->modelops.setup = 0;
  tcopflow->modelops.setnumvariablesandconstraints = 0;
  tcopflow->modelops.setvariablebounds = 0;
  tcopflow->modelops.setconstraintbounds = 0;
  tcopflow->modelops.setvariableandconstraintbounds = 0;
  tcopflow->modelops.setinitialguess = 0;
  tcopflow->modelops.computeconstraints = 0;
  tcopflow->modelops.computejacobian = 0;
  tcopflow->modelops.computehessian = 0;
  tcopflow->modelops.computeobjandgradient = 0;
  tcopflow->modelops.computeobjective = 0;
  tcopflow->modelops.computegradient = 0;

  ierr = PetscStrcpy(tcopflow->modelname, modelname);
  CHKERRQ(ierr);
  /* Call the underlying implementation constructor */
  ierr = (*r)(tcopflow);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  TCOPFLOWSetSolver - Sets the solver for TCOPFLOW

  Input Parameters:
+ tcopflow - tcopflow application object
- solvername - name of the solver
*/
PetscErrorCode TCOPFLOWSetSolver(TCOPFLOW tcopflow, const char *solvername) {
  PetscErrorCode ierr, (*r)(TCOPFLOW) = NULL;
  PetscInt i;
  PetscFunctionBegin;
  PetscBool match;
  for (i = 0; i < tcopflow->nsolversregistered; i++) {
    ierr =
        PetscStrcmp(tcopflow->TCOPFLOWSolverList[i].name, solvername, &match);
    CHKERRQ(ierr);
    if (match) {
      r = tcopflow->TCOPFLOWSolverList[i].create;
      break;
    }
  }

  if (!r)
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE,
             "Unknown type for TCOPFLOW Solver %s", solvername);

  /* Initialize (Null) the function pointers */
  tcopflow->solverops.destroy = 0;
  tcopflow->solverops.solve = 0;
  tcopflow->solverops.setup = 0;

  ierr = PetscStrcpy(tcopflow->solvername, solvername);
  CHKERRQ(ierr);
  /* Call the underlying implementation constructor */
  ierr = (*r)(tcopflow);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  TCOPFLOWSetNetworkData - Sets and reads the network data

  Input Parameter
+  tcopflow - The TCOPFLOW object
-  netfile - The name of the network file

  Notes: The input data must be in MATPOWER data format
*/
PetscErrorCode TCOPFLOWSetNetworkData(TCOPFLOW tcopflow, const char netfile[]) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscMemcpy(tcopflow->netfile, netfile,
                     PETSC_MAX_PATH_LEN * sizeof(char));
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode TCOPFLOWSetContingency(TCOPFLOW tcopflow, Contingency *ctgc) {
  PetscFunctionBegin;
  tcopflow->ctgc = ctgc;
  PetscFunctionReturn(0);
}

PetscErrorCode TCOPFLOWSetScenario(TCOPFLOW tcopflow, Scenario *scen) {
  PetscFunctionBegin;
  tcopflow->scen = scen;
  PetscFunctionReturn(0);
}

/*
  TCOPFLOWSetUp - Sets up an multi-period optimal power flow application object

  Input Parameters:
. tcopflow - the TCOPFLOW object

  Notes:
  This routine sets up the TCOPFLOW object and the underlying PS object. It
  also distributes the PS object when used in parallel.
*/
PetscErrorCode TCOPFLOWSetUp(TCOPFLOW tcopflow) {
  PetscErrorCode ierr;
  PetscBool solverset;
  char modelname[32] = "GENRAMP";
  char solvername[32] = "IPOPT";
  PS ps;
  PetscInt i;
  OPFLOW opflow;

  PetscFunctionBegin;

  ierr =
      PetscOptionsBegin(tcopflow->comm->type, NULL, "TCOPFLOW options", NULL);
  CHKERRQ(ierr);
  ierr = PetscOptionsString("-tcopflow_model", "TCOPFLOW model type", "",
                            modelname, modelname, 32, NULL);
  CHKERRQ(ierr);
  ierr = PetscOptionsString("-tcopflow_solver", "TCOPFLOW solver type", "",
                            solvername, solvername, 32, &solverset);
  CHKERRQ(ierr);
  ierr =
      PetscOptionsBool("-tcopflow_iscoupling",
                       "Include coupling between first stage and second stage",
                       "", tcopflow->iscoupling, &tcopflow->iscoupling, NULL);
  CHKERRQ(ierr);
  ierr = PetscOptionsReal("-tcopflow_dT", "Length of time-step (minutes)", "",
                          tcopflow->dT, &tcopflow->dT, NULL);
  CHKERRQ(ierr);
  ierr = PetscOptionsReal("-tcopflow_duration", "Time horizon (hours)", "",
                          tcopflow->duration, &tcopflow->duration, NULL);
  CHKERRQ(ierr);
  ierr = PetscOptionsReal("-tcopflow_tolerance", "optimization tolerance", "",
                          tcopflow->tolerance, &tcopflow->tolerance, NULL);
  CHKERRQ(ierr);
  PetscOptionsEnd();

  tcopflow->Nt = round(tcopflow->duration * 60.0 / tcopflow->dT) + 1;

  ExaGOLog(EXAGO_LOG_INFO,
           "TCOPFLOW: Duration = {:f} hours, timestep = {:f} minutes, number "
           "of time-steps = {:d}",
           tcopflow->duration, tcopflow->dT, tcopflow->Nt);

  /* Set Model */
  ierr = TCOPFLOWSetModel(tcopflow, modelname);
  CHKERRQ(ierr);

  /* Set solver */
  if (solverset) {
    if (tcopflow->solver)
      ierr = (*tcopflow->solverops.destroy)(tcopflow);
    ierr = TCOPFLOWSetSolver(tcopflow, solvername);
    CHKERRQ(ierr);
    ExaGOLog(EXAGO_LOG_INFO, "TCOPFLOW: Using {} solver", solvername);
  } else {
    if (!tcopflow->solver) {
      ierr = TCOPFLOWSetSolver(tcopflow, TCOPFLOWSOLVER_IPOPT);
      CHKERRQ(ierr);
      ExaGOLog(EXAGO_LOG_INFO, "TCOPFLOW: Using {} solver",
               TCOPFLOWSOLVER_IPOPT);
    }
  }

  /* Create OPFLOW objects */
  ierr = PetscCalloc1(tcopflow->Nt, &tcopflow->opflows);
  CHKERRQ(ierr);
  for (i = 0; i < tcopflow->Nt; i++) {
    ierr = OPFLOWCreate(PETSC_COMM_SELF, &tcopflow->opflows[i]);
    CHKERRQ(ierr);
    ierr = OPFLOWSetModel(tcopflow->opflows[i], OPFLOWMODEL_PBPOL);
    CHKERRQ(ierr);

    /* Read network data file */
    ierr = OPFLOWReadMatPowerData(tcopflow->opflows[i], tcopflow->netfile);
    CHKERRQ(ierr);
    /* Set up the PS object for opflow */
    ps = tcopflow->opflows[i]->ps;
    ierr = PSSetUp(ps);
    CHKERRQ(ierr);

    if (tcopflow->scen) {
      // Scenario is set, apply it */
      ierr = PSApplyScenario(ps, *tcopflow->scen);
      CHKERRQ(ierr);
    }

    if (tcopflow->ctgc) {
      // Contingency is set, apply it */
      ierr = PSApplyContingency(ps, *tcopflow->ctgc);
      CHKERRQ(ierr);
    }
    /* Set up OPFLOW object */
    ierr = OPFLOWSetUp(tcopflow->opflows[i]);
    CHKERRQ(ierr);
  }

  /* Read Pload profiles */
  if (tcopflow->ploadprofileset) {
    ierr = TCOPFLOWReadPloadProfile(tcopflow, tcopflow->ploadprofile);
  }

  /* Read Qload profiles */
  if (tcopflow->qloadprofileset) {
    ierr = TCOPFLOWReadQloadProfile(tcopflow, tcopflow->qloadprofile);
  }

  /* Read wind generation profiles */
  if (tcopflow->windgenprofileset) {
    ierr = TCOPFLOWReadWindGenProfile(tcopflow, tcopflow->windgenprofile);
  }

  ierr = PetscCalloc1(tcopflow->Nt, &tcopflow->nconineqcoup);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(tcopflow->Nt, &tcopflow->nxi);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(tcopflow->Nt, &tcopflow->ngi);
  CHKERRQ(ierr);

  /* Set number of variables and constraints */
  ierr = (*tcopflow->modelops.setnumvariablesandconstraints)(
      tcopflow, tcopflow->nxi, tcopflow->ngi, tcopflow->nconineqcoup);

  ierr = PetscCalloc1(tcopflow->Nt, &tcopflow->xstarti);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(tcopflow->Nt, &tcopflow->gstarti);
  CHKERRQ(ierr);

  tcopflow->nx = tcopflow->Nx = tcopflow->nxi[0];
  tcopflow->Ncon = tcopflow->ngi[0];
  tcopflow->Nconcoup = tcopflow->nconineqcoup[0];
  opflow = tcopflow->opflows[0];
  tcopflow->Nconeq = opflow->nconeq;
  tcopflow->Nconineq = opflow->nconineq;

  for (i = 1; i < tcopflow->Nt; i++) {
    tcopflow->xstarti[i] = tcopflow->xstarti[i - 1] + tcopflow->nxi[i - 1];
    tcopflow->gstarti[i] = tcopflow->gstarti[i - 1] + tcopflow->ngi[i - 1];
    tcopflow->Nx += tcopflow->nxi[i];
    tcopflow->Ncon += tcopflow->ngi[i];
    tcopflow->Nconcoup += tcopflow->nconineqcoup[i];
    opflow = tcopflow->opflows[i];
    tcopflow->Nconeq += opflow->nconeq;
    tcopflow->Nconineq += opflow->nconineq;
  }

  /* Create vector X */
  ierr = VecCreate(tcopflow->comm->type, &tcopflow->X);
  CHKERRQ(ierr);
  ierr = VecSetSizes(tcopflow->X, tcopflow->Nx, PETSC_DECIDE);
  CHKERRQ(ierr);
  ierr = VecSetFromOptions(tcopflow->X);
  CHKERRQ(ierr);

  ierr = VecDuplicate(tcopflow->X, &tcopflow->Xl);
  CHKERRQ(ierr);
  ierr = VecDuplicate(tcopflow->X, &tcopflow->Xu);
  CHKERRQ(ierr);
  ierr = VecDuplicate(tcopflow->X, &tcopflow->gradobj);
  CHKERRQ(ierr);

  /* vector for constraints */
  ierr = VecCreate(tcopflow->comm->type, &tcopflow->G);
  CHKERRQ(ierr);
  ierr = VecSetSizes(tcopflow->G, tcopflow->Ncon, PETSC_DECIDE);
  CHKERRQ(ierr);
  ierr = VecSetFromOptions(tcopflow->G);
  CHKERRQ(ierr);

  /* Constraint bounds vectors  */
  ierr = VecDuplicate(tcopflow->G, &tcopflow->Gl);
  CHKERRQ(ierr);
  ierr = VecDuplicate(tcopflow->G, &tcopflow->Gu);
  CHKERRQ(ierr);

  /* Constraint Jacobian */
  ierr = MatCreate(tcopflow->comm->type, &tcopflow->Jac);
  CHKERRQ(ierr);
  ierr = MatSetSizes(tcopflow->Jac, tcopflow->Ncon, tcopflow->Nx,
                     tcopflow->Ncon, tcopflow->Nx);
  CHKERRQ(ierr);
  ierr = MatSetUp(tcopflow->Jac);
  CHKERRQ(ierr);
  ierr = MatSetFromOptions(tcopflow->Jac);
  CHKERRQ(ierr);

  /* Hessian */
  ierr = MatCreate(tcopflow->comm->type, &tcopflow->Hes);
  CHKERRQ(ierr);
  ierr = MatSetSizes(tcopflow->Hes, tcopflow->Nx, tcopflow->Nx, tcopflow->Nx,
                     tcopflow->Nx);
  CHKERRQ(ierr);
  ierr = MatSetUp(tcopflow->Hes);
  CHKERRQ(ierr);
  ierr = MatSetFromOptions(tcopflow->Hes);
  CHKERRQ(ierr);

  /* Lagrangian multipliers */
  ierr = VecDuplicate(tcopflow->G, &tcopflow->Lambda);
  CHKERRQ(ierr);

  ierr = (*tcopflow->solverops.setup)(tcopflow);
  CHKERRQ(ierr);

  ExaGOLog(EXAGO_LOG_INFO, "{}", "TCOPFLOW: Setup completed");

  tcopflow->setupcalled = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/*
  TCOPFLOWSolve - Solves the AC multi-period optimal power flow

  Input Parameters:
. tcopflow - the multi-period optimal power flow application object
*/
PetscErrorCode TCOPFLOWSolve(TCOPFLOW tcopflow) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if (!tcopflow->setupcalled) {
    ierr = TCOPFLOWSetUp(tcopflow);
  }

  /* Set bounds on variables */
  if (tcopflow->modelops.setvariablebounds) {
    ierr = (*tcopflow->modelops.setvariablebounds)(tcopflow, tcopflow->Xl,
                                                   tcopflow->Xu);
    CHKERRQ(ierr);
  }

  /* Set bounds on constraints */
  if (tcopflow->modelops.setconstraintbounds) {
    ierr = (*tcopflow->modelops.setconstraintbounds)(tcopflow, tcopflow->Gl,
                                                     tcopflow->Gu);
    CHKERRQ(ierr);
  }

  /* Set initial guess */
  if (tcopflow->modelops.setinitialguess) {
    ierr = (*tcopflow->modelops.setinitialguess)(tcopflow, tcopflow->X);
    CHKERRQ(ierr);
  }

  ierr = VecSet(tcopflow->Lambda, 1.0);
  CHKERRQ(ierr);

  /* Solve */
  ierr = (*tcopflow->solverops.solve)(tcopflow);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  TCOPFLOWGetObjective - Returns the objective function value

  Input Parameters:
+ TCOPFLOW - the TCOPFLOW object
- obj    - the objective function value

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode TCOPFLOWGetObjective(TCOPFLOW tcopflow, PetscReal *obj) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*tcopflow->solverops.getobjective)(tcopflow, obj);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  TCOPFLOWGetSolution - Returns the TCOPFLOW solution for a given time-step

  Input Parameters:
+ TCOPFLOW - the TCOPFLOW object
. tnum     - time-step (0 is the first time-step)
- X        - the tcopflow solution for the given time-step

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode TCOPFLOWGetSolution(TCOPFLOW tcopflow, PetscInt tnum, Vec *X) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*tcopflow->solverops.getsolution)(tcopflow, tnum, X);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  TCOPFLOWGetConstraints - Returns the TCOPFLOW constraints for a given
time-step

  Input Parameters:
+ TCOPFLOW - the TCOPFLOW object
. tnum     - time-step (0 for the first time-step)
- G        - the tcopflow constraints

  Notes: Should be called after the optimization finishes.
         Equality constraints first followed by inequality constraints
*/
PetscErrorCode TCOPFLOWGetConstraints(TCOPFLOW tcopflow, PetscInt tnum,
                                      Vec *G) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*tcopflow->solverops.getconstraints)(tcopflow, tnum, G);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  TCOPFLOWGetConstraintMultipliers - Returns the TCOPFLOW constraint multipliers
for a given time-step

  Input Parameters:
+ TCOPFLOW - the TCOPFLOW object
. tnum     - time-step (0 for the first time-step)
- G    - the tcopflow constraint lagrange multipliers

  Notes: Should be called after the optimization finishes.
    Equality constraint multipliers first followed by inequality constraint
multipliers
*/
PetscErrorCode TCOPFLOWGetConstraintMultipliers(TCOPFLOW tcopflow,
                                                PetscInt tnum, Vec *Lambda) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr =
      (*tcopflow->solverops.getconstraintmultipliers)(tcopflow, tnum, Lambda);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  TCOPFLOWGetConvergenceStatus - Did TCOPFLOW converge?

  Input Parameters:
+ TCOPFLOW - the TCOPFLOW object
- status - PETSC_TRUE if converged, PETSC_FALSE otherwise

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode TCOPFLOWGetConvergenceStatus(TCOPFLOW tcopflow,
                                            PetscBool *status) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*tcopflow->solverops.getconvergencestatus)(tcopflow, status);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/**
 * @brief Sets the solver tolerance
 *
 * @param[out] tcopflow TCOPFLOW application object
 * @param[in]  tol      Solver tolerance
 *
 * Tolerance is specified during solve stage for each individual solver.
 *
 */
PetscErrorCode TCOPFLOWSetTolerance(TCOPFLOW tcopflow, PetscReal tol) {
  PetscFunctionBegin;
  tcopflow->tolerance = tol;
  PetscFunctionReturn(0);
}

/**
 * @brief Get the solver tolerance
 *
 * @param[out] tcopflow TCOPFLOW application object
 * @param[in]  tol      Solver tolerance
 */
PetscErrorCode TCOPFLOWGetTolerance(TCOPFLOW tcopflow, PetscReal *tol) {
  PetscFunctionBegin;
  *tol = tcopflow->tolerance;
  PetscFunctionReturn(0);
}

/**
 * @brief Returns the number of iterations for a given solver
 *
 * @param[in]  tcopflow TCOPFLOW application object
 * @param[out] iter     Number of iterations
 */
PetscErrorCode TCOPFLOWGetNumIterations(TCOPFLOW tcopflow, PetscInt *iter) {
  PetscErrorCode ierr;
  *iter = tcopflow->numiter;
  PetscFunctionReturn(0);
}
