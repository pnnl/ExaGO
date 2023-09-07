#include <private/opflowimpl.h>

/*
  OPFLOWPrintSolution - Prints OPFLOW details to stdout.

  Input Parameters:
+ opflow - the OPFLOW object

*/
PetscErrorCode OPFLOWPrintSolution(OPFLOW opflow) {
  PetscErrorCode ierr;
  PetscBool conv_status;
  PetscReal cost;

  PetscFunctionBegin;
  if (!opflow->solutiontops) {
    ierr = OPFLOWSolutionToPS(opflow);
    CHKERRQ(ierr);
  }

#ifndef EXAGO_DISABLE_LOGGING
  /* Print to stdout */
  ierr = PetscPrintf(
      opflow->comm->type,
      "=============================================================\n");
  CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type, "\t\tOptimal Power Flow\n");
  CHKERRQ(ierr);
  ierr = PetscPrintf(
      opflow->comm->type,
      "=============================================================\n");
  CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type, "%-35s %s\n", "Model",
                     opflow->modelname.c_str());
  CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type, "%-35s %s\n", "Solver",
                     opflow->solvername.c_str());
  CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type, "%-35s %s\n", "Objective",
                     OPFLOWObjectiveTypes[opflow->objectivetype]);
  CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type, "%-35s %s\n", "Initialization",
                     OPFLOWInitializationTypes[opflow->initializationtype]);
  CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type, "%-35s %s\n", "Gen. bus voltage mode",
                     OPFLOWGenBusVoltageTypes[opflow->genbusvoltagetype]);
  CHKERRQ(ierr);

  ierr = PetscPrintf(opflow->comm->type, "%-35s %s\n", "Load loss allowed",
                     opflow->include_loadloss_variables ? "YES" : "NO");
  CHKERRQ(ierr);
  if (opflow->include_loadloss_variables) {
    ierr = PetscPrintf(opflow->comm->type, "%-35s %g\n",
                       "Load loss penalty ($)", opflow->loadloss_penalty);
    CHKERRQ(ierr);
  }

  ierr =
      PetscPrintf(opflow->comm->type, "%-35s %s\n", "Power imbalance allowed",
                  opflow->include_powerimbalance_variables ? "YES" : "NO");
  CHKERRQ(ierr);
  if (opflow->include_powerimbalance_variables) {
    ierr = PetscPrintf(opflow->comm->type, "%-35s %g\n",
                       "Power imbalance penalty ($)",
                       opflow->powerimbalance_penalty);
    CHKERRQ(ierr);
  }

  ierr = PetscPrintf(opflow->comm->type, "%-35s %s\n",
                     "Ignore line flow constraints",
                     opflow->ignore_lineflow_constraints ? "YES" : "NO");
  CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type, "\n");
  CHKERRQ(ierr);
  MPI_Barrier(opflow->comm->type);

  ierr = PetscPrintf(opflow->comm->type, "%-35s %d\n", "Number of variables",
                     opflow->Nx);
  CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type, "%-35s %d\n",
                     "Number of equality constraints", opflow->Nconeq);
  CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type, "%-35s %d\n",
                     "Number of inequality constraints", opflow->Nconineq);
  CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type, "\n");
  CHKERRQ(ierr);

  ierr = OPFLOWGetConvergenceStatus(opflow, &conv_status);
  CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type, "%-35s %s\n", "Convergence status",
                     conv_status ? "CONVERGED" : "DID NOT CONVERGE");
  CHKERRQ(ierr);
  ierr = OPFLOWGetObjective(opflow, &cost);
  CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type, "%-35s %-7.2f\n", "Objective value",
                     cost);
  CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type, "\n");
  CHKERRQ(ierr);

  MPI_Barrier(opflow->comm->type);

  ierr = PSPrintSystemSummary(opflow->ps);
  CHKERRQ(ierr);
#endif

  PetscFunctionReturn(0);
}

/*
  OPFLOWSaveSolution - Saves the OPFLOW solution to file

  Input Parameters:
+ opflow - the OPFLOW object
. format - the output file format (csv, matpower)
- outfile  - Name of output file
*/
PetscErrorCode OPFLOWSaveSolution(OPFLOW opflow, OutputFormat format,
                                  const char outfile[]) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscLogEventBegin(opflow->outputlogger, 0, 0, 0, 0);
  CHKERRQ(ierr);

  if (!opflow->solutiontops) {
    ierr = OPFLOWSolutionToPS(opflow);
    CHKERRQ(ierr);
  }

  ierr = PSSaveSolution(opflow->ps, format, outfile);
  CHKERRQ(ierr);

  ierr = PetscLogEventEnd(opflow->outputlogger, 0, 0, 0, 0);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  OPFLOWSaveSolution - Saves the OPFLOW solution to file using
                       the output format set (default MATPOWER)

  Input Parameters:
+ opflow - the OPFLOW object
- outfile  - Name of output file
*/
PetscErrorCode OPFLOWSaveSolutionDefault(OPFLOW opflow, const char outfile[]) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = OPFLOWSaveSolution(opflow, opflow->outputformat, outfile);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
