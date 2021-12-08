#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>
#include <private/tcopflowimpl.h>

/*
  SCOPFLOWPrintSolution - Prints SCOPFLOW solution to stdout for the given
contingency

  Input Parameters:
+ scopflow - the SCOPFLOW object
- cont_num - the contingency number (0 for base)
  Notes:
   Only prints out the system summary for base case

*/
PetscErrorCode SCOPFLOWPrintSolution(SCOPFLOW scopflow, PetscInt cont_num) {
  PetscErrorCode ierr;
  PetscBool conv_status;
  PetscReal cost;
  PetscInt c;
  TCOPFLOW tcopflow;
  OPFLOW opflow;
  const PetscScalar *x, *lambda;

  PetscFunctionBegin;
  MPI_Barrier(scopflow->comm->type);
  c = cont_num - scopflow->cstart;
  if (cont_num >= scopflow->cstart && cont_num < scopflow->cend) {
    if (!scopflow->ismultiperiod) {
      opflow = scopflow->opflows[c];

      if (!opflow->solutiontops) {
        ierr = SCOPFLOWGetSolution(scopflow, cont_num, &opflow->X);
        CHKERRQ(ierr);
        ierr = SCOPFLOWGetConstraintMultipliers(scopflow, cont_num,
                                                &opflow->Lambda);
        CHKERRQ(ierr);
        ierr = OPFLOWSolutionToPS(opflow);
        CHKERRQ(ierr);
      }
    } else {
      tcopflow = scopflow->tcopflows[c];

      ierr = SCOPFLOWGetSolution(scopflow, cont_num, &tcopflow->X);
      CHKERRQ(ierr);
      ierr = SCOPFLOWGetConstraintMultipliers(scopflow, cont_num,
                                              &tcopflow->Lambda);
      CHKERRQ(ierr);

      ierr = VecGetArrayRead(tcopflow->X, &x);
      CHKERRQ(ierr);
      ierr = VecGetArrayRead(tcopflow->Lambda, &lambda);
      CHKERRQ(ierr);

      opflow = tcopflow->opflows[0];

      ierr = VecPlaceArray(opflow->X, x + tcopflow->xstarti[0]);
      ierr = VecPlaceArray(opflow->Lambda, x + tcopflow->gstarti[0]);

      if (!opflow->solutiontops) {
        ierr = OPFLOWSolutionToPS(opflow);
        CHKERRQ(ierr);
      }

      ierr = VecResetArray(opflow->X);
      ierr = VecResetArray(opflow->Lambda);

      ierr = VecRestoreArrayRead(tcopflow->X, &x);
      CHKERRQ(ierr);
      ierr = VecRestoreArrayRead(tcopflow->Lambda, &lambda);
      CHKERRQ(ierr);
    }

    /* Print to stdout */
    ierr = PetscPrintf(
        scopflow->comm->type,
        "=============================================================\n");
    CHKERRQ(ierr);
    ierr = PetscPrintf(scopflow->comm->type,
                       "\tSecurity-Constrained Optimal Power Flow\n");
    CHKERRQ(ierr);
    ierr = PetscPrintf(
        scopflow->comm->type,
        "=============================================================\n");
    CHKERRQ(ierr);

    ierr = PetscPrintf(scopflow->comm->type, "%-35s %s\n", "OPFLOW Formulation",
                       opflow->modelname);
    CHKERRQ(ierr);
    ierr = PetscPrintf(scopflow->comm->type, "%-35s %s\n", "Solver",
                       scopflow->solvername);
    CHKERRQ(ierr);
    ierr = PetscPrintf(scopflow->comm->type, "%-35s %s\n", "Initialization",
                       OPFLOWInitializationTypes[opflow->initializationtype]);
    CHKERRQ(ierr);
    ierr = PetscPrintf(scopflow->comm->type, "%-35s %d\n",
                       "Number of contingencies", scopflow->Nc - 1);
    CHKERRQ(ierr);

    ierr = PetscPrintf(scopflow->comm->type, "%-35s %s\n", "Load loss allowed",
                       opflow->include_loadloss_variables ? "YES" : "NO");
    CHKERRQ(ierr);
    if (opflow->include_loadloss_variables) {
      ierr = PetscPrintf(scopflow->comm->type, "%-35s %g\n",
                         "Load loss penalty ($)", opflow->loadloss_penalty);
      CHKERRQ(ierr);
    }

    ierr = PetscPrintf(scopflow->comm->type, "%-35s %s\n",
                       "Power imbalance allowed",
                       opflow->include_powerimbalance_variables ? "YES" : "NO");
    CHKERRQ(ierr);
    if (opflow->include_powerimbalance_variables) {
      ierr = PetscPrintf(scopflow->comm->type, "%-35s %g\n",
                         "Power imbalance penalty ($)",
                         opflow->powerimbalance_penalty);
      CHKERRQ(ierr);
    }

    ierr = PetscPrintf(scopflow->comm->type, "%-35s %s\n",
                       "Ignore line flow constraints",
                       opflow->ignore_lineflow_constraints ? "YES" : "NO");
    CHKERRQ(ierr);
    ierr = PetscPrintf(scopflow->comm->type, "\n");
    CHKERRQ(ierr);

    ierr = PetscPrintf(scopflow->comm->type, "%-35s %d\n",
                       "Number of variables", scopflow->Nx);
    CHKERRQ(ierr);
    ierr = PetscPrintf(scopflow->comm->type, "%-35s %d\n",
                       "Number of equality constraints", scopflow->Nconeq);
    CHKERRQ(ierr);
    ierr = PetscPrintf(scopflow->comm->type, "%-35s %d\n",
                       "Number of inequality constraints", scopflow->Nconineq);
    CHKERRQ(ierr);
    ierr = PetscPrintf(scopflow->comm->type, "%-35s %d\n",
                       "Number of coupling constraints", scopflow->Nconcoup);
    CHKERRQ(ierr);

    ierr = PetscPrintf(scopflow->comm->type, "\n");
    CHKERRQ(ierr);
  }
  MPI_Barrier(scopflow->comm->type);
  ierr = SCOPFLOWGetConvergenceStatus(scopflow, &conv_status);
  CHKERRQ(ierr);
  ierr = PetscPrintf(scopflow->comm->type, "%-35s %s\n", "Convergence status",
                     conv_status ? "CONVERGED" : "DID NOT CONVERGE");
  CHKERRQ(ierr);
  ierr = SCOPFLOWGetObjective(scopflow, &cost);
  CHKERRQ(ierr);
  ierr = PetscPrintf(scopflow->comm->type, "%-35s %-7.2f\n", "Objective value",
                     cost);
  CHKERRQ(ierr);
  ierr = PetscPrintf(scopflow->comm->type, "\n");
  CHKERRQ(ierr);

  if (cont_num >= scopflow->cstart && cont_num < scopflow->cend) {
    ierr = PSPrintSystemSummary(opflow->ps);
    CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/*
  SCOPFLOWSaveSolution - Saves the SCOPFLOW solution for the given contingency
to file

  Input Parameters:
+ scopflow - the OPFLOW object
. cont_num - the contingency number (0 for base case)
. format - the output file format (csv, matpower)
- outfile  - Name of output file
*/
PetscErrorCode SCOPFLOWSaveSolution(SCOPFLOW scopflow, PetscInt cont_num,
                                    OutputFormat format, const char outfile[]) {
  PetscErrorCode ierr;
  TCOPFLOW tcopflow;
  OPFLOW opflow;

  PetscFunctionBegin;
  if (!scopflow->ismultiperiod) {
    opflow = scopflow->opflows[cont_num];
  } else {
    tcopflow = scopflow->tcopflows[cont_num];
    opflow = tcopflow->opflows[0];
  }

  if (!opflow->solutiontops) {
    ierr = SCOPFLOWGetSolution(scopflow, cont_num, &opflow->X);
    CHKERRQ(ierr);
    ierr =
        SCOPFLOWGetConstraintMultipliers(scopflow, cont_num, &opflow->Lambda);
    CHKERRQ(ierr);
    ierr = OPFLOWSolutionToPS(opflow);
    CHKERRQ(ierr);
  }

  ierr = OPFLOWSaveSolution(opflow, format, outfile);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  SCOPFLOWSaveSolutionAll - Saves all SCOPFLOW solutions

  Input Parameters:
+ scopflow - the OPFLOW object
. format - the output file format (csv, matpower)
- outdir  - Name of output directory

  Notes:
   Save all SCOPFLOW solutions (one contingency per file)
*/
PetscErrorCode SCOPFLOWSaveSolutionAll(SCOPFLOW scopflow, OutputFormat format,
                                       const char outdir[]) {
  PetscErrorCode ierr;
  TCOPFLOW tcopflow;
  OPFLOW opflow;
  PetscInt c, t;
  char filename[64];
  char outfile[256];
  char tcopflowdirname[64];
  char tcopflowdir[256];
  const PetscScalar *x, *lambda;

  PetscFunctionBegin;

  if (!scopflow->comm->rank) {
    ierr = PetscMkdir(outdir);
    CHKERRQ(ierr);
  }
  MPI_Barrier(scopflow->comm->type);

  for (c = 0; c < scopflow->nc; c++) {
    if (scopflow->ismultiperiod) {
      tcopflow = scopflow->tcopflows[c];

      ierr =
          PetscSNPrintf(tcopflowdirname, 64, "cont_%d", scopflow->cstart + c);
      CHKERRQ(ierr);
      ierr = PetscStrncpy(tcopflowdir, outdir, 256);
      CHKERRQ(ierr);
      ierr = PetscStrlcat(tcopflowdir, "/", 256);
      CHKERRQ(ierr);
      ierr = PetscStrlcat(tcopflowdir, tcopflowdirname, 256);
      CHKERRQ(ierr);

      ierr = SCOPFLOWGetSolution(scopflow, scopflow->cstart + c, &tcopflow->X);
      CHKERRQ(ierr);
      ierr = SCOPFLOWGetConstraintMultipliers(scopflow, scopflow->cstart + c,
                                              &tcopflow->Lambda);
      CHKERRQ(ierr);

      ierr = VecGetArrayRead(tcopflow->X, &x);
      CHKERRQ(ierr);
      ierr = VecGetArrayRead(tcopflow->Lambda, &lambda);
      CHKERRQ(ierr);

      for (t = 0; t < tcopflow->Nt; t++) {
        opflow = tcopflow->opflows[t];

        ierr = PetscMkdir(tcopflowdir);
        CHKERRQ(ierr);

        ierr = PetscSNPrintf(filename, 64, "t%d.m", t);
        CHKERRQ(ierr);
        ierr = PetscStrncpy(outfile, tcopflowdir, 256);
        CHKERRQ(ierr);
        ierr = PetscStrlcat(outfile, "/", 256);
        CHKERRQ(ierr);
        ierr = PetscStrlcat(outfile, filename, 256);
        CHKERRQ(ierr);

        ierr = VecPlaceArray(opflow->X, x + tcopflow->xstarti[t]);
        ierr = VecPlaceArray(opflow->Lambda, x + tcopflow->gstarti[t]);

        if (!opflow->solutiontops) {
          ierr = OPFLOWSolutionToPS(opflow);
          CHKERRQ(ierr);
        }
        ierr = OPFLOWSaveSolution(opflow, format, outfile);
        CHKERRQ(ierr);

        ierr = VecResetArray(opflow->X);
        ierr = VecResetArray(opflow->Lambda);
      }
      ierr = VecRestoreArrayRead(tcopflow->X, &x);
      CHKERRQ(ierr);
      ierr = VecRestoreArrayRead(tcopflow->Lambda, &lambda);
      CHKERRQ(ierr);
    } else {
      opflow = scopflow->opflows[c];
      if (!opflow->solutiontops) {
        ierr = SCOPFLOWGetSolution(scopflow, scopflow->cstart + c, &opflow->X);
        CHKERRQ(ierr);
        ierr = SCOPFLOWGetConstraintMultipliers(scopflow, scopflow->cstart + c,
                                                &opflow->Lambda);
        CHKERRQ(ierr);
        ierr = OPFLOWSolutionToPS(opflow);
        CHKERRQ(ierr);
      }
      ierr = PetscSNPrintf(filename, 64, "cont_%d", scopflow->cstart + c);
      CHKERRQ(ierr);
      ierr = PetscStrncpy(outfile, outdir, 256);
      CHKERRQ(ierr);
      ierr = PetscStrlcat(outfile, "/", 256);
      CHKERRQ(ierr);
      ierr = PetscStrlcat(outfile, filename, 256);
      CHKERRQ(ierr);

      ierr = OPFLOWSaveSolution(opflow, format, outfile);
      CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}
