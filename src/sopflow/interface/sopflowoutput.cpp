#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>
#include <private/sopflowimpl.h>

/*
  SOPFLOWPrintSolution - Prints SOPFLOW solution to stdout for the given
scenario

  Input Parameters:
+ sopflow - the SOPFLOW object
- scen_num - the scenario number (0 for base)
  Notes:
   Only prints out the system summary for base case

*/
PetscErrorCode SOPFLOWPrintSolution(SOPFLOW sopflow, PetscInt scen_num) {
  PetscErrorCode ierr;
  PetscBool conv_status;
  PetscReal cost;
  PetscInt c;
  SCOPFLOW scopflow;
  OPFLOW opflow;
  const PetscScalar *x, *lambda;

  PetscFunctionBegin;
  MPI_Barrier(sopflow->comm->type);
  c = scen_num - sopflow->sstart;
  if (scen_num >= sopflow->sstart && scen_num < sopflow->send) {
    if (!sopflow->ismulticontingency) {
      opflow = sopflow->opflows[c];
      if (!opflow->solutiontops) {
        ierr = SOPFLOWGetSolution(sopflow, scen_num, &opflow->X);
        CHKERRQ(ierr);
        ierr =
            SOPFLOWGetConstraintMultipliers(sopflow, scen_num, &opflow->Lambda);
        CHKERRQ(ierr);
        ierr = OPFLOWSolutionToPS(opflow);
        CHKERRQ(ierr);
      }
    } else {
      scopflow = sopflow->scopflows[c];
      if (!scopflow->ismultiperiod) {
        ierr = SOPFLOWGetSolution(sopflow, scen_num, &scopflow->X);
        CHKERRQ(ierr);
        ierr = SOPFLOWGetConstraintMultipliers(sopflow, scen_num,
                                               &scopflow->Lambda);
        CHKERRQ(ierr);

        ierr = VecGetArrayRead(scopflow->X, &x);
        CHKERRQ(ierr);
        ierr = VecGetArrayRead(scopflow->Lambda, &lambda);
        CHKERRQ(ierr);

        opflow = scopflow->opflows[0];

        ierr = OPFLOWComputeObjective(opflow, opflow->X, &opflow->obj);

        ierr = VecPlaceArray(opflow->X, x);
        CHKERRQ(ierr);
        ierr = VecPlaceArray(opflow->Lambda, lambda);
        CHKERRQ(ierr);

        if (!opflow->solutiontops) {
          ierr = OPFLOWSolutionToPS(opflow);
        }

        ierr = VecResetArray(opflow->X);
        CHKERRQ(ierr);
        ierr = VecResetArray(opflow->Lambda);
        CHKERRQ(ierr);

        ierr = VecRestoreArrayRead(scopflow->X, &x);
        CHKERRQ(ierr);
        ierr = VecRestoreArrayRead(scopflow->Lambda, &lambda);
        CHKERRQ(ierr);
      } else {
        SETERRQ(
            PETSC_COMM_SELF, PETSC_ERR_SUP,
            "SOPFLOWPrintOutput: multiperiod scopf case not implemented yet\n");
      }
    }

    /* Print to stdout */
    ierr = PetscPrintf(
        sopflow->comm->type,
        "=============================================================\n");
    CHKERRQ(ierr);
    ierr =
        PetscPrintf(sopflow->comm->type, "\tStochastic Optimal Power Flow\n");
    CHKERRQ(ierr);
    ierr = PetscPrintf(
        sopflow->comm->type,
        "=============================================================\n");
    CHKERRQ(ierr);

    ierr = PetscPrintf(sopflow->comm->type, "%-35s %d\n", "Number of scenarios",
                       sopflow->Ns);
    CHKERRQ(ierr);
    ierr = PetscPrintf(sopflow->comm->type, "%-35s %s\n",
                       "Multi-contingency scenarios?",
                       (sopflow->flatten_contingencies ||
                        sopflow->ismulticontingency || (sopflow->Nc - 1))
                           ? "YES"
                           : "NO");
    CHKERRQ(ierr);
    if (sopflow->ismulticontingency || sopflow->flatten_contingencies) {
      int nc;
      if (sopflow->flatten_contingencies) {
        nc = sopflow->Nc - 1;
      } else {
        nc = sopflow->scopflows[0]->Nc;
      }
      ierr = PetscPrintf(sopflow->comm->type, "%-35s %d\n",
                         "Contingencies per scenario", nc);
      CHKERRQ(ierr);
    }

    ierr = PetscPrintf(sopflow->comm->type, "%-35s %s\n", "Solver",
                       sopflow->solvername);
    CHKERRQ(ierr);

    ierr = PetscPrintf(sopflow->comm->type, "%-35s %s\n", "Initialization",
                       OPFLOWInitializationTypes[opflow->initializationtype]);
    CHKERRQ(ierr);
    ierr = PetscPrintf(sopflow->comm->type, "%-35s %s\n", "Load loss allowed",
                       opflow->include_loadloss_variables ? "YES" : "NO");
    CHKERRQ(ierr);
    if (opflow->include_loadloss_variables) {
      ierr = PetscPrintf(sopflow->comm->type, "%-35s %g\n",
                         "Load loss penalty ($)", opflow->loadloss_penalty);
      CHKERRQ(ierr);
    }

    ierr = PetscPrintf(sopflow->comm->type, "%-35s %s\n",
                       "Power imbalance allowed",
                       opflow->include_powerimbalance_variables ? "YES" : "NO");
    CHKERRQ(ierr);
    if (opflow->include_powerimbalance_variables) {
      ierr = PetscPrintf(sopflow->comm->type, "%-35s %g\n",
                         "Power imbalance penalty ($)",
                         opflow->powerimbalance_penalty);
      CHKERRQ(ierr);
    }

    ierr = PetscPrintf(sopflow->comm->type, "%-35s %s\n",
                       "Ignore line flow constraints",
                       opflow->ignore_lineflow_constraints ? "YES" : "NO");
    CHKERRQ(ierr);
    ierr = PetscPrintf(sopflow->comm->type, "\n");
    CHKERRQ(ierr);

#if 0 
    ierr = PetscPrintf(sopflow->comm->type, "%-35s %d\n", "Number of variables",
                       sopflow->Nx);
    CHKERRQ(ierr);
    ierr = PetscPrintf(sopflow->comm->type, "%-35s %d\n",
                       "Number of equality constraints", sopflow->Nconeq);
    CHKERRQ(ierr);
    ierr = PetscPrintf(sopflow->comm->type, "%-35s %d\n",
                       "Number of inequality constraints", sopflow->Nconineq);
    CHKERRQ(ierr);
    ierr = PetscPrintf(sopflow->comm->type, "%-35s %d\n",
                       "Number of coupling constraints", sopflow->Nconcoup);
    CHKERRQ(ierr);

    ierr = PetscPrintf(sopflow->comm->type, "\n");
    CHKERRQ(ierr);
#endif
  }
  MPI_Barrier(sopflow->comm->type);
  ierr = SOPFLOWGetConvergenceStatus(sopflow, &conv_status);
  CHKERRQ(ierr);
  ierr = PetscPrintf(sopflow->comm->type, "%-35s %s\n", "Convergence status",
                     conv_status ? "CONVERGED" : "DID NOT CONVERGE");
  CHKERRQ(ierr);
  ierr = SOPFLOWGetBaseObjective(sopflow, &cost);
  CHKERRQ(ierr);
  ierr = PetscPrintf(sopflow->comm->type, "%-35s %-7.2f\n",
                     "Objective value (base)", cost);
  CHKERRQ(ierr);
  ierr = PetscPrintf(sopflow->comm->type, "\n");
  CHKERRQ(ierr);

  if (scen_num >= sopflow->sstart && scen_num < sopflow->send) {
    ierr = PSPrintSystemSummary(opflow->ps);
    CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/*
   SOPFLOWAppendScenario - Appends scenario information to MINIMAL output

*/
static PetscErrorCode SOPFLOWAppendScenario(int scen_num,
                                            const char outfile[]) {
  PetscErrorCode ierr = 0;
  FILE *fd;
  char filename[PETSC_MAX_PATH_LEN];
  char ext[] = ".txt";

  PetscFunctionBegin;

  /* FIXME: This repeats code in PSSaveSolution_MINIMAL */

  strcpy(filename, outfile);
  /* Add extension to file name */
  ierr = PetscStrlcat(filename, ext, 256);
  CHKERRQ(ierr);

  fd = fopen(filename, "a");
  if (fd == NULL) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
            "Cannot (re)open output file %s", filename);
    CHKERRQ(ierr);
  }
  fprintf(fd, "Scenario: %d\n", scen_num);
  CHKERRQ(ierr);
  fclose(fd);
  PetscFunctionReturn(0);
}

/*
  SOPFLOWSaveSolutionBase - Saves the SOPFLOW solution for the given scenario to
file

  Input Parameters:
+ sopflow - the OPFLOW object
. scen_num - the scenario number (0 for base case)
. use_default_format - ignore format and let OPFLOW determine format
. format - the output file format (csv, matpower)
- outfile  - Name of output file
*/
static PetscErrorCode SOPFLOWSaveSolutionBase(SOPFLOW sopflow,
                                              PetscInt scen_num,
                                              PetscBool use_default_format,
                                              OutputFormat format,
                                              const char outfile[]) {
  PetscErrorCode ierr;
  OPFLOW opflow = sopflow->opflows[scen_num];
  bool is_minimal_output;

  PetscFunctionBegin;

  if (!opflow->solutiontops) {
    ierr = SOPFLOWGetSolution(sopflow, scen_num, &opflow->X);
    CHKERRQ(ierr);
    ierr = SOPFLOWGetConstraintMultipliers(sopflow, scen_num, &opflow->Lambda);
    CHKERRQ(ierr);
    ierr = OPFLOWSolutionToPS(opflow);
    CHKERRQ(ierr);
  }

  if (use_default_format) {
    ierr = OPFLOWSaveSolutionDefault(opflow, outfile);
    is_minimal_output = (opflow->outputformat == MINIMAL);
  } else {
    ierr = OPFLOWSaveSolution(opflow, format, outfile);
    is_minimal_output = (format == MINIMAL);
  }
  CHKERRQ(ierr);

  /* For MINIMAL format only, *append* scenario information to
     OPFLOW solution output */
  if (is_minimal_output) {
    int snum;
    snum = sopflow->scen_num[scen_num];
    ierr = SOPFLOWAppendScenario(snum, outfile);
    CHKERRQ(ierr);
    if (sopflow->Nc > 1) {
      int cnum;
      Contingency *thecon;
      cnum = sopflow->cont_num[scen_num];
      thecon = &(sopflow->ctgclist->cont[cnum]);
      ierr = ContigencyAppendMinimal(thecon, cnum, PSSE, outfile);
      CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(0);
}

/*
  SOPFLOWSaveSolution - Saves the SOPFLOW solution for the given scenario to
file with the specified format

  Input Parameters:
+ sopflow - the OPFLOW object
. scen_num - the scenario number (0 for base case)
. format - the output file format (csv, matpower)
- outfile  - Name of output file
*/
PetscErrorCode SOPFLOWSaveSolution(SOPFLOW sopflow, PetscInt scen_num,
                                   OutputFormat format, const char outfile[]) {

  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr =
      SOPFLOWSaveSolutionBase(sopflow, scen_num, PETSC_FALSE, format, outfile);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  SOPFLOWSaveSolutionDefault - Saves the SOPFLOW solution for the given scenario
to file with the default format

  Input Parameters:
+ sopflow - the OPFLOW object
. scen_num - the scenario number (0 for base case)
- outfile  - Name of output file
*/
PetscErrorCode SOPFLOWSaveSolutionDefault(SOPFLOW sopflow, PetscInt scen_num,
                                          const char outfile[]) {

  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr =
      SOPFLOWSaveSolutionBase(sopflow, scen_num, PETSC_TRUE, MATPOWER, outfile);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SOPFLOWSaveSolutionAll - Saves all SOPFLOW solutions

  Input Parameters:
+ sopflow - the SOPFLOW object
. use_default_format - ignore format and use default OPFLOW format
. format - the output file format (csv, matpower)
- outdir  - Name of output directory

  Notes:
   Save all SOPFLOW solutions (one scenario per file)
*/
static PetscErrorCode SOPFLOWSaveSolutionAllBase(SOPFLOW sopflow,
                                                 PetscBool use_default_format,
                                                 OutputFormat format,
                                                 const char outdir[]) {
  PetscErrorCode ierr;
  SCOPFLOW scopflow;
  OPFLOW opflow;
  PetscInt s;
  char filename[64];
  char outfile[256];
  char scopflowdirname[64];
  char scopflowdir[256];
  PetscScalar *x, *lambda;
  bool is_minimal_output;

  PetscFunctionBegin;

  ierr = PetscLogEventBegin(sopflow->outputlogger, 0, 0, 0, 0);
  CHKERRQ(ierr);

  if (!sopflow->comm->rank) {
    ierr = PetscMkdir(outdir);
    CHKERRQ(ierr);
  }
  MPI_Barrier(sopflow->comm->type);

  for (s = 0; s < sopflow->ns; s++) {
    if (!sopflow->ismulticontingency) {
      opflow = sopflow->opflows[s];
      if (!opflow->solutiontops) {
        ierr = SOPFLOWGetSolution(sopflow, sopflow->sstart + s, &opflow->X);
        CHKERRQ(ierr);
        ierr = SOPFLOWGetConstraintMultipliers(sopflow, sopflow->sstart + s,
                                               &opflow->Lambda);
        CHKERRQ(ierr);
        ierr = OPFLOWSolutionToPS(opflow);
        CHKERRQ(ierr);
      }
      ierr = PetscSNPrintf(filename, 64, "scen_%d", sopflow->sstart + s);
      CHKERRQ(ierr);
      ierr = PetscStrncpy(outfile, outdir, 256);
      CHKERRQ(ierr);
      ierr = PetscStrlcat(outfile, "/", 256);
      CHKERRQ(ierr);
      ierr = PetscStrlcat(outfile, filename, 256);
      CHKERRQ(ierr);

      if (use_default_format) {
        ierr = OPFLOWSaveSolutionDefault(opflow, outfile);
        is_minimal_output = (opflow->outputformat == MINIMAL);
      } else {
        ierr = OPFLOWSaveSolution(opflow, format, outfile);
        is_minimal_output = (format == MINIMAL);
      }
      CHKERRQ(ierr);

      /* For MINIMAL format only, *append* scenario information to
         OPFLOW solution output */
      if (is_minimal_output) {
        int snum;
        snum = sopflow->scen_num[s];
        ierr = SOPFLOWAppendScenario(snum, outfile);
        CHKERRQ(ierr);
        if (sopflow->Nc > 1) {
          int cnum;
          Contingency *thecon;
          cnum = sopflow->cont_num[s];
          thecon = &(sopflow->ctgclist->cont[cnum]);
          ierr = ContigencyAppendMinimal(thecon, cnum, PSSE, outfile);
          CHKERRQ(ierr);
        }
      }
    } else {
      scopflow = sopflow->scopflows[s];

      ierr = PetscSNPrintf(scopflowdirname, 64, "scen_%d", sopflow->sstart + s);
      CHKERRQ(ierr);
      ierr = PetscStrncpy(scopflowdir, outdir, 256);
      CHKERRQ(ierr);
      ierr = PetscStrlcat(scopflowdir, "/", 256);
      CHKERRQ(ierr);
      ierr = PetscStrlcat(scopflowdir, scopflowdirname, 256);
      CHKERRQ(ierr);

      ierr = SOPFLOWGetSolution(sopflow, sopflow->sstart + s, &scopflow->X);
      CHKERRQ(ierr);
      ierr = SOPFLOWGetConstraintMultipliers(sopflow, sopflow->sstart + s,
                                             &scopflow->Lambda);
      CHKERRQ(ierr);

      if (use_default_format) {
        ierr = SCOPFLOWSaveSolutionAllDefault(scopflow, scopflowdir);
      } else {
        ierr = SCOPFLOWSaveSolutionAll(scopflow, format, scopflowdir);
      }
      CHKERRQ(ierr);
    }
  }

  ierr = PetscLogEventEnd(sopflow->outputlogger, 0, 0, 0, 0);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  SOPFLOWSaveSolutionAll - Saves all SOPFLOW solutions

  Input Parameters:
+ sopflow - the SOPFLOW object
. format - the output file format (csv, matpower)
- outdir  - Name of output directory

  Notes:
   Save all SOPFLOW solutions (one scenario per file)
*/
PetscErrorCode SOPFLOWSaveSolutionAll(SOPFLOW sopflow, OutputFormat format,
                                      const char outdir[]) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = SOPFLOWSaveSolutionAllBase(sopflow, PETSC_FALSE, format, outdir);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  SOPFLOWSaveSolutionAll - Saves all SOPFLOW solutions using the default OPFLOW
format

  Input Parameters:
+ sopflow - the SOPFLOW object
- outdir  - Name of output directory

  Notes:
   Save all SOPFLOW solutions (one scenario per file)
*/
PetscErrorCode SOPFLOWSaveSolutionAllDefault(SOPFLOW sopflow,
                                             const char outdir[]) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = SOPFLOWSaveSolutionAllBase(sopflow, PETSC_TRUE, MATPOWER, outdir);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
