static char help[] = "User example calling SCOPFLOW.\n\n";

#include <exago_config.h>
#include <scopflow.h>
#include <utils.h>

int main(int argc, char **argv) {
  PetscErrorCode ierr;
  SCOPFLOW scopflow;
  char file[PETSC_MAX_PATH_LEN];
  char ctgcfile[PETSC_MAX_PATH_LEN];
  PetscBool flg = PETSC_FALSE, flgctgc = PETSC_FALSE,
            ismultiperiod = PETSC_TRUE;
  PetscBool print_output = PETSC_FALSE, save_output = PETSC_FALSE;
  PetscLogStage stages[3];
  char appname[] = "scopflow";
  MPI_Comm comm = MPI_COMM_WORLD;

  /** Use `ExaGOLogSetLoggingFileName("opflow-logfile");` to log the output. */
  ierr = ExaGOInitialize(comm, &argc, &argv, appname, help);
  if (ierr) {
    fprintf(stderr, "Could not initialize ExaGO application %s.\n", appname);
    return ierr;
  }

  ierr = PetscOptionsGetBool(NULL, NULL, "-print_output", &print_output, NULL);
  CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL, NULL, "-save_output", &save_output, NULL);
  CHKERRQ(ierr);

  /* Register stages for profiling application code sections */
  ierr = PetscLogStageRegister("Reading Data", &stages[0]);
  CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Set up", &stages[1]);
  CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Solve", &stages[2]);
  CHKERRQ(ierr);

  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL, NULL, "-netfile", file, PETSC_MAX_PATH_LEN,
                               &flg);
  CHKERRQ(ierr);

  /* Get contingency data file from command line */
  ierr = PetscOptionsGetString(NULL, NULL, "-ctgcfile", ctgcfile,
                               PETSC_MAX_PATH_LEN, &flgctgc);
  CHKERRQ(ierr);

  /* Stage 1 - Application creation and reading data */
  ierr = PetscLogStagePush(stages[0]);
  CHKERRQ(ierr);

  /* Create SCOPFLOW object */
  ierr = SCOPFLOWCreate(MPI_COMM_WORLD, &scopflow);
  CHKERRQ(ierr);

  /* Set Network Data file */
  if (flg) {
    ierr = SCOPFLOWSetNetworkData(scopflow, file);
    CHKERRQ(ierr);
  } else {
    ierr = SCOPFLOWSetNetworkData(scopflow, "datafiles/case9/case9mod.m");
    CHKERRQ(ierr);
  }

  /* Set Contingency Data file */
  if (flgctgc) {
    std::string ext = FileNameExtension(ctgcfile);
    if (ext == "con") {
      ierr = SCOPFLOWSetContingencyData(scopflow, PSSE, ctgcfile);
    } else if (ext == "cont") {
      ierr = SCOPFLOWSetContingencyData(scopflow, NATIVE, ctgcfile);
    } else {
      std::stringstream errs;
      errs << "Unknown contingency data format: " << ctgcfile << ", (" << ext
           << ")";
      throw ExaGOError(errs.str().c_str());
    }
    CHKERRQ(ierr);
  }

  /* Set a subset of contingencies to be selected. Can use the option
   * -scopflow_Nc instead */
  /*   ierr = SCOPFLOWSetNumContingencies(scopflow,2);CHKERRQ(ierr); */

  ierr = PetscLogStagePop();
  CHKERRQ(ierr);

  /* Stage 2 - Set up SCOPFLOW application */
  ierr = PetscLogStagePush(stages[1]);
  CHKERRQ(ierr);

  /* Set up */
  ierr = SCOPFLOWSetUp(scopflow);
  CHKERRQ(ierr);

  /* Stage 3 - Solve */
  ierr = PetscLogStagePush(stages[2]);
  CHKERRQ(ierr);

  /* Solve */
  ierr = SCOPFLOWSolve(scopflow);
  CHKERRQ(ierr);

  ierr = PetscLogStagePop();
  CHKERRQ(ierr);

  /* Print solution */
  if (print_output) {
    ierr = SCOPFLOWPrintSolution(scopflow, 0);
    CHKERRQ(ierr);
  }

  /* Save solution */
  if (save_output) {
    ierr = SCOPFLOWSaveSolutionAll(scopflow, MATPOWER, "scopflowout");
    CHKERRQ(ierr);
  }

  /* Destroy SCOPFLOW object */
  ierr = SCOPFLOWDestroy(&scopflow);
  CHKERRQ(ierr);

  ExaGOFinalize();
  // PetscFinalize();
  return 0;
}
