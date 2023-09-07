static char help[] = "User example calling OPFLOW.\n\n";

#include <opflow.h>
#include <exago_config.h>
#include <utils.h>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  PetscErrorCode ierr;
  OPFLOW opflow;
  char file[PETSC_MAX_PATH_LEN], gicfile[PETSC_MAX_PATH_LEN],
      output_name[PETSC_MAX_PATH_LEN];
  PetscBool flg = PETSC_FALSE, print_output = PETSC_FALSE,
            save_output = PETSC_FALSE;
  PetscLogStage stages[3];
  char appname[] = "opflow";
  const char default_output_name[] = "opflowout";
  MPI_Comm comm = MPI_COMM_WORLD;
  PetscBool selfcheck = PETSC_FALSE;

  /** Use `ExaGOLogSetLoggingFileName("opflow-logfile");` to log the output. */
  ierr = ExaGOInitialize(comm, &argc, &argv, appname, help);
  if (ierr) {
    fprintf(stderr, "Could not initialize ExaGO application %s.\n", appname);
    return ierr;
  }

  ierr = PetscOptionsGetBool(NULL, NULL, "-print_output", &print_output, NULL);
  ExaGOCheckError(ierr);
  ierr = PetscOptionsGetString(NULL, NULL, "-save_output", &output_name[0],
                               PETSC_MAX_PATH_LEN, &save_output);
  ExaGOCheckError(ierr);

  /* Register stages for profiling application code sections */
  ierr = PetscLogStageRegister("Reading Data", &stages[0]);
  ExaGOCheckError(ierr);
  ierr = PetscLogStageRegister("Set up", &stages[1]);
  ExaGOCheckError(ierr);
  ierr = PetscLogStageRegister("Solve", &stages[2]);
  ExaGOCheckError(ierr);

  /* Stage 1 - Application creation and reading data */
  ierr = PetscLogStagePush(stages[0]);
  ExaGOCheckError(ierr);

  ExaGOLog(EXAGO_LOG_INFO, "{}", "Creating OPFlow\n");

  /* Create OPFLOW object */
  ierr = OPFLOWCreate(comm, &opflow);
  ExaGOCheckError(ierr);

  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL, NULL, "-netfile", file, PETSC_MAX_PATH_LEN,
                               &flg);
  ExaGOCheckError(ierr);
  if (flg) {
    ierr = OPFLOWReadMatPowerData(opflow, file);
    ExaGOCheckError(ierr);
  } else {
    ierr = OPFLOWReadMatPowerData(opflow, "datafiles/case9/case9mod.m");
    ExaGOCheckError(ierr);
  }

  /* Get gic data file from command line */
  ierr = PetscOptionsGetString(NULL, NULL, "-gicfile", gicfile,
                               PETSC_MAX_PATH_LEN, &flg);
  ExaGOCheckError(ierr);
  if (flg) {
    ierr = OPFLOWSetGICData(opflow, gicfile);
    ExaGOCheckError(ierr);
  }

  ierr = PetscLogStagePop();
  ExaGOCheckError(ierr);

  /* Stage 2 - Set up OPFLOW application */
  ierr = PetscLogStagePush(stages[1]);
  ExaGOCheckError(ierr);

  /* Set up */
  ierr = OPFLOWSetUp(opflow);
  ExaGOCheckError(ierr);

  ierr = PetscLogStagePop();
  ExaGOCheckError(ierr);

  /* Stage 3 - Solve */
  ierr = PetscLogStagePush(stages[2]);
  ExaGOCheckError(ierr);

  /* Solve */
  ierr = OPFLOWSolve(opflow);
  ExaGOCheckError(ierr);

  ierr = PetscLogStagePop();
  ExaGOCheckError(ierr);

  if (print_output) {
    ierr = OPFLOWPrintSolution(opflow);
    ExaGOCheckError(ierr);
  }

  if (save_output) {
    // deal with "-save_output 0/1" (the way it used to be)
    if (strlen(output_name) == 1 && output_name[0] == '1') {
      strncpy(output_name, default_output_name, PETSC_MAX_PATH_LEN);
    } else if (strlen(output_name) == 1 && output_name[0] == '0') {
      save_output = PETSC_FALSE;
    } else if (strlen(output_name) == 0) {
      strncpy(output_name, default_output_name, PETSC_MAX_PATH_LEN);
    }
    if (save_output) {
      ierr = OPFLOWSaveSolutionDefault(opflow, output_name);
      ExaGOCheckError(ierr);
    }
  }

  /* Destroy OPFLOW object */
  ierr = OPFLOWDestroy(&opflow);
  ExaGOCheckError(ierr);

  ExaGOFinalize();
  return 0;
}
