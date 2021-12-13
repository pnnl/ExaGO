static char help[] = "User example calling TCOPFLOW.\n\n";

#include <exago_config.h>
#include <tcopflow.h>
#include <tcopflowselfcheck.h>
#include <utils.hpp>

int main(int argc, char **argv) {
  PetscErrorCode ierr;
  TCOPFLOW tcopflow;
  char file[PETSC_MAX_PATH_LEN];
  char ploadprofile[PETSC_MAX_PATH_LEN];
  char qloadprofile[PETSC_MAX_PATH_LEN];
  char windgenprofile[PETSC_MAX_PATH_LEN];
  PetscBool flg = PETSC_FALSE, flg1 = PETSC_FALSE, flg2 = PETSC_FALSE,
            flg3 = PETSC_FALSE;
  PetscBool print_output = PETSC_FALSE, save_output = PETSC_FALSE;
#if defined(PETSC_USE_LOG)
  PetscLogStage stages[2];
#endif
  char filename[] = "/tcopflowoptions";
  MPI_Comm comm = MPI_COMM_WORLD;
  char appname[] = "tcopflow";

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

  /* PetscLogStageRegister */
  ierr = PetscLogStageRegister("Reading Data", &stages[0]);
  CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Solve", &stages[1]);
  CHKERRQ(ierr);

  /* Petsc Stage Push 1 */
  ierr = PetscLogStagePush(stages[0]);
  CHKERRQ(ierr);

  /* Create TCOPFLOW object */
  ierr = TCOPFLOWCreate(PETSC_COMM_WORLD, &tcopflow);
  CHKERRQ(ierr);

  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL, NULL, "-netfile", file, PETSC_MAX_PATH_LEN,
                               &flg);
  CHKERRQ(ierr);

  /* Set Network Data file */
  if (flg) {
    ierr = TCOPFLOWSetNetworkData(tcopflow, file);
    CHKERRQ(ierr);
  } else {
    ierr = TCOPFLOWSetNetworkData(tcopflow, "datafiles/case9/case9mod.m");
    CHKERRQ(ierr);
  }

  /* Set loadp and loadq files */
  ierr = PetscOptionsGetString(NULL, NULL, "-tcopflow_ploadprofile",
                               ploadprofile, PETSC_MAX_PATH_LEN, &flg1);
  CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL, NULL, "-tcopflow_qloadprofile",
                               qloadprofile, PETSC_MAX_PATH_LEN, &flg2);
  CHKERRQ(ierr);
  if (flg1 && flg2) {
    ierr = TCOPFLOWSetLoadProfiles(tcopflow, ploadprofile, qloadprofile);
    CHKERRQ(ierr);
  } else if (flg1) {
    ierr = TCOPFLOWSetLoadProfiles(tcopflow, ploadprofile, NULL);
    CHKERRQ(ierr);
  } else if (flg2) {
    ierr = TCOPFLOWSetLoadProfiles(tcopflow, NULL, qloadprofile);
    CHKERRQ(ierr);
  }

  ierr = PetscOptionsGetString(NULL, NULL, "-tcopflow_windgenprofile",
                               windgenprofile, PETSC_MAX_PATH_LEN, &flg3);
  CHKERRQ(ierr);
  if (flg3) {
    ierr = TCOPFLOWSetWindGenProfiles(tcopflow, windgenprofile);
    CHKERRQ(ierr);
  }

  /* End of First Stage and Start of Second */
  ierr = PetscLogStagePop();
  CHKERRQ(ierr);

  ierr = PetscLogStagePush(stages[1]);
  CHKERRQ(ierr);
  /* Solve */
  ierr = TCOPFLOWSolve(tcopflow);
  CHKERRQ(ierr);
  /*End of Final Stage */
  ierr = PetscLogStagePop();
  CHKERRQ(ierr);

  /* Print solution */
  if (print_output) {
    ierr = TCOPFLOWPrintSolution(tcopflow, 0);
    CHKERRQ(ierr);
  }

  /* Save solution */
  if (save_output) {
    ierr = TCOPFLOWSaveSolutionAll(tcopflow, MATPOWER, "tcopflowout");
    CHKERRQ(ierr);
  }

  int ret = ExaGOSelfcheckTCOPFLOW(tcopflow);

  /* Destroy OPFLOW object */
  ierr = TCOPFLOWDestroy(&tcopflow);
  CHKERRQ(ierr);

  ExaGOFinalize();
  return ret;
}
