static char help[] =
    "User example calling stochastic optimal power flow SOPFLOW.\n\n";

#include <exago_config.h>
#include <sopflow.h>
#include <utils.h>

int main(int argc, char **argv) {
  PetscErrorCode ierr;
  SOPFLOW sopflow;
  char file[PETSC_MAX_PATH_LEN];
  char scenfile[PETSC_MAX_PATH_LEN];
  char outputdir[PETSC_MAX_PATH_LEN];
  PetscBool outputdir_set;
  PetscBool flg = PETSC_FALSE, flgscen = PETSC_FALSE;
  PetscBool print_output = PETSC_FALSE, save_output = PETSC_FALSE;
  PetscLogStage stages[3];
  char appname[] = "sopflow";
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
  ierr = PetscOptionsGetString(NULL, NULL, "-sopflow_output_directory",
                               outputdir, PETSC_MAX_PATH_LEN, &outputdir_set);
  CHKERRQ(ierr);
  if (!outputdir_set) {
    strcpy(outputdir, "sopflowout");
  }

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

  /* Get scenario data file from command line */
  ierr = PetscOptionsGetString(NULL, NULL, "-scenfile", scenfile,
                               PETSC_MAX_PATH_LEN, &flgscen);
  CHKERRQ(ierr);

  /* Stage 1 - Application creation and reading data */
  ierr = PetscLogStagePush(stages[0]);
  CHKERRQ(ierr);

  /* Create SOPFLOW object */
  ierr = SOPFLOWCreate(PETSC_COMM_WORLD, &sopflow);
  CHKERRQ(ierr);

  /* Set Network Data file */
  if (flg) {
    ierr = SOPFLOWSetNetworkData(sopflow, file);
    CHKERRQ(ierr);
  } else {
    ierr = SOPFLOWSetNetworkData(sopflow, "datafiles/case9/case9mod.m");
    CHKERRQ(ierr);
  }

  /* Set Scenario Data file */
  if (flgscen) {
    ierr = SOPFLOWSetScenarioData(sopflow, SOPFLOW_NATIVE_SINGLEPERIOD, WIND,
                                  scenfile);
    CHKERRQ(ierr);
  }

  /* Set a subset of scenarios to be selected. Can use the option -sopflow_Ns
   * instead */
  /*   ierr = SOPFLOWSetNumScenarios(sopflow,2);CHKERRQ(ierr); */

  ierr = PetscLogStagePop();
  CHKERRQ(ierr);

  /* Stage 2 - Set up SOPFLOW application */
  ierr = PetscLogStagePush(stages[1]);
  CHKERRQ(ierr);

  /* Set up */
  ierr = SOPFLOWSetUp(sopflow);
  CHKERRQ(ierr);

  /* Stage 3 - Solve */
  ierr = PetscLogStagePush(stages[2]);
  CHKERRQ(ierr);

  /* Solve */
  ierr = SOPFLOWSolve(sopflow);
  CHKERRQ(ierr);

  ierr = PetscLogStagePop();
  CHKERRQ(ierr);

  /* Print solution */
  if (print_output) {
    ierr = SOPFLOWPrintSolution(sopflow, 0);
    CHKERRQ(ierr);
  }

  /* Save solution */
  if (save_output) {
    ierr = SOPFLOWSaveSolutionAllDefault(sopflow, outputdir);
    CHKERRQ(ierr);
  }

  /* Destroy SOPFLOW object */
  ierr = SOPFLOWDestroy(&sopflow);
  CHKERRQ(ierr);

  ExaGOFinalize();
  // PetscFinalize();
  return 0;
}
