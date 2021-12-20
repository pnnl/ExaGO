static char help[] =
    "User example calling PFLOW. Reads data in PSSE raw format\n\n";

#include <pflow.h>
#include <utils.h>
#include <exago_config.h>

int main(int argc, char **argv) {
  PetscErrorCode ierr;
  PFLOW pflow;
  char file[PETSC_MAX_PATH_LEN];
  PetscBool flg;
  PetscLogStage read, setup, solve;
  MPI_Comm comm = MPI_COMM_WORLD;
  char appname[] = "pflow";

  ierr = ExaGOInitialize(comm, &argc, &argv, appname, help);
  if (ierr) {
    fprintf(stderr, "Could not initialize ExaGO application %s\n", appname);
    return ierr;
  }

  ierr = PetscLogStageRegister("ReadData", &read);
  CHKERRQ(ierr);
  ierr = PetscLogStageRegister("SetUp", &setup);
  CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Solve", &solve);
  CHKERRQ(ierr);

  /* Create PFLOW object */
  ierr = PFLOWCreate(PETSC_COMM_WORLD, &pflow);
  CHKERRQ(ierr);

  ierr = PetscLogStagePush(read);
  CHKERRQ(ierr);
  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL, NULL, "-netfile", file, PETSC_MAX_PATH_LEN,
                               &flg);
  CHKERRQ(ierr);
  /* Read Network Data file */
  if (flg) {
    if (strstr(file, ".raw") != NULL) {
      ierr = PFLOWReadPSSERawData(pflow, file);
      CHKERRQ(ierr);
    } else {
      ierr = PFLOWReadMatPowerData(pflow, file);
      CHKERRQ(ierr);
    }
  } else {
    ierr = PFLOWReadMatPowerData(pflow, "datafiles/case9/case9mod.m");
    CHKERRQ(ierr);
  }
  PetscLogStagePop();

  ierr = PetscLogStagePush(setup);
  CHKERRQ(ierr);
  /* Set up */
  ierr = PFLOWSetUp(pflow);
  CHKERRQ(ierr);
  PetscLogStagePop();

  ierr = PetscLogStagePush(solve);
  CHKERRQ(ierr);
  /* Solve */
  ierr = PFLOWSolve(pflow);
  CHKERRQ(ierr);
  PetscLogStagePop();

  /* Update line flows, Pgen, Qgen, and other parameters */
  //  ierr = PFLOWPostSolve(pflow);CHKERRQ(ierr);

  /* Destroy PFLOW object */
  ierr = PFLOWDestroy(&pflow);
  CHKERRQ(ierr);

  ExaGOFinalize();
  return 0;
}
