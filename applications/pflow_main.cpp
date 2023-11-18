static char help[] =
    "User example calling PFLOW. Reads data in PSSE raw format\n\n";

#include <pflow.h>
#include <utils.h>
#include <exago_config.h>

int main(int argc, char **argv) {
  PetscErrorCode ierr;
  PFLOW pflow;
  char file[PETSC_MAX_PATH_LEN], gicfile[PETSC_MAX_PATH_LEN],
    output_name[PETSC_MAX_PATH_LEN];
  PetscBool flg = PETSC_FALSE, print_output = PETSC_FALSE,
    save_output = PETSC_FALSE;
  PetscLogStage read, setup, solve;
  MPI_Comm comm = MPI_COMM_WORLD;
  char appname[] = "pflow";
  const char default_output_name[] = "pflowout";

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

  ierr = PetscOptionsGetBool(NULL, NULL, "-print_output", &print_output, NULL);
  ExaGOCheckError(ierr);
  ierr = PetscOptionsGetString(NULL, NULL, "-output_file", output_name,
                               PETSC_MAX_PATH_LEN, &save_output);
  ExaGOCheckError(ierr);


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

  /* Get gic data file from command line */
  ierr = PetscOptionsGetString(NULL, NULL, "-gicfile", gicfile,
                               PETSC_MAX_PATH_LEN, &flg);
  ExaGOCheckError(ierr);
  if (flg) {
    ierr = PFLOWSetGICData(pflow, gicfile);
    ExaGOCheckError(ierr);
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

  if(print_output) {
    ierr = PFLOWPrintSolution(pflow);
    CHKERRQ(ierr);
  }

  if (save_output) {
    ierr = PFLOWSaveSolutionDefault(pflow, output_name);
    ExaGOCheckError(ierr);
  }

  
  /* Destroy PFLOW object */
  ierr = PFLOWDestroy(&pflow);
  CHKERRQ(ierr);

  ExaGOFinalize();
  return 0;
}
