static char help[] = "User example calling OPFLOW.\n\n";

#include <opflow.h>
#include <exago_config.h>
#include <utils.h>

int main(int argc,char **argv)
{
  PetscErrorCode     ierr;
  OPFLOW             opflow;
  OutputFormat       fmt=MATPOWER;
  char               file[PETSC_MAX_PATH_LEN];
  PetscBool          flg=PETSC_FALSE,print_output=PETSC_FALSE,save_output=PETSC_FALSE;
  PetscLogStage    stages[3];
  char options_pathname[200] = EXAGO_OPTIONS_DIR;
  char filename[] = "/opflowoptions";
  strcat(options_pathname, filename);

  const int npths = 3;
  char** pths = malloc(PETSC_MAX_PATH_LEN * npths);
  pths[0] = strdup(options_pathname);
  pths[1] = strdup("./opflowoptions");
  pths[2] = strdup("./options/opflowoptions");
  const int idx = doFilesExist(pths, npths);
  if (idx < 0)
  {
    fputs("Could not options file", stderr);
    return 1;
  }

  PetscInitialize(&argc,&argv,pths[idx],help);

  ierr = PetscOptionsGetBool(NULL,NULL,"-print_output",&print_output,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-save_output",&save_output,NULL);CHKERRQ(ierr);
  
  /* Register stages for profiling application code sections */
  ierr = PetscLogStageRegister("Reading Data",&stages[0]);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Set up",&stages[1]);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Solve",&stages[2]);CHKERRQ(ierr);

  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL,NULL,"-netfile",file,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);

  /* Stage 1 - Application creation and reading data */
  ierr = PetscLogStagePush(stages[0]);CHKERRQ(ierr);

  /* Create OPFLOW object */
  ierr = OPFLOWCreate(PETSC_COMM_WORLD,&opflow);CHKERRQ(ierr);
  
  /* Read Network Data file */
  if(flg) {
    ierr = OPFLOWReadMatPowerData(opflow,file);CHKERRQ(ierr);
  } else {
    ierr = OPFLOWReadMatPowerData(opflow,"datafiles/case9/case9mod.m");CHKERRQ(ierr);
  }

  ierr = PetscLogStagePop();CHKERRQ(ierr);

  /* Stage 2 - Set up OPFLOW application */
  ierr = PetscLogStagePush(stages[1]);CHKERRQ(ierr);

  /* Set up */
  ierr = OPFLOWSetUp(opflow);CHKERRQ(ierr);

  ierr = PetscLogStagePop();CHKERRQ(ierr);

  /* Stage 3 - Solve */
  ierr = PetscLogStagePush(stages[2]);CHKERRQ(ierr);

  /* Solve */
  ierr = OPFLOWSolve(opflow);CHKERRQ(ierr);

  ierr = PetscLogStagePop();CHKERRQ(ierr);
  
  if(print_output) {
    ierr = OPFLOWPrintSolution(opflow);CHKERRQ(ierr);
  }

  if(save_output) {
    ierr = OPFLOWSaveSolution(opflow,fmt,"opflowout");CHKERRQ(ierr);
  }

  /* Destroy OPFLOW object */
  ierr = OPFLOWDestroy(&opflow);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}
