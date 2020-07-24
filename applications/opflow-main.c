static char help[] = "User example calling OPFLOW.\n\n";

#include <opflow.h>
#include <scopflow_config.h>
#include <utils.h>

int main(int argc,char **argv)
{
  PetscErrorCode     ierr;
  OPFLOW             opflow;
  OutputFormat       fmt=MATPOWER;
  char               file[PETSC_MAX_PATH_LEN];
  PetscBool          flg;
  #if defined(PETSC_USE_LOG)
    PetscLogStage    stages[2];
  #endif
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

  /* PetscLogStageRegister */
  ierr = PetscLogStageRegister("Reading Data",&stages[0]);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Allocation & Setting",&stages[1]);CHKERRQ(ierr);

  /* Petsc Stage Push 1 */
  ierr = PetscLogStagePush(stages[0]);CHKERRQ(ierr);
  /* Create OPFLOW object */
  ierr = OPFLOWCreate(PETSC_COMM_WORLD,&opflow);CHKERRQ(ierr);

  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL,NULL,"-netfile",file,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  /* Read Network Data file */
  if(flg) {
    ierr = OPFLOWReadMatPowerData(opflow,file);CHKERRQ(ierr);
  } else {
    ierr = OPFLOWReadMatPowerData(opflow,"datafiles/case9mod.m");CHKERRQ(ierr);
  }
  /* End of First Stage and Start of Second */
  ierr = PetscLogStagePop();CHKERRQ(ierr);
  ierr = PetscLogStagePush(stages[1]);CHKERRQ(ierr);
  /* Solve */
  ierr = OPFLOWSolve(opflow);CHKERRQ(ierr);

  ierr = OPFLOWPrintSolution(opflow);CHKERRQ(ierr);

  ierr = OPFLOWSaveSolution(opflow,fmt,"opflowout");CHKERRQ(ierr);

  /*End of Final Stage */
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  /* Destroy OPFLOW object */
  ierr = OPFLOWDestroy(&opflow);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}
