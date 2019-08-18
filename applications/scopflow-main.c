static char help[] = "User example calling SCOPFLOW.\n\n";

#include <scopflow.h>


int main(int argc,char **argv)
{
  PetscErrorCode    ierr;
  SCOPFLOW          scopflow;
  char              file[PETSC_MAX_PATH_LEN];
  PetscBool         flg;
  #if defined(PETSC_USE_LOG)
    PetscLogStage stages[2];
  #endif

  PetscInitialize(&argc,&argv,"scopflowoptions",help);

  /* PetscLogStageRegister */
  ierr = PetscLogStageRegister("Reading Data",&stages[0]);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Solve",&stages[1]);CHKERRQ(ierr);

  /* Petsc Stage Push 1 */
  ierr = PetscLogStagePush(stages[0]);CHKERRQ(ierr);
  /* Create OPFLOW object */
  ierr = SCOPFLOWCreate(PETSC_COMM_WORLD,&scopflow);CHKERRQ(ierr);

  ierr = SCOPFLOWSetNumScenarios(scopflow,2);CHKERRQ(ierr);

  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL,NULL,"-netfile",file,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);

  /* Read Network Data file */
  if(flg) {
    ierr = SCOPFLOWSetNetworkData(scopflow,file);CHKERRQ(ierr);
  } else {
    ierr = SCOPFLOWSetNetworkData(scopflow,"datafiles/case9mod.m");CHKERRQ(ierr);
  }
  /* End of First Stage and Start of Second */
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  ierr = PetscLogStagePush(stages[1]);CHKERRQ(ierr);
  /* Solve */
  ierr = SCOPFLOWSolve(scopflow);CHKERRQ(ierr);
  /*End of Final Stage */
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  /* Destroy OPFLOW object */
  ierr = SCOPFLOWDestroy(&scopflow);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}
