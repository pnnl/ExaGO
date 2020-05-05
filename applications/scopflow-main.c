static char help[] = "User example calling SCOPFLOW.\n\n";

#include <scopflow_config.h>
#include <scopflow.h>


int main(int argc,char **argv)
{
  PetscErrorCode    ierr;
  SCOPFLOW          scopflow;
  Vec               X;
  PetscBool         conv_status;
  PetscReal         obj=0.0;
  char              file[PETSC_MAX_PATH_LEN];
  char              ctgcfile[PETSC_MAX_PATH_LEN];
  PetscBool         flg=PETSC_FALSE,flgctgc=PETSC_FALSE;
  #if defined(PETSC_USE_LOG)
    PetscLogStage stages[2];
  #endif
  char options_pathname[200] = SCOPFLOW_OPTIONS_DIR;
  char* filename = "/scopflowoptions";
  printf("%s\n", options_pathname);
  printf("%s\n", filename);
  strcat(options_pathname, filename);
  printf("%s\n", options_pathname);

  PetscInitialize(&argc,&argv,options_pathname,help);

  /* PetscLogStageRegister */
  ierr = PetscLogStageRegister("Reading Data",&stages[0]);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Solve",&stages[1]);CHKERRQ(ierr);

  /* Petsc Stage Push 1 */
  ierr = PetscLogStagePush(stages[0]);CHKERRQ(ierr);

  /* Create SCOPFLOW object */
  ierr = SCOPFLOWCreate(PETSC_COMM_WORLD,&scopflow);CHKERRQ(ierr);

  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL,NULL,"-netfile",file,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);

  /* Get contingency data file from command line */
  ierr = PetscOptionsGetString(NULL,NULL,"-ctgcfile",ctgcfile,PETSC_MAX_PATH_LEN,&flgctgc);CHKERRQ(ierr);

  /* Set Network Data file */
  if(flg) {
    ierr = SCOPFLOWSetNetworkData(scopflow,file);CHKERRQ(ierr);
  } else {
    ierr = SCOPFLOWSetNetworkData(scopflow,"datafiles/case9mod.m");CHKERRQ(ierr);
  }

  /* Set Contingency Data file */
  if(flgctgc) {
    ierr = SCOPFLOWSetContingencyData(scopflow,ctgcfile);CHKERRQ(ierr);
  }
  
  /* Set a subset of scenarios to be selected. Can use the option -scopflow_Ns instead */
  /*   ierr = SCOPFLOWSetNumScenarios(scopflow,2);CHKERRQ(ierr); */

  /* End of First Stage and Start of Second */
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  ierr = PetscLogStagePush(stages[1]);CHKERRQ(ierr);
  /* Solve */
  ierr = SCOPFLOWSolve(scopflow);CHKERRQ(ierr);
  /*End of Final Stage */
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  /* Get convergence status */
  ierr = SCOPFLOWGetConvergenceStatus(scopflow,&conv_status);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"SCOPFLOW %s\n",conv_status?"converged":"did not converge");

  /* Get objective function */
  ierr = SCOPFLOWGetObjective(scopflow,&obj);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"SCOPFLOW objective = %lf\n",obj);CHKERRQ(ierr);

  /* Get solution */
  ierr = SCOPFLOWGetBaseCaseSolution(scopflow,&X);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"SCOPFLOW base case solution\n");
  ierr = VecView(X,0);CHKERRQ(ierr);

  /* Destroy SCOPFLOW object */
  ierr = SCOPFLOWDestroy(&scopflow);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}
