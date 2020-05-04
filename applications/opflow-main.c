static char help[] = "User example calling OPFLOW.\n\n";

#include <opflow.h>
#include <scopflow_config.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  OPFLOW             opflow;
  char              file[PETSC_MAX_PATH_LEN];
  PetscBool         flg;
  PetscReal         obj;
  Vec               X,G,Lambda;
  PetscBool         conv_status;
  #if defined(PETSC_USE_LOG)
    PetscLogStage stages[2];
  #endif
  char options_pathname[200] = SCOPFLOW_OPTIONS_DIR;
  char* filename = "/opflowoptions";
  printf("%s\n", options_pathname);
  printf("%s\n", filename);
  strcat(options_pathname, filename);
  printf("%s\n", options_pathname);

  PetscInitialize(&argc,&argv,options_pathname,help);

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

  /* Get convergence status */
  ierr = OPFLOWGetConvergenceStatus(opflow,&conv_status);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"OPFLOW %s\n",conv_status?"converged":"did not converge");

  /* Get objective function */
  ierr = OPFLOWGetObjective(opflow,&obj);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"OPFLOW objective = %lf\n",obj);CHKERRQ(ierr);

  /* Get solution */
  ierr = OPFLOWGetSolution(opflow,&X);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"OPFLOW solution\n");
  ierr = VecView(X,0);CHKERRQ(ierr);

  /* Get constraints */
  ierr = OPFLOWGetConstraints(opflow,&G);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"OPFLOW constraints\n");
  ierr = VecView(G,0);CHKERRQ(ierr);

  /* Get constraint multipliers */
  ierr = OPFLOWGetConstraintMultipliers(opflow,&Lambda);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"OPFLOW constraint lagrange multipliers\n");
  ierr = VecView(Lambda,0);CHKERRQ(ierr);

  /*End of Final Stage */
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  /* Destroy OPFLOW object */
  ierr = OPFLOWDestroy(&opflow);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}
