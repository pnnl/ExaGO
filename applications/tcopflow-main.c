static char help[] = "User example calling TCOPFLOW.\n\n";

#include <scopflow_config.h>
#include <tcopflow.h>


int main(int argc,char **argv)
{
  PetscErrorCode    ierr;
  TCOPFLOW          tcopflow;
  char              file[PETSC_MAX_PATH_LEN];
  char              ploadprofile[PETSC_MAX_PATH_LEN];
  char              qloadprofile[PETSC_MAX_PATH_LEN];
  char              windgenprofile[PETSC_MAX_PATH_LEN];
  PetscBool         flg=PETSC_FALSE,flg1=PETSC_FALSE,flg2=PETSC_FALSE,flg3=PETSC_FALSE;
  #if defined(PETSC_USE_LOG)
    PetscLogStage stages[2];
  #endif
  char options_pathname[200] = SCOPFLOW_OPTIONS_DIR;
  char filename[] = "/tcopflowoptions";
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

  /* Create TCOPFLOW object */
  ierr = TCOPFLOWCreate(PETSC_COMM_WORLD,&tcopflow);CHKERRQ(ierr);

  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL,NULL,"-netfile",file,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);

  /* Set Network Data file */
  if(flg) {
    ierr = TCOPFLOWSetNetworkData(tcopflow,file);CHKERRQ(ierr);
  } else {
    ierr = TCOPFLOWSetNetworkData(tcopflow,"datafiles/case9mod.m");CHKERRQ(ierr);
  }


  /* Set loadp and loadq files */
  ierr = PetscOptionsGetString(NULL,NULL,"-tcopflow_ploadprofile",ploadprofile,PETSC_MAX_PATH_LEN,&flg1);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL,NULL,"-tcopflow_qloadprofile",qloadprofile,PETSC_MAX_PATH_LEN,&flg2);CHKERRQ(ierr);
  if(flg1 && flg2) {
    ierr = TCOPFLOWSetLoadProfiles(tcopflow,ploadprofile,qloadprofile);CHKERRQ(ierr);
  } else if(flg1) {
    ierr = TCOPFLOWSetLoadProfiles(tcopflow,ploadprofile,NULL);CHKERRQ(ierr);
  } else if(flg2) {
    ierr = TCOPFLOWSetLoadProfiles(tcopflow,NULL,qloadprofile);CHKERRQ(ierr);
  }

  ierr = PetscOptionsGetString(NULL,NULL,"-tcopflow_windgenprofile",windgenprofile,PETSC_MAX_PATH_LEN,&flg3);CHKERRQ(ierr);
  if(flg3) {
    ierr = TCOPFLOWSetWindGenProfiles(tcopflow,windgenprofile);CHKERRQ(ierr);
  }

  /* End of First Stage and Start of Second */
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  ierr = PetscLogStagePush(stages[1]);CHKERRQ(ierr);
  /* Solve */
  ierr = TCOPFLOWSolve(tcopflow);CHKERRQ(ierr);
  /*End of Final Stage */
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  /* Print solution */
  ierr = TCOPFLOWPrintSolution(tcopflow,0);CHKERRQ(ierr);

  ierr = TCOPFLOWSaveSolutionAll(tcopflow,CSV,"tcopflowout");CHKERRQ(ierr);

  /* Destroy OPFLOW object */
  ierr = TCOPFLOWDestroy(&tcopflow);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}
