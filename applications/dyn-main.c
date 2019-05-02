static char help[] = "User example calling DYN.\n";

#include <dyn.h>
#include <private/dynimpl.h>

int main(int argc,char **argv)
{
  PetscErrorCode  ierr;
  DYN             dyn;
  PetscBool      flg;
  char           netfile[PETSC_MAX_PATH_LEN],dyrfile[PETSC_MAX_PATH_LEN],eventfile[PETSC_MAX_PATH_LEN];

  ierr = PetscInitialize(&argc,&argv,"dynoptions",help);CHKERRQ(ierr);

  /* Create DYN object */
  ierr = DYNCreate(PETSC_COMM_WORLD,&dyn);CHKERRQ(ierr);

  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL,NULL,"-netfile",netfile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);

  /* Read Network Data file */
  if(flg) {
    ierr = DYNReadMatPowerData(dyn,netfile);CHKERRQ(ierr);
  } else {
    ierr = DYNReadMatPowerData(dyn,"datafiles/case9mod.m");CHKERRQ(ierr);
  }

  /* Get dyr data file from command line */
  ierr = PetscOptionsGetString(NULL,NULL,"-dyrfile",dyrfile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  /* Read Machine Data file */
  if(flg) {
    ierr = DYNReadDyrData(dyn,dyrfile);CHKERRQ(ierr);
  } else {
    ierr = DYNReadDyrData(dyn,"datafiles/case9mod1.dyr");CHKERRQ(ierr);
  }

  /* Get event data file from command line */
  ierr = PetscOptionsGetString(NULL,NULL,"-eventfile",eventfile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  /* Read event data file */
  if(flg) {
    ierr = DYNReadEventData(dyn,eventfile);CHKERRQ(ierr);
  } else {
    ierr = DYNReadEventData(dyn,"datafiles/case9mod.event");CHKERRQ(ierr);
  }

  /* Set time step */
  ierr = DYNSetStartTimeAndStep(dyn,0.0,0.01);CHKERRQ(ierr);

  /* Set Duration */
  ierr = DYNSetDuration(dyn,10000,1.0);CHKERRQ(ierr);

  /* Initialize to steady-state */
  ierr = DYNSetComputePrePflow(dyn,PETSC_TRUE);CHKERRQ(ierr);

  /* Solve */
  ierr = DYNSolve(dyn);CHKERRQ(ierr);

  ierr = DYNDestroy(&dyn);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}

