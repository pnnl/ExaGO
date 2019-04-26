static char help[] = "User example calling DYNOPFLOW.\n\n";

#include <dynopflow.h>


int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  DYNOPFLOW      dynopflow;

  ierr = PetscInitialize(&argc,&argv,NULL,help);CHKERRQ(ierr);
  
  /* Create DYNOPFLOW object */
  ierr = DYNOPFLOWCreate(PETSC_COMM_WORLD,&dynopflow);CHKERRQ(ierr);

  /* Read Network Data file */
  ierr = DYNOPFLOWReadMatPowerData(dynopflow,"datafiles/solved_2016_05_20T14_19_26.m");CHKERRQ(ierr);

  /* Read dynamic model data file */
  ierr = DYNOPFLOWReadDyrData(dynopflow,"/Users/frank/tmp/tsced/data/tsopf-data/case118/case118.dyn");CHKERRQ(ierr);

  /* Read event data file */
  ierr = DYNOPFLOWReadEventData(dynopflow,"/Users/frank/tmp/tsced/data/tsopf-data/case118/case118_159.event");CHKERRQ(ierr);

  /* Set time step */
  ierr = DYNOPFLOWSetStartTimeAndStep(dynopflow,0.0,0.01);CHKERRQ(ierr);

  /* Set Duration */
  ierr = DYNOPFLOWSetDuration(dynopflow,10000,1.0);CHKERRQ(ierr);
  
  /* Solve */
  ierr = DYNOPFLOWSolve(dynopflow);CHKERRQ(ierr);

  /* Destroy DYNOPF object */
  ierr = DYNOPFLOWDestroy(&dynopflow);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}
  
