static char help[] = "User example calling OPFLOW.\n\n";

#include <opflow.h>


int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  OPFLOW             opflow;
  char              file[PETSC_MAX_PATH_LEN];
  PetscBool         flg;

  PetscInitialize(&argc,&argv,"opflowoptions",help);
  
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

  /* Solve */
  ierr = OPFLOWSolve(opflow);CHKERRQ(ierr);

  /* Destroy OPFLOW object */
  ierr = OPFLOWDestroy(&opflow);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}
  
