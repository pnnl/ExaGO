static char help[] = "User example calling PFLOW. Reads data in PSSE raw format\n\n";

#include <pflow.h>


int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PFLOW             pflow;
  char              file[PETSC_MAX_PATH_LEN];
  PetscBool         flg;

  PetscInitialize(&argc,&argv,"pflowoptions",help);
  
  /* Create PFLOW object */
  ierr = PFLOWCreate(PETSC_COMM_WORLD,&pflow);CHKERRQ(ierr);

  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL,NULL,"-netfile",file,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  /* Read Network Data file */
  if(flg) {
    if(strstr(file,".raw") != NULL) {
      ierr = PFLOWReadPSSERawData(pflow,file);CHKERRQ(ierr);
    } else {
      ierr = PFLOWReadMatPowerData(pflow,file);CHKERRQ(ierr);
    }
  } else {
    ierr = PFLOWReadMatPowerData(pflow,"datafiles/case9mod.m");CHKERRQ(ierr);
  }

  /* Solve */
  ierr = PFLOWSolve(pflow);CHKERRQ(ierr);

  /* Update line flows, Pgen, Qgen, and other parameters */
  ierr = PFLOWPostSolve(pflow);CHKERRQ(ierr);

  /* Destroy PFLOW object */
  ierr = PFLOWDestroy(&pflow);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}
  
