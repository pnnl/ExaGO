static char help[] = "User example calling PFLOW and changing branch status.\n \
                      In this example, power flow on a 9-bus case is run first.\n\
                      This is followed by tripping two lines that result in the\n\
                      network splitting into two islands. A new reference bus is\n\
                      set for the second island, since the original reference bus\n\
                      is in island 1.\n\n";



#include <pflow.h>
#include <exago_config.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PFLOW             pflow;
  char options_pathname[200] = EXAGO_OPTIONS_DIR;
  char filename[] = "/pflowoptions";
  printf("%s\n", options_pathname);
  printf("%s\n", filename);
  strcat(options_pathname, filename);
  printf("%s\n", options_pathname);

  PetscInitialize(&argc,&argv,options_pathname,help);
  
  /* Create PFLOW object */
  ierr = PFLOWCreate(PETSC_COMM_WORLD,&pflow);CHKERRQ(ierr);

  ierr = PFLOWReadMatPowerData(pflow,"datafiles/case9mod.m");CHKERRQ(ierr);

  /* Solve */
  ierr = PFLOWSolve(pflow);CHKERRQ(ierr);

  /* Update line flows, Pgen, Qgen, and other parameters */
  ierr = PFLOWPostSolve(pflow);CHKERRQ(ierr);

  /* Trip lines 4-6 and 8-9 so that two islands are created
     Island 1: 1-4-5-7-2-8
     Island 2: 3-9-6
  */
  ierr = PetscPrintf(PETSC_COMM_SELF,"Tripping branch %d -- %d\n",8,9);CHKERRQ(ierr);
  /* Trip a line 4-6 */
   ierr = PFLOWSetLineStatus(pflow,8,9,"1 ",0);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_SELF,"Tripping branch %d -- %d\n",4,6);CHKERRQ(ierr);
  /* Trip a line 4-6 */
  ierr = PFLOWSetLineStatus(pflow,4,6,"1 ",0);CHKERRQ(ierr);

  /* Solve */
  ierr = PFLOWSolve(pflow);CHKERRQ(ierr);

  /* Update line flows, Pgen, Qgen, and other parameters */
  ierr = PFLOWPostSolve(pflow);CHKERRQ(ierr);

  /* Switch OFF generator 3, this will blackout island 2 */
  ierr = PFLOWSetGenStatus(pflow,3,"1 ",0);CHKERRQ(ierr);

  /* Solve */
  ierr = PFLOWSolve(pflow);CHKERRQ(ierr);

  /* Update line flows, Pgen, Qgen, and other parameters */
  ierr = PFLOWPostSolve(pflow);CHKERRQ(ierr);

  /* Destroy PFLOW object */
  ierr = PFLOWDestroy(&pflow);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}
  
