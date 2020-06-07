#include <iostream>
#include <string>

#include <opflow.h>
#include <scopflow_config.h>

#include "opflow/OpflowTests.hpp"

int main(int argc, char** argv)
{
  PetscErrorCode     ierr;
  OPFLOW             opflow,opflowtest;
  char               file[PETSC_MAX_PATH_LEN];
  PetscBool          flg;
  Vec                X;
  char options_pathname[200] = SCOPFLOW_OPTIONS_DIR;
  char* filename = "/opflowoptions";
  printf("%s\n", options_pathname);
  printf("%s\n", filename);
  strcat(options_pathname, filename);
  printf("%s\n", options_pathname);

  char help[] = "Unit tests for ACOPF\n";
  PetscInitialize(&argc,&argv,options_pathname,help);

  /* We create two OPFLOWs in this test. The first one is what we
     use as a reference. This will have the OPFLOW model type = OPFLOWMODEL_PBPOL, 
     i.e., the non-component-assembly version of the model.
     We will use the solution from this opflow in the computation of the
     various kernels of the test opflow (opflowtest). opflowtest will use the
     model type=OPFLOWMODEL
  */

  /* Create reference opflow */
  ierr = OPFLOWCreate(PETSC_COMM_WORLD,&opflow);CHKERRQ(ierr);

  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL,NULL,"-netfile",file,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);

  if(!flg){
    strcpy(file,"../datafiles/case9mod.m");
  }
  /* Read Network data */
  ierr = OPFLOWReadMatPowerData(opflow,file);CHKERRQ(ierr);

  /* Set opflow model type to power balance polar (no component assembly) */
  ierr = OPFLOWSetModel(opflow,OPFLOWMODEL_PBPOL);CHKERRQ(ierr);

  /* Solve OPFLOW to get the reference solution */
  ierr = OPFLOWSolve(opflow);

  /* Get the solution */
  ierr = OPFLOWGetSolution(opflow,&X);CHKERRQ(ierr);

  /* Create the test opflow */
  ierr = OPFLOWCreate(PETSC_COMM_WORLD,&opflowtest);CHKERRQ(ierr);

  /* Read Network data */
  ierr = OPFLOWReadMatPowerData(opflowtest,file);CHKERRQ(ierr);

  /* Set opflow model type to power balance polar2 (component assembly) */
  ierr = OPFLOWSetModel(opflowtest,OPFLOWMODEL_PBPOL2);CHKERRQ(ierr);

  /* Set up */
  ierr = OPFLOWSetUp(opflowtest);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_SELF,"\n\nStarting OPFLOW tests\n");CHKERRQ(ierr);
  int fail = 0;
  /* Test OPFLOW methods */

  /* Test for Variable bounds */
  ierr = PetscPrintf(PETSC_COMM_SELF,"%-35s","Variable bounds test");CHKERRQ(ierr);
  Vec Xl,Xu;
  ierr = OPFLOWGetVariableBounds(opflow,&Xl,&Xu);CHKERRQ(ierr); /* Reference */
  fail = exago::tests::TestOPFLOWComputeVariableBounds(opflowtest,Xl,Xu);
  ierr = PetscPrintf(PETSC_COMM_SELF,"%s\n",fail?"Fail":"Pass");CHKERRQ(ierr);

  /* Test for objective value */
  ierr = PetscPrintf(PETSC_COMM_SELF,"%-35s","Objective function test");CHKERRQ(ierr);
  double obj;
  ierr = OPFLOWGetObjective(opflow,&obj);CHKERRQ(ierr); /* Reference */
  fail = exago::tests::TestOPFLOWComputeObjective(opflowtest,X,obj);
  ierr = PetscPrintf(PETSC_COMM_SELF,"%s\n",fail?"Fail":"Pass");CHKERRQ(ierr);

  /* Gradient */
  ierr = PetscPrintf(PETSC_COMM_SELF,"%-35s","Gradient test");CHKERRQ(ierr);
  Vec grad;
  ierr = VecDuplicate(X,&grad);CHKERRQ(ierr);
  ierr = OPFLOWComputeGradient(opflow,X,grad);CHKERRQ(ierr); /* Reference */
  fail = exago::tests::TestOPFLOWComputeGradient(opflowtest,X,grad);
  ierr = VecDestroy(&grad);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF,"%s\n",fail?"Fail":"Pass");CHKERRQ(ierr);

  /* Constraints */
  ierr = PetscPrintf(PETSC_COMM_SELF,"%-35s","Constraints test");CHKERRQ(ierr);
  Vec G;
  ierr = OPFLOWGetConstraints(opflow,&G);CHKERRQ(ierr); /* Reference */
  fail = exago::tests::TestOPFLOWComputeConstraints(opflowtest,X,G);
  ierr = PetscPrintf(PETSC_COMM_SELF,"%s\n",fail?"Fail":"Pass");CHKERRQ(ierr);

  /* Constraint bounds */
  ierr = PetscPrintf(PETSC_COMM_SELF,"%-35s","Constraint bounds test");CHKERRQ(ierr);
  Vec Gl,Gu;
  ierr = OPFLOWGetConstraintBounds(opflow,&Gl,&Gu);CHKERRQ(ierr); /* Reference */
  fail = exago::tests::TestOPFLOWComputeConstraintBounds(opflowtest,Gl,Gu);
  ierr = PetscPrintf(PETSC_COMM_SELF,"%s\n",fail?"Fail":"Pass");CHKERRQ(ierr);

  /* Destroy OPFLOW objects */
  ierr = OPFLOWDestroy(&opflow);CHKERRQ(ierr);
  ierr = OPFLOWDestroy(&opflowtest);CHKERRQ(ierr);

  PetscFinalize();

  return 0;    
}
