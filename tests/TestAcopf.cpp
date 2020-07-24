#include <iostream>
#include <cstdio>
#include <string>

#include <opflow.h>
#include <scopflow_config.h>

#include "opflow/OpflowTests.hpp"

int main(int argc, char** argv)
{
  const bool     isTestOpflowModelPBPOL2    = true;
  const bool     isTestOpflowModelPBCAR     = false;
  const bool     isTestOpflowModelIBCAR2    = false;
  const bool     isTestOpflowModelIBCAR     = false;
  PetscErrorCode ierr;
  PetscBool      flg;
  OPFLOW         opflow;
  char           file_c_str[PETSC_MAX_PATH_LEN];
  std::string file;
  std::string options_pathname(EXAGO_OPTIONS_DIR);
  std::string filename("/opflowoptions");
  std::cout << options_pathname << "\n";
  std::cout << filename << "\n";
  options_pathname += filename;
  std::cout << options_pathname << "\n";

  char help[] = "Unit tests for ACOPF\n";
  PetscInitialize(&argc, &argv, options_pathname.c_str(), help);

  /* We create two OPFLOWs in this test. The first one is what we
     use as a reference. This will have the OPFLOW model type = OPFLOWMODEL_PBPOL, 
     i.e., the non-component-assembly version of the model.
     We will use the solution from this opflow in the computation of the
     various kernels of the test opflow (opflowtest). opflowtest will use the
     model type=OPFLOWMODEL
  */

  // Create reference data for testing 

  /* Create reference opflow */
  ierr = OPFLOWCreate(PETSC_COMM_WORLD,&opflow);CHKERRQ(ierr);

  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL,NULL,"-netfile",file_c_str,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);

  if(!flg){
    file = "../datafiles/case9mod.m";
  } 
  else{
    file = file_c_str;
  }
  
  /* Read Network data */
  ierr = OPFLOWReadMatPowerData(opflow,file.c_str());CHKERRQ(ierr);

  /* Set opflow model type to power balance polar (no component assembly) */
  ierr = OPFLOWSetModel(opflow,OPFLOWMODEL_PBPOL);CHKERRQ(ierr);

  /* Solve OPFLOW to get the reference solution */
  ierr = OPFLOWSolve(opflow);

  /* Get the reference solution */
  Vec X;
  ierr = OPFLOWGetSolution(opflow,&X);CHKERRQ(ierr);

  /* Get reference variable bounds */
  Vec Xl, Xu;
  ierr = OPFLOWGetVariableBounds(opflow,&Xl,&Xu);CHKERRQ(ierr);

  /* Get reference objective value */
  double obj;
  ierr = OPFLOWGetObjective(opflow,&obj);CHKERRQ(ierr);

  /* Get reference value for the objective radient */
  Vec grad;
  ierr = VecDuplicate(X,&grad);CHKERRQ(ierr);
  ierr = OPFLOWComputeGradient(opflow,X,grad);CHKERRQ(ierr);

  /* Get reference constraint residuals */
  Vec G;
  ierr = OPFLOWGetConstraints(opflow,&G);CHKERRQ(ierr);

  /* Get reference constraint bounds */
  Vec Gl, Gu;
  ierr = OPFLOWGetConstraintBounds(opflow,&Gl,&Gu);CHKERRQ(ierr);

  // Failure counter
  int fail = 0;

  if(isTestOpflowModelPBPOL2)
  {
    OPFLOW opflowtest;
    exago::tests::TestOpflow test;

    std::cout << "\nTesting power balance model in polar coordinates "
              << "(componentwise assembly) ... \n";

    // Create optimal power flow model
    ierr = OPFLOWCreate(PETSC_COMM_WORLD,&opflowtest);CHKERRQ(ierr);

    /* Read Network data */
    ierr = OPFLOWReadMatPowerData(opflowtest,file.c_str());CHKERRQ(ierr);

    /* Set opflow model type to power balance polar2 (component assembly) */
    ierr = OPFLOWSetModel(opflowtest,OPFLOWMODEL_PBPOL2);CHKERRQ(ierr);

    /* Set up */
    ierr = OPFLOWSetUp(opflowtest);CHKERRQ(ierr);

    fail += test.computeVariableBounds(opflowtest,Xl,Xu);
    fail += test.computeObjective(opflowtest,X,obj);
    fail += test.computeGradient(opflowtest,X,grad);
    fail += test.computeConstraints(opflowtest,X,G);
    fail += test.computeConstraintBounds(opflowtest,Gl,Gu);

    ierr = OPFLOWDestroy(&opflowtest);CHKERRQ(ierr);
  }

  if (isTestOpflowModelPBCAR)
  {
    OPFLOW                   opflowtest;
    exago::tests::TestOpflow test;

    std::cout << "\nTesting power balance model in cartesian coordinates "
              << "(nodewise assembly) ... \n";

    /* Set up test opflow */
    ierr = OPFLOWCreate(PETSC_COMM_WORLD,&opflowtest);CHKERRQ(ierr);
    ierr = OPFLOWReadMatPowerData(opflowtest,file.c_str());CHKERRQ(ierr);
    ierr = OPFLOWSetModel(opflowtest,OPFLOWMODEL_PBCAR);CHKERRQ(ierr);
    ierr = OPFLOWSetUp(opflowtest);CHKERRQ(ierr);

    fail += test.computeVariableBounds(opflowtest,Xl,Xu);
    fail += test.computeObjective(opflowtest,X,obj);
    fail += test.computeGradient(opflowtest,X,grad);
    fail += test.computeConstraints(opflowtest,X,G);
    fail += test.computeConstraintBounds(opflowtest,Gl,Gu);

    ierr = OPFLOWDestroy(&opflowtest);CHKERRQ(ierr);
  }

  if (isTestOpflowModelIBCAR)
  {
    OPFLOW opflowtest;
    exago::tests::TestOpflow test;

    std::cout << "\nTesting current balance model in cartesian coordinates "
              << "(nodewise assembly) ... \n";

    ierr = OPFLOWCreate(PETSC_COMM_WORLD,&opflowtest);CHKERRQ(ierr);
    ierr = OPFLOWReadMatPowerData(opflowtest,file.c_str());CHKERRQ(ierr);
    ierr = OPFLOWSetModel(opflowtest,OPFLOWMODEL_IBCAR);CHKERRQ(ierr);
    ierr = OPFLOWSetUp(opflowtest);CHKERRQ(ierr);

    fail += test.computeVariableBounds(opflowtest,Xl,Xu);
    fail += test.computeObjective(opflowtest,X,obj);
    fail += test.computeGradient(opflowtest,X,grad);
    fail += test.computeConstraints(opflowtest,X,G);
    fail += test.computeConstraintBounds(opflowtest,Gl,Gu);

    ierr = OPFLOWDestroy(&opflowtest);CHKERRQ(ierr);
  }

  if (isTestOpflowModelIBCAR2)
  {
    OPFLOW opflowtest;
    exago::tests::TestOpflow test;

    std::cout << "\nTesting current balance model in cartesian coordinates "
              << "(?) ... \n";

    ierr = OPFLOWCreate(PETSC_COMM_WORLD,&opflowtest);CHKERRQ(ierr);
    ierr = OPFLOWReadMatPowerData(opflowtest,file.c_str());CHKERRQ(ierr);
    ierr = OPFLOWSetModel(opflowtest,OPFLOWMODEL_IBCAR2);CHKERRQ(ierr);
    ierr = OPFLOWSetUp(opflowtest);CHKERRQ(ierr);

    fail += test.computeVariableBounds(opflowtest,Xl,Xu);
    fail += test.computeObjective(opflowtest,X,obj);
    fail += test.computeGradient(opflowtest,X,grad);
    fail += test.computeConstraints(opflowtest,X,G);
    fail += test.computeConstraintBounds(opflowtest,Gl,Gu);

    ierr = OPFLOWDestroy(&opflowtest);CHKERRQ(ierr);
  }

  /* Destroy OPFLOW objects */
  ierr = OPFLOWDestroy(&opflow);CHKERRQ(ierr);
  ierr = VecDestroy(&grad);CHKERRQ(ierr);
  PetscFinalize();

  return fail;
}
