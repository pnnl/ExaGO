#include <iostream>
#include <cstdio>
#include <string>

#include <opflow.h>
#include <scopflow_config.h>

#include "opflow/OpflowTests.hpp"

int main(int argc, char** argv)
{
  const bool     isTestOpflowModelPBPOL2 = true;
  const bool     isTestOpflowModelPBCAR  = false;
  const bool     isTestOpflowModelIBCAR2 = false;
  const bool     isTestOpflowModelIBCAR  = false;
  int            globalFail              = 0;
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

  /*
   * This test block has more verbose comments explaining what each step is
   * doing. The latter tests will reuse much of this, so refer to this test
   * for comments.
   */
  if(isTestOpflowModelPBPOL2)
  {
    Vec                      X;
    Vec                      Xl, Xu;
    Vec                      G;
    Vec                      Gl, Gu;
    Vec                      grad;
    double                   obj;
    OPFLOW                   opflowtest;
    exago::tests::TestOpflow test;
    int                      fail = 0;

    /* Get the solution */
    ierr = OPFLOWGetSolution(opflow,&X);CHKERRQ(ierr);

    /* Test for Variable bounds */
    ierr = OPFLOWGetVariableBounds(opflow,&Xl,&Xu);CHKERRQ(ierr); /* Reference */

    /* Test for objective value */
    ierr = OPFLOWGetObjective(opflow,&obj);CHKERRQ(ierr); /* Reference */

    /* Gradient */
    ierr = VecDuplicate(X,&grad);CHKERRQ(ierr);
    ierr = OPFLOWComputeGradient(opflow,X,grad);CHKERRQ(ierr); /* Reference */
    /* Constraints */

    ierr = OPFLOWGetConstraints(opflow,&G);CHKERRQ(ierr); /* Reference */

    /* Constraint bounds */
    ierr = OPFLOWGetConstraintBounds(opflow,&Gl,&Gu);CHKERRQ(ierr); /* Reference */

    /* Create the test opflow */
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

    ierr = VecDestroy(&G);CHKERRQ(ierr);
    ierr = VecDestroy(&Gl);CHKERRQ(ierr);
    ierr = VecDestroy(&Gu);CHKERRQ(ierr);
    ierr = VecDestroy(&X);CHKERRQ(ierr);
    ierr = VecDestroy(&Xl);CHKERRQ(ierr);
    ierr = VecDestroy(&Xu);CHKERRQ(ierr);
    ierr = VecDestroy(&grad);CHKERRQ(ierr);
    ierr = OPFLOWDestroy(&opflowtest);CHKERRQ(ierr);
    globalFail += fail;
  }

  if (isTestOpflowModelPBCAR)
  {
    Vec                      X;
    Vec                      Xl, Xu;
    Vec                      G;
    Vec                      Gl, Gu;
    Vec                      grad;
    double                   obj;
    OPFLOW                   opflowtest;
    exago::tests::TestOpflow test;
    int                      fail = 0;

    /* Recreate references */
    ierr = OPFLOWGetSolution(opflow,&X);CHKERRQ(ierr);
    ierr = OPFLOWGetVariableBounds(opflow,&Xl,&Xu);CHKERRQ(ierr); 
    ierr = OPFLOWGetObjective(opflow,&obj);CHKERRQ(ierr); 
    ierr = VecDuplicate(X,&grad);CHKERRQ(ierr);
    ierr = OPFLOWComputeGradient(opflow,X,grad);CHKERRQ(ierr); 
    ierr = OPFLOWGetConstraints(opflow,&G);CHKERRQ(ierr); 
    ierr = OPFLOWGetConstraintBounds(opflow,&Gl,&Gu);CHKERRQ(ierr); 

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

    ierr = VecDestroy(&G);CHKERRQ(ierr);
    ierr = VecDestroy(&Gl);CHKERRQ(ierr);
    ierr = VecDestroy(&Gu);CHKERRQ(ierr);
    ierr = VecDestroy(&X);CHKERRQ(ierr);
    ierr = VecDestroy(&Xl);CHKERRQ(ierr);
    ierr = VecDestroy(&Xu);CHKERRQ(ierr);
    ierr = VecDestroy(&grad);CHKERRQ(ierr);
    ierr = OPFLOWDestroy(&opflowtest);CHKERRQ(ierr);
    globalFail += fail;
  }

  if (isTestOpflowModelIBCAR)
  {
    Vec                      X;
    Vec                      Xl, Xu;
    Vec                      G;
    Vec                      Gl, Gu;
    Vec                      grad;
    double                   obj;
    OPFLOW                   opflowtest;
    exago::tests::TestOpflow test;
    int                      fail = 0;

    ierr = OPFLOWGetSolution(opflow,&X);CHKERRQ(ierr);
    ierr = OPFLOWGetVariableBounds(opflow,&Xl,&Xu);CHKERRQ(ierr); 
    ierr = OPFLOWGetObjective(opflow,&obj);CHKERRQ(ierr); 
    ierr = VecDuplicate(X,&grad);CHKERRQ(ierr);
    ierr = OPFLOWComputeGradient(opflow,X,grad);CHKERRQ(ierr); 
    ierr = OPFLOWGetConstraints(opflow,&G);CHKERRQ(ierr); 
    ierr = OPFLOWGetConstraintBounds(opflow,&Gl,&Gu);CHKERRQ(ierr); 

    ierr = OPFLOWCreate(PETSC_COMM_WORLD,&opflowtest);CHKERRQ(ierr);
    ierr = OPFLOWReadMatPowerData(opflowtest,file.c_str());CHKERRQ(ierr);
    ierr = OPFLOWSetModel(opflowtest,OPFLOWMODEL_IBCAR);CHKERRQ(ierr);
    ierr = OPFLOWSetUp(opflowtest);CHKERRQ(ierr);

    fail += test.computeVariableBounds(opflowtest,Xl,Xu);
    fail += test.computeObjective(opflowtest,X,obj);
    fail += test.computeGradient(opflowtest,X,grad);
    fail += test.computeConstraints(opflowtest,X,G);
    fail += test.computeConstraintBounds(opflowtest,Gl,Gu);

    ierr = VecDestroy(&G);CHKERRQ(ierr);
    ierr = VecDestroy(&Gl);CHKERRQ(ierr);
    ierr = VecDestroy(&Gu);CHKERRQ(ierr);
    ierr = VecDestroy(&X);CHKERRQ(ierr);
    ierr = VecDestroy(&Xl);CHKERRQ(ierr);
    ierr = VecDestroy(&Xu);CHKERRQ(ierr);
    ierr = VecDestroy(&grad);CHKERRQ(ierr);
    ierr = OPFLOWDestroy(&opflowtest);CHKERRQ(ierr);
    globalFail += fail;
  }

  if (isTestOpflowModelIBCAR2)
  {
    Vec                      X;
    Vec                      Xl, Xu;
    Vec                      G;
    Vec                      Gl, Gu;
    Vec                      grad;
    double                   obj;
    OPFLOW                   opflowtest;
    exago::tests::TestOpflow test;
    int                      fail = 0;

    ierr = OPFLOWGetSolution(opflow,&X);CHKERRQ(ierr);
    ierr = OPFLOWGetVariableBounds(opflow,&Xl,&Xu);CHKERRQ(ierr); 
    ierr = OPFLOWGetObjective(opflow,&obj);CHKERRQ(ierr); 
    ierr = VecDuplicate(X,&grad);CHKERRQ(ierr);
    ierr = OPFLOWComputeGradient(opflow,X,grad);CHKERRQ(ierr); 
    ierr = OPFLOWGetConstraints(opflow,&G);CHKERRQ(ierr); 
    ierr = OPFLOWGetConstraintBounds(opflow,&Gl,&Gu);CHKERRQ(ierr); 

    ierr = OPFLOWCreate(PETSC_COMM_WORLD,&opflowtest);CHKERRQ(ierr);
    ierr = OPFLOWReadMatPowerData(opflowtest,file.c_str());CHKERRQ(ierr);
    ierr = OPFLOWSetModel(opflowtest,OPFLOWMODEL_IBCAR2);CHKERRQ(ierr);
    ierr = OPFLOWSetUp(opflowtest);CHKERRQ(ierr);

    fail += test.computeVariableBounds(opflowtest,Xl,Xu);
    fail += test.computeObjective(opflowtest,X,obj);
    fail += test.computeGradient(opflowtest,X,grad);
    fail += test.computeConstraints(opflowtest,X,G);
    fail += test.computeConstraintBounds(opflowtest,Gl,Gu);

    ierr = VecDestroy(&G);CHKERRQ(ierr);
    ierr = VecDestroy(&Gl);CHKERRQ(ierr);
    ierr = VecDestroy(&Gu);CHKERRQ(ierr);
    ierr = VecDestroy(&X);CHKERRQ(ierr);
    ierr = VecDestroy(&Xl);CHKERRQ(ierr);
    ierr = VecDestroy(&Xu);CHKERRQ(ierr);
    ierr = VecDestroy(&grad);CHKERRQ(ierr);
    ierr = OPFLOWDestroy(&opflowtest);CHKERRQ(ierr);
    globalFail += fail;
  }

  // Why does this cause a SEGV?
  // ierr = OPFLOWDestroy(&opflow);CHKERRQ(ierr);
  PetscFinalize();

  return globalFail;    
}
