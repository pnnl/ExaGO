#include <iostream>
#include <cstdio>
#include <string>

#include <opflow.h>
#include <scopflow_config.h>

#include "opflow/OpflowTests.hpp"

/* Converts an array xin in natural ordering to an array xout in sparse-dense                                                                   
   ordering                                                                                                                                     
*/
void naturaltospdense(const double *xin,double *xout,int *idxn2sd_map,int nx)
{
  int i;

  for(i=0; i < nx; i++) {
    xout[idxn2sd_map[i]] = xin[i];
  }
}

/* Converts an array xin in sparse dense ordering to an array xout in natural                                                                   
   ordering                                                                                                                                     
*/
void spdensetonatural(const double *xin,double *xout,int *idxn2sd_map,int nx)
{
  int i;

  for(i=0; i < nx; i++) {
    xout[i] = xin[idxn2sd_map[i]];
  }
}

int main(int argc, char** argv)
{
  const bool     isTestOpflowModelPBPOLHIOP = true;
  const bool     isTestOpflowModelPBPOL2    = false;
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

  /* Set opflow solver type to IPOPT */
  ierr = OPFLOWSetSolver(opflow,OPFLOWSOLVER_TAO);CHKERRQ(ierr);

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

  /* Get reference equality constraint and inequality constraint Jacobian */
  Mat Jeq,Jineq;
  ierr = OPFLOWGetConstraintJacobian(opflow,&Jeq,&Jineq);CHKERRQ(ierr);

  Vec Lambda;
  Mat Hess;
  PetscScalar obj_factor;
  ierr = OPFLOWGetConstraintMultipliers(opflow,&Lambda);CHKERRQ(ierr);
  ierr = OPFLOWGetHessian(opflow,&Hess,&obj_factor);CHKERRQ(ierr);

  // Failure counter
  int fail = 0;

  if(isTestOpflowModelPBPOLHIOP)
  {
    OPFLOW opflowtest;
    exago::tests::TestOpflow test;

    std::cout << "\nTesting custom power balance model in polar coordinates for HIOP"
              << "(componentwise assembly) ... \n";

    // Create optimal power flow model
    ierr = OPFLOWCreate(PETSC_COMM_WORLD,&opflowtest);CHKERRQ(ierr);

    /* Read Network data */
    ierr = OPFLOWReadMatPowerData(opflowtest,file.c_str());CHKERRQ(ierr);

    /* Set opflow model type to power balance polar2 (component assembly) */
    ierr = OPFLOWSetModel(opflowtest,OPFLOWMODEL_PBPOLHIOP);CHKERRQ(ierr);

    /* Set solver to HIOP */
    ierr = OPFLOWSetSolver(opflowtest,OPFLOWSOLVER_HIOPNEW);CHKERRQ(ierr);

    /* Set up */
    ierr = OPFLOWSetUp(opflowtest);CHKERRQ(ierr);

    int nx,nconeq,nconineq,*idxn2sd_map;
    ierr = OPFLOWGetSizes(opflowtest,&nx,&nconeq,&nconineq);CHKERRQ(ierr);
    ierr = OPFLOWGetVariableOrdering(opflowtest,&idxn2sd_map);CHKERRQ(ierr);

    std::cout << "nx = " << nx << std::endl  << "nconeq = " << nconeq << std::endl
              << "nconineq = " << nconineq << std::endl;

    double *x_vec, *xl_vec, *xu_vec, *grad_vec, \
           *x_ref, *xl_ref, *xu_ref, *grad_ref, *g_ref, *gl_ref, *gu_ref;

    ierr = PetscMalloc1(nx,&x_ref);CHKERRQ(ierr);
    ierr = PetscMalloc1(nx,&xl_ref);CHKERRQ(ierr);
    ierr = PetscMalloc1(nx,&xu_ref);CHKERRQ(ierr);
    ierr = PetscMalloc1(nx,&grad_ref);CHKERRQ(ierr);

    // Petsc Documentation mentions "You MUST call VecRestoreArray() when you no longer need access to the array."...
    ierr = VecGetArray(X,&x_vec);CHKERRQ(ierr);
    ierr = VecGetArray(Xl,&xl_vec);CHKERRQ(ierr);
    ierr = VecGetArray(Xu,&xu_vec);CHKERRQ(ierr);
    ierr = VecGetArray(grad,&grad_vec);CHKERRQ(ierr);

    ierr = VecGetArray(G,&g_ref);CHKERRQ(ierr);
    ierr = VecGetArray(Gl,&gl_ref);CHKERRQ(ierr);
    ierr = VecGetArray(Gu,&gu_ref);CHKERRQ(ierr);

    /* Convert from natural to sparse dense ordering */
    naturaltospdense(x_vec,x_ref,idxn2sd_map,nx);CHKERRQ(ierr);
    naturaltospdense(xl_vec,xl_ref,idxn2sd_map,nx);CHKERRQ(ierr);
    naturaltospdense(xu_vec,xu_ref,idxn2sd_map,nx);CHKERRQ(ierr);
    naturaltospdense(grad_vec,grad_ref,idxn2sd_map,nx);CHKERRQ(ierr);
    /* _ref pointers are now in sparse-dense ordering */
    
    fail += test.computeVariableBounds(opflowtest,xl_ref,xu_ref);
    fail += test.computeObjective(opflowtest,x_ref,obj);
    fail += test.computeGradient(opflowtest,x_ref,grad_ref);
    fail += test.computeConstraints(opflowtest,x_ref,g_ref);
    fail += test.computeConstraintBounds(opflowtest,gl_ref,gu_ref);

    //    fail += test.computeConstraintJacobian(opflowtest,X,Jeq,Jineq);
    //    fail += test.computeHessian(opflowtest,X,Lambda,obj_factor,Hess);

    ierr = PetscFree(x_ref);CHKERRQ(ierr);
    ierr = PetscFree(xl_ref);CHKERRQ(ierr);
    ierr = PetscFree(xu_ref);CHKERRQ(ierr);
    ierr = PetscFree(grad_ref);CHKERRQ(ierr);

    ierr = OPFLOWDestroy(&opflowtest);CHKERRQ(ierr);
  }

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
    fail += test.computeConstraintJacobian(opflowtest,X,Jeq,Jineq);
    fail += test.computeHessian(opflowtest,X,Lambda,obj_factor,Hess);

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
