#include <iostream>
#include <cstdio>
#include <string>

#include<unistd.h>

#include <private/opflowimpl.h>
#include <exago_config.h>
#include <utils.h>

#include "opflow_tests.h"
#include "test_acopf_utils.h"

/**
 * @brief Unit test driver for gradient function
 * @see /tests/unit/opflow/opflow_tests.h for kernel tested by this driver
 *
 * You can pass the following option to the gradientAcopf executatable through
 * the command line (implemented using PETSc options):
 *
 *    ~ -netfile <data_file> : Specifies the input data file to test against.
 * Default value is `/<exago_dir>/datafiles/case9/case9mod.m`. See directory
 * datafiles for other potential inputs.
 *
 *
 */
int main(int argc, char **argv) {
  PetscErrorCode ierr;
  PetscBool flg;
  PetscBool flg2;
  Vec X, grad;
  Vec Xl;
  Vec Xu;
  Vec computeXl;
  Vec computeXu;
  IS x_is;
  int nx, nconeq, nconineq;
  PetscScalar *x_arr, *grad_arr;
  int fail = 0;
  char file_c_str[PETSC_MAX_PATH_LEN];
  std::string file;
  char appname[] = "opflow";
  MPI_Comm comm = MPI_COMM_WORLD;

  char help[] = "Unit tests for setting variable bounds running opflow\n";

  /** Use `ExaGOLogSetLoggingFileName("opflow-logfile");` to log the output. */
  //ExaGOLogSetLoggingFileName("opflow-logfile");
  ierr = ExaGOInitialize(comm, &argc, &argv, appname, help);
  if (ierr) {
    fprintf(stderr, "Could not initialize ExaGO application %s.\n", appname);
    return ierr;
  }

  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL, NULL, "-netfile", file_c_str,
                               PETSC_MAX_PATH_LEN, &flg);
  CHKERRQ(ierr);

  if (!flg) {
    file = "../datafiles/case9/case9mod.m";
  } else {
    file.assign(file_c_str);
  }

  OPFLOW opflowtest;
  exago::tests::TestOpflow test;

  /* Set up test opflow */
  ierr = OPFLOWCreate(PETSC_COMM_WORLD, &opflowtest);
  CHKERRQ(ierr);
  ierr = OPFLOWReadMatPowerData(opflowtest, file.c_str());
  CHKERRQ(ierr);
  ierr = OPFLOWSetInitializationType(opflowtest, OPFLOWINIT_FROMFILE);
  CHKERRQ(ierr);
  ierr = OPFLOWSetUp(opflowtest);
  CHKERRQ(ierr);
  
  
  ierr = OPFLOWGetVariableBounds(opflowtest, &Xl, &Xu);
  CHKERRQ(ierr);
  printf("xlxu\n");
  VecView(Xl,PETSC_VIEWER_STDOUT_SELF);
  VecView(Xu, PETSC_VIEWER_STDOUT_SELF);

  ierr = OPFLOWGetVariableBounds(opflowtest, &computeXl, &computeXu);
  CHKERRQ(ierr);
  printf("compute\n");
  VecView(computeXl, PETSC_VIEWER_STDOUT_SELF);
  VecView(computeXu, PETSC_VIEWER_STDOUT_SELF);

  ierr = OPFLOWComputeVariableBounds(opflowtest, computeXl, computeXu);
  CHKERRQ(ierr);
  printf("computexlxu\n");
  VecView(computeXl, PETSC_VIEWER_STDOUT_SELF);
  VecView(computeXu, PETSC_VIEWER_STDOUT_SELF);

/*  volatile int i = 0;
  char hostname[256];
  gethostname(hostname, sizeof(hostname));
  printf("PID %d on %s ready for attach\n", getpid(), hostname);
  fflush(stdout);
  while (0 == i)
      sleep(5);
*/
  
/*    
  VecView(Xl,PETSC_VIEWER_STDOUT_SELF);
  VecView(Xu, PETSC_VIEWER_STDOUT_SELF);
  VecView(computeXl, PETSC_VIEWER_STDOUT_SELF);
  VecView(computeXu, PETSC_VIEWER_STDOUT_SELF); */
  ierr = VecEqual(Xl,computeXl, &flg);
  CHKERRQ(ierr);
  printf("flag1 : %d \n", flg);
  ierr = VecEqual(Xu, computeXu, &flg2);
  CHKERRQ(ierr);
  printf("flag2 : %d \n", flg2);
  
  ierr = VecWhichEqual(Xl, computeXl, &x_is);
  CHKERRQ(ierr);
  ISView(x_is, PETSC_VIEWER_STDOUT_SELF);
  ierr = VecWhichEqual(Xu, computeXu, &x_is);
  CHKERRQ(ierr);
  ISView(x_is, PETSC_VIEWER_STDOUT_SELF);

  if(!(flg && flg2))
  {
    return fail + 1;
  }

  ierr = OPFLOWGetSolution(opflowtest, &X);
  CHKERRQ(ierr);

  printf("solution: \n");
  VecView(X, PETSC_VIEWER_STDOUT_SELF);
  
  ierr = OPFLOWGetSizes(opflowtest, &nx, &nconeq, &nconineq);
  CHKERRQ(ierr);

  printf("sizes got: \n");
  // Set Grad Based on X
  ierr = VecDuplicate(X, &grad);
  CHKERRQ(ierr);

  printf("vec dup: \n");
  ierr = VecGetArray(X, &x_arr);
  CHKERRQ(ierr);
  ierr = VecGetArray(grad, &grad_arr);
  CHKERRQ(ierr);

  printf("vec get: \n");
  // Gradient for all generator costs should be 1
  // 10 corresponds to locations of non-zero values
  for (int i = 0; i < nx; i++) {
    if (x_arr[i] == 10) {
      grad_arr[i] = 1.0;
    }
  }

  printf("for loop: \n");
  ierr = VecRestoreArray(X, &x_arr);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(grad, &grad_arr);
  CHKERRQ(ierr);

  printf("vec restore \n");
  
  fail += test.computeGradient(opflowtest, X, grad);

  printf("compute grad \n");
  printf("gradient:  \n");
  
  printf("fail: %d \n", fail);
  
  
  ierr = VecDestroy(&grad);
  CHKERRQ(ierr);

  printf("vec destroy \n");
  ierr = OPFLOWDestroy(&opflowtest);
  CHKERRQ(ierr);

  printf("opflow destroy \n");
  ExaGOFinalize();
  printf("exago final \n");
  return fail;
}

