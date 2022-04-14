#include <iostream>
#include <cstdio>
#include <string>

#include <private/opflowimpl.h>
#include <exago_config.h>
#include <utils.h>

#include "opflow_tests.h"
#include "test_acopf_utils.h"

#include <unistd.h>

/**
 * @brief Unit test driver for objective function
 * @see opflow/OpflowTests.hpp for kernel tested by this driver
 *
 * You can pass two options to the objectiveAcopf executatable through the
 * command line (implemented using PETSc options):
 *
 *    ~ -netfile <data_file> : Specifies the input data file to test against.
 * Default value is `/<exago_dir>/datafiles/case9/case9mod.m`. See directory
 * datafiles for other potential inputs.
 *
 *    ~ -num_copies <number> : Specifies the number of replications of the
 * network given through `-netfile`. If this is not set properly, test may fail
 *
 */
int main(int argc, char **argv) {
  PetscErrorCode ierr;
  PetscBool flg;
  Vec X, grad;
  int fail = 0;
  char file_c_str[PETSC_MAX_PATH_LEN];
  // char validation_c_str[PETSC_MAX_PATH_LEN];
  std::string file;
  char appname[] = "opflow";
  MPI_Comm comm = MPI_COMM_WORLD;
  // int num_copies = 0;

  char help[] = "Unit tests for gradient function running opflow\n";

  // volatile int i = 0;
  // char hostname[256];
  // gethostname(hostname, sizeof(hostname));
  // printf("PID %d on %s ready for attach\n", getpid(), hostname);
  // fflush(stdout);
  // while (0 == i)
  //   sleep(5);

  /** Use `ExaGOLogSetLoggingFileName("opflow-logfile");` to log the output. */
  ierr = ExaGOInitialize(comm, &argc, &argv, appname, help);
  if (ierr) {
    fprintf(stderr, "Could not initialize ExaGO application %s.\n", appname);
    return ierr;
  }

  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL, NULL, "-netfile", file_c_str,
                               PETSC_MAX_PATH_LEN, &flg);
  CHKERRQ(ierr);

  /* Get num_copies from command line */
  //ierr = PetscOptionsGetInt(NULL, NULL, "-num_copies", &num_copies, &flg);
  //CHKERRQ(ierr);

  if (!flg) {
    file = "../datafiles/case9/case9mod.m";
  } else {
    file.assign(file_c_str);
  }

  // Set obj_value as reference solution, and run as usual
  /* Get validation data file from command line */
  // ierr = PetscOptionsGetString(NULL, NULL, "-validation", validation_c_str,
  //                              PETSC_MAX_PATH_LEN, &flg);
  // CHKERRQ(ierr);

  OPFLOW opflowtest;
  exago::tests::TestOpflow test;

  /* Set up test opflow */
  ierr = OPFLOWCreate(PETSC_COMM_WORLD, &opflowtest);
  CHKERRQ(ierr);

  ierr = OPFLOWReadMatPowerData(opflowtest, "/ccsopen/home/rcrutherford/exago/exago-git/datafiles/unit/opflow/gradient/OFG_unittest1.m");
  CHKERRQ(ierr);

  ierr = OPFLOWSetInitializationType(opflowtest, OPFLOWINIT_FROMFILE);
  CHKERRQ(ierr);

  ierr = OPFLOWSetUp(opflowtest);
  CHKERRQ(ierr);

  ierr = OPFLOWGetSolution(opflowtest, &X);
  CHKERRQ(ierr);

  VecView(X, 	PETSC_VIEWER_STDOUT_SELF);

  // Set X From configuration
  int nx, nconeq, nconineq;
  ierr = OPFLOWGetSizes(opflowtest, &nx, &nconeq, &nconineq);
  CHKERRQ(ierr);

  for(int i= 0; i< nx; i++)
  {
    VecSetValue(X, i, 10, INSERT_VALUES);
  }

  ierr = VecDuplicate(X, &grad);
  CHKERRQ(ierr);
  // readFromFile(&grad, validation_c_str);

  // fail += test.computeGradient(opflowtest, X, grad);
  Vec grad_tmp;
  
  ierr = VecDuplicate(grad, &grad_tmp);
  CHKERRQ(ierr);
  ierr = VecSet(grad_tmp, 0.0);
  CHKERRQ(ierr);
  ierr = OPFLOWComputeGradient(opflowtest, X, grad_tmp);
  CHKERRQ(ierr);

  //fail += verifyAnswer(grad_tmp, grad);
  VecView(grad, PETSC_VIEWER_STDOUT_SELF);
  VecView(grad_tmp, PETSC_VIEWER_STDOUT_SELF);

  ierr = VecDestroy(&grad_tmp);
  CHKERRQ(ierr); 

  ierr = VecDestroy(&grad);
  CHKERRQ(ierr);

  ierr = OPFLOWDestroy(&opflowtest);
  CHKERRQ(ierr);

  ExaGOFinalize();
  return fail;
}
