#include <iostream>
#include <cstdio>
#include <string>

#include <private/opflowimpl.h>
#include <exago_config.h>
#include <utils.h>

#include "opflow_tests.h"
#include "test_acopf_utils.h"

// Debugging
#include <unistd.h>

/**
 * @brief Unit test driver for gradient function
 * @see opflow/opflow_tests.h for kernel tested by this driver
 *
 * You can pass two options to the gradient_acopf executatable through the
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
  std::string file;
  char appname[] = "opflow";
  MPI_Comm comm = MPI_COMM_WORLD;
  int num_copies = 1;

  // Debugging
  /*
  volatile int i = 0;
  char hostname[256];
  gethostname(hostname, sizeof(hostname));
  printf("PID %d on %s ready for attach\n", getpid(), hostname);
  fflush(stdout);
  while (0 == i)
    sleep(5);
    */

  char help[] = "Unit tests for gradient function running opflow\n";

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

  /* Get network data file from command line */
  ierr = PetscOptionsGetInt(NULL, NULL, "-num_copies", &num_copies, &flg);
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
  ierr = OPFLOWSetModel(opflowtest, OPFLOWMODEL_PBPOL);
  CHKERRQ(ierr);
  ierr = OPFLOWSetInitializationType(opflowtest, OPFLOWINIT_FROMFILE);
  CHKERRQ(ierr);
  ierr = OPFLOWSetUp(opflowtest);
  CHKERRQ(ierr);

  ierr = OPFLOWGetSolution(opflowtest, &X);
  CHKERRQ(ierr);

  // Set the grad array to be the solution
  int nx, nconeq, nconineq;
  ierr = OPFLOWGetSizes(opflowtest, &nx, &nconeq, &nconineq);
  CHKERRQ(ierr);

  ierr = VecDuplicate(X, &grad);
  CHKERRQ(ierr);

  ierr = VecAssemblyBegin(grad);
  CHKERRQ(ierr);

  // Start at 2 as first two elements are always 0
  for(int i=2; i<nx; i++)
  {
    if((i - 2) % 12 == 4)
    {
      ierr = VecSetValue(grad, i, 4.5, INSERT_VALUES);
      CHKERRQ(ierr);
    }
    else if ((i - 2) % 12 == 8 || (i - 2) % 12 == 9)
    {
      ierr = VecSetValue(grad, i, 100000, INSERT_VALUES);
      CHKERRQ(ierr);
    }
  }

  ierr = VecAssemblyEnd(grad);
  CHKERRQ(ierr);

  // Call test.computeGradient

  ierr = VecView(X, 0);
  CHKERRQ(ierr);
  ierr = VecView(grad, 0);
  CHKERRQ(ierr);

  fail += test.computeGradient(opflowtest, X, grad);

  /*
  ierr = VecRestoreArray(X, &x_ref);
  CHKERRQ(ierr);
  */

  // This vector is destroyed by opflow internally...
  // ierr = VecDestroy(&X);
  // CHKERRQ(ierr);

  ierr = VecDestroy(&grad);
  CHKERRQ(ierr);

  ierr = OPFLOWDestroy(&opflowtest);
  CHKERRQ(ierr);

  ExaGOFinalize();
  return fail;
}
