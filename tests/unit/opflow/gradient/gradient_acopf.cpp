#include <iostream>
#include <cstdio>
#include <string>

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
  Vec X, grad;
  int nx, nconeq, nconineq;
  PetscScalar *x_arr, *grad_arr;
  int fail = 0;
  char file_c_str[PETSC_MAX_PATH_LEN];
  std::string file;
  char appname[] = "opflow";
  MPI_Comm comm = MPI_COMM_WORLD;

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
  ierr = OPFLOWGetSolution(opflowtest, &X);
  CHKERRQ(ierr);

  ierr = OPFLOWGetSizes(opflowtest, &nx, &nconeq, &nconineq);
  CHKERRQ(ierr);

  // Set Grad Based on X
  ierr = VecDuplicate(X, &grad);
  CHKERRQ(ierr);

  ierr = VecGetArray(X, &x_arr);
  CHKERRQ(ierr);
  ierr = VecGetArray(grad, &grad_arr);
  CHKERRQ(ierr);

  // Gradient for all generator costs should be 1
  // 10 corresponds to locations of non-zero values
  for (int i = 0; i < nx; i++) {
    if (x_arr[i] == 10) {
      grad_arr[i] = 1.0;
    }
  }

  ierr = VecRestoreArray(X, &x_arr);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(grad, &grad_arr);
  CHKERRQ(ierr);

  fail += test.computeGradient(opflowtest, X, grad);

  ierr = VecDestroy(&X);
  CHKERRQ(ierr);
  ierr = VecDestroy(&grad);
  CHKERRQ(ierr);

  ierr = OPFLOWDestroy(&opflowtest);
  CHKERRQ(ierr);

  ExaGOFinalize();
  return fail;
}
