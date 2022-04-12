#include <iostream>
#include <cstdio>
#include <string>

#include <private/opflowimpl.h>
#include <exago_config.h>
#include <utils.h>

#include "opflow_tests.h"
#include "test_acopf_utils.h"

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
  char validation_c_str[PETSC_MAX_PATH_LEN];
  std::string file;
  char appname[] = "opflow";
  MPI_Comm comm = MPI_COMM_WORLD;
  int num_copies = 0;

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

  /* Get num_copies from command line */
  ierr = PetscOptionsGetInt(NULL, NULL, "-num_copies", &num_copies, &flg);
  CHKERRQ(ierr);

  if (!flg) {
    file = "../datafiles/case9/case9mod.m";
  } else {
    file.assign(file_c_str);
  }

  // Set obj_value as reference solution, and run as usual
  /* Get validation data file from command line */
  ierr = PetscOptionsGetString(NULL, NULL, "-validation", validation_c_str,
                               PETSC_MAX_PATH_LEN, &flg);
  CHKERRQ(ierr);

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
  ierr = VecDuplicate(X, &grad);
  CHKERRQ(ierr);
  readFromFile(&grad, validation_c_str);

  // If we are using HIOP, need to convert X
  // The string lengths must be 65
  char modelname[64];
  char solvername[64];
  ierr = OPFLOWGetModel(opflowtest, modelname);
  ierr = OPFLOWGetSolver(opflowtest, solvername);

  if (strncmp(solvername, "HIOP", 64) == 0) {
    double *x_vec, *x_ref, *grad_vec, *grad_ref;
    ierr = VecGetArray(X, &x_vec);
    CHKERRQ(ierr);

    ierr = VecGetArray(grad, &grad_vec);
    CHKERRQ(ierr);

    int nx, nconeq, nconineq;
    ierr = OPFLOWGetSizes(opflowtest, &nx, &nconeq, &nconineq);
    CHKERRQ(ierr);

    int *idxn2sd_map;
    ierr = OPFLOWGetVariableOrdering(opflowtest, &idxn2sd_map);
    CHKERRQ(ierr);

    ierr = PetscMalloc1(nx, &x_ref);
    CHKERRQ(ierr);
    ierr = PetscMalloc1(nx, &grad_ref);
    CHKERRQ(ierr);

    /* Convert from natural to sparse dense ordering */
    naturaltospdense(x_vec, x_ref, idxn2sd_map, nx);
    CHKERRQ(ierr);
    naturaltospdense(grad_vec, grad_ref, idxn2sd_map, nx);
    CHKERRQ(ierr);
    for (int i = 0; i < nx; i++) {
      std::cout << "grad_ref[" << i << "] == " << grad_ref[i] << std::endl;
    }

    /* _ref pointers are now in sparse-dense ordering */

    // If we are running using the CPU model, nothing needs to be done
    if (strncmp(modelname, "POWER_BALANCE_HIOP", 64) == 0) {
      fail += test.computeGradient(opflowtest, x_ref, grad_ref);
    } else // Using model PBPOLRAJAHIOP
    {

      // Get resource manager instance
      auto &resmgr = umpire::ResourceManager::getInstance();

      // Get Allocator
      umpire::Allocator h_allocator = resmgr.getAllocator("HOST");

      // Register array xref with umpire
      umpire::util::AllocationRecord record_x{
          x_ref, sizeof(double) * nx, h_allocator.getAllocationStrategy()};
      resmgr.registerAllocation(x_ref, record_x);
      // Allocate and copy xref and lambdaref to device
      double *x_ref_dev;

#ifdef EXAGO_ENABLE_GPU

      ierr = OPFLOWSetHIOPComputeMode(opflowtest, "GPU");
      CHKERRQ(ierr);

      umpire::Allocator d_allocator = resmgr.getAllocator("DEVICE");
      x_ref_dev =
          static_cast<double *>(d_allocator.allocate(nx * sizeof(double)));
#else
      ierr = OPFLOWSetHIOPComputeMode(opflowtest, "CPU");
      CHKERRQ(ierr);
      x_ref_dev = x_ref;
#endif
      resmgr.copy(x_ref_dev, x_ref);

      // grad_ref is copied onto device by the test
      fail += test.computeGradient(opflowtest, x_ref_dev, grad_ref, resmgr);

#ifdef EXAGO_ENABLE_GPU
      d_allocator.deallocate(x_ref_dev);
#endif
    }

    ierr = VecRestoreArray(X, &x_ref);
    CHKERRQ(ierr);

    ierr = VecRestoreArray(grad, &grad_ref);
    CHKERRQ(ierr);

    ierr = PetscFree(x_ref);
    CHKERRQ(ierr);
    ierr = PetscFree(grad_ref);
    CHKERRQ(ierr);
  } else // model POWER_BALANCE_POLAR
  {
    fail += test.computeGradient(opflowtest, X, grad);
  }

  ierr = VecDestroy(&grad);
  CHKERRQ(ierr);

  ierr = OPFLOWDestroy(&opflowtest);
  CHKERRQ(ierr);

  ExaGOFinalize();
  return fail;
}
