#include <cstdio>
#include <iostream>
#include <string>

#include <exago_config.h>
#include <private/pflowimpl.h>
#include <utils.h>

#include <pflow.h>
#include "pflow/pflow_tests.h"
#include "test_acopf_utils.h"

//#if defined(EXAGO_ENABLE_RAJA)
//#include <RAJA/RAJA.hpp>
//#include <private/raja_exec_config.h>
//#include <umpire/Allocator.hpp>
//#include <umpire/ResourceManager.hpp>
//#endif

/**
 * TODO: Update this section
 * @brief Unit test driver for ACOPF models
 * @see opflow/opflow_tests.h for kernels tested by this driver
 *
 * You can pass several options to the TestAcopf executatable through the
 * command line (implemented using PETSc options):
 *
 *    ~ -netfile <data_file> : Specifies the input data file to test against.
 * Default value is `/<exago_dir>/datafiles/case9/case9mod.m`. See directory
 * datafiles for other potential inputs.
 *
 *    ~ -gen_test_data       : If used, generates an answer key using IPOPT to
 * test against. If not used, uses existing answer keys located in
 * `/<exago_dir>/datafiles/test_validation`.
 *
 *    ~ -write_test_data     : If used, generates test data, and then writes out
 * the results to `/<install_dir>/tests/datafiles/test_validation/<data_file>/`.
 *                             In order to save generated results, you should
 * copy them into `/<exago_dir>/datafiles/test_validation/<data_file>/` and add
 * them to git.
 *
 */
int main(int argc, char **argv) {
  //  const bool isTestPflowModelPBPOL = false;
  //#if defined(EXAGO_ENABLE_RAJA)
  //  const bool isTestPflowModelPBPOLRAJAHIOP = true;
  //  const bool isTestPflowModelPBPOLHIOP = false;
  //#else
  //  const bool isTestPflowModelPBPOLHIOP = true;
  //#endif
  PetscErrorCode ierr;
  PetscBool flg, gen_test_data, write_test_data;
  PetscInt iter;
  PetscBool converged = PETSC_FALSE;
  // Vec X, Xl, Xu, G, Gl, Gu, grad, Lambda;
  Mat Jac;
  int fail = 0;
  PetscLogStage stages[2];
  double obj_value, obj_factor;
  char file_c_str[PETSC_MAX_PATH_LEN];
  char validation_dir_c_str[PETSC_MAX_PATH_LEN];
  std::string file, validation_path;
  char appname[] = "pflow";
  MPI_Comm comm = MPI_COMM_WORLD;

  char help[] = "Unit tests for PFLOW\n";

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

  std::cout << file << std::endl;

  /* Place to store/read reference solutions */
  ierr = PetscOptionsGetString(NULL, NULL, "-validation_dir",
                               validation_dir_c_str, PETSC_MAX_PATH_LEN, &flg);
  CHKERRQ(ierr);

  /* This defaults to within the install directory */
  if (!flg) {
    validation_path = std::string(EXAGO_OPTIONS_DIR) +
                      std::string("../datafiles/test_validation/") +
                      getFileName(file) + "/";
  } else {
    validation_path =
        std::string(validation_dir_c_str) + "/" + getFileName(file) + "/";
  }

  ierr =
      PetscOptionsGetBool(NULL, NULL, "-gen_test_data", NULL, &gen_test_data);
  CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL, NULL, "-write_test_data", NULL,
                             &write_test_data);
  CHKERRQ(ierr);

  ierr = PetscLogStageRegister("Solution key", &stages[0]);
  CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Test stage", &stages[1]);
  CHKERRQ(ierr);

  ierr = PetscLogStagePush(stages[0]);

  if (gen_test_data || write_test_data) {
    std::cout << "Generating test_validation data from existing model."
              << std::endl;

    PFLOW pflow;

    /* Create reference pflow */
    ierr = PFLOWCreate(PETSC_COMM_WORLD, &pflow);
    CHKERRQ(ierr);

    /* Read Network data */
    ierr = PFLOWReadMatPowerData(pflow, file.c_str());
    CHKERRQ(ierr);

    /* Solve PFLOW to get the reference solution */
    ierr = PFLOWSolve(pflow);

    /* Get the reference solution */
    ierr = PFLOWConverged(pflow, &converged);
    CHKERRQ(ierr);

    ierr = PFLOWGetNumIterations(pflow, &iter);
    CHKERRQ(ierr);

    ierr = PFLOWCreateMatrix(pflow, &Jac);
    CHKERRQ(ierr);

    ierr = PFLOWGetJacobian(pflow, &Jac);
    CHKERRQ(ierr);

    if (write_test_data) {
      std::cout << "Writing answer keys to " << validation_path << std::endl;

      validate_directory(validation_path);
      saveToFile(Jac, validation_path + "Jacpf_valid.bin");
    }
  } else {
    std::cout << "Loading answer key from " << validation_path << std::endl;
    readFromFile(&Jac, validation_path + "Jacpf_valid.bin");
  }

  ierr = PetscLogStagePop();
  CHKERRQ(ierr);

  ierr = PetscLogStagePush(stages[1]);
  CHKERRQ(ierr);

  exago::tests::TestPflow test;
  PFLOW pflowtest;

  ierr = PFLOWCreate(PETSC_COMM_WORLD, &pflowtest);
  CHKERRQ(ierr);

  ierr = PFLOWReadMatPowerData(pflowtest, file.c_str());
  CHKERRQ(ierr);

  ierr = PFLOWSetUp(pflowtest);
  CHKERRQ(ierr);

  ierr = PFLOWSolve(pflowtest);
  CHKERRQ(ierr);

  ierr = PFLOWGetConvergenceStatus(pflowtest, &converged);
  CHKERRQ(ierr);

  std::cout << "Convergence status: " << (bool)converged << std::endl;

  std::cout << "testing Jacobian" << std::endl;
  fail += test.computeJacobian(pflowtest, Jac);

  // ierr = PFLOWDestroy(&pflowtest);
  // CHKERRQ(ierr);

  /* Destroy PFLOW objects */
  ierr = MatDestroy(&Jac);
  CHKERRQ(ierr);

  ExaGOFinalize();
  // PetscFinalize();
  return fail;
}
