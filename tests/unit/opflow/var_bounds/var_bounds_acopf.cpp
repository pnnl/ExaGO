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
 * @brief Unit test driver for variable bounds
 * @see /tests/unit/opflow/opflow_tests.h for kernel tested by this driver
 *
 * You can pass the following option to the varBoundsAcopf executatable through
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
  Vec X;
  Vec Xl;
  Vec Xu;
  int nx, nconeq, nconineq;
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
  //VecView(Xl,PETSC_VIEWER_STDOUT_SELF);
  //VecView(Xu, PETSC_VIEWER_STDOUT_SELF);

  ierr = OPFLOWGetSolution(opflowtest, &X);
  CHKERRQ(ierr);

  ierr = OPFLOWGetSizes(opflowtest, &nx, &nconeq, &nconineq);
  CHKERRQ(ierr);

  fail += test.computeVariableBounds(opflowtest, Xl, Xu); 
  
  ierr = OPFLOWDestroy(&opflowtest);
  CHKERRQ(ierr);
  ExaGOFinalize();
  return fail;
}

