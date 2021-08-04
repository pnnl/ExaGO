static char help[] = "SCOPFLOW Functionality Tests.\n\n";

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <toml.hpp>
#include <scopflow.h>
#include "toml_utils.hpp"

// Ensure that options for a given test case are internally consistent and that
// all required parameters are defined.
void ensure_options_are_consistent(toml::value testcase,
    toml::value presets = toml::value{})
{
  auto ensure_option_available = [&] (const std::string& opt) {
    bool is_available = testcase.contains(opt) || presets.contains(opt);
    if (!is_available)
    {
      std::stringstream errs;
      errs << "SCOPFLOW Test suite expected option '" << opt
        << "' to be available, but it was not found in this testsuite"
        << " configuration:\n";
      errs << testcase << "\nwith these presets:\n" << presets;
      throw ExaGOError(errs.str().c_str());
    }
  };

  for (const auto& opt : {"solver", "model", "network", "contingencies", "tolerance"})
    ensure_option_available(opt);

  bool is_multiperiod = false;
  set_if_found(is_multiperiod, presets, "multiperiod");
  set_if_found(is_multiperiod, testcase, "multiperiod");

  if (is_multiperiod)
    for (const auto& opt : {"qload", "pload", "dT", "windgen", "duration"})
      ensure_option_available(opt);
}

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::cerr << "Pass path to test cases TOML file as the first argument to "
      << "this test driver.\n";
    std::exit(1);
  }

  PetscErrorCode ierr;
  OutputFormat fmt=MATPOWER;
  MPI_Comm comm=MPI_COMM_WORLD;
  char appname[]="scopflow";
  int _argc = 0;
  ierr = ExaGOInitialize(comm,&argc,&argv,appname,help);

  ExaGOLog(EXAGO_LOG_INFO,"%s\n","Creating OPFlow Functionality Test");
  const auto testsuite = toml::parse(argv[1]);
  const auto testcases = toml::find(testsuite, "testcase");

  const auto testsuite_name = toml::find<std::string>(testsuite, "testsuite_name");
  ExaGOLog(EXAGO_LOG_INFO, "Running test suite '%s'\n", testsuite_name.c_str());

  int total = 0;
  int fail = 0;

  /* Stream used for formatting results */
  std::stringstream summary;
  summary.precision(12);
  summary << std::fixed;
  static constexpr int col_width = 35;

  /* Parameters used to set up an SCOPFLOW */
  std::string solver;
  std::string model;
  std::string network;
  std::string contingencies;
  std::string pload;
  std::string qload;
  std::string windgen;
  int num_contingencies;
  double tolerance;
  double duration;
  double dT;
  int mode;
  bool multiperiod;

  /* Parameters used to modify underlying opflow */
  int opflow_initialization;
  int opflow_genbusvoltage;

  /* Parameters used to determine success or failure of functionality test */
  int expected_num_iters;
  double expected_obj_value;

  /* Iterate over all test cases and run individual functionality tests */
  for (auto& testcase : testcases.as_array()) {
    total++;
    double obj_value, error;
    int numiter = 0;

    /* If presets are defined for the test suite, use these initially.
     * The fallback values in this scope are pretty arbitrary - better defaults
     * should be used. */
    if (testsuite.contains("presets"))
    {
      auto presets = toml::find(testsuite, "presets");

      ensure_options_are_consistent(testcase, presets);

      set_if_found(solver, presets, "solver");
      set_if_found(model, presets, "model");
      set_if_found(network, presets, "network");
      set_if_found(contingencies, presets, "contingencies");
      set_if_found(pload, presets, "pload");
      set_if_found(qload, presets, "qload");
      set_if_found(windgen, presets, "windgen");
      set_if_found(num_contingencies, presets, "num_contingencies");
      set_if_found(expected_num_iters, presets, "num_iters");
      set_if_found(expected_obj_value, presets, "obj_value");
      set_if_found(tolerance, presets, "tolerance");
      set_if_found(opflow_initialization, presets, "opflow_initialization");
      set_if_found(opflow_genbusvoltage, presets, "opflow_genbusvoltage");
      set_if_found(mode, presets, "mode");
      set_if_found(multiperiod, presets, "multiperiod");
      set_if_found(duration, presets, "duration");
      set_if_found(dT, presets, "dT");
    }
    else
    {
      // If there are no presets, ensure that each test case is internally
      // consistent
      ensure_options_are_consistent(testcase);
    }

    /* If values are found for an individual test case, let these overwrite
     * the global values */
    set_if_found(solver, testcase, "solver");
    set_if_found(model, testcase, "model");
    set_if_found(network, testcase, "network");
    set_if_found(contingencies, testcase, "contingencies");
    set_if_found(pload, testcase, "pload");
    set_if_found(qload, testcase, "qload");
    set_if_found(windgen, testcase, "windgen");
    set_if_found(num_contingencies, testcase, "num_contingencies");
    set_if_found(expected_num_iters, testcase, "num_iters");
    set_if_found(expected_obj_value, testcase, "obj_value");
    set_if_found(tolerance, testcase, "tolerance");
    set_if_found(opflow_initialization, testcase, "opflow_initialization");
    set_if_found(opflow_genbusvoltage, testcase, "opflow_genbusvoltage");
    set_if_found(mode, testcase, "mode");
    set_if_found(multiperiod, testcase, "multiperiod");
    set_if_found(duration, testcase, "duration");
    set_if_found(dT, testcase, "dT");

    SCOPFLOW scopflow;
    ierr = SCOPFLOWCreate(comm,&scopflow);CHKERRQ(ierr);

    ierr = SCOPFLOWSetTolerance(scopflow, tolerance);CHKERRQ(ierr);

    // Prepend installation directory to network path
    resolve_datafiles_path(network);
    ierr = SCOPFLOWSetNetworkData(scopflow,network.c_str());CHKERRQ(ierr);

    // Prepend installation directory to contingency file
    resolve_datafiles_path(contingencies);
    ierr = SCOPFLOWSetContingencyData(scopflow, NATIVE, contingencies.c_str());CHKERRQ(ierr);

    // Prepend installation directory to network path
    resolve_datafiles_path(pload);
    ierr = SCOPFLOWSetPLoadData(scopflow,pload.c_str());CHKERRQ(ierr);

    // Prepend installation directory to network path
    resolve_datafiles_path(qload);
    ierr = SCOPFLOWSetQLoadData(scopflow,qload.c_str());CHKERRQ(ierr);

    // Prepend installation directory to network path
    resolve_datafiles_path(windgen);
    ierr = SCOPFLOWSetWindGenProfile(scopflow,windgen.c_str());CHKERRQ(ierr);

    // Set number of contingencies
    ierr = SCOPFLOWSetNumContingencies(scopflow, num_contingencies);CHKERRQ(ierr);

    ierr = SCOPFLOWEnableMultiPeriod(scopflow, (PetscBool)multiperiod);CHKERRQ(ierr);
    ierr = SCOPFLOWSetTimeStep(scopflow, dT);CHKERRQ(ierr);
    ierr = SCOPFLOWSetDuration(scopflow, duration);CHKERRQ(ierr);

    ierr = SCOPFLOWSetSolver(scopflow, solver.c_str());CHKERRQ(ierr);

    ierr = SCOPFLOWSetInitilizationType(scopflow, (OPFLOWInitializationType)opflow_initialization);CHKERRQ(ierr);
    ierr = SCOPFLOWSetGenBusVoltageType(scopflow, (OPFLOWGenBusVoltageType)opflow_genbusvoltage);CHKERRQ(ierr);

    ierr = SCOPFLOWSetModel(scopflow, model.c_str());CHKERRQ(ierr);
    
    ierr = SCOPFLOWSetMode(scopflow, (PetscInt)mode);CHKERRQ(ierr);

    ierr = SCOPFLOWSetUp(scopflow);CHKERRQ(ierr);

    ierr = SCOPFLOWSolve(scopflow);CHKERRQ(ierr);

    /* Possible ways for a funcitonality test to fail */
    bool converge_failed = false;
    bool obj_failed = false;
    bool num_iter_failed = false;

    /* Test convergence status */
    PetscBool conv_status=PETSC_FALSE;
    ierr = SCOPFLOWGetConvergenceStatus(scopflow,&conv_status);CHKERRQ(ierr);
    if (conv_status==PETSC_FALSE)
    {
      converge_failed = true;
    }

    /* Test objective value */
    ierr = SCOPFLOWGetObjective(scopflow,&obj_value);CHKERRQ(ierr);
    if (!IsEqual(obj_value,expected_obj_value,tolerance,error))
    {
      obj_failed = true;
    }

    /* Test num iterations */
    ierr = SCOPFLOWGetNumIterations(scopflow,&numiter);CHKERRQ(ierr);
    if (numiter!=expected_num_iters)
    {
      num_iter_failed = true;
    }

    /* Did the current functionality test fail in any way? */
    bool local_fail = converge_failed || obj_failed || num_iter_failed;

    if (local_fail)
    {
      /* Print summary of configuration */
      summary << "[[testcase]]\n";
      fmt_row(summary, col_width, "solver", solver);
      fmt_row(summary, col_width, "model", model);
      fmt_row(summary, col_width, "network", network);
      fmt_row(summary, col_width, "pload", pload);
      fmt_row(summary, col_width, "qload", qload);
      fmt_row(summary, col_width, "windgen", windgen);
      fmt_row(summary, col_width, "contingencies", contingencies);
      fmt_row(summary, col_width, "num_contingencies", num_contingencies);
      fmt_row(summary, col_width, "opflow_intitialization", opflow_initialization);
      fmt_row(summary, col_width, "opflow_genbusvoltage", opflow_genbusvoltage);
      fmt_row(summary, col_width, "multiperiod", bool2str(multiperiod));
      fmt_row(summary, col_width, "duration", duration);
      fmt_row(summary, col_width, "dT", dT);
      fmt_row(summary, col_width, "mode", mode);
      fmt_row(summary, col_width, "num_iters", expected_num_iters);
      fmt_comment(summary, col_width, "Actual Number of Iterations", numiter);
      fmt_row(summary, col_width, "obj_value", expected_obj_value);
      fmt_comment(summary, col_width, "Actual Objective Value", obj_value);
      fmt_comment(summary, col_width, "Scaled Objective Value Error", error);
      fmt_row(summary, col_width, "tolerance", tolerance);
      fmt_comment(summary, col_width, "Did SCOPFLOW Converge", bool2str(conv_status));
    }

    ExaGOLog(EXAGO_LOG_INFO, "-- %s #%d: %s", testsuite_name.c_str(), total,
        (local_fail?"FAIL":"PASS"));

    ierr = SCOPFLOWDestroy(&scopflow);CHKERRQ(ierr);
    fail += (local_fail ? 1 : 0);
  }

  ExaGOLog(EXAGO_LOG_INFO,"%s%d/%d%s","Selfcheck SCOPFLOW got ", fail, total, " failures.");

  if (fail)
  {
    ExaGOLog(EXAGO_LOG_ERROR, "%s%s",
        testsuite_name.c_str(),
        ":\ntests with the following configurations failed to match expected values:\n");
    std::cerr << "# begin autogenerated TOML test suite\n";
    fmt_row(std::cerr, col_width, "testsuite_name",
        "'ExaGO Test Suite Automatically Generated from Failed Tests'");
    std::cerr << summary.str();
    std::cerr << "# end autogenerated TOML test suite\n";
  }
  ExaGOFinalize();
  return fail;
}
