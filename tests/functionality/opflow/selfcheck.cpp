static char help[] = "OPFLOW Functionality Tests.\n\n";

#include <iostream>
#include <iomanip>
#include <string>
#include <unordered_map>
#include <sstream>
#include <toml.hpp>
#include <opflow.h>
#include <exago_config.h>
#include <utils.hpp>
#include "toml_utils.hpp"
#include <version.hpp>

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

  for (const auto& opt : {"solver", "model", "network", "gen_bus_voltage_type",
      "tolerance", "has_gen_set_point"})
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
  char appname[]="opflow";
  int _argc = 0;
  ierr = ExaGOInitialize(comm,&argc,&argv,appname,help);

  ExaGOLog(EXAGO_LOG_INFO,"%s\n","Creating OPFlow Functionality Test");
  ExaGOLog(EXAGO_LOG_INFO,"Using TOML testsuite file %s\n", argv[1]);
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

  /* Parameters used to set up an OPFLOW */
  std::string solver;
  std::string model;
  std::string network;
  bool is_loadloss_active;
  bool is_powerimb_active;
  int gen_bus_voltage_type;
  double tolerance;
  bool has_gen_set_point;
  std::string description;
  std::string hiop_compute_mode;

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
      set_if_found(solver, presets, "solver");
      set_if_found(model, presets, "model");
      set_if_found(network, presets, "network");
      set_if_found(description, presets, "description");
      set_if_found(is_loadloss_active, presets, "is_loadloss_active");
      set_if_found(is_powerimb_active, presets, "is_powerimb_active");
      set_if_found(expected_num_iters, presets, "num_iters");
      set_if_found(gen_bus_voltage_type, presets, "gen_bus_voltage_type");
      set_if_found(expected_obj_value, presets, "obj_value");
      set_if_found(hiop_compute_mode, presets, "hiop_compute_mode");
      set_if_found(tolerance, presets, "tolerance");
      set_if_found(has_gen_set_point, presets, "has_gen_set_point");
      ensure_options_are_consistent(testcase, presets);
    }
    else
    {
      // If there are no presets, ensure that each test case is internally
      // consistent on its own
      ensure_options_are_consistent(testcase);
    }

    /* If values are found for an individual test case, let these overwrite
     * the global values */
    set_if_found(solver, testcase, "solver");
    set_if_found(model, testcase, "model");
    set_if_found(network, testcase, "network");
    set_if_found(description, testcase, "description");
    set_if_found(is_loadloss_active, testcase, "is_loadloss_active");
    set_if_found(is_powerimb_active, testcase, "is_powerimb_active");
    set_if_found(expected_num_iters, testcase, "num_iters");
    set_if_found(gen_bus_voltage_type, testcase, "gen_bus_voltage_type");
    set_if_found(expected_obj_value, testcase, "obj_value");
    set_if_found(hiop_compute_mode, testcase, "hiop_compute_mode");
    set_if_found(tolerance, testcase, "tolerance");
    set_if_found(has_gen_set_point, testcase, "has_gen_set_point");

    /* Possible ways for a funcitonality test to fail */
    bool converge_failed = false;
    bool obj_failed = false;
    bool num_iter_failed = false;

    PetscBool conv_status=PETSC_FALSE;
    OPFLOW opflow;
    ierr = OPFLOWCreate(comm,&opflow);CHKERRQ(ierr);

    ierr = OPFLOWSetGenBusVoltageType(opflow,
        static_cast<OPFLOWGenBusVoltageType>(gen_bus_voltage_type));CHKERRQ(ierr);

    ierr = OPFLOWSetTolerance(opflow, tolerance);CHKERRQ(ierr);

    // Prepend installation directory to network path
    resolve_datafiles_path(network);
    ierr = OPFLOWReadMatPowerData(opflow,network.c_str());CHKERRQ(ierr);

    ierr = OPFLOWSetSolver(opflow, solver.c_str());CHKERRQ(ierr);

    ierr = OPFLOWSetModel(opflow, model.c_str());CHKERRQ(ierr);

    ierr = OPFLOWHasBusPowerImbalance(opflow, (PetscBool)is_powerimb_active);CHKERRQ(ierr);

    ierr = OPFLOWHasLoadLoss(opflow, (PetscBool)is_loadloss_active);CHKERRQ(ierr);

    ierr = OPFLOWHasGenSetPoint(opflow,(PetscBool) has_gen_set_point);

    ierr = OPFLOWSetUp(opflow);CHKERRQ(ierr);

    if (solver == "HIOP") {
      ierr = OPFLOWSetHIOPComputeMode(opflow, hiop_compute_mode.c_str());CHKERRQ(ierr);
    }

    std::cout<<"Test Description: "<<description<<std::endl;
    ierr = OPFLOWSolve(opflow);CHKERRQ(ierr);

    /* Test convergence status */
    ierr = OPFLOWGetConvergenceStatus(opflow,&conv_status);CHKERRQ(ierr);
    if (conv_status==PETSC_FALSE)
    {
      converge_failed = true;
    }

    /* Test objective value */
    ierr = OPFLOWGetObjective(opflow,&obj_value);CHKERRQ(ierr);
    if (!IsEqual(obj_value,expected_obj_value,tolerance,error))
    {
      obj_failed = true;
    }

    /* Test num iterations */
    ierr = OPFLOWGetNumIterations(opflow,&numiter);CHKERRQ(ierr);
    if (numiter!=expected_num_iters)
    {
      num_iter_failed = true;
    }

    ierr = OPFLOWDestroy(&opflow);CHKERRQ(ierr);

    /* Did the current functionality test fail in any way? */
    bool local_fail = converge_failed || obj_failed || num_iter_failed;

    if (local_fail)
    {
      /* Print summary of configuration */
      summary << "[[testcase]]\n";
      fmt_row(summary, col_width, "solver", solver);
      fmt_row(summary, col_width, "model", model);
      fmt_row(summary, col_width, "network", network);
      fmt_row(summary, col_width, "is_powerimb_active", bool2str(is_powerimb_active));
      fmt_row(summary, col_width, "is_loadloss_active", bool2str(is_loadloss_active));
      fmt_row(summary, col_width, "num_iters", expected_num_iters);
      fmt_comment(summary, col_width, "Actual Number of Iterations", numiter);
      fmt_row(summary, col_width, "obj_value", expected_obj_value);
      fmt_comment(summary, col_width, "Actual Objective Value", obj_value);
      fmt_comment(summary, col_width, "Scaled Objective Value Error", error);
      fmt_row(summary, col_width, "tolerance", tolerance);
      fmt_row(summary, col_width, "gen_bus_voltage_type", gen_bus_voltage_type);
      fmt_comment(summary, col_width, "Did OPFLOW Converge", bool2str(conv_status));
      fmt_row(summary, col_width, "hiop_compute_mode", hiop_compute_mode);
      fmt_row(summary, col_width, "has_gen_set_point", bool2str(has_gen_set_point));
    }

    ExaGOLog(EXAGO_LOG_INFO, "-- %s #%d: %s", testsuite_name.c_str(), total,
        (local_fail?"FAIL":"PASS"));

    if (local_fail)
      fail++;
  }

  ExaGOLog(EXAGO_LOG_INFO,"%s%d/%d%s","Selfcheck OPFLOW got ", fail, total, " failures.");

  if (fail)
  {
    ExaGOLog(EXAGO_LOG_ERROR, "%s%s",
        testsuite_name.c_str(),
        ":\ntests with the following configurations failed to match expected values:\n");
    fmt_row(summary, col_width, "testsuite_name",
        "ExaGO Test Suite Automatically Generated from Failed Tests");
    std::cerr << summary.str() << "\n";
  }

  ExaGOFinalize();
  return fail;
}
