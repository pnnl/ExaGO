static char help[] = "SOPFLOW Functionality Tests.\n\n";

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <toml.hpp>
#include <sopflow.h>
#include "toml_utils.hpp"

struct SopflowFunctionalityTestParameters
{
  /* Communicator required to run funcitonality test */
  MPI_Comm comm=MPI_COMM_WORLD;

  /* Parameters required to set up and test a SCOPFlow */
  std::string solver="";
  std::string network="";
  std::string scenfile="";
  std::string contingencies="";
  std::string pload="";
  std::string qload="";
  std::string windgen="";
  std::string description = "";
  int num_contingencies;
  int num_scenarios;
  double tolerance;
  double duration;
  double dT;
  int mode;
  bool multiperiod;
  bool multicontingency;

  /* Parameters used to determine success or failure of functionality test */
  int expected_num_iters;
  double expected_obj_value;

  /* Actual values observed from the system-under-test. */
  double obj_value;
  double error;
  int numiter;
  PetscBool conv_status=PETSC_FALSE;

  /* Assign all member variables from a toml map if values are found */
  void assign_from(toml::value values)
  {
    // Default SOPFLOW settings
    set_if_found(solver, values, "solver");
    set_if_found(network, values, "network");
    set_if_found(scenfile, values, "scenfile");
    set_if_found(num_scenarios, values, "num_scenarios");
    set_if_found(tolerance, values, "tolerance");

    // Multi-contingency SOPFLOW settings
    set_if_found(multicontingency, values, "multicontingency");
    set_if_found(contingencies, values, "contingencies");
    set_if_found(num_contingencies, values, "num_contingencies");

    // Multi-contingency + multi-period SOPFLOW settings
    set_if_found(multiperiod, values, "multiperiod");
    set_if_found(pload, values, "pload");
    set_if_found(qload, values, "qload");
    set_if_found(windgen, values, "windgen");
    set_if_found(duration, values, "duration");
    set_if_found(dT, values, "dT");

    // Expected result + description
    set_if_found(description, values, "description");
    set_if_found(expected_num_iters, values, "num_iters");
    set_if_found(expected_obj_value, values, "obj_value");
  }
};

struct SopflowFunctionalityTests
  : public FunctionalityTestContext<SopflowFunctionalityTestParameters>
{
  SopflowFunctionalityTests(std::string testsuite_filename,
      ExaGOVerbosityLevel logging_verbosity=EXAGO_LOG_INFO)
    : FunctionalityTestContext(testsuite_filename, logging_verbosity)
  {}

  using Params = SopflowFunctionalityTestParameters;
  void ensure_options_are_consistent(toml::value testcase,
      toml::value presets = toml::value{}) override
  {
    auto ensure_option_available = [&] (const std::string& opt) {
      bool is_available = testcase.contains(opt) || presets.contains(opt);
      if (!is_available)
      {
        std::stringstream errs;
        errs << "SOPFLOW Test suite expected option '" << opt
          << "' to be available, but it was not found in this testsuite"
          << " configuration:\n";
        errs << testcase << "\nwith these presets:\n" << presets;
        throw ExaGOError(errs.str().c_str());
      }
    };

    for (const auto& opt : {"solver", "network", "scenfile", "num_scenarios", "tolerance"})
      ensure_option_available(opt);

    bool is_multicontingency = false;
    set_if_found(is_multicontingency, presets, "multicontingency");
    set_if_found(is_multicontingency, testcase, "multicontingency");

    if (is_multicontingency)
      for (const auto& opt : {"contingencies", "num_contingencies"})
        ensure_option_available(opt);

    bool is_multiperiod = false;
    set_if_found(is_multiperiod, presets, "multiperiod");
    set_if_found(is_multiperiod, testcase, "multiperiod");

    if (is_multiperiod)
      for (const auto& opt : {"qload", "pload", "dT", "windgen", "duration"})
        ensure_option_available(opt);
  }

  void initialize_test_parameters(
      Params& params,
      const toml::value& testcase, const toml::value& presets) override
  {
    params.assign_from(presets);
    params.assign_from(testcase);
  }

  toml::value create_failing_testcase(const Params& params) override
  {
    toml::value testcase;

    testcase["description"] = params.description;
    testcase["solver"] = params.solver;
    testcase["network"] = params.network;
    testcase["scenfile"] = params.scenfile;
    testcase["num_scenarios"] = params.num_scenarios;

    testcase["multicontingency"] = params.multicontingency;
    if (params.multicontingency)
    {
      testcase["contingencies"] = params.contingencies;
      testcase["num_contingencies"] = params.num_contingencies;
    }

    testcase["num_iters"] = params.expected_num_iters;
    testcase["observed_num_iters"] = params.numiter;
    testcase["obj_value"] = params.expected_obj_value;
    testcase["observed_obj_value"] = params.obj_value;
    testcase["scaled_objective_value_error"] = params.error;
    testcase["tolerance"] = params.tolerance;
    testcase["did_sopflow_converge"] = params.conv_status;

    testcase["multiperiod"] = params.multiperiod;
    if (params.multiperiod)
    {
      testcase["duration"] = params.duration;
      testcase["dT"] = params.dT;
      testcase["pload"] = params.pload;
      testcase["qload"] = params.qload;
      testcase["windgen"] = params.windgen;
    }

    return testcase;
  }

  void run_test_case(Params& params) override
  {
    PetscErrorCode ierr;
    SOPFLOW sopflow;

    std::cout<<"Test Description: "<<params.description<<std::endl;
    ierr = SOPFLOWCreate(params.comm,&sopflow);ExaGOCheckError(ierr);

    ierr = SOPFLOWSetTolerance(sopflow, params.tolerance);ExaGOCheckError(ierr);

    // Prepend installation directory to network path
    resolve_datafiles_path(params.network);
    ierr = SOPFLOWSetNetworkData(sopflow,params.network.c_str());ExaGOCheckError(ierr);

    // Prepend installation directory to scenario data
    resolve_datafiles_path(params.scenfile);
    ierr = SOPFLOWSetScenarioData(sopflow, SOPFLOW_NATIVE,WIND,params.scenfile.c_str());ExaGOCheckError(ierr);

    ierr = SOPFLOWSetNumScenarios(sopflow, params.num_scenarios);ExaGOCheckError(ierr);

    // TODO:
    // - Implement multi contingency
    // - Implement multi contingency + multiperiod

    ierr = SOPFLOWSetSolver(sopflow, params.solver.c_str());ExaGOCheckError(ierr);
    ierr = SOPFLOWSetUp(sopflow);ExaGOCheckError(ierr);
    ierr = SOPFLOWSolve(sopflow);ExaGOCheckError(ierr);

    /* Possible ways for a funcitonality test to fail */
    bool converge_failed = false;
    bool obj_failed = false;
    bool num_iter_failed = false;

    /* Test convergence status */
    ierr = SOPFLOWGetConvergenceStatus(sopflow,&params.conv_status);ExaGOCheckError(ierr);
    if (params.conv_status==PETSC_FALSE)
    {
      converge_failed = true;
    }

    /* Test objective value */
    ierr = SOPFLOWGetObjective(sopflow,&params.obj_value);ExaGOCheckError(ierr);
    if (!IsEqual(params.obj_value,params.expected_obj_value,params.tolerance,params.error))
    {
      obj_failed = true;
    }

    /* Test num iterations */
    ierr = SOPFLOWGetNumIterations(sopflow,&params.numiter);ExaGOCheckError(ierr);
    if (params.numiter!=params.expected_num_iters)
    {
      num_iter_failed = true;
    }

    /* Did the current functionality test fail in any way? */
    bool local_fail = converge_failed || obj_failed || num_iter_failed;

    if (local_fail)
      fail();
    else
      pass();

    ierr = SOPFLOWDestroy(&sopflow);ExaGOCheckError(ierr);
  }
};

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
  char appname[]="sopflow";
  int _argc = 0;
  ierr = ExaGOInitialize(comm,&argc,&argv,appname,help);

  ExaGOLog(EXAGO_LOG_INFO,"%s\n","Creating SOPFLOW Functionality Test");

  SopflowFunctionalityTests test{std::string(argv[1])};
  test.run_all_test_cases();
  test.print_report();

  ExaGOFinalize();
  return test.failures();
}
