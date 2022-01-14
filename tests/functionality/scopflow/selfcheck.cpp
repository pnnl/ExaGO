static char help[] = "SCOPFLOW Functionality Tests.\n\n";

#include "toml_utils.h"
#include <scopflow.h>

struct ScopflowFunctionalityTestParameters {
  /* Communicator required to run funcitonality test */
  MPI_Comm comm = MPI_COMM_WORLD;

  /* Parameters required to set up and test a SCOPFlow */
  std::string solver = "";
  std::string model = "";
  std::string network = "";
  std::string contingencies = "";
  std::string pload = "";
  std::string qload = "";
  std::string windgen = "";
  std::string description = "";
  int num_contingencies;
  double tolerance;
  double duration;
  double dT;
  int mode;
  bool multiperiod;

  /* Parameters used to modify underlying opflow */
  std::string opflow_initialization_string;
  int opflow_initialization;
  std::string opflow_genbusvoltage_string;
  int opflow_genbusvoltage;

  /* Parameters used to determine success or failure of functionality test */
  int expected_num_iters;
  double expected_obj_value;
  std::vector<std::string> reasons_for_failure;

  /* Actual values observed from the system-under-test. */
  double obj_value;
  double error;
  int numiter;
  PetscBool conv_status = PETSC_FALSE;

  /* Assign all member variables from a toml map if values are found */
  void assign_from(toml::value values) {
    set_if_found(solver, values, "solver");
    set_if_found(model, values, "model");
    set_if_found(network, values, "network");
    set_if_found(contingencies, values, "contingencies");
    set_if_found(pload, values, "pload");
    set_if_found(qload, values, "qload");
    set_if_found(windgen, values, "windgen");
    set_if_found(description, values, "description");
    set_if_found(num_contingencies, values, "num_contingencies");
    set_if_found(expected_num_iters, values, "num_iters");
    set_if_found(expected_obj_value, values, "obj_value");
    set_if_found(tolerance, values, "tolerance");
    set_if_found(opflow_initialization_string, values, "opflow_initialization");
    set_if_found(opflow_genbusvoltage_string, values, "opflow_genbusvoltage");
    set_if_found(mode, values, "mode");
    set_if_found(multiperiod, values, "multiperiod");
    set_if_found(duration, values, "duration");
    set_if_found(dT, values, "dT");

    if (opflow_genbusvoltage_string == "VARIABLE_WITHIN_BOUNDS") {
      opflow_genbusvoltage = VARIABLE_WITHIN_BOUNDS;
    } else if (opflow_genbusvoltage_string == "FIXED_WITHIN_QBOUNDS") {
      opflow_genbusvoltage = FIXED_WITHIN_QBOUNDS;
    } else if (opflow_genbusvoltage_string == "FIXED_AT_SETPOINT") {
      opflow_genbusvoltage = FIXED_AT_SETPOINT;
    }

    if (opflow_initialization_string == "MIDPOINT") {
      opflow_initialization = OPFLOWINIT_MIDPOINT;
    } else if (opflow_initialization_string == "FROMFILE") {
      opflow_initialization = OPFLOWINIT_FROMFILE;
    } else if (opflow_initialization_string == "ACPF") {
      opflow_initialization = OPFLOWINIT_ACPF;
    } else if (opflow_initialization_string == "FLATSTART") {
      opflow_initialization = OPFLOWINIT_FLATSTART;
    }
  }
};

struct ScopflowFunctionalityTests
    : public FunctionalityTestContext<ScopflowFunctionalityTestParameters> {
  ScopflowFunctionalityTests(std::string testsuite_filename,
                             int logging_verbosity = EXAGO_LOG_INFO)
      : FunctionalityTestContext(testsuite_filename, logging_verbosity) {}

  using Params = ScopflowFunctionalityTestParameters;
  void
  ensure_options_are_consistent(toml::value testcase,
                                toml::value presets = toml::value{}) override {
    auto ensure_option_available = [&](const std::string &opt) {
      bool is_available = testcase.contains(opt) || presets.contains(opt);
      if (!is_available) {
        std::stringstream errs;
        errs << "SCOPFLOW Test suite expected option '" << opt
             << "' to be available, but it was not found in this testsuite"
             << " configuration:\n";
        errs << testcase << "\nwith these presets:\n" << presets;
        throw ExaGOError(errs.str().c_str());
      }
    };

    for (const auto &opt :
         {"solver", "model", "network", "contingencies", "tolerance"})
      ensure_option_available(opt);

    bool is_multiperiod = false;
    set_if_found(is_multiperiod, presets, "multiperiod");
    set_if_found(is_multiperiod, testcase, "multiperiod");

    if (is_multiperiod)
      for (const auto &opt : {"qload", "pload", "dT", "windgen", "duration"})
        ensure_option_available(opt);
  }

  void initialize_test_parameters(Params &params, const toml::value &testcase,
                                  const toml::value &presets) override {
    params.assign_from(presets);
    params.assign_from(testcase);
  }

  toml::value create_failing_testcase(const Params &params) override {
    toml::value testcase;

    testcase["solver"] = params.solver;
    testcase["model"] = params.model;
    testcase["network"] = params.network;
    testcase["contingencies"] = params.contingencies;
    testcase["num_contingencies"] = params.num_contingencies;
    testcase["opflow_initialization"] = params.opflow_initialization;
    testcase["opflow_genbusvoltage"] = params.opflow_genbusvoltage;
    testcase["num_iters"] = params.expected_num_iters;
    testcase["observed_num_iters"] = params.numiter;
    testcase["obj_value"] = params.expected_obj_value;
    testcase["observed_obj_value"] = params.obj_value;
    testcase["scaled_objective_value_error"] = params.error;
    testcase["tolerance"] = params.tolerance;
    testcase["did_scopflow_converge"] = params.conv_status;
    testcase["mode"] = params.mode;
    testcase["reasons_for_failure"] = params.reasons_for_failure;

    testcase["multiperiod"] = params.multiperiod;
    if (params.multiperiod) {
      testcase["duration"] = params.duration;
      testcase["dT"] = params.dT;
      testcase["pload"] = params.pload;
      testcase["qload"] = params.qload;
      testcase["windgen"] = params.windgen;
    }

    return testcase;
  }

  void run_test_case(Params &params) override {
    PetscErrorCode ierr;
    SCOPFLOW scopflow;

    std::cout << "Test Description: " << params.description << std::endl;
    ierr = SCOPFLOWCreate(params.comm, &scopflow);
    ExaGOCheckError(ierr);

    ierr = SCOPFLOWSetTolerance(scopflow, params.tolerance);
    ExaGOCheckError(ierr);

    // Prepend installation directory to network path
    resolve_datafiles_path(params.network);
    ierr = SCOPFLOWSetNetworkData(scopflow, params.network.c_str());
    ExaGOCheckError(ierr);

    // Prepend installation directory to contingency file
    resolve_datafiles_path(params.contingencies);
    ierr = SCOPFLOWSetContingencyData(scopflow, NATIVE,
                                      params.contingencies.c_str());
    ExaGOCheckError(ierr);

    // Prepend installation directory to network path
    resolve_datafiles_path(params.pload);
    ierr = SCOPFLOWSetPLoadData(scopflow, params.pload.c_str());
    ExaGOCheckError(ierr);

    // Prepend installation directory to network path
    resolve_datafiles_path(params.qload);
    ierr = SCOPFLOWSetQLoadData(scopflow, params.qload.c_str());
    ExaGOCheckError(ierr);

    // Prepend installation directory to network path
    resolve_datafiles_path(params.windgen);
    ierr = SCOPFLOWSetWindGenProfile(scopflow, params.windgen.c_str());
    ExaGOCheckError(ierr);

    // Set number of contingencies
    ierr = SCOPFLOWSetNumContingencies(scopflow, params.num_contingencies);
    ExaGOCheckError(ierr);

    ierr = SCOPFLOWEnableMultiPeriod(scopflow, (PetscBool)params.multiperiod);
    ExaGOCheckError(ierr);
    ierr = SCOPFLOWSetTimeStep(scopflow, params.dT);
    ExaGOCheckError(ierr);
    ierr = SCOPFLOWSetDuration(scopflow, params.duration);
    ExaGOCheckError(ierr);

    ierr = SCOPFLOWSetSolver(scopflow, params.solver.c_str());
    ExaGOCheckError(ierr);

    ierr = SCOPFLOWSetInitilizationType(
        scopflow, (OPFLOWInitializationType)params.opflow_initialization);
    ExaGOCheckError(ierr);
    ierr = SCOPFLOWSetGenBusVoltageType(
        scopflow, (OPFLOWGenBusVoltageType)params.opflow_genbusvoltage);
    ExaGOCheckError(ierr);

    ierr = SCOPFLOWSetModel(scopflow, params.model.c_str());
    ExaGOCheckError(ierr);

    ierr = SCOPFLOWSetMode(scopflow, (PetscInt)params.mode);
    ExaGOCheckError(ierr);

    ierr = SCOPFLOWSetUp(scopflow);
    ExaGOCheckError(ierr);

    ierr = SCOPFLOWSolve(scopflow);
    ExaGOCheckError(ierr);

    /* Possible ways for a funcitonality test to fail */
    bool converge_failed = false;
    bool obj_failed = false;
    bool num_iter_failed = false;
    params.reasons_for_failure.clear();

    /* Test convergence status */
    ierr = SCOPFLOWGetConvergenceStatus(scopflow, &params.conv_status);
    ExaGOCheckError(ierr);
    if (params.conv_status == PETSC_FALSE) {
      converge_failed = true;
      params.reasons_for_failure.push_back("failed to converge");
    }

    /* Test objective value */
    ierr = SCOPFLOWGetBaseObjective(scopflow, &params.obj_value);
    ExaGOCheckError(ierr);
    if (!IsEqual(params.obj_value, params.expected_obj_value, params.tolerance,
                 params.error)) {
      obj_failed = true;
      params.reasons_for_failure.push_back(fmt::format(
          "expected objective value={} actual objective value={} tol={} err={}",
          params.expected_obj_value, params.obj_value, params.tolerance,
          params.error));
    }

    /* Test num iterations */
    ierr = SCOPFLOWGetNumIterations(scopflow, &params.numiter);
    ExaGOCheckError(ierr);
    if (params.numiter != params.expected_num_iters) {
      num_iter_failed = true;
      params.reasons_for_failure.push_back(
          fmt::format("expected {} num iters, got {}",
                      params.expected_num_iters, params.numiter));
    }

    /* Did the current functionality test fail in any way? */
    bool local_fail = converge_failed || obj_failed || num_iter_failed;

    if (local_fail)
      fail();
    else
      pass();

    ierr = SCOPFLOWDestroy(&scopflow);
    ExaGOCheckError(ierr);
  }
};

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Pass path to test cases TOML file as the first argument to "
              << "this test driver.\n";
    std::exit(1);
  }

  PetscErrorCode ierr;
  OutputFormat fmt = MATPOWER;
  MPI_Comm comm = MPI_COMM_WORLD;
  char appname[] = "scopflow";
  int _argc = 0;
  ierr = ExaGOInitialize(comm, &argc, &argv, appname, help);

  ExaGOLog(EXAGO_LOG_INFO, "{}", "Creating SCOPFlow Functionality Test");

  ScopflowFunctionalityTests test{std::string(argv[1])};
  test.run_all_test_cases();
  test.print_report();

  ExaGOFinalize();
  return test.failures();
}
