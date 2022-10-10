static char help[] = "SOPFLOW Functionality Tests.\n\n";

#include "toml_utils.h"
#include <opflow.h>
#include <sopflow.h>

struct SopflowFunctionalityTestParameters {
  /* Communicator required to run funcitonality test */
  MPI_Comm comm = MPI_COMM_WORLD;

  /* Parameters required to set up and test a SCOPFlow */
  std::string solver = "";
  std::string network = "";
  std::string scenfile = "";
  std::string contingencies = "";
  std::string pload = "";
  std::string qload = "";
  std::string windgen = "";
  std::string description = "";
  int num_contingencies = 0;
  int num_scenarios;
  double tolerance;
  double duration;
  double dT;
  int mode;
  bool multiperiod;
  bool multicontingency;
  std::string initialization_string = "MIDPOINT";
  OPFLOWInitializationType initialization_type;
  std::string gen_bus_voltage_string = "FIXED_WITHIN_QBOUNDS";
  OPFLOWGenBusVoltageType gen_bus_voltage_type;
  std::string subproblem_solver = "IPOPT";
  std::string subproblem_model = "POWER_BALANCE_POLAR";
  bool ignore_lineflow_constraints = false;
  std::string compute_mode = "auto";
  int verbosity_level = 0;
  bool flatten = true;

  /* Parameters used to determine success or failure of functionality test */
  int expected_num_iters;
  double expected_obj_value;

  /* Actual values observed from the system-under-test. */
  double obj_value;
  double error;
  int numiter;
  PetscBool conv_status = PETSC_FALSE;
  std::vector<std::string> reasons_for_failure;

  /* Assign all member variables from a toml map if values are found */
  void assign_from(toml::value values) {
    // Default SOPFLOW settings
    set_if_found(solver, values, "solver");
    set_if_found(network, values, "network");
    set_if_found(scenfile, values, "scenfile");
    set_if_found(num_scenarios, values, "num_scenarios");
    set_if_found(tolerance, values, "tolerance");
    set_if_found(initialization_string, values, "initialization_type");
    set_if_found(gen_bus_voltage_string, values, "gen_bus_voltage_type");
    set_if_found(subproblem_solver, values, "subproblem_solver");
    set_if_found(subproblem_model, values, "subproblem_model");
    set_if_found(ignore_lineflow_constraints, values,
                 "ignore_lineflow_constraints");
    set_if_found(verbosity_level, values, "verbosity_level");
    set_if_found(compute_mode, values, "compute_mode");

    // Multi-contingency SOPFLOW settings
    set_if_found(multicontingency, values, "multicontingency");
    set_if_found(contingencies, values, "contingencies");
    set_if_found(num_contingencies, values, "num_contingencies");
    set_if_found(flatten, values, "flatten");

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

    if (initialization_string == "MIDPOINT") {
      initialization_type = OPFLOWINIT_MIDPOINT;
    } else if (initialization_string == "FROMFILE") {
      initialization_type = OPFLOWINIT_FROMFILE;
    } else if (initialization_string == "ACPF") {
      initialization_type = OPFLOWINIT_ACPF;
    } else if (initialization_string == "FLATSTART") {
      initialization_type = OPFLOWINIT_FLATSTART;
    }

    if (gen_bus_voltage_string == "VARIABLE_WITHIN_BOUNDS") {
      gen_bus_voltage_type = VARIABLE_WITHIN_BOUNDS;
    } else if (gen_bus_voltage_string == "FIXED_WITHIN_QBOUNDS") {
      gen_bus_voltage_type = FIXED_WITHIN_QBOUNDS;
    } else if (gen_bus_voltage_string == "FIXED_AT_SETPOINT") {
      gen_bus_voltage_type = FIXED_AT_SETPOINT;
    }
  }
};

struct SopflowFunctionalityTests
    : public FunctionalityTestContext<SopflowFunctionalityTestParameters> {
  SopflowFunctionalityTests(std::string testsuite_filename,
                            int logging_verbosity = EXAGO_LOG_INFO)
      : FunctionalityTestContext(testsuite_filename, logging_verbosity) {}

  using Params = SopflowFunctionalityTestParameters;
  void
  ensure_options_are_consistent(toml::value testcase,
                                toml::value presets = toml::value{}) override {
    auto ensure_option_available = [&](const std::string &opt) {
      bool is_available = testcase.contains(opt) || presets.contains(opt);
      if (!is_available) {
        std::stringstream errs;
        errs << "SOPFLOW Test suite expected option '" << opt
             << "' to be available, but it was not found in this testsuite"
             << " configuration:\n";
        errs << testcase << "\nwith these presets:\n" << presets;
        throw ExaGOError(errs.str().c_str());
      }
    };

    for (const auto &opt :
         {"solver", "network", "scenfile", "num_scenarios", "tolerance"})
      ensure_option_available(opt);

    bool is_multicontingency = false;
    set_if_found(is_multicontingency, presets, "multicontingency");
    set_if_found(is_multicontingency, testcase, "multicontingency");

    if (is_multicontingency)
      for (const auto &opt : {"contingencies", "num_contingencies"})
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

    testcase["description"] = params.description;
    testcase["solver"] = params.solver;
    testcase["network"] = params.network;
    testcase["scenfile"] = params.scenfile;
    testcase["num_scenarios"] = params.num_scenarios;
    testcase["initialization_type"] = params.initialization_type;

    testcase["multicontingency"] = params.multicontingency;
    if (params.multicontingency) {
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
    if (params.multiperiod) {
      testcase["duration"] = params.duration;
      testcase["dT"] = params.dT;
      testcase["pload"] = params.pload;
      testcase["qload"] = params.qload;
      testcase["windgen"] = params.windgen;
    }
    testcase["reasons_for_failure"] = params.reasons_for_failure;

    return testcase;
  }

  void run_test_case(Params &params) override {
    PetscErrorCode ierr;
    SOPFLOW sopflow;

    std::cout << "Test Description: " << params.description << std::endl;
    ierr = SOPFLOWCreate(params.comm, &sopflow);
    ExaGOCheckError(ierr);

    ierr = SOPFLOWSetTolerance(sopflow, params.tolerance);
    ExaGOCheckError(ierr);

    // Prepend installation directory to network path
    resolve_datafiles_path(params.network);
    ierr = SOPFLOWSetNetworkData(sopflow, params.network.c_str());
    ExaGOCheckError(ierr);

    ierr = SOPFLOWSetNumScenarios(sopflow, params.num_scenarios);
    ExaGOCheckError(ierr);

    // Prepend installation directory to scenario data
    resolve_datafiles_path(params.scenfile);
    ierr = SOPFLOWSetScenarioData(sopflow, SOPFLOW_NATIVE_SINGLEPERIOD, WIND,
                                  params.scenfile.c_str());
    ExaGOCheckError(ierr);

    ierr = SOPFLOWSetInitializationType(sopflow, params.initialization_type);
    ExaGOCheckError(ierr);

    ierr = SOPFLOWSetGenBusVoltageType(sopflow, params.gen_bus_voltage_type);
    ExaGOCheckError(ierr);

    ierr = SOPFLOWSetSubproblemModel(sopflow, params.subproblem_model.c_str());
    ExaGOCheckError(ierr);

    ierr =
        SOPFLOWSetSubproblemSolver(sopflow, params.subproblem_solver.c_str());
    ExaGOCheckError(ierr);

    ierr = SOPFLOWSetIgnoreLineflowConstraints(
        sopflow, (PetscBool)params.ignore_lineflow_constraints);
    ExaGOCheckError(ierr);

    ierr =
        SOPFLOWSetSubproblemComputeMode(sopflow, params.compute_mode.c_str());
    ExaGOCheckError(ierr);

    ierr = SOPFLOWSetSubproblemVerbosityLevel(sopflow, params.verbosity_level);
    ExaGOCheckError(ierr);

    // TODO:
    // - Implement multi contingency
    ierr = SOPFLOWEnableMultiContingency(sopflow,
                                         (PetscBool)params.multicontingency);
    ExaGOCheckError(ierr);
    if (params.multicontingency) {
      resolve_datafiles_path(params.contingencies);
      ierr = SOPFLOWSetContingencyData(sopflow, NATIVE,
                                       params.contingencies.c_str());
      ExaGOCheckError(ierr);
      ierr = SOPFLOWSetNumContingencies(sopflow, params.num_contingencies);
      ExaGOCheckError(ierr);
      ierr = SOPFLOWFlattenContingencies(sopflow, (PetscBool)params.flatten);
      ExaGOCheckError(ierr);
      //   ierr =
      //   SOPFLOWEnableMultiPeriod(sopflow,(PetscBool)params.multiperiod);ExaGOCheckError(ierr);
      if (params.multiperiod) {
        resolve_datafiles_path(params.windgen);
        resolve_datafiles_path(params.pload);
        resolve_datafiles_path(params.qload);
        ierr = SOPFLOWSetWindGenProfile(sopflow, params.windgen.c_str());
        ExaGOCheckError(ierr);
        ierr =
            SOPFLOWSetTimeStepandDuration(sopflow, params.dT, params.duration);
        ExaGOCheckError(ierr);
        ierr = SOPFLOWSetLoadProfiles(sopflow, params.pload.c_str(),
                                      params.qload.c_str());
        ExaGOCheckError(ierr);
      }
    } else {
      ierr = SOPFLOWSetNumContingencies(sopflow, 1);
    }
    // - Implement multi contingency + multiperiod

    ierr = SOPFLOWSetSolver(sopflow, params.solver.c_str());
    ExaGOCheckError(ierr);
    ierr = SOPFLOWSetUp(sopflow);
    ExaGOCheckError(ierr);
    ierr = SOPFLOWSolve(sopflow);
    ExaGOCheckError(ierr);

    /* Possible ways for a funcitonality test to fail */
    bool converge_failed = false;
    bool obj_failed = false;
    bool num_iter_failed = false;
    params.reasons_for_failure.clear();

    /* Test convergence status */
    ierr = SOPFLOWGetConvergenceStatus(sopflow, &params.conv_status);
    ExaGOCheckError(ierr);
    if (params.conv_status == PETSC_FALSE) {
      converge_failed = true;
      params.reasons_for_failure.push_back("failed to converge");
    }

    /* Test objective value */
    ierr = SOPFLOWGetBaseObjective(sopflow, &params.obj_value);
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
    ierr = SOPFLOWGetNumIterations(sopflow, &params.numiter);
    ExaGOCheckError(ierr);
    if (params.expected_num_iters != -1 && params.numiter != params.expected_num_iters) {
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

    ierr = SOPFLOWDestroy(&sopflow);
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
  char appname[] = "sopflow";
  int _argc = 0;
  ierr = ExaGOInitialize(comm, &argc, &argv, appname, help);

  ExaGOLog(EXAGO_LOG_INFO, "{}", "Creating SOPFLOW Functionality Test");

  SopflowFunctionalityTests test{std::string(argv[1])};
  test.run_all_test_cases();
  test.print_report();

  ExaGOFinalize();
  return test.failures();
}
