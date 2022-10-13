static char help[] = "OPFLOW Functionality Tests.\n\n";

#include "toml_utils.h"
#include <opflow.h>

struct OpflowFunctionalityTestParameters {
  /* Communicator required to run funcitonality test */
  MPI_Comm comm = MPI_COMM_WORLD;

  /* Parameters required to set up and test a OPFlow */
  std::string solver = "";
  std::string model = "";
  std::string network = "";
  bool is_loadloss_active;
  bool is_powerimb_active;
  std::string gen_bus_voltage_string;
  OPFLOWGenBusVoltageType gen_bus_voltage_type;
  double tolerance;
  double has_gen_set_point;
  bool use_agc = false;
  double load_loss_penalty = 1000.0;
  double power_imbalance_penalty = 1000.0;
  std::string description = "";
  int hiop_verbosity_level = 0;
  std::string hiop_compute_mode;
  std::string hiop_mem_space;
  std::string initialization_string = "MIDPOINT";
  OPFLOWInitializationType initialization_type;

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
    set_if_found(solver, values, "solver");
    set_if_found(model, values, "model");
    set_if_found(network, values, "network");
    set_if_found(is_loadloss_active, values, "is_loadloss_active");
    set_if_found(is_powerimb_active, values, "is_powerimb_active");
    set_if_found(gen_bus_voltage_string, values, "gen_bus_voltage_type");
    set_if_found(description, values, "description");
    set_if_found(expected_num_iters, values, "num_iters");
    set_if_found(expected_obj_value, values, "obj_value");
    set_if_found(tolerance, values, "tolerance");
    set_if_found(use_agc, values, "use_agc");
    set_if_found(load_loss_penalty, values, "load_loss_penalty");
    set_if_found(power_imbalance_penalty, values, "power_imbalance_penalty");
    set_if_found(initialization_string, values, "initialization_type");
    set_if_found(hiop_compute_mode, values, "hiop_compute_mode");
    set_if_found(hiop_mem_space, values, "hiop_mem_space");

    if (gen_bus_voltage_string == "VARIABLE_WITHIN_BOUNDS") {
      gen_bus_voltage_type = VARIABLE_WITHIN_BOUNDS;
    } else if (gen_bus_voltage_string == "FIXED_WITHIN_QBOUNDS") {
      gen_bus_voltage_type = FIXED_WITHIN_QBOUNDS;
    } else if (gen_bus_voltage_string == "FIXED_AT_SETPOINT") {
      gen_bus_voltage_type = FIXED_AT_SETPOINT;
    }

    if (initialization_string == "MIDPOINT") {
      initialization_type = OPFLOWINIT_MIDPOINT;
    } else if (initialization_string == "FROMFILE") {
      initialization_type = OPFLOWINIT_FROMFILE;
    } else if (initialization_string == "ACPF") {
      initialization_type = OPFLOWINIT_ACPF;
    } else if (initialization_string == "FLATSTART") {
      initialization_type = OPFLOWINIT_FLATSTART;
    } else if (initialization_string == "DCOPF") {
      initialization_type = OPFLOWINIT_DCOPF;
    }
  }
};

struct OpflowFunctionalityTests
    : public FunctionalityTestContext<OpflowFunctionalityTestParameters> {
  OpflowFunctionalityTests(std::string testsuite_filename,
                           int logging_verbosity = EXAGO_LOG_INFO)
      : FunctionalityTestContext(testsuite_filename, logging_verbosity) {}

  using Params = OpflowFunctionalityTestParameters;
  void
  ensure_options_are_consistent(toml::value testcase,
                                toml::value presets = toml::value{}) override {
    auto ensure_option_available = [&](const std::string &opt) {
      bool is_available = testcase.contains(opt) || presets.contains(opt);
      if (!is_available) {
        std::stringstream errs;
        errs << "OPFLOW Test suite expected option '" << opt
             << "' to be available, but it was not found in this testsuite"
             << " configuration:\n";
        errs << testcase << "\nwith these presets:\n" << presets;
        throw ExaGOError(errs.str().c_str());
      }
    };

    for (const auto &opt :
         {"solver", "model", "network", "gen_bus_voltage_type", "tolerance",
          "has_gen_set_point"})
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
    testcase["is_loadloss_active"] = params.is_loadloss_active;
    testcase["is_powerimb_active"] = params.is_powerimb_active;
    testcase["gen_bus_voltage_type"] = params.gen_bus_voltage_type;
    testcase["description"] = params.description;
    testcase["num_iters"] = params.expected_num_iters;
    testcase["observed_num_iters"] = params.numiter;
    testcase["tolerance"] = params.tolerance;
    testcase["use_agc"] = params.use_agc;
    testcase["load_loss_penalty"] = params.load_loss_penalty;
    testcase["power_imbalance_penalty"] = params.power_imbalance_penalty;
    testcase["initialization_type"] = params.initialization_type;
    testcase["hiop_compute_mode"] = params.hiop_compute_mode;
    testcase["hiop_verbosity_level"] = params.hiop_verbosity_level;
    testcase["hiop_mem_space"] = params.hiop_mem_space;
    testcase["obj_value"] = params.expected_obj_value;
    testcase["observed_obj_value"] = params.obj_value;
    testcase["scaled_objective_value_error"] = params.error;
    testcase["did_opflow_converge"] = params.conv_status;
    testcase["reasons_for_failure"] = params.reasons_for_failure;

    return testcase;
  }

  void run_test_case(Params &params) override {
    PetscErrorCode ierr;
    OPFLOW opflow;

    std::cout << "Test Description: " << params.description << std::endl;
    ierr = OPFLOWCreate(params.comm, &opflow);
    ExaGOCheckError(ierr);

    ierr = OPFLOWSetHIOPVerbosityLevel(opflow, params.hiop_verbosity_level);
    ExaGOCheckError(ierr);

    ierr = OPFLOWSetHIOPMemSpace(opflow, params.hiop_mem_space.c_str());
    ExaGOCheckError(ierr);

   ierr = OPFLOWSetTolerance(opflow, params.tolerance);
    ExaGOCheckError(ierr);

    // Prepend installation directory to network path
    resolve_datafiles_path(params.network);
    ierr = OPFLOWReadMatPowerData(opflow, params.network.c_str());
    ExaGOCheckError(ierr);

    ierr = OPFLOWSetGenBusVoltageType(opflow, params.gen_bus_voltage_type);
    ExaGOCheckError(ierr);

    ierr = OPFLOWSetSolver(opflow, params.solver.c_str());
    ExaGOCheckError(ierr);

    ierr = OPFLOWSetModel(opflow, params.model.c_str());
    ExaGOCheckError(ierr);

    ierr = OPFLOWHasBusPowerImbalance(opflow,
                                      (PetscBool)params.is_powerimb_active);
    ExaGOCheckError(ierr);

    ierr = OPFLOWHasLoadLoss(opflow, (PetscBool)params.is_loadloss_active);
    ExaGOCheckError(ierr);

    ierr = OPFLOWHasGenSetPoint(opflow, (PetscBool)params.has_gen_set_point);
    ExaGOCheckError(ierr);

    ierr = OPFLOWUseAGC(opflow, (PetscBool)params.use_agc);
    ExaGOCheckError(ierr);
    ExaGOCheckError(ierr);

    ierr = OPFLOWSetLoadLossPenalty(opflow, params.load_loss_penalty);
    ExaGOCheckError(ierr);

    ierr = OPFLOWSetBusPowerImbalancePenalty(opflow,
                                             params.power_imbalance_penalty);
    ExaGOCheckError(ierr);

    ierr = OPFLOWSetInitializationType(opflow, params.initialization_type);
    ExaGOCheckError(ierr);

    ierr = OPFLOWSetUp(opflow);
    ExaGOCheckError(ierr);

    if (params.solver == "HIOP") {
      ierr = OPFLOWSetHIOPComputeMode(opflow, params.hiop_compute_mode.c_str());
      ExaGOCheckError(ierr);
    }

    ierr = OPFLOWSolve(opflow);
    ExaGOCheckError(ierr);

    /* Possible ways for a functionality test to fail */
    bool converge_failed = false;
    bool obj_failed = false;
    bool num_iter_failed = false;
    params.reasons_for_failure.clear();

    /* Test convergence status */
    ierr = OPFLOWGetConvergenceStatus(opflow, &params.conv_status);
    ExaGOCheckError(ierr);
    if (params.conv_status == PETSC_FALSE) {
      converge_failed = true;
      params.reasons_for_failure.push_back("failed to converge");
    }

    /* Test objective value */
    ierr = OPFLOWGetObjective(opflow, &params.obj_value);
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
    ierr = OPFLOWGetNumIterations(opflow, &params.numiter);
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

    ierr = OPFLOWDestroy(&opflow);
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
  char appname[] = "opflow";
  int _argc = 0;
  ierr = ExaGOInitialize(comm, &argc, &argv, appname, help);

  ExaGOLog(EXAGO_LOG_INFO, "{}", "Creating OPFlow Functionality Test");

  OpflowFunctionalityTests test{std::string(argv[1])};
  test.run_all_test_cases();
  test.print_report();

  ExaGOFinalize();
  return test.failures();
}
