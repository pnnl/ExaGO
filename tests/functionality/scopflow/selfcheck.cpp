static char help[] = "SCOPFLOW Functionality Tests.\n\n";

#include "toml_utils.h"
#include <scopflow.h>

struct ScopflowFunctionalityTestParameters {
  /* Communicator required to run funcitonality test */
  MPI_Comm comm = MPI_COMM_WORLD;

  /* Parameters required to set up and test a SCOPFlow */
  std::string solver = SCOPFLOWOptions::solver.default_value;
  std::string model = SCOPFLOWOptions::model.default_value;
  std::string network = "";
  std::string contingencies = "";
  std::string pload = "";
  std::string qload = "";
  std::string windgen = "";
  std::string description = "";
  int num_contingencies = 0;
  double tolerance = SCOPFLOWOptions::tolerance.default_value;
  double warning_tolerance = 0.01;
  int iter_range = 0;
  double duration = SCOPFLOWOptions::duration.default_value;
  double dT = SCOPFLOWOptions::dT.default_value;
  int mode = SCOPFLOWOptions::mode.default_value;
  bool multiperiod = SCOPFLOWOptions::enable_multiperiod.default_value;

  /* Parameters used to modify underlying opflow */
  std::string opflow_initialization_string =
      OPFLOWOptions::initialization.default_value;
  int opflow_initialization;
  std::string opflow_genbusvoltage_string =
      OPFLOWOptions::genbusvoltage.default_value;
  int opflow_genbusvoltage;
  std::string subproblem_solver = OPFLOWOptions::solver.default_value;
  std::string subproblem_model = OPFLOWOptions::model.default_value;
  bool enable_powerimbalance_variables =
      OPFLOWOptions::include_powerimbalance_variables.default_value;
  bool ignore_lineflow_constraints =
      OPFLOWOptions::ignore_lineflow_constraints.default_value;
  int verbosity_level = 3;
  std::string mem_space = "DEFAULT";
  std::string compute_mode = "cpu";

  /* Parameters used to determine success or failure of functionality test */
  int expected_num_iters;
  double expected_obj_value;
  std::vector<std::string> reasons_for_failure;
  std::vector<std::string> warnings;

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
    set_if_found(warning_tolerance, values, "warning_tolerance");
    set_if_found(opflow_initialization_string, values, "opflow_initialization");
    set_if_found(opflow_genbusvoltage_string, values, "opflow_genbusvoltage");
    set_if_found(mode, values, "mode");
    set_if_found(multiperiod, values, "multiperiod");
    set_if_found(duration, values, "duration");
    set_if_found(dT, values, "dT");
    set_if_found(subproblem_solver, values, "subproblem_solver");
    set_if_found(subproblem_model, values, "subproblem_model");
    set_if_found(compute_mode, values, "compute_mode");
    set_if_found(verbosity_level, values, "verbosity_level");
    set_if_found(mem_space, values, "mem_space");
    set_if_found(iter_range, values, "iter_range");
    set_if_found(enable_powerimbalance_variables, values,
                 "enable_powerimbalance_variables");
    set_if_found(ignore_lineflow_constraints, values,
                 "ignore_lineflow_constraints");

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

  using Params = ScopflowFunctionalityTestParameters;
  MPI_Comm comm;
  int nprocs;

  ScopflowFunctionalityTests(std::string testsuite_filename, MPI_Comm comm,
                             int logging_verbosity = EXAGO_LOG_INFO)
      : FunctionalityTestContext(testsuite_filename, logging_verbosity),
        comm{comm} {
    auto err = MPI_Comm_size(comm, &nprocs);
    if (err != MPI_SUCCESS) {
      throw ExaGOError("Error getting MPI num ranks");
    }
  }

  void
  ensure_options_are_consistent(toml::value testcase,
                                toml::value presets = toml::value{}) override {

    auto ensure_option_available = [&](const std::string &opt) {
      bool is_available = testcase.contains(opt) || presets.contains(opt);
      if (is_true_somewhere(!is_available,comm)) {
        std::stringstream errs;
        errs << "SCOPFLOW Test suite expected option '" << opt
             << "' to be available, but it was not found in this testsuite"
             << " configuration:\n";
        errs << testcase << "\nwith these presets:\n" << presets;
        throw ExaGOError(errs.str().c_str());
      }
    };

    for (const auto &opt :
         {"solver", "model", "network", "contingencies", "tolerance",
          "subproblem_model", "subproblem_solver", "compute_mode",
          "verbosity_level"})
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
    testcase["expected_num_iters"] = params.numiter;
    testcase["iter_range"] = params.iter_range;
    testcase["observed_num_iters"] = params.numiter;
    testcase["obj_value"] = params.expected_obj_value;
    testcase["observed_obj_value"] = params.obj_value;
    testcase["scaled_objective_value_error"] = params.error;
    testcase["tolerance"] = params.tolerance;
    testcase["warning_tolerance"] = params.warning_tolerance;
    testcase["did_scopflow_converge"] = params.conv_status;
    testcase["mode"] = params.mode;
    testcase["subproblem_solver"] = params.subproblem_solver;
    testcase["subproblem_model"] = params.subproblem_model;
    testcase["compute_mode"] = params.compute_mode;
    testcase["verbosity_level"] = params.verbosity_level;
    testcase["enable_powerimbalance_variables"] =
        params.enable_powerimbalance_variables;
    testcase["ignore_lineflow_constraints"] =
        params.ignore_lineflow_constraints;
    testcase["reasons_for_failure"] = params.reasons_for_failure;
    testcase["warnings"] = params.warnings;

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
    int rank;
    auto err = MPI_Comm_rank(comm, &rank);
    if (err != MPI_SUCCESS)
      throw ExaGOError("Error getting MPI rank number");

    if (rank == 0)
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
    std::string ext = FileNameExtension(params.contingencies);
    resolve_datafiles_path(params.contingencies);
    if (ext == "con") {
      ierr = SCOPFLOWSetContingencyData(scopflow, PSSE,
          params.contingencies.c_str());
    } else {
      ierr = SCOPFLOWSetContingencyData(scopflow, NATIVE,
          params.contingencies.c_str());
    }
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

    ierr = SCOPFLOWSetSolver(scopflow, params.solver);
    ExaGOCheckError(ierr);

    ierr = SCOPFLOWSetInitilizationType(
        scopflow, (OPFLOWInitializationType)params.opflow_initialization);
    ExaGOCheckError(ierr);
    ierr = SCOPFLOWSetGenBusVoltageType(
        scopflow, (OPFLOWGenBusVoltageType)params.opflow_genbusvoltage);
    ExaGOCheckError(ierr);

    ierr = SCOPFLOWSetModel(scopflow, params.model);
    ExaGOCheckError(ierr);

    ierr = SCOPFLOWSetMode(scopflow, (PetscInt)params.mode);
    ExaGOCheckError(ierr);

    ierr = SCOPFLOWSetSubproblemSolver(scopflow, params.subproblem_solver);
    ExaGOCheckError(ierr);

    ierr = SCOPFLOWSetSubproblemModel(scopflow, params.subproblem_model);
    ExaGOCheckError(ierr);

    ierr = SCOPFLOWSetSubproblemComputeMode(scopflow, params.compute_mode);
    ExaGOCheckError(ierr);

    ierr = SCOPFLOWSetSubproblemVerbosityLevel(
        scopflow, (PetscInt)params.verbosity_level);
    ExaGOCheckError(ierr);

    ierr = SCOPFLOWSetSubproblemMemSpace(scopflow, params.mem_space);
    ExaGOCheckError(ierr);

    ierr = SCOPFLOWEnablePowerImbalanceVariables(
        scopflow, (PetscBool)params.enable_powerimbalance_variables);
    ExaGOCheckError(ierr);

    ierr = SCOPFLOWIgnoreLineflowConstraints(
        scopflow, (PetscBool)params.ignore_lineflow_constraints);
    ExaGOCheckError(ierr);

    ierr = SCOPFLOWSetUp(scopflow);
    ExaGOCheckError(ierr);

    ierr = SCOPFLOWSolve(scopflow);
    ExaGOCheckError(ierr);

    /* Possible ways for a funcitonality test to fail */
    bool converge_failed = false;
    bool obj_failed = false;
    bool obj_warning = false;
    bool num_iter_failed = false;
    bool num_iter_warning = false;
    params.reasons_for_failure.clear();
    params.warnings.clear();

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
      if (!IsEqual(params.obj_value, params.expected_obj_value,
                   params.warning_tolerance, params.error)) {
        obj_failed = true;
#ifdef EXAGO_ENABLE_LOGGING
        params.reasons_for_failure.push_back(
            fmt::format("expected objective value={} actual objective value={} "
                        "tol={} err={}",
                        params.expected_obj_value, params.obj_value,
                        params.tolerance, params.error));
#else
        char sbuf[256];
        sprintf(sbuf,
                "expected objective value=%e actual objective value=%e tol=%e "
                "err=%e",
                params.expected_obj_value, params.obj_value, params.tolerance,
                params.error);
        params.reasons_for_failure.push_back(std::string(sbuf));
#endif
      } else {
        obj_warning = true;
#ifdef EXAGO_ENABLE_LOGGING
        params.warnings.push_back(fmt::format(
            "expected objective value={} actual objective value={} tol={}"
            " warning_to={} err={}",
            params.expected_obj_value, params.obj_value, params.tolerance,
            params.warning_tolerance, params.error));
#else
        char sbuf[256];
        sprintf(sbuf,
                "expected objective value=%e actual objective value=%e tol=%e"
                " warning_to=%e err=%e",
                params.expected_obj_value, params.obj_value, params.tolerance,
                params.warning_tolerance, params.error);
        params.warnings.push_back(std::string(sbuf));
#endif
      }
    }

    /* Test num iterations */
    ierr = SCOPFLOWGetNumIterations(scopflow, &params.numiter);
    ExaGOCheckError(ierr);
    if (params.expected_num_iters != -1 &&
        params.numiter != params.expected_num_iters) {
      int diff = abs(params.numiter - params.expected_num_iters);
      if (diff > params.iter_range) {
        num_iter_failed = true;
#ifdef EXAGO_ENABLE_LOGGING
        params.reasons_for_failure.push_back(
            fmt::format("expected {} num iters, got {}",
                        params.expected_num_iters, params.numiter));
#else
        char sbuf[256];
        sprintf(sbuf, "expected %d num iters, got %d",
                params.expected_num_iters, params.numiter);
        params.reasons_for_failure.push_back(std::string(sbuf));
#endif
      } else {
        num_iter_warning = true;
#ifdef EXAGO_ENABLE_LOGGING
        params.warnings.push_back(fmt::format("expected {} num iters, got {}",
                                              params.expected_num_iters,
                                              params.numiter));
#else
        char sbuf[256];
        sprintf(sbuf, "expected %d num iters, got %d",
                params.expected_num_iters, params.numiter);
        params.warnings.push_back(std::string(sbuf));
#endif
      }
    }

    /* Did the current functionality test fail in any way? */
    bool local_fail = converge_failed || obj_failed || num_iter_failed;

    if (local_fail)
      fail();
    else {
      if (num_iter_warning || obj_warning) {
        warning();
      } else {
        pass();
      }
    }

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
  std::string name;
  if (argc > 1) {
    name = "Test from TOML file ";
    std::string str(argv[1]);
    int idx = str.find_last_of('/');
    str.erase(0, idx + 1);
    name = "Test from TOML file scopflow/";
    name.append(str);
  } else {
    name = "UNAMED SCOPFLOW TEST";
  }

  PetscErrorCode ierr;
  MPI_Comm comm = MPI_COMM_WORLD;
  char appname[] = "scopflow";
  ierr = ExaGOInitialize(comm, &argc, &argv, appname, help);
  ExaGOCheckError(ierr);
  ExaGOLog(EXAGO_LOG_INFO, "{}", "Creating SCOPFlow Functionality Test");

  int my_rank;
  auto err = MPI_Comm_rank(comm, &my_rank);
  if (err != MPI_SUCCESS)
    throw ExaGOError("Error getting MPI rank number");

  if (my_rank == 0)
    ExaGOLog(EXAGO_LOG_INFO, "{}", "Creating SCOPFLOW Functionality Test");
  ScopflowFunctionalityTests test{std::string(argv[1]), comm};
  test.run_all_test_cases();
  test.print_report();
  std::string filename = test.set_file_name(argv[1]);
  filename.append(".warning");
  test.print_warning(filename, name);

  ExaGOFinalize();
  return test.failures();
}
