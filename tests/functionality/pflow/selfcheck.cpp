static char help[] = "PFLOW Functionality Tests.\n\n";

#include "toml_utils.h"
#include <pflow.h>

struct PflowFunctionalityTestParameters {
  /* Communicator required to run funcitonality test */
  MPI_Comm comm = MPI_COMM_WORLD;

  /* Parameters required to set up and test a PFlow */
  std::string network = "";
  std::string description = "";

  /* Parameters used to determine success or failure of functionality test */
  int expected_num_iters;
  std::vector<std::string> reasons_for_failure;

  /* Actual values observed from the system-under-test. */
  double error;
  int numiter;
  PetscBool conv_status = PETSC_FALSE;

  /* Assign all member variables from a toml map if values are found */
  void assign_from(toml::value values) {
    set_if_found(network, values, "network");
    set_if_found(description, values, "description");
    set_if_found(expected_num_iters, values, "num_iters");
  }
};

struct PflowFunctionalityTests
    : public FunctionalityTestContext<PflowFunctionalityTestParameters> {

  MPI_Comm comm;
  int nprocs;

  PflowFunctionalityTests(std::string testsuite_filename, MPI_Comm comm,
                          int logging_verbosity = EXAGO_LOG_INFO)
      : comm{comm},
        FunctionalityTestContext(testsuite_filename, logging_verbosity) {
    auto err = MPI_Comm_size(comm, &nprocs);
    if (err)
      throw ExaGOError("Error getting MPI num ranks");
  }

  using Params = PflowFunctionalityTestParameters;
  void
  ensure_options_are_consistent(toml::value testcase,
                                toml::value presets = toml::value{}) override {

    int n_preset_procs;
    set_if_found(n_preset_procs, presets, "n_procs");
    int n_testcase_procs = -1;
    set_if_found(n_testcase_procs, testcase, "n_procs");

    if (-1 != n_testcase_procs) {
      std::stringstream errs;
      errs << "Number of processes should be declared globally in the preset "
              "area of the test suite TOML file, not inside each testcase.\n"
           << "Testcase: " << testcase << "\nWith presets:\n"
           << presets;
      throw ExaGOError(errs.str().c_str());
    } else if (nprocs != n_preset_procs) {
      std::stringstream errs;
      errs << "PFLOW Functionality test suite found " << n_preset_procs
           << " processes specified in the presets of the test suite TOML "
              "file, but this test is being run with "
           << nprocs << " processes.\nTestcase: " << testcase
           << "\nWith presets:\n"
           << presets;
      throw ExaGOError(errs.str().c_str());
    }

    auto ensure_option_available = [&](const std::string &opt) {
      bool is_available = testcase.contains(opt) || presets.contains(opt);
      if (!is_available) {
        std::stringstream errs;
        errs << "PFLOW Test suite expected option '" << opt
             << "' to be available, but it was not found in this testsuite"
             << " configuration:\n";
        errs << testcase << "\nwith these presets:\n" << presets;
        throw ExaGOError(errs.str().c_str());
      }
    };

    for (const auto &opt : {"network", "n_procs"})
      ensure_option_available(opt);
  }

  void initialize_test_parameters(Params &params, const toml::value &testcase,
                                  const toml::value &presets) override {
    params.assign_from(presets);
    params.assign_from(testcase);
  }

  toml::value create_failing_testcase(const Params &params) override {
    toml::value testcase;

    testcase["network"] = params.network;
    testcase["num_iters"] = params.expected_num_iters;
    testcase["observed_num_iters"] = params.numiter;
    testcase["did_pflow_converge"] = params.conv_status;
    testcase["reasons_for_failure"] = params.reasons_for_failure;

    return testcase;
  }

  void run_test_case(Params &params) override {
    PetscErrorCode ierr;
    PFLOW pflow;

    std::cout << "Test Description: " << params.description << std::endl;
    ierr = PFLOWCreate(params.comm, &pflow);
    ExaGOCheckError(ierr);

    // Prepend installation directory to network path
    resolve_datafiles_path(params.network);
    if (strstr(params.network.c_str(), ".raw") != NULL) {
      ierr = PFLOWReadPSSERawData(pflow, params.network.c_str());
      ExaGOCheckError(ierr);
    } else {
      ierr = PFLOWReadMatPowerData(pflow, params.network.c_str());
      ExaGOCheckError(ierr);
    }

    ierr = PFLOWSetUp(pflow);
    ExaGOCheckError(ierr);

    ierr = PFLOWSolve(pflow);
    ExaGOCheckError(ierr);

    /* Possible ways for a funcitonality test to fail */
    bool converge_failed = false;
    bool num_iter_failed = false;

    /* Test convergence status */
    ierr = PFLOWConverged(pflow, &params.conv_status);
    ExaGOCheckError(ierr);
    if (params.conv_status == PETSC_FALSE) {
      converge_failed = true;
      params.reasons_for_failure.push_back("failed to converge");
    }

    /* Test num iterations */
    ierr = PFLOWGetNumIterations(pflow, &params.numiter);
    ExaGOCheckError(ierr);
    if (params.numiter != params.expected_num_iters) {
      num_iter_failed = true;
      params.reasons_for_failure.push_back(
          fmt::format("expected {} num iters, got {}",
                      params.expected_num_iters, params.numiter));
    }

    /* Did the current functionality test fail in any way? */
    bool local_fail = converge_failed || num_iter_failed;

    if (local_fail)
      fail();
    else
      pass();

    ierr = PFLOWDestroy(&pflow);
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
  char appname[] = "pflow";
  int _argc = 0;
  ierr = ExaGOInitialize(comm, &argc, &argv, appname, help);

  ExaGOLog(EXAGO_LOG_INFO, "{}", "Creating PFlow Functionality Test");

  PflowFunctionalityTests test{std::string(argv[1]), comm};
  test.run_all_test_cases();
  test.print_report();

  ExaGOFinalize();
  return test.failures();
}
