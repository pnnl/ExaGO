static char help[] = "SCOPFLOW Functionality Tests.\n\n";

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <toml.hpp>
#include <scopflow.h>
#include <exago_config.h>
#include <utils.hpp>
#include <version.hpp>

void resolve_datafiles_path(std::string &path)
{
  std::vector<std::string> prefix;
  prefix.push_back("./");
  prefix.push_back( std::string(EXAGO_OPTIONS_DIR) + "/../");
  for (int i=0; i<prefix.size(); i++) {
    std::ifstream f{(prefix[i]+path).c_str()};
    if (f.is_open()) {
      path = prefix[i]+path;
      f.close();
      break;
    }
  }
}

/* For formatting reports as TOML */
void fmt_row(std::stringstream& summary, int col_width, std::string key,
    std::string value)
{
  std::stringstream value_fmt;
  value_fmt << "'" << value << "'";

  summary
    << std::setw(col_width-1) << std::left << key << std::right << "="
    << std::setw(col_width-1) << std::left << value_fmt.str() << std::right
    << "\n";
}

template<typename T>
void fmt_row(std::stringstream& summary, int col_width, std::string key, T value)
{
  summary
    << std::setw(col_width-1) << std::left << key << std::right << "="
    << std::setw(col_width-1) << std::left << value << std::right << "\n";
}

template<typename T>
void fmt_comment(std::stringstream& summary, int col_width, std::string key, T value)
{
  summary
    << "#"
    << std::setw(col_width-2) << std::left << key << std::right << "="
    << std::setw(col_width-1) << std::left << value << std::right << "\n";
}

static bool is_enabled(std::string dependency)
{
  const auto& deps = ExaGOGetDependencies();
  auto it = deps.find(dependency);
  if (it == deps.end())
    return false;
  return true;
}

static std::string bool2str(bool b)
{
  return (b?"true":"false");
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
  int num_contingencies;
  double tolerance;
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
      solver = toml::find_or(presets, "solver", "not set");
      model = toml::find_or(presets, "model", "not set");
      network = toml::find_or(presets, "network", "not set");
      contingencies = toml::find_or(presets, "contingencies", "not set");
      num_contingencies = toml::find_or(presets, "num_contingencies", 0);
      expected_num_iters = toml::find_or(presets, "num_iters", 0);
      expected_obj_value = toml::find_or(presets, "obj_value", 0.0);
      tolerance = toml::find_or(presets, "tolerance", 0.0);
      opflow_initialization = toml::find_or(presets, "opflow_intitialization", 0);
      opflow_genbusvoltage = toml::find_or(presets, "opflow_genbusvoltage", 0);
      multiperiod = toml::find_or(presets, "multiperiod", false);
    }

    /* If values are found for an individual test case, let these overwrite
     * the global values */
    solver = toml::find_or(testcase, "solver", solver);
    model = toml::find_or(testcase, "model", model);
    network = toml::find_or(testcase, "network", network);
    contingencies = toml::find_or(testcase, "contingencies", contingencies);
    num_contingencies = toml::find_or(testcase, "num_contingencies", num_contingencies);
    expected_num_iters = toml::find_or(testcase, "num_iters", expected_num_iters);
    expected_obj_value = toml::find_or(testcase, "obj_value", expected_obj_value);
    tolerance = toml::find_or(testcase, "tolerance", tolerance);
    opflow_initialization = toml::find_or(testcase, "opflow_intitialization", opflow_initialization);
    opflow_genbusvoltage = toml::find_or(testcase, "opflow_genbusvoltage", opflow_genbusvoltage);
    multiperiod = toml::find_or(testcase, "multiperiod", multiperiod);

    SCOPFLOW scopflow;
    ierr = SCOPFLOWCreate(comm,&scopflow);CHKERRQ(ierr);

    ierr = SCOPFLOWSetTolerance(scopflow, tolerance);CHKERRQ(ierr);

    // Prepend installation directory to network path
    resolve_datafiles_path(network);
    ierr = SCOPFLOWSetNetworkData(scopflow,network.c_str());CHKERRQ(ierr);

    // Prepend installation directory to contingency file
    resolve_datafiles_path(contingencies);
    ierr = SCOPFLOWSetContingencyData(scopflow, NATIVE, contingencies.c_str());CHKERRQ(ierr);

    // Set number of contingencies
    ierr = SCOPFLOWSetNumContingencies(scopflow, num_contingencies);CHKERRQ(ierr);

    ierr = SCOPFLOWEnableMultiPeriod(scopflow, (PetscBool)multiperiod);CHKERRQ(ierr);

    ierr = SCOPFLOWSetSolver(scopflow, solver.c_str());CHKERRQ(ierr);

    ierr = SCOPFLOWSetInitilizationType(scopflow, (OPFLOWInitializationType)opflow_initialization);CHKERRQ(ierr);
    ierr = SCOPFLOWSetGenBusVoltageType(scopflow, (OPFLOWGenBusVoltageType)opflow_genbusvoltage);CHKERRQ(ierr);

    ierr = SCOPFLOWSetModel(scopflow, model.c_str());CHKERRQ(ierr);

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
      fmt_row(summary, col_width, "contingencies", contingencies);
      fmt_row(summary, col_width, "num_contingencies", num_contingencies);
      fmt_row(summary, col_width, "opflow_intitialization", opflow_initialization);
      fmt_row(summary, col_width, "opflow_genbusvoltage", opflow_genbusvoltage);
      fmt_row(summary, col_width, "multiperiod", bool2str(multiperiod));
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
    fmt_row(summary, col_width, "testsuite_name",
        "ExaGO Test Suite Automatically Generated from Failed Tests");
    ExaGOLog(EXAGO_LOG_ERROR, "%s", summary.str().c_str());
  }
  ExaGOFinalize();
  return fail;
}
