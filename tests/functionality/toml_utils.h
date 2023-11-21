#pragma once
#ifdef EXAGO_ENABLE_LOGGING
#include <spdlog/fmt/fmt.h>
#endif
#include <exago_config.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <toml.hpp>
#include <utils.h>
#include <version.h>

/**
 *
 * @brief Utilities for working with TOML files which describe test suites for
 * ExaGO's functionality tests.
 *
 * @author Asher Mancinelli <asher.mancinelli@pnnl.gov>
 *
 */

template <typename ValueType>
void set_if_found(ValueType &value, toml::value config,
                  const std::string &key) {
  if (config.contains(key) and config.count(key) > 0) {
    value = toml::find<ValueType>(config, key);
  }
}

void resolve_datafiles_path(std::string &path) {
  std::vector<std::string> prefixes;
  prefixes.push_back("./");
  prefixes.push_back(std::string(EXAGO_OPTIONS_DIR) + "/../");
  for (const auto &prefix : prefixes) {
    std::ifstream f{prefix + path};
    if (f.is_open()) {
      path = prefix + path;
      f.close();
      return;
    }
  }
}

/* For formatting reports as TOML */
void fmt_row(std::ostream &summary, int col_width, std::string key,
             std::string value) {
  std::stringstream value_fmt;
  value_fmt << "'" << value << "'";

  summary << std::setw(col_width - 1) << std::left << key << std::right << "="
          << std::setw(col_width - 1) << std::left << value_fmt.str()
          << std::right << "\n";
}

template <typename T>
void fmt_row(std::ostream &summary, int col_width, std::string key, T value) {
  summary << std::setw(col_width - 1) << std::left << key << std::right << "="
          << std::setw(col_width - 1) << std::left << value << std::right
          << "\n";
}

template <typename T>
void fmt_comment(std::ostream &summary, int col_width, std::string key,
                 T value) {
  summary << "#" << std::setw(col_width - 2) << std::left << key << std::right
          << "=" << std::setw(col_width - 1) << std::left << value << std::right
          << "\n";
}

bool is_true_somewhere(bool flag, MPI_Comm comm) {
  bool ret;
  int err = MPI_Allreduce(&flag, &ret, 1, MPI_CXX_BOOL, MPI_LOR, comm);
  if (err != MPI_SUCCESS) {
    throw ExaGOError("Error in is_true_somewhere for MPI_Allreduce");
  }
  return ret;
}

// static std::string bool2str(bool b) { return (b ? "true" : "false"); }

/**
 *  Context manager for functionality tests for ExaGO applications.
 *
 * Each test context has an internal struct to encapsulate the test parameters.
 * These parameters will be passed to each test runner. An example type for the
 * test parameters might be the following:
 *
 * struct TestParameterType
 * {
 *   std::string solver, model, network;
 *   std::size_t expected_num_iters;
 *   double expected_obj_value;
 * };
 *
 * A concrete structure for a functionality test would then be defined like so:
 *
 * struct ConcreteFunctionalityTest
 *   : FunctionalityTestContext<TestParameterType>
 * {
 *   // ...
 * };
 *
 */
template <typename TestParameters> struct FunctionalityTestContext {
private:
  TestParameters test_parameters_;
  inline TestParameters &test_parameters() { return test_parameters_; }

public:
  /* Column width used when formatting failing test suite */
  static constexpr int col_width = 35;

  FunctionalityTestContext(std::string testsuite_filename, MPI_Comm comm,
                           int logging_verbosity = EXAGO_LOG_INFO)
      : logging_verbosity_{logging_verbosity}, comm_{comm} {
    testsuite_ = toml::parse(testsuite_filename);
    testcases_ = toml::find(testsuite(), "testcase");
    testsuite_name_ = toml::find<std::string>(testsuite(), "testsuite_name");
    auto err = MPI_Comm_rank(comm, &rank);
    if (err != MPI_SUCCESS)
      throw ExaGOError("Error getting MPI rank number");
    err = MPI_Comm_size(comm, &nprocs);
    if (err != MPI_SUCCESS) {
      throw ExaGOError("Error getting MPI num ranks");
    }
  }

  void run_all_test_cases() {
    for (auto &testcase : testcases().as_array()) {
      toml::value presets{};
      set_if_found(presets, testsuite(), "presets");
      ensure_options_are_consistent(testcase, presets);
      initialize_test_parameters(test_parameters(), testcase, presets);
      run_test_case(test_parameters());
    }
  }

  std::string set_file_name(char *name) {
    std::string str(name);
    int idx = str.find_last_of('/');
    int len = str.length();
    str.erase(0, idx + 1);
    idx = str.find_last_of('.');
    len = str.length();
    str.erase(idx, len - idx);
    return str;
  }

  /** For each testcase in a given test suite, ensure the options are internally
   * consistent. Throw if invalid - we don't want to attempt to keep running
   * if the test suite is ill-formed.
   */
  virtual void ensure_options_are_consistent(toml::value testcase,
                                             toml::value presets) = 0;

  /**
   * Initialize test parameters for a given testcase and set of presets if they
   * exist.
   */
  virtual void initialize_test_parameters(TestParameters &test_parameters,
                                          const toml::value &testcase,
                                          const toml::value &presets) = 0;

  /**
   * Create summary of configuration. This toml value should be a valid
   * testcase, and is allowed to have extra unused values which hint to the
   * user what went wrong in the functionality test run.
   */
  virtual toml::value create_failing_testcase(const TestParameters &) = 0;

  /* Callback for each failing test */
  virtual inline void fail() {
    ExaGOLog(verbosity(), "{}", "-- FAIL");
    auto testcase = create_failing_testcase(test_parameters());
    failing_testcases_.push_back(testcase);
    failures_++;
    total_num_tests_++;
  }

  virtual inline void warning() {
    ExaGOLog(verbosity(), "{}", "-- WARNING");
    auto testcase = create_failing_testcase(test_parameters());
    warning_testcases_.push_back(testcase);
    warnings_++;
    total_num_tests_++;
  }

  /* Callback for each passing test */
  virtual inline void pass() {
    total_num_tests_++;
    ExaGOLog(verbosity(), "{}", "-- PASS");
  }

  void print_report() {
    if (rank != 0)
      return;
    ExaGOLog(verbosity(), "{:d} / {:d} tests failed. {:d} warnings \n",
             failures(), total_tests(), warnings());
    if (failures()) {
      ExaGOLog(verbosity(), "{}", "Summary of failing functionality tests:");
      std::cout.precision(12);
      std::cout << std::flush
                << "###########################################################"
                   "#####################\n"
                << "# Begin auto-generated TOML test suite\n"
                << failing_testsuite()
                << "# End auto-generated TOML test suite\n"
                << "###########################################################"
                   "#####################\n"
                << "\n";
    }
  }

  void print_warning(std::string filename, std::string name) {
    if (warnings()) {
      warning_stream.precision(12);
      warning_stream
          << std::flush
          << "[ExaGO] Summary of functionality tests with warnings for\n"
          << name << ":\n"
          << "###########################################################"
             "#####################\n"
          << "# Begin auto-generated TOML test suite\n"
          << warning_testsuite() << "# End auto-generated TOML test suite\n"
          << "###########################################################"
             "#####################\n"
          << "\n";
      // extract final filename from string
      std::string file = filename;
      filename = ("../");
      filename.append(file);
      if (rank == 0) {
        std::ofstream fout;
        fout.open(filename.c_str());
        fout << warning_stream.str() << std::endl;
        fout.close();
      }
    }
  }

  /* Return number of failures encountered thus far */
  inline const std::size_t &failures() const { return failures_; }
  inline const std::size_t &warnings() const { return warnings_; }
  inline const std::size_t &total_tests() const { return total_num_tests_; }

private:
  /* Execute a single test case in a test suite. Test case runs should use
   * the `fail` and `pass` methods to report the status of each test. */
  virtual void run_test_case(TestParameters &test_parameters) = 0;

  /* Output stream to write the test suite of failures to */
  toml::value failing_testsuite() {
    toml::table testsuite;
    std::stringstream desc;
    desc << "Auto-generated test suite based on testsuite with name '"
         << testsuite_name() << "'";
    testsuite["testsuite_name"] = desc.str();
    testsuite["testcase"] = failing_testcases_;
    return testsuite;
  }
  /* Output stream to write the test suite of warnings to */
  toml::value warning_testsuite() {
    toml::table testsuite;
    std::stringstream desc;
    desc << "Auto-generated test suite based on testsuite with name '"
         << testsuite_name() << "'";
    testsuite["testsuite_name"] = desc.str();
    testsuite["testcase"] = warning_testcases_;
    return testsuite;
  }

  /* Verbosity used by calls to ExaGOLog in methods of this structure. */
  inline int verbosity() const { return logging_verbosity_; }

  /* Const accessors for constant private members */
  inline const toml::value &testcases() const { return testcases_; }
  inline const toml::value &testsuite() const { return testsuite_; }
  inline const std::string &testsuite_name() const { return testsuite_name_; }

private:
  toml::value failing_testsuite_;
  std::vector<toml::value> failing_testcases_;
  toml::value warning_testsuite_;
  std::vector<toml::value> warning_testcases_;
  std::size_t failures_{0}, warnings_{0}, total_num_tests_{0};
  int logging_verbosity_;
  std::string testsuite_name_{""};
  toml::value testcases_;
  toml::value testsuite_;
  std::stringstream warning_stream;
  MPI_Comm comm_;

protected:
  int rank;
  int nprocs;
};
