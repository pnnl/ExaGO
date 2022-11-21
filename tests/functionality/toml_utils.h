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

static std::string bool2str(bool b) { return (b ? "true" : "false"); };

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

  FunctionalityTestContext(std::string testsuite_filename,
                           int logging_verbosity = EXAGO_LOG_INFO)
      : logging_verbosity_{logging_verbosity} {
    testsuite_ = toml::parse(testsuite_filename);
    testcases_ = toml::find(testsuite(), "testcase");
    testsuite_name_ = toml::find<std::string>(testsuite(), "testsuite_name");
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

  /* Callback for each passing test */
  virtual inline void pass() {
    total_num_tests_++;
    ExaGOLog(verbosity(), "{}", "-- PASS");
  }

  void print_report() {
    ExaGOLog(verbosity(), "{:d} / {:d} tests failed.\n", failures(),
             total_tests());
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

  /* Return number of failures encountered thus far */
  inline const std::size_t &failures() const { return failures_; }
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

  /* Verbosity used by calls to ExaGOLog in methods of this structure. */
  inline int verbosity() const { return logging_verbosity_; }

  /* Const accessors for constant private members */
  inline const toml::value &testcases() const { return testcases_; }
  inline const toml::value &testsuite() const { return testsuite_; }
  inline const std::string &testsuite_name() const { return testsuite_name_; }

private:
  toml::value failing_testsuite_;
  std::vector<toml::value> failing_testcases_;
  std::size_t failures_{0}, total_num_tests_{0};
  int logging_verbosity_;
  std::string testsuite_name_{""};
  toml::value testcases_;
  toml::value testsuite_;
};
