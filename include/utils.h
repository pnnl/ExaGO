#ifndef EXAGO_UTILS_H
#define EXAGO_UTILS_H
#include "exago_config.h"
#include <petsc.h>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#ifdef EXAGO_ENABLE_LOGGING
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/fmt.h>
#include <spdlog/fmt/bundled/color.h>
#endif


/**
 *
 * @file utils.hpp
 *
 * The key APIs this header exposes are as follows:
 *
 * ## ExaGOError
 *
 * This error class may be constructed from an error code, a message, or nothing
 * at all (discouraged, as it's ambiguous). If constructed from an error code,
 * ExaGO will assume the error code comes from PETSc, and will query PETSc to
 * determine the error message. This error message will then be propagated
 * upwards.
 *
 * ## ExaGOLog
 *
 * ExaGO log messages are conditionally logged depending on the minimum
 * verbosity level.
 *
 * ExaGO verbosity levels are between [0, 10], however setting the minimum
 * verbosity level to values outside this range is allowed. For example, to log
 * all log messages, the user may set the log level to -1. Or, to disable all
 * logging at runtime, the user may set the log level to 99 or any value above
 * 10.
 *
 * The utility enum ExaGOLogLevel is provided to make calls to ExaGOLog more
 * readable, however passing ints directly is perfectly allowed as well.
 *
 * The minimum verbosity level a message must be to be logged may be set/queried
 * with the ExaGOLogSetMinLogLevel/ExaGOLogGetMinLogLevel functions.
 *
 * This logging facility will be entirely disabled at build-time if the CMake
 * variable EXAGO_ENABLE_LOGGING is set to OFF at configure-time.
 *
 */

enum ExaGOLogLevel : signed int {
  EXAGO_LOG_INFO = 0,
  EXAGO_LOG_WARN = 5,
  EXAGO_LOG_ERROR = 10,
};

/**
 * @brief ExaGO Error interfaces between error codes encountered in ExaGO
 * proper as well as error codes recieved from PETSc.
 */
struct ExaGOError : public std::exception {
  ExaGOError() : message{"ExaGO Error"} {}
  ExaGOError(const char *message) : message{message} {}
  ExaGOError(PetscErrorCode);

  /* The name _what_ is not in PascalCase like the rest of ExaGO because
   * ExaGOError inherits from the standard library exception which defines
   * _what_. */
  virtual const char *what() const noexcept { return message.c_str(); };
  virtual bool IsPetscError() const noexcept { return is_petsc_error; }

protected:
  std::string message;
  bool is_petsc_error = false;
};

/* Used to interface Petsc error return codes with ExaGO errors */
extern void ExaGOCheckError(int e);

/* Get the name of ExaGO's logger within spdlog */
extern PetscErrorCode ExaGOLogGetLoggerName(std::string &s);

/** Get minimum loglevel for a log to be printed */
extern PetscErrorCode ExaGOLogGetMinLogLevel(int &);

/** Set minimum loglevel for a log to be printed */
extern PetscErrorCode ExaGOLogSetMinLogLevel(int);

/**
 * @brief Implementation to log string according to ExaGO build configuration.
 * @note To log on every rank, you may use the overload which does not take a
 * communicator.
 */
template <typename... Args>
inline void ExaGOLog(MPI_Comm comm, int level, std::string fmt, Args... args) {
#ifdef EXAGO_ENABLE_LOGGING

  /* Check that the rank is 0 before logging */
  int rank;
  int ierr = MPI_Comm_rank(comm, &rank);
  ExaGOCheckError(ierr);
  if (0 != rank) {
    return;
  }

  /* Check that the current log level is greater than or equal to the current
   * minimum required log level */
  int loglevel;
  ierr = ExaGOLogGetMinLogLevel(loglevel);
  ExaGOCheckError(ierr);
  if (level >= loglevel) {

    /* Perform the actual logging */
    std::string logname;
    ierr = ExaGOLogGetLoggerName(logname);
    ExaGOCheckError(ierr);
    auto logger = spdlog::get(logname);

    /*
     * Because we handler our own verbosity levels, we just use `info` for all
     * log messages.
     */
    logger->info(fmt, args...);
  }
#endif
}

/**
 * @brief Overload which does not take an MPI communicator and will log on every
 * rank.
 */
template <typename... Args>
void inline ExaGOLog(int level, std::string fmt, Args... args) {
#ifdef EXAGO_ENABLE_LOGGING
  ExaGOLog(PETSC_COMM_SELF, level, fmt, args...);
#endif
}

template <typename T> struct ExaGOOption {
  std::string opt;
  std::string desc;
  T default_value;
  ExaGOOption(std::string const &opt, std::string const &desc, T const &value)
      : opt{opt}, desc{desc}, default_value{value} {}
};

template <> struct ExaGOOption<void> {
  std::string opt;
  std::string desc;
  ExaGOOption(std::string const &opt, std::string const &desc)
      : opt{opt}, desc{desc} {}
};

template <> struct ExaGOOption<std::string> {
  std::string default_value;
  std::string desc;
  std::string opt;
  std::vector<std::string> possible_values;

  ExaGOOption(std::string const &opt, std::string const &desc,
              std::string const &default_value,
              std::vector<std::string> const &possible_values)
      : opt{opt}, desc{desc}, possible_values{possible_values},
        default_value{default_value} {}

  /**
   * \brief Get the enum value of a stringy enum option.
   *
   * \param names array of names to be searched for the default value
   * \param num_items number of items in `names` to be searched
   * \param prefix prefix to add to elements of `names`. PETSc enum options take
   *        a similar argument, so we handle it here to minimize the number of
   *        changes required to effectively handle command line options.
   */
  std::size_t ToEnum(const char *const names[], std::size_t num_items,
                     std::string const &prefix = "") const {
    const auto it =
        std::find_if(names, names + num_items, [=](const char *name) {
          return (prefix + std::string(name)) == default_value;
        });
    std::size_t enumval = std::distance(names, it);
    if (enumval == num_items) {
      throw ExaGOError((std::string{"could not find default value '"} +
                        default_value +
                        "' when attempting to convert string option '" + opt +
                        "' into an enum value")
                           .c_str());
    }
    return enumval;
  }
};

using ExaGOStringOption = ExaGOOption<std::string>;
using ExaGOBoolOption = ExaGOOption<PetscBool>;
using ExaGOIntOption = ExaGOOption<int>;
using ExaGORealOption = ExaGOOption<double>;
using ExaGOFlagOption = ExaGOOption<void>;

template <typename T>
inline std::string ExaGOFormatOption(ExaGOOption<T> const &opt, std::size_t indent = 0,
                              std::string indentstr = "\t") {
  using std::is_same;
  using U =
      typename std::remove_reference<typename std::remove_cv<T>::type>::type;
  std::string typestr =
      (is_same<U, bool>::value or is_same<U, PetscBool>::value)
          ? "bool"
          : is_same<U, double>::value
                ? "real"
                : is_same<U, int>::value ? "int" : "unknown_type";
  std::string tab = "";
  for (int i = 0; i < indent; i++)
    tab += "\t";
#ifdef EXAGO_ENABLE_LOGGING
  return fmt::format("{0}{2} {3}\n{0}{1}{4} (type: {5})\n", tab, indentstr,
                     opt.opt,
                     fmt::format(fmt::emphasis::bold, "{}", opt.default_value),
                     opt.desc, typestr);
#else
  char sbuf[256];
  sprintf(sbuf, "%s%s\n%s%s%s (type: %s)\n", tab.c_str(), opt.opt.c_str(),
          tab.c_str(), indentstr.c_str(), opt.desc.c_str(), typestr.c_str());
  return std::string(sbuf);
#endif
}

template <>
inline std::string ExaGOFormatOption(ExaGOOption<void> const &opt, std::size_t indent,
                              std::string indentstr) {
  std::string tab = "";
  for (int i = 0; i < indent; i++)
    tab += "\t";
#ifdef EXAGO_ENABLE_LOGGING
  return fmt::format("{0}{2}\n{0}{1}{3} (type: flag)\n", tab, indentstr,
                     opt.opt, opt.desc);
#else
  char sbuf[256];
  sprintf(sbuf, "%s%s\n%s%s%s (type: flag)\n", tab.c_str(), opt.opt.c_str(),
          tab.c_str(), indentstr.c_str(), opt.desc.c_str());
  return std::string(sbuf);
#endif
}

template <>
inline std::string ExaGOFormatOption(ExaGOStringOption const &opt, std::size_t indent,
                              std::string indentstr) {
  std::string tab = "";
  for (int i = 0; i < indent; i++)
    tab += "\t";

  /* If there are no other possible values but just a single default value,
   * we won't try to use parenthesis. Parens only make sense when we're
   * demonstrating a range of possible values */
  std::string values = (opt.possible_values.size() ? "(" : "");

  /* If the argument is optional, we have a default value which we will
   * embolden to make clear */
#ifdef EXAGO_ENABLE_LOGGING
  values += fmt::format(fmt::emphasis::bold, "{}{}", opt.default_value,
                        (opt.possible_values.size() ? "|" : ""));
#else
  char sbuf[256];
  sprintf(sbuf, "%s%s", opt.default_value.c_str(),
          (opt.possible_values.size() ? "|" : ""));
  values += std::string(sbuf);
#endif
  auto it = opt.possible_values.begin();
  while (it != opt.possible_values.end()) {
    values += *it;
    if (++it != opt.possible_values.end())
      values += "|";
  }
  values += (opt.possible_values.size() ? ")" : "");
#ifdef EXAGO_ENABLE_LOGGING
  return fmt::format("{0}{2} {3}\n{1}{0}{4} (type: string)\n", tab, indentstr,
                     opt.opt, values, opt.desc);
#else
  sprintf(sbuf, "%s%s\n%s%s%s (type: flag)\n", tab.c_str(), opt.opt.c_str(),
          indentstr.c_str(), tab.c_str(), opt.desc.c_str());
  return std::string(sbuf);
#endif
}

template <typename T> void print(T const &opt) {
  static constexpr std::size_t indent = 1;
  static const std::string indentstr = "\t";
#ifdef EXAGO_ENABLE_LOGGING
  fmt::print(ExaGOFormatOption(opt, indent, indentstr));
  fmt::print("\n");
#else
  printf("%s", ExaGOFormatOption(opt, indent, indentstr).c_str());
  printf("\n");
#endif
};

/*
 * printHelp - Global flag that tells applications to print help menu and exit
 */
extern PetscBool printHelp;

/**
 * Initialize an ExaGO application.
 *
 * @note this takes care of Petsc initialization, so don't this function in
 * conjunction with `PetscInitialize.`
 */
extern "C" PetscErrorCode ExaGOInitialize(MPI_Comm, int *argc, char ***argv,
                                          char *appname, char *help);

/**
 * Teardown for an ExaGO application.
 *
 * @note this takes care of Petsc finalization, so don't this function in
 * conjunction with `PetscFinalize`.
 */
extern "C" PetscErrorCode ExaGOFinalize();

/** Returns 1 if files exists, else 0 */
bool DoesFileExist(const char *);

/** Returns 1 if directory exists, else 0 */
bool DoesDirExist(const char *);

/** First path in _files_ to be statable as a regular file */
std::vector<std::string>::const_iterator
FirstExistingFile(const std::vector<std::string> &files);

/** Find file name extension */
std::string FileNameExtension(const std::string &filename);

inline bool IsEqual(double value, double ref, double tolerance, double &error) {
  error = std::fabs(value - ref) / (1. + std::fabs(ref));
  return (error < tolerance);
}

inline bool IsEqual(int value, int ref, int tolerance, int &error) {
  error = std::abs(value - ref);
  return (error < tolerance);
}

/** Used by python interface to get self communicattor */
PETSC_EXTERN PetscErrorCode ExaGOGetSelfCommunicator(MPI_Comm *);

#endif
