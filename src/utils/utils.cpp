#include <common.h>
#include <dirent.h>
#include <exago_config.h>
#include <exago_build_options.h>
#include <iostream>
#include <string>
#include <regex>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/fmt.h>
#include <spdlog/fmt/bundled/color.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <sys/utsname.h>
#include <version.h>
#include <utils.h>
#if (EXAGO_ENABLE_IPOPT || EXAGO_ENABLE_HIOP)
#include <opflow.h>
#include <sopflow.h>
#include <scopflow.h>
#endif

namespace ExaGODefaultOptions {
const auto help = ExaGOFlagOption("-help", "Print help message");
const auto version = ExaGOFlagOption("-version", "Print version information");
const auto config = ExaGOFlagOption(
    "-config", "Print configuration options used to build ExaGO");
const auto options_file = ExaGOStringOption(
    "-options_file",
    "Path to options file used to load additional ExaGO configuration options",
    "/path/to/options_file", {});

/* Right now, these are handled in the application driver. If we want these
 * printed with the rest of the application help options, we'll have to move
 * the network file and print/save options into the applciation structures
 * and handle them in the solve/finalize functions. For now, they will be
 * omitted from the help message. */
// const auto network = ExaGOStringOption("-netfile", "Path to network file",
//                                        "/path/to/netfile", {});
} // namespace ExaGODefaultOptions

namespace {

/*
 * ExaGOOption formatting rationale:
 *
 * 0: the indentation level passed to the function
 * 1: indentation string, so the description is one level above the option
 * 2: the option name
 * 3: the default value or a regex-looking string describing the possible
 *    options, with any default value bolded
 * 4: the description
 * 5: the type of the option
 */
template <typename T>
std::string ExaGOFormatOption(ExaGOOption<T> const &opt, std::size_t indent = 0,
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
  return fmt::format("{0}{2} {3}\n{0}{1}{4} (type: {5})\n", tab, indentstr,
                     opt.opt,
                     fmt::format(fmt::emphasis::bold, "{}", opt.default_value),
                     opt.desc, typestr);
}

template <>
std::string ExaGOFormatOption(ExaGOOption<void> const &opt, std::size_t indent,
                              std::string indentstr) {
  std::string tab = "";
  for (int i = 0; i < indent; i++)
    tab += "\t";
  return fmt::format("{0}{2}\n{0}{1}{3} (type: flag)\n", tab, indentstr,
                     opt.opt, opt.desc);
}

template <>
std::string ExaGOFormatOption(ExaGOStringOption const &opt, std::size_t indent,
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
  values += fmt::format(fmt::emphasis::bold, "{}{}", opt.default_value,
                        (opt.possible_values.size() ? "|" : ""));
  auto it = opt.possible_values.begin();
  while (it != opt.possible_values.end()) {
    values += *it;
    if (++it != opt.possible_values.end())
      values += "|";
  }
  values += (opt.possible_values.size() ? ")" : "");
  return fmt::format("{0}{2} {3}\n{1}{0}{4} (type: string)\n", tab, indentstr,
                     opt.opt, values, opt.desc);
}

} // namespace

static char ExaGOCurrentAppName[128];

ExaGOError::ExaGOError(PetscErrorCode ierr) : is_petsc_error{true} {
  const char *error_message;
  char *specific_error_message;
  PetscErrorMessage(ierr, &error_message, &specific_error_message);
  message = error_message;
  message += specific_error_message;
}

void ExaGOCheckError(int e) {
  if (static_cast<int>(e) > PETSC_ERR_MIN_VALUE &&
      static_cast<int>(e) < PETSC_ERR_MAX_VALUE)
    throw ExaGOError(e);
}

static std::string ExaGOLoggerFileName;
static bool ExaGOLogUseFile = false;
static bool ExaGOLogIsLoggerInitialized = false;
static const std::string ExaGOLoggerName = "exago_logger";

static int ExaGOLogMinLogLevel = 0;

/** ExaGO logging communicator to determine rank*/
static MPI_Comm ExaGOLogComm = MPI_COMM_SELF;

PetscErrorCode ExaGOLogIsUsingLogFile(bool *flg) {
  *flg = ExaGOLogUseFile;
  return 0;
}

PetscErrorCode ExaGOLogGetLoggerName(std::string &s) {
  s = ExaGOLoggerName;
  return 0;
}

PetscErrorCode ExaGOGetLoggerPointer(std::shared_ptr<spdlog::logger> &logger) {
  logger = spdlog::get(ExaGOLoggerName);
  return 0;
}

/** Set ExaGO logging communicator */
PetscErrorCode ExaGOLogSetComm(MPI_Comm c) {
  ExaGOLogComm = c;
  return 0;
}

PetscErrorCode ExaGOLogInitialize() {
  bool use_file = false;
  PetscErrorCode ierr = 0;

  ierr = ExaGOLogIsUsingLogFile(&use_file);

  /** Normally we would check the error, but in here we don't even have logging
   * set up to handle the error so we have to just return a code. */
  if (ierr != 0)
    return 1;

  if (use_file) {
    std::string filename;
    ierr = ExaGOLogGetLoggerName(filename);
    /** We still don't have proper logging set up, so just return the error */
    if (ierr != 0)
      return 1;
    auto logger = spdlog::basic_logger_mt(ExaGOLoggerName, filename);
    logger->set_pattern("%^[ExaGO]%$ %v");
  } else {
    auto logger = spdlog::stdout_color_mt(ExaGOLoggerName);
    logger->set_pattern("%^[ExaGO]%$ %v");
  }
  ExaGOLogIsLoggerInitialized = true;
  return 0;
}

PetscErrorCode ExaGOLogGetMinLogLevel(int &l) {
  l = ExaGOLogMinLogLevel;
  return 0;
}

PetscErrorCode ExaGOLogSetMinLogLevel(int l) {
  ExaGOLogMinLogLevel = l;
  return 0;
}

/**
 * @brief Checks whether or not a directory entry can be opened at the given
 * path using only POSIX standard C system calls
 *
 * @param[in] pth filepath to stat
 * @see stat
 **/
bool DoesFileExist(const char *pth) {
  struct stat path_stat;
  stat(pth, &path_stat);
  return S_ISREG(path_stat.st_mode);
}

/**
 * @brief Checks whether or not a directory entry can be opened at the given
 * path using only POSIX standard C system calls.
 *
 * @param[in] path path to verify is statable
 * @return 0 if pth==nullptr or pth cannot be opened as directory with syscall
 *         1 else
 *
 * @see opendir
 **/
bool DoesDirExist(const char *pth) {
  if (pth == nullptr)
    return false;
  DIR *dp;
  dp = opendir(pth);
  if (dp == nullptr)
    return false;
  return true;
}

/**
 * @brief Finds first path to be statable as a regular file.
 *
 * @param[in] files Iterable collection of files to be searched
 * @return first iterator _it_ to point to an element of _files_ such that _it_
 * is statable as a regular file.
 *
 * @see DoesFileExist
 **/
std::vector<std::string>::const_iterator
FirstExistingFile(const std::vector<std::string> &files) {
  auto it = files.begin();
  for (; it != files.end(); it++) {
    auto message = std::string{"-- Checking "} + *it + " exists: ";
    if (DoesFileExist(it->c_str())) {
      message += "yes";
      ExaGOLog(EXAGO_LOG_INFO, "{}", message.c_str());
      return it;
    }
    message += "no";
    ExaGOLog(EXAGO_LOG_INFO, "{}", message.c_str());
  }

  /* No files were found to exist! Just return iterator to files.end() to follow
   * conventions in the STL. */
  return it;
}

/**
 * @brief Finds and returns a file name extension
 *
 * @param[filename] string containing file name
 *
 * @return the part of @c filename after the last ".", if found, empty
 * string otherwise
 *
 **/
std::string FileNameExtension(const std::string &filename) {
  std::string ext;
  std::smatch m;
  bool found = std::regex_search(
      filename, m,
      std::regex("\\.\\([^.]*\\)$", std::regex::icase | std::regex::grep));
  if (found) {
    ext = m.str(1);
  } else {
    ext = "";
  }

  return ext;
}

namespace {
template <typename T> static inline void print(T const &opt) {
  static constexpr std::size_t indent = 1;
  static const std::string indentstr = "\t";
  fmt::print(ExaGOFormatOption(opt, indent, indentstr));
  fmt::print("\n");
};
} // namespace

/**
 * @brief Print ExaGO command-line options.
 *
 * @pre Any explicitly managed resources are cleaned up before this is called -
 * this function will std::exit for you.
 *
 * @note the return type is PetscErrorCode even though ExaGOHelpPrintf is
 *       [[noreturn]]. This is to conform to PetscHelpPrintf, so ExaGOHelpPrintf
 *       may be set to the default help function pointer in Petsc.
 */
[[noreturn]] PetscErrorCode ExaGOHelpPrintf(MPI_Comm comm,
                                            const char _appname[], ...) {
  PetscErrorCode ierr;

  std::string versionstr;
  ExaGOVersionGetVersionStr(versionstr);
  fmt::print("ExaGO {} built on {}\n", versionstr, __DATE__);

  const std::string appname = _appname;
  fmt::print("\nArguments for {}:\n", appname);
  print(ExaGODefaultOptions::help);
  print(ExaGODefaultOptions::version);
  print(ExaGODefaultOptions::config);
  print(ExaGODefaultOptions::options_file);

  // See comment at ExaGODefaultOptions
  // print(O::network);

#if (EXAGO_ENABLE_IPOPT || EXAGO_ENABLE_HIOP)
  if (appname == "opflow") {
    print(OPFLOWOptions::model);
    print(OPFLOWOptions::solver);
    print(OPFLOWOptions::initialization);
    print(OPFLOWOptions::objective);
    print(OPFLOWOptions::genbusvoltage);
    print(OPFLOWOptions::has_gensetpoint);
    print(OPFLOWOptions::use_agc);
    print(OPFLOWOptions::tolerance);
    print(OPFLOWOptions::ignore_lineflow_constraints);
    print(OPFLOWOptions::include_loadloss_variables);
    print(OPFLOWOptions::loadloss_penalty);
    print(OPFLOWOptions::include_powerimbalance_variables);
    print(OPFLOWOptions::powerimbalance_penalty);

#ifdef EXAGO_ENABLE_HIOP
    print(OPFLOWOptions::hiop_compute_mode);
    print(OPFLOWOptions::hiop_verbosity_level);

#ifdef EXAGO_ENABLE_IPOPT
    print(OPFLOWOptions::hiop_ipopt_debug);
#endif

#endif

    // These options are temporarily ommitted. See comment under
    // ExaGODefaultOptions.
    //
    // fprintf(stderr, "\t -netfile ...\n");
    // fprintf(stderr, "\t -print_output <0|1>\n");
    // fprintf(stderr, "\t -save_output <0|1>\n");
    // fprintf(stderr, "\n");
  } else if (appname == "pflow") {
    /* pflow application driver does not take any additional arguments */
  } else if (appname == "sopflow") {
    print(SOPFLOWOptions::sopflow_model);
    print(SOPFLOWOptions::sopflow_solver);
    print(SOPFLOWOptions::opflow_model);
    print(SOPFLOWOptions::subproblem_model);
    print(SOPFLOWOptions::subproblem_solver);
    print(SOPFLOWOptions::ctgcfile);
    print(SOPFLOWOptions::iscoupling);
    print(SOPFLOWOptions::Ns);
    print(SOPFLOWOptions::Nc);
    print(SOPFLOWOptions::mode);
    print(SOPFLOWOptions::enable_multicontingency);
    print(SOPFLOWOptions::enable_multicontingency);
    print(SOPFLOWOptions::tolerance);
    print(SOPFLOWOptions::windgen);
  } else if (appname == "scopflow") {
    print(SCOPFLOWOptions::model);
    print(SCOPFLOWOptions::solver);
    print(SCOPFLOWOptions::subproblem_model);
    print(SCOPFLOWOptions::subproblem_solver);
    print(SCOPFLOWOptions::Nc);
    print(SCOPFLOWOptions::mode);
    print(SCOPFLOWOptions::enable_multiperiod);

    print(SCOPFLOWOptions::tolerance);
    print(SCOPFLOWOptions::dT);
    print(SCOPFLOWOptions::duration);
    print(SCOPFLOWOptions::windgenprofile);
    print(SCOPFLOWOptions::ploadprofile);
    print(SCOPFLOWOptions::qloadprofile);

  } else if (appname == "tcopflow") {
    // TODO: Update to use `ExaGOOption`s for options handling command line
    // options
    fmt::print("WARNING: TCOPFLOW application is not stable and potentially "
               "very out of date.");
    fprintf(stderr, "============== Help Options for application tcopflow "
                    "==============\n");
    fprintf(stderr, " General usage: ./%s <options>\n", appname.c_str());
    fprintf(stderr, " Options:\n");
    fprintf(stderr, "\t -netfile <netfilename>\n");
    fprintf(stderr, "\t -tcopflow_windgenprofile <windgenproffilename>\n");
    fprintf(stderr, "\t -tcopflow_ploadprofile <ploadprofile_filename>\n");
    fprintf(stderr, "\t -tcopflow_qloadprofile <qloadprofile_filename>\n");
    fprintf(stderr, "\t -tcopflow_iscoupling <0|1>\n");
    fprintf(stderr, "\t -opflow_ignore_lineflow_constraints <0|1>\n");
    fprintf(stderr, "\t -tcopflow_duration <Hours>\n");
    fprintf(stderr, "\t -tcopflow_dT <Minutes>\n");
    fprintf(stderr, "\t -tcopflow_tolerance <1e-6|...>\n");
    fprintf(stderr, "\n");
  } else {
    fmt::print("Please enter a valid application name.\n");
  }
#endif

  std::exit(EXIT_SUCCESS);
}

/**
 * Print help, version info if options are set
 */
PetscErrorCode ExaGOPrintHelpVersionInfo(int *argc, char ***argv,
                                         char *appname) {
  PetscErrorCode ierr;
  // If ExaGO is called in a unique way
  // Currently called by:
  //  - Logger Unit test
  //  - Python wrapper
  if (argc == NULL || argv == NULL) {
    return 0;
  }
  // *argv, *argv + *argc segfaults if NULL
  std::vector<std::string> args(*argv, *argv + *argc);
  auto args_it = args.begin();
  args_it++; /* Skip first argument, eg binary name */

  for (; args_it != args.end(); args_it++) {
    /* Help message */
    if (*args_it == "--help" or *args_it == "-help" or *args_it == "-h") {
      (*PetscHelpPrintf)(MPI_COMM_NULL, (appname ? appname : ""));
    }

    if (*args_it == "-c" or *args_it == "--config" or *args_it == "-config") {
      std::string vi;
      ierr = ExaGOVersionGetFullVersionInfo(vi);
      ExaGOCheckError(ierr);
      fmt::print("{}\n", vi);
      std::exit(EXIT_SUCCESS);
    }

    /* Version message */
    else if (*args_it == "-v" or *args_it == "-version" or
             *args_it == "--version") {
      std::string versionstr;
      ierr = ExaGOVersionGetVersionStr(versionstr);
      ExaGOCheckError(ierr);
      fmt::print("ExaGO {} built on {}\n", versionstr, __DATE__);
      std::exit(EXIT_SUCCESS);
    }
  }
  return 0;
}

static int initialized;

PetscErrorCode ExaGOInitialize(MPI_Comm comm, int *argc, char ***argv,
                               char *appname, char *help) {
  int i;
  PetscErrorCode ierr;
  strcpy(ExaGOCurrentAppName, appname);
  PetscHelpPrintf = &ExaGOHelpPrintf;

  /* This will print help or version info when -h or -version options are set */
  ierr = ExaGOPrintHelpVersionInfo(argc, argv, appname);

  /* Initialize MPI communicator, if not already */

  MPI_Initialized(&initialized);
  if (!initialized) {
    comm = MPI_COMM_WORLD;
    MPI_Init(argc, argv);
  }

  /* Set up ExaGO logger */
  ExaGOLogSetComm(comm);
  ierr = ExaGOLogInitialize();
  if (ierr) {
    fmt::print("Could not initialize ExaGO logger.\n");
    return ierr;
  }

  /* Call PetscInitialize without any options file */
  PETSC_COMM_WORLD = comm;
  ierr = PetscInitialize(argc, argv, nullptr, help);
  CHKERRQ(ierr);

  /* Insert options file */
  char optfile_name[PETSC_MAX_PATH_LEN];
  PetscBool flg;
  ierr = PetscOptionsGetString(NULL, NULL, "-options_file", optfile_name,
                               sizeof(optfile_name), &flg);
  CHKERRQ(ierr);
  if (flg) {
    ierr = PetscOptionsInsertFile(comm, NULL, optfile_name, PETSC_TRUE);
    CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode ExaGOFinalize() {
  PetscFunctionBegin;
  ExaGOLog(EXAGO_LOG_INFO, "Finalizing {} application.", ExaGOCurrentAppName);
  bool flg;
  ExaGOLogIsUsingLogFile(&flg);
  if (flg) {
    std::string filename;
    ExaGOLogGetLoggerName(filename);
    ExaGOLog(EXAGO_LOG_INFO, "See logfile {} for output.", filename);
  }
  PetscFinalize();

  if (!initialized) {
    MPI_Finalize();
  }
  PetscFunctionReturn(0);
}

#undef EXAGO_LOG_ENSURE_INITIALIZED

/* Used by python interface to get communicator */
PetscErrorCode ExaGOGetSelfCommunicator(MPI_Comm *comm) {
  PetscFunctionBegin;
  *comm = PETSC_COMM_SELF;
  PetscFunctionReturn(0);
}
