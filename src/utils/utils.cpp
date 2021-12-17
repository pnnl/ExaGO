#include <common.h>
#include <dirent.h>
#include <exago_config.h>
#include <iostream>
#include <string>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <sys/utsname.h>
#include <utils.h>
#include <version.h>

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
 * @brief Print ExaGO command-line options.
 *
 * @pre Any explicitly managed resources are cleaned up before this is called -
 * this function will std::exit for you.
 */
PetscErrorCode ExagoHelpPrintf(MPI_Comm comm, const char appname[], ...) {
  PetscErrorCode ierr;

  char *versionstr;
  ierr = ExaGOVersionGetFullVersionInfo(&versionstr);
  CHKERRQ(ierr);
  fprintf(
      stderr,
      "===================================================================\n");
  fprintf(stderr, "ExaGO Version Info:\n\n%s", versionstr);
  if (strcmp(appname, "opflow") == 0) {
    fprintf(stderr,
            "============== Help Options for application %s ==============\n\n",
            appname);
    fprintf(stderr, " General usage: mpiexec -n <N> ./%s <options>\n", appname);
    fprintf(stderr, " Options:\n");
    fprintf(stderr, "\t -netfile <netfilename>\n");
    fprintf(stderr, "\t -opflow_model <POWER_BALANCE_POLAR|...>\n");
    fprintf(stderr, "\t -opflow_solver <IPOPT|...>\n");
    fprintf(stderr, "\t -opflow_initialization <MIDPOINT|...>\n");
    fprintf(stderr, "\t -opflow_ignore_lineflow_constraints <0|1>\n");
    fprintf(stderr, "\t -opflow_include_loadloss_variables <0|1>\n");
    fprintf(stderr, "\t -opflow_include_powerimbalance_variables <0|1>\n");
    fprintf(stderr, "\t -opflow_loadloss_penalty <Penalty ($)>\n");
    fprintf(stderr, "\t -opflow_powerimbalance_penalty <Penalty ($)>\n");
    fprintf(stderr, "\t -opflow_genbusvoltage <FIXED_WITHIN_QBOUNDS|...>\n");
    fprintf(stderr, "\t -opflow_has_gensetpoint <0|1>\n");
    fprintf(stderr, "\t -opflow_objective <MIN_GEN_COST|...>\n");
    fprintf(stderr, "\t -opflow_use_agc <0|1>\n");
    fprintf(stderr, "\t -opflow_tolerance <1e-6|...>\n");
    fprintf(stderr, "\t -hiop_compute_mode <hybrid|...>\n");
    fprintf(stderr, "\t -hiop_verbosity_level <0-10>\n");
    fprintf(stderr, "\t -hiop_tolerance <1e-6|...>\n");
    fprintf(stderr, "\t -print_output <0|1>\n");
    fprintf(stderr, "\t -save_output <0|1>\n");
    fprintf(stderr, "\n");
  } else if (strcmp(appname, "pflow") == 0) {
    fprintf(stderr,
            "============== Help Options for application %s ==============\n",
            appname);
  } else if (strcmp(appname, "sopflow") == 0) {
    fprintf(stderr,
            "============== Help Options for application %s ==============\n",
            appname);
    fprintf(stderr, " General usage: mpiexec -n <N> ./%s <options>\n", appname);
    fprintf(stderr, " Options:\n");
    fprintf(stderr, "\t -netfile <netfilename>\n");
    fprintf(stderr, "\t -ctgcfile <ctgcfilename>\n");
    fprintf(stderr, "\t -scenfile <scenfilename>\n");
    fprintf(stderr, "\t -opflow_model <POWER_BALANCE_POLAR|...>\n");
    fprintf(stderr, "\t -sopflow_solver <IPOPT|...>\n");
    fprintf(stderr, "\t -sopflow_mode <0|1>\n");
    fprintf(stderr, "\t -sopflow_Ns <Ns>\n");
    fprintf(stderr, "\t -sopflow_enable_multicontingency <0|1>\n");
    fprintf(stderr, "\t -scopflow_enable_multiperiod <0|1>\n");
    fprintf(stderr, "\t -sopflow_tolerance <1e-6|...>\n");
    fprintf(stderr, "\n");
  } else if (strcmp(appname, "scopflow") == 0) {
    fprintf(stderr,
            "============== Help Options for application %s ==============\n",
            appname);
    fprintf(stderr, " General usage: mpiexec -n <N> ./%s <options>\n", appname);
    fprintf(stderr, " Options:\n");
    fprintf(stderr, "\t -netfile <netfilename>\n");
    fprintf(stderr, "\t -ctgcfile <ctgcfilename>\n");
    fprintf(stderr, "\t -scopflow_windgenprofile <windgenproffilename>\n");
    fprintf(stderr, "\t -scopflow_ploadprofile <ploadprofile_filename>\n");
    fprintf(stderr, "\t -scopflow_qloadprofile <qloadprofile_filename>\n");
    fprintf(stderr, "\t -opflow_model <POWER_BALANCE_POLAR|...>\n");
    fprintf(stderr, "\t -scopflow_solver <IPOPT|...>\n");
    fprintf(stderr, "\t -scopflow_mode <0|1>\n");
    fprintf(stderr, "\t -sopflow_Nc <Nc>\n");
    fprintf(stderr, "\t -scopflow_enable_multiperiod <0|1>\n");
    fprintf(stderr, "\t -scopflow_duration <Hours>\n");
    fprintf(stderr, "\t -scopflow_dT <Minutes>\n");
    fprintf(stderr, "\t -scopflow_tolerance <1e-6|...>\n");
    fprintf(stderr, "\n");

  } else if (strcmp(appname, "tcopflow") == 0) {
    fprintf(stderr,
            "============== Help Options for application %s ==============\n",
            appname);
    fprintf(stderr, " General usage: mpiexec -n <N> ./%s <options>\n", appname);
    fprintf(stderr, " Options:\n");
    fprintf(stderr, "\t -netfile <netfilename>\n");
    fprintf(stderr, "\t -options_file <tcopflowoptionsfilename>\n");
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
    fprintf(stderr, "Please enter a valid application name.\n");
  }
  std::exit(EXIT_SUCCESS);
}

/**
 * Parse arguments for initializing an ExaGO application.
 */
PetscErrorCode ExaGOInitializeParseArgs(int *argc, char ***argv, char *appname,
                                        char *help,
                                        std::vector<std::string> optfiles,
                                        bool *no_optfile) {
  PetscErrorCode ierr;
  std::vector<std::string> args(*argv, *argv + *argc);
  auto args_it = args.begin();
  args_it++; /* Skip first argument, eg binary name */

  for (; args_it != args.end(); args_it++) {
    /* Help message */
    if (*args_it == "-help" or *args_it == "-h") {
      (*PetscHelpPrintf)(MPI_COMM_NULL, appname);
    }

    /* Version message */
    else if (*args_it == "-v" or *args_it == "-version") {
      char *versionstr;
      ierr = ExaGOVersionGetFullVersionInfo(&versionstr);
      ExaGOCheckError(ierr);
      std::cerr << "ExaGO Version Info:\n\n"
                << versionstr << "ExaGO Help Message:\n"
                << help;
      free(versionstr);
      MPI_Finalize();
      std::exit(EXIT_SUCCESS);
    }

    /* Paths to options files */
    else if (*args_it == "-options_file") {
      /* Ensure the '-options_file' flag was not the last flag! If it was,
       * give the help message. */
      if (args_it + 1 == args.end())
        (*PetscHelpPrintf)(MPI_COMM_NULL, appname);

      args_it++;
      optfiles.push_back(*args_it);
    }
    /* If we need to skip options file */
    else if (*args_it == "-no_optfile") {
      *no_optfile = true;
    }
  }
  return 0;
}

PetscErrorCode ExaGOInitialize(MPI_Comm comm, int *argc, char ***argv,
                               char *appname, char *help) {
  int i;
  PetscErrorCode ierr;
  strcpy(ExaGOCurrentAppName, appname);
  PetscHelpPrintf = &ExagoHelpPrintf;

  /* Skip options file if necessary */
  bool no_optfile = false;

  /* Prefix checked for application options */
  std::string options_pathname = EXAGO_OPTIONS_DIR;

  /* Filename used when searching for application-specific options */
  std::string filename = std::string{ExaGOCurrentAppName} + "options";

  /* Paths to option files for given application */
  std::vector<std::string> optfiles;

  /* Skip argument parsing. Used when called via external programs such as
   * python */
  bool skip_args = false;

  if (argc == nullptr || argv == nullptr || *argc == 0) {
    skip_args = true;
    argc = nullptr;
    argv = nullptr;
    no_optfile = true;
    comm = MPI_COMM_WORLD;
  }

  if (!skip_args)
    ExaGOInitializeParseArgs(argc, argv, appname, help, optfiles, &no_optfile);

  int initialized;
  MPI_Initialized(&initialized);
  if (!initialized)
    MPI_Init(argc, argv);

  ExaGOLogSetComm(comm);
  ierr = ExaGOLogInitialize();
  if (ierr) {
    fprintf(stderr, "Could not initialize ExaGO logger.\n");
    return ierr;
  }

  // Skip options file setting if -no_optfile is passed
  if (no_optfile) {
    PetscInitialize(argc, argv, nullptr, help);
    return 0;
  }

  optfiles.push_back(std::string{options_pathname} + "/" + filename);
  optfiles.push_back("./" + std::string{filename});
  optfiles.push_back("./options/" + std::string{filename});

  auto file_it = FirstExistingFile(optfiles);

  if (file_it == optfiles.end()) {
    ExaGOLog(EXAGO_LOG_ERROR, "Could not find options file for application {}.",
             ExaGOCurrentAppName);
    throw std::runtime_error{"No options files were found"};
  }

  PetscInitialize(argc, argv, file_it->c_str(), help);
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
  MPI_Finalize();
  PetscFunctionReturn(0);
}

#undef EXAGO_LOG_ENSURE_INITIALIZED

/* Used by python interface to get communicator */
PetscErrorCode ExaGOGetSelfCommunicator(MPI_Comm *comm) {
  PetscFunctionBegin;
  *comm = PETSC_COMM_SELF;
  PetscFunctionReturn(0);
}
