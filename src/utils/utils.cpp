#include <common.h>
#include <dirent.h>
#include <exago_config.h>
#include <exago_build_options.h>
#include <iostream>
#include <string>
#include <regex>
#include <stdarg.h>
#include <sys/stat.h>
#include <sys/utsname.h>
#include <version.h>
#include <utils.h>

PetscBool printHelp;

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
 * the network file and print/save options into the application structures
 * and handle them in the solve/finalize functions. For now, they will be
 * omitted from the help message. */
// const auto network = ExaGOStringOption("-netfile", "Path to network file",
//                                        "/path/to/netfile", {});
} // namespace ExaGODefaultOptions

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
#ifdef EXAGO_ENABLE_LOGGING
  *flg = ExaGOLogUseFile;
#endif
  return 0;
}

PetscErrorCode ExaGOLogGetLoggerName(std::string &s) {
  s = ExaGOLoggerName;
  return 0;
}

#ifdef EXAGO_ENABLE_LOGGING
PetscErrorCode ExaGOGetLoggerPointer(std::shared_ptr<spdlog::logger> &logger) {
  logger = spdlog::get(ExaGOLoggerName);
  return 0;
}
#endif

/** Set ExaGO logging communicator */
PetscErrorCode ExaGOLogSetComm(MPI_Comm c) {
  ExaGOLogComm = c;
  return 0;
}

PetscErrorCode ExaGOLogInitialize() {
#ifdef EXAGO_ENABLE_LOGGING
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
#endif
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


/**
 * @brief Print ExaGO command-line options.
 *
 * @pre Any explicitly managed resources are cleaned up before this is called -
 * this function will std::exit for you.
 *
 */
PetscErrorCode ExaGOHelpPrintf() {
  PetscErrorCode ierr;

  std::string versionstr;
  ExaGOVersionGetVersionStr(versionstr);
#ifdef EXAGO_ENABLE_LOGGING
  fmt::print("ExaGO {} built on {}\n", versionstr, __DATE__);
#else
  printf("ExaGO %s built on %s\n", versionstr.c_str(), __DATE__);
#endif

#ifdef EXAGO_ENABLE_LOGGING
  fmt::print("\nExaGO options\n");
#else
  printf("\nExaGO options\n");
#endif
  print(ExaGODefaultOptions::help);
  print(ExaGODefaultOptions::version);
  print(ExaGODefaultOptions::config);
  print(ExaGODefaultOptions::options_file);

  return 0;

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
      printHelp = PETSC_TRUE;
      ExaGOHelpPrintf();
    }

    if (*args_it == "-c" or *args_it == "--config" or *args_it == "-config") {
      std::string vi;
      ierr = ExaGOVersionGetFullVersionInfo(vi);
      ExaGOCheckError(ierr);
#ifdef EXAGO_ENABLE_LOGGING
      fmt::print("{}\n", vi);
#else
      printf("%s\n", vi.c_str());
#endif
      std::exit(EXIT_SUCCESS);
    }

    /* Version message */
    else if (*args_it == "-v" or *args_it == "-version" or
             *args_it == "--version") {
      std::string versionstr;
      ierr = ExaGOVersionGetVersionStr(versionstr);
      ExaGOCheckError(ierr);
#ifdef EXAGO_ENABLE_LOGGING
      fmt::print("ExaGO {} built on {}\n", versionstr, __DATE__);
#else
      printf("ExaGO %s built on %s\n", versionstr.c_str(), __DATE__);
#endif
      std::exit(EXIT_SUCCESS);
    }
  }
  return 0;
}

// This function is used to suppress printing of PETSc help
PetscErrorCode PetscHelpPrintfNone(MPI_Comm comm,const char empty[],...)
{
  return 0;
}

static int initialized;

PetscErrorCode ExaGOInitialize(MPI_Comm comm, int *argc, char ***argv,
                               char *appname, char *help) {
  int i;
  PetscErrorCode ierr;
  strcpy(ExaGOCurrentAppName, appname);

  // Suppress printing PETSc help
  PetscHelpPrintf = &PetscHelpPrintfNone;

  printHelp = PETSC_FALSE;
  
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
#ifdef EXAGO_ENABLE_LOGGING
    fmt::print("Could not initialize ExaGO logger.\n");
#else
    printf("Could not initialize ExaGO logger.\n");
#endif
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
