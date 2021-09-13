#include <iostream>
#include <common.h>
#include <sys/utsname.h>
#include <sys/stat.h>
#include <version.hpp>
#include <stdarg.h>
#include <dirent.h>
#include <exago_config.h>
#include <utils.hpp>

const char* ExaGOVerbosityNames[] = {"INFO","WARNING","ERROR"};

static char ExaGOCurrentAppName[128];

ExaGOError::ExaGOError(PetscErrorCode ierr) : is_petsc_error{true} {
  const char *error_message;
  char *specific_error_message;
  PetscErrorMessage(ierr, &error_message, &specific_error_message);
  message = error_message;
  message += specific_error_message;
}

void ExaGOCheckError(int e)
{
  if (static_cast<int>(e) > PETSC_ERR_MIN_VALUE
      && static_cast<int>(e) < PETSC_ERR_MAX_VALUE)
    throw ExaGOError(e);
}

static PetscBool ExaGOLogIsLoggerInitialized=PETSC_FALSE;

/** This macro is undef'ed at the end of this header, don't use elsewhere.
    Only for use within PetscErrorCode-returning ExaGOLog-related functions. */
#define EXAGO_LOG_ENSURE_INITIALIZED() \
  do { \
    if (!ExaGOLogIsLoggerInitialized) \
    { \
      fprintf(stderr,"An `ExaGOLog` function has been called, but ExaGOLog" \
          " is not initialized!\n"); \
      exit(1); \
    } \
  } while(0)

static FILE* ExaGOLogFilePointer=NULL;

static char ExaGOLogFileName[PETSC_MAX_PATH_LEN]="filename_unset";

static PetscBool ExaGOLogUseLogFile=PETSC_FALSE;

static ExaGOVerbosityLevel ExaGOLogMinLogLevel = EXAGO_LOG_INFO;

/** ExaGO logging communicator to determine rank*/
static MPI_Comm ExaGOLogComm=MPI_COMM_SELF;

PetscErrorCode ExaGOLogGetLoggingFileName(char** name)
{
  if(ExaGOLogFileName==NULL)
    return 1;
  *name = strdup(ExaGOLogFileName);
  return 0;
}

PetscErrorCode ExaGOLogSetLoggingFileName(char* name)
{
  strcpy(ExaGOLogFileName,name);
  ExaGOLogUseLogFile = PETSC_TRUE;
  return 0;
}

PetscErrorCode ExaGOLogIsUsingLogFile(PetscBool *flg)
{
  *flg = ExaGOLogUseLogFile;
  return 0;
}

PetscErrorCode ExaGOLogGetLoggingFilePointer(FILE** fp)
{
  if(ExaGOLogFilePointer==NULL)
    return 1;
  *fp = ExaGOLogFilePointer;
  return 0;
}

PetscErrorCode ExaGOLogSetLoggingFilePointer(FILE* fp)
{
  ExaGOLogFilePointer = fp;
  return 0;
}

/** Set ExaGO logging communicator */
PetscErrorCode ExaGOLogSetComm(MPI_Comm c)
{
  ExaGOLogComm = c;
  return 0;
}

PetscErrorCode ExaGOLogIsInitialized(PetscBool*flg)
{
  *flg = ExaGOLogIsLoggerInitialized;
  return 0;
}

PetscErrorCode ExaGOLogInitialize()
{
  PetscBool flg=PETSC_FALSE;
  PetscErrorCode ierr=0;
  char* filename;
  FILE* fp;

  ierr=ExaGOLogIsUsingLogFile(&flg);
  /** Normally we would check the error, but in here we don't even have logging
   * set up to handle the error so we have to just return a code. */
  if(ierr != 0)
    return 1;

  /** If we're not using a logfile, there's no additional setup to perform */
  if (!flg)
    goto lbl_success; // so sue me

  ierr=ExaGOLogGetLoggingFileName(&filename);
  /** We still don't have proper logging set up, so just return the error */
  if(ierr != 0)
    return 1;

  fp = fopen(ExaGOLogFileName, "w");
  ExaGOLogSetLoggingFilePointer(fp);

lbl_success:
  ExaGOLogSetLoggingFilePointer(stdout);
  ExaGOLogIsLoggerInitialized = PETSC_TRUE;
  return 0;
}

PetscErrorCode ExaGOLogFinalize()
{
  EXAGO_LOG_ENSURE_INITIALIZED();

  PetscBool flg=PETSC_FALSE;
  PetscErrorCode ierr=0;
  ierr=ExaGOLogIsUsingLogFile(&flg);
  /** If we're not using a logfile, there's no teardown to perform */
  if (!flg)
    return 0;

  FILE *fp;
  ExaGOLogGetLoggingFilePointer(&fp);
  fclose(fp);

  ExaGOLogIsLoggerInitialized = PETSC_FALSE;
  return 0;
}

PetscErrorCode ExaGOLogGetMinLogLevel(ExaGOVerbosityLevel *l)
{
  *l = ExaGOLogMinLogLevel;
  return 0;
}

PetscErrorCode ExaGOLogSetMinLogLevel(ExaGOVerbosityLevel l)
{
  ExaGOLogMinLogLevel = l;
  return 0;
}

/** ExaGO logging function with format string */
void ExaGOLogImpl(ExaGOVerbosityLevel verb, const char *fmt, ...)
{
  EXAGO_LOG_ENSURE_INITIALIZED();
  int ierr;

  if(verb<ExaGOLogMinLogLevel)
    return;

  FILE *fp;
  ierr=ExaGOLogGetLoggingFilePointer(&fp);

  static int log_bufsize = 32768;
  char buf[log_bufsize];
  va_list args;
  va_start(args, fmt);
  vsprintf(buf, fmt, args);
  char *endl = strtok(buf, "\n");
  while (endl!=NULL)
  {
    fprintf(fp,"[ExaGO %s]: %s\n",ExaGOVerbosityNames[verb],endl);
    endl = strtok(NULL, "\n");
  }
  va_end(args);
  fflush(fp);
}

/**
 * @brief Checks whether or not a directory entry can be opened at the given
 * path using only POSIX standard C system calls
 *
 * @param[in] pth filepath to stat
 * @see stat
 **/
bool DoesFileExist(const char* pth)
{
  struct stat path_stat;
  stat(pth, &path_stat);
  return S_ISREG(path_stat.st_mode);
}

/**
 * @brief Checks whether or not a directory entry can be opened at the given
 * path using only POSIX standard C system calls.
 *
 * @param[in] path path to verify is statable
 * @return 0 if pth==NULL or pth cannot be opened as directory with syscall
 *         1 else
 *
 * @see opendir
 **/
bool DoesDirExist(const char* pth)
{
  if (pth==NULL) return false;
  DIR *dp;
  dp = opendir(pth);
  if (dp == NULL) return false;
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
std::vector<std::string>::const_iterator FirstExistingFile(
    const std::vector<std::string> &files)
{
  auto it = files.begin();
  for (; it != files.end(); it++)
  {
    auto message = std::string{"-- Checking "} + *it + " exists: ";
    if (DoesFileExist(it->c_str()))
    {
      message += "yes";
      ExaGOLog(EXAGO_LOG_INFO, "%s", message.c_str());
      return it;
    }
    message += "no";
    ExaGOLog(EXAGO_LOG_INFO, "%s", message.c_str());
  }

  /* No files were found to exist! Just return iterator to files.end() to follow
   * conventions in the STL. */
  return it;
}

#undef EXAGO_LOG_ENSURE_INITIALIZED

/**
 * @brief Print ExaGO command-line options.
 *
 * @pre Any explicitly managed resources are cleaned up before this is called -
 * this function will std::exit for you.
 */
PetscErrorCode ExagoHelpPrintf(MPI_Comm comm, const char appname[],...)
{
  PetscErrorCode ierr;

  char* versionstr;
  ierr=ExaGOVersionGetFullVersionInfo(&versionstr);CHKERRQ(ierr);
  fprintf(stderr, "===================================================================\n");
  fprintf(stderr, "ExaGO Version Info:\n\n%s", versionstr); 
  if (strcmp(appname,"opflow") == 0)
  {
    fprintf(stderr, "============== Help Options for application %s ==============\n\n", appname); 
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
  } 
  else if (strcmp(appname,"pflow")==0) 
  {
    fprintf(stderr, "============== Help Options for application %s ==============\n", appname);
  } 
  else if (strcmp(appname,"sopflow") == 0) 
  {
    fprintf(stderr, "============== Help Options for application %s ==============\n", appname);
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
  }
  else if (strcmp(appname,"scopflow") == 0) 
  {
    fprintf(stderr, "============== Help Options for application %s ==============\n", appname);
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

  }
  else if (strcmp(appname,"tcopflow") == 0) 
  {
    fprintf(stderr, "============== Help Options for application %s ==============\n", appname);
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
  }
  else 
  {
    fprintf(stderr, "Please enter a valid application name.\n");
  }
  std::exit(EXIT_SUCCESS);
}

/**
 * Parse arguments for initializing an ExaGO application.
 */
PetscErrorCode ExaGOInitializeParseArgs(int* argc, char*** argv, char* appname,
    char *help, std::vector<std::string> optfiles, bool* no_optfile)
{
  PetscErrorCode ierr;
  std::vector<std::string> args(*argv, *argv + *argc);
  auto args_it = args.begin();
  args_it++; /* Skip first argument, eg binary name */

  for(; args_it != args.end(); args_it++)
  {
    /* Help message */
    if (*args_it == "-help" or *args_it == "-h")
    {
      (*PetscHelpPrintf)(MPI_COMM_NULL, appname);
    }

    /* Version message */
    else if (*args_it == "-v" or *args_it == "-version")
    {
      char* versionstr;
      ierr = ExaGOVersionGetFullVersionInfo(&versionstr);
      ExaGOCheckError(ierr);
      std::cerr << "ExaGO Version Info:\n\n" << versionstr
        << "ExaGO Help Message:\n" << help;
      free(versionstr);
      MPI_Finalize();
      std::exit(EXIT_SUCCESS);
    }

    /* Paths to options files */
    else if (*args_it == "-options_file")
    {
      /* Ensure the '-options_file' flag was not the last flag! If it was,
       * give the help message. */
      if (args_it + 1 == args.end())
        (*PetscHelpPrintf)(MPI_COMM_NULL, appname);

      args_it++;
      optfiles.push_back(*args_it);
    }
    /* If we need to skip options file */
    else if (*args_it == "-no_optfile")
    {
      *no_optfile = true;
    }
  }
  return 0;
}

PetscErrorCode ExaGOInitialize(MPI_Comm comm, int* argc, char*** argv,
    char* appname, char *help)
{
  int i;
  PetscErrorCode ierr;
  strcpy(ExaGOCurrentAppName,appname);
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

  if (argc == NULL || argv == NULL || *argc == 0)
  {
    skip_args = true;
    argc = NULL;
    argv = NULL;
    no_optfile = true;
    comm = MPI_COMM_WORLD;
  }

  if (!skip_args)
    ExaGOInitializeParseArgs(argc, argv, appname, help, optfiles, &no_optfile);

  int initialized;
  MPI_Initialized(&initialized);
  if(!initialized)
    MPI_Init(argc,argv);

  ExaGOLogSetComm(comm);
  ierr = ExaGOLogInitialize();
  if(ierr)
  {
    fprintf(stderr, "Could not initialize ExaGO logger.\n");
    return ierr;
  }

  FILE *fp;
  ExaGOLogGetLoggingFilePointer(&fp);
  if(fp==NULL)
    ExaGOLogSetLoggingFilePointer(stderr);

  // Skip options file setting if -no_optfile is passed
  if (no_optfile)
  {
    PetscInitialize(argc,argv,NULL,help);
    return 0;
  }

  optfiles.push_back(std::string{options_pathname} + "/" + filename);
  optfiles.push_back("./" + std::string{filename});
  optfiles.push_back("./options/" + std::string{filename});

  auto file_it = FirstExistingFile(optfiles);

  if (file_it == optfiles.end())
  {
    ExaGOLog(EXAGO_LOG_ERROR,
        "Could not find options file for application %s.\n",
        ExaGOCurrentAppName);
    throw std::runtime_error{"No options files were found"};
  }
 
  PetscInitialize(argc,argv,file_it->c_str(),help);
  return 0;
}


PetscErrorCode ExaGOFinalize()
{
  PetscFunctionBegin;
  ExaGOLog(EXAGO_LOG_INFO,"Finalizing %s application.\n",ExaGOCurrentAppName);
  PetscBool flg;
  ExaGOLogIsUsingLogFile(&flg);
  if (flg)
  {
    char* filename;
    ExaGOLogGetLoggingFileName(&filename);
    ExaGOLog(EXAGO_LOG_INFO,"See logfile %s for output.\n",filename);
  }
  ExaGOLogFinalize();
  PetscFinalize();
  MPI_Finalize();
  PetscFunctionReturn(0);
}

#undef EXAGO_LOG_ENSURE_INITIALIZED

/* Used by python interface to get communicator */
PetscErrorCode ExaGOGetSelfCommunicator(MPI_Comm *comm)
{
  PetscFunctionBegin;
  *comm = PETSC_COMM_SELF;
  PetscFunctionReturn(0);
}
