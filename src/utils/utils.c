#include <common.h>
#include <private/utilsimpl.h>
#include <sys/utsname.h>
#include <sys/stat.h>
#include <stdarg.h>
#include <dirent.h>
#include <exago_config.h>

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

/** Rank from which logging statements will be written */
static int ExaGOLogLoggingRank=0;

/** Current MPI rank */
static int ExaGOLogCurrentRank=0;

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

  /** If we're not on the logging rank, there's no setup to perform */
  ierr=MPI_Comm_rank(ExaGOLogComm,&ExaGOLogCurrentRank);CHKERRQ(ierr);
  if (ExaGOLogCurrentRank!=ExaGOLogLoggingRank)
  {
    ExaGOLogIsLoggerInitialized = PETSC_TRUE;
    return 0;
  }

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

  /** If we're not on the logging rank, there's no teardown to perform */
  if (ExaGOLogCurrentRank!=ExaGOLogLoggingRank)
  {
    ExaGOLogIsLoggerInitialized = PETSC_FALSE;
    return 0;
  }
 

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

PetscErrorCode ExaGOLogUseEveryRank(PetscBool use)
{
  ExaGOLogLoggingRank = use ? -1 : 0;
  return 0;
}

/** ExaGO logging function with format string */
void ExaGOLogImpl(ExaGOVerbosityLevel verb, char *fmt, ...)
{
  EXAGO_LOG_ENSURE_INITIALIZED();
  int ierr;
  if (ExaGOLogLoggingRank>=0 && ExaGOLogCurrentRank!=ExaGOLogLoggingRank)
    return;

  if(verb<ExaGOLogMinLogLevel)
    return;

  FILE *fp;
  ierr=ExaGOLogGetLoggingFilePointer(&fp);

  char buf[1024];
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
}

/**
 * @brief Checks whether or not a directory entry can be opened at the given
 * path using only POSIX standard C system calls
 *
 * @param[in] pth filepath to stat
 * @see stat
 **/
int doesFileExist(char* pth)
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
int doesDirExist(char* pth)
{
  if (pth==NULL) return 0;
  DIR *dp;
  dp = opendir(pth);
  if (dp == NULL) return 0;
  return 1;
}

/**
 * @brief Verifies that all paths passed are statable as regular files
 *
 * @param[in] pths  array of paths to check exist
 * @param[in] npths number of paths to check
 * @return first i in [0, npths) such that pths[i] is statable as a regular file
 *         -2 if pths==NULL or pths[i]==NULL for any i in [0, npths),
 *         -1 if !doesFileExist(pths[i]) for i in [0, npths)
 * 
 * @see doesFileExist
 **/
int anyFileExist(char** pths, int npths)
{
  char buf[1024];
  if (pths==NULL)
  {
    ExaGOLog(EXAGO_LOG_INFO,"Function `%s` was called with argument"
        " pths==NULL.\n",__func__);
    return -2;
  }

  for (int i=0; i<npths; i++)
  {
    if (pths[i]==NULL)
    {
      continue;
    }
    sprintf(buf,"-- Checking %-70s exists: ",pths[i]);
    if (doesFileExist(pths[i]))
    {
      strcat(buf,"yes");
      ExaGOLog(EXAGO_LOG_INFO,"%s",buf);
      return i;
    }
    else
    {
      strcat(buf,"no");
      ExaGOLog(EXAGO_LOG_INFO,"%s",buf);
    }
  }
  return -1;
}

int isEqual(double value,double reference,double tol,double*err)
{
  *err = fabs(value-reference)/(1. + fabs(reference));
  return (*err < tol);
}

int isEqualInt(int a,int b,int tol,int*err)
{
  *err = abs(a-b);
  return (*err < tol);
}

#undef EXAGO_LOG_ENSURE_INITIALIZED
