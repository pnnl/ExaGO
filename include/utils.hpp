#ifndef EXAGO_UTILS_H
#define EXAGO_UTILS_H
#include <petsc.h>
#include <string>
#include <vector>
#include "exago_config.h"
#include <string>
#include <stdexcept>

enum ExaGOVerbosityLevel {
  EXAGO_LOG_INFO=0,
  EXAGO_LOG_WARN,
  EXAGO_LOG_ERROR,
  EXAGO_LOG_DISABLE,
  EXAGO_LOG_NUM_LOG_LEVELS,
};

extern const char* ExaGOVerbosityNames[EXAGO_LOG_NUM_LOG_LEVELS];

/**
 * @brief ExaGO Error interfaces between error codes encountered in ExaGO
 * proper as well as error codes recieved from PETSc.
 */
struct ExaGOError : public std::exception {
  ExaGOError() : message{"ExaGO Error"} {}
  ExaGOError(const char* message) : message{message} {}
  ExaGOError(PetscErrorCode);

  /* The name _what_ is not in PascalCase like the rest of ExaGO because
   * ExaGOError inherits from the standard library exception which defines
   * _what_. */
  virtual const char* what() const noexcept { return message.c_str(); };
  virtual bool IsPetscError() const noexcept { return is_petsc_error; }
protected:
  std::string message;
  bool is_petsc_error = false;
};

/* Used to interface Petsc error return codes with ExaGO errors */
extern void ExaGOCheckError(int e);

/**
 * Set the name for the logfile to be used; must be set before
 * `ExaGOInitialize` is called.
 */
extern PetscErrorCode ExaGOLogSetLoggingFileName(char*);

/** Retrieve the filename being used for the logfile */
extern PetscErrorCode ExaGOLogGetLoggingFileName(char**);

/** Indicates whether a logfile is being used or stdout/stderr */
extern PetscErrorCode ExaGOLogIsUsingLogFile(PetscBool*);

/** Get minimum loglevel for a log to be printed */
extern PetscErrorCode ExaGOLogGetMinLogLevel(ExaGOVerbosityLevel*);

/** Set minimum loglevel for a log to be printed */
extern PetscErrorCode ExaGOLogSetMinLogLevel(ExaGOVerbosityLevel);

/** Parameter indicates whether each MPI rank should print out each log
 * message. ExaGOLog is not currently able to log to every rank
 * /and/ log to a file. */
extern PetscErrorCode ExaGOLogUseEveryRank(PetscBool);

#if !defined(EXAGO_DISABLE_LOGGING)
#include <spdlog/spdlog.h>

/** 
 * @brief Implementation to log string according to ExaGO build configuration.
 * @note Users should use `ExaGOLog` instead, as this can be disabled for
 * optimization.
 */
/*
extern void ExaGOLogImpl(ExaGOVerbosityLevel,const char*,...);
*/
template<typename... Args> void ExaGOLogImpl(ExaGOVerbosityLevel level, std::string fmt, Args... args)
{
  auto logger = spdlog::get("exago_logger");
  switch (level) {
    case EXAGO_LOG_INFO:
      logger->info(fmt,args...);
      break;
    case EXAGO_LOG_WARN:
      logger->warn(fmt,args...);
      break;
    case EXAGO_LOG_ERROR:
      logger->error(fmt,args...);
      break;
    case EXAGO_LOG_DISABLE:
      logger->info(fmt,args...);
      break;
    case EXAGO_LOG_NUM_LOG_LEVELS:
      logger->info(fmt,args...);
      break;
    default:
      break;
  }
}

/**
 * @brief `ExaGOLog` is the user-facing logging function.
 * @warning Cannot take only two arguments - must provide a format and
 * parameters. For example, if you only need to log a single string, prefer:
 * `ExaGOLog(EXAGO_LOG_INFO,"%s","My interesting log message");`.
 */
#define ExaGOLog(level,fmt,...) ExaGOLogImpl(level,fmt,__VA_ARGS__)

#else

/** If logging is disabled, logging is a no-op */
#define ExaGOLog(x,y,...) (void)(x)

#endif

extern PetscErrorCode ExaGOLogIsInitialized(PetscBool*);

/**
 * Initialize an ExaGO application.
 *
 * @note this takes care of Petsc initialization, so don't this function in
 * conjunction with `PetscInitialize.`
 */
extern "C" PetscErrorCode ExaGOInitialize(MPI_Comm,int*argc,char***argv,
    char*appname,char*help);

/**
 * Teardown for an ExaGO application.
 *
 * @note this takes care of Petsc finalization, so don't this function in
 * conjunction with `PetscFinalize`.
 */
extern "C" PetscErrorCode ExaGOFinalize();

/** Returns 1 if files exists, else 0 */
bool DoesFileExist(const char*);

/** Returns 1 if directory exists, else 0 */
bool DoesDirExist(const char*);

/** First path in _files_ to be statable as a regular file */
std::vector<std::string>::const_iterator FirstExistingFile(
    const std::vector<std::string> &files);

inline bool IsEqual(double value, double ref, double tolerance, double& error)
{
  error = std::fabs(value-ref)/(1. + std::fabs(ref));
  return (error < tolerance);
}

inline bool IsEqual(int value, int ref, int tolerance, int& error)
{
  error = std::abs(value - ref);
  return (error < tolerance);
}

/** Used by python interface to get self communicattor */
PETSC_EXTERN PetscErrorCode ExaGOGetSelfCommunicator(MPI_Comm*);

#endif
