#include <version.h>
#include <common.h>
#include <utils.h>
#include <exago_config.h>

const char* ExaGODependencyNames[ExaGONumDependencies]
  = {"PETSC","MPI","Ipopt","HiOp","GPU","RAJA"};

const PetscBool ExaGOIsDependencyEnabled[ExaGONumDependencies] =
{
#ifdef EXAGO_ENABLE_PETSC
    PETSC_TRUE,
#else
    PETSC_FALSE,
#endif

#ifdef EXAGO_ENABLE_MPI
    PETSC_TRUE,
#else
    PETSC_FALSE,
#endif

#ifdef EXAGO_ENABLE_IPOPT
    PETSC_TRUE,
#else
    PETSC_FALSE,
#endif

#ifdef EXAGO_ENABLE_HIOP
    PETSC_TRUE,
#else
    PETSC_FALSE,
#endif

#ifdef EXAGO_ENABLE_GPU
    PETSC_TRUE,
#else
    PETSC_FALSE,
#endif

#ifdef EXAGO_ENABLE_RAJA
    PETSC_TRUE,
#else
    PETSC_FALSE,
#endif
};

PetscErrorCode ExaGOVersionGetFullVersionInfo(char** str)
{
  *str = malloc(2048);
  strcat(*str, "ExaGO version ");
  strcat(*str, EXAGO_VERSION);
  strcat(*str, " released on ");
  strcat(*str, EXAGO_RELEASE_DATE);
  strcat(*str, "\nbuilt with:\n");
  int i;
  for(i=0; i<ExaGONumDependencies; i++)
  {
    char buf[1024];
    sprintf(buf, "\t%-20s%20s\n", ExaGODependencyNames[i], ExaGOIsDependencyEnabled[i]?"YES":"NO");
    strcat(*str, buf);
  }
  return 0;
}

PetscErrorCode ExaGOVersionGetReleaseDate(char** str)
{
  *str = strdup(EXAGO_RELEASE_DATE);
  return 0;
}

PetscErrorCode ExaGOVersionGetVersion(int *major,int *minor,int *patch)
{
  *major = atoi(EXAGO_VERSION_MAJOR);
  *minor = atoi(EXAGO_VERSION_MINOR);
  *patch = atoi(EXAGO_VERSION_PATCH);
  return 0;
}

PetscErrorCode ExaGOVersionGetVersionStr(char **str)
{
  *str = strdup(EXAGO_VERSION);
  return 0;
}
