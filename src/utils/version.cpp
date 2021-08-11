#include <version.hpp>
#include <common.h>
#include <utils.hpp>
#include <exago_config.h>
#include <unordered_map>

namespace {
  static const std::unordered_map<std::string, bool> ExaGODependencies = {
    { 
      "petsc",
#ifdef EXAGO_ENABLE_PETSC
      true,
#else
      false,
#endif
    },
    { 
      "mpi",
#ifdef EXAGO_ENABLE_MPI
      true,
#else
      false,
#endif
    },
    { 
      "ipopt",
#ifdef EXAGO_ENABLE_IPOPT
      true,
#else
      false,
#endif
    },
    { 
      "hiop",
#ifdef EXAGO_ENABLE_HIOP
      true,
#else
      false,
#endif
    },
    { 
      "hiop_distributed",
#ifdef EXAGO_ENABLE_HIOP_DISTRIBUTED
      true,
#else
      false,
#endif
    },
    { 
      "hiop_sparse",
#ifdef EXAGO_ENABLE_HIOP_SPARSE
      true,
#else
      false,
#endif
    },
    { 
      "gpu",
#ifdef EXAGO_ENABLE_GPU
      true,
#else
      false,
#endif
    },
    { 
      "raja",
#ifdef EXAGO_ENABLE_RAJA
      true,
#else
      false,
#endif
    },
  };
}

PetscErrorCode ExaGOVersionGetFullVersionInfo(char** str)
{
  *str = (char*)malloc(2048);
  strcat(*str, "ExaGO version ");
  strcat(*str, EXAGO_VERSION);
  strcat(*str, " released on ");
  strcat(*str, EXAGO_RELEASE_DATE);
  strcat(*str, "\nbuilt with:\n");

  for(const auto& dep : ExaGODependencies)
  {
    const auto& name = dep.first;
    const auto& is_enabled = dep.second;
    char buf[1024];
    sprintf(buf, "\t%-20s%20s\n", name.c_str(), is_enabled?"YES":"NO");
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

const std::unordered_map<std::string, bool>& ExaGOGetDependencies()
{
  return ExaGODependencies;
}
