#include <common.h>
#include <exago_config.h>
#include <unordered_map>
#include <utils.h>
#if (EXAGO_ENABLE_IPOPT || EXAGO_ENABLE_HIOP)
#include <opflow.h>
#include <scopflow.h>
#include <sopflow.h>
#endif
#include <version.h>

namespace {
static const std::unordered_map<std::string, bool> ExaGOBoolConfigOptions = {
    {
        "EXAGO_BUILD_SHARED",
#ifdef EXAGO_BUILD_SHARED
        true,
#else
        false,
#endif
    },
    {
        "EXAGO_BUILD_STATIC",
#ifdef EXAGO_BUILD_STATIC
        true,
#else
        false,
#endif
    },
    {
        "EXAGO_ENABLE_PETSC",
#ifdef EXAGO_ENABLE_PETSC
        true,
#else
        false,
#endif
    },
    {
        "EXAGO_ENABLE_MPI",
#ifdef EXAGO_ENABLE_MPI
        true,
#else
        false,
#endif
    },
    {
        "EXAGO_ENABLE_IPOPT",
#ifdef EXAGO_ENABLE_IPOPT
        true,
#else
        false,
#endif
    },
    {
        "EXAGO_ENABLE_HIOP",
#ifdef EXAGO_ENABLE_HIOP
        true,
#else
        false,
#endif
    },
    {
        "EXAGO_ENABLE_HIOP_DISTRIBUTED",
#ifdef EXAGO_ENABLE_HIOP_DISTRIBUTED
        true,
#else
        false,
#endif
    },
    {
        "EXAGO_ENABLE_HIOP_SPARSE",
#ifdef EXAGO_ENABLE_HIOP_SPARSE
        true,
#else
        false,
#endif
    },
    {
        "EXAGO_ENABLE_GPU",
#ifdef EXAGO_ENABLE_GPU
        true,
#else
        false,
#endif
    },
    {
        "EXAGO_ENABLE_CUDA",
#ifdef EXAGO_ENABLE_CUDA
        true,
#else
        false,
#endif
    },
    {
        "EXAGO_ENABLE_HIP",
#ifdef EXAGO_ENABLE_HIP
        true,
#else
        false,
#endif
    },
    {
        "EXAGO_ENABLE_PYTHON",
#ifdef EXAGO_ENABLE_PYTHON
        true,
#else
        false,
#endif
    },
    {
        "EXAGO_ENABLE_RAJA",
#ifdef EXAGO_ENABLE_RAJA
        true,
#else
        false,
#endif
    },
};
static const std::unordered_map<std::string, std::string>
    ExaGOStringConfigOptions = {
        {"EXAGO_OPTIONS_DIR", EXAGO_OPTIONS_DIR},
        {"CMAKE_INSTALL_PREFIX", CMAKE_INSTALL_PREFIX},
        {"CMAKE_BUILD_TYPE", CMAKE_BUILD_TYPE},
        {"CMAKE_C_COMPILER", CMAKE_C_COMPILER},
        {"CMAKE_CXX_COMPILER", CMAKE_CXX_COMPILER},
        {"PETSC_DIR", EXAGO_PETSC_DIR},
        {"PETSC_INCLUDES", EXAGO_PETSC_INCLUDES},
        {"PETSC_LIBRARIES", EXAGO_PETSC_LIBRARIES},
#ifdef EXAGO_ENABLE_HIOP
        {"HiOp_DIR", EXAGO_HiOp_DIR},
#endif
#ifdef EXAGO_ENABLE_IPOPT
        {"IPOPT_ROOT_DIR", EXAGO_IPOPT_ROOT_DIR},
        {"IPOPT_INCLUDES", EXAGO_IPOPT_INCLUDES},
        {"IPOPT_LIBRARIES", EXAGO_IPOPT_LIBRARIES},
#endif
};
} // namespace

PetscErrorCode ExaGOVersionGetFullVersionInfo(std::string &str) {
  str = fmt::format("ExaGO version {} built on {}\n", EXAGO_VERSION, __DATE__);

  /* Recreate configure command */
  str +=
      fmt::format("Built with the following command:\n$ cmake -B {} -S {} \\\n",
                  CMAKE_BUILD_DIR, CMAKE_SOURCE_DIR);
  for (auto const &conf_opt : ExaGOBoolConfigOptions) {
    const auto &name = conf_opt.first;
    const auto &is_enabled = conf_opt.second;
    str +=
        fmt::format("\t-D{}:BOOL={} \\\n", name, (is_enabled ? "ON" : "OFF"));
  }
  auto confit = ExaGOStringConfigOptions.begin();
  while (confit != ExaGOStringConfigOptions.end()) {
    const auto &name = (*confit).first;
    const auto &value = (*confit).second;
    /* If it's not the final option, add a newline escape */
    str +=
        fmt::format("\t-D{}:STRING=\"{}\" {}\n", name, value,
                    (++confit == ExaGOStringConfigOptions.end() ? "" : "\\"));
  }
  return 0;
}

PetscErrorCode ExaGOVersionGetVersion(int *major, int *minor, int *patch) {
  *major = atoi(EXAGO_VERSION_MAJOR);
  *minor = atoi(EXAGO_VERSION_MINOR);
  *patch = atoi(EXAGO_VERSION_PATCH);
  return 0;
}

PetscErrorCode ExaGOVersionGetVersionStr(std::string &str) {
  str = EXAGO_VERSION;
  return 0;
}

const std::unordered_map<std::string, bool> &ExaGOGetBoolConfigOptions() {
  return ExaGOBoolConfigOptions;
}
const std::unordered_map<std::string, std::string> &
ExaGOGetStringConfigOptions() {
  return ExaGOStringConfigOptions;
}
