message(STATUS "Loading CMake cache for a GCC+HIP+MPI build")

set(prefix ${CMAKE_SOURCE_DIR}/install)
message(STATUS "Setting initial installation prefix to ${prefix}")
set(CMAKE_INSTALL_PREFIX
    ${prefix}
    CACHE PATH ""
)

set(arch gfx90a)
message(STATUS "Using initial AMD GPU target architecture of ${arch}")
set(AMDGPU_TARGETS
    ${arch}
    CACHE STRING ""
)

message(STATUS "Just building static libraries for now.")
set(EXAGO_BUILD_SHARED
    OFF
    CACHE BOOL ""
)
set(EXAGO_BUILD_STATIC
    ON
    CACHE BOOL ""
)

message(STATUS "Building in debug mode")
set(CMAKE_BUILD_TYPE
    Debug
    CACHE STRING ""
)

message(STATUS "Enabling GPU (HIP), HiOp, MPI, PETSC, and RAJA")
set(EXAGO_ENABLE_GPU
    ON
    CACHE BOOL ""
)
set(EXAGO_ENABLE_CUDA
    OFF
    CACHE BOOL ""
)
set(EXAGO_ENABLE_HIP
    ON
    CACHE BOOL ""
)
set(EXAGO_ENABLE_HIOP
    ON
    CACHE BOOL ""
)
set(EXAGO_ENABLE_IPOPT
    ON
    CACHE BOOL ""
)
set(EXAGO_ENABLE_MPI
    ON
    CACHE BOOL ""
)
set(EXAGO_ENABLE_PETSC
    ON
    CACHE BOOL ""
)
set(EXAGO_ENABLE_RAJA
    ON
    CACHE BOOL ""
)
set(PETSC_DIR
    $ENV{PETSC_DIR}
    CACHE PATH ""
)
set(EXAGO_RUN_TESTS
    ON
    CACHE BOOL ""
)

message(STATUS "Disabling Python when building with Ipopt on Crusher")
set(EXAGO_ENABLE_PYTHON
    OFF
    CACHE BOOL ""
)

message(STATUS "Enabling Logging")
set(EXAGO_ENABLE_LOGGING
    ON
    CACHE BOOL ""
)

message(STATUS "Done setting initial CMake cache")
