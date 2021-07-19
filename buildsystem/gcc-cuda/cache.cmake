message(STATUS "Loading CMake cache for a GCC+CUDA+MPI build")

set(prefix ${CMAKE_SOURCE_DIR}/install)
message(STATUS "Setting initial installation prefix to ${prefix}")
set(CMAKE_INSTALL_PREFIX
    ${prefix}
    CACHE PATH ""
)

set(arch 60)
message(STATUS "Using initial cuda architecture of ${arch}")
set(CMAKE_CUDA_ARCHITECTURES
    ${arch}
    CACHE STRING ""
)

message(STATUS "Building both shared and static libraries")
set(EXAGO_BUILD_SHARED
    ON
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

message(STATUS "Enabling GPU, HiOp, MPI, PETSC, and RAJA")
set(EXAGO_ENABLE_GPU
    ON
    CACHE BOOL ""
)
set(EXAGO_ENABLE_CUDA
    ON
    CACHE BOOL ""
)
set(EXAGO_ENABLE_HIOP
    ON
    CACHE BOOL ""
)
set(EXAGO_ENABLE_IPOPT
    OFF
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

message(STATUS "Enablying Python")
set(EXAGO_ENABLE_PYTHON
    ON
    CACHE BOOL ""
)

message(STATUS "Done setting initial CMake cache")
