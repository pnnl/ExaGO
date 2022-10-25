execute_process(COMMAND which gcc-9 OUTPUT_VARIABLE C_COMPILER)
execute_process(COMMAND which g++-9 OUTPUT_VARIABLE CXX_COMPILER)

string(STRIP ${C_COMPILER} C_COMPILER)
string(STRIP ${CXX_COMPILER} CXX_COMPILER)

message(STATUS "Using C_COMPILER ${C_COMPILER}")
message(STATUS "Using CXX_COMPILER ${CXX_COMPILER}")
set(CMAKE_C_COMPILER ${C_COMPILER} CACHE PATH "" FORCE)
set(CMAKE_CXX_COMPILER ${CXX_COMPILER} CACHE PATH "" FORCE)

message(STATUS "Building in Debug mode")
set(CMAKE_BUILD_TYPE
    Debug
    CACHE STRING ""
)

message(STATUS "Enabling HiOp, Ipopt and RAJA")
set(EXAGO_ENABLE_HIOP
    ON
    CACHE BOOL ""
)
set(EXAGO_ENABLE_IPOPT
    ON
    CACHE BOOL ""
)
set(EXAGO_ENABLE_RAJA
    ON
    CACHE BOOL ""
)

set(prefix ${CMAKE_SOURCE_DIR}/install)
message(STATUS "Setting initial installation prefix to ${prefix}")
set(CMAKE_INSTALL_PREFIX
    ${prefix}
    CACHE PATH ""
)

