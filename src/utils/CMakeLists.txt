# Sources for utility functions shared across libraries

# Specify sources for utility functions shared across libraries
set(UTILS_SRC utils.cpp version.cpp scopflowutils.cpp comm.cpp)

# Headers
set(UTILS_HEADERS ${CMAKE_SOURCE_DIR}/include/common.h
                  ${CMAKE_SOURCE_DIR}/include/utils.h
)

set_source_files_properties(${UTILS_SRC} PROPERTIES LANGUAGE CXX)

exago_add_library(
  UTILS
  SOURCES
  ${UTILS_SRC}
  HEADERS
  ${UTILS_HEADERS}
  LINK_LIBRARIES
  PUBLIC
  PkgConfig::PETSC
  PUBLIC
  $<$<BOOL:${EXAGO_ENABLE_LOGGING}>:spdlog::spdlog>
  PUBLIC
  $<$<BOOL:${EXAGO_ENABLE_MPI}>:MPI::MPI_C>
  PUBLIC
  $<$<BOOL:${EXAGO_ENABLE_MPI}>:MPI::MPI_CXX>
  OUTPUT_NAME
  exago_utils
)

# Going to have this target supply/propogate compiler flags and see if that
# works. -Wpedantic for CUDA compile options added a significant amount of new
# errors, so leaving out for now. -Werror is problematic on MacOS
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "AppleClang")
  set(CUDA_COMPILE_OPTIONS
      "$<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CONFIG:DEBUG>>:SHELL:--compiler-options -Wall,-Wextra>"
  )
  set(CXX_COMPILE_OPTIONS
      "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:DEBUG>>:-Wall;-Wextra;-Wpedantic>"
  )
else()
  set(CUDA_COMPILE_OPTIONS
      "$<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CONFIG:DEBUG>>:SHELL:--compiler-options -Werror,-Wall,-Wextra>"
  )
  set(CXX_COMPILE_OPTIONS
      "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:DEBUG>>:-Werror;-Wall;-Wextra;-Wpedantic>"
  )
endif()

if(EXAGO_BUILD_SHARED)
  target_compile_options(
    UTILS_shared PUBLIC ${CUDA_COMPILE_OPTIONS} ${CXX_COMPILE_OPTIONS}
  )
endif()
if(EXAGO_BUILD_STATIC)
  target_compile_options(
    UTILS_static PUBLIC ${CUDA_COMPILE_OPTIONS} ${CXX_COMPILE_OPTIONS}
  )
endif()
