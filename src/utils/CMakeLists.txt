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

# Going to have this target supply/propogate compiler flags and see if that works
set(COMPILE_OPTIONS "$<$<CONFIG:DEBUG>:-Wall;-Wextra;-Wpedantic;-Werror>")
if(EXAGO_BUILD_SHARED)
  target_compile_options(UTILS_shared PUBLIC ${COMPILE_OPTIONS})
endif()
if(EXAGO_BUILD_STATIC)
  target_compile_options(UTILS_static PUBLIC ${COMPILE_OPTIONS})
endif()

