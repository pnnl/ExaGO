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
target_compile_options(UTILS PUBLIC "$<$<CONFIG:DEBUG>:-Wall;-Wextra;-Wpedantic;-Werror>")
