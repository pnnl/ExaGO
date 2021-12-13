# Sources for utility functions shared across libraries

# Specify sources for utility functions shared across libraries
set(UTILS_SRC utils.cpp version.cpp scopflowutils.cpp comm.cpp)

# Headers
set(UTILS_HEADERS ${CMAKE_SOURCE_DIR}/include/common.h
                  ${CMAKE_SOURCE_DIR}/include/utils.hpp
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
  PETSC::SNES
  PUBLIC
  $<$<BOOL:${EXAGO_ENABLE_MPI}>:MPI::MPI_C>
  PUBLIC
  $<$<BOOL:${EXAGO_ENABLE_MPI}>:MPI::MPI_CXX>
  OUTPUT_NAME
  exago_utils
)
