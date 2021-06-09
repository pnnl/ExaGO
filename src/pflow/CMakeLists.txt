# Sources for power flow functions shared across libraries

# Specify PowerSystem module sources (shared across libraries)
set(PFLOW_SRC pflow.c)

# Header files for pflow
set(PFLOW_INCLUDE ${CMAKE_SOURCE_DIR}/include/pflow.h)

set_source_files_properties(${PFLOW_SRC} PROPERTIES LANGUAGE CXX)

exago_add_library(
  PFLOW
  SOURCES
  ${PFLOW_SRC}
  HEADERS
  ${PFLOW_INCLUDE}
  LINK_LIBRARIES
  PUBLIC
  ExaGO::UTILS
  PUBLIC
  ExaGO::PS
  PUBLIC
  PETSC::ALL
  INTERFACE
  $<$<BOOL:${EXAGO_ENABLE_MPI}>:MPI::MPI_C>
  OUTPUT_NAME
  exago_pflow
)
