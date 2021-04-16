# Sources for power flow functions shared across libraries

# Specify PowerSystem module sources (shared across libraries)
set(PFLOW_SRC pflow.c)

# Header files for pflow
set(PFLOW_INCLUDE ${CMAKE_SOURCE_DIR}/include/pflow.h)

exago_add_library(
  PFLOW
  SOURCES
  ${PFLOW_SRC}
  HEADERS
  ${PFLOW_INCLUDE}
  LINK_LIBRARIES
  PRIVATE
  ExaGO::UTILS
  PRIVATE
  ExaGO::PS
  PRIVATE
  PETSC::SNES
  $<$<BOOL:EXAGO_ENABLE_MPI>:MPI::MPI_C>
  OUTPUT_NAME
  exago_pflow
)
