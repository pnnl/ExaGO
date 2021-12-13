# Sources for power flow functions shared across libraries

# Specify PowerSystem module sources (shared across libraries)
set(PFLOW_SRC pflow.cpp)

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
  OUTPUT_NAME
  exago_pflow
)
