set(TCOPFLOW_SOLVER_SRC solver/ipopt/tcopflow_ipopt.cpp)

set(TCOPFLOW_MODEL_SRC model/genramp/genramp.cpp)

set(TCOPFLOW_SRC
    interface/tcopflow.cpp interface/tcopflowregi.cpp
    interface/tcopflowoutput.cpp interface/tcopflowreadprofiles.cpp
    ${TCOPFLOW_MODEL_SRC} ${TCOPFLOW_SOLVER_SRC}
)

# Header files for tcopflow
set(TCOPFLOW_INCLUDE ${CMAKE_SOURCE_DIR}/include/tcopflow.h)

set_source_files_properties(${TCOPFLOW_SRC} PROPERTIES LANGUAGE CXX)

exago_add_library(
  TCOPFLOW
  SOURCES
  ${TCOPFLOW_SRC}
  HEADERS
  ${TCOPFLOW_INCLUDE}
  LINK_LIBRARIES
  PUBLIC
  ExaGO::OPFLOW
  PRIVATE
  $<$<BOOL:${EXAGO_ENABLE_IPOPT}>:Ipopt>
  OUTPUT_NAME
  exago_tcopflow
)
