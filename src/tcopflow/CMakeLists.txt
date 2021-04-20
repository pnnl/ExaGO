set(TCOPFLOW_SOLVER_SRC solver/ipopt/tcopflow-ipopt.c)

set(TCOPFLOW_MODEL_SRC model/genramp/genramp.c)

set(TCOPFLOW_SRC
    interface/tcopflow.c interface/tcopflowregi.c interface/tcopflowoutput.c
    interface/tcopflowreadprofiles.c ${TCOPFLOW_MODEL_SRC}
    ${TCOPFLOW_SOLVER_SRC}
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
  PRIVATE
  PETSC::ALL
  ExaGO::UTILS
  ExaGO::PS
  ${EXAGO_PFLOW_LIBS}
  ExaGO::OPFLOW
  ExaGO::SCOPFLOW
  $<$<BOOL:${EXAGO_ENABLE_MPI}>:MPI::MPI_C>
  PUBLIC
  $<$<BOOL:${EXAGO_ENABLE_IPOPT}>:Ipopt>
  OUTPUT_NAME
  exago_tcopflow
)
