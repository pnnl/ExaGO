set(SCOPFLOW_SOLVER_SRC solver/ipopt/scopflow-ipopt.c
                        solver/empar/scopflow-empar.c
)

set(SCOPFLOW_MODEL_SRC model/genramp/genramp.c model/genrampt/genrampt.c)

set(SCOPFLOW_SRC
    interface/scopflow.c interface/scopflowregi.c interface/scopflowoutput.c
    interface/scopflowctgc.c ${SCOPFLOW_MODEL_SRC} ${SCOPFLOW_SOLVER_SRC}
)

# Header files for scopflow
set(SCOPFLOW_INCLUDE ${CMAKE_SOURCE_DIR}/include/scopflow.h)

set_source_files_properties(${SCOPFLOW_SRC} PROPERTIES LANGUAGE CXX)

exago_add_library(
  SCOPFLOW
  SOURCES
  ${SCOPFLOW_SRC}
  HEADERS
  ${SCOPFLOW_INCLUDE}
  LINK_LIBRARIES
  PUBLIC
  PETSC::SNES
  ExaGO::UTILS
  ExaGO::PS
  ${EXAGO_PFLOW_LIBS}
  ExaGO::OPFLOW
  $<$<BOOL:${EXAGO_ENABLE_IPOPT}>:Ipopt>
  PRIVATE
  $<$<BOOL:${EXAGO_ENABLE_MPI}>:MPI::MPI_C>
  OUTPUT_NAME
  exago_scopflow
)
