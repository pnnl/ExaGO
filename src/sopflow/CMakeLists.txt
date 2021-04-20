set(SOPFLOW_SOLVER_SRC solver/ipopt/sopflow-ipopt.c
                       solver/empar/sopflow-empar.c
)

set(SOPFLOW_MODEL_SRC model/genramp/sopflow-genramp.c
                      model/genrampc/sopflow-genrampc.c
)

set(SOPFLOW_SRC
    interface/sopflow.c interface/sopflowregi.c interface/sopflowoutput.c
    interface/sopflowscen.c ${SOPFLOW_MODEL_SRC} ${SOPFLOW_SOLVER_SRC}
)

# Header files for sopflow
set(SOPFLOW_INCLUDE ${CMAKE_SOURCE_DIR}/include/sopflow.h)

set_source_files_properties(${SOPFLOW_SRC} PROPERTIES LANGUAGE CXX)

exago_add_library(
  SOPFLOW
  SOURCES
  ${SOPFLOW_SRC}
  HEADERS
  ${SOPFLOW_INCLUDE}
  LINK_LIBRARIES
  PRIVATE
  PETSC::SNES
  ExaGO::UTILS
  ExaGO::PS
  ${EXAGO_PFLOW_LIBS}
  ExaGO::OPFLOW
  ExaGO::TCOPFLOW
  ${EXAGO_SCOPFLOW_LIBS}
  $<$<BOOL:${EXAGO_ENABLE_MPI}>:MPI::MPI_C>
  PUBLIC
  $<$<BOOL:${EXAGO_ENABLE_IPOPT}>:Ipopt>
  OUTPUT_NAME
  exago_sopflow
)
