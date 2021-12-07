set(SOPFLOW_SOLVER_SRC
    solver/ipopt/sopflow-ipopt.c solver/empar/sopflow-empar.c
    solver/hiop/sopflow-hiop.cpp
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
  PUBLIC
  ExaGO::SCOPFLOW
  PRIVATE
  $<$<BOOL:${EXAGO_ENABLE_IPOPT}>:Ipopt>
  OUTPUT_NAME
  exago_sopflow
)
