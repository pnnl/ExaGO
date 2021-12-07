set(SCOPFLOW_SOLVER_SRC
    solver/ipopt/scopflow-ipopt.c solver/empar/scopflow-empar.c
    solver/hiop/scopflow-hiop.cpp
)

set(SCOPFLOW_MODEL_SRC model/genramp/genramp.c model/genrampt/genrampt.c)

set(SCOPFLOW_SRC
    interface/scopflow.c interface/scopflowregi.c interface/scopflowoutput.c
    ${SCOPFLOW_MODEL_SRC} ${SCOPFLOW_SOLVER_SRC}
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
  ExaGO::CTGCLIST
  PUBLIC
  ExaGO::TCOPFLOW
  PRIVATE
  $<$<BOOL:${EXAGO_ENABLE_IPOPT}>:Ipopt>
  OUTPUT_NAME
  exago_scopflow
)
