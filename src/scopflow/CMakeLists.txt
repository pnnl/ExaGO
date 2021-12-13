set(SCOPFLOW_SOLVER_SRC
    solver/ipopt/scopflow_ipopt.cpp solver/empar/scopflow_empar.cpp
    solver/hiop/scopflow_hiop.cpp
)

set(SCOPFLOW_MODEL_SRC model/genramp/genramp.cpp model/genrampt/genrampt.cpp)

set(SCOPFLOW_SRC
    interface/scopflow.cpp interface/scopflowregi.cpp
    interface/scopflowoutput.cpp ${SCOPFLOW_MODEL_SRC} ${SCOPFLOW_SOLVER_SRC}
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
