set(SOPFLOW_SOLVER_SRC
    solver/ipopt/sopflow_ipopt.cpp solver/empar/sopflow_empar.cpp
    solver/hiop/sopflow_hiop.cpp
)

set(SOPFLOW_MODEL_SRC model/genramp/sopflow_genramp.cpp
                      model/genrampc/sopflow_genrampc.cpp
)

set(SOPFLOW_SRC
    interface/sopflow.cpp interface/sopflowregi.cpp interface/sopflowoutput.cpp
    interface/sopflowscen.cpp ${SOPFLOW_MODEL_SRC} ${SOPFLOW_SOLVER_SRC}
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
