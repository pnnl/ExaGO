# Models
set(OPFLOW_FORM_SRC
    model/current-bal-cartesian/ibcar.c model/current-bal-cartesian/ibcar2.c
    model/power-bal-polar/pbpol.c model/power-bal-cartesian/pbcar.c
    model/power-bal-hiop/pbpolhiop.c model/power-bal-hiop/pbpolhiopkernels.c
)

# Build RAJA/Umpire modules if RAJA and Umpire are enabled
if(EXAGO_ENABLE_RAJA)
  set(OPFLOW_FORM_SRC ${OPFLOW_FORM_SRC} model/power-bal-hiop/pbpolrajahiop.cpp
                      model/power-bal-hiop/pbpolrajahiopkernels.cpp
  )
endif()

set(OPFLOW_SOLVER_SRC
    solver/tao/opflow-tao.c solver/ipopt/opflow-ipopt.c
    solver/hiop/opflow-hiop.cpp solver/hiop/opflow-hiopsparse.cpp
)

set(OPFLOW_SRC interface/opflow.c interface/opflowregi.c
               interface/opflowoutput.c ${OPFLOW_FORM_SRC} ${OPFLOW_SOLVER_SRC}
)

set_source_files_properties(${OPFLOW_SRC} PROPERTIES LANGUAGE CXX)

if(EXAGO_ENABLE_RAJA AND EXAGO_ENABLE_CUDA)
  set_source_files_properties(
    model/power-bal-hiop/pbpolrajahiopkernels.cpp PROPERTIES LANGUAGE CUDA
  )
endif()

# Header files for opflow
set(OPFLOW_INCLUDE ${CMAKE_SOURCE_DIR}/include/opflow.h)

exago_add_library(
  OPFLOW
  SOURCES
  ${OPFLOW_SRC}
  HEADERS
  ${OPFLOW_INCLUDE}
  LINK_LIBRARIES
  PUBLIC
  $<$<BOOL:${EXAGO_ENABLE_RAJA}>:umpire>
  PUBLIC
  $<$<BOOL:${EXAGO_ENABLE_RAJA}>:RAJA>
  PUBLIC
  $<$<BOOL:${EXAGO_ENABLE_OMP}>:OpenMP::OpenMP_CXX>
  $<$<BOOL:${EXAGO_ENABLE_HIOP}>:HiOp>
  PUBLIC
  $<$<BOOL:${EXAGO_ENABLE_IPOPT}>:Ipopt>
  PUBLIC
  ExaGO::PFLOW
  PUBLIC
  PETSC::ALL
  ${LAPACK_LIRBARIES}
  OUTPUT_NAME
  exago_opflow
)
