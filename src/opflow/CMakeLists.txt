# Models
set(OPFLOW_FORM_SRC
    model/current_bal_cartesian/ibcar.cpp
    model/current_bal_cartesian/ibcar2.cpp
    model/power_bal_polar/pbpol.cpp
    model/power_bal_cartesian/pbcar.cpp
    model/power_bal_hiop/pbpolhiop.cpp
    model/power_bal_hiop/pbpolhiopkernels.cpp
    model/dcopf/dcopf.cpp
)

# Build RAJA/Umpire modules if RAJA and Umpire are enabled
if(EXAGO_ENABLE_RAJA)
  set(OPFLOW_FORM_SRC ${OPFLOW_FORM_SRC} model/power_bal_hiop/pbpolrajahiop.cpp
                      model/power_bal_hiop/pbpolrajahiopkernels.cpp
  )
endif()

set(OPFLOW_SOLVER_SRC
    solver/tao/opflow_tao.cpp solver/ipopt/opflow_ipopt.cpp
    solver/hiop/opflow_hiop.cpp solver/hiop/opflow_hiopsparse.cpp
)

set(OPFLOW_SRC
    interface/opflow.cpp interface/opflowregi.cpp interface/opflowoutput.cpp
    ${OPFLOW_FORM_SRC} ${OPFLOW_SOLVER_SRC}
)

set_source_files_properties(${OPFLOW_SRC} PROPERTIES LANGUAGE CXX)

if(EXAGO_ENABLE_RAJA AND EXAGO_ENABLE_CUDA)
  set_source_files_properties(
    model/power_bal_hiop/pbpolrajahiopkernels.cpp PROPERTIES LANGUAGE CUDA
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
  PUBLIC
  $<$<BOOL:${EXAGO_ENABLE_HIOP}>:HiOp::HiOp>
  PUBLIC
  $<$<BOOL:${EXAGO_ENABLE_IPOPT}>:Ipopt>
  PUBLIC
  ExaGO::PFLOW
  PUBLIC
  ${LAPACK_LIBRARIES}
  OUTPUT_NAME
  exago_opflow
)
