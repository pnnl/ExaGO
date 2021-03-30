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

# Compile library
add_library(sopflow ${SOPFLOW_SRC})

# link libraries
target_link_libraries(
  sopflow
  PRIVATE PETSC::SNES
          ${EXAGO_UTIL_LIBS}
          ${EXAGO_PS_LIBS}
          ${EXAGO_PFLOW_LIBS}
          ${EXAGO_OPFLOW_LIBS}
          ${EXAGO_TCOPFLOW_LIBS}
          ${EXAGO_SCOPFLOW_LIBS}
)

if(EXAGO_ENABLE_MPI)
  target_link_libraries(sopflow PRIVATE MPI::MPI_C)
endif()

if(EXAGO_ENABLE_IPOPT)
  target_link_libraries(
    sopflow
    PUBLIC
    INTERFACE Ipopt
  )
  target_include_directories(sopflow PUBLIC ${IPOPT_INCLUDE_DIR})
endif()

# Install library
install(TARGETS sopflow DESTINATION lib)

# Install header
install(FILES ${SOPFLOW_INCLUDE} DESTINATION include)

set(EXAGO_SOPFLOW_LIBS
    ${EXAGO_SOPFLOW_LIBS} sopflow
    PARENT_SCOPE
)
