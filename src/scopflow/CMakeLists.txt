set(SCOPFLOW_SOLVER_SRC
    solver/ipopt/scopflow-ipopt.c
    solver/ipopt/scopflow-ipopt2.c
    solver/empar/scopflow-empar.c
)

set(SCOPFLOW_MODEL_SRC
  model/genramp/genramp.c
  model/genrampt/genrampt.c
)

set(SCOPFLOW_SRC
	interface/scopflow.c
	interface/scopflowregi.c
	interface/scopflowoutput.c
	interface/scopflowctgc.c
	${SCOPFLOW_MODEL_SRC}
	${SCOPFLOW_SOLVER_SRC}
)

# Header files for scopflow
set(SCOPFLOW_INCLUDE
 ${CMAKE_SOURCE_DIR}/include/scopflow.h
)

#Compile library
add_library(scopflow ${SCOPFLOW_SRC})

# link libraries
target_link_libraries(scopflow PRIVATE
			     PETSC::SNES
			     ${EXAGO_UTIL_LIBS}
			     ${EXAGO_PS_LIBS}
			     ${EXAGO_PFLOW_LIBS}
			     ${EXAGO_OPFLOW_LIBS}
			     ${EXAGO_TCOPFLOW_LIBS}
)

if(EXAGO_HAVE_MPI)
  target_link_libraries(scopflow PRIVATE MPI::MPI_C)
endif()

if(EXAGO_HAVE_IPOPT)
	target_link_libraries(scopflow PUBLIC INTERFACE Ipopt)
	target_include_directories(scopflow PUBLIC ${IPOPT_INCLUDE_DIR})
endif()

# Install library
install(TARGETS scopflow DESTINATION lib)

# Install header
install(FILES ${SCOPFLOW_INCLUDE} DESTINATION include)

set(EXAGO_SCOPFLOW_LIBS ${EXAGO_SCOPFLOW_LIBS} scopflow PARENT_SCOPE)
