set(SCOPFLOW_SOLVER_SRC
    solver/ipopt/scopflow-ipopt.c
)

set(SCOPFLOW_SRC
	interface/scopflow.c
	interface/scopflowregi.c
	${SCOPFLOW_SOLVER_SRC}
)

#Compile library
add_library(scopflow SHARED ${SCOPFLOW_SRC})

# link libraries
target_link_libraries(scopflow PRIVATE
			     PETSC::SNES
			     MPI::MPI_C
			     ${IPOPT_LIBRARY}
			     ${SCOPFLOW_UTIL_LIBS}
			     ${SCOPFLOW_PS_LIBS}
			     ${SCOPFLOW_PFLOW_LIBS}
			     ${SCOPFLOW_OPFLOW_LIBS}
)

# Install shared objects
install(TARGETS scopflow DESTINATION lib)

set(SCOPFLOW_LIBS ${SCOPFLOW_LIBS} scopflow PARENT_SCOPE)
