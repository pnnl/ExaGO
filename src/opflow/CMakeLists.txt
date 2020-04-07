# Formulations
set(OPFLOW_FORM_SRC
    formulation/current-bal-cartesian/ibcar.c
    formulation/current-bal-cartesian/ibcar2.c
    formulation/power-bal-polar/pbpol.c
    formulation/power-bal-cartesian/pbcar.c
)

set(OPFLOW_SOLVER_SRC
    solver/tao/opflow-tao.c
)

set(OPFLOW_SRC
	interface/opflow.c
	interface/opflowregi.c
	${OPFLOW_FORM_SRC}
	${OPFLOW_SOLVER_SRC}
)

#Compile library
add_library(opflow SHARED ${OPFLOW_SRC})

# link libraries
target_link_libraries(opflow PRIVATE
			     PETSC::SNES
			     MPI::MPI_C
			     ${SCOPFLOW_UTIL_LIBS}
			     ${SCOPFLOW_PS_LIBS}
			     ${SCOPFLOW_PFLOW_LIBS}
)

# Install shared objects
install(TARGETS opflow DESTINATION lib)

set(SCOPFLOW_OPFLOW_LIBS ${SCOPFLOW_OPFLOW_LIBS} opflow PARENT_SCOPE)
