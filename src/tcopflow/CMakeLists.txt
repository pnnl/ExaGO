set(TCOPFLOW_SOLVER_SRC
    solver/ipopt/tcopflow-ipopt.c
)

set(TCOPFLOW_SRC
	interface/tcopflow.c
	interface/tcopflowregi.c
	interface/tcopflowoutput.c
	interface/tcopflowreadprofiles.c
	${TCOPFLOW_SOLVER_SRC}
)

# Header files for tcopflow
set(TCOPFLOW_INCLUDE
 ${CMAKE_SOURCE_DIR}/include/tcopflow.h
)

#Compile library
add_library(tcopflow SHARED ${TCOPFLOW_SRC})

# link libraries
target_link_libraries(tcopflow PRIVATE
			     PETSC::SNES
			     MPI::MPI_C
			     ${IPOPT_LIBRARY}
			     ${SCOPFLOW_UTIL_LIBS}
			     ${SCOPFLOW_PS_LIBS}
			     ${SCOPFLOW_PFLOW_LIBS}
			     ${SCOPFLOW_OPFLOW_LIBS}
)

# Install library
install(TARGETS tcopflow DESTINATION lib)

# Install header
install(FILES ${TCOPFLOW_INCLUDE} DESTINATION include)

set(TCOPFLOW_LIBS ${TCOPFLOW_LIBS} tcopflow PARENT_SCOPE)
