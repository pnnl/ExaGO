set(SCOPFLOW_SOLVER_SRC
    solver/ipopt/scopflow-ipopt.c
    solver/empar/scopflow-empar.c
)

set(SCOPFLOW_SRC
	interface/scopflow.c
	interface/scopflowregi.c
	interface/scopflowoutput.c
	interface/scopflowctgc.c
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
			     MPI::MPI_C
			     ${SCOPFLOW_UTIL_LIBS}
			     ${SCOPFLOW_PS_LIBS}
			     ${SCOPFLOW_PFLOW_LIBS}
			     ${SCOPFLOW_OPFLOW_LIBS}
)

if(SCOPFLOW_HAVE_IPOPT)
	target_link_libraries(scopflow PUBLIC INTERFACE Ipopt)
	target_include_directories(scopflow PUBLIC ${IPOPT_INCLUDE_DIR})
endif()

if(SCOPFLOW_HAVE_HIOP)
	target_link_libraries(scopflow PUBLIC INTERFACE HiOp)
	target_include_directories(scopflow PUBLIC ${HIOP_INCLUDE_DIR})
endif()

# Install library
install(TARGETS scopflow DESTINATION lib)

# Install header
install(FILES ${SCOPFLOW_INCLUDE} DESTINATION include)

set(SCOPFLOW_LIBS ${SCOPFLOW_LIBS} scopflow PARENT_SCOPE)
