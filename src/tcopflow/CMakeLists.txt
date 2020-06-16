set(TCOPFLOW_SOLVER_SRC
    solver/ipopt/tcopflow-ipopt.c
)

set(TCOPFLOW_MODEL_SRC
  model/genramp/genramp.c
)

set(TCOPFLOW_SRC
	interface/tcopflow.c
	interface/tcopflowregi.c
	interface/tcopflowoutput.c
	interface/tcopflowreadprofiles.c
	${TCOPFLOW_MODEL_SRC}
	${TCOPFLOW_SOLVER_SRC}
)

# Header files for tcopflow
set(TCOPFLOW_INCLUDE
 ${CMAKE_SOURCE_DIR}/include/tcopflow.h
)

#Compile library
add_library(tcopflow ${TCOPFLOW_SRC})

# link libraries
target_link_libraries(tcopflow PRIVATE
			     PETSC::SNES
			     MPI::MPI_C
			     ${SCOPFLOW_UTIL_LIBS}
			     ${SCOPFLOW_PS_LIBS}
			     ${SCOPFLOW_PFLOW_LIBS}
			     ${SCOPFLOW_OPFLOW_LIBS}
)

if(SCOPFLOW_HAVE_IPOPT)
	target_link_libraries(tcopflow PUBLIC INTERFACE Ipopt)
	target_include_directories(tcopflow PUBLIC ${IPOPT_INCLUDE_DIR})
endif()

if(SCOPFLOW_HAVE_HIOP)
	target_link_libraries(tcopflow PUBLIC INTERFACE HiOp)
	target_include_directories(tcopflow PUBLIC ${HIOP_INCLUDE_DIR})
endif()

# Install library
install(TARGETS tcopflow DESTINATION lib)

# Install header
install(FILES ${TCOPFLOW_INCLUDE} DESTINATION include)

set(TCOPFLOW_LIBS ${TCOPFLOW_LIBS} tcopflow PARENT_SCOPE)
