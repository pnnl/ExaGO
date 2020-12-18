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
			     ${EXAGO_UTIL_LIBS}
			     ${EXAGO_PS_LIBS}
			     ${EXAGO_PFLOW_LIBS}
			     ${EXAGO_OPFLOW_LIBS}
)

if(EXAGO_HAVE_MPI)
  target_link_libraries(tcopflow PRIVATE MPI::MPI_C)
endif()

if(EXAGO_HAVE_IPOPT)
	target_link_libraries(tcopflow PUBLIC INTERFACE Ipopt)
	target_include_directories(tcopflow PUBLIC ${IPOPT_INCLUDE_DIR})
endif()

# Install library
install(TARGETS tcopflow DESTINATION lib)

# Install header
install(FILES ${TCOPFLOW_INCLUDE} DESTINATION include)

set(EXAGO_TCOPFLOW_LIBS ${EXAGO_TCOPFLOW_LIBS} tcopflow PARENT_SCOPE)
