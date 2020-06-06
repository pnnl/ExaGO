# Models
set(OPFLOW_FORM_SRC
    model/current-bal-cartesian/ibcar.c
    model/current-bal-cartesian/ibcar2.c
    model/power-bal-polar/pbpol.c
    model/power-bal-cartesian/pbcar.c
)

set(OPFLOW_SOLVER_SRC
    solver/tao/opflow-tao.c
    solver/ipopt/opflow-ipopt.c
    solver/hiop/opflow-hiop.cpp
)

set(OPFLOW_SRC
	interface/opflow.c
	interface/opflowregi.c
	interface/opflowoutput.c
	${OPFLOW_FORM_SRC}
	${OPFLOW_SOLVER_SRC}
)


# Header files for opflow
set(OPFLOW_INCLUDE
 ${CMAKE_SOURCE_DIR}/include/opflow.h
)

#Compile library
add_library(opflow ${OPFLOW_SRC})

# link libraries
target_link_libraries(opflow PRIVATE
			     PETSC::ALL
					 MPI::MPI_C
			     ${SCOPFLOW_UTIL_LIBS}
			     ${SCOPFLOW_PS_LIBS}
			     ${SCOPFLOW_PFLOW_LIBS}
           ${MAGMA_LIBRARY}
			     -lblas -llapack
)

if(SCOPFLOW_HAVE_IPOPT)
	target_link_libraries(opflow PUBLIC INTERFACE Ipopt)
	target_include_directories(opflow PUBLIC ${IPOPT_INCLUDE_DIR})
endif()

if(SCOPFLOW_HAVE_HIOP)
	target_link_libraries(opflow PUBLIC INTERFACE HiOp)
	target_include_directories(opflow PUBLIC ${HIOP_INCLUDE_DIR})
endif()

# Install shared object
install(TARGETS opflow DESTINATION lib)

# Install header
install(FILES ${OPFLOW_INCLUDE} DESTINATION include)

# Install options file
if(SCOPFLOW_HAVE_HIOP)
  # Need to figure out where should options file go eventually, putting
  # it in bin because HIOP needs to have the options file in the same directory
  # as executable
  install(FILES ${CMAKE_SOURCE_DIR}/hiop.options DESTINATION bin)
  install(FILES ${CMAKE_SOURCE_DIR}/hiop.options DESTINATION ${PROJECT_BINARY_DIR})
endif()

set(SCOPFLOW_OPFLOW_LIBS ${SCOPFLOW_OPFLOW_LIBS} opflow PARENT_SCOPE)
