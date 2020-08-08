# Models
set(OPFLOW_FORM_SRC
    model/current-bal-cartesian/ibcar.c
    model/current-bal-cartesian/ibcar2.c
    model/power-bal-polar/pbpol.c
    model/power-bal-cartesian/pbcar.c
    model/power-bal-polar2/pbpol2.c
    model/power-bal-polar2/pbpol2params.c
    model/power-bal-polar2/pbpol2kernels.c
    model/power-bal-hiop/pbpolhiop.c
    model/power-bal-hiop/pbpolhiopkernels.c
)

set(OPFLOW_SOLVER_SRC
    solver/tao/opflow-tao.c
    solver/ipopt/opflow-ipopt.c
    solver/hiop/opflow-hiop.cpp
    solver/hiop/opflow-hiopnew.cpp
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
			     ${EXAGO_UTIL_LIBS}
			     ${EXAGO_PS_LIBS}
			     ${EXAGO_PFLOW_LIBS}
			     -lblas -llapack
)

# If Ipopt interface is enabled, link to Ipopt libs
if(EXAGO_HAVE_IPOPT)
	target_link_libraries(opflow PUBLIC INTERFACE Ipopt)
	target_include_directories(opflow PUBLIC ${IPOPT_INCLUDE_DIR})
endif()

# If HiOp interface is enabled, link to HiOp libs
if(EXAGO_HAVE_HIOP)
	target_link_libraries(opflow PUBLIC INTERFACE HiOp)
  target_include_directories(opflow PUBLIC ${HIOP_INCLUDE_DIR})
  # If HiOp is using GPU add CUDA includes
  if(EXAGO_ENABLE_GPU)
    target_include_directories(opflow PUBLIC ${CUDA_TOOLKIT_INCLUDE})
  endif()
endif()

# Install shared object
install(TARGETS opflow DESTINATION lib)

# Install header
install(FILES ${OPFLOW_INCLUDE} DESTINATION include)

# Install options file
if(EXAGO_HAVE_IPOPT)
  # Need to figure out where should options file go eventually, putting
  # it in bin because IPOPT needs to have the options file in the same directory
  # as executable
  install(FILES ${CMAKE_SOURCE_DIR}/ipopt.opt DESTINATION bin)
  install(FILES ${CMAKE_SOURCE_DIR}/ipopt.opt DESTINATION ${PROJECT_BINARY_DIR})
  file(COPY ${CMAKE_SOURCE_DIR}/ipopt.opt DESTINATION bin)
  file(COPY ${CMAKE_SOURCE_DIR}/ipopt.opt DESTINATION ${PROJECT_BINARY_DIR})
endif()

if(EXAGO_HAVE_HIOP)
  # Need to figure out where should options file go eventually, putting
  # it in bin because HIOP needs to have the options file in the same directory
  # as executable
  install(FILES ${CMAKE_SOURCE_DIR}/hiop.options DESTINATION bin)
  install(FILES ${CMAKE_SOURCE_DIR}/hiop.options DESTINATION ${PROJECT_BINARY_DIR})
  file(COPY ${CMAKE_SOURCE_DIR}/hiop.options DESTINATION bin)
  file(COPY ${CMAKE_SOURCE_DIR}/hiop.options DESTINATION ${PROJECT_BINARY_DIR})
endif()

set(EXAGO_OPFLOW_LIBS ${EXAGO_OPFLOW_LIBS} opflow PARENT_SCOPE)
