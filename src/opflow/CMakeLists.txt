# Models
set(OPFLOW_FORM_SRC
    model/current-bal-cartesian/ibcar.cpp
    model/current-bal-cartesian/ibcar2.cpp
    model/power-bal-polar/pbpol.cpp
    model/power-bal-cartesian/pbcar.cpp
    model/power-bal-hiop/pbpolhiop.cpp
    model/power-bal-hiop/pbpolhiopkernels.cpp
)

# Build RAJA/Umpire modules if RAJA and Umpire are enabled
if(EXAGO_HAVE_RAJA)
  set(OPFLOW_FORM_SRC ${OPFLOW_FORM_SRC} 
    model/power-bal-hiop/pbpolrajahiop.cpp
    model/power-bal-hiop/pbpolrajahiopkernels.cpp
  )
endif()

if(EXAGO_HAVE_RAJA AND EXAGO_ENABLE_GPU)
  set_source_files_properties(model/power-bal-hiop/pbpolrajahiopkernels.cpp PROPERTIES LANGUAGE CUDA)
endif()


set(OPFLOW_SOLVER_SRC
    solver/tao/opflow-tao.cpp
    solver/ipopt/opflow-ipopt.cpp
    solver/hiop/opflow-hiop.cpp
)

set(OPFLOW_SRC
	interface/opflow.cpp
	interface/opflowregi.cpp
	interface/opflowoutput.cpp
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
			     ${EXAGO_UTIL_LIBS}
			     ${EXAGO_PS_LIBS}
			     ${EXAGO_PFLOW_LIBS}
			     ${LAPACK_LIBRARIES}
)

if(EXAGO_HAVE_MPI)
  target_link_libraries(opflow PRIVATE MPI::MPI_C)
endif()

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

if(EXAGO_HAVE_RAJA)
  target_link_libraries(opflow PUBLIC umpire RAJA OpenMP::OpenMP_CXX)
  target_include_directories(opflow PUBLIC ${UMPIRE_INCLUDE_DIRS} ${camp_INSTALL_PREFIX}/include)
endif()


# Install shared object
install(TARGETS opflow DESTINATION lib)

# Install header
install(FILES ${OPFLOW_INCLUDE} DESTINATION include)

set(EXAGO_OPFLOW_LIBS ${EXAGO_OPFLOW_LIBS} opflow PARENT_SCOPE)