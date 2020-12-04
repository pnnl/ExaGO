set(SOPFLOW_SOLVER_SRC
    solver/ipopt/sopflow-ipopt.c
    solver/ipopt/sopflow-ipopt2.c
    solver/empar/sopflow-empar.c
)

set(SOPFLOW_MODEL_SRC
  model/genramp/sopflow-genramp.c
)

set(SOPFLOW_SRC
	interface/sopflow.c
	interface/sopflowregi.c
	interface/sopflowoutput.c
	interface/sopflowscen.c
	${SOPFLOW_MODEL_SRC}
	${SOPFLOW_SOLVER_SRC}
)

# Header files for sopflow
set(SOPFLOW_INCLUDE
 ${CMAKE_SOURCE_DIR}/include/sopflow.h
)

#Compile library
add_library(sopflow ${SOPFLOW_SRC})

# link libraries
target_link_libraries(sopflow PRIVATE
			     PETSC::SNES
			     MPI::MPI_C
			     ${EXAGO_UTIL_LIBS}
			     ${EXAGO_PS_LIBS}
			     ${EXAGO_PFLOW_LIBS}
			     ${EXAGO_OPFLOW_LIBS}
)

if(EXAGO_HAVE_IPOPT)
	target_link_libraries(sopflow PUBLIC INTERFACE Ipopt)
	target_include_directories(sopflow PUBLIC ${IPOPT_INCLUDE_DIR})
endif()

if(EXAGO_HAVE_HIOP)
	target_link_libraries(sopflow PUBLIC INTERFACE HiOp)
	target_include_directories(sopflow PUBLIC ${HIOP_INCLUDE_DIR})
endif()

# Install library
install(TARGETS sopflow DESTINATION lib)

# Install header
install(FILES ${SOPFLOW_INCLUDE} DESTINATION include)

set(SOPFLOW_LIBS ${SOPFLOW_LIBS} sopflow PARENT_SCOPE)
