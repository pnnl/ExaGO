# Specify PowerSystem module sources (shared across libraries)
set(PFLOW_SRC
  pflow.c
)

# Header files for pflow
set(PFLOW_INCLUDE
  ${CMAKE_SOURCE_DIR}/include/pflow.h
)

# Compile library
add_library(pflow ${PFLOW_SRC})

target_link_libraries(pflow PRIVATE
			 PETSC::SNES 
                         MPI::MPI_C
			 ${EXAGO_UTIL_LIBS}
			 ${EXAGO_PS_LIBS})

# Install libraries
install(TARGETS pflow DESTINATION lib)

#Install headers
install(FILES ${PFLOW_INCLUDE} DESTINATION include)

set(EXAGO_PFLOW_LIBS ${EXAGO_PFLOW_LIBS} pflow PARENT_SCOPE)
