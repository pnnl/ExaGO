# Specify PowerSystem module sources (shared across libraries)
set(PFLOW_SRC
  pflow.c
)

# Compile library
add_library(pflow SHARED ${PFLOW_SRC})

target_link_libraries(pflow PRIVATE
			 PETSC::SNES 
                         MPI::MPI_C
			 ${SCOPFLOW_UTIL_LIBS}
			 ${SCOPFLOW_PS_LIBS})

# Install shared objects
install(TARGETS pflow DESTINATION lib)

set(SCOPFLOW_PFLOW_LIBS ${SCOPFLOW_PFLOW_LIBS} pflow PARENT_SCOPE)
