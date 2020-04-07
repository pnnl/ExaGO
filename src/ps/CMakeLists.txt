
# Specify PowerSystem module sources (shared across libraries)
set(PS_SRC
  ps.c
  psreaddata.c
  psislanding.c
)

# Compile library
add_library(ps SHARED ${PS_SRC})

target_link_libraries(ps PRIVATE
			 PETSC::SNES 
                         MPI::MPI_C
			 ${SCOPFLOW_UTIL_LIBS})

# Install shared objects
install(TARGETS ps DESTINATION lib)

set(SCOPFLOW_PS_LIBS ${SCOPFLOW_PS_LIBS} ps PARENT_SCOPE)
