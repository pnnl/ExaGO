
# Specify PowerSystem module sources (shared across libraries)
set(PS_SRC
  ps.c
  psreaddata.c
  psislanding.c
)

set(PS_INCLUDE
  ${CMAKE_SOURCE_DIR}/include/ps.h
)

# Compile library
add_library(ps SHARED ${PS_SRC})

target_link_libraries(ps PRIVATE
			 PETSC::SNES 
                         MPI::MPI_C
			 ${SCOPFLOW_UTIL_LIBS})

# Install library
install(TARGETS ps DESTINATION lib)

# Install header
install(FILES ${PS_INCLUDE} DESTINATION include)

set(SCOPFLOW_PS_LIBS ${SCOPFLOW_PS_LIBS} ps PARENT_SCOPE)
