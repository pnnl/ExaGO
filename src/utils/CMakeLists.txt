
# Sources for utility functions shared across libraries

# Specify sources for utility functions shared across libraries
set(UTILS_SRC
  utils.c
  comm.c
)

# Add library
add_library(util SHARED ${UTILS_SRC})

# Link libraries
target_link_libraries(util PRIVATE PETSC::SNES MPI::MPI_C)

# Install shared objects
install(TARGETS util DESTINATION lib)

set(SCOPFLOW_UTIL_LIBS ${SCOPFLOW_UTIL_LIBS} util PARENT_SCOPE)
