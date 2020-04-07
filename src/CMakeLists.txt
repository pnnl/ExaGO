
# Debugging messages
# message("I am in the src directory!")

# Specify PowerSystem module sources (shared across libraries)
set(POWERSYS_SRC
  ps/ps.c
  ps/psreaddata.c
  ps/psislanding.c
)

# Specify sources for utility functions shared across libraries
set(UTILS_SRC
  utils/utils.c
  utils/comm.c
)

# Specify sources for PowerFlow library
set(PFLOW_SRC
  pflow/pflow.c
  ${POWERSYS_SRC}
  ${UTILS_SRC}
)

# Compile pflow library
add_library(pflow SHARED ${PFLOW_SRC})

# Link pflow library
target_link_libraries(pflow PRIVATE PETSC::SNES MPI::MPI_C)

# Add pflow library to SCOPFLOW_PFLOW_LIBS variable in the parent scope
set(SCOPFLOW_PFLOW_LIBS ${SCOPFLOW_PFLOW_LIBS} pflow PARENT_SCOPE)

# Header files for pflow library to install
set(PFLOW_INCLUDE
  ../include/pflow.h
)

# Install shared objects
install(TARGETS pflow DESTINATION lib)

# Install headers
install(FILES ${PFLOW_INCLUDE} DESTINATION include)
