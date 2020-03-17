# Debugging messages
# message("I am in the applications directory!")

# Define sources for pflow-proto application
set(PFLOWPROTO_SRC
  pflow-proto.c
)

# Compile pflowproto application
add_executable(pflowproto ${PFLOWPROTO_SRC})

# Link pflowproto application
target_link_libraries(pflowproto 
  ${SCOPFLOW_PFLOW_LIBS}
  ${SCOPFLOW_MATH_LIB}
  PETSC::SNES
  MPI::MPI_C
)

# Install the executable
install(TARGETS pflowproto RUNTIME DESTINATION bin)
