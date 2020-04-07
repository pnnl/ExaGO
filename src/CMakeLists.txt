
add_subdirectory(utils)
add_subdirectory(ps)
add_subdirectory(pflow)

# Add variables in the parent scope
set(SCOPFLOW_UTIL_LIBS ${SCOPFLOW_UTIL_LIBS} PARENT_SCOPE)
set(SCOPFLOW_PS_LIBS ${SCOPFLOW_PS_LIBS} PARENT_SCOPE)
set(SCOPFLOW_PFLOW_LIBS ${SCOPFLOW_PFLOW_LIBS} pflow PARENT_SCOPE)

# Header files for pflow library to install
set(PFLOW_INCLUDE
  ../include/pflow.h
)

# Install headers
install(FILES ${PFLOW_INCLUDE} DESTINATION include)
