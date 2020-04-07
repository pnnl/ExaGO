
add_subdirectory(utils)
add_subdirectory(ps)
add_subdirectory(pflow)
add_subdirectory(opflow)

# Add variables in the parent scope
set(SCOPFLOW_UTIL_LIBS ${SCOPFLOW_UTIL_LIBS} PARENT_SCOPE)
set(SCOPFLOW_PS_LIBS ${SCOPFLOW_PS_LIBS} PARENT_SCOPE)
set(SCOPFLOW_PFLOW_LIBS ${SCOPFLOW_PFLOW_LIBS} pflow PARENT_SCOPE)
set(SCOPFLOW_OPFLOW_LIBS ${SCOPFLOW_OPFLOW_LIBS} opflow PARENT_SCOPE)

# Header files for pflow
set(PFLOW_INCLUDE
  ../include/pflow.h
)

# Header files for opflow
set(OPFLOW_INCLUDE
 ../include/opflow.h
)

# Install headers
install(FILES ${PFLOW_INCLUDE} DESTINATION include)
install(FILES ${OPFLOW_INCLUDE} DESTINATION include)
