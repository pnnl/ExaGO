
add_subdirectory(utils)
add_subdirectory(ps)
add_subdirectory(pflow)
add_subdirectory(opflow)
add_subdirectory(scopflow)

# Add variables in the parent scope
set(SCOPFLOW_UTIL_LIBS ${SCOPFLOW_UTIL_LIBS} PARENT_SCOPE)
set(SCOPFLOW_PS_LIBS ${SCOPFLOW_PS_LIBS} PARENT_SCOPE)
set(SCOPFLOW_PFLOW_LIBS ${SCOPFLOW_PFLOW_LIBS} pflow PARENT_SCOPE)
set(SCOPFLOW_OPFLOW_LIBS ${SCOPFLOW_OPFLOW_LIBS} opflow PARENT_SCOPE)
set(SCOPFLOW_LIBS ${SCOPFLOW_LIBS} scopflow PARENT_SCOPE)
