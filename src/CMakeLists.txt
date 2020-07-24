
add_subdirectory(utils)
add_subdirectory(ps)
add_subdirectory(pflow)
add_subdirectory(opflow)
add_subdirectory(scopflow)
add_subdirectory(tcopflow)

# Add variables in the parent scope
set(EXAGO_UTIL_LIBS ${EXAGO_UTIL_LIBS} PARENT_SCOPE)
set(EXAGO_PS_LIBS ${EXAGO_PS_LIBS} PARENT_SCOPE)
set(EXAGO_PFLOW_LIBS ${EXAGO_PFLOW_LIBS} pflow PARENT_SCOPE)
set(EXAGO_OPFLOW_LIBS ${EXAGO_OPFLOW_LIBS} opflow PARENT_SCOPE)
set(SCOPFLOW_LIBS ${SCOPFLOW_LIBS} scopflow PARENT_SCOPE)
set(TCOPFLOW_LIBS ${TCOPFLOW_LIBS} tcopflow PARENT_SCOPE)

