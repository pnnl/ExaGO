add_subdirectory(utils)
add_subdirectory(ps)
add_subdirectory(pflow)
add_subdirectory(opflow)
add_subdirectory(scopflow)
add_subdirectory(sopflow)
add_subdirectory(tcopflow)

# Add variables in the parent scope
set(EXAGO_PFLOW_LIBS
    ${EXAGO_PFLOW_LIBS} pflow
    PARENT_SCOPE
)
set(EXAGO_TCOPFLOW_LIBS
    ${EXAGO_TCOPFLOW_LIBS} tcopflow
    PARENT_SCOPE
)
