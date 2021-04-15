#[[ Functionality tests for ExaGO applications ]]

set(EXAGO_APP_LIBS
    ExaGO::SOPFLOW
    ExaGO::SCOPFLOW
    ${EXAGO_TCOPFLOW_LIBS}
    ExaGO::OPFLOW
    ${EXAGO_PFLOW_LIBS}
    ExaGO::PS
    ExaGO::UTILS
    ${EXAGO_MATH_LIB}
    PETSC::SNES
)

if(EXAGO_ENABLE_MPI)
  list(APPEND EXAGO_APP_LIBS MPI::MPI_C)
endif()

add_subdirectory(opflow)
add_subdirectory(pflow)
add_subdirectory(tcopflow)
add_subdirectory(scopflow)
add_subdirectory(sopflow)
