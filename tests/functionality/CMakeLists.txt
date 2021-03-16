#[[ Functionality tests for ExaGO applications ]]

set(EXAGO_APP_LIBS
  ${EXAGO_SOPFLOW_LIBS}
  ${EXAGO_SCOPFLOW_LIBS}
  ${EXAGO_TCOPFLOW_LIBS}
  ${EXAGO_OPFLOW_LIBS}
  ${EXAGO_PFLOW_LIBS}
  ${EXAGO_PS_LIBS}
  ${EXAGO_UTIL_LIBS}
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
