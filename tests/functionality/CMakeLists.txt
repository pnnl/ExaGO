#[[ Functionality tests for ExaGO applications ]]

add_subdirectory(utils)

set(EXAGO_APP_LIBS
  ${EXAGO_UTIL_LIBS}
  ${EXAGO_PS_LIBS}
  ${EXAGO_PFLOW_LIBS}
  ${EXAGO_OPFLOW_LIBS}
	${EXAGO_SOPFLOW_LIBS}
	${EXAGO_SCOPFLOW_LIBS}
	${EXAGO_TCOPFLOW_LIBS}
  ${EXAGO_MATH_LIB}
  testFunctionality
  PETSC::SNES
)

if(EXAGO_HAVE_MPI)
  list(APPEND EXAGO_APP_LIBS MPI::MPI_C)
endif()

add_executable(testPflowFunctionality    pflow.c)
add_executable(testOpflowFunctionality   opflow.c)
add_executable(testScopflowFunctionality scopflow.c)
add_executable(testSopflowFunctionality  sopflow.c)
add_executable(testTcopflowFunctionality tcopflow.c)

#[[ Configure functionality tests ]]
foreach(target testPflowFunctionality 
    testOpflowFunctionality 
    testScopflowFunctionality 
    testSopflowFunctionality 
    testTcopflowFunctionality)
  target_link_libraries(${target} ${EXAGO_APP_LIBS})
  install(TARGETS ${target} RUNTIME DESTINATION bin)
endforeach()
