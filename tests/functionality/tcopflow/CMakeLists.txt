#[[ Functionality tests for ExaGO TCOPFLOW applications ]]

set_source_files_properties(
  tcopflow.c tcopflowselfcheck.c PROPERTIES LANGUAGE CXX
)

add_executable(testTcopflowFunctionality tcopflow.c tcopflowselfcheck.c)
target_include_directories(
  testTcopflowFunctionality PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

#[[ Configure TCOPFLOW functionality tests ]]
target_link_libraries(testTcopflowFunctionality ${EXAGO_APP_LIBS})
if(EXAGO_INSTALL_TESTS)
  install(TARGETS testTcopflowFunctionality
          RUNTIME DESTINATION tests/functionality
  )

  message(STATUS "Configuring TCOPFLOW functionality tests")
  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_TCOPFLOW
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:testTcopflowFunctionality>
    -netfile
    ${CMAKE_SOURCE_DIR}/datafiles/case9/case9mod_gen3_wind.m
    -tcopflow_windgenprofile
    ${CMAKE_SOURCE_DIR}/datafiles/case9/scenarios_9bus.csv
    -tcopflow_qloadprofile
    ${CMAKE_SOURCE_DIR}/datafiles/case9/load_Q.csv
    -tcopflow_ploadprofile
    ${CMAKE_SOURCE_DIR}/datafiles/case9/load_P.csv
    -tcopflow_iscoupling
    0
    -tcopflow_dT
    5.0
    -tcopflow_duration
    0.5
    -tcopflow_tolerance
    ${tolerance}
    -opflow_ignore_lineflow_constraints
    0
    -opflow_initialization
    ACPF
    -opflow_genbusvoltage
    VARIABLE_WITHIN_BOUNDS
    DEPENDS
    IPOPT
  )
endif()
