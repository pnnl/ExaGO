#[[ Functionality tests for ExaGO SOPFLOW applications ]]

set_source_files_properties(selfcheck.cpp PROPERTIES LANGUAGE CXX)
add_executable(testSopflowFunctionalityTOML selfcheck.cpp)
target_include_directories(
  testSopflowFunctionalityTOML PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                       ${CMAKE_CURRENT_SOURCE_DIR}/..
)

set(tolerance "1e-3")

#[[ Configure SOPFLOW functionality tests ]]
target_link_libraries(testSopflowFunctionalityTOML ExaGO::SOPFLOW)
if(EXAGO_INSTALL_TESTS)
  message(STATUS "Configuring SOPFLOW functionality tests")

  install(TARGETS testSopflowFunctionalityTOML
          RUNTIME DESTINATION tests/functionality
  )

  if(EXAGO_ENABLE_ALL_TESTS)
  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SOPFLOW_SCENARIO_TOML
    DEPENDS
    IPOPT
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:testSopflowFunctionalityTOML>
    ${CMAKE_CURRENT_SOURCE_DIR}/sopflow_multiscenario.toml
    -no_optfile
  )
  endif()

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SOPFLOW_CONTINGENCY_TOML
    DEPENDS
    IPOPT
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:testSopflowFunctionalityTOML>
    ${CMAKE_CURRENT_SOURCE_DIR}/sopflow_multicontingency.toml
    -no_optfile
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SOPFLOW_CONTINGENCY_MULTIPERIOD_TOML
    DEPENDS
    IPOPT
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:testSopflowFunctionalityTOML>
    ${CMAKE_CURRENT_SOURCE_DIR}/sopflow_multiperiod.toml
    -no_optfile
  )

endif()
