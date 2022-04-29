#[[ Functionality tests for ExaGO SOPFLOW applications ]]

set_source_files_properties(selfcheck.cpp PROPERTIES LANGUAGE CXX)
add_executable(test_sopflow_functionality selfcheck.cpp)
target_include_directories(
  test_sopflow_functionality PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
                                    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

set(tolerance "1e-3")

#[[ Configure SOPFLOW functionality tests ]]
target_link_libraries(test_sopflow_functionality ExaGO::SOPFLOW)
if(EXAGO_INSTALL_TESTS)
  message(STATUS "Configuring SOPFLOW functionality tests")

  install(TARGETS test_sopflow_functionality
          RUNTIME DESTINATION tests/functionality
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SOPFLOW_SCENARIO_TOML
    DEPENDS
    IPOPT
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:test_sopflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/sopflow_multiscenario.toml
    -no_optfile
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SOPFLOW_SCENARIO2_TOML
    DEPENDS
    IPOPT
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:test_sopflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/sopflow_multiscenario2.toml
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SOPFLOW_SCENARIO3_TOML
    DEPENDS
    IPOPT
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:test_sopflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/sopflow_multiscenario3.toml
    -opflow_ignore_lineflow_constraints
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SOPFLOW_CONTINGENCY_TOML
    DEPENDS
    IPOPT
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:test_sopflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/sopflow_multicontingency.toml
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SOPFLOW_SCENARIO_SERIAL_TOML
    DEPENDS
    IPOPT
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:test_sopflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/sopflow_multiscenario_serial.toml
    -no_optfile
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SOPFLOW_SCENARIO_PARALLEL_TOML
    DEPENDS
    IPOPT
    COMMAND
    ${MPICMD}
    "-n"
    3
    $<TARGET_FILE:test_sopflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/sopflow_multiscenario_parallel.toml
    -no_optfile
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SOPFLOW_SCENARIO_RAJA_GPU_TOML
    DEPENDS
    IPOPT
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:test_sopflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/sopflow_multiscenario_raja.toml
    -no_optfile
    -opflow_include_powerimbalance_variables
    -opflow_ignore_lineflow_constraints
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SOPFLOW_SCENARIO_MPI_RAJA_GPU_TOML
    DEPENDS
    IPOPT
    COMMAND
    ${MPICMD}
    "-n"
    3
    $<TARGET_FILE:test_sopflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/sopflow_multiscenario_mpi_raja.toml
    -no_optfile
    -opflow_include_powerimbalance_variables
    -opflow_ignore_lineflow_constraints
  )
endif()

# Disabled test for now, feature not ready yet
if(EXAGO_ENABLE_ALL_TESTS)
  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SOPFLOW_CONTINGENCY_MULTIPERIOD_TOML
    DEPENDS
    IPOPT
    COMMAND
    ${MPICMD}
    "-n"
    3
    $<TARGET_FILE:test_sopflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/sopflow_multiperiod.toml
  )
endif()
