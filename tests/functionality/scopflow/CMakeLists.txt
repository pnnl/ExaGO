#[[ Functionality tests for ExaGO SCOPFLOW applications ]]

add_executable(test_scopflow_functionality selfcheck.cpp)
target_include_directories(
  test_scopflow_functionality PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..
)

#[[ Configure SCOPFLOW functionality tests ]]
target_link_libraries(test_scopflow_functionality ExaGO::SCOPFLOW)
if(EXAGO_INSTALL_TESTS)
  install(TARGETS test_scopflow_functionality
          RUNTIME DESTINATION tests/functionality
  )
  # SCOPFLOW Functionality Tests
  message(STATUS "Configuring SCOPFLOW Functionality Tests.")

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SCOPFLOW_IPOPT_CONT_TESTSUITE
    DEPENDS
    IPOPT
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:test_scopflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/ipopt_cont.toml
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SCOPFLOW_EMPAR_CONT_TESTSUITE
    DEPENDS
    IPOPT
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:test_scopflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/empar_cont.toml
  )

  if(EXAGO_ENABLE_ALL_TESTS)
    exago_add_test(
      NAME
      FUNCTIONALITY_TEST_SCOPFLOW_IPOPT_MULTIPERIOD_TESTSUITE
      DEPENDS
      IPOPT
      COMMAND
      ${RUNCMD}
      $<TARGET_FILE:test_scopflow_functionality>
      ${CMAKE_CURRENT_SOURCE_DIR}/ipopt_multiperiod.toml
    )
  endif()

  # Tests for HIOP pridecomp
  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SCOPFLOW_HIOP_TESTSUITE
    DEPENDS
    IPOPT
    HIOP
    HIOP_DISTRIBUTED
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:test_scopflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/hiop_cont.toml
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SCOPFLOW_HIOP_MPI_TESTSUITE
    DEPENDS
    IPOPT
    HIOP
    HIOP_DISTRIBUTED
    COMMAND
    ${MPICMD}
    "-n"
    3
    $<TARGET_FILE:test_scopflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/hiop_cont_mpi.toml
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SCOPFLOW_HIOP_SERIAL_TESTSUITE
    DEPENDS
    IPOPT
    HIOP
    HIOP_DISTRIBUTED
    COMMAND
    ${MPICMD}
    "-n"
    3
    $<TARGET_FILE:test_scopflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/hiop_serial.toml
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SCOPFLOW_HIOP_RAJA_TESTSUITE
    DEPENDS
    IPOPT
    HIOP
    HIOP_DISTRIBUTED
    COMMAND
    ${MPICMD}
    "-n"
    3
    $<TARGET_FILE:test_scopflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/hiop_raja.toml
  )
endif()
