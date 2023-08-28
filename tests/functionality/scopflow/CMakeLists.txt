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
  set_tests_properties(
    FUNCTIONALITY_TEST_SCOPFLOW_EMPAR_CONT_TESTSUITE PROPERTIES LABELS
                                                                "incline-skip"
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
    ${CMAKE_CURRENT_SOURCE_DIR}/hiop_cont_pbpol.toml
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SCOPFLOW_HIOP_PBHIOP_TESTSUITE
    DEPENDS
    IPOPT
    HIOP
    HIOP_DISTRIBUTED
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:test_scopflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/hiop_cont_pbhiop.toml
    -hiop_mem_space
    HOST
  )
  set_tests_properties(
    FUNCTIONALITY_TEST_SCOPFLOW_HIOP_PBHIOP_TESTSUITE PROPERTIES LABELS
                                                              "newell-skip"
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SCOPFLOW_HIOP_PBHIOP_MPI_TESTSUITE
    DEPENDS
    IPOPT
    HIOP
    HIOP_DISTRIBUTED
    COMMAND
    ${MPICMD}
    "-n"
    3
    $<TARGET_FILE:test_scopflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/hiop_cont_mpi_pbhiop.toml
    -hiop_mem_space
    HOST
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
  set_tests_properties(
    FUNCTIONALITY_TEST_SCOPFLOW_HIOP_MPI_TESTSUITE PROPERTIES LABELS
                                                              "incline-skip"
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
  set_tests_properties(
    FUNCTIONALITY_TEST_SCOPFLOW_HIOP_SERIAL_TESTSUITE PROPERTIES LABELS
                                                                 "incline-skip"
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
  set_tests_properties(
    FUNCTIONALITY_TEST_SCOPFLOW_HIOP_RAJA_TESTSUITE PROPERTIES LABELS
                                                               "incline-skip"
  )

  foreach(fmt "MATPOWER" "CSV" "JSON" "MINIMAL")
    set(tname "SCOPFLOW_SOLUTION_OUTPUT_${fmt}")
    exago_add_test(
      NAME
      "${tname}"
      COMMAND
      ${RUNCMD}
      $<TARGET_FILE:app_scopflow>
      -netfile
      ${CMAKE_SOURCE_DIR}/datafiles/case9/case9mod.m
      -ctgcfile
      ${CMAKE_SOURCE_DIR}/datafiles/case9/case9.cont
      -scopflow_solver
      IPOPT
      -scopflow_Nc
      -1
      -save_output
      1
      -opflow_output_format
      "${fmt}"
      -scopflow_output_directory
      "scopflowout_${fmt}"
      DEPENDS
      IPOPT
    )
  endforeach()

endif()
