#[[ Functionality tests for ExaGO OPFLOW applications ]]

add_executable(test_opflow_functionality selfcheck.cpp)

target_include_directories(
  test_opflow_functionality PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
                                   ${CMAKE_CURRENT_SOURCE_DIR}/..
)

#[[ Configure OPFLOW functionality tests ]]
target_link_libraries(test_opflow_functionality ExaGO::OPFLOW)
if(EXAGO_INSTALL_TESTS)

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_OPFLOW_HIOP_TOML_TESTSUITE
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:test_opflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/hiop_pbpolhiop.toml
    DEPENDS
    HIOP
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_OPFLOW_RAJAHIOP_TOML_TESTSUITE
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:test_opflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/hiop_pbpolrajahiop.toml
    DEPENDS
    HIOP
    RAJA
  )

  if(EXAGO_ENABLE_HIOP_SPARSE)
    exago_add_test(
      NAME
      FUNCTIONALITY_TEST_OPFLOW_HIOPSPARSE_TOML_TESTSUITE
      COMMAND
      ${RUNCMD}
      $<TARGET_FILE:test_opflow_functionality>
      ${CMAKE_CURRENT_SOURCE_DIR}/hiopsparse_pbpol.toml
      DEPENDS
      HIOP_SPARSE
    )
  endif()

  # Note: All cartesian and/or current balance models are way behind. These
  # models should not be used and their tests should be disabled.
  if(EXAGO_FEATURE_CARTESIAN_ACTIVE)
    exago_add_test(
      NAME
      FUNCTIONALITY_TEST_OPFLOW_IPOPT_CARTESIAN_TOML_TESTSUITE
      COMMAND
      ${RUNCMD}
      $<TARGET_FILE:test_opflow_functionality>
      ${CMAKE_CURRENT_SOURCE_DIR}/ipopt_pbcart.toml
      DEPENDS
      IPOPT
    )
  endif()

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_OPFLOW_IPOPT_POLAR_TOML_TESTSUITE
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:test_opflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/ipopt_pbpol.toml
    DEPENDS
    IPOPT
  )

  install(TARGETS test_opflow_functionality
          RUNTIME DESTINATION tests/functionality
  )

  set(test_suites hiop_pbpolhiop.toml hiop_pbpolrajahiop.toml
                  hiopsparse_pbpol.toml ipopt_pbpol.toml ipopt_pbcart.toml
  )

  install(FILES ${test_suites} DESTINATION share/exago/testsuites)

endif()
