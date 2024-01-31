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

  if(EXAGO_ENABLE_GPU)
    exago_add_test(
      NAME
      FUNCTIONALITY_TEST_OPFLOW_RAJAHIOP_SPARSE_GPU_TOML_TESTSUITE
      COMMAND
      ${RUNCMD}
      $<TARGET_FILE:test_opflow_functionality>
      ${CMAKE_CURRENT_SOURCE_DIR}/hiop_pbpolrajahiopsparse_gpu.toml
      DEPENDS
      HIOP
      RAJA
    )
    if(EXAGO_ENABLE_HIOP AND EXAGO_ENABLE_RAJA)
      # This test passes on Ascent...
      set_tests_properties(
        FUNCTIONALITY_TEST_OPFLOW_RAJAHIOP_SPARSE_GPU_TOML_TESTSUITE
        PROPERTIES LABELS "incline-skip;"
      )
    endif()
  endif(EXAGO_ENABLE_GPU)

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

  exago_add_test(
    NAME
    OPFLOW_MONITOR_KVLEVEL_IPOPT
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:app_opflow>
    -netfile
    ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.m
    -opflow_solver
    IPOPT
    -opflow_model
    POWER_BALANCE_POLAR
    -opflow_initialization
    ACPF
    -opflow_monitor_line_kvlevels
    115.0
    DEPENDS
    IPOPT
    HIOP
  )

  exago_add_test(
    NAME
    OPFLOW_LINEFLOW_PENALTY
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:app_opflow>
    -netfile
    ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.m
    -opflow_solver
    IPOPT
    -opflow_model
    POWER_BALANCE_POLAR
    -opflow_initialization
    ACPF
    -opflow_allow_lineflow_violation
    -opflow_lineflowviolation_penalty
    100
    -no_optfile
    DEPENDS
    IPOPT
    HIOP
  )

  exago_add_test(
    NAME
    OPFLOW_MONITOR_KVLEVEL_HIOP_CPU
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:app_opflow>
    -netfile
    ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.m
    -opflow_solver
    HIOP
    -opflow_model
    POWER_BALANCE_HIOP
    -hiop_compute_mode
    CPU
    -hiop_mem_space
    HOST
    -hiop_verbosity_level
    3
    -opflow_initialization
    ACPF
    -opflow_monitor_line_kvlevels
    115.0
    DEPENDS
    IPOPT
    HIOP
  )

  foreach(fmt "MATPOWER" "CSV" "JSON" "MINIMAL")
    set(tname "OPFLOW_SOLUTION_OUTPUT_${fmt}")
    exago_add_test(
      NAME
      "${tname}"
      COMMAND
      ${RUNCMD}
      $<TARGET_FILE:app_opflow>
      -netfile
      ${CMAKE_SOURCE_DIR}/datafiles/case9/case9mod.m
      -opflow_solver
      IPOPT
      -opflow_output_format
      "${fmt}"
      -save_output
      1
      DEPENDS
      IPOPT
    )
  endforeach()

  # exago_add_test( NAME OPFLOW_MONITOR_KVLEVEL_HIOP_GPU COMMAND ${RUNCMD}
  # $<TARGET_FILE:app_opflow> -netfile
  # ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.m -opflow_solver HIOP
  # -opflow_model PBPOLRAJAHIOP -hiop_compute_mode GPU -hiop_mem_space DEVICE
  # -hiop_verbosity_level 3 -opflow_initialization ACPF
  # -opflow_monitor_line_kvlevels 115.0 DEPENDS IPOPT HIOP )

  install(TARGETS test_opflow_functionality
          RUNTIME DESTINATION tests/functionality
  )

  set(test_suites hiop_pbpolhiop.toml hiop_pbpolrajahiop.toml
                  hiopsparse_pbpol.toml ipopt_pbpol.toml ipopt_pbcart.toml
  )

  install(FILES ${test_suites} DESTINATION share/exago/testsuites)

endif()
