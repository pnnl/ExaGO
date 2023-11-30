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
    FUNCTIONALITY_TEST_SOPFLOW_SCENARIO3_PBHIOP_TOML
    DEPENDS
    IPOPT
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:test_sopflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/sopflow_multiscenario3_pbhiop.toml
    -opflow_ignore_lineflow_constraints
    -hiop_mem_space
    HOST
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SOPFLOW_SCENARIO3_RAJA_TOML
    DEPENDS
    IPOPT
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:test_sopflow_functionality>
    ${CMAKE_CURRENT_SOURCE_DIR}/sopflow_multiscenario3_raja.toml
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
    -hiop_mem_space
    HOST
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
  if(EXAGO_ENABLE_IPOPT)
    set_tests_properties(
      FUNCTIONALITY_TEST_SOPFLOW_SCENARIO_RAJA_GPU_TOML
      PROPERTIES LABELS "incline-skip"
    )
  endif()

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

  foreach(fmt "MATPOWER" "CSV" "JSON" "MINIMAL")
    foreach(flat 1 0)
      if(flat)
        set(tname "SOPFLOW_SOLUTION_OUTPUT_${fmt}_FLATTENED")
        set(outdir "sopflowout_${fmt}_flat")
      else()
        set(tname "SOPFLOW_SOLUTION_OUTPUT_${fmt}")
        set(outdir "sopflowout_${fmt}")
      endif()
      exago_add_test(
        NAME
        "${tname}"
        COMMAND
        ${RUNCMD}
        $<TARGET_FILE:app_sopflow>
        -netfile
        ${CMAKE_SOURCE_DIR}/datafiles/case9/case9mod.m
        -ctgcfile
        ${CMAKE_SOURCE_DIR}/datafiles/case9/case9.cont
        -scenfile
        ${CMAKE_SOURCE_DIR}/datafiles/case9/10_scenarios_9bus.csv
        -sopflow_solver
        IPOPT
        -sopflow_Ns
        8
        -sopflow_Nc
        3
        -sopflow_enable_multicontingency
        1
        -sopflow_flatten_contingencies
        ${flat}
        -save_output
        1
        -opflow_output_format
        "${fmt}"
        -sopflow_output_directory
        "${outdir}"
        DEPENDS
        IPOPT
      )
    endforeach()
  endforeach()

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
