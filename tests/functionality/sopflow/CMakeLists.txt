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

  if(EXAGO_INSTALL_TESTS)
    exago_add_test(
      NAME
      FUNCTIONALITY_TEST_SOPFLOW_SCENARIO_TOML
      DEPENDS
      IPOPT
      COMMAND
      ${RUNCMD}
      $<TARGET_FILE:test_sopflow_functionality>
      ${CMAKE_CURRENT_SOURCE_DIR}/sopflow_multiscenario.toml
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
  endif()

  # exago_add_test( NAME FUNCTIONALITY_TEST_SOPFLOW_SCENARIO_SERIAL_TOML DEPENDS
  # IPOPT COMMAND ${RUNCMD} $<TARGET_FILE:test_sopflow_functionality>
  # ${CMAKE_CURRENT_SOURCE_DIR}/sopflow_multiscenario_serial.toml )

  # exago_add_test( NAME FUNCTIONALITY_TEST_SOPFLOW_SCENARIO_PARALLEL_TOML
  # DEPENDS IPOPT COMMAND ${RUNCMD} $<TARGET_FILE:test_sopflow_functionality>
  # ${CMAKE_CURRENT_SOURCE_DIR}/sopflow_multiscenario_parallel.toml )

  # exago_add_test( NAME FUNCTIONALITY_TEST_SOPFLOW_SCENARIO_SERIAL_TOML DEPENDS
  # IPOPT COMMAND ${RUNCMD} $<TARGET_FILE:test_sopflow_functionality>
  # ${CMAKE_CURRENT_SOURCE_DIR}/sopflow_multiscenario_serial.toml -no_optfile
  # -sopflow_flatten_contingencies )

  # exago_add_test( NAME FUNCTIONALITY_TEST_SOPFLOW_SCENARIO_PARALLEL_TOML
  # DEPENDS IPOPT COMMAND ${RUNCMD} $<TARGET_FILE:test_sopflow_functionality>
  # ${CMAKE_CURRENT_SOURCE_DIR}/sopflow_multiscenario_parallel.toml -no_optfile
  # -sopflow_flatten_contingencies )

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

  # if (FALSE) exago_add_test( NAME SOPFLOW_200bus_HIOP_SERIAL COMMAND ${RUNCMD}
  # $<TARGET_FILE:app_sopflow> -netfile
  # ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.m -scenfile
  # ${CMAKE_SOURCE_DIR}/datafiles/10_scenarios_ACTIVSg200.csv -sopflow_solver
  # HIOP -sopflow_Ns -1 -opflow_initialization ACPF
  # -opflow_ignore_lineflow_constraints 1 -sopflow_subproblem_solver HIOP
  # -sopflow_subproblem_model POWER_BALANCE_HIOP -hiop_compute_mode CPU
  # -hiop_verbosity_level 3 DEPENDS IPOPT HIOP )

  # exago_add_test( NAME SOPFLOW_200bus_HIOP_MPI COMMAND ${MPICMD} "-n" 3
  # $<TARGET_FILE:app_sopflow> -netfile
  # ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.m -scenfile
  # ${CMAKE_SOURCE_DIR}/datafiles/10_scenarios_ACTIVSg200.csv -sopflow_solver
  # HIOP -sopflow_Ns -1 -sopflow_subproblem_solver HIOP
  # -sopflow_subproblem_model POWER_BALANCE_HIOP -hiop_compute_mode CPU
  # -hiop_verbosity_level 3 -opflow_initialization ACPF
  # -opflow_ignore_lineflow_constraints 1 DEPENDS IPOPT HIOP ) endif()

  # exago_add_test( NAME SOPFLOW_9bus_HIOP_SERIAL_GPU COMMAND ${RUNCMD}
  # $<TARGET_FILE:app_sopflow> -netfile
  # ${CMAKE_SOURCE_DIR}/datafiles/case9/case9mod_gen3_wind.m -scenfile
  # ${CMAKE_SOURCE_DIR}/datafiles/case9/10_scenarios_9bus.csv -sopflow_solver
  # HIOP -sopflow_Ns -1 -opflow_initialization ACPF -sopflow_subproblem_solver
  # HIOP -sopflow_subproblem_model PBPOLRAJAHIOP -hiop_compute_mode GPU
  # -hiop_verbosity_level 3 DEPENDS IPOPT HIOP )

  # if (FALSE) exago_add_test( NAME SOPFLOW_200bus_HIOP_SERIAL_GPU COMMAND
  # ${RUNCMD} $<TARGET_FILE:app_sopflow> -netfile
  # ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.m -scenfile
  # ${CMAKE_SOURCE_DIR}/datafiles/10_scenarios_ACTIVSg200.csv -sopflow_solver
  # HIOP -sopflow_Ns -1 -opflow_initialization ACPF
  # -opflow_ignore_lineflow_constraints 1 -sopflow_subproblem_solver HIOP
  # -sopflow_subproblem_model PBPOLRAJAHIOP -hiop_compute_mode CPU
  # -hiop_verbosity_level 3 DEPENDS IPOPT HIOP )

  # exago_add_test( NAME SOPFLOW_200bus_HIOP_MPI_GPU COMMAND ${MPICMD} "-n" 3
  # $<TARGET_FILE:app_sopflow> -netfile
  # ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.m -scenfile
  # ${CMAKE_SOURCE_DIR}/datafiles/10_scenarios_ACTIVSg200.csv -sopflow_solver
  # HIOP -sopflow_Ns -1 -opflow_initialization ACPF
  # -opflow_ignore_lineflow_constraints 1 -sopflow_subproblem_solver HIOP
  # -sopflow_subproblem_model PBPOLRAJAHIOP -hiop_compute_mode GPU
  # -hiop_verbosity_level 3 DEPENDS IPOPT HIOP ) endif()
endif()
