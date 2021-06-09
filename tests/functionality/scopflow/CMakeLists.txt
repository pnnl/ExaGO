#[[ Functionality tests for ExaGO SCOPFLOW applications ]]

set_source_files_properties(
  scopflow.c scopflowselfcheck.c PROPERTIES LANGUAGE CXX
)

add_executable(testScopflowFunctionality scopflow.c scopflowselfcheck.c)
target_include_directories(
  testScopflowFunctionality PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

set(tolerance "1e-3")

#[[ Configure SCOPFLOW functionality tests ]]
target_link_libraries(testScopflowFunctionality ExaGO::SCOPFLOW)
if(EXAGO_INSTALL_TESTS)
  install(TARGETS testScopflowFunctionality
          RUNTIME DESTINATION tests/functionality
  )

  foreach(nc RANGE 1 9)
    # SCOPFLOW single-period tests
    exago_add_test(
      NAME
      SCOPFLOW_9bus_${nc}cont
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
      ${nc}
      DEPENDS
      IPOPT
    )
    # SCOPFLOW multi-period tests
    exago_add_test(
      NAME
      SCOPFLOW_9bus_multiperiod_${nc}cont
      COMMAND
      ${RUNCMD}
      $<TARGET_FILE:app_scopflow>
      -netfile
      ${CMAKE_SOURCE_DIR}/datafiles/case9/case9mod_gen3_wind.m
      -ctgcfile
      ${CMAKE_SOURCE_DIR}/datafiles/case9/case9.cont
      -scopflow_enable_multiperiod
      1
      -scopflow_ploadprofile
      ${CMAKE_SOURCE_DIR}/datafiles/case9/load_P.csv
      -scopflow_qloadprofile
      ${CMAKE_SOURCE_DIR}/datafiles/case9/load_Q.csv
      -scopflow_windgenprofile
      ${CMAKE_SOURCE_DIR}/datafiles/case9/scenarios_9bus.csv
      -scopflow_solver
      IPOPT
      -scopflow_Nc
      ${nc}
      DEPENDS
      IPOPT
    )
  endforeach(nc)

  # Tests for HIOP pridecomp
  if(EXAGO_ENABLE_HIOP_DISTRIBUTED)
    exago_add_test(
      NAME
      SCOPFLOW_9bus_HIOP_SERIAL
      COMMAND
      ${RUNCMD}
      $<TARGET_FILE:app_scopflow>
      -netfile
      ${CMAKE_SOURCE_DIR}/datafiles/case9/case9mod.m
      -ctgcfile
      ${CMAKE_SOURCE_DIR}/datafiles/case9/case9.cont
      -scopflow_solver
      HIOP
      -scopflow_Nc
      -1
      -opflow_initialization
      ACPF
      DEPENDS
      IPOPT
      HIOP
    )

    exago_add_test(
      NAME
      SCOPFLOW_200bus_HIOP_SERIAL
      COMMAND
      ${RUNCMD}
      $<TARGET_FILE:app_scopflow>
      -netfile
      ${CMAKE_SOURCE_DIR}/datafiles/case9/case9mod.m
      -ctgcfile
      ${CMAKE_SOURCE_DIR}/datafiles/case9/case9.cont
      -scopflow_solver
      HIOP
      -scopflow_Nc
      20
      -opflow_initialization
      ACPF
      DEPENDS
      IPOPT
      HIOP
    )

    exago_add_test(
      NAME
      SCOPFLOW_200bus_HIOP_MPI
      COMMAND
      ${MPICMD}
      "-n"
      4
      $<TARGET_FILE:app_scopflow>
      -netfile
      ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.m
      -ctgcfile
      ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.cont
      -scopflow_solver
      HIOP
      -scopflow_Nc
      -1
      -opflow_initialization
      ACPF
      DEPENDS
      IPOPT
      HIOP
    )
  endif()

  # SCOPFLOW Functionality Tests
  message(STATUS "Configuring SCOPFLOW Functionality Tests.")
  file(GLOB network_file ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.m)
  get_filename_component(input_file ${network_file} NAME)

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SCOPFLOW_200bus_IPOPT
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:testScopflowFunctionality>
    -netfile
    ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.m
    -ctgcfile
    ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.cont
    -scopflow_solver
    IPOPT
    -scopflow_Nc
    0
    -scopflow_enable_multiperiod
    0
    -scopflow_tolerance
    ${tolerance}
    -opflow_initialization
    ACPF
    -opflow_genbusvoltage
    VARIABLE_WITHIN_BOUNDS
    DEPENDS
    IPOPT
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SCOPFLOW_200bus_IPOPT_15nCont
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:testScopflowFunctionality>
    -netfile
    ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.m
    -ctgcfile
    ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.cont
    -scopflow_solver
    IPOPT
    -scopflow_Nc
    15
    -scopflow_enable_multiperiod
    0
    -scopflow_tolerance
    ${tolerance}
    -opflow_initialization
    ACPF
    -opflow_genbusvoltage
    VARIABLE_WITHIN_BOUNDS
    DEPENDS
    IPOPT
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SCOPFLOW_200bus_IPOPT_mode_Corrective
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:testScopflowFunctionality>
    -netfile
    ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.m
    -ctgcfile
    ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.cont
    -scopflow_solver
    IPOPT
    -scopflow_Nc
    15
    -scopflow_enable_multiperiod
    0
    -scopflow_mode
    1
    -scopflow_tolerance
    ${tolerance}
    -opflow_initialization
    ACPF
    -opflow_genbusvoltage
    VARIABLE_WITHIN_BOUNDS
    DEPENDS
    IPOPT
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SCOPFLOW_9bus_IPOPT_9contingencies
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:testScopflowFunctionality>
    -netfile
    ${CMAKE_SOURCE_DIR}/datafiles/case9/case9mod_gen3_wind.m
    -ctgcfile
    ${CMAKE_SOURCE_DIR}/datafiles/case9/case9.cont
    -scopflow_solver
    IPOPT
    -scopflow_Nc
    9
    -scopflow_enable_multiperiod
    0
    -scopflow_tolerance
    ${tolerance}
    -opflow_initialization
    ACPF
    -opflow_genbusvoltage
    VARIABLE_WITHIN_BOUNDS
    DEPENDS
    IPOPT
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SCOPFLOW_9bus_IPOPT_9contingencies_load_pq_scene_wind
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:testScopflowFunctionality>
    -netfile
    ${CMAKE_SOURCE_DIR}/datafiles/case9/case9mod_gen3_wind.m
    -ctgcfile
    ${CMAKE_SOURCE_DIR}/datafiles/case9/case9.cont
    -scopflow_ploadprofile
    ${CMAKE_SOURCE_DIR}/datafiles/case9/load_P.csv
    -scopflow_qloadprofile
    ${CMAKE_SOURCE_DIR}/datafiles/case9/load_Q.csv
    -scopflow_windgenprofile
    ${CMAKE_SOURCE_DIR}/datafiles/case9/scenarios_9bus.csv
    -scopflow_solver
    IPOPT
    -scopflow_Nc
    9
    -scopflow_mode
    1
    -scopflow_enable_multiperiod
    1
    -scopflow_dT
    5.0
    -scopflow_duration
    0.16666667
    -scopflow_tolerance
    ${tolerance}
    -opflow_initialization
    ACPF
    -opflow_genbusvoltage
    VARIABLE_WITHIN_BOUNDS
    DEPENDS
    IPOPT
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SCOPFLOW_200bus_IPOPT_10contingencies_load_pq_scene_wind
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:testScopflowFunctionality>
    -netfile
    ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.m
    -ctgcfile
    ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.cont
    -scopflow_ploadprofile
    ${CMAKE_SOURCE_DIR}/datafiles/TAMU200_scenarios/load_P.csv
    -scopflow_qloadprofile
    ${CMAKE_SOURCE_DIR}/datafiles/TAMU200_scenarios/load_Q.csv
    -scopflow_windgenprofile
    ${CMAKE_SOURCE_DIR}/datafiles/TAMU200_scenarios/scenarios_200bus.csv
    -scopflow_solver
    IPOPT
    -scopflow_Nc
    10
    -scopflow_enable_multiperiod
    1
    -scopflow_dT
    5.0
    -scopflow_duration
    0.16666667
    -scopflow_tolerance
    ${tolerance}
    -opflow_initialization
    ACPF
    -opflow_genbusvoltage
    VARIABLE_WITHIN_BOUNDS
    DEPENDS
    IPOPT
  )
endif()
