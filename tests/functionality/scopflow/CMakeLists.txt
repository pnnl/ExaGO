#[[ Functionality tests for ExaGO SCOPFLOW applications ]]

add_executable(testScopflowFunctionalityTOML selfcheck.cpp)

set(tolerance "1e-3")

#[[ Configure SCOPFLOW functionality tests ]]
target_link_libraries(testScopflowFunctionalityTOML ExaGO::SCOPFLOW)
if(EXAGO_INSTALL_TESTS)
  install(TARGETS testScopflowFunctionalityTOML
          RUNTIME DESTINATION tests/functionality
  )
  # SCOPFLOW Functionality Tests
  message(STATUS "Configuring SCOPFLOW Functionality Tests.")
  exago_add_test(
    NAME
    testScopflowFunctionalityTOML
    DEPENDS
    IPOPT
    COMMAND
    $<TARGET_FILE:testScopflowFunctionalityTOML>
    ${CMAKE_CURRENT_SOURCE_DIR}/testcases.toml
    -no_optfile
  )
  # --- Old smoke tests
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
endif()
