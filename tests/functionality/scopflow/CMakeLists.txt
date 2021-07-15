#[[ Functionality tests for ExaGO SCOPFLOW applications ]]

add_executable(testScopflowFunctionalityTOML selfcheck.cpp)
target_include_directories(
  testScopflowFunctionalityTOML PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..
)

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
    FUNCTIONALITY_TEST_SCOPFLOW_IPOPT_CONT_TESTSUITE
    DEPENDS
    IPOPT
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:testScopflowFunctionalityTOML>
    ${CMAKE_CURRENT_SOURCE_DIR}/ipopt_cont.toml
    -no_optfile
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SCOPFLOW_EMPAR_CONT_TESTSUITE
    DEPENDS
    IPOPT
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:testScopflowFunctionalityTOML>
    ${CMAKE_CURRENT_SOURCE_DIR}/empar_cont.toml
    -no_optfile
  )

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_SCOPFLOW_IPOPT_MULTIPERIOD_TESTSUITE
    DEPENDS
    IPOPT
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:testScopflowFunctionalityTOML>
    ${CMAKE_CURRENT_SOURCE_DIR}/ipopt_multiperiod.toml
    -no_optfile
  )

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
    $<TARGET_FILE:testScopflowFunctionalityTOML>
    ${CMAKE_CURRENT_SOURCE_DIR}/hiop_cont.toml
    -no_optfile
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
    4
    $<TARGET_FILE:testScopflowFunctionalityTOML>
    ${CMAKE_CURRENT_SOURCE_DIR}/hiop_nocont_mpi.toml
    -no_optfile
  )

  # Tests for HIOP pridecomp (these tests are duplicated above and should be
  # removed after verifying that they work
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
      -scopflow_subproblem_solver
      HIOP
      -scopflow_subproblem_model
      POWER_BALANCE_HIOP
      -hiop_compute_mode
      CPU
      -hiop_verbosity_level
      3
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
      ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.m
      -ctgcfile
      ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.cont
      -scopflow_solver
      HIOP
      -scopflow_Nc
      20
      -opflow_initialization
      ACPF
      -scopflow_subproblem_solver
      HIOP
      -scopflow_subproblem_model
      POWER_BALANCE_HIOP
      -hiop_compute_mode
      CPU
      -hiop_verbosity_level
      3
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
      -scopflow_subproblem_solver
      HIOP
      -scopflow_subproblem_model
      POWER_BALANCE_HIOP
      -hiop_compute_mode
      CPU
      -hiop_verbosity_level
      3
      -opflow_initialization
      ACPF
      DEPENDS
      IPOPT
      HIOP
    )

    # These GPU tests will fail because of the bad_cast bug in HIOP
    exago_add_test(
      NAME
      SCOPFLOW_9bus_HIOP_SERIAL_GPU
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
      -scopflow_subproblem_solver
      HIOP
      -scopflow_subproblem_model
      PBPOLRAJAHIOP
      -hiop_compute_mode
      GPU
      -hiop_verbosity_level
      3
      DEPENDS
      IPOPT
      HIOP
    )

    exago_add_test(
      NAME
      SCOPFLOW_200bus_HIOP_SERIAL_GPU
      COMMAND
      ${RUNCMD}
      $<TARGET_FILE:app_scopflow>
      -netfile
      ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.m
      -ctgcfile
      ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.cont
      -scopflow_solver
      HIOP
      -scopflow_Nc
      20
      -opflow_initialization
      ACPF
      -scopflow_subproblem_solver
      HIOP
      -scopflow_subproblem_model
      PBPOLRAJAHIOP
      -hiop_compute_mode
      CPU
      -hiop_verbosity_level
      3
      DEPENDS
      IPOPT
      HIOP
    )

    exago_add_test(
      NAME
      SCOPFLOW_200bus_HIOP_MPI_GPU
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
      -scopflow_subproblem_solver 
      HIOP
      -scopflow_subproblem_model
      PBPOLRAJAHIOP
      -hiop_compute_mode
      GPU
      -hiop_verbosity_level
      3
      DEPENDS
      IPOPT
      HIOP
    )
  endif()
endif()
