#[[ Functionality tests for ExaGO OPFLOW applications ]]

set_source_files_properties(opflow.c opflowselfcheck.c PROPERTIES LANGUAGE CXX)

add_executable(testOpflowFunctionality opflow.c opflowselfcheck.c)
target_include_directories(
  testOpflowFunctionality PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

set(tolerance "1e-3")

#[[ Configure OPFLOW functionality tests ]]
target_link_libraries(testOpflowFunctionality ${EXAGO_APP_LIBS})
if(EXAGO_INSTALL_TESTS)
  install(TARGETS testOpflowFunctionality
          RUNTIME DESTINATION tests/functionality
  )

  message(STATUS "Running OPFLOW model tests")
  foreach(opflow_init FROMFILE FLATSTART ACPF MIDPOINT)
    foreach(opflow_obj MIN_GEN_COST MIN_GENSETPOINT_DEVIATION)
      exago_add_test(
        NAME
        OPFLOW_POWER_BALANCE_POLAR_${opflow_init}_${opflow_obj}
        COMMAND
        ${RUNCMD}
        $<TARGET_FILE:app_opflow>
        -opflow_model
        POWER_BALANCE_POLAR
        -opflow_initialization
        ${opflow_init}
        -opflow_solver
        IPOPT
        -opflow_objective
        ${opflow_obj}
        -opflow_genbusvoltage
        VARIABLE_WITHIN_BOUNDS
        DEPENDS
        IPOPT
        NETFILES
        ${network_files}
      )
    endforeach()
  endforeach()

  foreach(opflow_init ACPF MIDPOINT)
    exago_add_test(
      NAME
      OPFLOW_HIOP_CPU_${opflow_init}
      COMMAND
      ${RUNCMD}
      $<TARGET_FILE:app_opflow>
      -opflow_initialization
      ${opflow_init}
      -hiop_compute_mode
      CPU
      -opflow_solver
      HIOP
      -opflow_ignore_lineflow_constraints
      0
      -opflow_model
      POWER_BALANCE_HIOP
      -opflow_genbusvoltage
      VARIABLE_WITHIN_BOUNDS
      DEPENDS
      HIOP
      NETFILES
      ${network_files}
    )

    if(EXAGO_ENABLE_HIOP_SPARSE)
      exago_add_test(
        NAME
        OPFLOW_HIOP_SPARSE_${opflow_init}
        COMMAND
        ${RUNCMD}
        $<TARGET_FILE:app_opflow>
        -opflow_initialization
        ${opflow_init}
        -hiop_compute_mode
        CPU
        -opflow_solver
        HIOPSPARSE
        -opflow_ignore_lineflow_constraints
        0
        -opflow_model
        POWER_BALANCE_POLAR
        -opflow_genbusvoltage
        VARIABLE_WITHIN_BOUNDS
        DEPENDS
        HIOP
        NETFILES
        ${network_files}
      )
    endif()

    if(EXAGO_ENABLE_GPU)
      set(hiop_compute_mode GPU)
    else()
      set(hiop_compute_mode HYBRID)
    endif()
    exago_add_test(
      NAME
      OPFLOW_HIOP_GPU_${opflow_init}
      COMMAND
      ${RUNCMD}
      $<TARGET_FILE:app_opflow>
      -opflow_initialization
      ${opflow_init}
      -hiop_compute_mode
      ${hiop_compute_mode}
      -opflow_solver
      HIOP
      -opflow_ignore_lineflow_constraints
      0
      -opflow_model
      PBPOLRAJAHIOP
      -opflow_genbusvoltage
      VARIABLE_WITHIN_BOUNDS
      DEPENDS
      GPU
      HIOP
      RAJA
      NETFILES
      ${network_files}
    )
  endforeach(opflow_init ACPF MIDPOINT)

  exago_add_test(
    NAME
    OPFLOW_TAO
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:app_opflow>
    -opflow_solver
    TAO
    -opflow_model
    POWER_BALANCE_CARTESIAN
    NETFILES
    ${network_files}
  )

  message(STATUS "Configuring OPFLOW functionality tests")
  message(STATUS "Skipping test `OPFLOW_HIOP_POWER_BALANCE_HIOP_case_ACTIVSg200"
                 ".m_FUNCTIONALITY_TEST`."
  )
  message(STATUS "See issue https://gitlab.pnnl.gov/exasgd/frameworks/exago/-/i"
                 "ssues/157 for more information."
  )

  foreach(model POWER_BALANCE_HIOP PBPOLRAJAHIOP)
    exago_add_test(
      NAME
      FUNCTIONALITY_TEST_OPFLOW_HIOP_${model}
      DEPENDS
      HIOP
      RAJA
      NETFILES
      ${network_files_no200}
      COMMAND
      ${RUNCMD}
      $<TARGET_FILE:testOpflowFunctionality>
      -opflow_solver
      HIOP
      -opflow_model
      ${model}
      -opflow_tolerance
      ${tolerance}
      -opflow_genbusvoltage
      VARIABLE_WITHIN_BOUNDS
    )
  endforeach()

  if(EXAGO_ENABLE_HIOP_SPARSE)
    exago_add_test(
      NAME
      FUNCTIONALITY_TEST_OPFLOW_HIOPSPARSE_${model}
      DEPENDS
      HIOP
      RAJA
      NETFILES
      ${network_files_no200}
      COMMAND
      ${RUNCMD}
      $<TARGET_FILE:testOpflowFunctionality>
      -opflow_solver
      HIOPSPARSE
      -opflow_model
      POWER_BALANCE_POLAR
      -opflow_tolerance
      ${tolerance}
      -opflow_genbusvoltage
      VARIABLE_WITHIN_BOUNDS
    )
  endif()

  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_OPFLOW_TAO_POWER_BALANCE_CARTESIAN
    NETFILES
    ${network_files}
    DEPENDS
    IPOPT
    HIOP
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:testOpflowFunctionality>
    -opflow_solver
    TAO
    -opflow_model
    POWER_BALANCE_CARTESIAN
    -opflow_tolerance
    ${tolerance}
  )

  foreach(model POWER_BALANCE_POLAR POWER_BALANCE_CARTESIAN
                CURRENT_BALANCE_CARTESIAN CURRENT_BALANCE_CARTESIAN2
  )
    exago_add_test(
      NAME
      FUNCTIONALITY_TEST_OPFLOW_IPOPT_${model}
      DEPENDS
      IPOPT
      NETFILES
      ${network_files}
      COMMAND
      ${RUNCMD}
      $<TARGET_FILE:testOpflowFunctionality>
      -opflow_solver
      IPOPT
      -opflow_model
      ${model}
      -netfile
      ${test_net}
      -opflow_tolerance
      ${tolerance}
    )
  endforeach()

  # Functionality test using IPOPT Adapter for 2nd order derivative test
  exago_add_test(
    NAME
    FUNCTIONALITY_TEST_OPFLOW_IPOPT_DEBUG
    DEPENDS
    IPOPT
    HIOP
    RAJA
    NETFILES
    ${CMAKE_SOURCE_DIR}/datafiles/case9/case9mod.m
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:app_opflow>
    -opflow_solver
    HIOP
    -opflow_model
    POWER_BALANCE_HIOP
    -hiop_compute_mode
    CPU
    -hiop_ipopt_debug
    -opflow_genbusvoltage
    VARIABLE_WITHIN_BOUNDS
  )

endif()
