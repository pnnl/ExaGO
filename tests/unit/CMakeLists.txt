add_subdirectory(opflow)

if(EXAGO_ENABLE_RAJA AND EXAGO_ENABLE_CUDA)
  set_source_files_properties(test_acopf.cpp PROPERTIES LANGUAGE CUDA)
endif()

add_executable(test_acopf test_acopf.cpp utils/test_acopf_utils.cpp)
target_link_libraries(test_acopf ExaGO::OPFLOW)
target_include_directories(test_acopf PRIVATE ./utils)

add_executable(test_error_handler test_error_handler.cpp)
target_link_libraries(test_error_handler ExaGO::UTILS)

add_executable(test_logger test_logger.cpp)
target_link_libraries(test_logger ExaGO::UTILS)

if(EXAGO_INSTALL_TESTS)
  install(TARGETS test_logger test_acopf test_error_handler
          DESTINATION tests/unit
  )

  # Unit tests
  exago_add_test(
    NAME
    "UNIT_TESTS_OPFLOW"
    DEPENDS
    HIOP
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:test_acopf>
    -opflow_genbusvoltage
    VARIABLE_WITHIN_BOUNDS
    -validation_dir
    ${CMAKE_SOURCE_DIR}/datafiles/test_validation
    NETFILES
    ${network_files}
  )

  exago_add_test(
    NAME UNIT_TESTS_UTILS COMMAND $<TARGET_FILE:test_error_handler>
  )

  string(REPLACE ";" " " MPICMD_STR "${MPICMD}")
  exago_add_test(
    NAME
    UNIT_TESTS_LOGGER
    COMMAND
    ${CMAKE_CURRENT_SOURCE_DIR}/run_test_logger.pl
    ${MPICMD_STR}
    $<TARGET_FILE:test_logger>
    DEPENDS
    EXAGO_ENABLE_LOGGING
    )
endif()
