#[[ Functionality tests for ExaGO PFLOW applications ]]

set_source_files_properties(pflow.c selfcheck.cpp PROPERTIES LANGUAGE CXX)

add_executable(test_pflow_functionality selfcheck.cpp)

target_include_directories(
  test_pflow_functionality PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
                                  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

#[[ Configure PFLOW functionality tests ]]
target_link_libraries(test_pflow_functionality ExaGO::PFLOW)
if(EXAGO_INSTALL_TESTS)
  message(STATUS "Configuring PFLOW Functionality Tests.")

  # PFLOW Functionality Tests
  if(EXAGO_ENABLE_MPI)
    set(np_max, 3)
  else()
    set(np_max, 1)
  endif()
  foreach(np RANGE 1 np_max)
    exago_add_test(
      NAME
      FUNCTIONALITY_TEST_PFLOW_TESTSUITE_${np}_proc
      COMMAND
      ${MPICMD}
      "-n"
      ${np}
      ${EXAGO_EXTRA_MPI_FLAGS}
      $<TARGET_FILE:test_pflow_functionality>
      ${CMAKE_CURRENT_SOURCE_DIR}/pflow${np}procs.toml
    )
  endforeach(np)

  install(TARGETS test_pflow_functionality
          RUNTIME DESTINATION tests/functionality
  )

endif()
