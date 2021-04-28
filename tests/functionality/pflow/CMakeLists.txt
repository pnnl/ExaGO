#[[ Functionality tests for ExaGO PFLOW applications ]]

set_source_files_properties(pflow.c pflowselfcheck.c PROPERTIES LANGUAGE CXX)

add_executable(testPflowFunctionality pflow.c pflowselfcheck.c)
target_include_directories(
  testPflowFunctionality PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

#[[ Configure PFLOW functionality tests ]]
target_link_libraries(testPflowFunctionality ${EXAGO_APP_LIBS})
if(EXAGO_INSTALL_TESTS)
  install(TARGETS testPflowFunctionality
          RUNTIME DESTINATION tests/functionality
  )

  message(STATUS "Configuring PFLOW functionality tests")
  foreach(np RANGE 1 2)
    exago_add_test(
      NAME
      FUNCTIONALITY_TEST_PFLOW_${np}proc
      COMMAND
      ${MPICMD}
      "-n"
      ${np}
      ${EXAGO_EXTRA_MPI_FLAGS}
      $<TARGET_FILE:testPflowFunctionality>
      NETFILES
      ${CMAKE_SOURCE_DIR}/datafiles/case9/case9mod.m
    )
  endforeach(np)

endif()
