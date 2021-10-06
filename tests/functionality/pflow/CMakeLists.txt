#[[ Functionality tests for ExaGO PFLOW applications ]]

set_source_files_properties(
  pflow.c pflowselfcheck.c selfcheck.cpp PROPERTIES LANGUAGE CXX
)

add_executable(testPflowFunctionality pflow.c pflowselfcheck.c)

target_include_directories(
  testPflowFunctionality PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(testPflowFunctionalityTOML selfcheck.cpp)
target_include_directories(
  testPflowFunctionalityTOML PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
                                    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

#[[ Configure PFLOW functionality tests ]]
target_link_libraries(testPflowFunctionality ExaGO::PFLOW)
target_link_libraries(testPflowFunctionalityTOML ExaGO::PFLOW)
if(EXAGO_INSTALL_TESTS)
  message(STATUS "Configuring PFLOW Functionality Tests.")

  # PFLOW Functionality Tests
  foreach(np RANGE 1 2)
    exago_add_test(
      NAME
      FUNCTIONALITY_TEST_PFLOW_TESTSUITE_${np}_proc
      COMMAND
      ${MPICMD}
      "-n"
      ${np}
      ${EXAGO_EXTRA_MPI_FLAGS}
      $<TARGET_FILE:testPflowFunctionalityTOML>
      ${CMAKE_CURRENT_SOURCE_DIR}/pflow.toml
      -no_optfile
    )
  endforeach(np)

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

  install(TARGETS testPflowFunctionality
          RUNTIME DESTINATION tests/functionality
  )
  install(TARGETS testPflowFunctionalityTOML
          RUNTIME DESTINATION tests/functionality
  )

endif()
