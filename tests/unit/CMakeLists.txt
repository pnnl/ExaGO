if(EXAGO_ENABLE_RAJA AND EXAGO_ENABLE_GPU)
  set_source_files_properties(TestAcopf.cpp PROPERTIES LANGUAGE CUDA)
endif()

add_executable(testAcopf TestAcopf.cpp utils/TestAcopfUtils.cpp)
target_link_libraries(testAcopf ExaGO::OPFLOW)
target_include_directories(testAcopf PRIVATE ./utils)

add_executable(testErrorHandler TestErrorHandler.cpp)
target_link_libraries(testErrorHandler ExaGO::UTILS)

if(EXAGO_INSTALL_TESTS)
  install(TARGETS testAcopf testErrorHandler DESTINATION ./)

  # Unit tests
  exago_add_test(
    NAME
    "UNIT_TESTS_OPFLOW"
    DEPENDS
    HIOP
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:testAcopf>
    -opflow_genbusvoltage
    VARIABLE_WITHIN_BOUNDS
    -validation_dir
    ${CMAKE_SOURCE_DIR}/datafiles/test_validation
    NETFILES
    ${network_files}
  )

  exago_add_test(NAME UNIT_TESTS_UTILS COMMAND $<TARGET_FILE:testErrorHandler>)
endif()
