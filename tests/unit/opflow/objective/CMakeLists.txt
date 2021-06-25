if(EXAGO_ENABLE_RAJA AND EXAGO_ENABLE_GPU)
  set_source_files_properties(objectiveAcopf.cpp PROPERTIES LANGUAGE CUDA)
endif()

add_executable(objectiveAcopf objectiveAcopf.cpp ../../utils/TestAcopfUtils.cpp)
target_link_libraries(objectiveAcopf ExaGO::OPFLOW)
target_include_directories(objectiveAcopf PRIVATE ../../utils)

if(EXAGO_INSTALL_TESTS)
  install(TARGETS objectiveAcopf DESTINATION ./)
  # TODO:
  # - Test against the different solvers
  # - Use the 3 and 600 bus networks, not defaults
  exago_add_test(
    NAME
    "UNIT_TESTS_OBJECTIVE"
    DEPENDS
    HIOP
    COMMAND
    ${RUNCMD}
    $<TARGET_FILE:objectiveAcopf>
    -opflow_genbusvoltage
    VARIABLE_WITHIN_BOUNDS
    NETFILES
    ${network_files}
  )

endif()

