if(EXAGO_ENABLE_RAJA AND EXAGO_ENABLE_GPU)
  set_source_files_properties(objectiveAcopf.cpp PROPERTIES LANGUAGE CUDA)
endif()

add_executable(objectiveAcopf objectiveAcopf.cpp ../../utils/TestAcopfUtils.cpp)
target_link_libraries(objectiveAcopf ExaGO::OPFLOW)
target_include_directories(objectiveAcopf PRIVATE ../../utils)

# Network files to run on - doing 3 and 600 bus examples.
set(prefix ${PROJECT_SOURCE_DIR}/tests/unit/opflow/objective/)
set(obj_network_files 
  OF-unittestx3.m 
  OF-unittestx600.m)
# Map num_copies to a specific netfile
set(num_copies
  3
  600)

list(TRANSFORM obj_network_files PREPEND ${prefix})
install(FILES ${obj_network_files} DESTINATION ${EXAGO_DATAFILES_DIR})

if(EXAGO_INSTALL_TESTS)
  install(TARGETS objectiveAcopf DESTINATION ./)
  # Iterate over networks, matching network file to num_copies
  foreach(network num IN ZIP_LISTS obj_network_files num_copies)
    get_filename_component(net ${network} NAME)
    set(testname "UNIT_TESTS_OBJECTIVE_SIZE_${net}")
    exago_add_test(
      NAME
      ${testname}
      DEPENDS
      HIOP
      COMMAND
      ${RUNCMD}
      $<TARGET_FILE:objectiveAcopf>
      -opflow_genbusvoltage
      VARIABLE_WITHIN_BOUNDS
      -opflow_initialization 
      FROMFILE
      -netfile
      ${network}
      -num_copies
      ${num}
    )
  endforeach()
endif()

