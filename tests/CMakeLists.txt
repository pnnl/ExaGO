if(EXAGO_RUN_TESTS)

  # Test cases
  set(network_files
      ${CMAKE_SOURCE_DIR}/datafiles/case9/case9mod.m
      ${CMAKE_SOURCE_DIR}/datafiles/case118.m
      ${CMAKE_SOURCE_DIR}/datafiles/case_ACTIVSg200.m
  )
  set(network_files_no200 ${CMAKE_SOURCE_DIR}/datafiles/case9/case9mod.m
                          ${CMAKE_SOURCE_DIR}/datafiles/case118.m
  )

  message(STATUS "Configuring tests")

  # Unit tests
endif(EXAGO_RUN_TESTS)

add_subdirectory(unit)
add_subdirectory(functionality)

if(EXAGO_ENABLE_PYTHON)
  add_subdirectory(interfaces/python)
endif()
