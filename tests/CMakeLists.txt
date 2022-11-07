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
endif(EXAGO_RUN_TESTS)

include_directories(${PROJECT_SOURCE_DIR}/tpl/toml11)
# TODO - ensure we can build and test with different subsets of ExaGO
if(EXAGO_ENABLE_IPOPT OR EXAGO_ENABLE_HIOP)
  add_subdirectory(unit)
endif()

add_subdirectory(functionality)

if(EXAGO_ENABLE_PYTHON)
  add_subdirectory(interfaces/python)
endif()
