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

# Add toml11 subproject for selfcheck tests
include(ExternalProject)
ExternalProject_Add(
  toml11
  GIT_REPOSITORY https://github.com/ToruNiina/toml11.git
  GIT_TAG v3.7.0
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}
             -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
             -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
             -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
)
include_directories(${CMAKE_BINARY_DIR}/include)

add_subdirectory(unit)
add_subdirectory(functionality)

if(EXAGO_ENABLE_PYTHON)
  add_subdirectory(interfaces/python)
endif()
