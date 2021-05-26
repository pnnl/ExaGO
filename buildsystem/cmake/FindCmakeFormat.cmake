#[[

@file FindClangFormat.cmake
@brief Finds a cmake-format executable and creates a custom target to check or
format ExaGO's source.

@note Adapted from this repository: https://github.com/zemasoft/cmakeformat-cmake

#]]

if(NOT CMAKEFORMAT_EXECUTABLE)
  set(CMAKEFORMAT_EXECUTABLE cmake-format)
endif()

if(NOT EXISTS ${CMAKEFORMAT_EXECUTABLE})
  find_program(cmakeformat_executable_tmp ${CMAKEFORMAT_EXECUTABLE})
  if(cmakeformat_executable_tmp)
    set(CMAKEFORMAT_EXECUTABLE ${cmakeformat_executable_tmp})
    unset(cmakeformat_executable_tmp)
  else()
    message(
      FATAL_ERROR
        "ExaGO Clang Format executable ${CMAKEFORMAT_EXECUTABLE} not found!"
    )
  endif()
endif()

add_custom_target(
  cmake-lint
  COMMAND
    CMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR}
    CMAKEFORMAT_EXECUTABLE=${CMAKEFORMAT_EXECUTABLE}
    ${CMAKE_SOURCE_DIR}/buildsystem/cmake-lint/build.sh --run
  COMMENT "Running cmake-format linting..."
)

add_custom_target(
  cmake-lint-apply
  COMMAND
    CMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR}
    CMAKEFORMAT_EXECUTABLE=${CMAKEFORMAT_EXECUTABLE}
    ${CMAKE_SOURCE_DIR}/buildsystem/cmake-lint-apply/build.sh --run
  COMMENT "Applying cmake-format linting..."
)
