#[[

@file FindClangFormat.cmake
@brief Finds a clang-format executable and creates a custom target to check or
format ExaGO's source.

@note Adapted from this repository: https://github.com/zemasoft/clangformat-cmake

#]]

if(NOT CLANGFORMAT_EXECUTABLE)
  set(CLANGFORMAT_EXECUTABLE clang-format)
endif()

if(NOT EXISTS ${CLANGFORMAT_EXECUTABLE})
  find_program(clangformat_executable_tmp ${CLANGFORMAT_EXECUTABLE})
  if(clangformat_executable_tmp)
    set(CLANGFORMAT_EXECUTABLE ${clangformat_executable_tmp})
    unset(clangformat_executable_tmp)
  else()
    message(
      FATAL_ERROR
        "ExaGO Clang Format executable ${CLANGFORMAT_EXECUTABLE} not found!"
    )
  endif()
endif()

add_custom_target(
  lint
  COMMAND
    CMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR}
    CLANGFORMAT_EXECUTABLE=${CLANGFORMAT_EXECUTABLE}
    ${CMAKE_SOURCE_DIR}/buildsystem/lint/build.sh
  COMMENT "Running clang-format linting..."
)

add_custom_target(
  lint-apply
  COMMAND
    CMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR}
    CLANGFORMAT_EXECUTABLE=${CLANGFORMAT_EXECUTABLE}
    ${CMAKE_SOURCE_DIR}/buildsystem/lint-apply/build.sh
  COMMENT "Applying clang-format linting..."
)
