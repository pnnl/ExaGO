find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  # Update submodules as needed
  message(STATUS "Submodule update")
  execute_process(
    COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    RESULT_VARIABLE GIT_SUBMOD_RESULT
  )
  if(NOT GIT_SUBMOD_RESULT EQUAL "0")
    message(
      FATAL_ERROR
        "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules"
    )
  endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/tpl/toml11/CMakeLists.txt")
  message(
    FATAL_ERROR
      "The toml submodule was not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again."
  )
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/tpl/spdlog/CMakeLists.txt")
  message(
    FATAL_ERROR
      "The spdlog submodule was not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again."
  )
endif()
