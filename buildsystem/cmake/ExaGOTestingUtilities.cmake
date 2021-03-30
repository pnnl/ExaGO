#[[#############################################################################

ExaGOTestingUtilities.cmake

@author Asher Mancinelli <asher.mancinelli@pnnl.gov>
@remark List variable `network_files` must be defined in enclosing scope when
  this macro is invoked.
@remark This file exports the following: exago_add_test. Ideally, the functions
  described here should allow all tests in top-level CMakeLists to be written
  declaratively.

#]]

#[[

@brief Adds test to ExaGO's testing framework if all dependencies are truthy.
@param NAME: String for the test name to be created.
@param DEPENDS: Only add the test if the DEPENDS variables are truthy.
@param COMMAND: Command to be invoked for the test. For tests that define
  NETFILES (see below), the command will have `-netfile ${netfile}` appended to
  it for each netfile in the list.
@param NETFILES: List of netfiles this test will run with. If list is ommitted,
  the `-netfile` flag will not be appended to the test command.

#]]
function(exago_add_test)
  set(flags "")
  set(one_value_args NAME)
  set(multi_value_args DEPENDS COMMAND NETFILES)

  cmake_parse_arguments(
    "EXAGO" "${flags}" "${one_value_args}" "${multi_value_args}" ${ARGN}
  )

  set(has_all_dependencies ON)

  # Search the dependencies to ensure they are all found before adding the test.
  foreach(dependency ${EXAGO_DEPENDS})
    if(NOT ${EXAGO_ENABLE_${dependency}})
      set(has_all_dependencies OFF)
      break()
    endif()
  endforeach()

  # Create a unique name for each test based on the netfile it uses
  if(${has_all_dependencies})

    # If length of netfiles is 0, we're likely just using this function for
    # dependency checking. In this case, we just run with the command and forget
    # about netfiles.
    list(LENGTH EXAGO_NETFILES len)
    if(${len} EQUAL 0)
      add_test(NAME ${EXAGO_NAME} COMMAND ${EXAGO_COMMAND})
      return()
    endif()

    # Otherwise, add the -netfile command and generate len(EXAGO_NETFILES) tests
    foreach(netfile ${EXAGO_NETFILES})
      # Add -netfile <path> to test command
      set(command ${EXAGO_COMMAND})
      list(APPEND command -netfile ${netfile})

      # Create name unique to input name and datafile
      get_filename_component(netfile_name ${netfile} NAME)
      add_test(NAME "${EXAGO_NAME}_${netfile_name}" COMMAND ${command})
    endforeach()
  endif()

endfunction(exago_add_test)
