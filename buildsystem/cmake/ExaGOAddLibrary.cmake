#[[#############################################################################

ExaGOAddLibrary.cmake

@author Cameron Rutherford <cameron.rutherford@pnnl.gov>
@remark This file exports the following: exago_add_library. The macro
  described here should allow all libraries in any CMakeLists.txt to be
  written identically.
@remark Using a marco instead of a function as none of the inputs to
  the macro should be modified in any capacity.

#]]
# ##############################################################################

#[[

@brief Creates an object library + shared/static library
@param STATIC_ONLY : Only create static lib (if EXAGO_BUILD_STATIC is enabled)
@param SHARED_ONLY : Only create shared lib (if EXAGO_BUILD_SHARED is enabled)
@param OUTPUT_NAME : Specifies the output name of the library to be created
@param SOURCES : List of source files used to create the library
@param HEADERS : List of header files to install with the library
@param INCLUDE_SUBDIR : optional subdirectory for headers to be installed
@param LINK_LIBRARIES : List of libraries that are linked against using
  target_link_libraries for each target within the macro

#]]

macro(exago_add_library target)
  set(options STATIC_ONLY SHARED_ONLY)
  set(oneValueArgs OUTPUT_NAME INCLUDE_SUBDIR)
  set(multiValueArgs SOURCES HEADERS LINK_LIBRARIES)

  # Parse arguments
  cmake_parse_arguments(
    exago_add_library "${options}" "${oneValueArgs}" "${multiValueArgs}"
    ${ARGN}
  )

  # Library types that we want to create
  set(_libtypes "")
  if(EXAGO_BUILD_STATIC AND (NOT exago_add_library_SHARED_ONLY))
    set(_libtypes "STATIC")
  endif()
  if(EXAGO_BUILD_SHARED AND (NOT exago_add_library_STATIC_ONLY))
    set(_libtypes "${_libtypes};SHARED")
  endif()

  # Build Libraries
  foreach(_libtype ${_libtypes})
    # Add library suffix so internal library names are unique
    if(${_libtype} MATCHES "STATIC")
      set(_lib_suffix "_static")
    else()
      set(_lib_suffix "_shared")
    endif()

    # Source files for target
    set(sources ${exago_add_library_SOURCES})

    # Obj target also needs a unique name
    set(obj_target ${target}_obj${_lib_suffix})

    # -- Create object library --

    add_library(${obj_target} OBJECT ${sources})

    if(exago_add_library_LINK_LIBRARIES)
      if(${_lib_type} MATCHES "STATIC")
        append_static_suffix(exago_add_library_LINK_LIBRARIES _all_libs)
      else()
        set(_all_libs ${exago_add_library_LINK_LIBRARIES})
      endif()
      target_link_libraries(${obj_target} ${_all_libs})
    endif()

    # Object libraries need PIC code enabled
    set_target_properties(
      ${obj_target} PROPERTIES POSITION_INDEPENDENT_CODE TRUE
    )

    set(_actual_target_name ${target}${_lib_suffix})

    add_library(
      ${_actual_target_name} ${_libtype} $<TARGET_OBJECTS:${obj_target}>
    )

    if(exago_add_library_LINK_LIBRARIES)
      target_link_libraries(
        ${_actual_target_name} ${exago_add_library_LINK_LIBRARIES}
      )
    endif()

    # Generic Include directories to be added Bulding : public, config/export
    # and shared/private headers Installing: installed include directory Can
    # also add macro option INCLUDE_DIRECTORIES to customize this
    target_include_directories(
      ${_actual_target_name}
      PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
             $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
             $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/private>
             $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )

    # Enabling namespace usage with exported targets
    add_library(ExaGO::${target} ALIAS ${_actual_target_name})

    if(exago_add_library_OUTPUT_NAME)
      set_target_properties(
        ${_actual_target_name}
        PROPERTIES OUTPUT_NAME ${exago_add_library_OUTPUT_NAME}
                   CLEAN_DIRECT_OUTPUT 1
      )
    else()
      set_target_properties(
        ${_actual_target_name} PROPERTIES OUTPUT_NAME ${target}
                                          CLEAN_DIRECT_OUTPUT 1
      )
    endif()

    install(
      TARGETS ${_actual_target_name}
      DESTINATION lib
      EXPORT exago-targets
    )

    # Install header files
    if(exago_add_library_HEADERS)
      install(FILES ${exago_add_library_HEADERS}
              DESTINATION "include/${exago_add_library_INCLUDE_SUBDIR}"
      )
    endif()

  endforeach()
endmacro()

# Macro to append static suffix to library names Currently using hard coded
# "_shared" but this can be changed
macro(append_static_suffix libs_in libs_out)
  set(_STATIC_LIB_SUFFIX "_static")
  set(${libs_out} "")
  foreach(_lib ${${libs_in}})
    if(TARGET ${_lib}${_STATIC_LIB_SUFFIX})
      list(APPEND ${libs_out} ${_lib}${_STATIC_LIB_SUFFIX})
    else()
      list(APPEND ${libs_out} ${_lib})
    endif()
  endforeach()
endmacro()
