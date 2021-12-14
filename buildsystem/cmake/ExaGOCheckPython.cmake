set(EXAGO_MIN_PYTHON_VERSION "3.6")

find_package(Python COMPONENTS Interpreter Development)

if(${Python_VERSION} VERSION_LESS ${EXAGO_MIN_PYTHON_VERSION})
  message(
    STATUS
      "Found Python ${Python_VERSION}, but ExaGO requires "
      "${EXAGO_MIN_PYTHON_VERSION}. Option EXAGO_ENABLE_PYTHON will be disabled."
  )
  set(EXAGO_ENABLE_PYTHON OFF)
else()
  set(EXAGO_PYTHON_SITELIB
      "lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages"
  )
  message(STATUS "Found Python interpreter: ${Python_EXECUTABLE}")
  message(STATUS "Using Python site library prefix: ${EXAGO_PYTHON_SITELIB}")
  message(STATUS "ExaGO Python bindings will be installed in: "
                 "${CMAKE_INSTALL_PREFIX}/${EXAGO_PYTHON_SITELIB}"
  )
endif()
