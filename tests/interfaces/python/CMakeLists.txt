set(pythonpath ${CMAKE_INSTALL_PREFIX}/${EXAGO_PYTHON_SITELIB}:$ENV{PYTHONPATH})
set(pythonpath ${pythonpath}:${PROJECT_SOURCE_DIR}/tests/interfaces/python)

find_program(PYTEST_EXE "pytest")
if("${PYTEST_EXE}" STREQUAL "PYTEST_EXE-NOTFOUND")
  message(FATAL_ERROR "Pytest module path not found")
endif()

list(TRANSFORM RUNCMD APPEND " ")
string(CONCAT CMD ${RUNCMD} " " ${PYTEST_EXE})

if(EXAGO_ENABLE_HIOP)
  exago_add_test(
    NAME
    Python_pytest
    COMMAND
    bash
    -c
    "
    export PYTHONPATH=${pythonpath}
    cd ${PROJECT_SOURCE_DIR}/tests/interfaces/python
    ${CMD} -m nocomm -v .
    "
    DEPENDS
    HIOP
  )
  if(EXAGO_ENABLE_MPI)
    exago_add_test(
      NAME
      Python_pytest_MPI
      COMMAND
      bash
      -c
      "
      export PYTHONPATH=${pythonpath}
      cd ${PROJECT_SOURCE_DIR}/tests/interfaces/python
      ${CMD} -m MPI -v .
      "
      DEPENDS
      HIOP
    )
  endif()
endif()
