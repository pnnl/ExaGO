set(pythonpath ${CMAKE_INSTALL_PREFIX}/${EXAGO_PYTHON_SITELIB}:$ENV{PYTHONPATH})
set(pythonpath ${pythonpath}:${PROJECT_SOURCE_DIR}/tests/interfaces/python)

find_program(PYTEST_EXE "pytest")
if("${PYTEST_EXE}" STREQUAL "PYTEST_EXE-NOTFOUND")
  message(FATAL_ERROR "Pytest module path not found")
endif()

list(TRANSFORM RUNCMD APPEND " ")
string(CONCAT CMD ${RUNCMD} " " ${PYTEST_EXE})

# Set up a list of test to run based on current configuration
set(TEST_FILES "test_0_initialize.py")

# The PFLOW API should always be available, so let's test it
string(APPEND TEST_FILES " " "test_1_pflow.py")

# The OPFLOW API tests are currently written only for HiOP
if(EXAGO_ENABLE_HIOP)
  string(APPEND TEST_FILES " " "test_2_opflow.py")
endif()

# The SCOPFLOW and SOPFLOW API tests are currently written only for IPOPT
if(EXAGO_ENABLE_IPOPT)
  string(APPEND TEST_FILES " " "test_3_scopflow.py" " " "test_4_sopflow.py" " " "test_5_tcopflow.py")
endif()
string(APPEND TEST_FILES " " "test_9_finalize.py")

exago_add_test(
  NAME
  Python_pytest
  COMMAND
  bash
  -c
  "
  export PYTHONPATH=${pythonpath}
  cd ${PROJECT_SOURCE_DIR}/tests/interfaces/python
  ${CMD} -m nocomm -v ${TEST_FILES}
  "
)

exago_add_test(
  NAME
  Python_pytest_MPI
  COMMAND
  bash
  -c
  "
  export PYTHONPATH=${pythonpath}
  cd ${PROJECT_SOURCE_DIR}/tests/interfaces/python
  ${CMD} -m MPI -v ${TEST_FILES}
  "
  DEPENDS
  MPI
)
