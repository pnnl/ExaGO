set(pythonpath ${CMAKE_INSTALL_PREFIX}/${EXAGO_PYTHON_SITELIB}:$ENV{PYTHONPATH})

# Want a space separated string of the RUNCMD, so we append a space to each
# element for concatentation later
list(TRANSFORM RUNCMD APPEND " ")

foreach(testcase test_pflow test_opflow)
  string(CONCAT CMD ${RUNCMD} ${Python_EXECUTABLE} " " ${testcase}.py)
  exago_add_test(
    NAME
    Python_${testcase}
    COMMAND
    bash
    -c
    "
    export PYTHONPATH=${pythonpath}
    cd ${PROJECT_SOURCE_DIR}/tests/interfaces/python
    ${CMD}
    "
  )
endforeach()
