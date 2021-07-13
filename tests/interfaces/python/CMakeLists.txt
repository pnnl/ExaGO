set(pythonpath ${CMAKE_INSTALL_PREFIX}/${EXAGO_PYTHON_SITELIB}:$ENV{PYTHONPATH})

foreach(testcase test_config test_pflow test_opflow)
  exago_add_test(
    NAME
    Python_${testcase}
    COMMAND
    bash
    -c
    "
    export PYTHONPATH=${pythonpath}
    cd ${PROJECT_SOURCE_DIR}/tests/interfaces/python
    ${Python_EXECUTABLE} ${testcase}.py
    "
  )
endforeach()
