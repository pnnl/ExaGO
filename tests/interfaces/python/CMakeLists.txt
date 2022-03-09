set(pythonpath ${CMAKE_INSTALL_PREFIX}/${EXAGO_PYTHON_SITELIB}:$ENV{PYTHONPATH})

# Want a space separated string of the RUNCMD, so we append a space to each
# element for concatentation later
list(TRANSFORM RUNCMD APPEND " ")

foreach(testcase test_pflow test_opflow)
  # Run MPI tests depending MPI being enabled
  string(CONCAT CMD ${RUNCMD} ${Python_EXECUTABLE} " "
                $<$<BOOL:${EXAGO_ENABLE_MPI}>:MPI_> ${testcase}.py
  )
  if(EXAGO_ENABLE_MPI)
    set(testname Python_${testcase}_MPI)
  else()
    set(testname Python_${testcase})
  endif()
  exago_add_test(
    NAME
    ${testname}
    COMMAND
    bash
    -c
    "
    export PYTHONPATH=${pythonpath}
    cd ${PROJECT_SOURCE_DIR}/tests/interfaces/python
    ${CMD}
    "
    DEPENDS
    HIOP
  )
endforeach()

# Test initializing without communicator
string(CONCAT CMD ${RUNCMD} ${Python_EXECUTABLE} " "
              $<$<BOOL:${EXAGO_ENABLE_MPI}>:MPI_> test_initialize.py
)
if(EXAGO_ENABLE_MPI)
  set(testname Python_test_initialize_MPI)
else()
  set(testname Python_test_initialize)
endif()
exago_add_test(
  NAME
  ${testname}
  COMMAND
  bash
  -c
  "
  export PYTHONPATH=${pythonpath}
  cd ${PROJECT_SOURCE_DIR}/tests/interfaces/python
  ${CMD}
  "
)
