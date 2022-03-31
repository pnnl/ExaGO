if(EXAGO_ENABLE_RAJA AND EXAGO_ENABLE_CUDA)
  set_source_files_properties(objective_acopf.cpp PROPERTIES LANGUAGE CUDA)
endif()

add_executable(
  objective_acopf objective_acopf.cpp ../../utils/test_acopf_utils.cpp
)
target_link_libraries(objective_acopf ExaGO::OPFLOW)
target_include_directories(objective_acopf PRIVATE ../../utils)
target_include_directories(objective_acopf PRIVATE ../)

# Network files to run on - doing 3 and 600 bus examples.
set(prefix ${EXAGO_DATAFILES_DIR}/unit/opflow/objective/)
set(obj_network_files OF_unittestx3.m OF_unittestx600.m)
# Map num_copies to a specific netfile
set(num_copies 3 600)

list(TRANSFORM obj_network_files PREPEND ${prefix})

set(solvers IPOPT HIOP HIOP HIOPSPARSE)
set(models POWER_BALANCE_POLAR POWER_BALANCE_HIOP PBPOLRAJAHIOP
           POWER_BALANCE_POLAR
)
set(dependencies IPOPT HIOP HIOP HIOP_SPARSE)

if(EXAGO_INSTALL_TESTS)
  install(TARGETS objective_acopf DESTINATION tests/unit/opflow/objective)
  foreach(
    model
    solver
    dependency
    IN
    ZIP_LISTS
    models
    solvers
    dependencies
  )
    # Iterate over networks, matching network file to num_copies
    foreach(network num IN ZIP_LISTS obj_network_files num_copies)
      get_filename_component(net ${network} NAME)
      set(testname "UNIT_TESTS_OBJECTIVE_SIZE_${net}_${solver}_${model}")
      exago_add_test(
        NAME
        ${testname}
        DEPENDS
        ${dependency}
        COMMAND
        ${RUNCMD}
        $<TARGET_FILE:objective_acopf>
        -opflow_model
        ${model}
        -opflow_solver
        ${solver}
        -netfile
        ${network}
        -num_copies
        ${num}
        -validation
        ${CMAKE_CURRENT_SOURCE_DIR}/of.csv
      )
    endforeach()
  endforeach()
endif()
