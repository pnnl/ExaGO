if(EXAGO_ENABLE_RAJA AND EXAGO_ENABLE_CUDA)
  set_source_files_properties(gradient_acopf.cpp PROPERTIES LANGUAGE CUDA)
endif()

add_executable(
  gradient_acopf gradient_acopf.cpp ../../utils/test_acopf_utils.cpp
)
target_link_libraries(gradient_acopf ExaGO::OPFLOW)
target_include_directories(gradient_acopf PRIVATE ../../utils)
target_include_directories(gradient_acopf PRIVATE ../)

# Network files to run on
set(prefix ${EXAGO_DATAFILES_DIR}/unit/opflow/gradient/)
set(obj_network_files OFG_unittest1.m OFG_unittestx3.m OFG_unittestx600.m)

list(TRANSFORM obj_network_files PREPEND ${prefix})

set(models POWER_BALANCE_POLAR POWER_BALANCE_POLAR)
set(solvers IPOPT HIOPSPARSE)
set(dependencies IPOPT HIOP_SPARSE)

if(EXAGO_INSTALL_TESTS)
  install(TARGETS gradient_acopf DESTINATION tests/unit/opflow/gradient)
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
    # Iterate over networks
    foreach(network ${obj_network_files})
      get_filename_component(net ${network} NAME)
      set(testname "UNIT_TESTS_GRADIENT_${net}_${solver}_${model}")
      exago_add_test(
        NAME
        ${testname}
        DEPENDS
        ${dependency}
        COMMAND
        ${RUNCMD}
        $<TARGET_FILE:gradient_acopf>
        -opflow_model
        ${model}
        -opflow_solver
        ${solver}
        -netfile
        ${network}
      )
    endforeach()
  endforeach()
endif()
