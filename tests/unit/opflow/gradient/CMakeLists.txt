if(EXAGO_ENABLE_RAJA AND EXAGO_ENABLE_GPU)
  set_source_files_properties(gradient_acopf.cpp PROPERTIES LANGUAGE CUDA)
endif()

set(GRADIENT_SRC
  gradient_acopf.cpp
  ../../utils/test_acopf_utils.cpp
  ../../../../../src/opflow/model/power_bal_polar/pbpol.cpp)

add_executable(gradient_acopf ${GRADIENT_SRC})

target_link_libraries(gradient_acopf ExaGO::OPFLOW)
target_include_directories(gradient_acopf PRIVATE ../../utils)
target_include_directories(gradient_acopf PRIVATE ../)

# Network files to run on - doing 3 and 600 bus examples.
set(prefix ${EXAGO_DATAFILES_DIR}/unit/opflow/gradient/)
set(obj_network_files OFG_unittestx3.m OFG_unittestx600.m)
# Map num_copies to a specific netfile
set(num_copies 3 600)

list(TRANSFORM obj_network_files PREPEND ${prefix})

set(solvers IPOPT HIOP HIOP HIOPSPARSE)
set(models POWER_BALANCE_POLAR POWER_BALANCE_HIOP PBPOLRAJAHIOP POWER_BALANCE_POLAR)
# Logic in CMake is HiOp::Sparse/HIOP_SPARSE, while opflow options are HIOPSPARSE
set(dependencies IPOPT HIOP HIOP HIOP_SPARSE)

if(EXAGO_INSTALL_TESTS)
  install(TARGETS gradient_acopf DESTINATION tests/unit/opflow/gradient)
  # Iterate over networks, matching network file to num_copies
  foreach(network num IN ZIP_LISTS obj_network_files num_copies)
    get_filename_component(net ${network} NAME)
    foreach(model solver dependency IN ZIP_LISTS models solvers dependencies)
      set(testname "UNIT_TESTS_GRADIENT_SIZE_${num}_${solver}_${model}")
      exago_add_test(
        NAME
        ${testname}
        DEPENDS
        ${dependency}
        COMMAND
        ${RUNCMD}
        $<TARGET_FILE:gradient_acopf>
        -opflow_include_loadloss_variables 1
        -opflow_include_powerimbalance_variables 1
        -netfile
        ${network}
        -opflow_model
        ${model}
        -opflow_solver
        ${solver}
        -num_copies
        ${num}
      )
    endforeach()
  endforeach()
endif()

