set(EXAGO_PFLOW_APP_LIBS ${EXAGO_UTIL_LIBS} ${EXAGO_PS_LIBS}
                         ${EXAGO_PFLOW_LIBS} ${EXAGO_MATH_LIB} PETSC::SNES
)

if(EXAGO_ENABLE_MPI)
  list(APPEND EXAGO_PFLOW_APP_LIBS MPI::MPI_CXX)
endif()

if(EXAGO_ENABLE_RAJA AND EXAGO_ENABLE_GPU)
  set_source_files_properties(TestAcopf.cpp PROPERTIES LANGUAGE CUDA)
endif()

add_executable(testAcopf TestAcopf.cpp utils/TestAcopfUtils.cpp)
target_link_libraries(testAcopf ${EXAGO_OPFLOW_LIBS} ${EXAGO_PFLOW_APP_LIBS})
target_include_directories(testAcopf PRIVATE ./utils)

if(EXAGO_INSTALL_TESTS)
  install(TARGETS testAcopf DESTINATION tests/unit)
endif()
