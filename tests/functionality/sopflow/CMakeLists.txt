#[[ Functionality tests for ExaGO SOPFLOW applications ]]

add_executable(testSopflowFunctionality sopflow.c sopflowselfcheck.c)
target_include_directories(testSopflowFunctionality PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

#[[ Configure SOPFLOW functionality tests ]]
target_link_libraries(testSopflowFunctionality ${EXAGO_APP_LIBS})
if(EXAGO_INSTALL_TESTS) 
  install(TARGETS testSopflowFunctionality RUNTIME DESTINATION bin/tests/functionality)
endif()
