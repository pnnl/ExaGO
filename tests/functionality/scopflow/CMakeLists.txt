#[[ Functionality tests for ExaGO SCOPFLOW applications ]]

add_executable(testScopflowFunctionality scopflow.c scopflowselfcheck.c)
target_include_directories(testScopflowFunctionality PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

#[[ Configure SCOPFLOW functionality tests ]]
target_link_libraries(testScopflowFunctionality ${EXAGO_APP_LIBS})
if(EXAGO_INSTALL_TESTS) 
  install(TARGETS testScopflowFunctionality RUNTIME DESTINATION tests/functionality)
endif()
