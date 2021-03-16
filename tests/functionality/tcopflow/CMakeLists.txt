#[[ Functionality tests for ExaGO TCOPFLOW applications ]]

add_executable(testTcopflowFunctionality tcopflow.c tcopflowselfcheck.c)
target_include_directories(testTcopflowFunctionality PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

#[[ Configure TCOPFLOW functionality tests ]]
target_link_libraries(testTcopflowFunctionality ${EXAGO_APP_LIBS})
if(EXAGO_INSTALL_TESTS) 
  install(TARGETS testTcopflowFunctionality RUNTIME DESTINATION bin/tests/functionality)
endif()
