#[[ Functionality tests for ExaGO OPFLOW applications ]]

add_executable(testOpflowFunctionality opflow.c opflowselfcheck.c)
target_include_directories(
  testOpflowFunctionality PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

#[[ Configure OPFLOW functionality tests ]]
target_link_libraries(testOpflowFunctionality ${EXAGO_APP_LIBS})
if(EXAGO_INSTALL_TESTS)
  install(TARGETS testOpflowFunctionality
          RUNTIME DESTINATION tests/functionality
  )
endif()
