#[[ Functionality tests for ExaGO PFLOW applications ]]

add_executable(testPflowFunctionality pflow.c pflowselfcheck.c)
target_include_directories(
  testPflowFunctionality PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

#[[ Configure PFLOW functionality tests ]]
target_link_libraries(testPflowFunctionality ${EXAGO_APP_LIBS})
if(EXAGO_INSTALL_TESTS)
  install(TARGETS testPflowFunctionality
          RUNTIME DESTINATION tests/functionality
  )
endif()
