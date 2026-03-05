add_executable(test_psse_parser test_psse_parser.cpp)
target_link_libraries(test_psse_parser ExaGO::PS)
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data)
add_test(NAME psse_parser
  COMMAND test_psse_parser
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data
)
