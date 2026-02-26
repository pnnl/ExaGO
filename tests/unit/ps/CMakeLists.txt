configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/ieee9bus_v33.raw
  ${CMAKE_CURRENT_BINARY_DIR}/ieee9bus_v33.raw COPYONLY
)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/ieee9bus_v34_shunts.raw
  ${CMAKE_CURRENT_BINARY_DIR}/ieee9bus_v34_shunts.raw COPYONLY
)

add_executable(test_psse_parser test_psse_parser.cpp)
target_link_libraries(test_psse_parser ExaGO::PS)
add_test(NAME psse_parser COMMAND test_psse_parser)
