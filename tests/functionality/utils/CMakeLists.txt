add_library(testFunctionality selfcheckopflow.c selfcheckpflow.c selfchecksopflow.c)
target_link_libraries(testFunctionality opflow pflow sopflow)
target_include_directories(testFunctionality PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
