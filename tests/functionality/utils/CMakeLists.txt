add_library(testFunctionality selfcheckopflow.c selfcheckpflow.c selfchecksopflow.c selfchecktcopflow.c)
target_link_libraries(testFunctionality opflow pflow sopflow tcopflow)
target_include_directories(testFunctionality PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
