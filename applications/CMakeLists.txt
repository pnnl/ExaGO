set_source_files_properties(
  pflow_main.cpp opflow_main.cpp scopflow_main.cpp sopflow_main.cpp
  tcopflow_main.cpp PROPERTIES LANGUAGE CXX
)

# Define sources for applications PFLOW
set(PFLOW_SRC pflow_main.cpp)

# OPFLOW
set(OPFLOW_SRC opflow_main.cpp)

# SCOPFLOW
set(SCOPFLOW_SRC scopflow_main.cpp)

# SOPFLOW
set(SOPFLOW_SRC sopflow_main.cpp)

# TCOPFLOW
set(TCOPFLOW_SRC tcopflow_main.cpp)

# Add executables PFLOW
add_executable(app_pflow ${PFLOW_SRC})
# OPFLOW
add_executable(app_opflow ${OPFLOW_SRC})
# SCOPFLOW
add_executable(app_scopflow ${SCOPFLOW_SRC})
# SOPFLOW
add_executable(app_sopflow ${SOPFLOW_SRC})
# TCOPFLOW
add_executable(app_tcopflow ${TCOPFLOW_SRC})

# Link libraries
target_link_libraries(app_pflow ExaGO::PFLOW)
target_link_libraries(app_opflow ExaGO::OPFLOW)
target_link_libraries(app_scopflow ExaGO::SCOPFLOW)
target_link_libraries(app_sopflow ExaGO::SOPFLOW)
target_link_libraries(app_tcopflow ExaGO::TCOPFLOW)

set_target_properties(app_pflow PROPERTIES OUTPUT_NAME pflow)
set_target_properties(app_opflow PROPERTIES OUTPUT_NAME opflow)
set_target_properties(app_scopflow PROPERTIES OUTPUT_NAME scopflow)
set_target_properties(app_sopflow PROPERTIES OUTPUT_NAME sopflow)
set_target_properties(app_tcopflow PROPERTIES OUTPUT_NAME tcopflow)

# Install the executables
install(TARGETS app_pflow RUNTIME DESTINATION bin)
install(TARGETS app_opflow RUNTIME DESTINATION bin)
install(TARGETS app_scopflow RUNTIME DESTINATION bin)
install(TARGETS app_sopflow RUNTIME DESTINATION bin)
install(TARGETS app_tcopflow RUNTIME DESTINATION bin)
