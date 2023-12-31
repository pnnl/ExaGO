set_source_files_properties(pflow_main.cpp PROPERTIES LANGUAGE CXX)

# Define sources for applications PFLOW
set(PFLOW_SRC pflow_main.cpp)
# Add executables PFLOW
add_executable(app_pflow ${PFLOW_SRC})
# Link libraries
target_link_libraries(app_pflow ExaGO::PFLOW)
set_target_properties(app_pflow PROPERTIES OUTPUT_NAME pflow)
install(TARGETS app_pflow RUNTIME DESTINATION bin)

if(EXAGO_ENABLE_IPOPT OR EXAGO_ENABLE_HIOP)

  set_source_files_properties(
    opflow_main.cpp scopflow_main.cpp sopflow_main.cpp tcopflow_main.cpp
    PROPERTIES LANGUAGE CXX
  )

  # OPFLOW
  set(OPFLOW_SRC opflow_main.cpp)
  add_executable(app_opflow ${OPFLOW_SRC})
  target_link_libraries(app_opflow ExaGO::OPFLOW)
  set_target_properties(app_opflow PROPERTIES OUTPUT_NAME opflow)
  install(TARGETS app_opflow RUNTIME DESTINATION bin)

  if(EXAGO_ENABLE_IPOPT)

    # SCOPFLOW
    set(SCOPFLOW_SRC scopflow_main.cpp)
    add_executable(app_scopflow ${SCOPFLOW_SRC})
    target_link_libraries(app_scopflow ExaGO::SCOPFLOW)
    set_target_properties(app_scopflow PROPERTIES OUTPUT_NAME scopflow)
    install(TARGETS app_scopflow RUNTIME DESTINATION bin)

    # SOPFLOW
    set(SOPFLOW_SRC sopflow_main.cpp)
    add_executable(app_sopflow ${SOPFLOW_SRC})
    target_link_libraries(app_sopflow ExaGO::SOPFLOW)
    set_target_properties(app_sopflow PROPERTIES OUTPUT_NAME sopflow)
    install(TARGETS app_sopflow RUNTIME DESTINATION bin)

    # TCOPFLOW
    set(TCOPFLOW_SRC tcopflow_main.cpp)
    add_executable(app_tcopflow ${TCOPFLOW_SRC})
    target_link_libraries(app_tcopflow ExaGO::TCOPFLOW)
    set_target_properties(app_tcopflow PROPERTIES OUTPUT_NAME tcopflow)
    install(TARGETS app_tcopflow RUNTIME DESTINATION bin)
  endif()
endif()
