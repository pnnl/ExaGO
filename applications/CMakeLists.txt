# Define sources for applications PFLOW
set(PFLOW_SRC pflow-main.c)

# OPFLOW
set(OPFLOW_SRC opflow-main.c)

# SCOPFLOW
set(SCOPFLOW_SRC scopflow-main.c)

# SOPFLOW
set(SOPFLOW_SRC sopflow-main.c)

# TCOPFLOW
set(TCOPFLOW_SRC tcopflow-main.c)

set(EXAGO_PFLOW_APP_LIBS ExaGO::UTILS ExaGO::PS ${EXAGO_PFLOW_LIBS}
                         ${EXAGO_MATH_LIB} PETSC::SNES
)

if(EXAGO_ENABLE_MPI)
  list(APPEND EXAGO_PFLOW_APP_LIBS MPI::MPI_C)
endif()

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
target_link_libraries(app_pflow ${EXAGO_PFLOW_APP_LIBS})

target_link_libraries(app_opflow ExaGO::OPFLOW ${EXAGO_PFLOW_APP_LIBS})

target_link_libraries(
  app_scopflow ExaGO::SCOPFLOW ExaGO::TCOPFLOW ExaGO::OPFLOW
  ${EXAGO_PFLOW_APP_LIBS}
)

target_link_libraries(
  app_sopflow ExaGO::SOPFLOW ExaGO::SCOPFLOW ExaGO::TCOPFLOW ExaGO::OPFLOW
  ${EXAGO_PFLOW_APP_LIBS}
)

target_link_libraries(app_tcopflow ExaGO::TCOPFLOW ${EXAGO_PFLOW_APP_LIBS})

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
