

# Define sources for applications
# PFLOW
set(PFLOW_SRC
  pflow-main.cpp
)

# OPFLOW
set(OPFLOW_SRC
  opflow-main.cpp
)

#SCOPFLOW
set(SCOPFLOW_SRC
  scopflow-main.cpp
)

#SOPFLOW
set(SOPFLOW_SRC
  sopflow-main.cpp
)

#TCOPFLOW
set(TCOPFLOW_SRC
  tcopflow-main.cpp
)

set(EXAGO_PFLOW_APP_LIBS   
  ${EXAGO_UTIL_LIBS}	
  ${EXAGO_PS_LIBS}
  ${EXAGO_PFLOW_LIBS}	 
  ${EXAGO_MATH_LIB}
  PETSC::SNES
)

if(EXAGO_HAVE_MPI)
  list(APPEND EXAGO_PFLOW_APP_LIBS MPI::MPI_C)
endif()

# Add executables
# PFLOW
add_executable(app_pflow ${PFLOW_SRC})
# OPFLOW
add_executable(app_opflow ${OPFLOW_SRC})
#SCOPFLOW
add_executable(app_scopflow ${SCOPFLOW_SRC})
#SOPFLOW
add_executable(app_sopflow ${SOPFLOW_SRC})
#TCOPFLOW
add_executable(app_tcopflow ${TCOPFLOW_SRC})

# Link libraries
target_link_libraries(app_pflow 
	${EXAGO_PFLOW_APP_LIBS}
)

target_link_libraries(app_opflow
	${EXAGO_OPFLOW_LIBS} 
	${EXAGO_PFLOW_APP_LIBS}
)

target_link_libraries(app_scopflow
	${EXAGO_SCOPFLOW_LIBS}
	${EXAGO_TCOPFLOW_LIBS}
	${EXAGO_OPFLOW_LIBS}
	${EXAGO_PFLOW_APP_LIBS}
)

target_link_libraries(app_sopflow
	${EXAGO_SOPFLOW_LIBS}
	${EXAGO_SCOPFLOW_LIBS}
	${EXAGO_TCOPFLOW_LIBS}
	${EXAGO_OPFLOW_LIBS}
	${EXAGO_PFLOW_APP_LIBS}
)

target_link_libraries(app_tcopflow
	${EXAGO_TCOPFLOW_LIBS}
	${EXAGO_PFLOW_APP_LIBS}
)

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