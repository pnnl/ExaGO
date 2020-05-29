

# Define sources for applications
# PFLOWPROTO
set(PFLOWPROTO_SRC
  pflow-proto.c
)

#OPFLOWPROTO
set(OPFLOWPROTO_SRC
  opflow-proto3.c
)

# PFLOW
set(PFLOW_SRC
  pflow-main.c
)

# PFLOW2
set(PFLOW2_SRC
  pflow-main2.c
)

# OPFLOW
set(OPFLOW_SRC
  opflow-main.c
)

#SCOPFLOW
set(SCOPFLOW_SRC
  scopflow-main.c
)

#TCOPFLOW
set(TCOPFLOW_SRC
  tcopflow-main.c
)

set(SCOPFLOW_PFLOW_APP_LIBS   
  ${SCOPFLOW_UTIL_LIBS}	
  ${SCOPFLOW_PS_LIBS}
  ${SCOPFLOW_PFLOW_LIBS}	 
  ${SCOPFLOW_MATH_LIB}
  PETSC::SNES
  MPI::MPI_C
)

# Add executables
# PFLOWPROTO
add_executable(app_pflowproto ${PFLOWPROTO_SRC})
# PFLOW
add_executable(app_pflow ${PFLOW_SRC})
# PFLOW2
add_executable(app_pflow2 ${PFLOW2_SRC})
# OPFLOW
add_executable(app_opflow ${OPFLOW_SRC})
# OPFLOWPROTO
add_executable(app_opflowproto ${OPFLOWPROTO_SRC})
#SCOPFLOW
add_executable(app_scopflow ${SCOPFLOW_SRC})
#TCOPFLOW
add_executable(app_tcopflow ${TCOPFLOW_SRC})

# Link libraries
target_link_libraries(app_pflowproto 
	${SCOPFLOW_PFLOW_APP_LIBS}
)

target_link_libraries(app_pflow 
	${SCOPFLOW_PFLOW_APP_LIBS}
)

target_link_libraries(app_pflow2 
	${SCOPFLOW_PFLOW_APP_LIBS}
)

target_link_libraries(app_opflow
	${SCOPFLOW_OPFLOW_LIBS} 
	${SCOPFLOW_PFLOW_APP_LIBS}
)

target_link_libraries(app_opflowproto
	${SCOPFLOW_OPFLOW_LIBS} 
	${SCOPFLOW_PFLOW_APP_LIBS}
)

target_link_libraries(app_scopflow
	${SCOPFLOW_LIBS}
	${SCOPFLOW_OPFLOW_LIBS}
	${SCOPFLOW_PFLOW_APP_LIBS}
)

target_link_libraries(app_tcopflow
	${TCOPFLOW_LIBS}
	${TCOPFLOW_OPFLOW_LIBS}
	${SCOPFLOW_PFLOW_APP_LIBS}
)

# Install the executables
install(TARGETS app_pflowproto RUNTIME DESTINATION bin)
install(TARGETS app_pflow RUNTIME DESTINATION bin)
install(TARGETS app_pflow2 RUNTIME DESTINATION bin)
install(TARGETS app_opflow RUNTIME DESTINATION bin)
install(TARGETS app_opflowproto RUNTIME DESTINATION bin)
install(TARGETS app_scopflow RUNTIME DESTINATION bin)
install(TARGETS app_tcopflow RUNTIME DESTINATION bin)




