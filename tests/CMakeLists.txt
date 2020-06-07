
set(SCOPFLOW_PFLOW_APP_LIBS   
  ${SCOPFLOW_UTIL_LIBS}	
  ${SCOPFLOW_PS_LIBS}
  ${SCOPFLOW_PFLOW_LIBS}	 
  ${SCOPFLOW_MATH_LIB}
  PETSC::SNES
  MPI::MPI_C
)


add_executable(testAcopf TestAcopf.cpp)
target_link_libraries(testAcopf
	${SCOPFLOW_OPFLOW_LIBS} 
	${SCOPFLOW_PFLOW_APP_LIBS}
)

install(TARGETS testAcopf RUNTIME DESTINATION tests)
