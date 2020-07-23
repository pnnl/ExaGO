
set(EXAGO_PFLOW_APP_LIBS   
  ${EXAGO_UTIL_LIBS}	
  ${EXAGO_PS_LIBS}
  ${EXAGO_PFLOW_LIBS}	 
  ${EXAGO_MATH_LIB}
  PETSC::SNES
  MPI::MPI_CXX
)

add_executable(testAcopf TestAcopf.cpp)
target_link_libraries(testAcopf
	${EXAGO_OPFLOW_LIBS} 
	${EXAGO_PFLOW_APP_LIBS}
)
target_include_directories(testAcopf PRIVATE ./utils)

install(TARGETS testAcopf DESTINATION tests)
