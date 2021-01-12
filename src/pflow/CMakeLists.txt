# Specify PowerSystem module sources (shared across libraries)
set(PFLOW_SRC
  pflow.cpp
)

# Header files for pflow
set(PFLOW_INCLUDE
  ${CMAKE_SOURCE_DIR}/include/pflow.h
)

# Compile library
add_library(pflow ${PFLOW_SRC})

target_link_libraries(pflow PRIVATE
			 PETSC::SNES 
			 ${EXAGO_UTIL_LIBS}
			 ${EXAGO_PS_LIBS})

if(EXAGO_HAVE_MPI)
  target_link_libraries(pflow PRIVATE MPI::MPI_C)
endif()

# Install libraries
install(TARGETS pflow DESTINATION lib)

#Install headers
install(FILES ${PFLOW_INCLUDE} DESTINATION include)

set(EXAGO_PFLOW_LIBS ${EXAGO_PFLOW_LIBS} pflow PARENT_SCOPE)
