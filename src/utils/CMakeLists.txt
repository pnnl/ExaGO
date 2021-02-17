
# Sources for utility functions shared across libraries

# Specify sources for utility functions shared across libraries
set(UTILS_SRC
  utils.c
  version.c
  driverutils.c
  scopflowutils.c
  comm.c
)

# Headers
set(UTILS_INCLUDE
  ${CMAKE_SOURCE_DIR}/include/common.h
  ${CMAKE_SOURCE_DIR}/include/utils.h
)

# Add library
add_library(psutil ${UTILS_SRC})

# Link libraries
target_link_libraries(psutil PRIVATE PETSC::SNES)

if(EXAGO_ENABLE_MPI)
  target_link_libraries(psutil PRIVATE MPI::MPI_C)
endif()

# Install library
install(TARGETS psutil DESTINATION lib)

# Install header
install(FILES ${UTILS_INCLUDE} DESTINATION include)
#EXAGO_EXAGO_
set(EXAGO_UTIL_LIBS ${EXAGO_UTIL_LIBS} psutil PARENT_SCOPE)
