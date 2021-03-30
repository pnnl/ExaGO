# Specify PowerSystem module sources (shared across libraries)
set(PS_SRC ps.c psreaddata.c psoutput.c psislanding.c)

set(PS_INCLUDE ${CMAKE_SOURCE_DIR}/include/ps.h)

# Compile library
add_library(ps ${PS_SRC})

target_link_libraries(ps PRIVATE PETSC::SNES ${EXAGO_UTIL_LIBS})

if(EXAGO_ENABLE_MPI)
  target_link_libraries(ps PRIVATE MPI::MPI_C)
endif()

# Install library
install(TARGETS ps DESTINATION lib)

# Install header
install(FILES ${PS_INCLUDE} DESTINATION include)

set(EXAGO_PS_LIBS
    ${EXAGO_PS_LIBS} ps
    PARENT_SCOPE
)
