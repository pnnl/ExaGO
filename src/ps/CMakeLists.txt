# Specify PowerSystem module sources (shared across libraries)
set(PS_SRC ps.c psreaddata.c psoutput.c psislanding.c)

set(PS_HEADERS ${CMAKE_SOURCE_DIR}/include/ps.h)

exago_add_library(
  PS
  SOURCES
  ${PS_SRC}
  HEADERS
  ${PS_INCLUDE}
  LINK_LIBRARIES
  PRIVATE
  PETSC::SNES
  PRIVATE
  ExaGO::UTILS
  PRIVATE
  $<$<BOOL:${EXAGO_ENABLE_MPI}>:MPI::MPI_C>
  OUTPUT_NAME
  exago_ps
)
