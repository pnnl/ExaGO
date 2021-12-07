# Specify PowerSystem module sources (shared across libraries)
set(PS_SRC ps.c psreaddata.c psoutput.c psislanding.c)

set(PS_HEADERS ${CMAKE_SOURCE_DIR}/include/ps.h ${CMAKE_SOURCE_DIR})

set_source_files_properties(${PS_SRC} PROPERTIES LANGUAGE CXX)

exago_add_library(
  PS
  SOURCES
  ${PS_SRC}
  HEADERS
  ${PS_INCLUDE}
  LINK_LIBRARIES
  PUBLIC
  ExaGO::UTILS
  OUTPUT_NAME
  exago_ps
)
