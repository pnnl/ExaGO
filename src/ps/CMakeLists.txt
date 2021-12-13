# Specify PowerSystem module sources (shared across libraries)
set(PS_SRC ps.cpp psreaddata.cpp psoutput.cpp psislanding.cpp)

set(PS_HEADERS ${CMAKE_SOURCE_DIR}/include/ps.h)

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
