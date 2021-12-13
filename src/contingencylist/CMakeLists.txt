# Specify PowerSystem module sources (shared across libraries)
set(CTGCLIST_SRC contingencylist.cpp)

set(CTGCLIST_INCLUDE ${CMAKE_SOURCE_DIR}/include/private/contingencylist.h)

set_source_files_properties(${CTGCLIST_SRC} PROPERTIES LANGUAGE CXX)

exago_add_library(
  CTGCLIST
  SOURCES
  ${CTGCLIST_SRC}
  HEADERS
  ${CTGCLIST_INCLUDE}
  LINK_LIBRARIES
  PUBLIC
  ExaGO::UTILS
  OUTPUT_NAME
  exago_ctgclist
)
