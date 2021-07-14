#[[

Exports target `hiop_hip` which finds all hip libraries needed by hiop.

]]

find_library(
  HIPBLAS_LIBRARY
  NAMES hipblas
  PATHS ENV LD_LIBRARY_PATH ENV DYLD_LIBRARY_PATH
  PATH_SUFFIXES lib64 lib
)

if(HIPBLAS_LIBRARY)
  get_filename_component(HIPBLAS_LIBRARY_DIR ${HIPBLAS_LIBRARY} DIRECTORY)
endif()

find_path(
  HIPBLAS_INCLUDE_DIR
  NAMES hipblas.h
  PATHS ${HIPBLAS_LIBRARY_DIR}/..
  PATH_SUFFIXES include
)

if(HIPBLAS_LIBRARY)
  include_directories(${HIPBLAS_INCLUDE_DIR})
  add_library(Hipblas INTERFACE)
  target_link_libraries(Hipblas INTERFACE ${HIPBLAS_LIBRARY})
  target_include_directories(Hipblas INTERFACE ${HIPBLAS_INCLUDE_DIR})
  message(STATUS "Found Hipblas include: ${HIPBLAS_INCLUDE_DIR}")
  message(STATUS "Found Hipblas library: ${HIPBLAS_LIBRARY}")
else()
  message(STATUS "Hipblas was not found.")
endif()

set(HIPBLAS_INCLUDE_DIR CACHE PATH "Path to hipblas.h")
set(HIPBLAS_LIBRARY CACHE PATH "Path to hipblas library")
