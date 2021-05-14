#[[

Finds HiOp include directory and libraries and exports target `HiOp`

User may set:
- HIOP_ROOT_DIR

]]

include(CheckCXXSymbolExists)

find_library(HIOP_LIBRARY
  NAMES
  libhiop.so hiop # Prefer shared over static since static causes issues
                  # linking with mpi on some platforms
  PATHS
  ${HIOP_DIR} $ENV{HIOP_DIR} ${HIOP_ROOT_DIR}
  ENV LD_LIBRARY_PATH ENV DYLD_LIBRARY_PATH
  PATH_SUFFIXES
  lib64 lib)

if(HIOP_LIBRARY)
  get_filename_component(HIOP_LIBRARY_DIR ${HIOP_LIBRARY} DIRECTORY)
endif()

find_path(
  HIOP_INCLUDE_DIR
  NAMES hiopInterface.hpp
  PATHS ${HIOP_DIR} ${HIOP_ROOT_DIR} $ENV{HIOP_DIR} ${HIOP_LIBRARY_DIR}/..
  PATH_SUFFIXES include
)

if(HIOP_LIBRARY AND HIOP_INCLUDE_DIR)
  message(STATUS "Found HiOp library: " ${HIOP_LIBRARY})
  message(STATUS "Found HiOp include directory: " ${HIOP_INCLUDE_DIR})

  add_library(HiOp INTERFACE)
  target_link_libraries(HiOp INTERFACE ${HIOP_LIBRARY})
  target_include_directories(HiOp INTERFACE ${HIOP_INCLUDE_DIR})

  if(EXAGO_ENABLE_GPU)
    include(FindMagma)
    target_link_libraries(HiOp INTERFACE Magma)
    target_include_directories(HiOp INTERFACE ${CUDA_TOOLKIT_INCLUDE})
  endif()

  check_cxx_symbol_exists(HIOP_SPARSE ${HIOP_INCLUDE_DIR}/hiop_defs.hpp EXAGO_ENABLE_HIOP_SPARSE)
  check_cxx_symbol_exists(HIOP_USE_COINHSL ${HIOP_INCLUDE_DIR}/hiop_defs.hpp EXAGO_HIOP_USE_COINHSL)  

  if(EXAGO_ENABLE_HIOP_SPARSE AND EXAGO_HIOP_USE_COINHSL)
    include(FindExaGOCOINHSL)
    if(NOT COINHSL_LIBRARY)
      message(
        FATAL_ERROR "HIOP_SPARSE is enabled, but COINHSL could not be found."
      )
    endif()
    target_link_libraries(HiOp INTERFACE COINHSL)
  endif()
else()
  if(NOT HIOP_LIB)
    message(STATUS "HiOp library not found! Please provide correct filepath.")
  endif()
  if(NOT HIOP_INCLUDE_DIR)
    message(
      STATUS "HiOp include directory  not found! Please provide correct path."
    )
  endif()
endif()

set(HIOP_INCLUDE_DIR CACHE PATH "Path to HiOp header files")
set(HIOP_LIBRARY CACHE FILEPATH "Path to HiOp library")
