#[[

Finds HiOp include directory and libraries and exports target `HiOp`

User may set:
- HIOP_ROOT_DIR

]]

include(CheckCXXSymbolExists)

find_library(
  HIOP_LIBRARY
  NAMES libhiop.so hiop # Prefer shared over static since static causes issues
                        # linking with mpi on some platforms
  PATHS ${HIOP_DIR}
        $ENV{HIOP_DIR}
        ${HIOP_ROOT_DIR}
        ENV
        LD_LIBRARY_PATH
        ENV
        DYLD_LIBRARY_PATH
  PATH_SUFFIXES lib64 lib
)

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

  set(HIOP_DEFINITIONS ${HIOP_INCLUDE_DIR}/hiop_defs.hpp)
  message(STATUS "HiOp definitions: ${HIOP_DEFINITIONS}")

  # Force cmake to use mpi flags if we're building with mpi. Otherwise, the
  # underlying call to `try_compile` will fail and the the symbol will be
  # reported as missing even if it's defined in the header.
  if(EXAGO_ENABLE_MPI)
    message(STATUS "Adding MPI compile flags for HiOp sparse/coinhsl detection")
    set(CMAKE_REQUIRED_FLAGS -I${MPI_C_INCLUDE_DIRS} ${MPI_C_COMPILE_OPTIONS}
                             -I${MPI_CXX_INCLUDE_DIRS} ${MPI_CXX_COMPILE_OPTIONS})
    set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
  endif()

  check_cxx_symbol_exists(
    HIOP_SPARSE ${HIOP_DEFINITIONS} EXAGO_ENABLE_HIOP_SPARSE
  )

  check_cxx_symbol_exists(
    HIOP_USE_COINHSL ${HIOP_DEFINITIONS} EXAGO_HIOP_USE_COINHSL
  )

  if(EXAGO_ENABLE_GPU)
    include(FindMagma)

    if(EXAGO_ENABLE_CUDA)
      target_include_directories(HiOp INTERFACE ${CUDA_TOOLKIT_INCLUDE})
    elseif(EXAGO_ENABLE_HIP)
      include(FindHipblas)
      target_link_libraries(Magma INTERFACE Hipblas)
    endif()

    target_link_libraries(HiOp INTERFACE Magma)
  endif()

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
