# If HIOP_DIR is not set then try to find it at default locations
if(NOT HIOP_DIR)
   find_path(HIOP_DIR 
     NAMES 
     _dist-default-build
     lib/libhiop.a 
     HINTS
     /usr
     /usr/local
     /opt/local
     ~/local/ipopt
     # PNNL specific paths
     /qfs/projects/exasgd/newell/hiop
     /qfs/projects/exasgd/newell/hiop_shared_cpu
     /qfs/projects/exasgd/newell/hiop_master-cuda_10.2-gcc_7.4.0-metis_5.1.0-magma_2.5.2_cuda10.2-openmpi_3.1.5
     /qfs/projects/exasgd/newell/hiop_PIC
     /qfs/projects/exasgd/newell/hiop-cpu # PNNL power9 system-specific
   )
endif()

# Exit if HIOP_DIR is not set and cannot be found in default locations
if(NOT HIOP_DIR)
   message(FATAL_ERROR "HiOP directory could not be found. Please specify the directory using the flag -DHIOP_DIR=<location_of_HiOP_install>")
endif(NOT HIOP_DIR)

# Ipopt library location
find_library(HIOP_LIBRARY NAME libhiop.so libhiop.dylib libhiop.a
  HINTS 
  ${HIOP_DIR}/lib
  ${HIOP_DIR}/lib64
)
message(STATUS "HIOP library directory=${HIOP_LIBRARY}")

find_library(MAGMA_LIBRARY NAME libmagma.so libmagma.a libmagam.dylib magma.pc
  HINTS
  ${MAGMA_DIR}/lib
  ${MAGMA_DIR}/lib/pkgconfig
)
message(STATUS "MAGMA DIR=${MAGMA_DIR}")
message(STATUS "MAGMA library directory=${MAGMA_LIBRARY}")


# HiOP include directories
# Find HiOP header path and ensure all needed files are there
find_path(HIOP_INCLUDE_DIR NAME hiopNlpFormulation.hpp 
  HINTS 
  ${HIOP_DIR}/include
)
message(STATUS "HIOP Include directory = ${HIOP_INCLUDE_DIR}")

if(HIOP_INCLUDE_DIR AND HIOP_LIBRARY)
  include_directories(${HIOP_INCLUDE_DIR})
  set(SCOPFLOW_HAVE_HIOP 1)
else()
  message(FATAL_ERROR "HIOP not found!")
endif()
