# If IPOPT_DIR is not set then try to find it at default locations
if(NOT IPOPT_DIR)
   find_path(IPOPT_DIR NAMES share/coin/doc/Ipopt/ipopt_addlibs_cpp.txt HINTS ~/local/ipopt)
endif()

# Exit if IPOPT_DIR is not set and cannot be found in default locations
if(NOT IPOPT_DIR)
   message(FATAL_ERROR "Ipopt directory could not be found. Please specify the directory using the flag -DIPOPT_DIR=<location_of_Ipopt_install>")
endif(NOT IPOPT_DIR)

# Ipopt library location
find_library(IPOPT_LIBRARY NAME libipopt.a libipopt.so libipopt.dylib HINTS ${IPOPT_DIR}/lib)
message(STATUS "IPOPT library directory=${IPOPT_LIBRARY}")

# Ipopt include directories
# Find Ipopt header path and ensure all needed files are there
find_path(IPOPT_INCLUDE_DIR NAME IpIpoptNLP.hpp HINTS ${IPOPT_DIR}/include/coin)
message(STATUS "IPOPT Include directory = ${IPOPT_INCLUDE_DIR}")

if(IPOPT_INCLUDE_DIR AND IPOPT_LIBRARY)
  include_directories(${IPOPT_INCLUDE_DIR})
  set(SCOPFLOW_HAVE_IPOPT TRUE)
else()
  message(FATAL_ERROR "IPOPT not found!")
endif()