# Find Ipopt installation path
find_path (IPOPT_DIR NAMES share/coin/doc/Ipopt/ipopt_addlibs_cpp.txt HINTS ~/local/ipopt)
message (STATUS "Found Ipopt in ${IPOPT_DIR}")

# Find Ipopt header path and ensure all needed files are there
find_path(IPOPT_INCLUDE_DIR
  IpTNLP.hpp
  HINTS ${IPOPT_DIR}/include/coin
)
message (STATUS "Found Ipopt headers in ${IPOPT_INCLUDE_DIR}")

# Assume Ipopt lib directory is in the same place as the include directory
set(IPOPT_LIBRARY_DIR ${IPOPT_DIR}/lib)

# Ipopt modules needed for the build
# The order matters in case of static build!
set(IPOPT_MODULES 
  ipopt 
  coinmetis
  coinmumps
)

# Find each Ipopt module and add it to the list of libraries to link
set(IPOPT_LIBRARY)
foreach(mod ${IPOPT_MODULES})
  find_library(IPOPT_${mod}
    NAMES ${mod}
    HINTS ${IPOPT_LIBRARY_DIR}
  )
  if(IPOPT_${mod})
    set(IPOPT_LIBRARY ${IPOPT_LIBRARY} ${IPOPT_${mod}})
  else()
    # unset ${IPOPT_LIBRARY_DIR} and ask user to supply it
  endif()
endforeach()
message (STATUS "Found Ipopt libraries ${IPOPT_LIBRARY}")
