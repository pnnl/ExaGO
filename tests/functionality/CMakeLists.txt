#[[ Functionality tests for ExaGO applications ]]
add_subdirectory(pflow)
if(EXAGO_ENABLE_IPOPT OR EXAGO_ENABLE_HIOP)
  add_subdirectory(opflow)
endif()

# Disable tests for these applications when IPOPT is not installed because
# TCOPFLOW is dependent on IPOPT, SCOPFLOW and SOPFLOW also depend on IPOPT
# even when using HIOP solver.
if(EXAGO_ENABLE_IPOPT)
  add_subdirectory(tcopflow)
  add_subdirectory(scopflow)
  add_subdirectory(sopflow)
endif()

