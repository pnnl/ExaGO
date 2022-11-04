#[[ Functionality tests for ExaGO applications ]]
add_subdirectory(pflow)
if(EXAGO_ENABLE_IPOPT OR EXAGO_ENABLE_HIOP)
  add_subdirectory(opflow)
  add_subdirectory(tcopflow)
  add_subdirectory(scopflow)
  add_subdirectory(sopflow)
endif()
