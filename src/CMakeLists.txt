add_subdirectory(utils)
add_subdirectory(ps)
add_subdirectory(pflow)

if(EXAGO_ENABLE_HIOP OR EXAGO_ENABLE_IPOPT)
  add_subdirectory(contingencylist)
  add_subdirectory(opflow)
endif()

# Disable these applications when IPOPT is not installed since these
# applications have a dependency on IPOPT. So, if the user installs HIOP and
# does not install IPOPT then these applications will be disabled
if(EXAGO_ENABLE_IPOPT)
  add_subdirectory(scopflow)
  add_subdirectory(sopflow)
  add_subdirectory(tcopflow)
endif()
