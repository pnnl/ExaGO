#[[

Calls find_package to search for HiOp. Requires HiOp >= 0.5.3

User may set:
- HiOp_DIR

]]

find_package(HiOp REQUIRED)

if(TARGET HiOp::HiOp)
  if(HiOp::SPARSE AND TARGET HiOp::COINHSL)
    set(EXAGO_ENABLE_HIOP_SPARSE
        ON
        CACHE BOOL "Enable HiOp Sparse" FORCE
    )
  else()
    set(EXAGO_ENABLE_HIOP_SPARSE
        OFF
        CACHE BOOL "Enable HiOp Sparse" FORCE
    )
  endif()
  mark_as_advanced(FORCE HiOp::SPARSE)
else()
  message(FATAL_ERROR "Find_package could not load HiOp")
endif()
