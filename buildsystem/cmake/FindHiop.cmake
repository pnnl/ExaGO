#[[

Finds HiOp include directory and libraries and exports target `HiOp`

User may set:
- HIOP_ROOT_DIR

]]

include(CheckCXXSymbolExists)

find_package(HiOp REQUIRED HINTS ${HIOP_ROOT_DIR})

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
else()
  message(FATAL_ERROR "Find_package could not load HiOp")
endif()
