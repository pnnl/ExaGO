#include <exago_config.h>

#if defined(EXAGO_ENABLE_RAJA)
#if defined(EXAGO_ENABLE_HIOP_SPARSE)

#pragma once

#include <umpire/Allocator.hpp>
#include <umpire/ResourceManager.hpp>

#include <opflow.h>
#include "pbpolrajahiopsparse.hpp"
#include "paramsrajahiop.h"

/**

  @brief Takes a raw pointer to some data and registers it with the
  allocator and resource manager.

  TODO: Ensure that we do not have to deregister the allocation with
  the resource manager. This is done with `resmgr.deregisterAllocation(T* ptr)`,
  but must we do this with ever allocation? Will we incurr a memory leak if not?

*/
template <typename T, typename SizeType>
void registerWith(T *ptr, SizeType N, umpire::ResourceManager &resmgr,
                  umpire::Allocator &allocator) {
  umpire::util::AllocationRecord record{ptr, sizeof(T) * N,
                                        allocator.getAllocationStrategy()};
  resmgr.registerAllocation(ptr, record);
}

#endif // EXAGO_ENABLE_HIOP_SPARSE
#endif // EXAGO_ENABLE_HIOP_RAJA
