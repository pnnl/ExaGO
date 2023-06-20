#include <exago_config.h>

#if defined(EXAGO_ENABLE_RAJA)

#pragma once

#include "pbpolrajahiop.h"

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

#endif // EXAGO_ENABLE_HIOP
