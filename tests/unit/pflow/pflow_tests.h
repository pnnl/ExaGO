
#pragma once
#include <cassert>
#include <test_base.h>
#include <type_traits>

#include <common.h>
#include <exago_config.h>
#include <pflow.h>
#include <private/pflowimpl.h>

//#if defined(EXAGO_ENABLE_RAJA)
//#include <RAJA/RAJA.hpp>
//#include <umpire/Allocator.hpp>
//#include <umpire/ResourceManager.hpp>
//#endif

//#define cleanup(fail, opflow)                                                  \
//  printMessage(fail, __func__, getRank(opflow));                               \
//  return reduceReturn(fail, opflow);

namespace exago {
namespace tests {

class TestPflow : public TestBase {
public:
  TestPflow() = default;

  LocalOrdinalType trivialTest(PFLOW pflow) {
    PetscErrorCode ierr;
    LocalOrdinalType fail = 0;
    RealType obj_val;
    return(fail);  
}

}; // class TestOpflow : public TestBase

} // namespace tests
} // namespace exago
