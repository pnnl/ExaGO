
#pragma once
#include "petscsys.h"
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

  LocalOrdinalType computeJacobian(PFLOW pflow, Mat JRef) {

    PetscErrorCode ierr;
    LocalOrdinalType fail = 0;
    Mat J;

    ierr = PFLOWCreateMatrix(pflow, &J);
    CHKERRQ(ierr);

    ierr = PFLOWGetJacobian(pflow, &J);
    CHKERRQ(ierr);

    fail += verifyAnswer(J, JRef);

    return (fail);
  }

  virtual int verifyAnswer(Mat a, Mat b, const RealType &tol = eps) const {
    int ncols, ncolsref;
    int fail = 0;
    const int *cols, *colsref;
    const double *vals, *valsref;
    PetscInt nrow, ncol, nrowref, ncolref;
    PetscErrorCode ierr;
    auto idx = [&ncol](double *mat, int r, int c) {
      return mat[(r * ncol) + c];
    };

    // MatView(a, 0);
    // MatView(b, 0);

    ierr = MatGetSize(a, &nrow, &ncol);
    CHKERRQ(ierr);
    ierr = MatGetSize(b, &nrowref, &ncolref);
    CHKERRQ(ierr);

    if (nrow != nrowref) {
      std::cout << "Failed due to row count: J: " << nrow
                << " JRef: " << nrowref << std::endl;
      fail++;
    }
    if (ncol != ncolref) {
      std::cout << "Failed due to column count: J: " << ncol
                << " JRef: " << ncolref << std::endl;
      fail++;
    }

    if (fail == 0) {
      for (int i = 0; i < nrow; i++) {
        ierr = MatGetRow(a, i, &ncols, &cols, &vals);
        CHKERRQ(ierr);
        ierr = MatGetRow(b, i, &ncolsref, &colsref, &valsref);
        for (int j = 0; j < ncols; j++) {
          // std::cout << "J: " << vals[j] << " JRef: " << valsref[j] <<
          // std::endl;
          if (!isEqual(vals[j], valsref[j], tol)) {
            std::cout << "Failed for index (" << i << ", " << cols[j]
                      << ") : " << vals[j] << " != " << valsref[j] << std::endl;
            fail++;
          }
        }
      }
    }

    return fail;
  }

}; // class TestOpflow : public TestBase

} // namespace tests
} // namespace exago
