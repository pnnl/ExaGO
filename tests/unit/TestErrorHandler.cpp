#include <petsc.h>
#include <utils.hpp>

#include "utils/testBase.hpp"
using namespace exago::tests;

/* These _Throw*_ functions test various constructors for ExaGOErrors */

/* Throw an ExaGOError with the default constructor */
static void ThrowExaGOErrorDefaultConstructor() noexcept(false)
{
  throw ExaGOError();
}

/* Throw an ExaGOError with a string */
static void ThrowExaGOErrorStringConstructor() noexcept(false)
{
  throw ExaGOError("Some error");
}

/* Throw an ExaGOError from a PetscErrorCode */
static void ThrowExaGOErrorFromPetscError() noexcept(false)
{
  /* The error doesn't actually matter here - we're just ensuring that
   * constructing an ExaGOError from a petsc error code is functional */
  PetscErrorCode ierr = PETSC_ERR_MEM;
  throw ExaGOError(ierr);
}

//------------------------------------------------------------------------------

struct TestExaGOErrorHandler : public TestBase {
  static int TestExaGOErrorDefaultConstructor()
  {
    static int fail = 0;
    try {
      ThrowExaGOErrorDefaultConstructor();
      /* If this line isn't reached, something went wrong with the ExaGOError */
      fail++;
    } catch (const ExaGOError& e) {
      /* Success! */
    }
    printMessage(fail, __func__, 0);
    return fail;
  }

  static int TestExaGOErrorStringConstructor()
  {
    static int fail = 0;
    try {
      ThrowExaGOErrorStringConstructor();
      /* If this line isn't reached, something went wrong with the ExaGOError */
      fail++;
    } catch (const ExaGOError& e) {
      /* Success! */
    }
    printMessage(fail, __func__, 0);
    return fail;
  }

  static int TestExaGOErrorPetscErrorCodeConstructor()
  {
    static int fail = 0;
    try {
      ThrowExaGOErrorFromPetscError();
      /* If this line isn't reached, something went wrong with the ExaGOError */
      fail++;
    } catch (const ExaGOError& e) {
      if (!e.IsPetscError()) fail++;
    }
    printMessage(fail, __func__, 0);
    return fail;
  }

  static int TestExaGOErrorCheckMacro()
  {
    static int fail = 0;
    try {
      /* A common idiom is to retrieve the return code from an ExaGO or PETSc
       * function call and run an error checking macro on the return code before
       * continuing execution. We provide a similar macro in _utils.hpp_ which
       * throws an ExaGOError via the PetscErrorCode initializer if the code
       * falls within PetscError's range of values. */

      PetscErrorCode ierr = PETSC_ERR_MEM; /* pretend this value is returned from
                                              a Petsc function call */
      ExaGOCheckError(ierr); /* Check the error code just like we already do with
                                the CHKERRQ family of Petsc macros */
    } catch (const ExaGOError& e) {
      if (!e.IsPetscError()) fail++;
    }
    printMessage(fail, __func__, 0);
    return fail;
  }

  static int RunAllTests() {
    static int fail = 0;
    fail += TestExaGOErrorDefaultConstructor();
    fail += TestExaGOErrorStringConstructor();
    fail += TestExaGOErrorPetscErrorCodeConstructor();
    fail += TestExaGOErrorCheckMacro();
    return fail;
  }
};

int main()
{
  return TestExaGOErrorHandler().RunAllTests();
}
