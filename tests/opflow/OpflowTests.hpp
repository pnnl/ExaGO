#pragma once
#include <type_traits>
#include <cassert>
#include <testBase.hpp>

#include <opflow.h>
#include <common.h>
#include <private/opflowimpl.h>
#include <scopflow_config.h>

#define cleanup(fail, opflow) \
  printMessage(fail, __func__, getRank(opflow)); \
  return reduceReturn(fail, opflow);

namespace exago { namespace tests {

class TestOpflow : public TestBase
{
public:
  TestOpflow() = default;

  LocalOrdinalType computeObjective(OPFLOW opflow, Vec X, RealType obj_ref)
  {
    PetscErrorCode   ierr;
    LocalOrdinalType fail = 0;
    RealType         obj_val;
    ierr = OPFLOWComputeObjective(opflow,X,&obj_val);CHKERRQ(ierr);
    fail += !isEqual(obj_val, obj_ref);
    cleanup(fail, opflow);
  }

  LocalOrdinalType computeGradient(OPFLOW opflow, Vec X, Vec grad_ref)
  {
    PetscErrorCode   ierr;
    Vec              grad;
    LocalOrdinalType fail = 0;

    ierr = VecDuplicate(grad_ref,&grad);CHKERRQ(ierr);
    ierr = OPFLOWComputeGradient(opflow,X,grad);CHKERRQ(ierr);
    fail += verifyAnswer(grad, grad_ref);
    ierr = VecDestroy(&grad);CHKERRQ(ierr);

    cleanup(fail, opflow);
  }

  LocalOrdinalType computeVariableBounds(OPFLOW opflow,Vec Xlref,Vec Xuref)
  {
    PetscErrorCode   ierr;
    LocalOrdinalType fail = 0;
    Vec              Xu;
    Vec              Xl;

    ierr = VecDuplicate(Xlref,&Xl);CHKERRQ(ierr);
    ierr = VecDuplicate(Xuref,&Xu);CHKERRQ(ierr);
    ierr = OPFLOWComputeVariableBounds(opflow,Xl,Xu);CHKERRQ(ierr);
    fail += verifyAnswer(Xl, Xlref);
    fail += verifyAnswer(Xu, Xuref);

    ierr = VecDestroy(&Xl);CHKERRQ(ierr);
    ierr = VecDestroy(&Xu);CHKERRQ(ierr);
    cleanup(fail, opflow);
  }

  LocalOrdinalType computeConstraints(OPFLOW opflow,Vec X,Vec Gref)
  {
    PetscErrorCode ierr;
    LocalOrdinalType fail = 0;
    Vec G;
    
    ierr = VecDuplicate(Gref,&G);CHKERRQ(ierr);
    ierr = OPFLOWComputeConstraints(opflow,X,G);CHKERRQ(ierr);
    fail += verifyAnswer(G, Gref);

    ierr = VecDestroy(&G);CHKERRQ(ierr);
    cleanup(fail, opflow);
  }

  LocalOrdinalType computeConstraintBounds(OPFLOW opflow,Vec Glref,Vec Guref)
  {
    PetscErrorCode ierr;
    LocalOrdinalType fail = 0;
    Vec Gl, Gu;
    RealType dglnorm, dgunorm;
    
    ierr = VecDuplicate(Glref,&Gl);CHKERRQ(ierr);
    ierr = VecDuplicate(Guref,&Gu);CHKERRQ(ierr);
    ierr = OPFLOWComputeConstraintBounds(opflow,Gl,Gu);CHKERRQ(ierr);
    fail += verifyAnswer(Gl, Glref);
    fail += verifyAnswer(Gu, Guref);

    ierr = VecDestroy(&Gl);CHKERRQ(ierr);
    ierr = VecDestroy(&Gu);CHKERRQ(ierr);
    cleanup(fail, opflow);
  }

private:
  virtual MPI_Comm getMPIComm(OPFLOW opflow) const
  {
    if(opflow == nullptr) THROW_NULL_DEREF;
    return opflow->comm->type;
  }

  virtual LocalOrdinalType getRank(OPFLOW opflow) const
  {
    if(opflow == nullptr) THROW_NULL_DEREF;
    LocalOrdinalType rank;
    MPI_Comm_rank(getMPIComm(opflow), &rank);
    return rank;
  }

  /**
   * @brief Base implementation for reduceReturn. Most other implementations of
   * this method will grab the MPI communicator from the respective type and
   * call this function to perform the reduction.
   *
   */
  LocalOrdinalType reduceReturn(LocalOrdinalType fail, MPI_Comm comm) const
  {
    LocalOrdinalType gFail;
    int ierr = MPI_Allreduce(&fail, &gFail, 1, MPI_INT, MPI_SUM, comm);assert(ierr==0);
    return static_cast<LocalOrdinalType>(gFail != 0);
  }

  virtual LocalOrdinalType reduceReturn(LocalOrdinalType fail, OPFLOW opflow) const
  {
    if(opflow == nullptr) THROW_NULL_DEREF;
    MPI_Comm comm = getMPIComm(opflow);
    return reduceReturn(fail, comm);
  }

  virtual LocalOrdinalType reduceReturn(LocalOrdinalType fail, PetscObject obj) const
  {
    MPI_Comm       comm;
    PetscErrorCode ierr;
    ierr = PetscObjectGetComm(obj, &comm);CHKERRQ(ierr);
    return reduceReturn(fail, comm);
  }

  /**
   * @brief Verifies that two vectors are equivilant within tolerance.
   *
   */
  virtual int verifyAnswer(Vec a, Vec b) const
  {
    LocalOrdinalType   size;
    LocalOrdinalType   fail = 0;
    PetscErrorCode     ierr;
    PetscScalar const* avals;
    PetscScalar const* bvals;
    ierr = VecGetSize(a, &size);CHKERRQ(ierr);
    ierr = VecGetArrayRead(a, &avals);CHKERRQ(ierr);
    ierr = VecGetArrayRead(b, &bvals);CHKERRQ(ierr);
    for (int i=0; i<size; i++)
    {
      if (!isEqual(avals[i], bvals[i]))
        fail++;
    }
    ierr = VecRestoreArrayRead(a, &avals);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(b, &bvals);CHKERRQ(ierr);
    return reduceReturn(fail, (PetscObject)a);
  }
};

}} // namespace exago::tests
