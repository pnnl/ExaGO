
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

  LocalOrdinalType computeObjective(OPFLOW opflow, double *xref, RealType obj_ref)
  {
    PetscErrorCode   ierr;
    LocalOrdinalType fail = 0;
    RealType         obj_val;
    ierr = (*opflow->modelops.computeobjectivearray)(opflow, xref, &obj_val);CHKERRQ(ierr);
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

  LocalOrdinalType computeGradient(OPFLOW opflow, double *xref, double *grad_ref)
  {
    PetscErrorCode   ierr;
    double*          grad;
    LocalOrdinalType fail = 0;
    int              nx,nconeq,nconineq;
    
    ierr = OPFLOWGetSizes(opflow,&nx,&nconeq,&nconineq);CHKERRQ(ierr);

    ierr = PetscMalloc1(nx,&grad);CHKERRQ(ierr);

    ierr = (*opflow->modelops.computegradientarray)(opflow, xref, grad);CHKERRQ(ierr);

    fail += verifyAnswer(grad_ref, grad, nx);

    ierr = PetscFree(grad);CHKERRQ(ierr);

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

    LocalOrdinalType computeVariableBounds(OPFLOW opflow,double *xlref,double *xuref)
  {
    PetscErrorCode   ierr;
    LocalOrdinalType fail = 0;
    double           *xu,*xl;
    int              nx,nconeq,nconineq;
    
    ierr = OPFLOWGetSizes(opflow,&nx,&nconeq,&nconineq);CHKERRQ(ierr);

    ierr = PetscMalloc1(nx,&xu);CHKERRQ(ierr);
    ierr = PetscMalloc1(nx,&xl);CHKERRQ(ierr);

    ierr = (*opflow->modelops.setvariableboundsarray)(opflow,xl,xu);CHKERRQ(ierr);

    fail += verifyAnswer(xlref,xl,nx);
    fail += verifyAnswer(xuref,xu,nx);

    ierr = PetscFree(xl);CHKERRQ(ierr);
    ierr = PetscFree(xu);CHKERRQ(ierr);

    cleanup(fail,opflow);
  }

  LocalOrdinalType computeConstraints(OPFLOW opflow,Vec X,Vec Gref)
  {
    PetscErrorCode ierr;
    LocalOrdinalType fail = 0;
    Vec G;
    
    ierr = VecDuplicate(Gref,&G);CHKERRQ(ierr);
    ierr = OPFLOWComputeConstraints(opflow,X,G);CHKERRQ(ierr);
    fail += verifyAnswer(G, Gref,1e-3);

    ierr = VecDestroy(&G);CHKERRQ(ierr);
    cleanup(fail, opflow);
  }

  LocalOrdinalType computeConstraints(OPFLOW opflow, double *xref, double *gref)
  {
    PetscErrorCode ierr;
    LocalOrdinalType fail = 0;
    double *g, *ge, *gi;

    int              nx,nconeq,nconineq;
    
    ierr = OPFLOWGetSizes(opflow,&nx,&nconeq,&nconineq);CHKERRQ(ierr);

    g = new double[nconeq + nconineq];
    ge = g;
    if(opflow->nconineq) gi = g + opflow->nconeq;

    ierr = (*opflow->modelops.computeequalityconstraintsarray)(opflow, xref, ge);CHKERRQ(ierr);
    if(opflow->nconineq) {
      ierr = (*opflow->modelops.computeinequalityconstraintsarray)(opflow, xref, gi);CHKERRQ(ierr);
    }

    fail += verifyAnswer(g, gref, nconeq + nconineq,1e-3);

    delete[] g;

    cleanup(fail, opflow);
  }

  LocalOrdinalType computeConstraintBounds(OPFLOW opflow,Vec Glref,Vec Guref)
  {
    PetscErrorCode ierr;
    LocalOrdinalType fail = 0;
    Vec Gl, Gu;
    
    ierr = VecDuplicate(Glref,&Gl);CHKERRQ(ierr);
    ierr = VecDuplicate(Guref,&Gu);CHKERRQ(ierr);
    ierr = OPFLOWComputeConstraintBounds(opflow,Gl,Gu);CHKERRQ(ierr);
    fail += verifyAnswer(Gl, Glref);
    fail += verifyAnswer(Gu, Guref);

    ierr = VecDestroy(&Gl);CHKERRQ(ierr);
    ierr = VecDestroy(&Gu);CHKERRQ(ierr);
    cleanup(fail, opflow);
  }

   LocalOrdinalType computeConstraintBounds(OPFLOW opflow,double *glref,double *guref)
  {
    PetscErrorCode ierr;
    LocalOrdinalType fail = 0;
    double *gl, *gu;
    int nx, nconeq, nconineq;

    ierr = OPFLOWGetSizes(opflow, &nx, &nconeq, &nconineq);CHKERRQ(ierr);

    gl = new double[nconeq + nconineq];
    gu = new double[nconeq + nconineq];

    ierr = (*opflow->modelops.setconstraintboundsarray)(opflow,gl,gu);CHKERRQ(ierr);
    fail += verifyAnswer(gl, glref, nconeq + nconineq);
    fail += verifyAnswer(gu, guref, nconeq + nconineq);

    delete[] gl;
    delete[] gu;

    cleanup(fail, opflow);
  }

  LocalOrdinalType computeConstraintJacobian(OPFLOW opflow,Vec X,Mat Jeqref, Mat Jineqref)
  {
    PetscErrorCode   ierr;
    LocalOrdinalType fail = 0;
    Mat              Jeq, Jineq;
    Vec              temp1, temp2;
    PetscScalar      alpha = -1.0;
    PetscInt         nrow, ncol;
    const RealType   local_tol = 1e-4; // tolerance specific to current test

    ierr = MatDuplicate(Jeqref,MAT_DO_NOT_COPY_VALUES,&Jeq);CHKERRQ(ierr);
    ierr = MatDuplicate(Jineqref,MAT_DO_NOT_COPY_VALUES,&Jineq);CHKERRQ(ierr);
    ierr = OPFLOWComputeConstraintJacobian(opflow,X,Jeq,Jineq);CHKERRQ(ierr);

    ierr = MatGetSize(Jeq,&nrow,&ncol);CHKERRQ(ierr);
    ierr = VecCreate(MPI_COMM_SELF,&temp1);CHKERRQ(ierr);
    ierr = VecSetSizes(temp1,nrow,nrow);CHKERRQ(ierr);
    ierr = VecSetFromOptions(temp1);

    ierr = VecDuplicate(temp1,&temp2);CHKERRQ(ierr);

    ierr = VecSet(temp1,0.0);
    ierr = VecSet(temp2,0.0);

    ierr = MatMult(Jeq,X,temp1);CHKERRQ(ierr);
    ierr = MatMult(Jeqref,X,temp2);CHKERRQ(ierr);

    fail += verifyAnswer(temp1, temp2, local_tol);

    ierr = VecDestroy(&temp1);CHKERRQ(ierr);
    ierr = VecDestroy(&temp2);CHKERRQ(ierr);

    ierr = MatGetSize(Jineq,&nrow,&ncol);CHKERRQ(ierr);
    ierr = VecCreate(MPI_COMM_SELF,&temp1);CHKERRQ(ierr);
    ierr = VecSetSizes(temp1,nrow,nrow);CHKERRQ(ierr);
    ierr = VecSetFromOptions(temp1);

    ierr = VecDuplicate(temp1,&temp2);CHKERRQ(ierr);

    ierr = VecSet(temp1,0.0);
    ierr = VecSet(temp2,0.0);

    ierr = MatMult(Jineq,X,temp1);CHKERRQ(ierr);
    ierr = MatMult(Jineqref,X,temp2);CHKERRQ(ierr);

    fail += verifyAnswer(temp1, temp2, local_tol);

    ierr = VecDestroy(&temp1);CHKERRQ(ierr);
    ierr = VecDestroy(&temp2);CHKERRQ(ierr);

    ierr = MatDestroy(&Jeq);CHKERRQ(ierr);
    ierr = MatDestroy(&Jineq);CHKERRQ(ierr);

    cleanup(fail, opflow);
  }

  LocalOrdinalType computeHessian(OPFLOW opflow,Vec X,Vec Lambda,PetscScalar obj_factor,Mat Hessref)
  {
    PetscErrorCode ierr;
    LocalOrdinalType fail = 0;
    Mat Hess;
    Vec temp1,temp2;
    
    ierr = MatDuplicate(Hessref,MAT_DO_NOT_COPY_VALUES,&Hess);CHKERRQ(ierr);

    ierr = OPFLOWComputeHessian(opflow,X,Lambda,obj_factor,Hess);CHKERRQ(ierr);

    ierr = VecDuplicate(X,&temp1);CHKERRQ(ierr);
    ierr = VecDuplicate(temp1,&temp2);CHKERRQ(ierr);

    ierr = VecSet(temp1,0.0);
    ierr = VecSet(temp2,0.0);

    ierr = MatMult(Hess,X,temp1);CHKERRQ(ierr);
    ierr = MatMult(Hessref,X,temp2);CHKERRQ(ierr);

    fail += verifyAnswer(temp1, temp2);

    ierr = VecDestroy(&temp1);CHKERRQ(ierr);
    ierr = VecDestroy(&temp2);CHKERRQ(ierr);

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
  virtual int verifyAnswer(const double* avals, const double* bvals, int size,const RealType& tol=eps) const
  {
    LocalOrdinalType   fail = 0;

    for (int i=0; i<size; i++)
    {
      if (!isEqual(avals[i], bvals[i], tol))
        fail++;
    }

    return(fail);
    //    return reduceReturn(fail, (PetscObject)a);
  }


  /**
   * @brief Verifies that two vectors are equivilant within tolerance.
   *
   */
  virtual int verifyAnswer(Vec a, Vec b, const RealType& tol=eps) const
  {
    LocalOrdinalType   size;
    LocalOrdinalType   fail = 0;
    PetscErrorCode     ierr;
    PetscScalar const* avals;
    PetscScalar const* bvals;
    ierr = VecGetSize(a, &size);CHKERRQ(ierr);
    ierr = VecGetArrayRead(a, &avals);CHKERRQ(ierr);
    ierr = VecGetArrayRead(b, &bvals);CHKERRQ(ierr);

    fail = verifyAnswer(avals,bvals,size,tol);

    ierr = VecRestoreArrayRead(a, &avals);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(b, &bvals);CHKERRQ(ierr);
    return reduceReturn(fail, (PetscObject)a);
  }
};

}} // namespace exago::tests
