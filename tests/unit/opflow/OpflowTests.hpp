
#pragma once
#include <type_traits>
#include <cassert>
#include <testBase.hpp>

#include <opflow.h>
#include <common.h>
#include <private/opflowimpl.h>
#include <exago_config.h>

#if defined(EXAGO_ENABLE_RAJA)
#include <umpire/Allocator.hpp>
#include <umpire/ResourceManager.hpp>
#include <RAJA/RAJA.hpp>
#endif

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

  // Note: xref is allocated on the GPU when running the GPU tests
  LocalOrdinalType computeObjective(OPFLOW opflow, double *xref, RealType obj_ref)
  {
    PetscErrorCode   ierr;
    LocalOrdinalType fail = 0;
    RealType         obj_val=0.0;
    ierr = OPFLOWComputeObjectiveArray(opflow, xref, &obj_val);CHKERRQ(ierr);
    fail += !isEqual(obj_val, obj_ref);
    cleanup(fail, opflow);
  }

  LocalOrdinalType computeGradient(OPFLOW opflow, Vec X, Vec grad_ref)
  {
    PetscErrorCode   ierr;
    Vec              grad;
    LocalOrdinalType fail = 0;

    ierr = VecDuplicate(grad_ref,&grad);CHKERRQ(ierr);
    ierr = VecSet(grad,0.0);CHKERRQ(ierr);
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

    ierr = OPFLOWComputeGradientArray(opflow,xref,grad);CHKERRQ(ierr);

    fail += verifyAnswer(grad_ref, grad, nx);

    ierr = PetscFree(grad);CHKERRQ(ierr);

    cleanup(fail, opflow);
  }

#if defined(EXAGO_ENABLE_RAJA)
  LocalOrdinalType computeGradient(OPFLOW opflow, double *xref_dev, double *grad_ref,umpire::ResourceManager& resmgr)
  {
    PetscErrorCode   ierr;
    double*          grad;
    double*          grad_dev;
    LocalOrdinalType fail = 0;
    int              nx,nconeq,nconineq;

    // Get allocator
    umpire::Allocator h_allocator = resmgr.getAllocator("HOST");

    ierr = OPFLOWGetSizes(opflow,&nx,&nconeq,&nconineq);CHKERRQ(ierr);

    grad = static_cast<double*>(h_allocator.allocate(nx*sizeof(double)));

#ifdef EXAGO_ENABLE_GPU
    umpire::Allocator d_allocator = resmgr.getAllocator("DEVICE");
    grad_dev = static_cast<double*>(d_allocator.allocate(nx*sizeof(double)));
#else
    grad_dev = grad;
#endif

    ierr = OPFLOWComputeGradientArray(opflow,xref_dev,grad_dev);CHKERRQ(ierr);

    // Copy back from device to host
    resmgr.copy(grad,grad_dev);

    fail += verifyAnswer(grad_ref, grad, nx);

    h_allocator.deallocate(grad);
#ifdef EXAGO_ENABLE_GPU
    d_allocator.deallocate(grad_dev);
#endif
    cleanup(fail, opflow);
  }
#endif


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

#if defined(EXAGO_ENABLE_RAJA) 
  LocalOrdinalType computeVariableBounds(OPFLOW opflow,double *xlref,double *xuref,umpire::ResourceManager& resmgr)
  {
    PetscErrorCode   ierr;
    LocalOrdinalType fail = 0;
    double           *xl,*xu;
    double           *xu_dev,*xl_dev;
    int              nx,nconeq,nconineq;
    
    // Get allocator
    umpire::Allocator h_allocator = resmgr.getAllocator("HOST");

    ierr = OPFLOWGetSizes(opflow,&nx,&nconeq,&nconineq);CHKERRQ(ierr);

    // Create arrays on the device
    xl = static_cast<double*>(h_allocator.allocate(nx*sizeof(double)));
    xu = static_cast<double*>(h_allocator.allocate(nx*sizeof(double)));
#ifdef EXAGO_ENABLE_GPU
    umpire::Allocator d_allocator = resmgr.getAllocator("DEVICE");
    xl_dev = static_cast<double*>(d_allocator.allocate(nx*sizeof(double)));
    xu_dev = static_cast<double*>(d_allocator.allocate(nx*sizeof(double)));
#else
    xl_dev = xl;
    xu_dev = xu;
#endif
    ierr = (*opflow->modelops.setvariableboundsarray)(opflow,xl_dev,xu_dev);CHKERRQ(ierr);

    // Copy back from device to host
    resmgr.copy(xl,xl_dev);
    resmgr.copy(xu,xu_dev);

    fail += verifyAnswer(xlref,xl,nx);
    fail += verifyAnswer(xuref,xu,nx);

    h_allocator.deallocate(xl);
    h_allocator.deallocate(xu);
#ifdef EXAGO_ENABLE_GPU
    d_allocator.deallocate(xl_dev);
    d_allocator.deallocate(xu_dev);
#endif
    cleanup(fail,opflow);
  }
#endif

  LocalOrdinalType computeConstraints(OPFLOW opflow,Vec X,Vec Gref)
  {
    PetscErrorCode ierr;
    LocalOrdinalType fail = 0;
    Vec G;
    
    ierr = VecDuplicate(Gref,&G);CHKERRQ(ierr);
    ierr = VecSet(G,0.0);CHKERRQ(ierr);
    ierr = OPFLOWComputeConstraints(opflow,X,G);CHKERRQ(ierr);
    fail += verifyAnswer(G, Gref);

    ierr = VecDestroy(&G);CHKERRQ(ierr);
    cleanup(fail, opflow);
  }

  LocalOrdinalType computeConstraints(OPFLOW opflow, double *xref, double *gref)
  {
    PetscErrorCode ierr;
    LocalOrdinalType fail = 0;
    double *g, *ge, *gi;
    int nx,nconeq,nconineq;
    
    ierr = OPFLOWGetSizes(opflow,&nx,&nconeq,&nconineq);CHKERRQ(ierr);

    g = new double[nconeq + nconineq];
    ge = g;
    if(opflow->nconineq) gi = g + opflow->nconeq;

    ierr = OPFLOWComputeEqualityConstraintsArray(opflow,xref,ge);CHKERRQ(ierr);

    if(opflow->nconineq) {
      ierr = OPFLOWComputeInequalityConstraintsArray(opflow,xref,gi);CHKERRQ(ierr);
    }

    fail += verifyAnswer(g, gref, nconeq + nconineq);

    delete[] g;

    cleanup(fail, opflow);
  }

#if defined(EXAGO_ENABLE_RAJA)
  LocalOrdinalType computeConstraints(OPFLOW opflow, double *xref_dev, double *gref,umpire::ResourceManager& resmgr)
  {
    PetscErrorCode   ierr;
    LocalOrdinalType fail = 0;
    double           *g;
    double           *g_dev,*ge_dev,*gi_dev;
    int              nx,nconeq,nconineq,i;

    // Get allocator
    umpire::Allocator h_allocator = resmgr.getAllocator("HOST");
    
    ierr = OPFLOWGetSizes(opflow,&nx,&nconeq,&nconineq);CHKERRQ(ierr);

    g     = static_cast<double*>(h_allocator.allocate(opflow->ncon*sizeof(double)));
#ifdef EXAGO_ENABLE_GPU
    umpire::Allocator d_allocator = resmgr.getAllocator("DEVICE");
    g_dev = static_cast<double*>(d_allocator.allocate(opflow->ncon*sizeof(double)));
#else
    g_dev = g;
#endif

    for(i=0; i < opflow->ncon; i++) g[i] = 0.0;

    ge_dev = g_dev;
    if(opflow->nconineq) {
      gi_dev = g_dev + opflow->nconeq;
    }

    ierr = OPFLOWComputeEqualityConstraintsArray(opflow,xref_dev,ge_dev);CHKERRQ(ierr);

    if(opflow->nconineq) {
      ierr = OPFLOWComputeInequalityConstraintsArray(opflow,xref_dev,gi_dev);CHKERRQ(ierr);
    }
   
    // Copy back from device to host
    resmgr.copy(g,g_dev);

    fail += verifyAnswer(g, gref, nconeq + nconineq);

    h_allocator.deallocate(g);
#ifdef EXAGO_ENABLE_GPU
    d_allocator.deallocate(g_dev);
#endif
    cleanup(fail, opflow);
  }
#endif

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

#if defined(EXAGO_ENABLE_RAJA)
  LocalOrdinalType computeConstraintBounds(OPFLOW opflow,double *glref,double *guref,umpire::ResourceManager& resmgr)
  {
    PetscErrorCode ierr;
    LocalOrdinalType fail = 0;
    double *gl, *gu;
    double *gl_dev,*gu_dev;
    int nx, nconeq, nconineq,i;

    // Get allocator
    umpire::Allocator h_allocator = resmgr.getAllocator("HOST");

    ierr = OPFLOWGetSizes(opflow, &nx, &nconeq, &nconineq);CHKERRQ(ierr);

    gl     = static_cast<double*>(h_allocator.allocate(opflow->ncon*sizeof(double)));
    gu     = static_cast<double*>(h_allocator.allocate(opflow->ncon*sizeof(double)));

#ifdef EXAGO_ENABLE_GPU
    umpire::Allocator d_allocator = resmgr.getAllocator("DEVICE");
    gl_dev = static_cast<double*>(d_allocator.allocate(opflow->ncon*sizeof(double)));
    gu_dev = static_cast<double*>(d_allocator.allocate(opflow->ncon*sizeof(double)));
#else
    gl_dev = gl;
    gu_dev = gu;
#endif

    for(i=0; i < opflow->ncon; i++) {
      gl[i] = 0.0;
      gu[i] = 0.0;
    }

    // Copy from host to device
    resmgr.copy(gl_dev, gl);
    resmgr.copy(gu_dev, gu);

    ierr = (*opflow->modelops.setconstraintboundsarray)(opflow,gl_dev,gu_dev);CHKERRQ(ierr);

    // Copy back from device to host
    resmgr.copy(gl,gl_dev);
    resmgr.copy(gu,gu_dev);

    fail += verifyAnswer(gl, glref, nconeq + nconineq);
    fail += verifyAnswer(gu, guref, nconeq + nconineq);

    h_allocator.deallocate(gl);
    h_allocator.deallocate(gu);
#ifdef EXAGO_ENABLE_GPU
    d_allocator.deallocate(gl_dev);
    d_allocator.deallocate(gu_dev);
#endif
    cleanup(fail, opflow);
  }
#endif

  LocalOrdinalType computeConstraintJacobian(OPFLOW opflow,Vec X,Mat Jeqref, Mat Jineqref)
  {
    PetscErrorCode   ierr;
    LocalOrdinalType fail = 0;
    Mat              Jeq, Jineq;
    Vec              temp1, temp2;
    PetscInt         nrow, ncol;

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

    fail += verifyAnswer(temp1, temp2);

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

    fail += verifyAnswer(temp1, temp2);

    ierr = VecDestroy(&temp1);CHKERRQ(ierr);
    ierr = VecDestroy(&temp2);CHKERRQ(ierr);

    ierr = MatDestroy(&Jeq);CHKERRQ(ierr);
    ierr = MatDestroy(&Jineq);CHKERRQ(ierr);

    cleanup(fail, opflow);
  }

  LocalOrdinalType computeConstraintJacobian(OPFLOW opflow,double *x,Mat Jeqref_nat, Mat Jineqref_nat)
  {
    PetscErrorCode   ierr;
    LocalOrdinalType fail = 0;
    Mat              Jeqref_spdn, Jeqref_sparse, Jeqref_dense;
    PetscInt         nrow, ncol;
    int              *row_map, *col_map, *spcol_map, *dncol_map, *idxn2sd_map;
    int              nxsparse = 2*opflow->ps->ngenON;
    int              nxdense = 2*opflow->ps->nbus;

    // Re-order the reference solution using IS (Index sets), and then slice for the sparse and dense components
    ierr = OPFLOWGetVariableOrdering(opflow,&idxn2sd_map);CHKERRQ(ierr);
    ierr = MatGetSize(Jeqref_nat, &nrow, &ncol);CHKERRQ(ierr);

    IS col_map_IS, row_map_IS, spcol_map_IS, dncol_map_IS;

    col_map = new int[ncol];

    for(int i = 0; i < ncol; i++) col_map[idxn2sd_map[i]] = i;

    row_map = new int[nrow];
    for (int i = 0; i < nrow; i++) row_map[i] = i;

    spcol_map = new int[nxsparse];
    for (int i = 0; i < nxsparse; i++) spcol_map[i] = i;

    dncol_map = new int[nxdense];
    for (int i = 0; i < nxdense; i++) dncol_map[i] = i + nxsparse;

    ierr = ISCreateGeneral(PETSC_COMM_WORLD, nrow, row_map, PETSC_COPY_VALUES, &row_map_IS);CHKERRQ(ierr);
    ierr = ISCreateGeneral(PETSC_COMM_WORLD, ncol, col_map, PETSC_COPY_VALUES, &col_map_IS);CHKERRQ(ierr);
    ierr = ISCreateGeneral(PETSC_COMM_WORLD, nxsparse, spcol_map, PETSC_COPY_VALUES, &spcol_map_IS);CHKERRQ(ierr);
    ierr = ISCreateGeneral(PETSC_COMM_WORLD, nxdense, dncol_map, PETSC_COPY_VALUES, &dncol_map_IS);CHKERRQ(ierr);

    ierr = MatPermute(Jeqref_nat, row_map_IS, col_map_IS, &Jeqref_spdn);CHKERRQ(ierr);

    ierr = MatCreateSubMatrix(Jeqref_spdn, row_map_IS, spcol_map_IS, MAT_INITIAL_MATRIX, &Jeqref_sparse);CHKERRQ(ierr);
    ierr = MatCreateSubMatrix(Jeqref_spdn, row_map_IS, dncol_map_IS, MAT_INITIAL_MATRIX, &Jeqref_dense);CHKERRQ(ierr);

    // Get the hiop sparse matrix solution to test against
    int *iRow, *jCol;
    double *values;
    int nnz = nxsparse;

    ierr   = PetscMalloc1(nnz,&iRow);CHKERRQ(ierr);
    ierr   = PetscMalloc1(nnz,&jCol);CHKERRQ(ierr);
    ierr = PetscMalloc1(nnz,&values);CHKERRQ(ierr);

    ierr = (*opflow->modelops.computesparsejacobianhiop)(opflow, iRow, jCol, values);CHKERRQ(ierr);

    fail += verifyAnswer(Jeqref_sparse, nxsparse,iRow, jCol, values);

    // Get the hiop dense matrix solution to test against
    double *Jeq_dense;

    ierr = PetscMalloc1(nrow*nxdense,&Jeq_dense);CHKERRQ(ierr);

    ierr = (*opflow->modelops.computedenseequalityconstraintjacobianhiop)(opflow, x, Jeq_dense);CHKERRQ(ierr);

    fail += verifyAnswer(Jeqref_dense, Jeq_dense);

    // Checking for the presence of inequality constraints on the given problem
    if(opflow->Nconineq) {
      // Since the inequality contraints have a different number of rows, we must re-create some things first
      IS ineq_row_map_IS;
      int *ineq_row_map;
      PetscInt nrow_ineq, ncol_ineq;

      ierr = MatGetSize(Jineqref_nat, &nrow_ineq, &ncol_ineq);CHKERRQ(ierr);

      ineq_row_map = new int[nrow_ineq];
      for (int i = 0; i < nrow_ineq; i++) {
	ineq_row_map[i] = i;
      }

      ierr = ISCreateGeneral(PETSC_COMM_WORLD, nrow_ineq, ineq_row_map, PETSC_COPY_VALUES, &ineq_row_map_IS);CHKERRQ(ierr);

      // Re-order and get dense components
      Mat Jineqref_spdn, Jineqref_dense;
      ierr = MatPermute(Jineqref_nat, ineq_row_map_IS, col_map_IS, &Jineqref_spdn);CHKERRQ(ierr);
      ierr = MatCreateSubMatrix(Jineqref_spdn, ineq_row_map_IS, dncol_map_IS, MAT_INITIAL_MATRIX, &Jineqref_dense);CHKERRQ(ierr);
      
      // Get the hiop dense matrix solution
      double *Jineq_dense;

      ierr = PetscMalloc1(nrow_ineq*nxdense,&Jineq_dense);CHKERRQ(ierr);

      ierr = (*opflow->modelops.computedenseinequalityconstraintjacobianhiop)(opflow, x, Jineq_dense);CHKERRQ(ierr);

      fail += verifyAnswer(Jineqref_dense, Jineq_dense);

      //Cleanup                                                                                                         
      delete[] ineq_row_map;
      ierr = MatDestroy(&Jineqref_spdn);CHKERRQ(ierr);
      ierr = MatDestroy(&Jineqref_dense);CHKERRQ(ierr);

      ierr = PetscFree(Jineq_dense);CHKERRQ(ierr);
    }

    // Cleanup                    
    ierr = MatDestroy(&Jeqref_spdn);CHKERRQ(ierr);
    ierr = MatDestroy(&Jeqref_sparse);CHKERRQ(ierr);
    ierr = MatDestroy(&Jeqref_dense);CHKERRQ(ierr);

    ierr = ISDestroy(&col_map_IS);CHKERRQ(ierr);
    ierr = ISDestroy(&row_map_IS);CHKERRQ(ierr);
    ierr = ISDestroy(&spcol_map_IS);CHKERRQ(ierr);
    ierr = ISDestroy(&dncol_map_IS);CHKERRQ(ierr);

    
    delete[] col_map;
    delete[] row_map;
    delete[] spcol_map;
    delete[] dncol_map;

    ierr = PetscFree(iRow);CHKERRQ(ierr);
    ierr = PetscFree(jCol);CHKERRQ(ierr);
    ierr = PetscFree(values);CHKERRQ(ierr);

    ierr = PetscFree(Jeq_dense);CHKERRQ(ierr);

    cleanup(fail, opflow);
  }

#if defined(EXAGO_ENABLE_RAJA)
  /**
   * @brief Specific test for computing the constraint Jacobian using HiOp
   * 
   * @param x is in sparse-dense ordering
   * @param Jeqref is in application ordering
   * @param Jineqref is in application ordering
   * 
   * @pre The flag _opflow_->nconineq determines the presence of inequality constraints.
   *      If the flag is set, checking against Jineqref_nat can be skipped entirely.
   * @pre _Jeqref_nat_ and _Jineqref_nat_ have the same number of columns, but not the same number of rows
   */
  LocalOrdinalType computeConstraintJacobian(OPFLOW opflow,double *x_dev,Mat Jeqref_nat, Mat Jineqref_nat,umpire::ResourceManager& resmgr)
  {
    PetscErrorCode   ierr;
    LocalOrdinalType fail = 0;
    Mat              Jeqref_spdn, Jeqref_sparse, Jeqref_dense;
    PetscInt         nrow, ncol;
    int              *row_map, *col_map, *spcol_map, *dncol_map, *idxn2sd_map;
    int              nxsparse = 2*opflow->ps->ngenON;
    int              nxdense = 2*opflow->ps->nbus;

    // Get allocators
    umpire::Allocator h_allocator = resmgr.getAllocator("HOST");

    // Re-order the reference solution using IS (Index sets), and then slice for the sparse and dense components
    ierr = OPFLOWGetVariableOrdering(opflow,&idxn2sd_map);CHKERRQ(ierr);
    ierr = MatGetSize(Jeqref_nat, &nrow, &ncol);CHKERRQ(ierr);

    IS col_map_IS, row_map_IS, spcol_map_IS, dncol_map_IS;

    col_map = new int[ncol];

    for(int i = 0; i < ncol; i++) col_map[idxn2sd_map[i]] = i;

    row_map = new int[nrow];
    for (int i = 0; i < nrow; i++) row_map[i] = i;

    spcol_map = new int[nxsparse];
    for (int i = 0; i < nxsparse; i++) spcol_map[i] = i;

    dncol_map = new int[nxdense];
    for (int i = 0; i < nxdense; i++) dncol_map[i] = i + nxsparse;

    ierr = ISCreateGeneral(PETSC_COMM_WORLD, nrow, row_map, PETSC_COPY_VALUES, &row_map_IS);CHKERRQ(ierr);
    ierr = ISCreateGeneral(PETSC_COMM_WORLD, ncol, col_map, PETSC_COPY_VALUES, &col_map_IS);CHKERRQ(ierr);
    ierr = ISCreateGeneral(PETSC_COMM_WORLD, nxsparse, spcol_map, PETSC_COPY_VALUES, &spcol_map_IS);CHKERRQ(ierr);
    ierr = ISCreateGeneral(PETSC_COMM_WORLD, nxdense, dncol_map, PETSC_COPY_VALUES, &dncol_map_IS);CHKERRQ(ierr);

    ierr = MatPermute(Jeqref_nat, row_map_IS, col_map_IS, &Jeqref_spdn);CHKERRQ(ierr);

    ierr = MatCreateSubMatrix(Jeqref_spdn, row_map_IS, spcol_map_IS, MAT_INITIAL_MATRIX, &Jeqref_sparse);CHKERRQ(ierr);
    ierr = MatCreateSubMatrix(Jeqref_spdn, row_map_IS, dncol_map_IS, MAT_INITIAL_MATRIX, &Jeqref_dense);CHKERRQ(ierr);

    // Get the hiop sparse matrix solution to test against
    int *iRow, *jCol,*iRow_dev,*jCol_dev;
    double *values,*values_dev;
    int nnz = nxsparse;

    iRow = static_cast<int*>(h_allocator.allocate(nnz*sizeof(int)));
    jCol = static_cast<int*>(h_allocator.allocate(nnz*sizeof(int)));
    values = static_cast<double*>(h_allocator.allocate(nnz*sizeof(double)));

#ifdef EXAGO_ENABLE_GPU
    umpire::Allocator d_allocator = resmgr.getAllocator("DEVICE");
    iRow_dev = static_cast<int*>(d_allocator.allocate(nnz*sizeof(int)));
    jCol_dev = static_cast<int*>(d_allocator.allocate(nnz*sizeof(int)));
    values_dev = static_cast<double*>(d_allocator.allocate(nnz*sizeof(double)));
#else
    iRow_dev = iRow;
    jCol_dev = jCol;
    values_dev = values;
#endif

    ierr = (*opflow->modelops.computesparsejacobianhiop)(opflow, iRow_dev, jCol_dev, values_dev);CHKERRQ(ierr);

    resmgr.copy(iRow,iRow_dev);
    resmgr.copy(jCol,jCol_dev);
    resmgr.copy(values,values_dev);

    fail += verifyAnswer(Jeqref_sparse, nnz, iRow, jCol, values);

    // Get the hiop dense matrix solution to test against
    double *Jeq_dense,*Jeq_dense_dev;

    Jeq_dense = static_cast<double*>(h_allocator.allocate(nrow*nxdense*sizeof(double*)));
#ifdef EXAGO_ENABLE_GPU
    Jeq_dense_dev  = static_cast<double*>(d_allocator.allocate(nrow*nxdense*sizeof(double*)));
#else
    Jeq_dense_dev = Jeq_dense;
#endif

    ierr = (*opflow->modelops.computedenseequalityconstraintjacobianhiop)(opflow, x_dev, Jeq_dense_dev);CHKERRQ(ierr);

    // Copy back from device
    resmgr.copy(Jeq_dense, Jeq_dense_dev);

    fail += verifyAnswer(Jeqref_dense, Jeq_dense);

    // Checking for the presence of inequality constraints on the given problem
    if(opflow->Nconineq) {
      // Since the inequality contraints have a different number of rows, we must re-create some things first
      IS ineq_row_map_IS;
      int *ineq_row_map;
      PetscInt nrow_ineq, ncol_ineq;

      ierr = MatGetSize(Jineqref_nat, &nrow_ineq, &ncol_ineq);CHKERRQ(ierr);

      ineq_row_map = new int[nrow_ineq];
      for (int i = 0; i < nrow_ineq; i++) 
      {
	      ineq_row_map[i] = i;
      }

      ierr = ISCreateGeneral(PETSC_COMM_WORLD, nrow_ineq, ineq_row_map, PETSC_COPY_VALUES, &ineq_row_map_IS);CHKERRQ(ierr);

      // Re-order and get dense components
      Mat Jineqref_spdn, Jineqref_dense;
      ierr = MatPermute(Jineqref_nat, ineq_row_map_IS, col_map_IS, &Jineqref_spdn);CHKERRQ(ierr);
      ierr = MatCreateSubMatrix(Jineqref_spdn, ineq_row_map_IS, dncol_map_IS, MAT_INITIAL_MATRIX, &Jineqref_dense);CHKERRQ(ierr);
      
      // Get the hiop dense matrix solution
      double *Jineq_dense,*Jineq_dense_dev;

      Jineq_dense = static_cast<double*>(h_allocator.allocate(nxdense*nrow_ineq*sizeof(double*)));
#ifdef EXAGO_ENABLE_GPU
      Jineq_dense_dev  = static_cast<double*>(d_allocator.allocate(nxdense*nrow_ineq*sizeof(double*)));
#else
      Jineq_dense_dev = Jineq_dense;
#endif

      ierr = (*opflow->modelops.computedenseinequalityconstraintjacobianhiop)(opflow, x_dev, Jineq_dense_dev);CHKERRQ(ierr);

	    // Copy back from the device
	    resmgr.copy(Jineq_dense, Jineq_dense_dev);

      fail += verifyAnswer(Jineqref_dense, Jineq_dense);

      //Cleanup                                                                                                         
      delete[] ineq_row_map;
      ierr = MatDestroy(&Jineqref_spdn);CHKERRQ(ierr);
      ierr = MatDestroy(&Jineqref_dense);CHKERRQ(ierr);

	    h_allocator.deallocate(Jineq_dense);
#ifdef EXAGO_ENABLE_GPU
	    d_allocator.deallocate(Jineq_dense_dev);
#endif
    }

    // Cleanup                    
    ierr = MatDestroy(&Jeqref_spdn);CHKERRQ(ierr);
    ierr = MatDestroy(&Jeqref_sparse);CHKERRQ(ierr);
    ierr = MatDestroy(&Jeqref_dense);CHKERRQ(ierr);

    ierr = ISDestroy(&col_map_IS);CHKERRQ(ierr);
    ierr = ISDestroy(&row_map_IS);CHKERRQ(ierr);
    ierr = ISDestroy(&spcol_map_IS);CHKERRQ(ierr);
    ierr = ISDestroy(&dncol_map_IS);CHKERRQ(ierr);
    delete[] col_map;
    delete[] row_map;
    delete[] spcol_map;
    delete[] dncol_map;

    h_allocator.deallocate(iRow);
    h_allocator.deallocate(jCol);
    h_allocator.deallocate(values);
    h_allocator.deallocate(Jeq_dense);
#ifdef EXAGO_ENABLE_GPU
    d_allocator.deallocate(iRow_dev);
    d_allocator.deallocate(jCol_dev);
    d_allocator.deallocate(values_dev);
    d_allocator.deallocate(Jeq_dense_dev);
#endif
    cleanup(fail, opflow);
  }
#endif

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
  /**
   * @brief Specific test for computing the constraint Hessian using HiOp
   * 
   * @param x is in sparse-dense ordering
   * @param Hessref is in application ordering
   * 
   * @pre _Hessref_ is a square matrix
   */
  LocalOrdinalType computeHessian(OPFLOW opflow,double *x_ref,double *lambda_ref,PetscScalar obj_factor,Mat Hessref)
  {
    PetscErrorCode   ierr;
    LocalOrdinalType fail = 0;
    PetscInt         nrow, ncol;
    int              nxsparse = 2*opflow->ps->ngenON;
    int              nxdense = 2*opflow->ps->nbus;
    int              *permutation_map, *spcol_map, *dncol_map, *idxn2sd_map;
    Mat              Hessref_spdn, Hessref_sparse, Hessref_dense;

    // Re-order the reference solution using IS (Index sets), and then slice for the sparse and dense components
    ierr = OPFLOWGetVariableOrdering(opflow,&idxn2sd_map);CHKERRQ(ierr);
    ierr = MatGetSize(Hessref, &nrow, &ncol);CHKERRQ(ierr);

    IS permutation_map_IS, spcol_map_IS, dncol_map_IS;

    permutation_map = new int[ncol];
    for(int i = 0; i < ncol; i++) {
      permutation_map[idxn2sd_map[i]] = i;
    }

    spcol_map = new int[nxsparse];
    for (int i = 0; i < nxsparse; i++) {
      spcol_map[i] = i;
    }

    dncol_map = new int[nxdense];
    for (int i = 0; i < nxdense; i++) {
      dncol_map[i] = i + nxsparse;
    }

    ierr = ISCreateGeneral(PETSC_COMM_WORLD, nrow, permutation_map, PETSC_COPY_VALUES, &permutation_map_IS);CHKERRQ(ierr);
    ierr = ISCreateGeneral(PETSC_COMM_WORLD, nxsparse, spcol_map, PETSC_COPY_VALUES, &spcol_map_IS);CHKERRQ(ierr);
    ierr = ISCreateGeneral(PETSC_COMM_WORLD, nxdense, dncol_map, PETSC_COPY_VALUES, &dncol_map_IS);CHKERRQ(ierr);

    ierr = MatPermute(Hessref, permutation_map_IS, permutation_map_IS, &Hessref_spdn);CHKERRQ(ierr);

    ierr = MatCreateSubMatrix(Hessref_spdn, spcol_map_IS, spcol_map_IS, MAT_INITIAL_MATRIX, &Hessref_sparse);CHKERRQ(ierr);
    ierr = MatCreateSubMatrix(Hessref_spdn, dncol_map_IS, dncol_map_IS, MAT_INITIAL_MATRIX, &Hessref_dense);CHKERRQ(ierr);

    // Test the dense component
    double *hess_dense;
    hess_dense = new double[nxdense*nxdense];
    ierr = (*opflow->modelops.computedensehessianhiop)(opflow, x_ref, lambda_ref, hess_dense);CHKERRQ(ierr);

    fail += verifyAnswer(Hessref_dense, hess_dense);

    // Test the sparse component
    int *iRow, *jCol;
    double *values;
    int nnz = nxsparse;
    iRow = new int[nnz]();
    jCol = new int[nnz]();
    values = new double[nnz]();

    opflow->obj_factor = obj_factor;
    ierr = (*opflow->modelops.computesparsehessianhiop)(opflow, x_ref, iRow, jCol, values);CHKERRQ(ierr);

    fail += verifyAnswer(Hessref_sparse, nnz, iRow, jCol, values);

    // Cleanup
    delete[] permutation_map;
    delete[] spcol_map;
    delete[] dncol_map;
    delete[] iRow;
    delete[] jCol;
    delete[] values;
    ierr = MatDestroy(&Hessref_spdn);CHKERRQ(ierr);
    ierr = MatDestroy(&Hessref_sparse);CHKERRQ(ierr);
    ierr = MatDestroy(&Hessref_dense);CHKERRQ(ierr);    
    ierr = ISDestroy(&permutation_map_IS);CHKERRQ(ierr);
    ierr = ISDestroy(&spcol_map_IS);CHKERRQ(ierr);
    ierr = ISDestroy(&dncol_map_IS);CHKERRQ(ierr);
    cleanup(fail, opflow);
  }

#if defined(EXAGO_ENABLE_RAJA)
  LocalOrdinalType computeHessian(OPFLOW opflow,double *x_ref_dev,double *lambda_ref_dev,PetscScalar obj_factor,Mat Hessref,umpire::ResourceManager& resmgr,double* hess_dense,double* hess_dense_dev)
  {
    PetscErrorCode   ierr;
    LocalOrdinalType fail = 0;
    PetscInt         nrow, ncol;
    int              nxsparse = 2*opflow->ps->ngenON;
    int              nxdense = 2*opflow->ps->nbus;
    int              *permutation_map, *spcol_map, *dncol_map, *idxn2sd_map;
    Mat              Hessref_spdn, Hessref_sparse, Hessref_dense;

    // Get allocator
    umpire::Allocator h_allocator = resmgr.getAllocator("HOST");

    // Re-order the reference solution using IS (Index sets), and then slice for the sparse and dense components
    ierr = OPFLOWGetVariableOrdering(opflow,&idxn2sd_map);CHKERRQ(ierr);
    ierr = MatGetSize(Hessref, &nrow, &ncol);CHKERRQ(ierr);

    IS permutation_map_IS, spcol_map_IS, dncol_map_IS;

    permutation_map = new int[ncol];
    for(int i = 0; i < ncol; i++) {
      permutation_map[idxn2sd_map[i]] = i;
    }

    spcol_map = new int[nxsparse];
    for (int i = 0; i < nxsparse; i++) {
      spcol_map[i] = i;
    }

    dncol_map = new int[nxdense];
    for (int i = 0; i < nxdense; i++) {
      dncol_map[i] = i + nxsparse;
    }

    ierr = ISCreateGeneral(PETSC_COMM_WORLD, nrow, permutation_map, PETSC_COPY_VALUES, &permutation_map_IS);CHKERRQ(ierr);
    ierr = ISCreateGeneral(PETSC_COMM_WORLD, nxsparse, spcol_map, PETSC_COPY_VALUES, &spcol_map_IS);CHKERRQ(ierr);
    ierr = ISCreateGeneral(PETSC_COMM_WORLD, nxdense, dncol_map, PETSC_COPY_VALUES, &dncol_map_IS);CHKERRQ(ierr);

    ierr = MatPermute(Hessref, permutation_map_IS, permutation_map_IS, &Hessref_spdn);CHKERRQ(ierr);

    ierr = MatCreateSubMatrix(Hessref_spdn, spcol_map_IS, spcol_map_IS, MAT_INITIAL_MATRIX, &Hessref_sparse);CHKERRQ(ierr);
    ierr = MatCreateSubMatrix(Hessref_spdn, dncol_map_IS, dncol_map_IS, MAT_INITIAL_MATRIX, &Hessref_dense);CHKERRQ(ierr);

    // Test sparse Hessian
    int *iRow, *jCol,*iRow_dev,*jCol_dev;
    double *values,*values_dev;
    int nnz = nxsparse;

    iRow = static_cast<int*>(h_allocator.allocate(nnz*sizeof(int)));
    jCol = static_cast<int*>(h_allocator.allocate(nnz*sizeof(int)));
    values = static_cast<double*>(h_allocator.allocate(nnz*sizeof(double)));
#ifdef EXAGO_ENABLE_GPU
    umpire::Allocator d_allocator = resmgr.getAllocator("DEVICE");
    iRow_dev = static_cast<int*>(d_allocator.allocate(nnz*sizeof(int)));
    jCol_dev = static_cast<int*>(d_allocator.allocate(nnz*sizeof(int)));
    values_dev = static_cast<double*>(d_allocator.allocate(nnz*sizeof(double)));
#else
    iRow_dev = iRow;
    jCol_dev = jCol;
    values_dev = values;
#endif

    opflow->obj_factor = obj_factor;
    ierr = (*opflow->modelops.computesparsehessianhiop)(opflow, x_ref_dev, iRow_dev, jCol_dev, values_dev);CHKERRQ(ierr);

    // Copy back from the device
    resmgr.copy(iRow,iRow_dev);
    resmgr.copy(jCol,jCol_dev);
    resmgr.copy(values,values_dev);

    fail += verifyAnswer(Hessref_sparse, nnz, iRow, jCol, values);

    // Test dense Hessian
    ierr = (*opflow->modelops.computedensehessianhiop)(opflow, x_ref_dev, lambda_ref_dev, hess_dense_dev);CHKERRQ(ierr);

    resmgr.copy(hess_dense,hess_dense_dev);

    fail += verifyAnswer(Hessref_dense, hess_dense);

    delete[] permutation_map;
    delete[] spcol_map;
    delete[] dncol_map;
    ierr = MatDestroy(&Hessref_spdn);CHKERRQ(ierr);
    ierr = MatDestroy(&Hessref_sparse);CHKERRQ(ierr);
    ierr = MatDestroy(&Hessref_dense);CHKERRQ(ierr);
    ierr = ISDestroy(&permutation_map_IS);CHKERRQ(ierr);
    ierr = ISDestroy(&spcol_map_IS);CHKERRQ(ierr);
    ierr = ISDestroy(&dncol_map_IS);CHKERRQ(ierr);

    h_allocator.deallocate(iRow);
    h_allocator.deallocate(jCol);
    h_allocator.deallocate(values);
#ifdef EXAGO_ENABLE_GPU
    d_allocator.deallocate(iRow_dev);
    d_allocator.deallocate(jCol_dev);
    d_allocator.deallocate(values_dev);
#endif
    cleanup(fail, opflow);
  }
#endif

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
   * @brief Verifies that two arrays are equivilant within tolerance.
   *
   */
  virtual int verifyAnswer(const double* avals, const double* bvals, int size,const RealType& tol=eps) const
  {
    LocalOrdinalType   fail = 0;

    for (int i=0; i<size; i++)
    {
      if (!isEqual(avals[i], bvals[i], tol)) {
	printf("avals[%d] = %18.16f not equal to bvals[%d] = %18.16f. Exceeds set tolerance = %g\n",i,avals[i],i,bvals[i],tol);
        fail++;
      }
    }

    return(fail);
    //    return reduceReturn(fail, (PetscObject)a);
  }


  /**
   * @brief Verifies that two PETSc vectors are equivilant within tolerance.
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

  /**
   * @brief Compare two matrices. This overload verifies against a dense matrix
   * with a contiguous layout.
   */
  virtual int verifyAnswer(Mat a, double* b, const RealType& tol=eps) const
  {
    int              ncols;
    int              fail = 0;
    const int        *cols;
    const double     *vals;
    int              count = 0;
    PetscInt         nrow, ncol;
    PetscErrorCode   ierr;
    auto idx = [&ncol] (double* mat, int r, int c)
      { return mat[(r*ncol)+c]; };

    ierr = MatGetSize(a, &nrow, &ncol);CHKERRQ(ierr);
    for(int i = 0; i < nrow; i++) {
      ierr = MatGetRow(a, i, &ncols, &cols, &vals);CHKERRQ(ierr);
      for(int j = 0; j < ncols; j++) {
        if(!isEqual(vals[j], idx(b, i, cols[j]), tol)) {
          std::cout << "Failed for index (" << i << ", " << cols[j] << ") : " << vals[j] << " != " << idx(b, i, cols[j]) << std::endl;
          fail++;
        }
        count++;
      }
      ierr = MatRestoreRow(a, i, &ncols, &cols, &vals);CHKERRQ(ierr);
    }
    return fail;
  }

  /**
   * @brief Compare two matrices. Matrix a is in PETSc sparse matrix format and b is a dense matrix.
   * @note @abhyshr will we still need this overload?
   */
  virtual int verifyAnswer(Mat a, double** b, const RealType& tol=eps) const
  {
    int              ncols;
    int              fail = 0;
    const int        *cols;
    const double     *vals;
    int              count = 0;
    PetscInt         nrow, ncol;
    PetscErrorCode   ierr;

    ierr = MatGetSize(a, &nrow, &ncol);CHKERRQ(ierr);
    for(int i = 0; i < nrow; i++) {
      ierr = MatGetRow(a, i, &ncols, &cols, &vals);CHKERRQ(ierr);
      for(int j = 0; j < ncols; j++) {
	if(!isEqual(vals[j], b[i][cols[j]], tol)) {
	  std::cout << "Failed for index (" << i << ", " << cols[j] << ") : " << vals[j] << " != " << b[i][cols[j]] << std::endl;
	  fail++;
	}
	count++;
      }
      ierr = MatRestoreRow(a, i, &ncols, &cols, &vals);CHKERRQ(ierr);
    }
    return fail;
  }
  
  /**
   * @brief Compares two matrices - a is in PETSc format and the other matrix is described by sparse triplet iRow, jCol, values
   *
   */
  virtual int verifyAnswer(Mat a, int nnz, int *iRow, int *jCol, double *values, const RealType &tol = eps) const
  {
    double         val;
    int            fail = 0;
    PetscErrorCode ierr;

    for(int i=0; i < nnz; i++) {
      ierr = MatGetValues(a,1,iRow+i,1,jCol+i,&val);CHKERRQ(ierr);

      if(!isEqual(val, values[i], tol)) {
	std::cout << "Failed for index (" << iRow[i] << ", " << jCol[i] << ") : " << val << " != " << values[i] << std::endl;
	fail++;
      }
    }
    return fail;
  }
}; // class TestOpflow : public TestBase
  
}} // namespace exago::tests
