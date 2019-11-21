#if defined(SCOPFLOW_HAVE_IPOPT)

#include <private/scopflowimpl.h>
#include "scopflow-ipopt.h"

static int CCMatrixToMatrixMarketValuesOnly(CCMatrix ccmatrix,PetscInt nz,PetscScalar *values)
{
  PetscErrorCode ierr;

  ierr = PetscMemcpy(values,ccmatrix->values,nz*sizeof(PetscScalar));

  return 0;
}

static int CCMatrixToMatrixMarketLocationsOnly(CCMatrix ccmatrix,PetscInt ncol,PetscInt *iRow,PetscInt *jCol,PetscInt roffset,PetscInt coffset,PetscInt nval)
{
  PetscInt *rowidx;
  PetscInt *colptr;
  PetscInt j,k,ctr=0;
  
  rowidx = ccmatrix->rowidx;
  colptr = ccmatrix->colptr;

  /* Copy from compressed column to (row,col,val) format */
  for(j=0; j < ncol; j++) {
    for(k=colptr[j]; k < colptr[j+1]; k++) {
      iRow[ctr] = rowidx[k] + roffset;
      jCol[ctr] = j + coffset;
      ctr++;
    }
  }
  if(ctr != nval) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,"Incorrect number of entries ctr = %d given = %d\n",ctr,nval);

  return 0;
}

/* IPOPT callback functions */
Bool eval_scopflow_f(PetscInt n, PetscScalar* x, Bool new_x,
            PetscScalar* obj_value, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  SCOPFLOW  scopflow=(SCOPFLOW)user_data;

  *obj_value = 0.0;
  ierr = VecPlaceArray(scopflow->X,x);CHKERRQ(ierr);
  ierr = (*scopflow->formops.computeobjective)(scopflow,scopflow->X,obj_value);CHKERRQ(ierr);
  ierr = VecResetArray(scopflow->X);CHKERRQ(ierr);
				
  return TRUE;
}

Bool eval_scopflow_grad_f(PetscInt n, PetscScalar* x, Bool new_x,
                 PetscScalar* grad_f, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  SCOPFLOW  scopflow=(SCOPFLOW)user_data;

  ierr = VecPlaceArray(scopflow->X,x);CHKERRQ(ierr);
  ierr = VecPlaceArray(scopflow->gradobj,grad_f);CHKERRQ(ierr);
  ierr = VecSet(scopflow->gradobj,0.0);CHKERRQ(ierr);
  ierr = (*scopflow->formops.computegradient)(scopflow,scopflow->X,scopflow->gradobj);CHKERRQ(ierr);
  ierr = VecResetArray(scopflow->X);CHKERRQ(ierr);
  ierr = VecResetArray(scopflow->gradobj);CHKERRQ(ierr);

  return TRUE;
}

Bool eval_scopflow_g(PetscInt n, PetscScalar* x, Bool new_x,
            PetscInt m, PetscScalar* g, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  SCOPFLOW         scopflow=(SCOPFLOW)user_data;

  ierr = VecPlaceArray(scopflow->X,x);CHKERRQ(ierr);

  /* Equality constraints */
  ierr = VecPlaceArray(scopflow->Ge,g);CHKERRQ(ierr);
  ierr = (*scopflow->formops.computeequalityconstraints)(scopflow,scopflow->X,scopflow->Ge);CHKERRQ(ierr);
  ierr = VecResetArray(scopflow->Ge);CHKERRQ(ierr);

  if(scopflow->Nconineq) {
    /* Inequality constraints */
    ierr = VecPlaceArray(scopflow->Gi,g+scopflow->nconeq);CHKERRQ(ierr);
    ierr = (*scopflow->formops.computeinequalityconstraints)(scopflow,scopflow->X,scopflow->Gi);CHKERRQ(ierr);
    ierr = VecResetArray(scopflow->Gi);CHKERRQ(ierr);
  }

  ierr = VecResetArray(scopflow->X);CHKERRQ(ierr);

  return TRUE;
}

Bool eval_scopflow_jac_g(PetscInt n, PetscScalar *x, Bool new_x,
                PetscInt m, PetscInt nele_jac,
                PetscInt *iRow, PetscInt *jCol, PetscScalar *values,
                UserDataPtr user_data)
{

  PetscErrorCode ierr;
  SCOPFLOW         scopflow=(SCOPFLOW)user_data;
  SCOPFLOWSolver_IPOPT ipopt=(SCOPFLOWSolver_IPOPT)scopflow->solver;
  PetscInt       *iRowstart = iRow,*jColstart=jCol;
  PetscInt       roffset,coffset;
  Mat_SeqAIJ     *aij;
  PetscInt       nrow,ncol;


  if(values == NULL) {
    /* Set locations only */

    roffset = 0;
    coffset = 0;

    /* Equality constrained Jacobian */
    ierr = (*scopflow->formops.computeequalityconstraintjacobian)(scopflow,scopflow->X,scopflow->Jac_Ge);CHKERRQ(ierr);
    /* Transpose the matrix to convert it to column compressed sparse format */
    ierr = MatTranspose(scopflow->Jac_Ge,MAT_REUSE_MATRIX,&ipopt->Jac_GeT);CHKERRQ(ierr);

    ierr = MatGetSize(ipopt->Jac_GeT,&nrow,&ncol);CHKERRQ(ierr);
    aij = (Mat_SeqAIJ*)ipopt->Jac_GeT->data;
    ierr = PetscMemcpy(ipopt->jac_ge->rowidx,aij->j,aij->nz*sizeof(PetscInt));CHKERRQ(ierr);
    ierr = PetscMemcpy(ipopt->jac_ge->colptr,aij->i,(nrow+1)*sizeof(PetscInt));CHKERRQ(ierr);
    ierr = PetscMemcpy(ipopt->jac_ge->values,aij->a,aij->nz*sizeof(PetscScalar));CHKERRQ(ierr);
    CCMatrixToMatrixMarketLocationsOnly(ipopt->jac_ge,scopflow->nx,iRowstart,jColstart,roffset,coffset,ipopt->nnz_jac_ge);

    /* Increment iRow,jCol pointers */
    iRowstart += ipopt->nnz_jac_ge;
    jColstart += ipopt->nnz_jac_ge;

    if(scopflow->Nconineq) {
      /* Inequality constrained Jacobian */
      roffset = scopflow->nconeq;

      ierr = (*scopflow->formops.computeinequalityconstraintjacobian)(scopflow,scopflow->X,scopflow->Jac_Gi);CHKERRQ(ierr);
      /* Transpose the matrix to convert it to column compressed sparse format */
      ierr = MatTranspose(scopflow->Jac_Gi,MAT_REUSE_MATRIX,&ipopt->Jac_GiT);CHKERRQ(ierr);

      ierr = MatGetSize(ipopt->Jac_GiT,&nrow,&ncol);CHKERRQ(ierr);
      aij = (Mat_SeqAIJ*)ipopt->Jac_GiT->data;
      ierr = PetscMemcpy(ipopt->jac_gi->rowidx,aij->j,aij->nz*sizeof(PetscInt));CHKERRQ(ierr);
      ierr = PetscMemcpy(ipopt->jac_gi->colptr,aij->i,(nrow+1)*sizeof(PetscInt));CHKERRQ(ierr);
      ierr = PetscMemcpy(ipopt->jac_gi->values,aij->a,aij->nz*sizeof(PetscScalar));CHKERRQ(ierr);
      
      CCMatrixToMatrixMarketLocationsOnly(ipopt->jac_gi,scopflow->nx,iRowstart,jColstart,roffset,coffset,ipopt->nnz_jac_gi);

    }
  } else {
    ierr = VecPlaceArray(scopflow->X,x);CHKERRQ(ierr);
    /* Compute equality constraint jacobian */
    ierr = (*scopflow->formops.computeequalityconstraintjacobian)(scopflow,scopflow->X,scopflow->Jac_Ge);CHKERRQ(ierr);

    ierr = MatTranspose(scopflow->Jac_Ge,MAT_REUSE_MATRIX,&ipopt->Jac_GeT);CHKERRQ(ierr);
    ierr = MatGetSize(ipopt->Jac_GeT,&nrow,&ncol);CHKERRQ(ierr);
    aij = (Mat_SeqAIJ*)ipopt->Jac_GeT->data;

    ierr = PetscMemcpy(ipopt->jac_ge->values,aij->a,aij->nz*sizeof(PetscScalar));CHKERRQ(ierr);

    CCMatrixToMatrixMarketValuesOnly(ipopt->jac_ge,ipopt->nnz_jac_ge,values);

    if(scopflow->Nconineq) {
      /* Compute inequality constraint jacobian */
      ierr = (*scopflow->formops.computeinequalityconstraintjacobian)(scopflow,scopflow->X,scopflow->Jac_Gi);CHKERRQ(ierr);

      ierr = MatTranspose(scopflow->Jac_Gi,MAT_REUSE_MATRIX,&ipopt->Jac_GiT);CHKERRQ(ierr);
      aij = (Mat_SeqAIJ*)ipopt->Jac_GiT->data;

      ierr = PetscMemcpy(ipopt->jac_gi->values,aij->a,aij->nz*sizeof(PetscScalar));CHKERRQ(ierr);
      
      CCMatrixToMatrixMarketValuesOnly(ipopt->jac_gi,ipopt->nnz_jac_gi,values+ipopt->nnz_jac_ge);
    }
    ierr = VecResetArray(scopflow->X);CHKERRQ(ierr);
  }

  return TRUE;
}

Bool eval_scopflow_h(PetscInt n, PetscScalar *x, Bool new_x, PetscScalar obj_factor,
            PetscInt m, PetscScalar *lambda, Bool new_lambda,
            PetscInt nele_hess, PetscInt *iRow, PetscInt *jCol,
            PetscScalar *values, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  PetscInt       nrow;
  Mat_SeqSBAIJ   *sbaij;
  SCOPFLOW         scopflow=(SCOPFLOW)user_data;
  SCOPFLOWSolver_IPOPT ipopt=(SCOPFLOWSolver_IPOPT)scopflow->solver;

  scopflow->obj_factor = obj_factor;

  if(values == NULL) {
    ierr = (*scopflow->formops.computehessian)(scopflow,scopflow->X,scopflow->Lambdae,scopflow->Lambdai,scopflow->Hes);CHKERRQ(ierr);
    ierr = MatGetSize(scopflow->Hes,&nrow,&nrow);CHKERRQ(ierr);
    ierr = MatConvert(scopflow->Hes,MATSEQSBAIJ,MAT_REUSE_MATRIX,&ipopt->Hes_sbaij);CHKERRQ(ierr);
    sbaij = (Mat_SeqSBAIJ*)ipopt->Hes_sbaij->data;
    ierr = PetscMemcpy(ipopt->hes->rowidx,sbaij->j,sbaij->nz*sizeof(PetscInt));CHKERRQ(ierr);
    ierr = PetscMemcpy(ipopt->hes->colptr,sbaij->i,(nrow+1)*sizeof(PetscInt));CHKERRQ(ierr);
    ierr = PetscMemcpy(ipopt->hes->values,sbaij->a,sbaij->nz*sizeof(PetscScalar));CHKERRQ(ierr);
    CCMatrixToMatrixMarketLocationsOnly(ipopt->hes,scopflow->nx,iRow,jCol,0,0,ipopt->nnz_hes);
  } else {
    ierr = VecPlaceArray(scopflow->X,x);CHKERRQ(ierr);
    ierr = VecPlaceArray(scopflow->Lambdae,lambda);CHKERRQ(ierr);
    if(scopflow->Nconineq) {
      ierr = VecPlaceArray(scopflow->Lambdai,lambda+scopflow->nconeq);CHKERRQ(ierr);
    }

    /* Compute non-zeros for Hessian */
    ierr = (*scopflow->formops.computehessian)(scopflow,scopflow->X,scopflow->Lambdae,scopflow->Lambdai,scopflow->Hes);CHKERRQ(ierr);
    ierr = MatConvert(scopflow->Hes,MATSEQSBAIJ,MAT_REUSE_MATRIX,&ipopt->Hes_sbaij);CHKERRQ(ierr);
    /* Since the Hessian is symmetric, we don't need to convert it to column compressed sparse format */
    sbaij = (Mat_SeqSBAIJ*)ipopt->Hes_sbaij->data;
    ipopt->nnz_hes = sbaij->nz;
    ierr = PetscMemcpy(ipopt->hes->values,sbaij->a,sbaij->nz*sizeof(PetscScalar));CHKERRQ(ierr);
    CCMatrixToMatrixMarketValuesOnly(ipopt->hes,ipopt->nnz_hes,values);

    ierr = VecResetArray(scopflow->X);CHKERRQ(ierr);
    ierr = VecResetArray(scopflow->Lambdae);CHKERRQ(ierr);
    if(scopflow->Nconineq) {
      ierr = VecResetArray(scopflow->Lambdai);CHKERRQ(ierr);
    }
  }

  return 1;
}

PetscErrorCode SCOPFLOWSolverSolve_IPOPT(SCOPFLOW scopflow)
{
  PetscErrorCode     ierr;
  SCOPFLOWSolver_IPOPT ipopt=scopflow->solver;
  Mat_SeqAIJ         *aij;
  Mat_SeqSBAIJ       *sbaij;
  PetscScalar        *x,*xl,*xu,*gl,*gu;

  PetscFunctionBegin;

  /* Compute nonzeros for the Jacobian */
  /* Equality constraint Jacobian */
  ierr = (*scopflow->formops.computeequalityconstraintjacobian)(scopflow,scopflow->X,scopflow->Jac_Ge);CHKERRQ(ierr);
  /* Transpose the matrix to convert it to column compressed sparse format */
  ierr = MatTranspose(scopflow->Jac_Ge,MAT_INITIAL_MATRIX,&ipopt->Jac_GeT);CHKERRQ(ierr);
  ierr = MatSetOption(ipopt->Jac_GeT,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);
  aij = (Mat_SeqAIJ*)ipopt->Jac_GeT->data;
  ipopt->nnz_jac_ge = aij->nz;

  /* Create ccmatrix object for equality constrained Jacobian */
  ierr = PetscCalloc1(1,&ipopt->jac_ge);CHKERRQ(ierr);
  ierr = PetscMalloc1(scopflow->nx+1,&ipopt->jac_ge->colptr);CHKERRQ(ierr);
  ierr = PetscMalloc1(ipopt->nnz_jac_ge,&ipopt->jac_ge->rowidx);CHKERRQ(ierr);
  ierr = PetscMalloc1(ipopt->nnz_jac_ge,&ipopt->jac_ge->values);CHKERRQ(ierr);

  ipopt->nnz_jac_gi = 0;
  if(scopflow->Nconineq) {
    ierr = (*scopflow->formops.computeinequalityconstraintjacobian)(scopflow,scopflow->X,scopflow->Jac_Gi);CHKERRQ(ierr);
    /* Transpose the matrix to convert it to column compressed sparse format */
    ierr = MatTranspose(scopflow->Jac_Gi,MAT_INITIAL_MATRIX,&ipopt->Jac_GiT);CHKERRQ(ierr);
    ierr = MatSetOption(ipopt->Jac_GiT,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);
    aij = (Mat_SeqAIJ*)ipopt->Jac_GiT->data;
    ipopt->nnz_jac_gi = aij->nz;

    /* Create ccmatrix object for inequality constrained Jacobian */
    ierr = PetscCalloc1(1,&ipopt->jac_gi);CHKERRQ(ierr);
    ierr = PetscMalloc1(scopflow->nx+1,&ipopt->jac_gi->colptr);CHKERRQ(ierr);
    ierr = PetscMalloc1(ipopt->nnz_jac_gi,&ipopt->jac_gi->rowidx);CHKERRQ(ierr);
    ierr = PetscMalloc1(ipopt->nnz_jac_gi,&ipopt->jac_gi->values);CHKERRQ(ierr);
  }   
  ipopt->nnz_jac_g = ipopt->nnz_jac_ge + ipopt->nnz_jac_gi;

  /* Compute non-zeros for Hessian */
  ierr = (*scopflow->formops.computehessian)(scopflow,scopflow->X,scopflow->Lambdae,scopflow->Lambdai,scopflow->Hes);CHKERRQ(ierr);
  /* Convert matrix to symmetric sbaij format needed for the IPOPT solver */
  ierr = MatSetOption(scopflow->Hes,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MatConvert(scopflow->Hes,MATSEQSBAIJ,MAT_INITIAL_MATRIX,&ipopt->Hes_sbaij);CHKERRQ(ierr);
  /* Since the Hessian is symmetric, we don't need to convert it to column compressed sparse format */
  sbaij = (Mat_SeqSBAIJ*)ipopt->Hes_sbaij->data;
  ipopt->nnz_hes = sbaij->nz;

  /* Create ccmatrix object for hessian */
  ierr = PetscCalloc1(1,&ipopt->hes);CHKERRQ(ierr);
  ierr = PetscMalloc1(scopflow->nx+1,&ipopt->hes->colptr);CHKERRQ(ierr);
  ierr = PetscMalloc1(ipopt->nnz_hes,&ipopt->hes->rowidx);CHKERRQ(ierr);
  ierr = PetscMalloc1(ipopt->nnz_hes,&ipopt->hes->values);CHKERRQ(ierr);

  /* Create IPOPT solver instance */
  ierr = VecGetArray(scopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Gu,&gu);CHKERRQ(ierr);

  /* Create IPOPT problem */
  ipopt->nlp = CreateIpoptProblem(scopflow->nx,xl,xu,scopflow->ncon,gl,gu,ipopt->nnz_jac_g,ipopt->nnz_hes,0,&eval_scopflow_f,
				   &eval_scopflow_g,&eval_scopflow_grad_f,
				   &eval_scopflow_jac_g,&eval_scopflow_h);

  ierr = VecRestoreArray(scopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Gu,&gu);CHKERRQ(ierr);

  ierr = VecGetArray(scopflow->X,&x);CHKERRQ(ierr);
  /* Solve */
  ipopt->solve_status = IpoptSolve(ipopt->nlp,x,NULL,&scopflow->obj,NULL,NULL,NULL,scopflow);

  ierr = VecRestoreArray(scopflow->X,&x);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverDestroy_IPOPT(SCOPFLOW scopflow)
{
  PetscErrorCode     ierr;
  SCOPFLOWSolver_IPOPT ipopt=scopflow->solver;

  PetscFunctionBegin;

  if(ipopt->nlp) {
    FreeIpoptProblem(ipopt->nlp);
  }
  ierr = MatDestroy(&ipopt->Jac_GeT);CHKERRQ(ierr);
  ierr = PetscFree(ipopt->jac_ge->colptr);CHKERRQ(ierr);
  ierr = PetscFree(ipopt->jac_ge->rowidx);CHKERRQ(ierr);
  ierr = PetscFree(ipopt->jac_ge->values);CHKERRQ(ierr);
  ierr = PetscFree(ipopt->jac_ge);CHKERRQ(ierr);

  if(scopflow->Nconineq) {
    ierr = MatDestroy(&ipopt->Jac_GiT);CHKERRQ(ierr);
    ierr = PetscFree(ipopt->jac_gi->colptr);CHKERRQ(ierr);
    ierr = PetscFree(ipopt->jac_gi->rowidx);CHKERRQ(ierr);
    ierr = PetscFree(ipopt->jac_gi->values);CHKERRQ(ierr);
    ierr = PetscFree(ipopt->jac_gi);CHKERRQ(ierr);
  }

  ierr = MatDestroy(&(ipopt->Hes_sbaij));CHKERRQ(ierr);
  ierr = PetscFree(ipopt->hes->colptr);CHKERRQ(ierr);
  ierr = PetscFree(ipopt->hes->rowidx);CHKERRQ(ierr);
  ierr = PetscFree(ipopt->hes->values);CHKERRQ(ierr);
  ierr = PetscFree(ipopt->hes);CHKERRQ(ierr);

  ierr = PetscFree(ipopt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverSetUp_IPOPT(SCOPFLOW scopflow)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverCreate_IPOPT(SCOPFLOW scopflow)
{
  PetscErrorCode ierr;
  SCOPFLOWSolver_IPOPT ipopt;
  
  PetscFunctionBegin;

  if(scopflow->comm->size > 1) SETERRQ1(PETSC_ERR_SUP,PETSC_COMM_SELF,"nrank = %d, IPOPT solver does not support execution in parallel\n",scopflow->comm->size); 
  ierr = PetscCalloc1(1,&ipopt);CHKERRQ(ierr);

  ipopt->nlp = NULL;
  ipopt->nnz_jac_g = 0;
  ipopt->nnz_hes = 0;
  scopflow->solver = ipopt;

  scopflow->solverops.setup = SCOPFLOWSolverSetUp_IPOPT;
  scopflow->solverops.solve = SCOPFLOWSolverSolve_IPOPT;
  scopflow->solverops.destroy = SCOPFLOWSolverDestroy_IPOPT;

  PetscFunctionReturn(0);
}

#endif
