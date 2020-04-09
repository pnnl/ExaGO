#include <scopflow_config.h>
#if defined(SCOPFLOW_HAVE_IPOPT)

#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>
#include "scopflow-ipopt.h"

static int CCMatrixToMatrixMarketValuesOnly(CCMatrix ccmatrix,PetscInt nz,PetscScalar *values)
{
  PetscErrorCode ierr;

  ierr = PetscMemcpy(values,ccmatrix->values,nz*sizeof(PetscScalar));CHKERRQ(ierr);

  return ierr;
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
  SCOPFLOWSolver_IPOPT scopflowipopt = (SCOPFLOWSolver_IPOPT)scopflow->solver;
  OPFLOW    opflow;
  PetscInt  i;
  PetscScalar *xi;
  PetscScalar opflowobj;
  PetscInt    k;

  *obj_value = 0.0;

  k = scopflow->first_stage_gen_cost_only?1:scopflow->Ns;

  for(i=0; i < k; i++) {
    opflowobj = 0.0;
    xi = x + scopflowipopt->xstarti[i];
    opflow = scopflow->opflows[i];
    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
    ierr = (*opflow->formops.computeobjective)(opflow,opflow->X,&opflowobj);CHKERRQ(ierr);
    *obj_value += opflowobj;
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  }
				
  return TRUE;
}

Bool eval_scopflow_grad_f(PetscInt n, PetscScalar* x, Bool new_x,
                 PetscScalar* grad_f, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  SCOPFLOW  scopflow=(SCOPFLOW)user_data;
  SCOPFLOWSolver_IPOPT scopflowipopt = (SCOPFLOWSolver_IPOPT)scopflow->solver;
  OPFLOW    opflow;
  PetscInt  i;
  PetscScalar *xi,*gradi;

  for(i=0; i < scopflow->Ns; i++) {
    opflow = scopflow->opflows[i];
    xi = x + scopflowipopt->xstarti[i];
    gradi = grad_f + scopflowipopt->xstarti[i];


    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->gradobj,gradi);CHKERRQ(ierr);
    ierr = VecSet(opflow->gradobj,0.0);CHKERRQ(ierr);

    if(!scopflow->first_stage_gen_cost_only) {
      ierr = (*opflow->formops.computegradient)(opflow,opflow->X,opflow->gradobj);CHKERRQ(ierr);
    } else {
      if(i == 0) {
	ierr = (*opflow->formops.computegradient)(opflow,opflow->X,opflow->gradobj);CHKERRQ(ierr);
      }
    }
    
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->gradobj);CHKERRQ(ierr);
  }
  return TRUE;
}

Bool eval_scopflow_g(PetscInt n, PetscScalar* x, Bool new_x,
            PetscInt m, PetscScalar* g, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  SCOPFLOW  scopflow=(SCOPFLOW)user_data;
  SCOPFLOWSolver_IPOPT scopflowipopt = (SCOPFLOWSolver_IPOPT)scopflow->solver;
  OPFLOW    opflow0,opflow;
  PetscInt  i,j,k,loc,loc0,ctr;
  PetscScalar *x0,*xi,*gi;
  PS        ps,ps0;
  PSBUS     bus,bus0;
  PSGEN     gen,gen0;

  x0 = x;

  opflow0 = scopflow->opflows[0];
  for(i=0; i < scopflow->Ns; i++) {
    xi   = x + scopflowipopt->xstarti[i];
    gi   = g + scopflowipopt->gstarti[i];

    opflow = scopflow->opflows[i];

    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);

    /* Equality constraints */
    ierr = VecPlaceArray(opflow->Ge,gi);CHKERRQ(ierr);
    ierr = (*opflow->formops.computeequalityconstraints)(opflow,opflow->X,opflow->Ge);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Ge);CHKERRQ(ierr);
    gi = gi + opflow->nconeq;
      
    if(opflow->Nconineq) {
      /* Inequality constraints */
      ierr = VecPlaceArray(opflow->Gi,gi);CHKERRQ(ierr);
      ierr = (*opflow->formops.computeinequalityconstraints)(opflow,opflow->X,opflow->Gi);CHKERRQ(ierr);
      ierr = VecResetArray(opflow->Gi);CHKERRQ(ierr);
      gi = gi + opflow->nconineq;
    }

    if(scopflow->nconineqcoup[i]) {
      ctr = 0;
      ps = opflow->ps;
      ps0 = opflow0->ps;
      for(j=0; j < ps->nbus; j++) {
	bus = &ps->bus[j];
	bus0 = &ps0->bus[j];
	ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
	ierr = PSBUSGetVariableLocation(bus0,&loc0);CHKERRQ(ierr);
	
	for(k=0; k < bus->ngen; k++) {
	  ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	  ierr = PSBUSGetGen(bus0,k,&gen0);CHKERRQ(ierr);

	  if(!gen->status) {
	    if(gen0->status) loc0 += 2;
	    continue;
	  } else {
	    loc += 2;
	    if(!gen0->status) continue;
	    loc0 += 2;
	  }

	  gi[ctr] = xi[loc] - x0[loc0]; /* PGi - PG0 */
	  ctr++;
	}
      }
    }
    
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  }

  return TRUE;
}

Bool eval_scopflow_jac_g(PetscInt n, PetscScalar *x, Bool new_x,
                PetscInt m, PetscInt nele_jac,
                PetscInt *iRow, PetscInt *jCol, PetscScalar *values,
                UserDataPtr user_data)
{

  PetscErrorCode ierr;
  SCOPFLOW         scopflow=(SCOPFLOW)user_data;
  SCOPFLOWSolver_IPOPT scopflowipopt=(SCOPFLOWSolver_IPOPT)scopflow->solver;
  OPFLOW         opflow,opflow0;
  OPFLOWSolver_IPOPT opflowipopt;
  PetscInt       *iRowstart = iRow,*jColstart=jCol;
  PetscInt       roffset,coffset;
  Mat_SeqAIJ     *aij;
  PetscInt       nrow,ncol;
  PetscScalar    *xi,*xarr,*valuesi;
  PetscInt       i,j,k,loc,loc0,x0loc,xiloc;
  PS             ps,ps0;
  PSBUS          bus,bus0;
  PSGEN          gen,gen0;

  if(values == NULL) {
    /* Set locations only */
    ierr = VecGetArray(scopflow->X,&xarr);CHKERRQ(ierr);
    opflow0 = scopflow->opflows[0];
    for(i=0; i < scopflow->Ns; i++) {
      opflow = scopflow->opflows[i];
      opflowipopt = (OPFLOWSolver_IPOPT)opflow->solver;

      roffset = scopflowipopt->gstarti[i];
      coffset = scopflowipopt->xstarti[i];

      xi = xarr + scopflowipopt->xstarti[i];
      ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);

      ierr = (*opflow->formops.computeequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Ge);CHKERRQ(ierr);
      /* Transpose the matrix to convert it to column compressed sparse format */
      ierr = MatTranspose(opflow->Jac_Ge,MAT_REUSE_MATRIX,&opflowipopt->Jac_GeT);CHKERRQ(ierr);

      ierr = MatGetSize(opflowipopt->Jac_GeT,&nrow,&ncol);CHKERRQ(ierr);
      aij = (Mat_SeqAIJ*)opflowipopt->Jac_GeT->data;
      ierr = PetscMemcpy(opflowipopt->jac_ge->rowidx,aij->j,aij->nz*sizeof(PetscInt));CHKERRQ(ierr);
      ierr = PetscMemcpy(opflowipopt->jac_ge->colptr,aij->i,(nrow+1)*sizeof(PetscInt));CHKERRQ(ierr);
      ierr = PetscMemcpy(opflowipopt->jac_ge->values,aij->a,aij->nz*sizeof(PetscScalar));CHKERRQ(ierr);
      CCMatrixToMatrixMarketLocationsOnly(opflowipopt->jac_ge,opflow->nx,iRowstart,jColstart,roffset,coffset,opflowipopt->nnz_jac_ge);

      /* Increment iRow,jCol pointers and roffset */
      iRowstart += opflowipopt->nnz_jac_ge;
      jColstart += opflowipopt->nnz_jac_ge;
      roffset += opflow->nconeq;

      if(opflow->Nconineq) {
	/* Inequality constrained Jacobian */
	ierr = (*opflow->formops.computeinequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Gi);CHKERRQ(ierr);
	/* Transpose the matrix to convert it to column compressed sparse format */
	ierr = MatTranspose(opflow->Jac_Gi,MAT_REUSE_MATRIX,&opflowipopt->Jac_GiT);CHKERRQ(ierr);

	ierr = MatGetSize(opflowipopt->Jac_GiT,&nrow,&ncol);CHKERRQ(ierr);
	aij = (Mat_SeqAIJ*)opflowipopt->Jac_GiT->data;
	ierr = PetscMemcpy(opflowipopt->jac_gi->rowidx,aij->j,aij->nz*sizeof(PetscInt));CHKERRQ(ierr);
	ierr = PetscMemcpy(opflowipopt->jac_gi->colptr,aij->i,(nrow+1)*sizeof(PetscInt));CHKERRQ(ierr);
	ierr = PetscMemcpy(opflowipopt->jac_gi->values,aij->a,aij->nz*sizeof(PetscScalar));CHKERRQ(ierr);
      
	CCMatrixToMatrixMarketLocationsOnly(opflowipopt->jac_gi,opflow->nx,iRowstart,jColstart,roffset,coffset,opflowipopt->nnz_jac_gi);

	iRowstart += opflowipopt->nnz_jac_gi;
	jColstart += opflowipopt->nnz_jac_gi;
	roffset += opflow->nconineq;
      }

      if(scopflow->nconineqcoup[i]) {
	ps = opflow->ps;
	ps0 = opflow0->ps;
	for(j=0; j < ps->nbus; j++) {
	  bus = &ps->bus[j];
	  bus0 = &ps0->bus[j];
	  ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
	  ierr = PSBUSGetVariableLocation(bus0,&loc0);CHKERRQ(ierr);
	  for(k=0; k < bus->ngen; k++) {
	    ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	    ierr = PSBUSGetGen(bus0,k,&gen0);CHKERRQ(ierr);

	    if(!gen->status) {
	      if(gen0->status) loc0 += 2;
	      continue;
	    } else {
	      loc += 2;
	      if(!gen0->status) continue;
	      loc0 += 2;
	    }

	    x0loc = scopflowipopt->xstarti[0] + loc0;
	    xiloc = scopflowipopt->xstarti[i] + loc;
	    iRowstart[0] = roffset;
	    jColstart[0] = x0loc;
	    iRowstart[1] = roffset;
	    jColstart[1] = xiloc;
	    iRowstart += 2;
	    jColstart += 2;
	    roffset += 1;
	  }
	}
      }

      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
    }
    ierr = VecRestoreArray(scopflow->X,&xarr);CHKERRQ(ierr);
  } else {

    opflow0 = scopflow->opflows[0];
    valuesi = values;
    for(i=0; i < scopflow->Ns; i++) {
      opflow = scopflow->opflows[i];
      opflowipopt = (OPFLOWSolver_IPOPT)opflow->solver;

      xi = x + scopflowipopt->xstarti[i];
      ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
      /* Compute equality constraint jacobian */
      ierr = (*opflow->formops.computeequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Ge);CHKERRQ(ierr);
      ierr = MatTranspose(opflow->Jac_Ge,MAT_REUSE_MATRIX,&opflowipopt->Jac_GeT);CHKERRQ(ierr);
      ierr = MatGetSize(opflowipopt->Jac_GeT,&nrow,&ncol);CHKERRQ(ierr);
      aij = (Mat_SeqAIJ*)opflowipopt->Jac_GeT->data;

      ierr = PetscMemcpy(opflowipopt->jac_ge->values,aij->a,aij->nz*sizeof(PetscScalar));CHKERRQ(ierr);

      CCMatrixToMatrixMarketValuesOnly(opflowipopt->jac_ge,opflowipopt->nnz_jac_ge,valuesi);
      valuesi += opflowipopt->nnz_jac_ge;

      if(opflow->Nconineq) {
	/* Compute inequality constraint jacobian */
	ierr = (*opflow->formops.computeinequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Gi);CHKERRQ(ierr);

	ierr = MatTranspose(opflow->Jac_Gi,MAT_REUSE_MATRIX,&opflowipopt->Jac_GiT);CHKERRQ(ierr);
	aij = (Mat_SeqAIJ*)opflowipopt->Jac_GiT->data;

	ierr = PetscMemcpy(opflowipopt->jac_gi->values,aij->a,aij->nz*sizeof(PetscScalar));CHKERRQ(ierr);
      
	CCMatrixToMatrixMarketValuesOnly(opflowipopt->jac_gi,opflowipopt->nnz_jac_gi,valuesi);
	valuesi += opflowipopt->nnz_jac_gi;
      }

      if(scopflow->nconineqcoup[i]) {
	ps = opflow->ps;
	ps0 = opflow0->ps;
	for(j=0; j < ps->nbus; j++) {
	  bus = &ps->bus[j];
	  bus0 = &ps0->bus[j];
	  ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
	  ierr = PSBUSGetVariableLocation(bus0,&loc0);CHKERRQ(ierr);
	  for(k=0; k < bus->ngen; k++) {
	    ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	    ierr = PSBUSGetGen(bus0,k,&gen0);CHKERRQ(ierr);
	    if(!gen->status) {
	      if(gen0->status) loc0 += 2;
	      continue;
	    } else {
	      loc += 2;
	      if(!gen0->status) continue;
	      loc0 += 2;
	    }

	    x0loc = scopflowipopt->xstarti[0] + loc0;
	    xiloc = scopflowipopt->xstarti[i] + loc;
	    valuesi[0] = -1;
	    valuesi[1] = 1;
	    valuesi += 2;
	  }
	}
      }

      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
    }
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
  SCOPFLOWSolver_IPOPT scopflowipopt=(SCOPFLOWSolver_IPOPT)scopflow->solver;
  OPFLOW          opflow;
  OPFLOWSolver_IPOPT opflowipopt;
  PetscScalar     *xi,*valuesi=values,*lameqi,*lamineqi;
  PetscInt         i;
  PetscInt        *iRowStart=iRow,*jColStart=jCol;
  PetscInt        roffset;

  scopflow->obj_factor = obj_factor;

  if(values == NULL) {

    for(i=0; i < scopflow->Ns; i++) {
      opflow = scopflow->opflows[i];
      opflowipopt = (OPFLOWSolver_IPOPT)opflow->solver;

      roffset = scopflowipopt->xstarti[i];

      ierr = MatGetSize(opflowipopt->Hes_sbaij,&nrow,NULL);CHKERRQ(ierr);
      sbaij = (Mat_SeqSBAIJ*)opflowipopt->Hes_sbaij->data;
      ierr = PetscMemcpy(opflowipopt->hes->rowidx,sbaij->j,sbaij->nz*sizeof(PetscInt));CHKERRQ(ierr);
      ierr = PetscMemcpy(opflowipopt->hes->colptr,sbaij->i,(nrow+1)*sizeof(PetscInt));CHKERRQ(ierr);
      ierr = PetscMemcpy(opflowipopt->hes->values,sbaij->a,sbaij->nz*sizeof(PetscScalar));CHKERRQ(ierr);
      CCMatrixToMatrixMarketLocationsOnly(opflowipopt->hes,opflow->nx,iRowStart,jColStart,roffset,roffset,opflowipopt->nnz_hes);

      iRowStart += opflowipopt->nnz_hes;
      jColStart += opflowipopt->nnz_hes;
    }
  } else {

    for(i=0; i < scopflow->Ns; i++) {
      opflow = scopflow->opflows[i];
      opflowipopt = (OPFLOWSolver_IPOPT)opflow->solver;
      opflow->obj_factor = obj_factor;
      if(!scopflow->first_stage_gen_cost_only) opflow->compute_obj_hessian = PETSC_TRUE;
      else {
	if(i == 0) opflow->compute_obj_hessian = PETSC_TRUE;
	else opflow->compute_obj_hessian = PETSC_FALSE;
      }

      roffset = scopflowipopt->xstarti[i];

      xi = x + roffset;
      lameqi = lambda + scopflowipopt->gstarti[i];

      ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
      ierr = VecPlaceArray(opflow->Lambdae,lameqi);CHKERRQ(ierr);
      if(opflow->Nconineq) {
	lamineqi = lameqi + opflow->nconeq;
	ierr = VecPlaceArray(opflow->Lambdai,lamineqi);CHKERRQ(ierr);
      }

      ierr = (*opflow->formops.computehessian)(opflow,opflow->X,opflow->Lambdae,opflow->Lambdai,opflow->Hes);CHKERRQ(ierr);
      ierr = MatConvert(opflow->Hes,MATSEQSBAIJ,MAT_REUSE_MATRIX,&opflowipopt->Hes_sbaij);CHKERRQ(ierr);
      /* Since the Hessian is symmetric, we don't need to convert it to column compressed sparse format */
      sbaij = (Mat_SeqSBAIJ*)opflowipopt->Hes_sbaij->data;
      ierr = PetscMemcpy(opflowipopt->hes->values,sbaij->a,sbaij->nz*sizeof(PetscScalar));CHKERRQ(ierr);
      CCMatrixToMatrixMarketValuesOnly(opflowipopt->hes,opflowipopt->nnz_hes,valuesi);
      valuesi += opflowipopt->nnz_hes;

      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
      ierr = VecResetArray(opflow->Lambdae);CHKERRQ(ierr);
      if(opflow->Nconineq) {
	ierr = VecResetArray(opflow->Lambdai);CHKERRQ(ierr);
      }
    }
  }

  return 1;
}

PetscErrorCode SCOPFLOWSolverSolve_IPOPT(SCOPFLOW scopflow)
{
  PetscErrorCode     ierr;
  SCOPFLOWSolver_IPOPT scopflowipopt = (SCOPFLOWSolver_IPOPT)scopflow->solver;
  OPFLOW             opflow;
  OPFLOWSolver_IPOPT opflowipopt;
  Mat_SeqAIJ         *aij;
  Mat_SeqSBAIJ       *sbaij;
  PetscScalar        *x,*xl,*xu,*gl,*gu,*xi,*lameqi,*lamineqi,*lam;
  PetscInt           i;

  PetscFunctionBegin;

  scopflowipopt->nnz_jac_ge = scopflowipopt->nnz_jac_gi = 0;

  ierr = VecGetArray(scopflow->X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Lambda,&lam);CHKERRQ(ierr);

  for(i=0; i < scopflow->Ns; i++) {
    opflow = scopflow->opflows[i];
    opflowipopt = (OPFLOWSolver_IPOPT)opflow->solver;
    xi = x + scopflowipopt->xstarti[i];

    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);

    /* Compute nonzeros for the Jacobian */
    /* Equality constraint Jacobian */
    ierr = (*opflow->formops.computeequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Ge);CHKERRQ(ierr);
    /* Transpose the matrix to convert it to column compressed sparse format */
    ierr = MatTranspose(opflow->Jac_Ge,MAT_INITIAL_MATRIX,&opflowipopt->Jac_GeT);CHKERRQ(ierr);
    ierr = MatSetOption(opflowipopt->Jac_GeT,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);
    aij = (Mat_SeqAIJ*)opflowipopt->Jac_GeT->data;

    opflowipopt->nnz_jac_ge = aij->nz;
    scopflowipopt->nnz_jac_ge += opflowipopt->nnz_jac_ge;

    /* Create ccmatrix object for equality constrained Jacobian */
    ierr = PetscCalloc1(1,&opflowipopt->jac_ge);CHKERRQ(ierr);
    ierr = PetscCalloc1(opflow->nx+1,&opflowipopt->jac_ge->colptr);CHKERRQ(ierr);
    ierr = PetscCalloc1(opflowipopt->nnz_jac_ge,&opflowipopt->jac_ge->rowidx);CHKERRQ(ierr);
    ierr = PetscCalloc1(opflowipopt->nnz_jac_ge,&opflowipopt->jac_ge->values);CHKERRQ(ierr);

    opflowipopt->nnz_jac_gi = 0;
    if(opflow->Nconineq) {
      ierr = (*opflow->formops.computeinequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Gi);CHKERRQ(ierr);
      /* Transpose the matrix to convert it to column compressed sparse format */
      ierr = MatTranspose(opflow->Jac_Gi,MAT_INITIAL_MATRIX,&opflowipopt->Jac_GiT);CHKERRQ(ierr);
      ierr = MatSetOption(opflowipopt->Jac_GiT,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);
      aij = (Mat_SeqAIJ*)opflowipopt->Jac_GiT->data;
      opflowipopt->nnz_jac_gi = aij->nz;
      scopflowipopt->nnz_jac_gi += opflowipopt->nnz_jac_gi;

      /* Create ccmatrix object for inequality constrained Jacobian */
      ierr = PetscCalloc1(1,&opflowipopt->jac_gi);CHKERRQ(ierr);
      ierr = PetscCalloc1(opflow->nx+1,&opflowipopt->jac_gi->colptr);CHKERRQ(ierr);
      ierr = PetscCalloc1(opflowipopt->nnz_jac_gi,&opflowipopt->jac_gi->rowidx);CHKERRQ(ierr);
      ierr = PetscCalloc1(opflowipopt->nnz_jac_gi,&opflowipopt->jac_gi->values);CHKERRQ(ierr);
    }   
    opflowipopt->nnz_jac_g = opflowipopt->nnz_jac_ge + opflowipopt->nnz_jac_gi;

    /* Add non-zeros for Jacobian of coupling constraints */
    if(scopflow->nconineqcoup[i]) scopflowipopt->nnz_jac_gi += 2*scopflow->nconineqcoup[i];


    /* Compute non-zeros for Hessian */

    lameqi = lam + scopflowipopt->gstarti[i];
    ierr = VecPlaceArray(opflow->Lambdae,lameqi);CHKERRQ(ierr);

    if(opflow->Nconineq) {
      lamineqi = lam + scopflowipopt->gstarti[i]+opflow->nconeq;
      ierr = VecPlaceArray(opflow->Lambdai,lamineqi);CHKERRQ(ierr);
    } else {
      opflow->Lambdai = NULL;
    }

    ierr = (*opflow->formops.computehessian)(opflow,opflow->X,opflow->Lambdae,opflow->Lambdai,opflow->Hes);CHKERRQ(ierr);
    /* Convert matrix to symmetric sbaij format needed for the IPOPT solver */
    ierr = MatSetOption(opflow->Hes,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
    ierr = MatConvert(opflow->Hes,MATSEQSBAIJ,MAT_INITIAL_MATRIX,&opflowipopt->Hes_sbaij);CHKERRQ(ierr);
    /* Since the Hessian is symmetric, we don't need to convert it to column compressed sparse format */
    sbaij = (Mat_SeqSBAIJ*)opflowipopt->Hes_sbaij->data;
    opflowipopt->nnz_hes = sbaij->nz;
    scopflowipopt->nnz_hes += opflowipopt->nnz_hes;

    /* Create ccmatrix object for hessian */
    ierr = PetscCalloc1(1,&opflowipopt->hes);CHKERRQ(ierr);
    ierr = PetscCalloc1(opflow->nx+1,&opflowipopt->hes->colptr);CHKERRQ(ierr);
    ierr = PetscCalloc1(opflowipopt->nnz_hes,&opflowipopt->hes->rowidx);CHKERRQ(ierr);
    ierr = PetscCalloc1(opflowipopt->nnz_hes,&opflowipopt->hes->values);CHKERRQ(ierr);

    ierr = VecResetArray(opflow->Lambdae);CHKERRQ(ierr);
    if(opflow->Nconineq) {
      ierr = VecResetArray(opflow->Lambdai);CHKERRQ(ierr);
    }
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(scopflow->X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Lambda,&lam);CHKERRQ(ierr);

  scopflowipopt->nnz_jac_g = scopflowipopt->nnz_jac_ge + scopflowipopt->nnz_jac_gi;

  /* Create IPOPT solver instance */
  ierr = VecGetArray(scopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Gu,&gu);CHKERRQ(ierr);

  /* Create IPOPT problem */
  scopflowipopt->nlp = CreateIpoptProblem(scopflow->Nx,xl,xu,scopflow->Ncon,gl,gu,scopflowipopt->nnz_jac_g,scopflowipopt->nnz_hes,0,&eval_scopflow_f,
				   &eval_scopflow_g,&eval_scopflow_grad_f,
				   &eval_scopflow_jac_g,&eval_scopflow_h);

  ierr = VecRestoreArray(scopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Gu,&gu);CHKERRQ(ierr);

  ierr = VecGetArray(scopflow->X,&x);CHKERRQ(ierr);

  /* Solve */
  scopflowipopt->solve_status = IpoptSolve(scopflowipopt->nlp,x,NULL,&scopflow->obj,NULL,NULL,NULL,scopflow);

  ierr = VecRestoreArray(scopflow->X,&x);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverDestroy_IPOPT(SCOPFLOW scopflow)
{
  PetscErrorCode     ierr;
  SCOPFLOWSolver_IPOPT ipopt = (SCOPFLOWSolver_IPOPT)scopflow->solver;

  PetscFunctionBegin;

  if(ipopt->nlp) {
    FreeIpoptProblem(ipopt->nlp);
    ipopt->nlp = NULL;
  }

  ierr = PetscFree(ipopt->xstarti);CHKERRQ(ierr);
  ierr = PetscFree(ipopt->gstarti);CHKERRQ(ierr);
  ierr = PetscFree(ipopt->nxi);CHKERRQ(ierr);
  ierr = PetscFree(ipopt->ngi);CHKERRQ(ierr);

  ierr = PetscFree(ipopt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverSetUp_IPOPT(SCOPFLOW scopflow)
{
  PetscErrorCode ierr;
  SCOPFLOWSolver_IPOPT ipopt = (SCOPFLOWSolver_IPOPT)scopflow->solver;
  PetscInt       i,j,k,ctr;
  OPFLOW         opflow,opflow0;
  PetscInt       ngenON;
  PetscScalar    *x,*xi,*xl,*xu,*xli,*xui,*gl,*gu,*gli,*gui;
  PS             ps,ps0;
  PSBUS          bus,bus0;
  PSGEN          gen,gen0;

  PetscFunctionBegin;

  ierr = PetscCalloc1(scopflow->Ns,&ipopt->xstarti);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->Ns,&ipopt->gstarti);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->Ns,&ipopt->nxi);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->Ns,&ipopt->ngi);CHKERRQ(ierr);

  scopflow->Nx = 0;
  scopflow->Ncon = 0;
  ipopt->xstarti[0] = 0;
  ipopt->gstarti[0] = 0;

  for(i=0; i < scopflow->Ns; i++) {
    opflow = scopflow->opflows[i];
    ierr = PSGetNumActiveGenerators(opflow->ps,&ngenON,NULL);CHKERRQ(ierr);
    ipopt->nxi[i] = opflow->nx;
    if(scopflow->iscoupling) scopflow->nconineqcoup[i] = (i == 0)?0:ngenON;
    else scopflow->nconineqcoup[i] = 0;

    ipopt->ngi[i] = opflow->ncon + scopflow->nconineqcoup[i];
    if(i < scopflow->Ns - 1) {
      ipopt->xstarti[i+1] = ipopt->xstarti[i] + ipopt->nxi[i];
      ipopt->gstarti[i+1] = ipopt->gstarti[i] + ipopt->ngi[i];
    }
    scopflow->Nx += ipopt->nxi[i];
    scopflow->Ncon += ipopt->ngi[i];
  }

  /* Create vector X */
  ierr = VecCreate(scopflow->comm->type,&scopflow->X);CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->X,scopflow->Nx,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(scopflow->X);CHKERRQ(ierr);

  ierr = VecDuplicate(scopflow->X,&scopflow->Xl);CHKERRQ(ierr);
  ierr = VecDuplicate(scopflow->X,&scopflow->Xu);CHKERRQ(ierr);

  /* vector for constraints */
  ierr = VecCreate(scopflow->comm->type,&scopflow->G);CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->G,scopflow->Ncon,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(scopflow->G);CHKERRQ(ierr);

  /* Constraint bounds vectors  */
  ierr = VecDuplicate(scopflow->G,&scopflow->Gl);CHKERRQ(ierr);
  ierr = VecDuplicate(scopflow->G,&scopflow->Gu);CHKERRQ(ierr);

  /* Lagrangian multipliers */
  ierr = VecDuplicate(scopflow->G,&scopflow->Lambda);CHKERRQ(ierr);

  /* Set Initial guess and Bounds on variables and constraints */
  ierr = VecGetArray(scopflow->X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Gu,&gu);CHKERRQ(ierr);

  opflow0 = scopflow->opflows[0];
  for(i=0; i < scopflow->Ns; i++) {
    opflow = scopflow->opflows[i];
    /* Set initial guess and bounds on variables */
    xi  = x  + ipopt->xstarti[i];
    xli = xl + ipopt->xstarti[i];
    xui = xu + ipopt->xstarti[i];

    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Xl,xli);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Xu,xui);CHKERRQ(ierr);

    /* Set bounds */
    ierr = (*opflow->formops.setvariablebounds)(opflow,opflow->Xl,opflow->Xu);CHKERRQ(ierr);

    /* Set initial guess */
    ierr = OPFLOWSetInitialGuess(opflow,opflow->X);CHKERRQ(ierr);

    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Xl);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Xu);CHKERRQ(ierr);

    /* Set bounds on constraints */
    gli = gl + ipopt->gstarti[i];
    gui = gu + ipopt->gstarti[i];
    
    ierr = VecPlaceArray(opflow->Gl,gli);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Gu,gui);CHKERRQ(ierr);

    ierr = (*opflow->formops.setconstraintbounds)(opflow,opflow->Gl,opflow->Gu);CHKERRQ(ierr);

    ierr = VecResetArray(opflow->Gl);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Gu);CHKERRQ(ierr);
    
    if(scopflow->nconineqcoup[i]) {
      ctr = 0;
      ps = opflow->ps;
      ps0 = opflow0->ps;
      /* Bounds on coupling constraints */
      for(j=0; j < ps->nbus; j++) {
	bus = &ps->bus[j];
	bus0 = &ps0->bus[j];

	for(k=0; k < bus->ngen; k++) {
	  ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	  ierr = PSBUSGetGen(bus0,k,&gen0);CHKERRQ(ierr);
	  if(!gen->status || !gen0->status) continue;
	  /* Generator can do a full ramp up to its max. capacity */
	  gli[opflow->ncon + ctr] = -gen->pt;
	  gui[opflow->ncon + ctr] =  gen->pt;
	  ctr++;
	}
      }
    }
  }
  
  ierr = VecRestoreArray(scopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Gu,&gu);CHKERRQ(ierr);

  /* Initialize Lagrange multiplier */
  ierr = VecSet(scopflow->Lambda,1.0);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverCreate_IPOPT(SCOPFLOW scopflow)
{
  PetscErrorCode ierr;
  SCOPFLOWSolver_IPOPT ipopt;
  
  PetscFunctionBegin;

  if(scopflow->comm->size > 1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"IPOPT solver does not support execution in parallel\n",scopflow->comm->size); 
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
