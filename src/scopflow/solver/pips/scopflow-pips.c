#if defined(SCOPFLOW_HAVE_PIPS)

#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>
#include "scopflow-pips.h"

int scopflow_init_x0(double* x0, CallBackDataPtr cbd) {
  int row = cbd->row_node_id;
  SCOPFLOW scopflow=cbd->prob;
  OPFLOW   opflow = scopflow->opflows[row];
  PetscErrorCode ierr;
  
  ierr = VecPlaceArray(opflow->X,x0);CHKERRQ(ierr);
  ierr = (*opflow->formops.setinitialguess)(opflow,opflow->X);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
	
  return 1;
}

int scopflow_prob_info(int* n, double* col_lb, double* col_ub, int* m,
		  double* row_lb, double* row_ub, CallBackDataPtr cbd) {
  int row = cbd->row_node_id;
  int col = cbd->col_node_id;
  SCOPFLOW scopflow = (SCOPFLOW)cbd->prob;
  OPFLOW   opflow = scopflow->opflows[row];
  PetscInt rank=scopflow->comm->rank;
  int type = cbd->typeflag;
  PetscErrorCode ierr;
  Vec            Xl,Xu,Gl,Gu;
  PS     ps;
  PSBUS  bus;
  PSGEN  gen;
  PetscInt ctr;
  PetscInt j,k;

  if(type == 1) {
    if(row_lb == NULL){
      *m = 0;
    }
    return 1;
  }
  
  /* Set sizes of variables and constraints */
  if(col_lb == NULL){
    *n = opflow->nx; 
    *m = opflow->ncon + scopflow->nconineqcoup[row];

    scopflow->ns++;

    ierr = PetscPrintf(PETSC_COMM_SELF,"Rank %d row %d col %d type %d Nvar=%d Ncon = %d\n",rank,row,col,type,*n,*m);CHKERRQ(ierr);
  } else {

    /* Bounds on variables and constraints */
    Xl = scopflow->opflows[row]->Xl;
    Xu = scopflow->opflows[row]->Xu;
    Gl = scopflow->opflows[row]->Gl;
    Gu = scopflow->opflows[row]->Gu;

    ierr = VecPlaceArray(Xl,col_lb);CHKERRQ(ierr);
    ierr = VecPlaceArray(Xu,col_ub);CHKERRQ(ierr);
    ierr = VecPlaceArray(Gl,row_lb);CHKERRQ(ierr);
    ierr = VecPlaceArray(Gu,row_ub);CHKERRQ(ierr);
    
    /* Set variable bounds */
    ierr = (*opflow->formops.setvariablebounds)(opflow,opflow->Xl,opflow->Xu);CHKERRQ(ierr);

    /* Set constraint bounds */
    ierr = (*opflow->formops.setconstraintbounds)(opflow,opflow->Gl,opflow->Gu);CHKERRQ(ierr);
    
    ierr = VecResetArray(Xl);CHKERRQ(ierr);
    ierr = VecResetArray(Xu);CHKERRQ(ierr);
    ierr = VecResetArray(Gl);CHKERRQ(ierr);
    ierr = VecResetArray(Gu);CHKERRQ(ierr);
  

    if(scopflow->nconineqcoup[row]) {
      /* Coupling constraint bounds */
      ctr = 0;
      ps = opflow->ps;
      /* Bounds on coupling constraints */
      for(j=0; j < ps->nbus; j++) {
	bus = &ps->bus[j];
	
	for(k=0; k < bus->ngen; k++) {
	  ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	  /* Generator can do a full ramp up to its max. capacity */
	  row_lb[opflow->ncon + ctr] = -gen->pt;
	  row_ub[opflow->ncon + ctr] =  gen->pt;
	  ctr++;
	}
      }
    }

    /* Copy over col_lb and col_ub to vectors Xl and Xu so
       that initialization of X can work correctly as
       Xinital = (Xl + Xu)/2
    */
    PetscScalar *xl,*xu;
    ierr = VecGetArray(Xl,&xl);CHKERRQ(ierr);
    ierr = VecGetArray(Xu,&xu);CHKERRQ(ierr);
    ierr = PetscMemcpy(xl,col_lb,opflow->nx*sizeof(PetscScalar));CHKERRQ(ierr);
    ierr = PetscMemcpy(xu,col_ub,opflow->nx*sizeof(PetscScalar));CHKERRQ(ierr);
    ierr = VecRestoreArray(Xl,&xl);CHKERRQ(ierr);
    ierr = VecRestoreArray(Xu,&xu);CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_SELF,"Rank %d: row = %d, col = %d set up problem info\n",rank,row,col);CHKERRQ(ierr);
  }
  return 1;
}

int scopflow_write_solution(double* x, double* lam_eq, double* lam_ieq,CallBackDataPtr cbd)
{
  int row = cbd->row_node_id;
  if(row == 0) {
  } else if(row == 1 || row == 2) {
  }
  return 1;
}

int scopflow_eval_f(double* x0, double* x1, double* obj, CallBackDataPtr cbd) 
{
  PetscErrorCode ierr;
  int row = cbd->row_node_id;
  int col = cbd->col_node_id;
  SCOPFLOW scopflow=(SCOPFLOW)cbd->prob;
  OPFLOW   opflow=scopflow->opflows[row];
  
  *obj = 0.0;
  if(row == col) {
    if(row == 0) {
      ierr = VecPlaceArray(opflow->X,x0);CHKERRQ(ierr);
      ierr = (*opflow->formops.computeobjective)(opflow,opflow->X,obj);CHKERRQ(ierr);
      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
    } else {
      if(!scopflow->first_stage_gen_cost_only) {
	ierr = VecPlaceArray(opflow->X,x1);CHKERRQ(ierr);
	ierr = (*opflow->formops.computeobjective)(opflow,opflow->X,obj);CHKERRQ(ierr);
	ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
      }
    }
  } else return 0;

  return 1;
}

int scopflow_eval_grad_f(double* x0, double* x1, double* grad, CallBackDataPtr cbd) 
{
  PetscErrorCode ierr;
  int row = cbd->row_node_id;
  int col = cbd->col_node_id;
  SCOPFLOW scopflow=(SCOPFLOW)cbd->prob;
  OPFLOW   opflow=scopflow->opflows[row];
  double   *x;
  
  if(row == col) {
    if(row == 0) {
      x = x0;
      ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
      ierr = VecPlaceArray(opflow->gradobj,grad);CHKERRQ(ierr);
      ierr = VecSet(opflow->gradobj,0.0);CHKERRQ(ierr);
      ierr = (*opflow->formops.computegradient)(opflow,opflow->X,opflow->gradobj);CHKERRQ(ierr);
      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
      ierr = VecResetArray(opflow->gradobj);CHKERRQ(ierr);
    } else {
      if(scopflow->first_stage_gen_cost_only) {
	ierr = VecPlaceArray(opflow->gradobj,grad);CHKERRQ(ierr);
	ierr = VecSet(opflow->gradobj,0.0);CHKERRQ(ierr);
	ierr = VecResetArray(opflow->gradobj);CHKERRQ(ierr);
	return 1;
      }
      x = x1;
      ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
      ierr = VecPlaceArray(opflow->gradobj,grad);CHKERRQ(ierr);
      ierr = VecSet(opflow->gradobj,0.0);CHKERRQ(ierr);
      ierr = (*opflow->formops.computegradient)(opflow,opflow->X,opflow->gradobj);CHKERRQ(ierr);
      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
      ierr = VecResetArray(opflow->gradobj);CHKERRQ(ierr);
    }
  } else {
    ierr = VecPlaceArray(opflow->gradobj,grad);CHKERRQ(ierr);
    ierr = VecSet(opflow->gradobj,0.0);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->gradobj);CHKERRQ(ierr);
  }
    
  return 1;
}

int scopflow_eval_h(double* x0, double* x1, double* lambda, int* nz, double* elts,
	       int* rowidx, int* colptr, CallBackDataPtr cbd) {

  PetscErrorCode ierr;
  int row = cbd->row_node_id;
  int col = cbd->col_node_id;
  SCOPFLOW scopflow=(SCOPFLOW)cbd->prob;
  Mat_SeqSBAIJ  *sbaij;
  OPFLOW   opflow = scopflow->opflows[row];
  OPFLOWSolver_IPOPT opflowipopt = opflow->solver;
  PetscInt nrow,ncol;
  PetscScalar *lameq,*lamineq;

  if(elts==NULL) {
    *nz = 0;
    if(row == col) {
      PetscScalar *x;
      if(row == 0) x = x0;
      else x = x1;

      *nz = opflowipopt->nnz_hes;
    }
  } else {
    if(row == col) {
      PetscScalar *x;
      if(row == 0) x = x0;
      else x = x1;

      ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);

      lameq = lambda;
      ierr = VecPlaceArray(opflow->Lambdae,lameq);CHKERRQ(ierr);
      if(opflow->Nconineq) {
	lamineq = lambda + opflow->nconeq;
	ierr = VecPlaceArray(opflow->Lambdai,lamineq);CHKERRQ(ierr);
      }

      /* Compute Hessian */
      ierr = (*opflow->formops.computehessian)(opflow,opflow->X,opflow->Lambdae,opflow->Lambdai,opflow->Hes);CHKERRQ(ierr);
      ierr = MatConvert(opflow->Hes,MATSEQSBAIJ,MAT_REUSE_MATRIX,&opflowipopt->Hes_sbaij);CHKERRQ(ierr);

      /* Since the Hessian is symmetric, we don't need to convert it to column compressed sparse format */
      sbaij = (Mat_SeqSBAIJ*)opflowipopt->Hes_sbaij->data;

      ierr = MatGetSize(opflow->Hes,&nrow,&ncol);CHKERRQ(ierr);

      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
      ierr = VecResetArray(opflow->Lambdae);CHKERRQ(ierr);
      if(opflow->Nconineq) {
	ierr = VecResetArray(opflow->Lambdai);CHKERRQ(ierr);
      }

      ierr = PetscMemcpy(rowidx,sbaij->j,sbaij->nz*sizeof(PetscInt));CHKERRQ(ierr);
      ierr = PetscMemcpy(colptr,sbaij->i,(nrow+1)*sizeof(PetscInt));CHKERRQ(ierr);
      ierr = PetscMemcpy(elts,sbaij->a,sbaij->nz*sizeof(PetscScalar));CHKERRQ(ierr);

      /*      ierr = MatView(opflow->Hes,0);CHKERRQ(ierr);
      exit(1);
      */
    }
  }
  
  return 1;
}

int scopflow_eval_g(double* x0, double* x1, double* eq_g, double* inq_g,
	       CallBackDataPtr cbd) 
{
  PetscErrorCode ierr;
  int row = cbd->row_node_id;
  SCOPFLOW scopflow=(SCOPFLOW)cbd->prob;
  OPFLOW   opflow=scopflow->opflows[row];
  double   *x;
  PS       ps;
  PSBUS    bus;
  PSGEN    gen;
  PetscScalar *inqcoup_g;
  PetscInt ctr,j,k,loc;

  if(row == 0) x = x0;
  else x = x1;


  ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);

  /* Equality constraints */
  ierr = VecPlaceArray(opflow->Ge,eq_g);CHKERRQ(ierr);
  ierr = (*opflow->formops.computeequalityconstraints)(opflow,opflow->X,opflow->Ge);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->Ge);CHKERRQ(ierr);

  if(opflow->Nconineq) {
    /* Inequality Constraints */
    ierr = VecPlaceArray(opflow->Gi,inq_g);CHKERRQ(ierr);
    ierr = (*opflow->formops.computeinequalityconstraints)(opflow,opflow->X,opflow->Gi);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Gi);CHKERRQ(ierr);
  }

  if(scopflow->nconineqcoup[row]) {
    inqcoup_g = inq_g + opflow->nconineq;
    ctr = 0;
    ps = opflow->ps;
    for(j=0; j < ps->nbus; j++) {
      bus = &ps->bus[j];
      ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
      
      for(k=0; k < bus->ngen; k++) {
	loc += 2;
	ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	inqcoup_g[ctr++] = x1[loc] - x0[loc]; /* PGi - PG0 */
      }
    }
  }

  ierr = VecResetArray(opflow->X);CHKERRQ(ierr);    

  return 1;
}


int scopflow_eval_jac_g(double* x0, double* x1, int* e_nz, double* e_elts,
		int* e_rowidx, int* e_colptr, int* i_nz, double* i_elts, int* i_rowidx,
		int* i_colptr, CallBackDataPtr cbd) {
  int row = cbd->row_node_id;
  int col = cbd->col_node_id;
  PetscErrorCode ierr;
  SCOPFLOW scopflow=(SCOPFLOW)cbd->prob;
  SCOPFLOWSolver_PIPS pips=scopflow->solver;
  /*  PetscInt rank=scopflow->comm->rank; */
  OPFLOW   opflow=scopflow->opflows[row];
  OPFLOWSolver_IPOPT opflowipopt = opflow->solver;
  Mat_SeqAIJ *aij;
  PetscInt    nrow,ncol;
  PS   ps;
  PSBUS bus;
  PetscInt ridx,cidx,roffset,j,k,loc;
  PetscScalar val;

  if(e_elts==NULL && i_elts == NULL) {
    /* Number of non-zeros in equality and inequality constraint Jacobian */
    opflow = scopflow->opflows[row];
    if(row == col) {
      PetscScalar *x;
      if(row == 0) x = x0;
      else x = x1;

      *e_nz = opflowipopt->nnz_jac_ge;
      *i_nz = opflowipopt->nnz_jac_gi; /* Includes coupling non-zeros */

    } else {
      if(col == 0) {
	*e_nz = 0;
	*i_nz = 0;
	if(scopflow->nconineqcoup[row])	*i_nz = scopflow->nconineqcoup[row];
      }
    }
  } else {
    if(row == col) {
      PetscScalar *x;
      if(row == 0) x = x0;
      else x = x1;

      ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);

      /* Equality constraints Jacobian */
      ierr = (*opflow->formops.computeequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Ge);CHKERRQ(ierr);
      /* Transpose the matrix to convert it to column compressed sparse format */
      ierr = MatTranspose(opflow->Jac_Ge,MAT_REUSE_MATRIX,&opflowipopt->Jac_GeT);CHKERRQ(ierr);

      ierr = MatGetSize(opflowipopt->Jac_GeT,&nrow,&ncol);CHKERRQ(ierr);
      aij = (Mat_SeqAIJ*)opflowipopt->Jac_GeT->data;

      ierr = PetscMemcpy(e_rowidx,aij->j,aij->nz*sizeof(PetscInt));CHKERRQ(ierr);
      ierr = PetscMemcpy(e_colptr,aij->i,(nrow+1)*sizeof(PetscInt));CHKERRQ(ierr);
      ierr = PetscMemcpy(e_elts,aij->a,aij->nz*sizeof(PetscScalar));CHKERRQ(ierr);

      if(opflow->Nconineq + scopflow->nconineqcoup[row]) {
	ierr = MatZeroEntries(opflow->Jac_Gi);CHKERRQ(ierr);
	if(opflow->Nconineq) {
	  /* Inequality constrained Jacobian */
	  ierr = (*opflow->formops.computeinequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Gi);CHKERRQ(ierr);
	}

	if(scopflow->nconineqcoup[row]) {
	  ps = opflow->ps;
	  roffset = opflow->Nconineq;
	  for(j=0; j < ps->nbus; j++) {
	    bus = &ps->bus[j];
	    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
	    for(k=0; k < bus->ngen; k++) {
	      loc += 2;
	      ridx = roffset; cidx = loc;
	      val  = 1.0;
	      ierr = MatSetValues(opflow->Jac_Gi,1,&ridx,1,&cidx,&val,INSERT_VALUES);CHKERRQ(ierr);
	      roffset += 1;
	    }
	  }
	  ierr = MatAssemblyBegin(opflow->Jac_Gi,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	  ierr = MatAssemblyEnd(opflow->Jac_Gi,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}

	/* Transpose the matrix to convert it to column compressed sparse format */
	ierr = MatTranspose(opflow->Jac_Gi,MAT_REUSE_MATRIX,&opflowipopt->Jac_GiT);CHKERRQ(ierr);	
	ierr = MatGetSize(opflowipopt->Jac_GiT,&nrow,&ncol);CHKERRQ(ierr);
	aij = (Mat_SeqAIJ*)opflowipopt->Jac_GiT->data;
	ierr = PetscMemcpy(i_rowidx,aij->j,aij->nz*sizeof(PetscInt));CHKERRQ(ierr);
	ierr = PetscMemcpy(i_colptr,aij->i,(nrow+1)*sizeof(PetscInt));CHKERRQ(ierr);
	ierr = PetscMemcpy(i_elts,aij->a,aij->nz*sizeof(PetscScalar));CHKERRQ(ierr);
      }
      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);

    } else {
      if(scopflow->nconineqcoup[row] && col == 0) {
	ierr = MatGetSize(pips->Jac_GicoupT[row],&nrow,&ncol);CHKERRQ(ierr);
	aij = (Mat_SeqAIJ*)pips->Jac_GicoupT[row]->data;
	ierr = PetscMemcpy(i_rowidx,aij->j,aij->nz*sizeof(PetscInt));CHKERRQ(ierr);
	ierr = PetscMemcpy(i_colptr,aij->i,(nrow+1)*sizeof(PetscInt));CHKERRQ(ierr);
	ierr = PetscMemcpy(i_elts,aij->a,aij->nz*sizeof(PetscScalar));CHKERRQ(ierr);
      }
    }
  }
  
  return 1;
}

PetscErrorCode SCOPFLOWSolverSolve_PIPS(SCOPFLOW scopflow)
{
  PetscErrorCode     ierr;
  SCOPFLOWSolver_PIPS pips=scopflow->solver;
  OPFLOW             opflow;
  OPFLOWSolver_IPOPT opflowipopt;
  Mat_SeqAIJ         *aij;
  Mat_SeqSBAIJ       *sbaij;
  PetscScalar        *x,*xi,*lameqi,*lamineqi,*lam;
  PetscInt           i;

  PetscFunctionBegin;

  pips->nnz_jac_ge = pips->nnz_jac_gi = 0;

  ierr = VecGetArray(scopflow->X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Lambda,&lam);CHKERRQ(ierr);

  ierr = PetscCalloc1(scopflow->Ns,&pips->Jac_GicoupT);CHKERRQ(ierr);
  for(i=0; i < scopflow->Ns; i++) {
    opflow = scopflow->opflows[i];
    opflowipopt = opflow->solver;
    xi = x + pips->xstarti[i];

    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);

    /* Compute nonzeros for the Jacobian */
    /* Equality constraint Jacobian */
    ierr = (*opflow->formops.computeequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Ge);CHKERRQ(ierr);
    /* Transpose the matrix to convert it to column compressed sparse format */
    ierr = MatTranspose(opflow->Jac_Ge,MAT_INITIAL_MATRIX,&opflowipopt->Jac_GeT);CHKERRQ(ierr);
    ierr = MatSetOption(opflowipopt->Jac_GeT,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);
    aij = (Mat_SeqAIJ*)opflowipopt->Jac_GeT->data;

    opflowipopt->nnz_jac_ge = aij->nz;
    pips->nnz_jac_ge += opflowipopt->nnz_jac_ge;

    /* Create ccmatrix object for equality constrained Jacobian */
    ierr = PetscCalloc1(1,&opflowipopt->jac_ge);CHKERRQ(ierr);
    ierr = PetscMalloc1(opflow->nx+1,&opflowipopt->jac_ge->colptr);CHKERRQ(ierr);
    ierr = PetscMalloc1(opflowipopt->nnz_jac_ge,&opflowipopt->jac_ge->rowidx);CHKERRQ(ierr);
    ierr = PetscMalloc1(opflowipopt->nnz_jac_ge,&opflowipopt->jac_ge->values);CHKERRQ(ierr);

    opflowipopt->nnz_jac_gi = 0;
    if(opflow->Nconineq + scopflow->nconineqcoup[i]) {
      /* Create a bigger inequality constrained Jacobian matrix that includes the coupling constraints */
      /* Also update the number of inequality constraints */
      ierr = MatDestroy(&opflow->Jac_Gi);CHKERRQ(ierr);

      ierr = MatCreate(opflow->comm->type,&opflow->Jac_Gi);CHKERRQ(ierr);
      ierr = MatSetSizes(opflow->Jac_Gi,opflow->Nconineq+scopflow->nconineqcoup[i],opflow->Nx,opflow->Nconineq+scopflow->nconineqcoup[i],opflow->Nx);CHKERRQ(ierr);
      ierr = MatSetUp(opflow->Jac_Gi);CHKERRQ(ierr);
      ierr = MatSetFromOptions(opflow->Jac_Gi);CHKERRQ(ierr);

      if(opflow->Nconineq) {
	ierr = (*opflow->formops.computeinequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Gi);CHKERRQ(ierr);
      }

      /* Add non-zeros for Jacobian of coupling constraints */
      if(scopflow->nconineqcoup[i]) {
	PS ps;
	PSBUS bus;
	PetscInt j,k,roffset,loc,ridx,cidx;
	PetscScalar val;
	ps = opflow->ps;
	roffset = opflow->Nconineq;
	for(j=0; j < ps->nbus; j++) {
	  bus = &ps->bus[j];
	  ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
	  for(k=0; k < bus->ngen; k++) {
	    loc += 2;
	    ridx = roffset; cidx = loc;
	    val  = 1.0;
	    ierr = MatSetValues(opflow->Jac_Gi,1,&ridx,1,&cidx,&val,INSERT_VALUES);CHKERRQ(ierr);
	    roffset += 1;
	  }
	}
	ierr = MatAssemblyBegin(opflow->Jac_Gi,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(opflow->Jac_Gi,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      }

      /* Transpose the matrix to convert it to column compressed sparse format */
      ierr = MatTranspose(opflow->Jac_Gi,MAT_INITIAL_MATRIX,&opflowipopt->Jac_GiT);CHKERRQ(ierr);
      ierr = MatSetOption(opflowipopt->Jac_GiT,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);
      aij = (Mat_SeqAIJ*)opflowipopt->Jac_GiT->data;
      opflowipopt->nnz_jac_gi = aij->nz;

      pips->nnz_jac_gi += opflowipopt->nnz_jac_gi;

      /* Create ccmatrix object for inequality constrained Jacobian */
      ierr = PetscCalloc1(1,&opflowipopt->jac_gi);CHKERRQ(ierr);
      ierr = PetscMalloc1(opflow->nx+1,&opflowipopt->jac_gi->colptr);CHKERRQ(ierr);
      ierr = PetscMalloc1(opflowipopt->nnz_jac_gi,&opflowipopt->jac_gi->rowidx);CHKERRQ(ierr);
      ierr = PetscMalloc1(opflowipopt->nnz_jac_gi,&opflowipopt->jac_gi->values);CHKERRQ(ierr);

      /* Create coupling matrix */
      ierr = MatCreate(PETSC_COMM_SELF,&pips->Jac_GicoupT[i]);CHKERRQ(ierr);
      ierr = MatSetSizes(pips->Jac_GicoupT[i],opflow->nx,opflow->nconineq+scopflow->nconineqcoup[i],opflow->nx,opflow->nconineq+scopflow->nconineqcoup[i]);CHKERRQ(ierr);
      ierr = MatSetFromOptions(pips->Jac_GicoupT[i]);CHKERRQ(ierr);
      ierr = MatSetUp(pips->Jac_GicoupT[i]);CHKERRQ(ierr);
      /* Add non-zeros for Jacobian of coupling constraints */
      if(scopflow->nconineqcoup[i]) {
	PS ps;
	PSBUS bus;
	PetscInt j,k,roffset,loc,ridx,cidx;
	PetscScalar val;
	ps = opflow->ps;
	roffset = opflow->Nconineq;
	for(j=0; j < ps->nbus; j++) {
	  bus = &ps->bus[j];
	  ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
	  for(k=0; k < bus->ngen; k++) {
	    loc += 2;
	    ridx = roffset; cidx = loc;
	    val  = -1.0;
	    ierr = MatSetValues(pips->Jac_GicoupT[i],1,&cidx,1,&ridx,&val,INSERT_VALUES);CHKERRQ(ierr);
	    roffset += 1;
	  }
	}
	ierr = MatAssemblyBegin(pips->Jac_GicoupT[i],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(pips->Jac_GicoupT[i],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

      }

    }
    opflowipopt->nnz_jac_g = opflowipopt->nnz_jac_ge + opflowipopt->nnz_jac_gi;

    /* Compute non-zeros for Hessian */

    lameqi = lam + pips->gstarti[i];
    ierr = VecPlaceArray(opflow->Lambdae,lameqi);CHKERRQ(ierr);

    if(opflow->Nconineq) {
      lamineqi = lam + pips->gstarti[i]+opflow->nconeq;
      ierr = VecPlaceArray(opflow->Lambdai,lamineqi);CHKERRQ(ierr);
    }
    ierr = (*opflow->formops.computehessian)(opflow,opflow->X,opflow->Lambdae,opflow->Lambdai,opflow->Hes);CHKERRQ(ierr);
    /* Convert matrix to symmetric sbaij format needed for the PIPS solver */
    ierr = MatSetOption(opflow->Hes,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
    ierr = MatConvert(opflow->Hes,MATSEQSBAIJ,MAT_INITIAL_MATRIX,&opflowipopt->Hes_sbaij);CHKERRQ(ierr);
    /* Since the Hessian is symmetric, we don't need to convert it to column compressed sparse format */
    sbaij = (Mat_SeqSBAIJ*)opflowipopt->Hes_sbaij->data;
    opflowipopt->nnz_hes = sbaij->nz;
    pips->nnz_hes += opflowipopt->nnz_hes;

    /* Create ccmatrix object for hessian */
    ierr = PetscCalloc1(1,&opflowipopt->hes);CHKERRQ(ierr);
    ierr = PetscMalloc1(opflow->nx+1,&opflowipopt->hes->colptr);CHKERRQ(ierr);
    ierr = PetscMalloc1(opflowipopt->nnz_hes,&opflowipopt->hes->rowidx);CHKERRQ(ierr);
    ierr = PetscMalloc1(opflowipopt->nnz_hes,&opflowipopt->hes->values);CHKERRQ(ierr);

    ierr = VecResetArray(opflow->Lambdae);CHKERRQ(ierr);
    if(opflow->Nconineq) {
      ierr = VecResetArray(opflow->Lambdai);CHKERRQ(ierr);
    }
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(scopflow->X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Lambda,&lam);CHKERRQ(ierr);

  pips->nnz_jac_g = pips->nnz_jac_ge + pips->nnz_jac_gi;

  /* Create PIPS solver instance */
  pips->nlp = CreatePipsNlpProblemStruct(scopflow->comm->type, scopflow->Ns-1,
							       scopflow_init_x0, scopflow_prob_info, scopflow_eval_f, scopflow_eval_g, scopflow_eval_grad_f, scopflow_eval_jac_g,
						               scopflow_eval_h, scopflow_write_solution, (UserDataPtr)scopflow);

  /* Solve */
  PipsNlpSolveStruct(pips->nlp);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverDestroy_PIPS(SCOPFLOW scopflow)
{
  PetscErrorCode     ierr;
  SCOPFLOWSolver_PIPS pips=scopflow->solver;
  OPFLOWSolver_IPOPT  opflowipopt;
  OPFLOW              opflow;

  PetscFunctionBegin;
  PetscInt i;

  if(pips->nlp) {
    /* Free Pips problem struct */
    FreePipsNlpProblemStruct(pips->nlp);
    pips->nlp = NULL;
  }

  ierr = PetscFree(pips->xstarti);CHKERRQ(ierr);
  ierr = PetscFree(pips->gstarti);CHKERRQ(ierr);
  ierr = PetscFree(pips->nxi);CHKERRQ(ierr);
  ierr = PetscFree(pips->ngi);CHKERRQ(ierr);

  for(i=0; i < scopflow->Ns; i++) {
    opflow = scopflow->opflows[i];
    if(opflow->Nconineq + scopflow->nconineqcoup[i]) {
      opflowipopt = opflow->solver;
      if(!opflow->Nconineq) {
	ierr = MatDestroy(&opflow->Jac_Gi);CHKERRQ(ierr);
	ierr = MatDestroy(&opflowipopt->Jac_GiT);CHKERRQ(ierr);
	ierr = PetscFree(opflowipopt->jac_gi->colptr);CHKERRQ(ierr);
	ierr = PetscFree(opflowipopt->jac_gi->rowidx);CHKERRQ(ierr);
	ierr = PetscFree(opflowipopt->jac_gi->values);CHKERRQ(ierr);
	ierr = PetscFree(opflowipopt->jac_gi);CHKERRQ(ierr);
      }
    }
    ierr = MatDestroy(&pips->Jac_GicoupT[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree(pips->Jac_GicoupT);CHKERRQ(ierr);

  ierr = PetscFree(pips);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverSetUp_PIPS(SCOPFLOW scopflow)
{
  PetscErrorCode ierr;
  SCOPFLOWSolver_PIPS pips = (SCOPFLOWSolver_PIPS)scopflow->solver;
  PetscInt       i,j,k,ctr;
  OPFLOW         opflow;
  PetscInt       ngen;
  PetscScalar    *x,*xi,*xl,*xu,*xli,*xui,*gl,*gu,*gli,*gui;
  PS             ps;
  PSBUS          bus;
  PSGEN          gen;

  PetscFunctionBegin;

  ierr = PetscCalloc1(scopflow->Ns,&pips->xstarti);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->Ns,&pips->gstarti);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->Ns,&pips->nxi);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->Ns,&pips->ngi);CHKERRQ(ierr);

  scopflow->Nx = 0;
  scopflow->Ncon = 0;
  pips->xstarti[0] = 0;
  pips->gstarti[0] = 0;

  ierr = PSGetNumGenerators(scopflow->opflows[0]->ps,&ngen,NULL);CHKERRQ(ierr);
  for(i=0; i < scopflow->Ns; i++) {
    opflow = scopflow->opflows[i];
    pips->nxi[i] = opflow->nx;
    if(scopflow->iscoupling) scopflow->nconineqcoup[i] = (i == 0)?0:ngen;
    else scopflow->nconineqcoup[i] = 0;

    pips->ngi[i] = opflow->ncon + scopflow->nconineqcoup[i];
    if(i < scopflow->Ns - 1) {
      pips->xstarti[i+1] = pips->xstarti[i] + pips->nxi[i];
      pips->gstarti[i+1] = pips->gstarti[i] + pips->ngi[i];
    }
    scopflow->Nx += pips->nxi[i];
    scopflow->Ncon += pips->ngi[i];
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

  for(i=0; i < scopflow->Ns; i++) {
    opflow = scopflow->opflows[i];
    /* Set initial guess and bounds on variables */
    xi  = x  + pips->xstarti[i];
    xli = xl + pips->xstarti[i];
    xui = xu + pips->xstarti[i];

    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Xl,xli);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Xu,xui);CHKERRQ(ierr);

    /* Set bounds */
    ierr = (*opflow->formops.setvariablebounds)(opflow,opflow->Xl,opflow->Xu);CHKERRQ(ierr);

    /* Set initial guess */
    ierr = (*opflow->formops.setinitialguess)(opflow,opflow->X);CHKERRQ(ierr);

    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Xl);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Xu);CHKERRQ(ierr);

    /* Set bounds on constraints */
    gli = gl + pips->gstarti[i];
    gui = gu + pips->gstarti[i];
    
    ierr = VecPlaceArray(opflow->Gl,gli);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Gu,gui);CHKERRQ(ierr);

    ierr = (*opflow->formops.setconstraintbounds)(opflow,opflow->Gl,opflow->Gu);CHKERRQ(ierr);

    ierr = VecResetArray(opflow->Gl);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Gu);CHKERRQ(ierr);
    
    if(scopflow->nconineqcoup[i]) {
      ctr = 0;
      ps = opflow->ps;
      /* Bounds on coupling constraints */
      for(j=0; j < ps->nbus; j++) {
	bus = &ps->bus[j];

	for(k=0; k < bus->ngen; k++) {
	  ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
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

PetscErrorCode SCOPFLOWSolverCreate_PIPS(SCOPFLOW scopflow)
{
  PetscErrorCode ierr;
  SCOPFLOWSolver_PIPS pips;
  
  PetscFunctionBegin;

  ierr = PetscCalloc1(1,&pips);CHKERRQ(ierr);

  pips->nlp = NULL;
  pips->nnz_jac_g = 0;
  pips->nnz_hes = 0;
  scopflow->solver = pips;

  scopflow->solverops.setup = SCOPFLOWSolverSetUp_PIPS;
  scopflow->solverops.solve = SCOPFLOWSolverSolve_PIPS;
  scopflow->solverops.destroy = SCOPFLOWSolverDestroy_PIPS;

  PetscFunctionReturn(0);
}

#endif
