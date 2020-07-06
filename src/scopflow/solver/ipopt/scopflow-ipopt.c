#include <scopflow_config.h>
#if defined(SCOPFLOW_HAVE_IPOPT)

#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>
#include "scopflow-ipopt.h"

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

  *obj_value = 0.0;

  for(i=0; i < scopflow->Nc; i++) {
    opflowobj = 0.0;
    xi = x + scopflowipopt->xstarti[i];
    opflow = scopflow->opflows[i];
    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
    ierr = (*opflow->modelops.computeobjective)(opflow,opflow->X,&opflowobj);CHKERRQ(ierr);
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

  for(i=0; i < scopflow->Nc; i++) {
    opflow = scopflow->opflows[i];
    xi = x + scopflowipopt->xstarti[i];
    gradi = grad_f + scopflowipopt->xstarti[i];


    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->gradobj,gradi);CHKERRQ(ierr);
    ierr = VecSet(opflow->gradobj,0.0);CHKERRQ(ierr);

    ierr = (*opflow->modelops.computegradient)(opflow,opflow->X,opflow->gradobj);CHKERRQ(ierr);
    
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
  for(i=0; i < scopflow->Nc; i++) {
    xi   = x + scopflowipopt->xstarti[i];
    gi   = g + scopflowipopt->gstarti[i];

    opflow = scopflow->opflows[i];

    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);

    /* Equality constraints */
    ierr = VecPlaceArray(opflow->Ge,gi);CHKERRQ(ierr);
    ierr = (*opflow->modelops.computeequalityconstraints)(opflow,opflow->X,opflow->Ge);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Ge);CHKERRQ(ierr);
    gi = gi + opflow->nconeq;
      
    if(opflow->Nconineq) {
      /* Inequality constraints */
      ierr = VecPlaceArray(opflow->Gi,gi);CHKERRQ(ierr);
      ierr = (*opflow->modelops.computeinequalityconstraints)(opflow,opflow->X,opflow->Gi);CHKERRQ(ierr);
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
  PetscInt       nrow,ncol;
  PetscScalar    *xi,*xarr;
  PetscInt       i,j,k,loc,loc0,x0loc,xiloc;
  PS             ps,ps0;
  PSBUS          bus,bus0;
  PSGEN          gen,gen0;
  PetscInt       nvals;
  const PetscInt *cols;
  const PetscScalar *vals;

  if(values == NULL) {
    /* Set locations only */
    ierr = VecGetArray(scopflow->X,&xarr);CHKERRQ(ierr);
    opflow0 = scopflow->opflows[0];
    for(i=0; i < scopflow->Nc; i++) {
      opflow = scopflow->opflows[i];
      opflowipopt = (OPFLOWSolver_IPOPT)opflow->solver;

      roffset = scopflowipopt->gstarti[i];
      coffset = scopflowipopt->xstarti[i];

      xi = xarr + scopflowipopt->xstarti[i];
      ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);

      ierr = (*opflow->modelops.computeequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Ge);CHKERRQ(ierr);

      ierr = MatGetSize(opflow->Jac_Ge,&nrow,&ncol);CHKERRQ(ierr);
      /* Copy over locations to triplet format */
      for(j=0; j < nrow; j++) {
	ierr = MatGetRow(opflow->Jac_Ge,j,&nvals,&cols,&vals);CHKERRQ(ierr);
	for(k=0; k < nvals; k++) {
	  iRowstart[k] = roffset + j;
	  jColstart[k] = coffset + cols[k];
	}
	/* Increment iRow,jCol pointers */
	iRowstart += nvals;
	jColstart += nvals;
	ierr = MatRestoreRow(opflow->Jac_Ge,j,&nvals,&cols,&vals);CHKERRQ(ierr);
      }

      roffset += opflow->nconeq;
      if(opflow->Nconineq) {
	/* Inequality constrained Jacobian */
	ierr = (*opflow->modelops.computeinequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Gi);CHKERRQ(ierr);

	ierr = MatGetSize(opflow->Jac_Gi,&nrow,&ncol);CHKERRQ(ierr);
	/* Copy over locations to triplet format */
	for(j=0; j < nrow; j++) {
	  ierr = MatGetRow(opflow->Jac_Gi,j,&nvals,&cols,&vals);CHKERRQ(ierr);
	  for(k=0; k < nvals; k++) {
	    iRowstart[k] = roffset + j;
	    jColstart[k] = coffset + cols[k];
	  }
	  /* Increment iRow,jCol pointers */
	  iRowstart += nvals;
	  jColstart += nvals;
	  ierr = MatRestoreRow(opflow->Jac_Gi,j,&nvals,&cols,&vals);CHKERRQ(ierr);
	}

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
    for(i=0; i < scopflow->Nc; i++) {
      opflow = scopflow->opflows[i];
      opflowipopt = (OPFLOWSolver_IPOPT)opflow->solver;

      xi = x + scopflowipopt->xstarti[i];
      ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
      /* Compute equality constraint jacobian */
      ierr = (*opflow->modelops.computeequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Ge);CHKERRQ(ierr);

      ierr = MatGetSize(opflow->Jac_Ge,&nrow,&ncol);CHKERRQ(ierr);
      /* Copy over values */
      for(j=0; j < nrow; j++) {
	ierr = MatGetRow(opflow->Jac_Ge,j,&nvals,&cols,&vals);CHKERRQ(ierr);
	for(k=0; k < nvals; k++) {
	  values[k] = vals[k];
	}
	values += nvals;
	ierr = MatRestoreRow(opflow->Jac_Ge,j,&nvals,&cols,&vals);CHKERRQ(ierr);
      }

      if(opflow->Nconineq) {
	/* Compute inequality constraint jacobian */
	ierr = (*opflow->modelops.computeinequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Gi);CHKERRQ(ierr);

	ierr = MatGetSize(opflow->Jac_Gi,&nrow,&ncol);CHKERRQ(ierr);
	/* Copy over values */
	for(j=0; j < nrow; j++) {
	  ierr = MatGetRow(opflow->Jac_Gi,j,&nvals,&cols,&vals);CHKERRQ(ierr);
	  for(k=0; k < nvals; k++) {
	    values[k] = vals[k];
	  }
	  values += nvals;
	  ierr = MatRestoreRow(opflow->Jac_Gi,j,&nvals,&cols,&vals);CHKERRQ(ierr);
	}

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
	    values[0] = -1;
	    values[1] = 1;
	    values += 2;
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
  PetscErrorCode       ierr;
  PetscInt             nrow;
  SCOPFLOW             scopflow=(SCOPFLOW)user_data;
  SCOPFLOWSolver_IPOPT scopflowipopt=(SCOPFLOWSolver_IPOPT)scopflow->solver;
  OPFLOW               opflow;
  OPFLOWSolver_IPOPT   opflowipopt;
  PetscScalar          *xi,*lameqi,*lamineqi;
  PetscInt             i;
  PetscInt             roffset;
  PetscInt             nvals;
  const PetscInt       *cols;
  const PetscScalar    *vals;
  PetscInt             j,k;
  PetscInt             ctr=0;

  scopflow->obj_factor = obj_factor;

  if(values == NULL) {

    for(i=0; i < scopflow->Nc; i++) {
      opflow = scopflow->opflows[i];
      opflowipopt = (OPFLOWSolver_IPOPT)opflow->solver;

      roffset = scopflowipopt->xstarti[i];

      ierr = MatGetSize(opflow->Hes,&nrow,NULL);CHKERRQ(ierr);

      /* Copy over locations to triplet format */
      /* Note that IPOPT
	 requires a lower diagonal Hessian (see note https://coin-or.github.io/Ipopt/INTERFACES.html#INTERFACE_CODE)
	 Hence, we only add lower diagonal locations
      */
      for(j=0; j < nrow; j++) {
	ierr = MatGetRow(opflow->Hes,j,&nvals,&cols,&vals);CHKERRQ(ierr);
	ctr = 0;
	for(k=0; k < nvals; k++) {
	  if(cols[k] >= j) { /* upper triangle */
	    /* save as lower triangle locations */
	    iRow[ctr] = roffset + cols[k];
	    jCol[ctr] = roffset + j;
	    ctr++;
	  }
	}
	iRow += ctr;
	jCol += ctr;
	ierr = MatRestoreRow(opflow->Hes,j,&nvals,&cols,&vals);CHKERRQ(ierr);
      }

    }
  } else {

    for(i=0; i < scopflow->Nc; i++) {
      opflow = scopflow->opflows[i];
      opflowipopt = (OPFLOWSolver_IPOPT)opflow->solver;
      opflow->obj_factor = obj_factor;

      roffset = scopflowipopt->xstarti[i];

      xi = x + roffset;
      lameqi = lambda + scopflowipopt->gstarti[i];

      ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
      ierr = VecPlaceArray(opflow->Lambdae,lameqi);CHKERRQ(ierr);
      if(opflow->Nconineq) {
	lamineqi = lameqi + opflow->nconeq;
	ierr = VecPlaceArray(opflow->Lambdai,lamineqi);CHKERRQ(ierr);
      }

      ierr = (*opflow->modelops.computehessian)(opflow,opflow->X,opflow->Lambdae,opflow->Lambdai,opflow->Hes);CHKERRQ(ierr);

      /* Copy over values */
      ierr = MatGetSize(opflow->Hes,&nrow,&nrow);CHKERRQ(ierr);
      for(j=0; j < nrow; j++) {
	ierr = MatGetRow(opflow->Hes,j,&nvals,&cols,&vals);CHKERRQ(ierr);
	ctr = 0;
	for(k=0; k < nvals; k++) {
	  if(cols[k] >= j) { /* Upper triangle values (same as lower triangle) */
	    values[ctr] = vals[k];
	    ctr++;
	  }
	}
	values += ctr;
	ierr = MatRestoreRow(opflow->Hes,j,&nvals,&cols,&vals);CHKERRQ(ierr);
      }

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
  PetscScalar        *x,*xl,*xu,*g,*gl,*gu,*xi,*lameqi,*lamineqi,*lam;
  PetscInt           i;
  MatInfo            info_eq,info_ineq,info_hes;

  PetscFunctionBegin;

  scopflowipopt->nnz_jac_ge = scopflowipopt->nnz_jac_gi = 0;

  ierr = VecGetArray(scopflow->X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Lambda,&lam);CHKERRQ(ierr);

  for(i=0; i < scopflow->Nc; i++) {
    opflow = scopflow->opflows[i];
    opflowipopt = (OPFLOWSolver_IPOPT)opflow->solver;
    xi = x + scopflowipopt->xstarti[i];

    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);

    /* Compute nonzeros for the Jacobian */
    /* Equality constraint Jacobian */
    ierr = (*opflow->modelops.computeequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Ge);CHKERRQ(ierr);
    ierr = MatGetInfo(opflow->Jac_Ge,MAT_LOCAL,&info_eq);CHKERRQ(ierr);

    opflowipopt->nnz_jac_ge = info_eq.nz_used;
    scopflowipopt->nnz_jac_ge += opflowipopt->nnz_jac_ge;

    opflowipopt->nnz_jac_gi = 0;
    if(opflow->Nconineq) {
      ierr = (*opflow->modelops.computeinequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Gi);CHKERRQ(ierr);
      ierr = MatGetInfo(opflow->Jac_Gi,MAT_LOCAL,&info_ineq);CHKERRQ(ierr);

      opflowipopt->nnz_jac_gi = info_ineq.nz_used;
      
      scopflowipopt->nnz_jac_gi += opflowipopt->nnz_jac_gi;
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

    ierr = (*opflow->modelops.computehessian)(opflow,opflow->X,opflow->Lambdae,opflow->Lambdai,opflow->Hes);CHKERRQ(ierr);
    ierr = MatGetInfo(opflow->Hes,MAT_LOCAL,&info_hes);CHKERRQ(ierr);

    opflowipopt->nnz_hes = (info_hes.nz_used  -opflow->nx)/2 + opflow->nx;

    scopflowipopt->nnz_hes += opflowipopt->nnz_hes;

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
  ierr = VecGetArray(scopflow->G,&g);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Lambda,&lam);CHKERRQ(ierr);

  /* Solve */
  scopflowipopt->solve_status = IpoptSolve(scopflowipopt->nlp,x,g,&scopflow->obj,lam,NULL,NULL,scopflow);

  ierr = VecRestoreArray(scopflow->X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->G,&g);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Lambda,&lam);CHKERRQ(ierr);

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

PetscErrorCode SCOPFLOWSolverGetObjective_IPOPT(SCOPFLOW scopflow,PetscReal *obj)
{
  PetscFunctionBegin;
  *obj = scopflow->obj;
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetSolution_IPOPT(SCOPFLOW scopflow,PetscInt cont_num,Vec *X)
{
  PetscErrorCode ierr;
  SCOPFLOWSolver_IPOPT scopflowipopt = (SCOPFLOWSolver_IPOPT)scopflow->solver;
  OPFLOW         opflow=scopflow->opflows[cont_num];
  Vec            Xi=opflow->X;
  PetscInt       nxi=opflow->nx;
  PetscScalar    *xi,*x;
  PetscInt       ix=scopflowipopt->xstarti[cont_num];

  PetscFunctionBegin;
  ierr = VecGetArray(Xi,&xi);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->X,&x);CHKERRQ(ierr);

  ierr = PetscArraycpy(xi,x+ix,nxi);CHKERRQ(ierr);

  ierr = VecRestoreArray(Xi,&xi);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->X,&x);CHKERRQ(ierr);

  *X = Xi;
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConstraints_IPOPT(SCOPFLOW scopflow,PetscInt cont_num,Vec *G)
{
  PetscErrorCode ierr;
  SCOPFLOWSolver_IPOPT scopflowipopt = (SCOPFLOWSolver_IPOPT)scopflow->solver;
  OPFLOW         opflow=scopflow->opflows[cont_num];
  Vec            Gi=opflow->G;
  PetscInt       ngi=opflow->ncon;
  PetscScalar    *gi,*g;
  PetscInt       ig=scopflowipopt->gstarti[cont_num];

  PetscFunctionBegin;
  ierr = VecGetArray(Gi,&gi);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->G,&g);CHKERRQ(ierr);

  ierr = PetscArraycpy(gi,g+ig,ngi);CHKERRQ(ierr);

  ierr = VecRestoreArray(Gi,&gi);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->G,&g);CHKERRQ(ierr);

  *G = Gi;

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConstraintMultipliers_IPOPT(SCOPFLOW scopflow,PetscInt cont_num,Vec *Lambda)
{
  PetscErrorCode ierr;
  SCOPFLOWSolver_IPOPT scopflowipopt = (SCOPFLOWSolver_IPOPT)scopflow->solver;
  OPFLOW         opflow=scopflow->opflows[cont_num];
  Vec            Lambdai=opflow->Lambda;
  PetscInt       ngi=opflow->ncon;
  PetscScalar    *lambdai,*lambda;
  PetscInt       ig=scopflowipopt->gstarti[cont_num];

  PetscFunctionBegin;
  ierr = VecGetArray(Lambdai,&lambdai);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Lambda,&lambda);CHKERRQ(ierr);

  ierr = PetscArraycpy(lambdai,lambda+ig,ngi);CHKERRQ(ierr);

  ierr = VecRestoreArray(Lambdai,&lambdai);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Lambda,&lambda);CHKERRQ(ierr);

  *Lambda = Lambdai;

  PetscFunctionReturn(0);
}


PetscErrorCode SCOPFLOWSolverGetConvergenceStatus_IPOPT(SCOPFLOW scopflow,PetscBool *status)
{
  SCOPFLOWSolver_IPOPT ipopt = (SCOPFLOWSolver_IPOPT)scopflow->solver;

  PetscFunctionBegin;
  if(ipopt->solve_status < 2) *status = PETSC_TRUE; /* See IpReturnCodes_inc.h in IPOPT. The first two denote convergence */
  else *status = PETSC_FALSE;

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

  ierr = PetscCalloc1(scopflow->Nc,&ipopt->xstarti);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->Nc,&ipopt->gstarti);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->Nc,&ipopt->nxi);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->Nc,&ipopt->ngi);CHKERRQ(ierr);

  scopflow->Nx = 0;
  scopflow->Ncon = 0;
  scopflow->Nconeq = 0;
  scopflow->Nconineq = 0;
  scopflow->Nconcoup = 0;
  ipopt->xstarti[0] = 0;
  ipopt->gstarti[0] = 0;

  for(i=0; i < scopflow->Nc; i++) {
    opflow = scopflow->opflows[i];
    ierr = PSGetNumActiveGenerators(opflow->ps,&ngenON,NULL);CHKERRQ(ierr);
    ipopt->nxi[i] = opflow->nx;
    if(scopflow->iscoupling) scopflow->nconineqcoup[i] = (i == 0)?0:ngenON;
    else scopflow->nconineqcoup[i] = 0;

    ipopt->ngi[i] = opflow->ncon + scopflow->nconineqcoup[i];
    if(i < scopflow->Nc - 1) {
      ipopt->xstarti[i+1] = ipopt->xstarti[i] + ipopt->nxi[i];
      ipopt->gstarti[i+1] = ipopt->gstarti[i] + ipopt->ngi[i];
    }
    scopflow->Nx +=        ipopt->nxi[i];
    scopflow->Ncon +=      ipopt->ngi[i];
    scopflow->Nconeq +=    opflow->nconeq;
    scopflow->Nconineq += opflow->nconineq;
    scopflow->Nconcoup +=  scopflow->nconineqcoup[i];
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
  for(i=0; i < scopflow->Nc; i++) {
    opflow = scopflow->opflows[i];
    /* Set initial guess and bounds on variables */
    xi  = x  + ipopt->xstarti[i];
    xli = xl + ipopt->xstarti[i];
    xui = xu + ipopt->xstarti[i];

    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Xl,xli);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Xu,xui);CHKERRQ(ierr);

    /* Set bounds */
    ierr = (*opflow->modelops.setvariablebounds)(opflow,opflow->Xl,opflow->Xu);CHKERRQ(ierr);

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

    ierr = (*opflow->modelops.setconstraintbounds)(opflow,opflow->Gl,opflow->Gu);CHKERRQ(ierr);

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

	  if(scopflow->makeup_power_source == 0) {
	    /* Only ref. bus responsible for make-up power for contingencies */
	    if(bus->ide == REF_BUS) {
	      gli[opflow->ncon + ctr] = -10000;
	      gui[opflow->ncon + ctr] =  10000;
	    } else {
	      gli[opflow->ncon + ctr] = 0.0;
	      gui[opflow->ncon + ctr] = 0.0;
	    }	    
	  } else {
	    gli[opflow->ncon + ctr] = -scopflow->mode*gen->ramp_rate_30min;
	    gui[opflow->ncon + ctr] =  scopflow->mode*gen->ramp_rate_30min;
	  }
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
  scopflow->solverops.getobjective = SCOPFLOWSolverGetObjective_IPOPT;
  scopflow->solverops.getsolution  = SCOPFLOWSolverGetSolution_IPOPT;
  scopflow->solverops.getconvergencestatus = SCOPFLOWSolverGetConvergenceStatus_IPOPT;
  scopflow->solverops.getconstraints = SCOPFLOWSolverGetConstraints_IPOPT;
  scopflow->solverops.getconstraintmultipliers = SCOPFLOWSolverGetConstraintMultipliers_IPOPT;

  PetscFunctionReturn(0);
}

#endif
