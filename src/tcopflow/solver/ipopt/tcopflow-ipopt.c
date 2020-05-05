#include <scopflow_config.h>
#if defined(SCOPFLOW_HAVE_IPOPT)

#include <private/opflowimpl.h>
#include <private/tcopflowimpl.h>
#include "tcopflow-ipopt.h"

/* IPOPT callback functions */
Bool eval_tcopflow_f(PetscInt n, PetscScalar* x, Bool new_x,
            PetscScalar* obj_value, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  TCOPFLOW  tcopflow=(TCOPFLOW)user_data;
  TCOPFLOWSolver_IPOPT tcopflowipopt = (TCOPFLOWSolver_IPOPT)tcopflow->solver;
  OPFLOW    opflow;
  PetscInt  i;
  PetscScalar *xi;
  PetscScalar opflowobj;
  PetscInt    k;

  *obj_value = 0.0;

  for(i=0; i < tcopflow->Ns; i++) {
    opflowobj = 0.0;
    xi = x + tcopflowipopt->xstarti[i];
    opflow = tcopflow->opflows[i];
    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
    ierr = (*opflow->formops.computeobjective)(opflow,opflow->X,&opflowobj);CHKERRQ(ierr);
    *obj_value += opflowobj;
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  }
				
  return TRUE;
}

Bool eval_tcopflow_grad_f(PetscInt n, PetscScalar* x, Bool new_x,
                 PetscScalar* grad_f, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  TCOPFLOW  tcopflow=(TCOPFLOW)user_data;
  TCOPFLOWSolver_IPOPT tcopflowipopt = (TCOPFLOWSolver_IPOPT)tcopflow->solver;
  OPFLOW    opflow;
  PetscInt  i;
  PetscScalar *xi,*gradi;

  for(i=0; i < tcopflow->Ns; i++) {
    opflow = tcopflow->opflows[i];
    xi = x + tcopflowipopt->xstarti[i];
    gradi = grad_f + tcopflowipopt->xstarti[i];


    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->gradobj,gradi);CHKERRQ(ierr);
    ierr = VecSet(opflow->gradobj,0.0);CHKERRQ(ierr);

    ierr = (*opflow->formops.computegradient)(opflow,opflow->X,opflow->gradobj);CHKERRQ(ierr);
    
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->gradobj);CHKERRQ(ierr);
  }
  return TRUE;
}

Bool eval_tcopflow_g(PetscInt n, PetscScalar* x, Bool new_x,
            PetscInt m, PetscScalar* g, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  TCOPFLOW  tcopflow=(TCOPFLOW)user_data;
  TCOPFLOWSolver_IPOPT tcopflowipopt = (TCOPFLOWSolver_IPOPT)tcopflow->solver;
  OPFLOW    opflowtpre,opflow;
  PetscInt  i,j,k,loc,loctpre,ctr;
  PetscScalar *xtpre,*xi,*gi;
  PS        ps,pstpre;
  PSBUS     bus,bustpre;
  PSGEN     gen,gentpre;

  for(i=0; i < tcopflow->Ns; i++) {
    xi   = x + tcopflowipopt->xstarti[i];
    gi   = g + tcopflowipopt->gstarti[i];

    opflow = tcopflow->opflows[i];

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

    if(tcopflow->nconineqcoup[i]) {
      ctr = 0;
      ps = opflow->ps;
      pstpre = opflowtpre->ps;
      for(j=0; j < ps->nbus; j++) {
	bus = &ps->bus[j];
	bustpre = &pstpre->bus[j];
	ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
	ierr = PSBUSGetVariableLocation(bustpre,&loctpre);CHKERRQ(ierr);
	
	for(k=0; k < bus->ngen; k++) {
	  ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	  ierr = PSBUSGetGen(bustpre,k,&gentpre);CHKERRQ(ierr);

	  if(!gen->status) {
	    if(gentpre->status) loctpre += 2;
	    continue;
	  } else {
	    loc += 2;
	    if(!gentpre->status) continue;
	    loctpre += 2;
	  }

	  gi[ctr] = xi[loc] - xtpre[loctpre]; /* PG(t) - PG(t-dT) */
	  ctr++;
	}
      }
    }
    
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
    opflowtpre = opflow;
    xtpre = xi;
  }

  return TRUE;
}

Bool eval_tcopflow_jac_g(PetscInt n, PetscScalar *x, Bool new_x,
                PetscInt m, PetscInt nele_jac,
                PetscInt *iRow, PetscInt *jCol, PetscScalar *values,
                UserDataPtr user_data)
{

  PetscErrorCode ierr;
  TCOPFLOW         tcopflow=(TCOPFLOW)user_data;
  TCOPFLOWSolver_IPOPT tcopflowipopt=(TCOPFLOWSolver_IPOPT)tcopflow->solver;
  OPFLOW         opflow,opflowtpre;
  OPFLOWSolver_IPOPT opflowipopt;
  PetscInt       *iRowstart = iRow,*jColstart=jCol;
  PetscInt       roffset,coffset;
  PetscInt       nrow,ncol;
  PetscScalar    *xi,*xarr;
  PetscInt       i,j,k,loc,loctpre,xtpreloc,xiloc;
  PS             ps,pstpre;
  PSBUS          bus,bustpre;
  PSGEN          gen,gentpre;
  PetscInt       nvals;
  const PetscInt *cols;
  const PetscScalar *vals;

  if(values == NULL) {
    /* Set locations only */
    ierr = VecGetArray(tcopflow->X,&xarr);CHKERRQ(ierr);
    for(i=0; i < tcopflow->Ns; i++) {
      opflow = tcopflow->opflows[i];
      opflowipopt = (OPFLOWSolver_IPOPT)opflow->solver;

      roffset = tcopflowipopt->gstarti[i];
      coffset = tcopflowipopt->xstarti[i];

      xi = xarr + tcopflowipopt->xstarti[i];
      ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);

      ierr = (*opflow->formops.computeequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Ge);CHKERRQ(ierr);

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
	ierr = (*opflow->formops.computeinequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Gi);CHKERRQ(ierr);

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

      if(tcopflow->nconineqcoup[i]) {
	ps = opflow->ps;
	pstpre = opflowtpre->ps;
	for(j=0; j < ps->nbus; j++) {
	  bus = &ps->bus[j];
	  bustpre = &pstpre->bus[j];
	  ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
	  ierr = PSBUSGetVariableLocation(bustpre,&loctpre);CHKERRQ(ierr);
	  for(k=0; k < bus->ngen; k++) {
	    ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	    ierr = PSBUSGetGen(bustpre,k,&gentpre);CHKERRQ(ierr);

	    if(!gen->status) {
	      if(gentpre->status) loctpre += 2;
	      continue;
	    } else {
	      loc += 2;
	      if(!gentpre->status) continue;
	      loctpre += 2;
	    }

	    xtpreloc = tcopflowipopt->xstarti[i-1] + loctpre;
	    xiloc = tcopflowipopt->xstarti[i] + loc;
	    iRowstart[0] = roffset;
	    jColstart[0] = xtpreloc;
	    iRowstart[1] = roffset;
	    jColstart[1] = xiloc;
	    iRowstart += 2;
	    jColstart += 2;
	    roffset += 1;
	  }
	}
      }

      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
      opflowtpre = opflow;
    }
    ierr = VecRestoreArray(tcopflow->X,&xarr);CHKERRQ(ierr);
  } else {

    for(i=0; i < tcopflow->Ns; i++) {
      opflow = tcopflow->opflows[i];
      opflowipopt = (OPFLOWSolver_IPOPT)opflow->solver;

      xi = x + tcopflowipopt->xstarti[i];
      ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
      /* Compute equality constraint jacobian */
      ierr = (*opflow->formops.computeequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Ge);CHKERRQ(ierr);

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
	ierr = (*opflow->formops.computeinequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Gi);CHKERRQ(ierr);

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

      if(tcopflow->nconineqcoup[i]) {
	ps = opflow->ps;
	pstpre = opflowtpre->ps;
	for(j=0; j < ps->nbus; j++) {
	  bus = &ps->bus[j];
	  bustpre = &pstpre->bus[j];
	  ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
	  ierr = PSBUSGetVariableLocation(bustpre,&loctpre);CHKERRQ(ierr);
	  for(k=0; k < bus->ngen; k++) {
	    ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	    ierr = PSBUSGetGen(bustpre,k,&gentpre);CHKERRQ(ierr);
	    if(!gen->status) {
	      if(gentpre->status) loctpre += 2;
	      continue;
	    } else {
	      loc += 2;
	      if(!gentpre->status) continue;
	      loctpre += 2;
	    }

	    xtpreloc = tcopflowipopt->xstarti[i-1] + loctpre;
	    xiloc = tcopflowipopt->xstarti[i] + loc;
	    values[0] = -1;
	    values[1] = 1;
	    values += 2;
	  }
	}
      }

      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
      opflowtpre = opflow;
    }
  }

  return TRUE;
}

Bool eval_tcopflow_h(PetscInt n, PetscScalar *x, Bool new_x, PetscScalar obj_factor,
            PetscInt m, PetscScalar *lambda, Bool new_lambda,
            PetscInt nele_hess, PetscInt *iRow, PetscInt *jCol,
            PetscScalar *values, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  PetscInt       nrow;
  TCOPFLOW         tcopflow=(TCOPFLOW)user_data;
  TCOPFLOWSolver_IPOPT tcopflowipopt=(TCOPFLOWSolver_IPOPT)tcopflow->solver;
  OPFLOW          opflow;
  OPFLOWSolver_IPOPT opflowipopt;
  PetscScalar     *xi,*valuesi=values,*lameqi,*lamineqi;
  PetscInt         i;
  PetscInt        *iRowStart=iRow,*jColStart=jCol;
  PetscInt        roffset;
  PetscInt       nvals;
  const PetscInt *cols;
  const PetscScalar *vals;
  PetscInt j,k;
  PetscInt ctr=0;

  tcopflow->obj_factor = obj_factor;

  if(values == NULL) {

    for(i=0; i < tcopflow->Ns; i++) {
      opflow = tcopflow->opflows[i];
      opflowipopt = (OPFLOWSolver_IPOPT)opflow->solver;

      roffset = tcopflowipopt->xstarti[i];

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

    for(i=0; i < tcopflow->Ns; i++) {
      opflow = tcopflow->opflows[i];
      opflowipopt = (OPFLOWSolver_IPOPT)opflow->solver;
      opflow->obj_factor = obj_factor;

      roffset = tcopflowipopt->xstarti[i];

      xi = x + roffset;
      lameqi = lambda + tcopflowipopt->gstarti[i];

      ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
      ierr = VecPlaceArray(opflow->Lambdae,lameqi);CHKERRQ(ierr);
      if(opflow->Nconineq) {
	lamineqi = lameqi + opflow->nconeq;
	ierr = VecPlaceArray(opflow->Lambdai,lamineqi);CHKERRQ(ierr);
      }

      ierr = (*opflow->formops.computehessian)(opflow,opflow->X,opflow->Lambdae,opflow->Lambdai,opflow->Hes);CHKERRQ(ierr);

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

PetscErrorCode TCOPFLOWSolverSolve_IPOPT(TCOPFLOW tcopflow)
{
  PetscErrorCode     ierr;
  TCOPFLOWSolver_IPOPT tcopflowipopt = (TCOPFLOWSolver_IPOPT)tcopflow->solver;
  OPFLOW             opflow;
  OPFLOWSolver_IPOPT opflowipopt;
  PetscScalar        *x,*xl,*xu,*gl,*gu,*xi,*lameqi,*lamineqi,*lam;
  PetscInt           i;
  MatInfo            info_eq,info_ineq,info_hes;

  PetscFunctionBegin;

  tcopflowipopt->nnz_jac_ge = tcopflowipopt->nnz_jac_gi = 0;

  ierr = VecGetArray(tcopflow->X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(tcopflow->Lambda,&lam);CHKERRQ(ierr);

  for(i=0; i < tcopflow->Ns; i++) {
    opflow = tcopflow->opflows[i];
    opflowipopt = (OPFLOWSolver_IPOPT)opflow->solver;
    xi = x + tcopflowipopt->xstarti[i];

    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);

    /* Compute nonzeros for the Jacobian */
    /* Equality constraint Jacobian */
    ierr = (*opflow->formops.computeequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Ge);CHKERRQ(ierr);
    ierr = MatGetInfo(opflow->Jac_Ge,MAT_LOCAL,&info_eq);CHKERRQ(ierr);

    opflowipopt->nnz_jac_ge = info_eq.nz_used;
    tcopflowipopt->nnz_jac_ge += opflowipopt->nnz_jac_ge;

    opflowipopt->nnz_jac_gi = 0;
    if(opflow->Nconineq) {
      ierr = (*opflow->formops.computeinequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Gi);CHKERRQ(ierr);
      ierr = MatGetInfo(opflow->Jac_Gi,MAT_LOCAL,&info_ineq);CHKERRQ(ierr);

      opflowipopt->nnz_jac_gi = info_ineq.nz_used;
      
      tcopflowipopt->nnz_jac_gi += opflowipopt->nnz_jac_gi;
    }   
    opflowipopt->nnz_jac_g = opflowipopt->nnz_jac_ge + opflowipopt->nnz_jac_gi;

    /* Add non-zeros for Jacobian of coupling constraints */
    if(tcopflow->nconineqcoup[i]) tcopflowipopt->nnz_jac_gi += 2*tcopflow->nconineqcoup[i];


    /* Compute non-zeros for Hessian */

    lameqi = lam + tcopflowipopt->gstarti[i];
    ierr = VecPlaceArray(opflow->Lambdae,lameqi);CHKERRQ(ierr);

    if(opflow->Nconineq) {
      lamineqi = lam + tcopflowipopt->gstarti[i]+opflow->nconeq;
      ierr = VecPlaceArray(opflow->Lambdai,lamineqi);CHKERRQ(ierr);
    } else {
      opflow->Lambdai = NULL;
    }

    ierr = (*opflow->formops.computehessian)(opflow,opflow->X,opflow->Lambdae,opflow->Lambdai,opflow->Hes);CHKERRQ(ierr);
    ierr = MatGetInfo(opflow->Hes,MAT_LOCAL,&info_hes);CHKERRQ(ierr);

    opflowipopt->nnz_hes = (info_hes.nz_used  -opflow->nx)/2 + opflow->nx;

    tcopflowipopt->nnz_hes += opflowipopt->nnz_hes;

    ierr = VecResetArray(opflow->Lambdae);CHKERRQ(ierr);
    if(opflow->Nconineq) {
      ierr = VecResetArray(opflow->Lambdai);CHKERRQ(ierr);
    }
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(tcopflow->X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(tcopflow->Lambda,&lam);CHKERRQ(ierr);

  tcopflowipopt->nnz_jac_g = tcopflowipopt->nnz_jac_ge + tcopflowipopt->nnz_jac_gi;

  /* Create IPOPT solver instance */
  ierr = VecGetArray(tcopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(tcopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecGetArray(tcopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(tcopflow->Gu,&gu);CHKERRQ(ierr);

  /* Create IPOPT problem */
  tcopflowipopt->nlp = CreateIpoptProblem(tcopflow->Nx,xl,xu,tcopflow->Ncon,gl,gu,tcopflowipopt->nnz_jac_g,tcopflowipopt->nnz_hes,0,&eval_tcopflow_f,
				   &eval_tcopflow_g,&eval_tcopflow_grad_f,
				   &eval_tcopflow_jac_g,&eval_tcopflow_h);

  ierr = VecRestoreArray(tcopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(tcopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecRestoreArray(tcopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(tcopflow->Gu,&gu);CHKERRQ(ierr);

  ierr = VecGetArray(tcopflow->X,&x);CHKERRQ(ierr);

  /* Solve */
  tcopflowipopt->solve_status = IpoptSolve(tcopflowipopt->nlp,x,NULL,&tcopflow->obj,NULL,NULL,NULL,tcopflow);

  ierr = VecRestoreArray(tcopflow->X,&x);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode TCOPFLOWSolverDestroy_IPOPT(TCOPFLOW tcopflow)
{
  PetscErrorCode     ierr;
  TCOPFLOWSolver_IPOPT ipopt = (TCOPFLOWSolver_IPOPT)tcopflow->solver;

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

PetscErrorCode TCOPFLOWSolverSetUp_IPOPT(TCOPFLOW tcopflow)
{
  PetscErrorCode ierr;
  TCOPFLOWSolver_IPOPT ipopt = (TCOPFLOWSolver_IPOPT)tcopflow->solver;
  PetscInt       i,j,k,ctr;
  OPFLOW         opflow,opflow0;
  PetscInt       ngenON;
  PetscScalar    *x,*xi,*xl,*xu,*xli,*xui,*gl,*gu,*gli,*gui;
  PS             ps,ps0;
  PSBUS          bus,bus0;
  PSGEN          gen,gen0;

  PetscFunctionBegin;

  ierr = PetscCalloc1(tcopflow->Ns,&ipopt->xstarti);CHKERRQ(ierr);
  ierr = PetscCalloc1(tcopflow->Ns,&ipopt->gstarti);CHKERRQ(ierr);
  ierr = PetscCalloc1(tcopflow->Ns,&ipopt->nxi);CHKERRQ(ierr);
  ierr = PetscCalloc1(tcopflow->Ns,&ipopt->ngi);CHKERRQ(ierr);

  tcopflow->Nx = 0;
  tcopflow->Ncon = 0;
  ipopt->xstarti[0] = 0;
  ipopt->gstarti[0] = 0;

  for(i=0; i < tcopflow->Ns; i++) {
    opflow = tcopflow->opflows[i];
    ierr = PSGetNumActiveGenerators(opflow->ps,&ngenON,NULL);CHKERRQ(ierr);
    ipopt->nxi[i] = opflow->nx;
    if(tcopflow->iscoupling) tcopflow->nconineqcoup[i] = (i == 0)?0:ngenON;
    else tcopflow->nconineqcoup[i] = 0;

    ipopt->ngi[i] = opflow->ncon + tcopflow->nconineqcoup[i];
    if(i < tcopflow->Ns - 1) {
      ipopt->xstarti[i+1] = ipopt->xstarti[i] + ipopt->nxi[i];
      ipopt->gstarti[i+1] = ipopt->gstarti[i] + ipopt->ngi[i];
    }
    tcopflow->Nx += ipopt->nxi[i];
    tcopflow->Ncon += ipopt->ngi[i];
  }

  /* Create vector X */
  ierr = VecCreate(tcopflow->comm->type,&tcopflow->X);CHKERRQ(ierr);
  ierr = VecSetSizes(tcopflow->X,tcopflow->Nx,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(tcopflow->X);CHKERRQ(ierr);

  ierr = VecDuplicate(tcopflow->X,&tcopflow->Xl);CHKERRQ(ierr);
  ierr = VecDuplicate(tcopflow->X,&tcopflow->Xu);CHKERRQ(ierr);

  /* vector for constraints */
  ierr = VecCreate(tcopflow->comm->type,&tcopflow->G);CHKERRQ(ierr);
  ierr = VecSetSizes(tcopflow->G,tcopflow->Ncon,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(tcopflow->G);CHKERRQ(ierr);

  /* Constraint bounds vectors  */
  ierr = VecDuplicate(tcopflow->G,&tcopflow->Gl);CHKERRQ(ierr);
  ierr = VecDuplicate(tcopflow->G,&tcopflow->Gu);CHKERRQ(ierr);

  /* Lagrangian multipliers */
  ierr = VecDuplicate(tcopflow->G,&tcopflow->Lambda);CHKERRQ(ierr);

  /* Set Initial guess and Bounds on variables and constraints */
  ierr = VecGetArray(tcopflow->X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(tcopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(tcopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecGetArray(tcopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(tcopflow->Gu,&gu);CHKERRQ(ierr);

  opflow0 = tcopflow->opflows[0];
  for(i=0; i < tcopflow->Ns; i++) {
    opflow = tcopflow->opflows[i];
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
    
    if(tcopflow->nconineqcoup[i]) {
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
	  gli[opflow->ncon + ctr] = -gen->ramp_rate_min*tcopflow->dT;
	  gui[opflow->ncon + ctr] =  gen->ramp_rate_min*tcopflow->dT;
	  ctr++;
	}
      }
    }
  }
  
  ierr = VecRestoreArray(tcopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(tcopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecRestoreArray(tcopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(tcopflow->Gu,&gu);CHKERRQ(ierr);

  /* Initialize Lagrange multiplier */
  ierr = VecSet(tcopflow->Lambda,1.0);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

PetscErrorCode TCOPFLOWSolverCreate_IPOPT(TCOPFLOW tcopflow)
{
  PetscErrorCode ierr;
  TCOPFLOWSolver_IPOPT ipopt;
  
  PetscFunctionBegin;

  if(tcopflow->comm->size > 1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"IPOPT solver does not support execution in parallel\n",tcopflow->comm->size); 
  ierr = PetscCalloc1(1,&ipopt);CHKERRQ(ierr);

  ipopt->nlp = NULL;
  ipopt->nnz_jac_g = 0;
  ipopt->nnz_hes = 0;
  tcopflow->solver = ipopt;

  tcopflow->solverops.setup = TCOPFLOWSolverSetUp_IPOPT;
  tcopflow->solverops.solve = TCOPFLOWSolverSolve_IPOPT;
  tcopflow->solverops.destroy = TCOPFLOWSolverDestroy_IPOPT;

  PetscFunctionReturn(0);
}

#endif
