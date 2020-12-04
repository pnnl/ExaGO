#include <exago_config.h>
#if defined(EXAGO_HAVE_IPOPT)

#include <private/opflowimpl.h>
#include <private/sopflowimpl.h>
#include "sopflow-ipopt.h"

/* IPOPT callback functions */
Bool eval_sopflow_f(PetscInt n, PetscScalar* x, Bool new_x,
            PetscScalar* obj_value, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  SOPFLOW  sopflow=(SOPFLOW)user_data;
  SOPFLOWSolver_IPOPT sopflowipopt = (SOPFLOWSolver_IPOPT)sopflow->solver;
  OPFLOW    opflow;
  PetscInt  i;
  PetscScalar *xi;
  PetscScalar opflowobj;

  *obj_value = 0.0;

  for(i=0; i < sopflow->ns; i++) {
    opflowobj = 0.0;
    xi = x + sopflowipopt->xstarti[i];
    opflow = sopflow->opflows[i];
    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
    ierr = (*opflow->modelops.computeobjective)(opflow,opflow->X,&opflowobj);CHKERRQ(ierr);
    *obj_value += opflowobj;
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  }
				
  return TRUE;
}

Bool eval_sopflow_grad_f(PetscInt n, PetscScalar* x, Bool new_x,
                 PetscScalar* grad_f, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  SOPFLOW  sopflow=(SOPFLOW)user_data;
  SOPFLOWSolver_IPOPT sopflowipopt = (SOPFLOWSolver_IPOPT)sopflow->solver;
  OPFLOW    opflow;
  PetscInt  i;
  PetscScalar *xi,*gradi;

  for(i=0; i < sopflow->ns; i++) {
    opflow = sopflow->opflows[i];
    xi = x + sopflowipopt->xstarti[i];
    gradi = grad_f + sopflowipopt->xstarti[i];


    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->gradobj,gradi);CHKERRQ(ierr);
    ierr = VecSet(opflow->gradobj,0.0);CHKERRQ(ierr);

    ierr = (*opflow->modelops.computegradient)(opflow,opflow->X,opflow->gradobj);CHKERRQ(ierr);
    
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->gradobj);CHKERRQ(ierr);
  }
  return TRUE;
}

Bool eval_sopflow_g(PetscInt n, PetscScalar* x, Bool new_x,
            PetscInt m, PetscScalar* g, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  SOPFLOW  sopflow=(SOPFLOW)user_data;
  SOPFLOWSolver_IPOPT sopflowipopt = (SOPFLOWSolver_IPOPT)sopflow->solver;
  OPFLOW    opflow0,opflow;
  PetscInt  i,j,k,loc,loc0,ctr;
  PetscScalar *x0,*xi,*gi;
  PS        ps,ps0;
  PSBUS     bus,bus0;
  PSGEN     gen,gen0;

  x0 = x;

  opflow0 = sopflow->opflows[0];
  for(i=0; i < sopflow->ns; i++) {
    xi   = x + sopflowipopt->xstarti[i];
    gi   = g + sopflowipopt->gstarti[i];

    opflow = sopflow->opflows[i];

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

    if(sopflow->nconineqcoup[i]) {
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

Bool eval_sopflow_jac_g(PetscInt n, PetscScalar *x, Bool new_x,
                PetscInt m, PetscInt nele_jac,
                PetscInt *iRow, PetscInt *jCol, PetscScalar *values,
                UserDataPtr user_data)
{

  PetscErrorCode ierr;
  SOPFLOW         sopflow=(SOPFLOW)user_data;
  SOPFLOWSolver_IPOPT sopflowipopt=(SOPFLOWSolver_IPOPT)sopflow->solver;
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
    ierr = VecGetArray(sopflow->X,&xarr);CHKERRQ(ierr);
    opflow0 = sopflow->opflows[0];
    for(i=0; i < sopflow->ns; i++) {
      opflow = sopflow->opflows[i];
      opflowipopt = (OPFLOWSolver_IPOPT)opflow->solver;

      roffset = sopflowipopt->gstarti[i];
      coffset = sopflowipopt->xstarti[i];

      xi = xarr + sopflowipopt->xstarti[i];
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

      if(sopflow->nconineqcoup[i]) {
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

	    x0loc = sopflowipopt->xstarti[0] + loc0;
	    xiloc = sopflowipopt->xstarti[i] + loc;
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
    ierr = VecRestoreArray(sopflow->X,&xarr);CHKERRQ(ierr);
  } else {

    opflow0 = sopflow->opflows[0];
    for(i=0; i < sopflow->ns; i++) {
      opflow = sopflow->opflows[i];
      opflowipopt = (OPFLOWSolver_IPOPT)opflow->solver;

      xi = x + sopflowipopt->xstarti[i];
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

      if(sopflow->nconineqcoup[i]) {
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

	    x0loc = sopflowipopt->xstarti[0] + loc0;
	    xiloc = sopflowipopt->xstarti[i] + loc;
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

Bool eval_sopflow_h(PetscInt n, PetscScalar *x, Bool new_x, PetscScalar obj_factor,
            PetscInt m, PetscScalar *lambda, Bool new_lambda,
            PetscInt nele_hess, PetscInt *iRow, PetscInt *jCol,
            PetscScalar *values, UserDataPtr user_data)
{
  PetscErrorCode       ierr;
  PetscInt             nrow;
  SOPFLOW             sopflow=(SOPFLOW)user_data;
  SOPFLOWSolver_IPOPT sopflowipopt=(SOPFLOWSolver_IPOPT)sopflow->solver;
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

  sopflow->obj_factor = obj_factor;

  if(values == NULL) {

    for(i=0; i < sopflow->ns; i++) {
      opflow = sopflow->opflows[i];
      opflowipopt = (OPFLOWSolver_IPOPT)opflow->solver;

      roffset = sopflowipopt->xstarti[i];

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

    for(i=0; i < sopflow->ns; i++) {
      opflow = sopflow->opflows[i];
      opflowipopt = (OPFLOWSolver_IPOPT)opflow->solver;
      opflow->obj_factor = obj_factor;

      roffset = sopflowipopt->xstarti[i];

      xi = x + roffset;
      lameqi = lambda + sopflowipopt->gstarti[i];

      ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
      ierr = VecPlaceArray(opflow->Lambdae,lameqi);CHKERRQ(ierr);
      if(opflow->nconineq) {
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

PetscErrorCode SOPFLOWSolverSolve_IPOPT(SOPFLOW sopflow)
{
  PetscErrorCode     ierr;
  SOPFLOWSolver_IPOPT sopflowipopt = (SOPFLOWSolver_IPOPT)sopflow->solver;
  OPFLOW             opflow;
  OPFLOWSolver_IPOPT opflowipopt;
  PetscScalar        *x,*xl,*xu,*g,*gl,*gu,*xi,*lameqi,*lamineqi,*lam;
  PetscInt           i;
  MatInfo            info_eq,info_ineq,info_hes;

  PetscFunctionBegin;

  sopflowipopt->nnz_jac_ge = sopflowipopt->nnz_jac_gi = 0;

  ierr = VecGetArray(sopflow->X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->Lambda,&lam);CHKERRQ(ierr);

  for(i=0; i < sopflow->ns; i++) {
    opflow = sopflow->opflows[i];
    opflowipopt = (OPFLOWSolver_IPOPT)opflow->solver;
    xi = x + sopflowipopt->xstarti[i];

    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);

    /* Compute nonzeros for the Jacobian */
    /* Equality constraint Jacobian */
    ierr = (*opflow->modelops.computeequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Ge);CHKERRQ(ierr);
    ierr = MatGetInfo(opflow->Jac_Ge,MAT_LOCAL,&info_eq);CHKERRQ(ierr);

    opflowipopt->nnz_jac_ge = info_eq.nz_used;
    sopflowipopt->nnz_jac_ge += opflowipopt->nnz_jac_ge;

    opflowipopt->nnz_jac_gi = 0;
    if(opflow->Nconineq) {
      ierr = (*opflow->modelops.computeinequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Gi);CHKERRQ(ierr);
      ierr = MatGetInfo(opflow->Jac_Gi,MAT_LOCAL,&info_ineq);CHKERRQ(ierr);

      opflowipopt->nnz_jac_gi = info_ineq.nz_used;
      
      sopflowipopt->nnz_jac_gi += opflowipopt->nnz_jac_gi;
    }   
    opflowipopt->nnz_jac_g = opflowipopt->nnz_jac_ge + opflowipopt->nnz_jac_gi;

    /* Add non-zeros for Jacobian of coupling constraints */
    if(sopflow->nconineqcoup[i]) sopflowipopt->nnz_jac_gi += 2*sopflow->nconineqcoup[i];

    /* Compute non-zeros for Hessian */

    lameqi = lam + sopflowipopt->gstarti[i];
    ierr = VecPlaceArray(opflow->Lambdae,lameqi);CHKERRQ(ierr);

    if(opflow->Nconineq) {
      lamineqi = lam + sopflowipopt->gstarti[i]+opflow->nconeq;
      ierr = VecPlaceArray(opflow->Lambdai,lamineqi);CHKERRQ(ierr);
    } else {
      opflow->Lambdai = NULL;
    }

    ierr = (*opflow->modelops.computehessian)(opflow,opflow->X,opflow->Lambdae,opflow->Lambdai,opflow->Hes);CHKERRQ(ierr);
    ierr = MatGetInfo(opflow->Hes,MAT_LOCAL,&info_hes);CHKERRQ(ierr);

    opflowipopt->nnz_hes = (info_hes.nz_used  -opflow->nx)/2 + opflow->nx;

    sopflowipopt->nnz_hes += opflowipopt->nnz_hes;

    ierr = VecResetArray(opflow->Lambdae);CHKERRQ(ierr);
    if(opflow->Nconineq) {
      ierr = VecResetArray(opflow->Lambdai);CHKERRQ(ierr);
    }
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(sopflow->X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->Lambda,&lam);CHKERRQ(ierr);

  sopflowipopt->nnz_jac_g = sopflowipopt->nnz_jac_ge + sopflowipopt->nnz_jac_gi;

  /* Create IPOPT solver instance */
  ierr = VecGetArray(sopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->Gu,&gu);CHKERRQ(ierr);

  /* Create IPOPT problem */
  sopflowipopt->nlp = CreateIpoptProblem(sopflow->Nx,xl,xu,sopflow->Ncon,gl,gu,sopflowipopt->nnz_jac_g,sopflowipopt->nnz_hes,0,&eval_sopflow_f,
				   &eval_sopflow_g,&eval_sopflow_grad_f,
				   &eval_sopflow_jac_g,&eval_sopflow_h);

  ierr = VecRestoreArray(sopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->Gu,&gu);CHKERRQ(ierr);

  ierr = VecGetArray(sopflow->X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->G,&g);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->Lambda,&lam);CHKERRQ(ierr);

  /* Solve */
  sopflowipopt->solve_status = IpoptSolve(sopflowipopt->nlp,x,g,&sopflow->obj,lam,NULL,NULL,sopflow);

  ierr = VecRestoreArray(sopflow->X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->G,&g);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->Lambda,&lam);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverDestroy_IPOPT(SOPFLOW sopflow)
{
  PetscErrorCode     ierr;
  SOPFLOWSolver_IPOPT ipopt = (SOPFLOWSolver_IPOPT)sopflow->solver;

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

PetscErrorCode SOPFLOWSolverGetObjective_IPOPT(SOPFLOW sopflow,PetscReal *obj)
{
  PetscFunctionBegin;
  *obj = sopflow->obj;
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetSolution_IPOPT(SOPFLOW sopflow,PetscInt scen_num,Vec *X)
{
  PetscErrorCode ierr;
  SOPFLOWSolver_IPOPT sopflowipopt = (SOPFLOWSolver_IPOPT)sopflow->solver;
  OPFLOW         opflow=sopflow->opflows[scen_num];
  Vec            Xi=opflow->X;
  PetscInt       nxi=opflow->nx;
  PetscScalar    *xi,*x;
  PetscInt       ix=sopflowipopt->xstarti[scen_num];

  PetscFunctionBegin;
  ierr = VecGetArray(Xi,&xi);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->X,&x);CHKERRQ(ierr);

  ierr = PetscArraycpy(xi,x+ix,nxi);CHKERRQ(ierr);

  ierr = VecRestoreArray(Xi,&xi);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->X,&x);CHKERRQ(ierr);

  *X = Xi;
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetConstraints_IPOPT(SOPFLOW sopflow,PetscInt scen_num,Vec *G)
{
  PetscErrorCode ierr;
  SOPFLOWSolver_IPOPT sopflowipopt = (SOPFLOWSolver_IPOPT)sopflow->solver;
  OPFLOW         opflow=sopflow->opflows[scen_num];
  Vec            Gi=opflow->G;
  PetscInt       ngi=opflow->ncon;
  PetscScalar    *gi,*g;
  PetscInt       ig=sopflowipopt->gstarti[scen_num];

  PetscFunctionBegin;
  ierr = VecGetArray(Gi,&gi);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->G,&g);CHKERRQ(ierr);

  ierr = PetscArraycpy(gi,g+ig,ngi);CHKERRQ(ierr);

  ierr = VecRestoreArray(Gi,&gi);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->G,&g);CHKERRQ(ierr);

  *G = Gi;

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetConstraintMultipliers_IPOPT(SOPFLOW sopflow,PetscInt scen_num,Vec *Lambda)
{
  PetscErrorCode ierr;
  SOPFLOWSolver_IPOPT sopflowipopt = (SOPFLOWSolver_IPOPT)sopflow->solver;
  OPFLOW         opflow=sopflow->opflows[scen_num];
  Vec            Lambdai=opflow->Lambda;
  PetscInt       ngi=opflow->ncon;
  PetscScalar    *lambdai,*lambda;
  PetscInt       ig=sopflowipopt->gstarti[scen_num];

  PetscFunctionBegin;
  ierr = VecGetArray(Lambdai,&lambdai);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->Lambda,&lambda);CHKERRQ(ierr);

  ierr = PetscArraycpy(lambdai,lambda+ig,ngi);CHKERRQ(ierr);

  ierr = VecRestoreArray(Lambdai,&lambdai);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->Lambda,&lambda);CHKERRQ(ierr);

  *Lambda = Lambdai;

  PetscFunctionReturn(0);
}


PetscErrorCode SOPFLOWSolverGetConvergenceStatus_IPOPT(SOPFLOW sopflow,PetscBool *status)
{
  SOPFLOWSolver_IPOPT ipopt = (SOPFLOWSolver_IPOPT)sopflow->solver;

  PetscFunctionBegin;
  if(ipopt->solve_status < 2) *status = PETSC_TRUE; /* See IpReturnCodes_inc.h in IPOPT. The first two denote convergence */
  else *status = PETSC_FALSE;

  PetscFunctionReturn(0);
}


PetscErrorCode SOPFLOWSolverSetUp_IPOPT(SOPFLOW sopflow)
{
  PetscErrorCode ierr;
  SOPFLOWSolver_IPOPT ipopt = (SOPFLOWSolver_IPOPT)sopflow->solver;
  PetscInt       i,j,k,ctr;
  OPFLOW         opflow,opflow0;
  PetscInt       ngenON;
  PetscScalar    *x,*xi,*xl,*xu,*xli,*xui,*gl,*gu,*gli,*gui;
  PS             ps,ps0;
  PSBUS          bus,bus0;
  PSGEN          gen,gen0;

  PetscFunctionBegin;

  ierr = PetscCalloc1(sopflow->ns,&ipopt->xstarti);CHKERRQ(ierr);
  ierr = PetscCalloc1(sopflow->ns,&ipopt->gstarti);CHKERRQ(ierr);
  ierr = PetscCalloc1(sopflow->ns,&ipopt->nxi);CHKERRQ(ierr);
  ierr = PetscCalloc1(sopflow->ns,&ipopt->ngi);CHKERRQ(ierr);

  sopflow->Nx = 0;
  sopflow->Ncon = 0;
  sopflow->Nconeq = 0;
  sopflow->Nconineq = 0;
  sopflow->Nconcoup = 0;
  ipopt->xstarti[0] = 0;
  ipopt->gstarti[0] = 0;

  for(i=0; i < sopflow->ns; i++) {
    opflow = sopflow->opflows[i];
    ierr = PSGetNumActiveGenerators(opflow->ps,&ngenON,NULL);CHKERRQ(ierr);
    ipopt->nxi[i] = opflow->nx;
    if(sopflow->iscoupling) sopflow->nconineqcoup[i] = (i == 0)?0:ngenON;
    else sopflow->nconineqcoup[i] = 0;

    ipopt->ngi[i] = opflow->ncon + sopflow->nconineqcoup[i];
    if(i < sopflow->ns - 1) {
      ipopt->xstarti[i+1] = ipopt->xstarti[i] + ipopt->nxi[i];
      ipopt->gstarti[i+1] = ipopt->gstarti[i] + ipopt->ngi[i];
    }
    sopflow->Nx +=        ipopt->nxi[i];
    sopflow->Ncon +=      ipopt->ngi[i];
    sopflow->Nconeq +=    opflow->nconeq;
    sopflow->Nconineq += opflow->nconineq;
    sopflow->Nconcoup +=  sopflow->nconineqcoup[i];
  }

  /* Create vector X */
  ierr = VecCreate(sopflow->comm->type,&sopflow->X);CHKERRQ(ierr);
  ierr = VecSetSizes(sopflow->X,sopflow->Nx,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(sopflow->X);CHKERRQ(ierr);

  ierr = VecDuplicate(sopflow->X,&sopflow->Xl);CHKERRQ(ierr);
  ierr = VecDuplicate(sopflow->X,&sopflow->Xu);CHKERRQ(ierr);

  /* vector for constraints */
  ierr = VecCreate(sopflow->comm->type,&sopflow->G);CHKERRQ(ierr);
  ierr = VecSetSizes(sopflow->G,sopflow->Ncon,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(sopflow->G);CHKERRQ(ierr);

  /* Constraint bounds vectors  */
  ierr = VecDuplicate(sopflow->G,&sopflow->Gl);CHKERRQ(ierr);
  ierr = VecDuplicate(sopflow->G,&sopflow->Gu);CHKERRQ(ierr);

  /* Lagrangian multipliers */
  ierr = VecDuplicate(sopflow->G,&sopflow->Lambda);CHKERRQ(ierr);

  /* Set Initial guess and Bounds on variables and constraints */
  ierr = VecGetArray(sopflow->X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->Gu,&gu);CHKERRQ(ierr);

  opflow0 = sopflow->opflows[0];
  for(i=0; i < sopflow->ns; i++) {
    opflow = sopflow->opflows[i];
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
    
    if(sopflow->nconineqcoup[i]) {
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

	  if(sopflow->makeup_power_source == 0) {
	    /* Only ref. bus responsible for make-up power for scenarios */
	    if(bus->ide == REF_BUS) {
	      gli[opflow->ncon + ctr] = -10000;
	      gui[opflow->ncon + ctr] =  10000;
	    } else {
	      gli[opflow->ncon + ctr] = 0.0;
	      gui[opflow->ncon + ctr] = 0.0;
	    }	    
	  } else {
	    gli[opflow->ncon + ctr] = -sopflow->mode*gen->ramp_rate_30min;
	    gui[opflow->ncon + ctr] =  sopflow->mode*gen->ramp_rate_30min;
	  }
	  ctr++;
	}
      }
    }
  }
  
  ierr = VecRestoreArray(sopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->Gu,&gu);CHKERRQ(ierr);

  /* Initialize Lagrange multiplier */
  ierr = VecSet(sopflow->Lambda,1.0);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverCreate_IPOPT(SOPFLOW sopflow)
{
  PetscErrorCode ierr;
  SOPFLOWSolver_IPOPT ipopt;
  
  PetscFunctionBegin;

  if(sopflow->comm->size > 1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"IPOPT solver does not support execution in parallel\n",sopflow->comm->size); 
  ierr = PetscCalloc1(1,&ipopt);CHKERRQ(ierr);

  ipopt->nlp = NULL;
  ipopt->nnz_jac_g = 0;
  ipopt->nnz_hes = 0;
  sopflow->solver = ipopt;

  sopflow->solverops.setup = SOPFLOWSolverSetUp_IPOPT;
  sopflow->solverops.solve = SOPFLOWSolverSolve_IPOPT;
  sopflow->solverops.destroy = SOPFLOWSolverDestroy_IPOPT;
  sopflow->solverops.getobjective = SOPFLOWSolverGetObjective_IPOPT;
  sopflow->solverops.getsolution  = SOPFLOWSolverGetSolution_IPOPT;
  sopflow->solverops.getconvergencestatus = SOPFLOWSolverGetConvergenceStatus_IPOPT;
  sopflow->solverops.getconstraints = SOPFLOWSolverGetConstraints_IPOPT;
  sopflow->solverops.getconstraintmultipliers = SOPFLOWSolverGetConstraintMultipliers_IPOPT;

  PetscFunctionReturn(0);
}

#endif
