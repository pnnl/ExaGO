
#include <private/sopflowimpl.h>
#include <private/opflowimpl.h>
#include "sopflow-genramp.h"

PetscErrorCode SOPFLOWModelDestroy_GENRAMP(SOPFLOW sopflow)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscFree(sopflow->model);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSetVariableBounds_GENRAMP(SOPFLOW sopflow,Vec Xl,Vec Xu)
{
  PetscErrorCode ierr;
  OPFLOW         opflow;
  PetscScalar    *xl,*xu,*xli,*xui;
  PetscInt       i;
  PetscFunctionBegin;

  ierr = VecGetArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(Xu,&xu);CHKERRQ(ierr);

  for(i=0; i < sopflow->Ns; i++) {
    opflow = sopflow->opflows[i];

    /* Set bounds on variables */
    xli = xl + sopflow->xstarti[i];
    xui = xu + sopflow->xstarti[i];

    ierr = VecPlaceArray(opflow->Xl,xli);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Xu,xui);CHKERRQ(ierr);

    /* Set bounds */
    ierr = (*opflow->modelops.setvariablebounds)(opflow,opflow->Xl,opflow->Xu);CHKERRQ(ierr);

    ierr = VecResetArray(opflow->Xl);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Xu);CHKERRQ(ierr);

    if(i > 0 && opflow->has_gensetpoint) {
      /* Modify the bounds on ramping variables */
      PetscInt       j,k;
      PS             ps = opflow->ps;
      PSBUS          bus;
      PSGEN          gen;

      for(j = 0; j < ps->nbus; j++) {
	bus = &ps->bus[j];
	for(k=0; k < bus->ngen; k++) {
	  ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	  if(!gen->status) continue;
	  if(sopflow->mode == 0) continue;
	  else {
	    xli[gen->startxpdevloc] = -gen->ramp_rate_30min;
	    xui[gen->startxpdevloc] = gen->ramp_rate_30min;
	  }
	}
      }
    } 
  }

  ierr = VecRestoreArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Xu,&xu);CHKERRQ(ierr);

  PetscFunctionReturn(0);

}

PetscErrorCode SOPFLOWSetConstraintBounds_GENRAMP(SOPFLOW sopflow,Vec Gl,Vec Gu)
{
  PetscInt       i,j,k,ctr;
  PetscErrorCode ierr;
  PetscScalar    *gl,*gu,*gli,*gui;
  OPFLOW         opflow,opflow0;
  PS             ps,ps0;
  PSBUS          bus,bus0;
  PSGEN          gen,gen0;
  
  PetscFunctionBegin;

  ierr = VecGetArray(sopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->Gu,&gu);CHKERRQ(ierr);

  opflow0 = sopflow->opflows[0];
  for(i=0; i < sopflow->Ns; i++) {
    opflow = sopflow->opflows[i];

    /* Set bounds on constraints */
    gli = gl + sopflow->gstarti[i];
    gui = gu + sopflow->gstarti[i];
    
    ierr = VecPlaceArray(opflow->Gl,gli);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Gu,gui);CHKERRQ(ierr);

    ierr = (*opflow->modelops.setconstraintbounds)(opflow,opflow->Gl,opflow->Gu);CHKERRQ(ierr);

    ierr = VecResetArray(opflow->Gl);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Gu);CHKERRQ(ierr);
    
    if(sopflow->nconineqcoup[i]) {
      ctr = 0;
      ps    = opflow->ps;
      ps0 = opflow0->ps;
      /* Bounds on inequality coupling constraints */
      for(j=0; j < ps->nbus; j++) {
	bus    = &ps->bus[j];
	bus0 = &ps0->bus[j];

	for(k=0; k < bus->ngen; k++) {
	  ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	  ierr = PSBUSGetGen(bus0,k,&gen0);CHKERRQ(ierr);
	  if(!gen->status || !gen0->status) continue;

	  if(sopflow->mode == 0) {
	    /* Only ref. bus and renewable generation can deviate from base-case set points */
	    if(bus->ide == REF_BUS || gen->genfuel_type == GENFUEL_WIND) {
	      gli[opflow->ncon + ctr] = -10000;
	      gui[opflow->ncon + ctr] =  10000;
	    } else {
	      gli[opflow->ncon + ctr] = 0.0;
	      gui[opflow->ncon + ctr] = 0.0;
	    }	    
	  } else {
	    gli[opflow->ncon + ctr] = -gen->ramp_rate_30min;
	    gui[opflow->ncon + ctr] =  gen->ramp_rate_30min;
	  }

	  ctr++;
	}
      }
    }
  }
  
  ierr = VecRestoreArray(sopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->Gu,&gu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSetVariableandConstraintBounds_GENRAMP(SOPFLOW sopflow,Vec Xl,Vec Xu, Vec Gl, Vec Gu)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = SOPFLOWSetVariableBounds_GENRAMP(sopflow,Xl,Xu);CHKERRQ(ierr);
  ierr = SOPFLOWSetConstraintBounds_GENRAMP(sopflow,Gl,Gu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSetInitialGuess_GENRAMP(SOPFLOW sopflow,Vec X)
{
  PetscErrorCode ierr;
  OPFLOW         opflow;
  PetscScalar    *x,*xi,*xl,*xu,*xli,*xui;
  PetscInt       i;
  PetscFunctionBegin;

  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->Xu,&xu);CHKERRQ(ierr);

  for(i=0; i < sopflow->Ns; i++) {
    opflow = sopflow->opflows[i];
    /* Set initial guess and bounds on variables */
    xi  = x  + sopflow->xstarti[i];
    xli = xl + sopflow->xstarti[i];
    xui = xu + sopflow->xstarti[i];

    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Xl,xli);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Xu,xui);CHKERRQ(ierr);

    /* Set initial guess */
    ierr = OPFLOWSetInitialGuess(opflow,opflow->X);CHKERRQ(ierr);

    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Xl);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Xu);CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->Xu,&xu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWComputeJacobian_GENRAMP(SOPFLOW sopflow,Vec X,Mat J)
{
  PetscErrorCode ierr;
  OPFLOW         opflow,opflow0;
  PetscInt       roffset,coffset;
  PetscInt       nrow,ncol;
  PetscScalar    *xi,*x;
  PetscInt       i,j,k,loc,loc0,x0loc,xiloc;
  PS             ps,ps0;
  PSBUS          bus,bus0;
  PSGEN          gen,gen0;
  PetscInt       nvals;
  const PetscInt *cols;
  const PetscScalar *vals;
  PetscInt       row,col,gloc;
  PetscScalar    val;


  PetscFunctionBegin;

  ierr = VecGetArray(X,&x);CHKERRQ(ierr);

  opflow0 = sopflow->opflows[0];
  for(i=0; i < sopflow->Ns; i++) {
    opflow = sopflow->opflows[i];

    roffset = sopflow->gstarti[i];
    coffset = sopflow->xstarti[i];

    xi = x + sopflow->xstarti[i];
    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);

    ierr = (*opflow->modelops.computeequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Ge);CHKERRQ(ierr);

    ierr = MatGetSize(opflow->Jac_Ge,&nrow,&ncol);CHKERRQ(ierr);
    /* Copy over locations and values to Jac */
    for(j=0; j < nrow; j++) {
      ierr = MatGetRow(opflow->Jac_Ge,j,&nvals,&cols,&vals);CHKERRQ(ierr);
      for(k=0; k < nvals; k++) {
	row = roffset + j;
	col = coffset + cols[k];
	val = vals[k];
	ierr = MatSetValues(J,1,&row,1,&col,&val,INSERT_VALUES);
      }
      ierr = MatRestoreRow(opflow->Jac_Ge,j,&nvals,&cols,&vals);CHKERRQ(ierr);
    }

    if(i > 0 && opflow->has_gensetpoint) {
      ps  = opflow->ps;
      ps0 = opflow0->ps;
      for(j=0; j < ps->nbus; j++) {
	bus = &ps->bus[j];
	bus0 = &ps0->bus[j];

	for(k=0; k < bus->ngen; k++) {
	  ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	  ierr = PSBUSGetGen(bus0,k,&gen0);CHKERRQ(ierr);
	    
	  if(!gen->status) continue;

	  loc0 = gen0->startxpsetloc;
	  gloc = gen->starteqloc+1;
	    
	  x0loc = sopflow->xstarti[0] + loc0;
	  row = roffset + gloc;
	  col = x0loc;
	  val = -1.;
	  ierr = MatSetValues(J,1,&row,1,&col,&val,INSERT_VALUES);CHKERRQ(ierr);
	}
      }
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
	  row = roffset + j;
	  col = coffset + cols[k];
	  val = vals[k];
	  ierr = MatSetValues(J,1,&row,1,&col,&val,INSERT_VALUES);
	}
	ierr = MatRestoreRow(opflow->Jac_Gi,j,&nvals,&cols,&vals);CHKERRQ(ierr);
      }

      roffset += opflow->nconineq;
    }

    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);

    if(sopflow->nconineqcoup[i]) {
      ps  = opflow->ps;
      ps0 = opflow0->ps;
      for(j=0; j < ps->nbus; j++) {
	bus = &ps->bus[j];
	bus0 = &ps0->bus[j];
	for(k=0; k < bus->ngen; k++) {
	  ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	  ierr = PSBUSGetGen(bus0,k,&gen0);CHKERRQ(ierr);
	  
	  if(gen->status && gen0->status) {
	    loc0 = gen0->startxpowloc;
	    loc  = gen->startxpowloc;

	    x0loc = sopflow->xstarti[0] + loc0;
	    xiloc = sopflow->xstarti[i] + loc;
	    row = roffset;
	    col = x0loc;
	    val = -1.;
	    ierr = MatSetValues(J,1,&row,1,&col,&val,INSERT_VALUES);CHKERRQ(ierr);
	    col = xiloc;
	    val = 1.;
	    ierr = MatSetValues(J,1,&row,1,&col,&val,INSERT_VALUES);CHKERRQ(ierr);
	  
	    roffset += 1;
	  }
	}
      }
    }
  }

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWComputeConstraints_GENRAMP(SOPFLOW sopflow,Vec X,Vec G)
{
  PetscErrorCode ierr;
  OPFLOW    opflow0,opflow;
  PetscInt  i,j,k,loc,loc0,ctr;
  PetscScalar *x0,*x,*xi,*g,*gi;
  PS        ps,ps0;
  PSBUS     bus,bus0;
  PSGEN     gen,gen0;

  PetscFunctionBegin;

  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(G,&g);CHKERRQ(ierr);
  
  opflow0 = sopflow->opflows[0];
  ps0 = opflow0->ps;
  x0 = x + sopflow->xstarti[0];

  for(i=0; i < sopflow->Ns; i++) {
    xi   = x + sopflow->xstarti[i];
    gi   = g + sopflow->gstarti[i];

    opflow = sopflow->opflows[i];

    if(i > 0 && opflow->has_gensetpoint) {
      ps = opflow->ps;
      for(j=0; j < ps->nbus; j++) {
	bus  = &ps->bus[j];
	bus0 = &ps0->bus[j];

	for(k=0; k < bus->ngen; k++) {
	  ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	  ierr = PSBUSGetGen(bus0,k,&gen0);CHKERRQ(ierr);
	  if(!gen->status) continue;
	  /* Update the generator set-point */
	  gen->pgs = x0[gen0->startxpsetloc];
	}
      }
    }
    
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
	
	for(k=0; k < bus->ngen; k++) {
	  ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	  ierr = PSBUSGetGen(bus0,k,&gen0);CHKERRQ(ierr);

	  if(gen->status && gen0->status) {
	    loc0 = gen0->startxpowloc;
	    loc  = gen->startxpowloc;

	    gi[ctr] = xi[loc] - x0[loc0]; /* PG(i) - PG(0) */
	    ctr++;
	  }
	}
      }
    }
    
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(G,&g);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWComputeObjective_GENRAMP(SOPFLOW sopflow,Vec X,PetscScalar *obj)
{
  PetscErrorCode ierr;
  OPFLOW         opflow;
  PetscInt       i;
  PetscScalar    *xi;
  PetscScalar    *x;

  PetscFunctionBegin;
  *obj = 0.0;

  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  for(i=0; i < sopflow->Ns; i++) {
    xi = x + sopflow->xstarti[i];
    opflow = sopflow->opflows[i];
    opflow->obj = 0.0;
    
    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
    ierr = (*opflow->modelops.computeobjective)(opflow,opflow->X,&opflow->obj);CHKERRQ(ierr);
    *obj += opflow->obj;
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  }
  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
				
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWComputeGradient_GENRAMP(SOPFLOW sopflow,Vec X,Vec Grad)
{
  PetscErrorCode ierr;
  OPFLOW    opflow;
  PetscInt  i;
  PetscScalar *x,*xi,*grad,*gradi;

  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Grad,&grad);CHKERRQ(ierr);

  for(i=0; i < sopflow->Ns; i++) {
    opflow = sopflow->opflows[i];
    xi    = x + sopflow->xstarti[i];
    gradi = grad + sopflow->xstarti[i];


    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->gradobj,gradi);CHKERRQ(ierr);
    ierr = VecSet(opflow->gradobj,0.0);CHKERRQ(ierr);

    ierr = (*opflow->modelops.computegradient)(opflow,opflow->X,opflow->gradobj);CHKERRQ(ierr);
    
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->gradobj);CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Grad,&grad);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWComputeObjandGradient_GENRAMP(SOPFLOW sopflow,Vec X,PetscScalar *obj,Vec Grad)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = SOPFLOWComputeObjective_GENRAMP(sopflow,X,obj);CHKERRQ(ierr);
  ierr = SOPFLOWComputeGradient_GENRAMP(sopflow,X,Grad);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWModelSetNumVariablesandConstraints_GENRAMP(SOPFLOW sopflow,PetscInt *nxi,PetscInt *ngi,PetscInt *nconeqcoup,PetscInt *nconineqcoup)
{
  PetscInt i,ngenON;
  OPFLOW   opflow;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  for(i=0; i < sopflow->ns; i++) {
    opflow = sopflow->opflows[i];
    ierr = PSGetNumActiveGenerators(opflow->ps,&ngenON,NULL);CHKERRQ(ierr);
    nxi[i] = opflow->nx;
    if(sopflow->iscoupling) {
      if(opflow->has_gensetpoint) nconeqcoup[i] = 0;
      else nconineqcoup[i] = (i == 0)?0:ngenON;
    }
    ngi[i] = opflow->ncon + nconeqcoup[i] + nconineqcoup[i];
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWComputeHessian_GENRAMP(SOPFLOW sopflow,Vec X,Vec Lambda,Mat H)
{
  PetscErrorCode    ierr;
  PetscInt          nrow;
  OPFLOW            opflow;
  PetscScalar       *x,*xi,*lambda,*lameqi,*lamineqi;
  PetscInt          i;
  PetscInt          roffset;
  PetscInt          nvals;
  const PetscInt    *cols;
  const PetscScalar *vals;
  PetscInt          j,k;
  PetscInt          row,col;
  PetscScalar       val;

  PetscFunctionBegin;

  ierr = MatZeroEntries(H);CHKERRQ(ierr);
  
  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Lambda,&lambda);CHKERRQ(ierr);
  for(i=0; i < sopflow->Ns; i++) {
    opflow = sopflow->opflows[i];
    opflow->obj_factor = sopflow->obj_factor;

    roffset = sopflow->xstarti[i];

    xi = x + roffset;
    lameqi = lambda + sopflow->gstarti[i];

    ierr = VecPlaceArray(opflow->X,xi);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Lambdae,lameqi);CHKERRQ(ierr);
    if(opflow->Nconineq) {
      lamineqi = lameqi + opflow->nconeq;
      ierr = VecPlaceArray(opflow->Lambdai,lamineqi);CHKERRQ(ierr);
    }

    ierr = (*opflow->modelops.computehessian)(opflow,opflow->X,opflow->Lambdae,opflow->Lambdai,opflow->Hes);CHKERRQ(ierr);

    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Lambdae);CHKERRQ(ierr);
    if(opflow->Nconineq) {
      ierr = VecResetArray(opflow->Lambdai);CHKERRQ(ierr);
    }

    /* Copy over values */
    ierr = MatGetSize(opflow->Hes,&nrow,&nrow);CHKERRQ(ierr);
    for(j=0; j < nrow; j++) {
      ierr = MatGetRow(opflow->Hes,j,&nvals,&cols,&vals);CHKERRQ(ierr);
      row = roffset + j;
      for(k=0; k < nvals; k++) {
	col = roffset + cols[k];
	val = vals[k];
	ierr = MatSetValues(H,1,&row,1,&col,&val,INSERT_VALUES);CHKERRQ(ierr);
      }
      ierr = MatRestoreRow(opflow->Hes,j,&nvals,&cols,&vals);CHKERRQ(ierr);
    }
  }

  ierr = MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWModelCreate_GENRAMP(SOPFLOW sopflow)
{
  GENRAMP genramp;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  
  ierr = PetscCalloc1(1,&genramp);CHKERRQ(ierr);

  sopflow->model = genramp;

  /* Inherit Ops */
  sopflow->modelops.destroy = SOPFLOWModelDestroy_GENRAMP;
  sopflow->modelops.setnumvariablesandconstraints = SOPFLOWModelSetNumVariablesandConstraints_GENRAMP;
  sopflow->modelops.setvariablebounds = SOPFLOWSetVariableBounds_GENRAMP;
  sopflow->modelops.setconstraintbounds = SOPFLOWSetConstraintBounds_GENRAMP;
  sopflow->modelops.setvariableandconstraintbounds = SOPFLOWSetVariableandConstraintBounds_GENRAMP;
  sopflow->modelops.setinitialguess = SOPFLOWSetInitialGuess_GENRAMP;
  sopflow->modelops.computeconstraints = SOPFLOWComputeConstraints_GENRAMP;
  sopflow->modelops.computejacobian = SOPFLOWComputeJacobian_GENRAMP;
  sopflow->modelops.computehessian = SOPFLOWComputeHessian_GENRAMP;
  sopflow->modelops.computeobjandgradient = SOPFLOWComputeObjandGradient_GENRAMP;
  sopflow->modelops.computeobjective = SOPFLOWComputeObjective_GENRAMP;
  sopflow->modelops.computegradient  = SOPFLOWComputeGradient_GENRAMP;
  
  PetscFunctionReturn(0);
}
