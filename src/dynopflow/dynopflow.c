#include <private/psimpl.h>
#include <private/dynopflowimpl.h>
#include <private/opflowimpl.h>
#include <private/dynimpl.h>
#include "IpStdCInterface.h"

/* External Functions from OPFLOW */
extern Bool eval_opflow_f(PetscInt n, PetscScalar* x, Bool new_x,
            PetscScalar* obj_value, UserDataPtr user_data);

extern Bool eval_opflow_grad_f(PetscInt n, PetscScalar* x, Bool new_x,
                 PetscScalar* grad_f, UserDataPtr user_data);

extern Bool eval_opflow_g(PetscInt n, PetscScalar* x, Bool new_x,
            PetscInt m, PetscScalar* g, UserDataPtr user_data);

extern Bool eval_opflow_jac_g(PetscInt n, PetscScalar *x, Bool new_x,
                PetscInt m, PetscInt nele_jac,
                PetscInt *iRow, PetscInt *jCol, PetscScalar *values,
                UserDataPtr user_data);

extern Bool eval_opflow_h(PetscInt n, PetscScalar *x, Bool new_x, PetscScalar obj_factor,
            PetscInt m, PetscScalar *lambda, Bool new_lambda,
            PetscInt nele_hess, PetscInt *iRow, PetscInt *jCol,
            PetscScalar *values, UserDataPtr user_data);

Bool eval_dynopflow_f(PetscInt n, PetscScalar* x, Bool new_x,
            PetscScalar* obj_value, UserDataPtr user_data)
{
  Bool           status;
  DYNOPFLOW      dynopflow=(DYNOPFLOW)user_data;
  OPFLOW         opflow=dynopflow->opflow;

  status = eval_opflow_f(n,x,new_x,obj_value,(UserDataPtr)opflow);

  return TRUE;
}

Bool eval_dynopflow_grad_f(PetscInt n, PetscScalar* x, Bool new_x,
                 PetscScalar* grad_f, UserDataPtr user_data)
{
  Bool        status;
  DYNOPFLOW   dynopflow=(DYNOPFLOW)user_data;
  OPFLOW      opflow=dynopflow->opflow;

  status = eval_opflow_grad_f(n,x,new_x,grad_f,(UserDataPtr)opflow);

  return TRUE;
}

PetscErrorCode DYNOPFLOWUpdateDYNParameters(DYNOPFLOW,Vec);
PetscErrorCode DYNOPFLOWComputeSensiP(DYNOPFLOW);
Bool eval_dynopflow_g(PetscInt n, PetscScalar* x, Bool new_x,
             PetscInt m, PetscScalar* g, UserDataPtr user_data)
{
  Bool status;
  DYNOPFLOW   dynopflow=(DYNOPFLOW)user_data;
  OPFLOW      opflow=dynopflow->opflow;
  DYN         dyn=dynopflow->dyn;
  PetscInt    i;
  PetscErrorCode ierr;
  Vec    C;

  status = eval_opflow_g(n,x,new_x,m,g,(UserDataPtr)opflow);

  if(!dynopflow->initsolve) {
    if(dynopflow->solvesimultaneous) {
      ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
    
      /* Initialize cost integrand vector */

      ierr = TSGetCostIntegral(dyn->ts,&C);CHKERRQ(ierr);
      ierr = VecSet(C,0.0);CHKERRQ(ierr);

      /* Compute dynamic constraints */
      ierr = DYNOPFLOWUpdateDYNParameters(dynopflow,opflow->X);CHKERRQ(ierr);
      ierr = DYNSetDuration(dyn,dyn->maxsteps,dyn->tmax);CHKERRQ(ierr);
      ierr = DYNSetStartTimeAndStep(dyn,dyn->t0,dyn->dt0);CHKERRQ(ierr);
      ierr = DYNSolve(dynopflow->dyn);CHKERRQ(ierr);

      PetscScalar *c;
      ierr  = TSGetCostIntegral(dyn->ts,&C);CHKERRQ(ierr);
      ierr = VecView(C,0);CHKERRQ(ierr);
      
      ierr = VecGetArray(C,&c);CHKERRQ(ierr);
      for(i=0; i < dynopflow->Ncondyn; i++) g[dynopflow->Nconopflow+i] = c[i];
      ierr = VecRestoreArray(C,&c);CHKERRQ(ierr);

      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
    } else {
      PetscScalar *gpre;
      PetscReal    vdot;
      ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
      ierr = VecWAXPY(dynopflow->Xtemp,-1.0,dynopflow->Xpre,opflow->X);CHKERRQ(ierr);
      ierr = VecGetArray(dynopflow->Gdynpre,&gpre);CHKERRQ(ierr);

      for(i=0; i < dynopflow->Ncondyn; i++) {
	ierr = VecDot(dynopflow->mu[i],dynopflow->Xtemp,&vdot);CHKERRQ(ierr);
	g[dynopflow->Nconopflow+i] = gpre[i] - vdot;
	//	ierr = PetscPrintf(PETSC_COMM_SELF,"g[%d] = %2.5f\n",dynopflow->Nconopflow+i,g[dynopflow->Nconopflow+i]);CHKERRQ(ierr);
      }
      ierr = VecRestoreArray(dynopflow->Gdynpre,&gpre);CHKERRQ(ierr);
      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
    }
  }
  return TRUE;
}

Bool eval_dynopflow_jac_g(PetscInt n, PetscScalar *x, Bool new_x,
                PetscInt m, PetscInt nele_jac,
                PetscInt *iRow, PetscInt *jCol, PetscScalar *values,
                UserDataPtr user_data)
{
  Bool status;
  DYNOPFLOW   dynopflow=(DYNOPFLOW)user_data;
  OPFLOW      opflow=dynopflow->opflow;
  DYN         dyn=dynopflow->dyn;
  PetscInt    i,j,ctr=dynopflow->nnz_opflow_jac_g,rowstart=dynopflow->Nconopflow;
  Vec         *dcdp=dynopflow->mu;
  PetscScalar *dcp;
  PetscErrorCode ierr;

  status = eval_opflow_jac_g(n,x,new_x,m,nele_jac,iRow,jCol,values,(UserDataPtr)opflow);

  if(values == NULL) {
    for(i=0; i < dynopflow->Ncondyn;i++) {
      for(j=0; j < opflow->n; j++) {
	iRow[ctr] = rowstart+i;
	jCol[ctr] = j;
	ctr++;
      }
    }
  } else {
    if(!dynopflow->initsolve) {
      if(dynopflow->solvesimultaneous) {
	if(new_x) {
	  ierr = PetscPrintf(PETSC_COMM_SELF,"New x is true in jac_g");CHKERRQ(ierr);
	}
	ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
	ierr = DYNOPFLOWUpdateDYNParameters(dynopflow,opflow->X);CHKERRQ(ierr);
	ierr = VecResetArray(opflow->X);CHKERRQ(ierr);

	ierr = DYNAdjointSolve(dyn);CHKERRQ(ierr);

	for(i=0; i < dynopflow->Ncondyn;i++) {
	  ierr = VecGetArray(dcdp[i],&dcp);CHKERRQ(ierr);
	  for(j=0; j < opflow->n; j++) {
	    values[ctr] = dcp[j];
	    ctr++;
	  }
	  ierr = VecRestoreArray(dcdp[i],&dcp);CHKERRQ(ierr);
	}
      } else {
	for(i=0; i < dynopflow->Ncondyn;i++) {
	  ierr = VecGetArray(dcdp[i],&dcp);CHKERRQ(ierr);
	  for(j=0; j < opflow->n; j++) {
	    values[ctr] = -dcp[j];
	    ctr++;
	  }
	  ierr = VecRestoreArray(dcdp[i],&dcp);CHKERRQ(ierr);
	}
      }
    }
  }    
      
  return TRUE;
}

Bool eval_dynopflow_h(Index n, Number *x, Bool new_x, Number obj_factor,
            Index m, Number *lambda, Bool new_lambda,
            Index nele_hess, Index *iRow, Index *jCol,
            Number *values, UserDataPtr user_data)
{
  return FALSE;
}

/*
  DYNOPFLOWComputeDYNJacobianP - Computes the Jacobian of the DYN equations w.r.t. OPFLOW variables (parameters)

  Input Parameters:
+ ts - the TS solver
. t  - the current time
. x  - the solution vector
- ctx - application context (dynopflow)

  Output Parameters
. jacP - Jacobian of the DAE equations w.r.t. parameters (OPFLOW variables)
*/
PetscErrorCode DYNOPFLOWComputeDYNJacobianP(TS ts,PetscReal t,Vec X,Mat jacP,void *ctx)
{
  PetscErrorCode ierr;
  DYNOPFLOW      dynopflow=(DYNOPFLOW)ctx;
  DYN            dyn=dynopflow->dyn;
  OPFLOW         opflow=dynopflow->opflow;
  const PetscScalar *x,*xopflow;

  PetscFunctionBegin;
  ierr = MatZeroEntries(jacP);CHKERRQ(ierr);

  Vec localX,localXopflow;
  PetscInt i;
  PS       dynps=dyn->ps,opflowps=opflow->ps;
  PetscBool ghostbus;

  /* Get DYN local vector */
  ierr = DMGetLocalVector(dynps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dynps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dynps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  /* Get OPFLOW local vector */
  ierr = DMGetLocalVector(opflowps->networkdm,&localXopflow);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(opflowps->networkdm,opflow->X,INSERT_VALUES,localXopflow);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(opflowps->networkdm,opflow->X,INSERT_VALUES,localXopflow);CHKERRQ(ierr);


  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(localXopflow,&xopflow);CHKERRQ(ierr);

  for(i=0; i < dynps->nbus; i++) {
    PetscScalar VD, VQ,Vm,Vm2;
    PetscInt    dynloc,opflowloc;
    PSBUS       dynbus,opflowbus;
    PetscInt    k,Pgloc,Qgloc,Valoc,Vmloc,ctr=0;
    PetscScalar Vmopflow;
    PetscInt    row[2],col[1];
    PetscScalar val[2];
    
    dynbus = &dynps->bus[i];
    opflowbus = &opflowps->bus[i];

    ierr = PSBUSGetVariableGlobalLocation(dynbus,&dynloc);CHKERRQ(ierr);
    ierr = PSBUSGetVariableLocation(opflowbus,&opflowloc);CHKERRQ(ierr);
    Valoc = opflowloc; Vmloc = opflowloc+1;

    ierr = PSBUSIsGhosted(dynbus,&ghostbus);CHKERRQ(ierr);
    
    if(dynbus->ide == ISOLATED_BUS) continue;
    
    VD = x[dynloc]; VQ = x[dynloc+1];
    Vm = PetscSqrtScalar(VD*VD + VQ*VQ);
    Vm2 = Vm*Vm;

    for(k=0; k < dynbus->nload; k++) {
      Vmopflow = xopflow[opflowloc+1];
      PSLOAD load;
      ierr = PSBUSGetLoad(dynbus,k,&load);CHKERRQ(ierr);
      if(!load->status) continue;
      PetscScalar Pd,Qd,dyp_dVmopflow,dyq_dVmopflow;
      Pd = load->pl;
      Qd = load->ql;
	
      dyp_dVmopflow = 2*Pd/(dynbus->vm*dynbus->vm*dynbus->vm);
      dyq_dVmopflow = 2*Qd/(dynbus->vm*dynbus->vm*dynbus->vm);

      row[0] = dynloc; row[1] = dynloc+1;
      col[0] = Vmloc;
      val[0] = dyq_dVmopflow*VD - dyp_dVmopflow*VQ;
      val[1] = -dyp_dVmopflow*VD - dyq_dVmopflow*VQ;

      ierr = SetMatrixValues(jacP,2,row,1,col,val);CHKERRQ(ierr);
    }

    /* Generator frequency deviation */
    for(k=0; k < dynbus->ngen; k++) {
      PSGEN gen;
      DYNGenModel dyngen;
      const PetscScalar *xdyn;
      PetscInt    startloc;
      PetscScalar frequency;
      PetscScalar freq_min = dynopflow->freq_min;
      PetscScalar freq_max = dynopflow->freq_max;
      PetscInt    dynlocglob;

      ierr = PSBUSGetGen(dynbus,k,&gen);CHKERRQ(ierr);
      ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);

      startloc = dynloc;
      xdyn = x+startloc;
      
      dynlocglob = startloc + dyngen->startloc;
      if(gen->status) {
	Pgloc = opflowloc + 2 + ctr;
	Qgloc = opflowloc + 2 + ctr + 1;
	/* Set partial derivatives of machine DAE equations w.r.t parameters */
	ierr = DYNGenModelDAERHSJacobianP(dyngen,t,xdyn,jacP,dynlocglob,Valoc,Vmloc,Pgloc,Qgloc);CHKERRQ(ierr);

	if(dyngen->dynexc) {
	  DYNExcModel dynexc;

	  ierr = PSGENGetDYNExc(gen,&dynexc);CHKERRQ(ierr);

	  dynlocglob = startloc + dynexc->startloc;
	  ierr = DYNExcModelDAERHSJacobianP(dynexc,t,xdyn,jacP,dynlocglob,Valoc,Vmloc);CHKERRQ(ierr);
	}
	ctr += 2;
      }
    }
  }

  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(localXopflow,&xopflow);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dynps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(opflowps->networkdm,&localXopflow);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(jacP,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jacP,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  //  ierr = MatView(jacP,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DYNOPFLOWCostIntegrand(TS ts,PetscReal t,Vec X,Vec R,void *ctx)
{
  PetscErrorCode    ierr;
  PetscScalar       *r;
  const PetscScalar *x;
  DYNOPFLOW         dynopflow=(DYNOPFLOW)ctx;
  DYN               dyn=dynopflow->dyn;
  PetscInt          ctr=0;

  PetscFunctionBegin;
  ierr = VecSet(R,0.0);CHKERRQ(ierr);
  ierr = VecGetArray(R,&r);CHKERRQ(ierr);

  Vec localX;
  PetscInt i;
  PS       ps=dyn->ps;
  PetscBool ghostbus;

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    PetscScalar VD, VQ,Vm,Vm2;
    PetscInt    loc;
    PSBUS       bus;
    PetscInt    k;
    
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    
    if(bus->ide == ISOLATED_BUS) continue;
    
    VD = x[loc]; VQ = x[loc+1];
    Vm = PetscSqrtScalar(VD*VD + VQ*VQ);
    Vm2 = Vm*Vm;

    /* Generator frequency deviation */
    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;
      DYNGenModel dyngen;
      const PetscScalar *xdyn;
      PetscInt    startloc;
      PetscScalar frequency;
      PetscScalar freq_min = dynopflow->freq_min;
      PetscScalar freq_max = dynopflow->freq_max;
      PetscScalar freq_viol;

      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);

      startloc = loc;
      xdyn = x+startloc;
      
      if(!gen->status) {
	if(!dynopflow->sum_freq_cons) r[ctr++] = 0.0;
	continue;
      }
      /* Get Machine frequency */
      ierr = DYNGenModelGetFrequency(dyngen,t,xdyn,&frequency);
      //      ierr = PetscPrintf(PETSC_COMM_SELF,"Time %3.2f Gen[%d]: Freq = %2.5f\n",t,i,frequency);CHKERRQ(ierr);
      freq_viol = dynopflow->scal*PetscPowScalarInt(PetscMax(0.0,PetscMax(frequency-freq_max,freq_min-frequency)),dynopflow->exp);CHKERRQ(ierr);
      if(dynopflow->sum_freq_cons) r[0] += freq_viol;
      else r[ctr++] = freq_viol;
    }
  }

  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = VecRestoreArray(R,&r);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

static PetscErrorCode DYNOPFLOWDRDYFunction(TS ts,PetscReal t,Vec X,Vec *drdy,void *ctx)
{
  PetscErrorCode    ierr;
  PetscScalar       *ry;
  const PetscScalar *x;
  DYNOPFLOW         dynopflow=(DYNOPFLOW)ctx;
  DYN               dyn=dynopflow->dyn;
  PetscInt          veci=0;
  Vec localX;
  PetscInt i;
  PS       ps=dyn->ps;
  PetscBool ghostbus;

  PetscFunctionBegin;
  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);

  for(i=0; i < dynopflow->Ncondyn; i++) {
    ierr = VecSet(drdy[i],0.0);CHKERRQ(ierr);
  }
  
  if(dynopflow->sum_freq_cons) {
    ierr = VecGetArray(drdy[0],&ry);CHKERRQ(ierr);
  }

  for(i=0; i < ps->nbus; i++) {
    PetscScalar VD, VQ,Vm,Vm2;
    PetscInt    loc;
    PSBUS       bus;
    PetscInt    k;
    
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    
    if(bus->ide == ISOLATED_BUS) continue;
    
    VD = x[loc]; VQ = x[loc+1];
    Vm = PetscSqrtScalar(VD*VD + VQ*VQ);
    Vm2 = Vm*Vm;

    /* Generator frequency deviation */
    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;
      DYNGenModel dyngen;
      const PetscScalar *xdyn;
      PetscInt    startloc;
      PetscScalar frequency;     /* Machine frequency */
      PetscScalar dfreq_dstate; /* Partial derivative of machine frequency w.r.t. to its state (dw, dn,w, others,delta) that governs its frequency */
      PetscInt    stateloc;     /* Location of the start relative to its bus */
      PetscScalar freq_min = dynopflow->freq_min;
      PetscScalar freq_max = dynopflow->freq_max;

      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);

      startloc = loc;
      xdyn = x+startloc;

      if(gen->status) {
	/* Get Machine frequency and its partial derivative w.r.t machine state*/
	ierr = DYNGenModelGetFrequency(dyngen,t,xdyn,&frequency);CHKERRQ(ierr);
	ierr = DYNGenModelGetdFreqdState(dyngen,t,xdyn,&dfreq_dstate,&stateloc);

	if(!dynopflow->sum_freq_cons) {
	  ierr = VecGetArray(drdy[veci],&ry);CHKERRQ(ierr);
	}
	if(frequency > freq_max) ry[startloc+stateloc] = dynopflow->exp*dynopflow->scal*PetscPowScalarInt(dfreq_dstate,dynopflow->exp-1);
	else if(frequency < freq_min) ry[startloc+stateloc] = dynopflow->exp*dynopflow->scal*PetscPowScalarInt(-dfreq_dstate,dynopflow->exp-1);

	if(!dynopflow->sum_freq_cons) {
	  ierr = VecRestoreArray(drdy[veci++],&ry);CHKERRQ(ierr);
	}
      } else veci++;
    }
  }

  if(dynopflow->sum_freq_cons) {
    ierr = VecRestoreArray(drdy[0],&ry);CHKERRQ(ierr);
  }
  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);


  PetscFunctionReturn(0);
}

static PetscErrorCode DYNOPFLOWDRDPFunction(TS ts,PetscReal t,Vec U,Vec *drdp,void *ctx)
{
  PetscErrorCode    ierr;
  PetscScalar       *rp;
  const PetscScalar *u;
  DYNOPFLOW         dynopflow=(DYNOPFLOW)ctx;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(U,&u);CHKERRQ(ierr);
  ierr = VecGetArray(drdp[0],&rp);CHKERRQ(ierr);
  rp[0] = 0.;
  ierr = VecRestoreArray(drdp[0],&rp);CHKERRQ(ierr);
  ierr = VecGetArrayRead(U,&u);CHKERRQ(ierr);  
  PetscFunctionReturn(0);
}

/*
  DYNOPFLOWUpdateDYNParameters - Updates generator Pg,Qg and bus Vm and Va for DYN

  Input Parameters:
+ dynopflow - the DYNOPFLOW object
- X         - OPFLOW solution vector
*/
PetscErrorCode DYNOPFLOWUpdateDYNParameters(DYNOPFLOW dynopflow,Vec X)
{
  PetscErrorCode ierr;
  OPFLOW        opflow=dynopflow->opflow;
  DYN           dyn=dynopflow->dyn;
  PS            dynps=dyn->ps; /* PS object in DYN */
  PS            opflowps=opflow->ps; /* PS object in OPFLOW */
  PSBUS         opflowbus,dynbus;
  PetscScalar   Vm,Va;
  PetscInt      nbus=opflowps->nbus,i,loc,k;
  PetscScalar   *xopflow;

  PetscFunctionBegin;

  /* Get the array from the OPFLOW vector */
  ierr = VecGetArray(X,&xopflow);CHKERRQ(ierr);

  for(i=0; i < nbus; i++) {
    opflowbus = &opflowps->bus[i]; /* OPFLOW bus */
    dynbus    = &dynps->bus[i]; /* DYN bus */

    /* Get the location of variables in xopflow array for this bus */
    ierr = PSBUSGetVariableLocation(opflowbus,&loc);CHKERRQ(ierr);
    /* Copy Va and Vm from xopflow to dynbus */
    dynbus->va = xopflow[loc]*180.0/PETSC_PI;
    dynbus->vm = xopflow[loc+1];
    for(k=0; k < dynbus->ngen; k++) {
      PSGEN gen;
      ierr = PSBUSGetGen(dynbus,k,&gen);CHKERRQ(ierr);
      loc += 2;
      /* Copy generator set point voltage, real and reactive power */
      gen->vs = dynbus->vm;
      gen->pg = xopflow[loc];
      gen->qg = xopflow[loc+1];
    }
  }

  ierr = VecRestoreArray(X,&xopflow);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNOPFLOWCreate - Creates a DYNOPFLOW application context.

  Input Parameters:
. mpicomm - MPI communicator

  Output Parmater:
. dynopflowout - the new DYNOPFLOW application context

*/
PetscErrorCode DYNOPFLOWCreate(MPI_Comm mpicomm,DYNOPFLOW *dynopflowout)
{
  PetscErrorCode ierr;
  DYNOPFLOW          dynopflow;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&dynopflow);CHKERRQ(ierr);

  /* Create the communicator context */
  ierr = COMMCreate(mpicomm,&dynopflow->comm);CHKERRQ(ierr);

  /* Create the DYN context */
  ierr = DYNCreate(PETSC_COMM_SELF,&dynopflow->dyn);CHKERRQ(ierr);

  /* Create the OPFLOW context */
  ierr = OPFLOWCreate(PETSC_COMM_SELF,&dynopflow->opflow);CHKERRQ(ierr);

  dynopflow->scal = 1.0;
  dynopflow->exp = 2;
  dynopflow->eta = 1e-2;
  dynopflow->freq_max = freq + 0.25;
  dynopflow->freq_min = freq - 0.25;
  dynopflow->cpmaxit = 100;
  dynopflow->solvesimultaneous = PETSC_TRUE;
  dynopflow->sum_freq_cons = PETSC_TRUE;

  dynopflow->setupcalled = PETSC_FALSE;
  *dynopflowout = dynopflow;

  PetscFunctionReturn(0);
}

/*
  DYNOPFLOWDestroy - Destroys the dynopflow application object that was created with DYNOPFLOWCreate().

  Input parameter
. dynopflow - the DYNOPFLOW application context pointer
*/
PetscErrorCode DYNOPFLOWDestroy(DYNOPFLOW *dynopflow)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = COMMDestroy(&(*dynopflow)->comm);CHKERRQ(ierr);

  if(!(*dynopflow)->solvesimultaneous) {
    ierr = VecDestroy(&(*dynopflow)->Xpre);CHKERRQ(ierr);
    ierr = VecDestroy(&(*dynopflow)->Xtemp);CHKERRQ(ierr);
    ierr = VecDestroy(&(*dynopflow)->Gdynpre);CHKERRQ(ierr);
  }

  ierr = DYNDestroy(&(*dynopflow)->dyn);CHKERRQ(ierr);
  ierr = OPFLOWDestroy(&(*dynopflow)->opflow);CHKERRQ(ierr);

  ierr = PetscFree(*dynopflow);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNOPFLOWReadMatPowerData - Reads the network data given in MATPOWER data format 

  Input Parameter
+  DYNOPFLOW - The DYNOPFLOW object
-  netfile - The name of the network file

   Notes:
   With multiple processes, each process reads the data file.
*/
PetscErrorCode DYNOPFLOWReadMatPowerData(DYNOPFLOW dynopflow,const char netfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* DYN and OPFLOW read network data file */
  ierr = DYNReadMatPowerData(dynopflow->dyn,netfile);CHKERRQ(ierr);
  ierr = OPFLOWReadMatPowerData(dynopflow->opflow,netfile);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNOPFLOWCreateGlobalVector - Returns a global vector of the appropriate size
  and distribution conforming to the distribution of the PS object.

  Input Paramereters:
. dynopflow - the optimal power flow application object

  Output Parameters:
. vec - the global vector

  Notes:
  DYNOPFLOWSetUp() must be called before calling this routine.
*/
PetscErrorCode DYNOPFLOWCreateGlobalVector(DYNOPFLOW dynopflow,Vec *vec)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if(!dynopflow->setupcalled) SETERRQ(dynopflow->comm->type,0,"DYNOPFLOWSetUp() must be called before calling PFLOWCreateGlobalVector");
  ierr = OPFLOWCreateGlobalVector(dynopflow->opflow,vec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNOPFLOWReadDyrData - Reads the machine dynamic data file

  Input Parameter
+  DYNOPFLOW - The DYNOPFLOW object
-  dyrfile - The name of the .dyr file

   Notes:
   With multiple processes, each process reads the data file.
*/
PetscErrorCode DYNOPFLOWReadDyrData(DYNOPFLOW dynopflow,const char dyrfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* DYN reads dyr file */
  ierr = DYNReadDyrData(dynopflow->dyn,dyrfile);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNOPFLOWReadEventData - Reads the event data file and sets up the events in the DYNOPFLOW object

  Input Parameters:
+ DYNOPFLOW - the dynopf object
- eventfile - the name of the event file

  Notes:
    Each processor reads the data file.
*/
PetscErrorCode DYNOPFLOWReadEventData(DYNOPFLOW dynopflow,const char eventfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* DYN reads event data file */
  ierr = DYNReadEventData(dynopflow->dyn,eventfile);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNOPFLOWSetDuration - Sets the duration of the dynamic simulation

  Input Parameters
+  DYNOPFLOW - the DYNOPFLOW object
.  max_steps  - the maximum number of steps that the time-stepping solver takes
-  tend - the end time
*/
PetscErrorCode DYNOPFLOWSetDuration(DYNOPFLOW dynopflow,PetscInt max_steps,PetscReal tend)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DYNSetDuration(dynopflow->dyn,max_steps,tend);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNOPFLOWSetInitialTimeAndStep - Sets the start time and the time step of the dynamic simulation

  Input Parameters
+  DYNOPFLOW - the DYNOPFLOW object
.  start_time - start time
.  time_step  - the step size (in seconds)

   Notes:
   For variable time-stepping methods, this step is used as the initial time step.
*/
PetscErrorCode DYNOPFLOWSetStartTimeAndStep(DYNOPFLOW dynopflow,PetscReal start_time, PetscReal time_step)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DYNSetStartTimeAndStep(dynopflow->dyn,start_time,time_step);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

extern PetscErrorCode OPFLOWSetVariableandConstraintBounds(OPFLOW,Vec,Vec,Vec,Vec);

/*
  DYNOPFLOWSetVariableandConstraintBounds - Sets the bounds on variables and constraints

  Input Parameters:
+ dynopflow - the DYNOPFLOW object
. Xl        - lower bound on X
. Xu        - upper bound on X
. Gl        - lower bound on g
- Gu        - upper bound on g
*/
PetscErrorCode DYNOPFLOWSetVariableandConstraintBounds(DYNOPFLOW dynopflow,Vec Xl,Vec Xu,Vec Gl,Vec Gu)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=dynopflow->opflow;
  PetscScalar    *gl,*gu;
  PetscInt       i;

  PetscFunctionBegin;
  /* Set OPFLOW variable and constraint bounds */
  ierr = OPFLOWSetVariableandConstraintBounds(opflow,Xl,Xu,Gl,Gu);CHKERRQ(ierr);
  /* Set the bounds on the dynamic constraints */
  ierr = VecGetArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(Gu,&gu);CHKERRQ(ierr);

  /* The dynamic constraints are located after the opflow constraints */
  gl += opflow->Ncon; gu+= opflow->Ncon;
  
  for(i=0; i < dynopflow->Ncondyn; i++) {
    gl[i] = -1e-10;
    gu[i] = dynopflow->eta;
  }

  ierr = VecRestoreArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gu,&gu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

extern PetscErrorCode OPFLOWGetConstraintJacobianNonzeros(OPFLOW,PetscInt*);

/*
  DYNOPFLOWGetConstraintJacobianNonzeros - Gets the number of nonzeros in the constraint Jacobian

  Input Parameters:
. dynopflow - the DYNOPFLOW object

  Output Parameters:
. nnz - the number of nonzeros
*/
PetscErrorCode DYNOPFLOWGetConstraintJacobianNonzeros(DYNOPFLOW dynopflow,PetscInt *nnz)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=dynopflow->opflow;

  PetscFunctionBegin;
  /* Get the number of nonzeros due to the OPFLOW Jacobian */
  ierr = OPFLOWGetConstraintJacobianNonzeros(opflow,&dynopflow->nnz_opflow_jac_g);CHKERRQ(ierr);
  /* Nonzeros from the dynamic constraints */
  dynopflow->nnz_dyn_jac_g = dynopflow->Ncondyn*opflow->n;
  *nnz += dynopflow->nnz_opflow_jac_g + dynopflow->nnz_dyn_jac_g;
  PetscFunctionReturn(0);
}

/*
  DYNOPFLOWComputeSensiP - Computes the sensitivity of the dynamic constraints w.r.t. the parameters

  Input Parameters:
. dynopflow - the DYNOPFLOW object
*/
PetscErrorCode DYNOPFLOWComputeSensiP(DYNOPFLOW dynopflow)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=dynopflow->opflow;
  DYN            dyn=dynopflow->dyn;

  PetscFunctionBegin;
  /* Update DYN parameters */
  ierr = DYNOPFLOWUpdateDYNParameters(dynopflow,opflow->X);CHKERRQ(ierr);

  /* Compute sensitivity w.r.t. parameters */
  ierr = DYNComputeSensiP(dyn);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNOPFLOWSetUp - Sets up a dynamics constrained optimal power flow object

  Input Parameters:
. dynopflow - the DYNOPFLOW object

  Notes:
  This routine sets up the DYNOPFLOW object and the underlying PS object.

*/
PetscErrorCode DYNOPFLOWSetUp(DYNOPFLOW dynopflow)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=dynopflow->opflow;
  DYN            dyn=dynopflow->dyn;
  PetscInt       ngen=0,i,xdynsize,xopflowsize;

  PetscFunctionBegin;

  ierr = PetscOptionsGetReal(NULL, NULL,"-dynopflow_cost_scale",&dynopflow->scal,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL, NULL,"-dynopflow_cost_exp",&dynopflow->exp,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL, NULL,"-dynopflow_cost_eta",&dynopflow->eta,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL, NULL,"-dynopflow_freq_max",&dynopflow->freq_max,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL, NULL,"-dynopflow_freq_min",&dynopflow->freq_min,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL, NULL,"-dynopflow_solvesimultaneous",&dynopflow->solvesimultaneous,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL, NULL,"-dynopflow_sum_frequency_constraints",&dynopflow->sum_freq_cons,NULL);CHKERRQ(ierr);

  /* Set up the PS object in OPFLOW */
  ierr = PSSetUp(opflow->ps);CHKERRQ(ierr);
  /* Set up the PS object in DYN */
  ierr = PSSetUp(dyn->ps);CHKERRQ(ierr);

  ierr = PSGetNumGenerators(dyn->ps,&ngen,NULL);CHKERRQ(ierr);
  
  opflow->Nconeq = 2*opflow->ps->Nbus;
  opflow->Nconineq = 2*opflow->ps->Nbranch;
  opflow->Ncon = opflow->Nconeq + opflow->Nconineq;

  dynopflow->Nconopflow = 2*opflow->ps->Nbus + 2*opflow->ps->Nbranch;
  dynopflow->Ncondyn    = dyn->ncostfcns;
  dynopflow->Ncon = dynopflow->Nconopflow + dynopflow->Ncondyn;

  /* Set OPFLOW part */
  /* Set up the PS object in OPFLOW */
  ierr = PSSetUp(opflow->ps);CHKERRQ(ierr);

  /* Create the solution vector and upper and lower bounds */
  ierr = PSCreateGlobalVector(opflow->ps,&opflow->X);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->Xl);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->Xu);CHKERRQ(ierr);

  /* Create the constraint function vector and its bounds */
  ierr = VecCreate(opflow->ps->comm->type,&opflow->G);CHKERRQ(ierr);
  ierr = VecSetSizes(opflow->G,dynopflow->Ncon,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(opflow->G);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->G,&opflow->Gl);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->G,&opflow->Gu);CHKERRQ(ierr);

  /* Create the gradient vector for OPFLOW */
  ierr = VecDuplicate(opflow->X,&opflow->gradobj);CHKERRQ(ierr);

  ierr = VecGetSize(opflow->X,&opflow->n);CHKERRQ(ierr);
  ierr = VecGetSize(opflow->G,&opflow->m);CHKERRQ(ierr);

  /* Create vectors for multipliers */
  ierr = VecDuplicate(opflow->G,&opflow->lambda_g);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->lambda_xl);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->lambda_xu);CHKERRQ(ierr);

  /* Constraint jacobian */
  ierr = DYNOPFLOWGetConstraintJacobianNonzeros(dynopflow,&dynopflow->nnz_jac_g);CHKERRQ(ierr);

  /* Set up DYN object */
  ierr = DYNSetUp(dyn);CHKERRQ(ierr);

  /* Set up vectors for sensitivities */
  ierr = VecDuplicateVecs(dyn->X,dynopflow->Ncondyn,&dynopflow->lambda);CHKERRQ(ierr);
  ierr = VecDuplicateVecs(opflow->X,dynopflow->Ncondyn,&dynopflow->mu);CHKERRQ(ierr);
  ierr = VecDuplicateVecs(dyn->X,opflow->n,&dynopflow->dy0dp);CHKERRQ(ierr);
  
  /* Set the number of cost gradients for dyn->ts */
  ierr = TSSetCostGradients(dyn->ts,dynopflow->Ncondyn,dynopflow->lambda,dynopflow->mu);CHKERRQ(ierr);

  /* Create the Jacobian of the DAE equations w.r.t. parameters */
  ierr = VecGetSize(opflow->X,&xopflowsize);CHKERRQ(ierr);
  ierr = VecGetSize(dyn->X,&xdynsize);CHKERRQ(ierr);
  ierr = MatCreate(dyn->ps->comm->type,&dynopflow->Jacp);CHKERRQ(ierr);
  ierr = MatSetSizes(dynopflow->Jacp,xdynsize,xopflowsize,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = MatSetFromOptions(dynopflow->Jacp);CHKERRQ(ierr);
  ierr = MatSetUp(dynopflow->Jacp);CHKERRQ(ierr); /* Skipping allocation for now */

  /* Set functions needed by TS for adjoint calculation */
  ierr = TSSetSaveTrajectory(dyn->ts);CHKERRQ(ierr);
  ierr = TSAdjointSetRHSJacobian(dyn->ts,dynopflow->Jacp,DYNOPFLOWComputeDYNJacobianP,(void*)dynopflow);CHKERRQ(ierr);

  ierr = TSSetCostIntegrand(dyn->ts,dynopflow->Ncondyn,(PetscErrorCode (*)(TS,PetscReal,Vec,Vec,void*))DYNOPFLOWCostIntegrand,
				   (PetscErrorCode (*)(TS,PetscReal,Vec,Vec*,void*))DYNOPFLOWDRDYFunction,
				   (PetscErrorCode (*)(TS,PetscReal,Vec,Vec*,void*))DYNOPFLOWDRDPFunction,PETSC_TRUE,(void*)dynopflow);CHKERRQ(ierr);

  if(!dynopflow->solvesimultaneous) {
    ierr = VecDuplicate(opflow->X,&dynopflow->Xpre);CHKERRQ(ierr);
    ierr = VecDuplicate(opflow->X,&dynopflow->Xtemp);CHKERRQ(ierr);
    ierr = VecCreate(opflow->ps->comm->type,&dynopflow->Gdynpre);CHKERRQ(ierr);
    ierr = VecSetSizes(dynopflow->Gdynpre,dynopflow->Ncondyn,PETSC_DECIDE);CHKERRQ(ierr);
    ierr = VecSetFromOptions(dynopflow->Gdynpre);CHKERRQ(ierr);
  }
  opflow->setupcalled = PETSC_TRUE;
  dyn->setupcalled = PETSC_TRUE;
  dynopflow->setupcalled = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/*
  DYNOPFLOWSolveInitialize - Solves the OPF without dynamic constraints.

  Also, sets up the constraints and the options for Ipopt
*/
PetscErrorCode DYNOPFLOWSolveInitialize(DYNOPFLOW dynopflow)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=dynopflow->opflow;
  DYN            dyn=dynopflow->dyn;
  PetscScalar    *xl,*xu,*gl,*gu;
  PetscScalar    *x,*g,*lg,*lxl,*lxu;

  PetscFunctionBegin;
  if(!dynopflow->setupcalled) {
    ierr = DYNOPFLOWSetUp(dynopflow);CHKERRQ(ierr);
  }

  if(opflow->nlp) {
    FreeIpoptProblem(opflow->nlp);
  }

  /* Set bounds on variables and constraints */
  ierr = DYNOPFLOWSetVariableandConstraintBounds(dynopflow,opflow->Xl,opflow->Xu,opflow->Gl,opflow->Gu);CHKERRQ(ierr);

  ierr = VecGetArray(opflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(opflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecGetArray(opflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(opflow->Gu,&gu);CHKERRQ(ierr);

  /* Create IPOPT problem */
  opflow->nlp = CreateIpoptProblem(opflow->n,xl,xu,opflow->m,gl,gu,dynopflow->nnz_jac_g,0,0,&eval_dynopflow_f,
				   &eval_dynopflow_g,&eval_dynopflow_grad_f,
				   &eval_dynopflow_jac_g,&eval_dynopflow_h);

  ierr = VecRestoreArray(opflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->Gu,&gu);CHKERRQ(ierr);

  /* Options for IPOPT. This need to go through PetscOptionsBegin later */

  AddIpoptNumOption(opflow->nlp, "tol", 1e-4);
  AddIpoptNumOption(opflow->nlp, "acceptable_tol", 1e-4);
  AddIpoptNumOption(opflow->nlp, "mu_init", 1e-1);
  /*  AddIpoptNumOption(opflow->nlp,"bound_frac",1e-4);
  AddIpoptNumOption(opflow->nlp,"bound_push",1e-4);
  AddIpoptNumOption(opflow->nlp, "dual_inf_tol", 1e-2);
  AddIpoptNumOption(opflow->nlp, "compl_inf_tol", 1e-2);
  AddIpoptNumOption(opflow->nlp, "constr_viol_tol", 1e-2);
  */
  AddIpoptStrOption(opflow->nlp, "mu_strategy", "monotone");
  AddIpoptStrOption(opflow->nlp, "print_user_options", "yes");
  AddIpoptStrOption(opflow->nlp, "output_file", "ipopt.out");
  AddIpoptStrOption(opflow->nlp, "hessian_approximation", "limited-memory");

  /*  AddIpoptStrOption(opflow->nlp, "derivative_test", "first-order"); */

  /* Set Initial Guess */
  ierr = OPFLOWSetInitialGuess(opflow,opflow->X);CHKERRQ(ierr);

  /* Initial solve without dynamic constraints for initialization */
  ierr = VecGetArray(opflow->X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(opflow->G,&g);CHKERRQ(ierr);
  ierr = VecGetArray(opflow->lambda_g,&lg);CHKERRQ(ierr);
  ierr = VecGetArray(opflow->lambda_xl,&lxl);CHKERRQ(ierr);
  ierr = VecGetArray(opflow->lambda_xu,&lxu);CHKERRQ(ierr);

  dynopflow->initsolve = PETSC_TRUE;
  /* Initial Solve without dynamic constraints for initialization */
  opflow->solve_status = IpoptSolve(opflow->nlp,x,g,&opflow->obj,lg,lxl,lxu,dynopflow);

  ierr = VecRestoreArray(opflow->X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->G,&g);CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->lambda_g,&lg);CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->lambda_xl,&lxl);CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->lambda_xu,&lxu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNOPFLOWSolve - Solves the dynamics constrained optimal power flow

  Input Parameters:
. dynopflow - the dynopf application object
*/
PetscErrorCode DYNOPFLOWSolve(DYNOPFLOW dynopflow)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=dynopflow->opflow;
  DYN            dyn=dynopflow->dyn;
  PetscScalar    *xl,*xu,*gl,*gu;
  PetscScalar    *x,*g,*lg,*lxl,*lxu;
  PetscInt       i,ctr;

  PetscFunctionBegin;
  /* Run OPF without dynamic constraints for initialization */
  ierr = DYNOPFLOWSolveInitialize(dynopflow);CHKERRQ(ierr);

  /* Solve with dynamic constraints */
  if(dynopflow->solvesimultaneous) {
    ierr = VecGetArray(opflow->X,&x);CHKERRQ(ierr);
    ierr = VecGetArray(opflow->G,&g);CHKERRQ(ierr);
    ierr = VecGetArray(opflow->lambda_g,&lg);CHKERRQ(ierr);
    ierr = VecGetArray(opflow->lambda_xl,&lxl);CHKERRQ(ierr);
    ierr = VecGetArray(opflow->lambda_xu,&lxu);CHKERRQ(ierr);


    AddIpoptStrOption(opflow->nlp,"warm_start_init_point","yes");
    AddIpoptNumOption(opflow->nlp,"warm_start_bound_push",1e-16);
    AddIpoptNumOption(opflow->nlp,"warm_start_bound_frac",1e-16);
    AddIpoptNumOption(opflow->nlp,"warm_start_slack_bound_push",1e-16);
    AddIpoptNumOption(opflow->nlp,"warm_start_slack_bound_frac",1e-16);
    AddIpoptNumOption(opflow->nlp,"warm_start_mult_bound_push",1e-16);
    /*     AddIpoptStrOption(opflow->nlp, "nlp_scaling_method","none"); */

    dynopflow->initsolve = PETSC_FALSE;
    opflow->solve_status = IpoptSolve(opflow->nlp,x,g,&opflow->obj,lg,lxl,lxu,dynopflow);

    ierr = VecRestoreArray(opflow->X,&x);CHKERRQ(ierr);
    ierr = VecRestoreArray(opflow->G,&g);CHKERRQ(ierr);
    ierr = VecRestoreArray(opflow->lambda_g,&lg);CHKERRQ(ierr);
    ierr = VecRestoreArray(opflow->lambda_xl,&lxl);CHKERRQ(ierr);
    ierr = VecRestoreArray(opflow->lambda_xu,&lxu);CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_SELF,"Objective function value = %lf\n",opflow->obj);CHKERRQ(ierr);
    ierr = VecView(opflow->X,0);CHKERRQ(ierr);
  } else { /* Cutting plane method */
    /* Copy X to Xpre */
    ierr = VecCopy(opflow->X,dynopflow->Xpre);CHKERRQ(ierr);

    for(ctr=0; ctr < dynopflow->cpmaxit; ctr++) {
      /* Compute dynamic constraints and sensitivities */
    
      /* Initialize cost integrand */
      Vec    C;
      ierr = TSGetCostIntegral(dyn->ts,&C);CHKERRQ(ierr);
      ierr = VecSet(C,0.0);CHKERRQ(ierr);

      /* Compute dynamic constraints */
      ierr = DYNOPFLOWUpdateDYNParameters(dynopflow,opflow->X);CHKERRQ(ierr);
      ierr = DYNSetDuration(dyn,dyn->maxsteps,dyn->tmax);CHKERRQ(ierr);
      ierr = DYNSetStartTimeAndStep(dyn,dyn->t0,dyn->dt0);CHKERRQ(ierr);
      ierr = DYNSolve(dynopflow->dyn);CHKERRQ(ierr);
      /* - - - - - - - - - - - -  Adjoint solve - - - - - - - - - - - - - */
      /* Initialize lambda and mu */
      for(i=0; i < dynopflow->Ncondyn; i++) {
	ierr = VecSet(dynopflow->lambda[i],0.0);CHKERRQ(ierr);
	ierr = VecSet(dynopflow->mu[i],0.0);CHKERRQ(ierr);
      }
      
      
      ierr = TSAdjointSolve(dyn->ts);CHKERRQ(ierr);
      
      PetscScalar *c;
      PetscBool    converged=PETSC_TRUE;
      ierr  = TSGetCostIntegral(dyn->ts,&C);CHKERRQ(ierr);
      ierr = VecView(C,0);CHKERRQ(ierr);
      ierr = VecGetArray(C,&c);CHKERRQ(ierr);
      for(i=0; i < dynopflow->Ncondyn; i++) {
	if(c[i] > dynopflow->eta) {
	  converged = PETSC_FALSE;
	  break;
	}
      }
      ierr = VecRestoreArray(C,&c);CHKERRQ(ierr);
      if(converged) {
	ierr = PetscPrintf(PETSC_COMM_SELF,"Objective function value = %lf\n",opflow->obj);CHKERRQ(ierr);
	ierr = VecView(opflow->X,0);CHKERRQ(ierr);
	break;
      }

      ierr = VecCopy(C,dynopflow->Gdynpre);CHKERRQ(ierr);
      
      ierr = DYNOPFLOWComputeSensiP(dynopflow);CHKERRQ(ierr);

      /* Do optimization */
      ierr = VecGetArray(opflow->X,&x);CHKERRQ(ierr);
      ierr = VecGetArray(opflow->G,&g);CHKERRQ(ierr);
      ierr = VecGetArray(opflow->lambda_g,&lg);CHKERRQ(ierr);
      ierr = VecGetArray(opflow->lambda_xl,&lxl);CHKERRQ(ierr);
      ierr = VecGetArray(opflow->lambda_xu,&lxu);CHKERRQ(ierr);
      
      dynopflow->initsolve = PETSC_FALSE;
      opflow->solve_status = IpoptSolve(opflow->nlp,x,g,&opflow->obj,lg,lxl,lxu,dynopflow);
      
      ierr = VecRestoreArray(opflow->X,&x);CHKERRQ(ierr);
      ierr = VecRestoreArray(opflow->G,&g);CHKERRQ(ierr);
      ierr = VecRestoreArray(opflow->lambda_g,&lg);CHKERRQ(ierr);
      ierr = VecRestoreArray(opflow->lambda_xl,&lxl);CHKERRQ(ierr);
      ierr = VecRestoreArray(opflow->lambda_xu,&lxu);CHKERRQ(ierr);
      
      ierr = PetscPrintf(PETSC_COMM_SELF,"Objective function value = %lf\n",opflow->obj);CHKERRQ(ierr);
      //      ierr = VecView(opflow->X,0);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

