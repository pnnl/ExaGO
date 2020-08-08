#include <scopflow_config.h>

#if defined(EXAGO_HAVE_HIOP)

#include <private/opflowimpl.h>
#include "pbpolhiop.h"

PetscErrorCode OPFLOWSetVariableandConstraintBounds_PBPOLHIOP(OPFLOW opflow,Vec Xl,Vec Xu, Vec Gl, Vec Gu)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = OPFLOWSetVariableBounds_PBPOLHIOP(opflow,Xl,Xu);CHKERRQ(ierr);
  ierr = OPFLOWSetConstraintBounds_PBPOLHIOP(opflow,Gl,Gu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetInitialGuess_PBPOLHIOP(OPFLOW opflow,Vec X)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  const PetscScalar    *xl,*xu;
  PetscScalar    *x;
  PetscInt       i;
  PSBUS          bus;
  PetscInt       loc;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(opflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArrayRead(opflow->Xu,&xu);CHKERRQ(ierr);
  
  for(i=0; i < ps->nbus; i++) {
    PetscInt k;

    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);

    if(bus->ide == ISOLATED_BUS) {
      x[loc] = bus->va*PETSC_PI/180.0;
      x[loc+1] = bus->vm;
    } else {
      if(opflow->initializationtype == OPFLOWINIT_MIDPOINT) {
	/* Initial guess for voltage angles and bounds on voltage magnitudes */
	x[loc]   = (xl[loc] + xu[loc])/2.0;
	x[loc+1] = (xl[loc+1] + xu[loc+1])/2.0;
      } else if(opflow->initializationtype == OPFLOWINIT_FROMFILE || opflow->initializationtype == OPFLOWINIT_ACPF) {
	x[loc] = bus->va*PETSC_PI/180.0;
	x[loc+1]   = PetscMax(bus->Vmin,PetscMin(bus->vm,bus->Vmax));
      } else if(opflow->initializationtype == OPFLOWINIT_FLATSTART) {
	x[loc] = 0.0;
	x[loc+1] = 1.0;
      }
    }

    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;

      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;

      loc = loc+2;

      if(opflow->initializationtype == OPFLOWINIT_MIDPOINT || opflow->initializationtype == OPFLOWINIT_FLATSTART) {
	x[loc]   = 0.5*(xl[loc] + xu[loc]);
	x[loc+1] = 0.5*(xl[loc+1] + xu[loc+1]);
      } else if(opflow->initializationtype == OPFLOWINIT_FROMFILE || opflow->initializationtype == OPFLOWINIT_ACPF) {
	x[loc] = PetscMax(gen->pb,PetscMin(gen->pg,gen->pt));
	x[loc+1] = PetscMax(gen->qb,PetscMin(gen->qg,gen->qt));
      }
    }

    if(opflow->include_loadloss_variables) {
      for(k=0; k < bus->nload; k++) {
	PSLOAD load;
	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
	loc += 2;
	/* Initial value for real and reactive power load loss */
	x[loc] = 0.0;
	x[loc+1] = 0.0;
      }
    } 

    if(opflow->include_powerimbalance_variables) {
      loc += 2;
      x[loc] = x[loc+1] = 0.0;
    }
  }

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(opflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(opflow->Xu,&xu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeConstraints_PBPOLHIOP(OPFLOW opflow,Vec X,Vec G)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeObjandGradient_PBPOLHIOP(OPFLOW opflow,Vec X,PetscScalar *obj,Vec Grad)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = OPFLOWComputeObjective_PBPOLHIOP(opflow,X,obj);CHKERRQ(ierr);
  ierr = OPFLOWComputeGradient_PBPOLHIOP(opflow,X,Grad);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelSetNumVariables_PBPOLHIOP(OPFLOW opflow,PetscInt *busnvar,PetscInt *branchnvar,PetscInt *nx)
{
  PetscInt i,ngen,nload,k;
  PS       ps=opflow->ps;
  PSBUS    bus;
  PSGEN    gen;
  PetscErrorCode ierr;
  PetscBool isghost;

  PetscFunctionBegin;
  
  *nx = 0;
  /* No variables for the branches */
  for(i=0; i < ps->nline; i++) {
    branchnvar[i] = 0;
    *nx += branchnvar[i];
  }

  /* Variables for the buses */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);
    //    if(isghost) continue;
    busnvar[i] = 2; /* 2 variables for the bus */
    ierr = PSBUSGetNGen(bus,&ngen);CHKERRQ(ierr);
    for(k=0; k < ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      busnvar[i] += 2; /* (2 variables for voltage + Pg, Qg for each gen) */
    }

    if(opflow->include_loadloss_variables) {
      ierr = PSBUSGetNLoad(bus,&nload);CHKERRQ(ierr);
      /* Load loss variables..Real and imaginary part of the load loss */
      busnvar[i] += 2*nload;
    }

    if(opflow->include_powerimbalance_variables) busnvar[i] += 2;
    if(!isghost) *nx += busnvar[i];
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelSetNumConstraints_PBPOLHIOP(OPFLOW opflow,PetscInt *branchnconeq,PetscInt *busnconeq,PetscInt *nconeq,PetscInt *nconineq)
{
  PetscInt i;
  PS       ps=opflow->ps;
  PSBUS    bus;
  PSLINE   line;
  PetscErrorCode ierr;
  PetscBool isghost;

  PetscFunctionBegin;
  *nconeq = *nconineq = 0;

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);
    if(isghost) continue;
    *nconeq += 2;
  }

  if(!opflow->ignore_lineflow_constraints) {
    for(i=0; i < ps->nline; i++) {
      line = &ps->line[i];
      if(line->status && line->rateA < 1e5) *nconineq += 2; /* Line flow constraints */
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeHessian_PBPOLHIOP(OPFLOW opflow,Vec X,Vec Lambdae,Vec Lambdai,Mat H)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = MatZeroEntries(H);CHKERRQ(ierr);

  /* Objective function Hessian */
  ierr = OPFLOWComputeObjectiveHessian_PBPOLHIOP(opflow,X,H);CHKERRQ(ierr);

  /* Equality constraints Hessian */
  ierr = OPFLOWComputeEqualityConstraintsHessian_PBPOLHIOP(opflow,X,Lambdae,H);CHKERRQ(ierr);
  
  /* Inequality constraints Hessian */
  if(opflow->nconineq) {
    ierr = OPFLOWComputeInequalityConstraintsHessian_PBPOLHIOP(opflow,X,Lambdai,H);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolutionToPS_PBPOLHIOP(OPFLOW opflow)
{
  PetscErrorCode ierr;
  PS             ps=(PS)opflow->ps;
  PetscInt       i,k;
  Vec            X,Lambda;
  PSBUS          bus;
  PSGEN          gen;
  PSLOAD         load;
  PSLINE         line;
  const PetscScalar *x,*lambda,*lambdae,*lambdai;
  PetscInt       loc,gloc=0;
  PetscScalar    Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  PetscScalar    Vmf,Vmt,thetaf,thetat,thetaft,thetatf;
  PetscScalar    Pf,Qf,Pt,Qt;
  PSBUS          busf,bust;
  const PSBUS    *connbuses;
  PetscInt       xlocf,xloct;

  PetscFunctionBegin;

  ierr = OPFLOWGetSolution(opflow,&X);CHKERRQ(ierr);
  ierr = OPFLOWGetConstraintMultipliers(opflow,&Lambda);CHKERRQ(ierr);

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(Lambda,&lambda);CHKERRQ(ierr);
  lambdae = lambda;
  if(opflow->Nconineq) {
    lambdai = lambdae + opflow->nconeq;
  }

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);

    bus->va = x[loc];
    bus->vm = x[loc+1];

    bus->mult_pmis = lambdae[gloc];
    bus->mult_qmis = lambdae[gloc+1];
    gloc += 2;

    for(k=0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;

      loc += 2;

      gen->pg = x[loc];
      gen->qg = x[loc+1];
    }

    if(opflow->include_loadloss_variables) {
      for(k=0; k < bus->nload; k++) {
	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
	loc += 2;
	load->pl = load->pl - x[loc];
	load->ql = load->ql - x[loc+1];
      }
    }

    if(opflow->include_powerimbalance_variables) {
      loc += 2;
      bus->pimb = x[loc];
      bus->qimb = x[loc+1];
    }
  }

  gloc = 0;

  if(!opflow->ignore_lineflow_constraints) {
    for(i=0; i<ps->nline; i++) {
      line = &ps->line[i];
      if(!line->status) {
	line->mult_sf = line->mult_st = 0.0;
	continue;
      }
      
      Gff = line->yff[0];
      Bff = line->yff[1];
      Gft = line->yft[0];
      Bft = line->yft[1];
      Gtf = line->ytf[0];
      Btf = line->ytf[1];
      Gtt = line->ytt[0];
      Btt = line->ytt[1];
      
      ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];
      
      ierr = PSBUSGetVariableLocation(busf,&xlocf);CHKERRQ(ierr);
      ierr = PSBUSGetVariableLocation(bust,&xloct);CHKERRQ(ierr);
      
      thetaf  = x[xlocf];
      Vmf     = x[xlocf+1];
      thetat  = x[xloct];
      Vmt     = x[xloct+1];
      thetaft = thetaf - thetat;
      thetatf = thetat - thetaf;
      
      Pf = Gff*Vmf*Vmf  + Vmf*Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
      Qf = -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
      
      Pt = Gtt*Vmt*Vmt  + Vmt*Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf));
      Qt = -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));

      line->pf = Pf;
      line->qf = Qf;
      line->pt = Pt;
      line->qt = Qt;
      line->sf = PetscSqrtScalar(Pf*Pf + Qf*Qf);
      line->st = PetscSqrtScalar(Pt*Pt + Qt*Qt);

      if(line->rateA > 1e5) {
	line->mult_sf = line->mult_st = 0.0;
      } else {
	line->mult_sf = lambdai[gloc];
	line->mult_st = lambdai[gloc+1];
	gloc += 2;
      }
    }
  }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Lambda,&lambda);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelSetUp_PBPOLHIOP(OPFLOW opflow)
{
  PetscErrorCode ierr;
  PBPOLHIOP      pbpolhiop=(PBPOLHIOP)opflow->model;
  PetscInt       *idxn2sd_map;
  PS             ps;
  PSBUS          bus;
  int            ngen,nxsparse,nxdense;
  PSGEN          gen;
  
  PetscFunctionBegin;

  /* Create natural to sparse dense variable mapping */
  
  idxn2sd_map = opflow->idxn2sd_map;
  
  ps = opflow->ps;
  nxsparse = 2*ps->ngenON;
  nxdense  = 2*ps->nbus;
  
  int i,k;
  int spct=0,dnct=0;
  int loc;
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    PSBUSGetVariableLocation(bus,&loc);

    idxn2sd_map[loc] = nxsparse + dnct;
    idxn2sd_map[loc+1] = nxsparse + dnct+1;

    dnct += 2;
    loc += 2;
    PSBUSGetNGen(bus,&ngen);
    for(k=0; k < ngen; k++) {
      PSBUSGetGen(bus,k,&gen);
      if(!gen->status) continue;

      idxn2sd_map[loc] = spct;
      idxn2sd_map[loc+1] = spct + 1;

      spct += 2;
      loc += 2;
    }
  }

  ierr = CreateBusParams(opflow,&pbpolhiop->busparams);
  ierr = CreateGenParams(opflow,&pbpolhiop->genparams);
  ierr = CreateLineParams(opflow,&pbpolhiop->lineparams);
  ierr = CreateLoadParams(opflow,&pbpolhiop->loadparams);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelDestroy_PBPOLHIOP(OPFLOW opflow)
{
  PetscErrorCode ierr;
  PBPOLHIOP pbpolhiop=(PBPOLHIOP)opflow->model;

  PetscFunctionBegin;

  ierr = DestroyBusParams(&pbpolhiop->busparams);
  ierr = DestroyGenParams(&pbpolhiop->genparams);
  ierr = DestroyLineParams(&pbpolhiop->lineparams);
  ierr = DestroyLoadParams(&pbpolhiop->loadparams);

  ierr = PetscFree(opflow->model);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelCreate_PBPOLHIOP(OPFLOW opflow)
{
  PBPOLHIOP pbpol;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  
  ierr = PetscCalloc1(1,&pbpol);CHKERRQ(ierr);

  opflow->model = pbpol;

  opflow->spdnordering = PETSC_TRUE;

  /* Inherit Ops */
  opflow->modelops.destroy                              = OPFLOWModelDestroy_PBPOLHIOP;
  opflow->modelops.setnumvariables                      = OPFLOWModelSetNumVariables_PBPOLHIOP;
  opflow->modelops.setnumconstraints                    = OPFLOWModelSetNumConstraints_PBPOLHIOP;
  opflow->modelops.setvariablebounds                    = OPFLOWSetVariableBounds_PBPOLHIOP;
  opflow->modelops.setvariableboundsarray               = OPFLOWSetVariableBoundsArray_PBPOLHIOP;
  opflow->modelops.setconstraintbounds                  = OPFLOWSetConstraintBounds_PBPOLHIOP;
  opflow->modelops.setconstraintboundsarray             = OPFLOWSetConstraintBoundsArray_PBPOLHIOP;
  opflow->modelops.setvariableandconstraintbounds       = OPFLOWSetVariableandConstraintBounds_PBPOLHIOP;
  opflow->modelops.setinitialguess                      = OPFLOWSetInitialGuess_PBPOLHIOP;
  opflow->modelops.computeequalityconstraints           = OPFLOWComputeEqualityConstraints_PBPOLHIOP;
  opflow->modelops.computeinequalityconstraints         = OPFLOWComputeInequalityConstraints_PBPOLHIOP;
  opflow->modelops.computeequalityconstraintsarray      = OPFLOWComputeEqualityConstraintsArray_PBPOLHIOP;
  opflow->modelops.computeinequalityconstraintsarray    = OPFLOWComputeInequalityConstraintsArray_PBPOLHIOP;
  opflow->modelops.computeequalityconstraintjacobian    = OPFLOWComputeEqualityConstraintJacobian_PBPOLHIOP;
  opflow->modelops.computeinequalityconstraintjacobian  = OPFLOWComputeInequalityConstraintJacobian_PBPOLHIOP;
  opflow->modelops.computehessian                       = OPFLOWComputeHessian_PBPOLHIOP;
  opflow->modelops.computeobjandgradient                = OPFLOWComputeObjandGradient_PBPOLHIOP;
  opflow->modelops.computeobjective                     = OPFLOWComputeObjective_PBPOLHIOP;
  opflow->modelops.computegradient                      = OPFLOWComputeGradient_PBPOLHIOP;
  opflow->modelops.computeobjectivearray                = OPFLOWComputeObjectiveArray_PBPOLHIOP;
  opflow->modelops.computegradientarray                 = OPFLOWComputeGradientArray_PBPOLHIOP;
  opflow->modelops.solutiontops                         = OPFLOWSolutionToPS_PBPOLHIOP;
  opflow->modelops.setup                                = OPFLOWModelSetUp_PBPOLHIOP;
  opflow->modelops.setsparsejacobianlocationshiop       = OPFLOWSetSparseJacobianLocations_PBPOLHIOP;
  
  PetscFunctionReturn(0);
}

#endif
