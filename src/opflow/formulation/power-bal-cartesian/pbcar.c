#include <private/opflowimpl.h>
#include "pbcar.h"

PetscErrorCode OPFLOWFormulationDestroy_PBCAR(OPFLOW opflow)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscFree(opflow->formulation);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetVariableBounds_PBCAR(OPFLOW opflow,Vec Xl,Vec Xu)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscScalar    *xl,*xu;
  PetscInt       i;
  PSBUS          bus;
  PetscInt       loc;
  PetscBool      isghost;
  Vec            localXl,localXu;

  PetscFunctionBegin;

  ierr = DMGetLocalVector(ps->networkdm,&localXl);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ps->networkdm,&localXu);CHKERRQ(ierr);
  ierr = VecSet(Xl,0.0);CHKERRQ(ierr);
  ierr = VecSet(Xu,0.0);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,Xl,INSERT_VALUES,localXl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,Xl,INSERT_VALUES,localXl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,Xu,INSERT_VALUES,localXu);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,Xu,INSERT_VALUES,localXu);CHKERRQ(ierr);

  /* Get array pointers */
  ierr = VecGetArray(localXl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(localXu,&xu);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    PetscInt k;

    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);
    if(isghost) continue;

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);

    /* Bounds on real and imaginary part of voltages */
    xl[loc] = -bus->Vmax; xu[loc] = bus->Vmax;
    xl[loc+1] = -bus->Vmax; xu[loc+1] = bus->Vmax;

    if(bus->ide == ISOLATED_BUS) {
      xl[loc] = xu[loc] = bus->vm*PetscCosScalar(bus->va*PETSC_PI/180.0);
      xl[loc+1] = xu[loc+1] = bus->vm*PetscSinScalar(bus->va*PETSC_PI/180.0);
    }

    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      loc = loc+2;
      xl[loc] = gen->pb; /* PGmin */
      xu[loc] = gen->pt; /* PGmax */
      xl[loc+1] = gen->qb; /* QGmin */
      xu[loc+1] = gen->qt; /* QGmax */
      /* pb, pt, qb, qt are converted in p.u. in ps.c */
    }

    if(opflow->include_loadloss_variables) {
      for(k=0; k < bus->nload; k++) {
	PSLOAD load;
	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
	loc += 2;
	if(!load->status) xl[loc] = xu[loc] = xl[loc+1] = xu[loc+1] = 0.0;
	else {
	  xl[loc] = PetscMin(0.0,load->pl);
	  xu[loc] = PetscMax(0.0,load->pl);
	  xl[loc+1] = PetscMin(0.0,load->ql);
	  xu[loc+1] = PetscMax(0.0,load->ql);
	}
      }
    }

    if(opflow->include_powerimbalance_variables) {
      loc += 2;
      xl[loc] = xl[loc+1] = PETSC_NINFINITY;
      xu[loc] = xu[loc+1] = PETSC_INFINITY;
    }
  }
  
  ierr = VecRestoreArray(localXl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(localXu,&xu);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(ps->networkdm,localXl,ADD_VALUES,Xl);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(ps->networkdm,localXl,ADD_VALUES,Xl);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(ps->networkdm,localXu,ADD_VALUES,Xu);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(ps->networkdm,localXu,ADD_VALUES,Xu);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(ps->networkdm,&localXl);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localXu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetConstraintBounds_PBCAR(OPFLOW opflow,Vec Gl,Vec Gu)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscScalar    *gl,*gu;
  PetscInt       i;
  PSLINE         line;
  PSBUS          bus;
  PetscInt       gloc=0;
  PetscBool      isghost;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(Gu,&gu);CHKERRQ(ierr);

  /* Equality constraint bounds for balance equations and ref. bus angle */
  for(i=0; i < ps->nbus; i++) {

    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);
    if(isghost) continue;

    gl[gloc]   = 0.0;   gu[gloc]   = 0.0;
    gl[gloc+1] = 0.0;   gu[gloc+1] = 0.0;

    gloc += 2;
    if(bus->ide == REF_BUS) {
      /* Equality constraint on reference bus angle */
      /* V_I -V_R*tan(Va) */
      gl[gloc] = gu[gloc] = 0;
      gloc++;
    }
  }

  if(opflow->nconineq) {
    /* Inequality constraints on voltage magnitude */
    for(i=0; i < ps->nbus; i++) {
      bus = &ps->bus[i]; 
      ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);
      if(isghost) continue;
     
      gl[gloc] = bus->Vmin*bus->Vmin;
      gu[gloc] = bus->Vmax*bus->Vmax;
      if(bus->ide == ISOLATED_BUS) {
	gl[gloc] = PETSC_NINFINITY;
	gu[gloc] = PETSC_INFINITY;
      }
      gloc++;
    }
    
    
    for(i=0; i < ps->nbranch; i++) {
      line = &ps->line[i];
      if(!line->status || line->rateA > 1e5) continue;
      
      /* Line flow inequality constraints */
      gl[gloc] = gl[gloc+1] = 0.0; 
      gu[gloc] = gu[gloc+1] = (line->rateA/ps->MVAbase)*(line->rateA/ps->MVAbase);
      gloc += 2;
    }
  }

  ierr = VecRestoreArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gu,&gu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetVariableandConstraintBounds_PBCAR(OPFLOW opflow,Vec Xl,Vec Xu, Vec Gl, Vec Gu)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = OPFLOWSetVariableBounds_PBCAR(opflow,Xl,Xu);CHKERRQ(ierr);
  ierr = OPFLOWSetConstraintBounds_PBCAR(opflow,Gl,Gu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetInitialGuess_PBCAR(OPFLOW opflow,Vec X)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  const PetscScalar *xl,*xu;
  PetscScalar    *x;
  PetscInt       i;
  PSBUS          bus;
  PetscInt       loc;
  PetscBool      isghost;
  Vec            localX,localXl,localXu;

  PetscFunctionBegin;

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ps->networkdm,&localXl);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ps->networkdm,&localXu);CHKERRQ(ierr);

  ierr = VecSet(X,0.0);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,opflow->Xl,INSERT_VALUES,localXl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,opflow->Xl,INSERT_VALUES,localXl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,opflow->Xu,INSERT_VALUES,localXu);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,opflow->Xu,INSERT_VALUES,localXu);CHKERRQ(ierr);


  /* Get array pointers */
  ierr = VecGetArray(localX,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(localXl,&xl);CHKERRQ(ierr);
  ierr = VecGetArrayRead(localXu,&xu);CHKERRQ(ierr);
  
  for(i=0; i < ps->nbus; i++) {
    PetscInt k;

    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);
    if(isghost) continue;

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);

    if(bus->ide == ISOLATED_BUS) {
      x[loc]   = bus->vm*PetscCosScalar(bus->va*PETSC_PI/180.0);
      x[loc+1] = bus->vm*PetscSinScalar(bus->va*PETSC_PI/180.0);
    } else {
      if(opflow->initializationtype == OPFLOWINIT_MIDPOINT) {
	x[loc]   = 0.5*(bus->Vmin + bus->Vmax)*PetscCosScalar(bus->va*PETSC_PI/180.0);
	x[loc+1] = 0; //0.5*(bus->Vmin + bus->Vmax)*PetscSinScalar(bus->va*PETSC_PI/180.0);
      } else if(opflow->initializationtype == OPFLOWINIT_FROMFILE || opflow->initializationtype == OPFLOWINIT_ACPF) {
	x[loc]   = PetscMax(bus->Vmin,PetscMin(bus->vm,bus->Vmax))*PetscCosScalar(bus->va*PETSC_PI/180.0);
	x[loc+1] = PetscMax(bus->Vmin,PetscMin(bus->vm,bus->Vmax))*PetscSinScalar(bus->va*PETSC_PI/180.0);
      } else if(opflow->initializationtype == OPFLOWINIT_FLATSTART) {
	x[loc] = 1.0;
	x[loc+1] = 0.0;
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

  ierr = VecRestoreArray(localX,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(localXl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(localXu,&xu);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(ps->networkdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(ps->networkdm,localX,ADD_VALUES,X);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localXl);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localXu);CHKERRQ(ierr);

  //  ierr = VecView(X,0);CHKERRQ(ierr);


  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeEqualityConstraints_PBCAR(OPFLOW opflow,Vec X,Vec Ge)
{
  PetscErrorCode ierr;
  PetscInt       i,k,nconnlines;
  PetscInt       xloc,xlocf,xloct;
  PetscScalar    Pg,Qg,Pd,Qd;
  PetscScalar    Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  PetscScalar    Vrf,Vif,Vrt,Vit;
  PetscScalar    Pf,Qf,Pt,Qt;
  PetscScalar    Vr,Vi,Vm;
  PS             ps=opflow->ps;
  PSLOAD         load;
  PSLINE         line;
  PSBUS          bus,busf,bust;
  PSGEN          gen;
  const PSBUS    *connbuses;
  const PSLINE   *connlines;
  PetscScalar    *ge;
  const PetscScalar *x;
  Vec            localX;
  PetscBool      isghost;
  PetscInt       gloc=0;

  PetscFunctionBegin;
  ierr = VecSet(Ge,0.0);CHKERRQ(ierr);
  ierr = VecSet(opflow->Gelocal,0.0);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = VecGetArray(opflow->Gelocal,&ge);CHKERRQ(ierr);
  
  for (i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);

    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);

    /* Real and imaginary voltages for the bus */
    Vr = x[xloc];
    Vi = x[xloc+1];
    Vm = PetscSqrtScalar(Vr*Vr + Vi*Vi);
    
    if(!isghost) {
      if (bus->ide == ISOLATED_BUS) {
	ge[gloc]   = Vr - bus->vm*PetscCosScalar(bus->va*PETSC_PI/180.0);
	ge[gloc+1] = Vi - bus->vm*PetscSinScalar(bus->va*PETSC_PI/180.0);
  
	gloc += 2;
	continue;
      }

      /* Shunt injections */
      ge[gloc]   += Vm*Vm*bus->gl;
      ge[gloc+1] += -Vm*Vm*bus->bl;

      /* Generation injection */
      for (k=0; k < bus->ngen; k++) {
	ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	if(!gen->status) continue;
	xloc += 2;
	Pg = x[xloc];
	Qg = x[xloc+1];
      
	ge[gloc]   += -Pg;
	ge[gloc+1] += -Qg;
      }

      /* Load injection */
      for (k=0; k < bus->nload; k++) {
	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
	if(opflow->include_loadloss_variables) {
	  xloc += 2;
	  Pd = load->pl - x[xloc];
	  Qd = load->ql - x[xloc+1];
	} else {
	  Pd = load->pl;
	  Qd = load->ql;
	}
	ge[gloc]   += Pd;
	ge[gloc+1] += Qd;
      }

      /* Power imbalance addition */
      if(opflow->include_powerimbalance_variables) {
	PetscScalar Pimb,Qimb;
	xloc += 2;
	Pimb = x[xloc];
	Qimb = x[xloc+1];
	ge[gloc]   += Pimb;
	ge[gloc+1] += Qimb;
      }
    }

    /* Branch flow injections */
    ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);
    for (k=0; k < nconnlines; k++) {
      line = connlines[k];
      if (!line->status) continue;

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
      
      Vrf  = x[xlocf];
      Vif  = x[xlocf+1];
      Vrt  = x[xloct];
      Vit  = x[xloct+1];

      if (bus == busf) {
	Pf =  Gff*(Vrf*Vrf + Vif*Vif) + Vrf*(Gft*Vrt - Bft*Vit) + Vif*(Bft*Vrt + Gft*Vit);
	Qf = -Bff*(Vrf*Vrf + Vif*Vif) + Vif*(Gft*Vrt - Bft*Vit) - Vrf*(Bft*Vrt + Gft*Vit);

        ge[gloc]   += Pf;
        ge[gloc+1] += Qf;
      } else {
	Pt =  Gtt*(Vrt*Vrt + Vit*Vit) + Vrt*(Gtf*Vrf - Btf*Vif) + Vit*(Btf*Vrf + Gtf*Vif);
	Qt = -Btt*(Vrt*Vrt + Vit*Vit) + Vit*(Gtf*Vrf - Btf*Vif) - Vrt*(Btf*Vrf + Gtf*Vif);

        ge[gloc]   += Pt;
        ge[gloc+1] += Qt;
      }
    }
    gloc += 2;
    if(bus->ide == REF_BUS) {
      if(!isghost) {
	/* Equality constraint for angle */
	ge[gloc] = Vi - Vr*PetscTanScalar(bus->va*PETSC_PI/180.0);
      }
      gloc++;
    }
  }

  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

  ierr = VecRestoreArray(opflow->Gelocal,&ge);CHKERRQ(ierr);

  ierr = VecScatterBegin(opflow->scattereqcon,opflow->Gelocal,Ge,ADD_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
  ierr = VecScatterEnd(opflow->scattereqcon,opflow->Gelocal,Ge,ADD_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeEqualityConstraintJacobian_PBCAR(OPFLOW opflow,Vec X,Mat Je)
{
  PetscErrorCode ierr;
  PetscInt       i,k,row[3],col[4],genctr,loadctr;
  PetscInt       nconnlines,locglob,loc,locglobf,locglobt,locf,loct;
  PetscScalar    Vr,Vi,val[8],Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  PetscScalar    Vrf,Vrt,Vif,Vit;
  PS             ps=opflow->ps;
  PSBUS          bus;
  PSLINE         line;
  PSBUS          busf,bust;
  PSGEN          gen;
  const PSLINE   *connlines;
  const PSBUS    *connbuses;
  const PetscScalar *xarr;
  PetscBool      isghost;
  Vec            localX;

  PetscFunctionBegin;
  ierr = MatZeroEntries(Je);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);
  for (i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);

    row[0] = opflow->eqconglobloc[i]; row[1] = row[0] + 1;

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bus,&locglob);CHKERRQ(ierr);

    Vr = xarr[loc];
    Vi = xarr[loc+1];

    col[0] = locglob; col[1] = locglob+1;

    if(!isghost) {
      /* Isolated and reference bus */
      if(bus->ide == ISOLATED_BUS) {
	val[0] = val[3] = 1.0;
	val[1] = val[2] = 0.0;
	ierr = MatSetValues(Je,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
	continue;
      }

      /* Partial derivative for shunt contribution */
      val[0] = 2*Vr*bus->gl;  val[1] = 2*Vi*bus->gl;
      val[2] = -2*Vr*bus->bl; val[3]= -2*Vi*bus->bl; 
      ierr = MatSetValues(Je,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
    
      genctr = 0;
      for (k=0; k < bus->ngen; k++) {
	ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	if(!gen->status) continue;
	val[0] = -1;
	col[0] = locglob + 2 + genctr;
	ierr = MatSetValues(Je,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
      
	col[0] = locglob + 2 + genctr+1;
	ierr = MatSetValues(Je,1,row+1,1,col,val,ADD_VALUES);CHKERRQ(ierr);
	genctr += 2;
      }

      if(opflow->include_loadloss_variables) {
	loadctr = 0;
	for(k=0; k < bus->nload; k++) {
	  val[0] = -1;
	  col[0] = locglob + 2 + 2*bus->ngen + loadctr;
	  ierr = MatSetValues(Je,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	  col[0] = locglob + 2 + 2*bus->ngen + loadctr + 1;
	  ierr = MatSetValues(Je,1,row+1,1,col,val,ADD_VALUES);CHKERRQ(ierr);
	  loadctr += 2;
	}
      }

      /* Power imbalance Jacobian terms */
      if(opflow->include_powerimbalance_variables) {
	val[0] = 1;
	col[0] = locglob + 2 + 2*bus->ngen + opflow->include_loadloss_variables*2*bus->nload;
	ierr = MatSetValues(Je,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	col[0] = locglob + 2 + 2*bus->ngen + opflow->include_loadloss_variables*2*bus->nload + 1;
	ierr = MatSetValues(Je,1,row+1,1,col,val,ADD_VALUES);CHKERRQ(ierr);
      }
    }

    /* Partial derivatives of network equations */
    /* Get the lines supporting the bus */
    ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);

    for (k=0; k < nconnlines; k++) {
      line = connlines[k];
      if (!line->status) continue;
      Gff = line->yff[0];
      Bff = line->yff[1];
      Gft = line->yft[0];
      Bft = line->yft[1];
      Gtf = line->ytf[0];
      Btf = line->ytf[1];
      Gtt = line->ytt[0];
      Btt = line->ytt[1];

      /* Get the connected buses to this line */
      ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];

      ierr = PSBUSGetVariableGlobalLocation(busf,&locglobf);CHKERRQ(ierr);
      ierr = PSBUSGetVariableGlobalLocation(bust,&locglobt);CHKERRQ(ierr);

      ierr = PSBUSGetVariableLocation(busf,&locf);CHKERRQ(ierr);
      ierr = PSBUSGetVariableLocation(bust,&loct);CHKERRQ(ierr);

      Vrf = xarr[locf];  Vif = xarr[locf+1];
      Vrt = xarr[loct];  Vit = xarr[loct+1];

      if (bus == busf) {
      	col[0] = locglobf; col[1] = locglobf+1; col[2] = locglobt; col[3] = locglobt+1;
	/* dPf_dVrf */
	val[0] = 2*Gff*Vrf + (Gft*Vrt - Bft*Vit);
	/* dPf_dVif */
	val[1] = 2*Gff*Vif + (Bft*Vrt + Gft*Vit);
	/* dPf_dVrt */
	val[2] = Vrf*Gft + Vif*Bft;
	/* dPf_dVit */
	val[3] = -Vrf*Bft + Vif*Gft;
	
	/* dQf_dVrf */
	val[4] = -2*Bff*Vrf - (Bft*Vrt + Gft*Vit);
	/* dQf_dVif */
	val[5] = -2*Bff*Vif + (Gft*Vrt - Bft*Vit);
	/* dQf_dVrt */
	val[6] =  Vif*Gft - Vrf*Bft;
	/* dQf_dVit */
	val[7] = -Vif*Bft - Vrf*Gft;
        ierr = MatSetValues(Je,2,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
      } else {
      	col[0] = locglobt; col[1] = locglobt+1; col[2] = locglobf; col[3] = locglobf+1;
	val[0] = 2*Gtt*Vrt + (Gtf*Vrf - Btf*Vif);
	val[1] = 2*Gtt*Vit + (Btf*Vrf + Gtf*Vif);
	val[2] = Vrt*Gtf  + Vit*Btf;
	val[3] = -Vrt*Btf + Vit*Gtf;

	val[4] = -2*Btt*Vrt - (Btf*Vrf + Gtf*Vif);
	val[5] = -2*Btt*Vit + (Gtf*Vrf - Btf*Vif);
	val[6] =  Vit*Gtf - Vrt*Btf;
	val[7] = -Vit*Btf - Vrt*Gtf;

        ierr = MatSetValues(Je,2,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
      }
    }

    if(!isghost) {
      if(bus->ide == REF_BUS) {
	/* Partial of equality constraint for ref. bus angle */
	row[2] = row[1] + 1;
	col[0] = locglob; col[1] = locglob + 1;
	val[0] = -PetscTanScalar(bus->va*PETSC_PI/180.0);
	val[1] = 1.0;
	ierr = MatSetValues(Je,1,row+2,2,col,val,ADD_VALUES);CHKERRQ(ierr);
      }
    }
  }

  ierr = VecRestoreArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(Je,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Je,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeInequalityConstraints_PBCAR(OPFLOW opflow,Vec X,Vec Gi)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscInt       gloc=0;
  PetscInt       xloc,xlocf,xloct;
  PetscScalar    *g;
  PetscScalar    Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  PetscScalar    Vr,Vi,Vrf,Vrt,Vif,Vit;
  PetscScalar    Pf,Qf,Pt,Qt,Sf2,St2;
  PS             ps=opflow->ps;
  PSBUS          bus;
  PSLINE         line;
  PSBUS          busf,bust;
  const PSBUS    *connbuses;
  const PetscScalar *x;
  PetscBool      isghost;
  Vec            localX;

  PetscFunctionBegin;
  ierr = VecSet(Gi,0.0);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Gi,&g);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);

    if(isghost) continue;

    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);

    Vr = x[xloc];
    Vi = x[xloc+1];
    
    g[gloc] = Vr*Vr + Vi*Vi;
    gloc++;
  }

  for(i=0; i<ps->nbranch; i++) {
    line = &ps->line[i];
    if(!line->status || line->rateA > 1e5) continue;

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

    Vrf  = x[xlocf];
    Vif  = x[xlocf+1];
    Vrt  = x[xloct];
    Vit  = x[xloct+1];

    Pf =  Gff*(Vrf*Vrf + Vif*Vif) + Vrf*(Gft*Vrt - Bft*Vit) + Vif*(Bft*Vrt + Gft*Vit);
    Qf = -Bff*(Vrf*Vrf + Vif*Vif) + Vif*(Gft*Vrt - Bft*Vit) - Vrf*(Bft*Vrt + Gft*Vit);

    Pt =  Gtt*(Vrt*Vrt + Vit*Vit) + Vrt*(Gtf*Vrf - Btf*Vif) + Vit*(Btf*Vrf + Gtf*Vif);
    Qt = -Btt*(Vrt*Vrt + Vit*Vit) + Vit*(Gtf*Vrf - Btf*Vif) - Vrt*(Btf*Vrf + Gtf*Vif);

    Sf2 = Pf*Pf + Qf*Qf;
    St2 = Pt*Pt + Qt*Qt;

    g[gloc]   = Sf2;
    g[gloc+1] = St2;

    gloc = gloc + 2;
  }
  
  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gi,&g);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeInequalityConstraintJacobian_PBCAR(OPFLOW opflow,Vec X,Mat Ji)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscInt       row[2],col[4];
  PetscInt       rstart,rend;
  PetscInt       gloc=0,xloc,xlocf,xloct,xlocglob,xlocfglob,xloctglob;
  PetscScalar    val[4];
  PetscScalar    Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  PetscScalar    Vr,Vi,Vrf,Vrt,Vif,Vit;
  PetscScalar    Pf,Qf,Pt,Qt;
  PetscScalar    dSf2_dPf, dSf2_dQf, dSt2_dPt, dSt2_dQt;
  PetscScalar    dPf_dVrf,dPf_dVif,dPf_dVrt,dPf_dVit;
  PetscScalar    dQf_dVrf,dQf_dVif,dQf_dVrt,dQf_dVit;
  PetscScalar    dPt_dVrf,dPt_dVif,dPt_dVrt,dPt_dVit;
  PetscScalar    dQt_dVrf,dQt_dVif,dQt_dVrt,dQt_dVit;
  PetscScalar    dSf2_dVrf,dSf2_dVif,dSf2_dVrt,dSf2_dVit;
  PetscScalar    dSt2_dVrf,dSt2_dVif,dSt2_dVrt,dSt2_dVit;
  PS             ps=opflow->ps;
  MPI_Comm       comm=opflow->comm->type;
  PSLINE         line;
  PSBUS          bus,busf,bust;
  const PSBUS    *connbuses;
  const PetscScalar *x;
  Vec            localX;
  PetscBool      isghost;

  PetscFunctionBegin;
  ierr = MatZeroEntries(Ji);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(Ji,&rstart,&rend);CHKERRQ(ierr);
  ierr = MPI_Barrier(comm);CHKERRQ(ierr);

  gloc = rstart;
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);
    if(isghost) continue;

    row[0] = gloc;

    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bus,&xlocglob);CHKERRQ(ierr);

    Vr = x[xloc]; Vi = x[xloc+1];

    col[0] = xlocglob; col[1] = xlocglob+1;
    val[0] = 2*Vr; val[1] = 2*Vi;

    ierr = MatSetValues(Ji,1,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);

    gloc += 1;
  }
    
  for (i=0; i < ps->nbranch; i++) {
    line = &ps->line[i];
    if(!line->status || line->rateA > 1e5) continue;

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
    ierr = PSBUSGetVariableGlobalLocation(busf,&xlocfglob);CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bust,&xloctglob);CHKERRQ(ierr);


    Vrf  = x[xlocf];
    Vif  = x[xlocf+1];
    Vrt  = x[xloct];
    Vit  = x[xloct+1];

    Pf =  Gff*(Vrf*Vrf + Vif*Vif) + Vrf*(Gft*Vrt - Bft*Vit) + Vif*(Bft*Vrt + Gft*Vit);
    Qf = -Bff*(Vrf*Vrf + Vif*Vif) + Vif*(Gft*Vrt - Bft*Vit) - Vrf*(Bft*Vrt + Gft*Vit);

    Pt =  Gtt*(Vrt*Vrt + Vit*Vit) + Vrt*(Gtf*Vrf - Btf*Vif) + Vit*(Btf*Vrf + Gtf*Vif);
    Qt = -Btt*(Vrt*Vrt + Vit*Vit) + Vit*(Gtf*Vrf - Btf*Vif) - Vrt*(Btf*Vrf + Gtf*Vif);

    dSf2_dPf = 2*Pf;
    dSf2_dQf = 2*Qf;
    dSt2_dPt = 2*Pt;
    dSt2_dQt = 2*Qt;

    dPf_dVrf = 2*Gff*Vrf + (Gft*Vrt - Bft*Vit);
    dPf_dVif = 2*Gff*Vif + (Bft*Vrt + Gft*Vit);
    dPf_dVrt = Vrf*Gft + Vif*Bft;
    dPf_dVit = -Vrf*Bft + Vif*Gft;

    dQf_dVrf = -2*Bff*Vrf - (Bft*Vrt + Gft*Vit);
    dQf_dVif = -2*Bff*Vif + (Gft*Vrt - Bft*Vit);
    dQf_dVrt =  Vif*Gft - Vrf*Bft;
    dQf_dVit = -Vif*Bft - Vrf*Gft;

    dPt_dVrt = 2*Gtt*Vrt + (Gtf*Vrf - Btf*Vif);
    dPt_dVit = 2*Gtt*Vit + (Btf*Vrf + Gtf*Vif);
    dPt_dVrf = Vrt*Gtf  + Vit*Btf;
    dPt_dVif = -Vrt*Btf + Vit*Gtf;

    dQt_dVrt = -2*Btt*Vrt - (Btf*Vrf + Gtf*Vif);
    dQt_dVit = -2*Btt*Vit + (Gtf*Vrf - Btf*Vif);
    dQt_dVrf =  Vit*Gtf - Vrt*Btf;
    dQt_dVif  = -Vit*Btf - Vrt*Gtf;

    dSf2_dVrf = dSf2_dPf*dPf_dVrf + dSf2_dQf*dQf_dVrf;
    dSf2_dVrt = dSf2_dPf*dPf_dVrt + dSf2_dQf*dQf_dVrt;
    dSf2_dVif = dSf2_dPf*dPf_dVif + dSf2_dQf*dQf_dVif;
    dSf2_dVit = dSf2_dPf*dPf_dVit + dSf2_dQf*dQf_dVit;

    row[0] = gloc;
    col[0] = xlocfglob; col[1] = xlocfglob+1; col[2] = xloctglob; col[3] = xloctglob+1;
    val[0] = dSf2_dVrf;
    val[1] = dSf2_dVif;
    val[2] = dSf2_dVrt;
    val[3] = dSf2_dVit;
    ierr = MatSetValues(Ji,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

    dSt2_dVrf = dSt2_dPt*dPt_dVrf + dSt2_dQt*dQt_dVrf;
    dSt2_dVrt = dSt2_dPt*dPt_dVrt + dSt2_dQt*dQt_dVrt;
    dSt2_dVif = dSt2_dPt*dPt_dVif + dSt2_dQt*dQt_dVif;
    dSt2_dVit = dSt2_dPt*dPt_dVit + dSt2_dQt*dQt_dVit;

    row[0] = gloc+1;
    col[0] = xloctglob; col[1] = xloctglob+1; col[2] = xlocfglob; col[3] = xlocfglob+1;
    val[0] = dSt2_dVrt;
    val[1] = dSt2_dVit;
    val[2] = dSt2_dVrf;
    val[3] = dSt2_dVif;
    ierr = MatSetValues(Ji,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

    gloc += 2;
  }

  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(Ji,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Ji,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeConstraints_PBCAR(OPFLOW opflow,Vec X,Vec G)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeObjective_PBCAR(OPFLOW opflow,Vec X,PetscScalar *obj)
{
  PetscErrorCode ierr;
  const PetscScalar *x;
  PS             ps=opflow->ps;
  PetscInt       i;
  PSBUS          bus;
  PSGEN          gen;
  PetscInt       loc;
  PetscBool      isghost;
  Vec            localX;
  PetscReal      localobj;

  PetscFunctionBegin;

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);

  localobj = 0.0;
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);
    if(isghost) continue;

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    
    PetscInt k;
    PetscScalar Pg;
    for(k=0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      loc = loc+2;
      Pg = x[loc]*ps->MVAbase;
      localobj += gen->cost_alpha*Pg*Pg + gen->cost_beta*Pg + gen->cost_gamma;
    }

    if(opflow->include_loadloss_variables) {
      PSLOAD load;
      PetscScalar Pdloss,Qdloss;
      for(k=0; k < bus->nload; k++) {
	loc += 2;
	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
	if(!load->status) continue;
	Pdloss = x[loc];
	Qdloss = x[loc+1];
	localobj += opflow->loadloss_penalty*ps->MVAbase*ps->MVAbase*(Pdloss*Pdloss + Qdloss*Qdloss);
      }
    }

    if(opflow->include_powerimbalance_variables) {
      PetscScalar Pimb,Qimb;
      loc += 2;
      Pimb = x[loc];
      Qimb = x[loc+1];
      localobj += opflow->powerimbalance_penalty*ps->MVAbase*ps->MVAbase*(Pimb*Pimb + Qimb*Qimb);
    }
  }
  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

  ierr = MPI_Allreduce(&localobj,obj,1,MPIU_REAL,MPI_SUM,opflow->comm->type);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeGradient_PBCAR(OPFLOW opflow,Vec X,Vec grad)
{
  PetscErrorCode ierr;
  const PetscScalar *x;
  PetscScalar    *df;
  PS             ps=opflow->ps;
  PetscInt       i;
  PSBUS          bus;
  PSGEN          gen;
  PetscInt       loc;
  Vec            localX,localgrad;
  PetscBool      isghost;

  PetscFunctionBegin;

  ierr = VecSet(grad,0.0);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ps->networkdm,&localgrad);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,grad,INSERT_VALUES,localgrad);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,grad,INSERT_VALUES,localgrad);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = VecGetArray(localgrad,&df);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);
    if(isghost) continue;

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    
    PetscInt k;
    PetscScalar Pg;
    for(k=0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      loc = loc+2;
      Pg = x[loc]*ps->MVAbase;
      df[loc] = ps->MVAbase*(2*gen->cost_alpha*Pg + gen->cost_beta);
    }

    if(opflow->include_loadloss_variables) {
      PSLOAD load;
      PetscScalar Pdloss,Qdloss;
      for(k=0; k < bus->nload; k++) {
	loc += 2;
	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
	if(!load->status) continue;
	Pdloss = x[loc];
	Qdloss = x[loc+1];
	df[loc]   = opflow->loadloss_penalty*ps->MVAbase*ps->MVAbase*2*Pdloss;
	df[loc+1] = opflow->loadloss_penalty*ps->MVAbase*ps->MVAbase*2*Qdloss; 
      }
    }

    if(opflow->include_powerimbalance_variables) {
      PetscScalar Pimb,Qimb;
      loc += 2;
      Pimb = x[loc];
      Qimb = x[loc+1];
      df[loc] = opflow->powerimbalance_penalty*ps->MVAbase*ps->MVAbase*2*Pimb;
      df[loc+1] = opflow->powerimbalance_penalty*ps->MVAbase*ps->MVAbase*2*Qimb;
    }
  }
  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(localgrad,&df);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(ps->networkdm,localgrad,INSERT_VALUES,grad);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(ps->networkdm,localgrad,INSERT_VALUES,grad);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localgrad);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeObjandGradient_PBCAR(OPFLOW opflow,Vec X,PetscScalar *obj,Vec Grad)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = OPFLOWComputeObjective_PBCAR(opflow,X,obj);CHKERRQ(ierr);
  ierr = OPFLOWComputeGradient_PBCAR(opflow,X,Grad);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWFormulationSetNumVariables_PBCAR(OPFLOW opflow,PetscInt *busnvar,PetscInt *branchnvar,PetscInt *nx)
{
  PetscInt i,ngen,nload,j;
  PS       ps=opflow->ps;
  PSBUS    bus;
  PSGEN    gen;
  PetscErrorCode ierr;
  PetscBool      isghost;

  PetscFunctionBegin;
  
  *nx = 0;
  /* No variables for the branches */
  for(i=0; i < ps->nbranch; i++) {
    branchnvar[i] = 0;
    *nx += branchnvar[i];
  }

  /* Variables for the buses */
  /* Real and imaginary part of voltages at each bus + Pg,Qg for each gen */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);
    busnvar[i] = 2; /* 2 variables for voltage */
    ierr = PSBUSGetNGen(bus,&ngen);CHKERRQ(ierr);
    for(j=0; j < ngen; j++) {
      ierr = PSBUSGetGen(bus,j,&gen);CHKERRQ(ierr);
      if(gen->status) busnvar[i] += 2; /* (2 variables Pg, Qg for each in-service gen) */
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

PetscErrorCode OPFLOWFormulationSetNumConstraints_PBCAR(OPFLOW opflow,PetscInt *branchnconeq,PetscInt *busnconeq,PetscInt *nconeq,PetscInt *nconineq)
{
  PetscInt i;
  PS       ps=opflow->ps;
  PSBUS    bus;
  PSLINE   line;
  PetscErrorCode ierr;
  PetscBool      isghost;

  PetscFunctionBegin;
  
  *nconeq = 0;
  *nconineq = 0;

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);
    busnconeq[i] = 2; /* Power balance constraints */
    if(bus->ide == REF_BUS) busnconeq[i] += 1;
    if(!isghost) {
      *nconeq += 2;
      if(bus->ide == REF_BUS) *nconeq += 1; /* Reference angle constraint */
      *nconineq += 1; /* Voltage magnitude constraint */
    }
  }

  for(i=0; i < ps->nbranch; i++) {
    line = &ps->line[i];
    if(line->status && line->rateA < 1e5) *nconineq += 2; /* Line flow constraints */
  }

  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeEqualityConstraintsHessian - Computes the Hessian for the equality constraints function part
  
  Input Parameters:
+ opflow   - the OPFLOW object
. X        - solution vector X
- Lambda   - Lagrangian multiplier vector

  Output Parameters:
. H - the Hessian part for the equality constraints

*/
PetscErrorCode OPFLOWComputeEqualityConstraintsHessian_PBCAR(OPFLOW opflow,Vec X,Vec Lambda,Mat H) 
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscInt       i,k;
  PSLINE         line;
  PSBUS          bus;
  PetscInt       nconnlines;
  const PSLINE   *connlines;
  const PSBUS    *connbuses;
  PetscInt       xloc,xlocf,xloct,xlocglob,xlocfglob,xloctglob;
  PSBUS          busf,bust;
  const PetscScalar *x;
  const PetscScalar *lambda;
  PetscInt       gloc=0;
  PetscInt       row[12],col[12];
  PetscScalar    val[12];
  Vec            localX;
  PetscBool      isghost;

  PetscFunctionBegin;

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecScatterBegin(opflow->scattereqcon,Lambda,opflow->Lambdaelocal,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd(opflow->scattereqcon,Lambda,opflow->Lambdaelocal,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(opflow->Lambdaelocal,&lambda);CHKERRQ(ierr);

  // For equality constraints (power flow) */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);
    
    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bus,&xlocglob);CHKERRQ(ierr);
    
    if(!isghost) {
      /* Shunt */
      row[0] = xlocglob; row[1] = xlocglob+1;
      col[0] = xlocglob; col[1] = xlocglob+1;
      val[0] = lambda[gloc]*2*bus->gl + lambda[gloc+1]*(-2*bus->bl);
      val[3] = val[0];
      val[1] = val[2] = 0.0;
      ierr = MatSetValues(H,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
    }

    ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);
    
    for(k=0; k < nconnlines; k++) {
      line = connlines[k];
      
      if(!line->status) continue;
      
      PetscScalar Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
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

      ierr = PSBUSGetVariableGlobalLocation(busf,&xlocfglob);CHKERRQ(ierr);
      ierr = PSBUSGetVariableGlobalLocation(bust,&xloctglob);CHKERRQ(ierr);
     
      if(bus == busf) {
	 
	PetscScalar dPf_dVrf_dVrf,dPf_dVrf_dVif,dPf_dVrf_dVrt,dPf_dVrf_dVit;
	PetscScalar dPf_dVif_dVrf,dPf_dVif_dVif,dPf_dVif_dVrt,dPf_dVif_dVit;
	PetscScalar dPf_dVrt_dVrf,dPf_dVrt_dVif,dPf_dVrt_dVrt,dPf_dVrt_dVit;
	PetscScalar dPf_dVit_dVrf,dPf_dVit_dVif,dPf_dVit_dVrt,dPf_dVit_dVit;

	/* Partial of     dPf_dVrf = 2*Gff*Vrf + (Gft*Vrt - Bft*Vit); */
	dPf_dVrf_dVrf = 2*Gff;
	dPf_dVrf_dVif = 0.0;
	dPf_dVrf_dVrt = Gft;
	dPf_dVrf_dVit = -Bft;

	/* Partial of dPf_dVif = 2*Gff*Vif + (Bft*Vrt + Gft*Vit); */
	dPf_dVif_dVrf = 0.0;
	dPf_dVif_dVif = 2*Gff;
	dPf_dVif_dVrt = Bft;
	dPf_dVif_dVit = Gft;

	/* Partial of dPf_dVrt = Vrf*Gft + Vif*Bft; */
	dPf_dVrt_dVrf = Gft;
	dPf_dVrt_dVif = Bft;
	dPf_dVrt_dVrt = 0.0;
	dPf_dVrt_dVit = 0.0;

	/* Partial of dPf_dVit = -Vrf*Bft + Vif*Gft; */
	dPf_dVit_dVrf = -Bft;
	dPf_dVit_dVif =  Gft;
	dPf_dVit_dVrt = 0.0;
	dPf_dVit_dVit = 0.0;

	PetscScalar dQf_dVrf_dVrf,dQf_dVrf_dVif,dQf_dVrf_dVrt,dQf_dVrf_dVit;
	PetscScalar dQf_dVif_dVrf,dQf_dVif_dVif,dQf_dVif_dVrt,dQf_dVif_dVit;
	PetscScalar dQf_dVrt_dVrf,dQf_dVrt_dVif,dQf_dVrt_dVrt,dQf_dVrt_dVit;
	PetscScalar dQf_dVit_dVrf,dQf_dVit_dVif,dQf_dVit_dVrt,dQf_dVit_dVit;

	/* Partial of     dQf_dVrf = -2*Bff*Vrf - (Bft*Vrt + Gft*Vit) */
	dQf_dVrf_dVrf = -2*Bff;
	dQf_dVrf_dVif = 0.0;
	dQf_dVrf_dVrt = -Bft;
	dQf_dVrf_dVit = -Gft;

	/* Partial of     dQf_dVif = -2*Bff*Vif + (Gft*Vrt - Bft*Vit); */
	dQf_dVif_dVrf = 0.0;
	dQf_dVif_dVif = -2*Bff;
	dQf_dVif_dVrt = Gft;
	dQf_dVif_dVit = -Bft;

	/* Partial of dQf_dVrt =  Vif*Gft - Vrf*Bft; */
	dQf_dVrt_dVrf = -Bft;
	dQf_dVrt_dVif = Gft;
	dQf_dVrt_dVrt = 0.0;
	dQf_dVrt_dVit = 0.0;

	/* Partial of dQf_dVit = -Vif*Bft - Vrf*Gft; */
	dQf_dVit_dVrf = -Gft;
	dQf_dVit_dVif = -Bft;
	dQf_dVit_dVrt = 0.0;
	dQf_dVit_dVit = 0.0;

	row[0] 	= xlocfglob; 	
	row[1] 	= xlocfglob+1;
	row[2]  = xloctglob;
	row[3]  = xloctglob+1;

	col[0] 	= xlocfglob;

	val[0] = val[1] = val[2] = val[3] = 0.0;
	val[0]  = lambda[gloc]*dPf_dVrf_dVrf + lambda[gloc+1]*dQf_dVrf_dVrf;
	val[1]  = lambda[gloc]*dPf_dVif_dVrf + lambda[gloc+1]*dQf_dVif_dVrf;
	val[2]  = lambda[gloc]*dPf_dVrt_dVrf + lambda[gloc+1]*dQf_dVrt_dVrf;
	val[3]  = lambda[gloc]*dPf_dVit_dVrf + lambda[gloc+1]*dQf_dVit_dVrf;
	ierr = MatSetValues(H,4,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	col[0] 	= xlocfglob+1;

	val[0] = val[1] = val[2] = val[3] = 0.0;
	val[0]  = lambda[gloc]*dPf_dVrf_dVif + lambda[gloc+1]*dQf_dVrf_dVif;
	val[1]  = lambda[gloc]*dPf_dVif_dVif + lambda[gloc+1]*dQf_dVif_dVif;
	val[2]  = lambda[gloc]*dPf_dVrt_dVif + lambda[gloc+1]*dQf_dVrt_dVif;
	val[3]  = lambda[gloc]*dPf_dVit_dVif + lambda[gloc+1]*dQf_dVit_dVif;

	ierr = MatSetValues(H,4,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	col[0] 	= xloctglob;

	val[0] = val[1] = val[2] = val[3] = 0.0;
	val[0]  = lambda[gloc]*dPf_dVrf_dVrt + lambda[gloc+1]*dQf_dVrf_dVrt;
	val[1]  = lambda[gloc]*dPf_dVif_dVrt + lambda[gloc+1]*dQf_dVif_dVrt;
	val[2]  = lambda[gloc]*dPf_dVrt_dVrt + lambda[gloc+1]*dQf_dVrt_dVrt;
	val[3]  = lambda[gloc]*dPf_dVit_dVrt + lambda[gloc+1]*dQf_dVit_dVrt;
	ierr = MatSetValues(H,4,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	col[0] 	= xloctglob+1;

	val[0] = val[1] = val[2] = val[3] = 0.0;
	val[0]  = lambda[gloc]*dPf_dVrf_dVit + lambda[gloc+1]*dQf_dVrf_dVit;
	val[1]  = lambda[gloc]*dPf_dVif_dVit + lambda[gloc+1]*dQf_dVif_dVit;
	val[2]  = lambda[gloc]*dPf_dVrt_dVit + lambda[gloc+1]*dQf_dVrt_dVit;
	val[3]  = lambda[gloc]*dPf_dVit_dVit + lambda[gloc+1]*dQf_dVit_dVit;
	ierr = MatSetValues(H,4,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
      } else {
	PetscScalar dPt_dVrf_dVrf,dPt_dVrf_dVif,dPt_dVrf_dVrt,dPt_dVrf_dVit;
	PetscScalar dPt_dVif_dVrf,dPt_dVif_dVif,dPt_dVif_dVrt,dPt_dVif_dVit;
	PetscScalar dPt_dVrt_dVrf,dPt_dVrt_dVif,dPt_dVrt_dVrt,dPt_dVrt_dVit;
	PetscScalar dPt_dVit_dVrf,dPt_dVit_dVif,dPt_dVit_dVrt,dPt_dVit_dVit;

	/* Partial of     dPt_dVrt = 2*Gtt*Vrt + (Gtf*Vrf - Btf*Vif); */
	dPt_dVrt_dVrt = 2*Gtt;
	dPt_dVrt_dVit = 0.0;
	dPt_dVrt_dVrf = Gtf;
	dPt_dVrt_dVif = -Btf;

	/* Partial of dPt_dVit = 2*Gtt*Vit + (Btf*Vrf + Gtf*Vif); */
	dPt_dVit_dVrt = 0.0;
	dPt_dVit_dVit = 2*Gtt;
	dPt_dVit_dVrf = Btf;
	dPt_dVit_dVif = Gtf;

	/* Partial of dPt_dVrf = Vrt*Gtf + Vit*Btf; */
	dPt_dVrf_dVrt = Gtf;
	dPt_dVrf_dVit = Btf;
	dPt_dVrf_dVrf = 0.0;
	dPt_dVrf_dVif = 0.0;

	/* Partial of dPt_dVif = -Vrt*Btf + Vit*Gtf; */
	dPt_dVif_dVrt = -Btf;
	dPt_dVif_dVit =  Gtf;
	dPt_dVif_dVrf = 0.0;
	dPt_dVif_dVif = 0.0;
	
	PetscScalar dQt_dVrf_dVrf,dQt_dVrf_dVif,dQt_dVrf_dVrt,dQt_dVrf_dVit;
	PetscScalar dQt_dVif_dVrf,dQt_dVif_dVif,dQt_dVif_dVrt,dQt_dVif_dVit;
	PetscScalar dQt_dVrt_dVrf,dQt_dVrt_dVif,dQt_dVrt_dVrt,dQt_dVrt_dVit;
	PetscScalar dQt_dVit_dVrf,dQt_dVit_dVif,dQt_dVit_dVrt,dQt_dVit_dVit;

	/* Partial of dQt_dVrt = -2*Btt*Vrt - (Btf*Vrf + Gtf*Vif) */
	dQt_dVrt_dVrt = -2*Btt;
	dQt_dVrt_dVit = 0.0;
	dQt_dVrt_dVrf = -Btf;
	dQt_dVrt_dVif = -Gtf;

	/* Partial of dQt_dVit = -2*Btt*Vit + (Gtf*Vrf - Btf*Vif); */
	dQt_dVit_dVrt = 0.0;
	dQt_dVit_dVit = -2*Btt;
	dQt_dVit_dVrf = Gtf;
	dQt_dVit_dVif = -Btf;

	/* Partial of dQt_dVrf =  Vit*Gtf - Vrt*Btf; */
	dQt_dVrf_dVrt = -Btf;
	dQt_dVrf_dVit = Gtf;
	dQt_dVrf_dVrf = 0.0;
	dQt_dVrf_dVif = 0.0;

	/* Partial of dQt_dVif = -Vit*Btf - Vrt*Gtf; */
	dQt_dVif_dVrt = -Gtf;
	dQt_dVif_dVit = -Btf;
	dQt_dVif_dVrf =  0.0;
	dQt_dVif_dVif =  0.0;

	row[0] 	= xloctglob; 	
	row[1] 	= xloctglob+1;
	row[2]  = xlocfglob;
	row[3]  = xlocfglob+1;

	col[0] 	= xloctglob;

	val[0] = val[1] = val[2] = val[3] = 0.0;
	val[0]  = lambda[gloc]*dPt_dVrt_dVrt + lambda[gloc+1]*dQt_dVrt_dVrt;
	val[1]  = lambda[gloc]*dPt_dVit_dVrt + lambda[gloc+1]*dQt_dVit_dVrt;
	val[2]  = lambda[gloc]*dPt_dVrf_dVrt + lambda[gloc+1]*dQt_dVrf_dVrt;
	val[3]  = lambda[gloc]*dPt_dVif_dVrt + lambda[gloc+1]*dQt_dVif_dVrt;
	ierr = MatSetValues(H,4,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	col[0] 	= xloctglob+1;

	val[0] = val[1] = val[2] = val[3] = 0.0;
	val[0]  = lambda[gloc]*dPt_dVrt_dVit + lambda[gloc+1]*dQt_dVrt_dVit;
	val[1]  = lambda[gloc]*dPt_dVit_dVit + lambda[gloc+1]*dQt_dVit_dVit;
	val[2]  = lambda[gloc]*dPt_dVrf_dVit + lambda[gloc+1]*dQt_dVrf_dVit;
	val[3]  = lambda[gloc]*dPt_dVif_dVit + lambda[gloc+1]*dQt_dVif_dVit;
	ierr = MatSetValues(H,4,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	col[0] 	= xlocfglob;

	val[0] = val[1] = val[2] = val[3] = 0.0;
	val[0]  = lambda[gloc]*dPt_dVrt_dVrf + lambda[gloc+1]*dQt_dVrt_dVrf;
	val[1]  = lambda[gloc]*dPt_dVit_dVrf + lambda[gloc+1]*dQt_dVit_dVrf;
	val[2]  = lambda[gloc]*dPt_dVrf_dVrf + lambda[gloc+1]*dQt_dVrf_dVrf;
	val[3]  = lambda[gloc]*dPt_dVif_dVrf + lambda[gloc+1]*dQt_dVif_dVrf;
	ierr = MatSetValues(H,4,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	col[0] 	= xlocfglob+1;

	val[0] = val[1] = val[2] = val[3] = 0.0;
	val[0]  = lambda[gloc]*dPt_dVrt_dVif + lambda[gloc+1]*dQt_dVrt_dVif;
	val[1]  = lambda[gloc]*dPt_dVit_dVif + lambda[gloc+1]*dQt_dVit_dVif;
	val[2]  = lambda[gloc]*dPt_dVrf_dVif + lambda[gloc+1]*dQt_dVrf_dVif;
	val[3]  = lambda[gloc]*dPt_dVif_dVif + lambda[gloc+1]*dQt_dVif_dVif;
	ierr = MatSetValues(H,4,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
      }
    }
    gloc += 2;
    if(bus->ide == REF_BUS) gloc++;
  }

  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(opflow->Lambdaelocal,&lambda);CHKERRQ(ierr);
  
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*
  OPFLOWComputeInequalityConstraintsHessian - Computes the Inequality Constraints Hessian

  Input Parameters:
+ opflow   - the OPFLOW object
. X        - the solution vector
- Lambda   - Lagrangian multipler vector

  Output Parameters:
+ H   - the Hessian matrix

*/
PetscErrorCode OPFLOWComputeInequalityConstraintsHessian_PBCAR(OPFLOW opflow, Vec X, Vec Lambda,Mat H)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscInt       i;
  PSLINE         line;
  const PSBUS    *connbuses;
  PetscInt       xloc,xlocf,xloct,xlocglob,xlocglobf,xlocglobt;
  PSBUS          bus,busf,bust;
  const PetscScalar *x;
  const PetscScalar *lambda;
  PetscInt       gloc=0;
  PetscInt       row[12],col[12];
  PetscScalar    val[12];
  PetscScalar    Vrf,Vif,Vrt,Vit;
  PetscScalar    dVm2_dVr2,dVm2_dVi2;
  PetscScalar    Pf,Qf,Pt,Qt;
  Vec            localX;
  PetscBool      isghost;
  PetscInt       rstart;

  PetscFunctionBegin;

  ierr = MatGetOwnershipRange(H,&rstart,NULL);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(Lambda,&lambda);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);

    if(isghost) continue;

    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bus,&xlocglob);CHKERRQ(ierr);

    dVm2_dVr2 = 2.0;
    dVm2_dVi2 = 2.0;
   
    row[0] = xlocglob; row[1] = xlocglob+1;
    col[0] = xlocglob; col[1] = xlocglob+1;
    val[0] = lambda[gloc]*dVm2_dVr2;
    val[1] = val[2] = 0.0;
    val[3] = lambda[gloc]*dVm2_dVi2;

    ierr = MatSetValues(H,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
    gloc++;
  }

  val[0] = val[1] = val[2] = val[3] = 0.0;
  // for the part of line constraints
  for(i=0; i < ps->nbranch; i++) {
    line = &ps->line[i];
    if(!line->status || line->rateA > 1e5) continue;

    PetscScalar Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
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

    ierr = PSBUSGetVariableGlobalLocation(busf,&xlocglobf);CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bust,&xlocglobt);CHKERRQ(ierr);

    Vrf  = x[xlocf];
    Vif  = x[xlocf+1];
    Vrt  = x[xloct];
    Vit  = x[xloct+1];

    Pf =  Gff*(Vrf*Vrf + Vif*Vif) + Vrf*(Gft*Vrt - Bft*Vit) + Vif*(Bft*Vrt + Gft*Vit);
    Qf = -Bff*(Vrf*Vrf + Vif*Vif) + Vif*(Gft*Vrt - Bft*Vit) - Vrf*(Bft*Vrt + Gft*Vit);

    Pt =  Gtt*(Vrt*Vrt + Vit*Vit) + Vrt*(Gtf*Vrf - Btf*Vif) + Vit*(Btf*Vrf + Gtf*Vif);
    Qt = -Btt*(Vrt*Vrt + Vit*Vit) + Vit*(Gtf*Vrf - Btf*Vif) - Vrt*(Btf*Vrf + Gtf*Vif);

    PetscScalar dSf2_dPf,dSf2_dQf,dSt2_dPt,dSt2_dQt;
    dSf2_dPf = 2*Pf;
    dSf2_dQf = 2*Qf;
    dSt2_dPt = 2*Pt;
    dSt2_dQt = 2*Qt;

    /* Need to calculat the partial derivatives w.r.t. these terms
    dSf2_dVrf = dSf2_dPf*dPf_dVrf + dSf2_dQf*dQf_dVrf;
    dSf2_dVrt = dSf2_dPf*dPf_dVrt + dSf2_dQf*dQf_dVrt;
    dSf2_dVif = dSf2_dPf*dPf_dVif + dSf2_dQf*dQf_dVif;
    dSf2_dVit = dSf2_dPf*dPf_dVit + dSf2_dQf*dQf_dVit;

    which would be

    dSf2_dVrf_dVrf = dSf2_dPf_dVrf*dPf_dVrf + dSf2_dPf*dPf_dVrf_Vrf + dSf2_dQf_dVrf*dQf_dVrf + dSf2_dQf*dQf_dVrf_dVrf
     = 2*dPf_dVr*dPf_dVrf + dSf2_dPf*dPf_dVrf_dVrf + 2*dQf_dVrf*dQf_dVrf + dSf2_dQf*dQf_dVrf_dVrf

    dSf2_dVrf_dVif = dSf2_dPf_dVif*dPf_dVrf + dSf2_dPf*dPf_dVrf_dVif
                     + dSf2_dPf_dVif*dPf_dVrf + dSf2_dPf*dPf_dVrf_dVif
    */
    PetscScalar    dPf_dVrf,dPf_dVif,dPf_dVrt,dPf_dVit;
    PetscScalar    dQf_dVrf,dQf_dVif,dQf_dVrt,dQf_dVit;
    PetscScalar    dPt_dVrf,dPt_dVif,dPt_dVrt,dPt_dVit;
    PetscScalar    dQt_dVrf,dQt_dVif,dQt_dVrt,dQt_dVit;

    dPf_dVrf = 2*Gff*Vrf + (Gft*Vrt - Bft*Vit);
    dPf_dVif = 2*Gff*Vif + (Bft*Vrt + Gft*Vit);
    dPf_dVrt = Vrf*Gft + Vif*Bft;
    dPf_dVit = -Vrf*Bft + Vif*Gft;

    dQf_dVrf = -2*Bff*Vrf - (Bft*Vrt + Gft*Vit);
    dQf_dVif = -2*Bff*Vif + (Gft*Vrt - Bft*Vit);
    dQf_dVrt =  Vif*Gft - Vrf*Bft;
    dQf_dVit = -Vif*Bft - Vrf*Gft;

    dPt_dVrt = 2*Gtt*Vrt + (Gtf*Vrf - Btf*Vif);
    dPt_dVit = 2*Gtt*Vit + (Btf*Vrf + Gtf*Vif);
    dPt_dVrf = Vrt*Gtf  + Vit*Btf;
    dPt_dVif = -Vrt*Btf + Vit*Gtf;

    dQt_dVrt = -2*Btt*Vrt - (Btf*Vrf + Gtf*Vif);
    dQt_dVit = -2*Btt*Vit + (Gtf*Vrf - Btf*Vif);
    dQt_dVrf =  Vit*Gtf - Vrt*Btf;
    dQt_dVif  = -Vit*Btf - Vrt*Gtf;

    PetscScalar d2Pf_dVrf_dVrf,d2Pf_dVrf_dVif,d2Pf_dVrf_dVrt,d2Pf_dVrf_dVit;
    PetscScalar d2Pf_dVif_dVrf,d2Pf_dVif_dVif,d2Pf_dVif_dVrt,d2Pf_dVif_dVit;
    PetscScalar d2Pf_dVrt_dVrf,d2Pf_dVrt_dVif,d2Pf_dVrt_dVrt,d2Pf_dVrt_dVit;
    PetscScalar d2Pf_dVit_dVrf,d2Pf_dVit_dVif,d2Pf_dVit_dVrt,d2Pf_dVit_dVit;

    /* Partial of     dPf_dVrf = 2*Gff*Vrf + (Gft*Vrt - Bft*Vit); */
    d2Pf_dVrf_dVrf = 2*Gff;
    d2Pf_dVrf_dVif = 0.0;
    d2Pf_dVrf_dVrt = Gft;
    d2Pf_dVrf_dVit = -Bft;

    /* Partial of dPf_dVif = 2*Gff*Vif + (Bft*Vrt + Gft*Vit); */
    d2Pf_dVif_dVrf = 0.0;
    d2Pf_dVif_dVif = 2*Gff;
    d2Pf_dVif_dVrt = Bft;
    d2Pf_dVif_dVit = Gft;
    
    /* Partial of dPf_dVrt = Vrf*Gft + Vif*Bft; */
    d2Pf_dVrt_dVrf = Gft;
    d2Pf_dVrt_dVif = Bft;
    d2Pf_dVrt_dVrt = 0.0;
    d2Pf_dVrt_dVit = 0.0;
    
    /* Partial of dPf_dVit = -Vrf*Bft + Vif*Gft; */
    d2Pf_dVit_dVrf = -Bft;
    d2Pf_dVit_dVif =  Gft;
    d2Pf_dVit_dVrt = 0.0;
    d2Pf_dVit_dVit = 0.0;
    
    PetscScalar d2Qf_dVrf_dVrf,d2Qf_dVrf_dVif,d2Qf_dVrf_dVrt,d2Qf_dVrf_dVit;
    PetscScalar d2Qf_dVif_dVrf,d2Qf_dVif_dVif,d2Qf_dVif_dVrt,d2Qf_dVif_dVit;
    PetscScalar d2Qf_dVrt_dVrf,d2Qf_dVrt_dVif,d2Qf_dVrt_dVrt,d2Qf_dVrt_dVit;
    PetscScalar d2Qf_dVit_dVrf,d2Qf_dVit_dVif,d2Qf_dVit_dVrt,d2Qf_dVit_dVit;
    
    /* Partial of     dQf_dVrf = -2*Bff*Vrf - (Bft*Vrt + Gft*Vit) */
    d2Qf_dVrf_dVrf = -2*Bff;
    d2Qf_dVrf_dVif = 0.0;
    d2Qf_dVrf_dVrt = -Bft;
    d2Qf_dVrf_dVit = -Gft;
    
    /* Partial of     dQf_dVif = -2*Bff*Vif + (Gft*Vrt - Bft*Vit); */
    d2Qf_dVif_dVrf = 0.0;
    d2Qf_dVif_dVif = -2*Bff;
    d2Qf_dVif_dVrt = Gft;
    d2Qf_dVif_dVit = -Bft;
    
    /* Partial of dQf_dVrt =  Vif*Gft - Vrf*Bft; */
    d2Qf_dVrt_dVrf = -Bft;
    d2Qf_dVrt_dVif = Gft;
    d2Qf_dVrt_dVrt = 0.0;
    d2Qf_dVrt_dVit = 0.0;
    
    /* Partial of dQf_dVit = -Vif*Bft - Vrf*Gft; */
    d2Qf_dVit_dVrf = -Gft;
    d2Qf_dVit_dVif = -Bft;
    d2Qf_dVit_dVrt = 0.0;
    d2Qf_dVit_dVit = 0.0;

    PetscScalar d2Pt_dVrf_dVrf,d2Pt_dVrf_dVif,d2Pt_dVrf_dVrt,d2Pt_dVrf_dVit;
    PetscScalar d2Pt_dVif_dVrf,d2Pt_dVif_dVif,d2Pt_dVif_dVrt,d2Pt_dVif_dVit;
    PetscScalar d2Pt_dVrt_dVrf,d2Pt_dVrt_dVif,d2Pt_dVrt_dVrt,d2Pt_dVrt_dVit;
    PetscScalar d2Pt_dVit_dVrf,d2Pt_dVit_dVif,d2Pt_dVit_dVrt,d2Pt_dVit_dVit;
    
    /* Partial of     dPt_dVrt = 2*Gtt*Vrt + (Gtf*Vrf - Btf*Vif); */
    d2Pt_dVrt_dVrt = 2*Gtt;
    d2Pt_dVrt_dVit = 0.0;
    d2Pt_dVrt_dVrf = Gtf;
    d2Pt_dVrt_dVif = -Btf;
    
    /* Partial of dPt_dVit = 2*Gtt*Vit + (Btf*Vrf + Gtf*Vif); */
    d2Pt_dVit_dVrt = 0.0;
    d2Pt_dVit_dVit = 2*Gtt;
    d2Pt_dVit_dVrf = Btf;
    d2Pt_dVit_dVif = Gtf;
    
    /* Partial of dPt_dVrf = Vrt*Gtf + Vit*Btf; */
    d2Pt_dVrf_dVrt = Gtf;
    d2Pt_dVrf_dVit = Btf;
    d2Pt_dVrf_dVrf = 0.0;
    d2Pt_dVrf_dVif = 0.0;
    
    /* Partial of dPt_dVif = -Vrt*Btf + Vit*Gtf; */
    d2Pt_dVif_dVrt = -Btf;
    d2Pt_dVif_dVit =  Gtf;
    d2Pt_dVif_dVrf = 0.0;
    d2Pt_dVif_dVif = 0.0;
    
    PetscScalar d2Qt_dVrf_dVrf,d2Qt_dVrf_dVif,d2Qt_dVrf_dVrt,d2Qt_dVrf_dVit;
    PetscScalar d2Qt_dVif_dVrf,d2Qt_dVif_dVif,d2Qt_dVif_dVrt,d2Qt_dVif_dVit;
    PetscScalar d2Qt_dVrt_dVrf,d2Qt_dVrt_dVif,d2Qt_dVrt_dVrt,d2Qt_dVrt_dVit;
    PetscScalar d2Qt_dVit_dVrf,d2Qt_dVit_dVif,d2Qt_dVit_dVrt,d2Qt_dVit_dVit;
    
    /* Partial of dQt_dVrt = -2*Btt*Vrt - (Btf*Vrf + Gtf*Vif) */
    d2Qt_dVrt_dVrt = -2*Btt;
    d2Qt_dVrt_dVit = 0.0;
    d2Qt_dVrt_dVrf = -Btf;
    d2Qt_dVrt_dVif = -Gtf;
    
    /* Partial of dQt_dVit = -2*Btt*Vit + (Gtf*Vrf - Btf*Vif); */
    d2Qt_dVit_dVrt = 0.0;
    d2Qt_dVit_dVit = -2*Btt;
    d2Qt_dVit_dVrf = Gtf;
    d2Qt_dVit_dVif = -Btf;
    
    /* Partial of dQt_dVrf =  Vit*Gtf - Vrt*Btf; */
    d2Qt_dVrf_dVrt = -Btf;
    d2Qt_dVrf_dVit = Gtf;
    d2Qt_dVrf_dVrf = 0.0;
    d2Qt_dVrf_dVif = 0.0;
    
    /* Partial of dQt_dVif = -Vit*Btf - Vrt*Gtf; */
    d2Qt_dVif_dVrt = -Gtf;
    d2Qt_dVif_dVit = -Btf;
    d2Qt_dVif_dVrf =  0.0;
    d2Qt_dVif_dVif =  0.0;
    
    PetscScalar d2Sf2_dVrf_dVrf=0.0,d2Sf2_dVrf_dVif=0.0,d2Sf2_dVrf_dVrt=0.0,d2Sf2_dVrf_dVit=0.0;
    PetscScalar d2St2_dVrf_dVrf=0.0,d2St2_dVrf_dVif=0.0,d2St2_dVrf_dVrt=0.0,d2St2_dVrf_dVit=0.0;

    d2Sf2_dVrf_dVrf = 2*dPf_dVrf*dPf_dVrf + dSf2_dPf*d2Pf_dVrf_dVrf +  2*dQf_dVrf*dQf_dVrf + dSf2_dQf*d2Qf_dVrf_dVrf;
    d2Sf2_dVrf_dVif = 2*dPf_dVif*dPf_dVrf + dSf2_dPf*d2Pf_dVrf_dVif +  2*dQf_dVif*dQf_dVrf + dSf2_dQf*d2Qf_dVrf_dVif;
    d2Sf2_dVrf_dVrt = 2*dPf_dVrt*dPf_dVrf + dSf2_dPf*d2Pf_dVrf_dVrt +  2*dQf_dVrt*dQf_dVrf + dSf2_dQf*d2Qf_dVrf_dVrt;
    d2Sf2_dVrf_dVit = 2*dPf_dVit*dPf_dVrf + dSf2_dPf*d2Pf_dVrf_dVit +  2*dQf_dVit*dQf_dVrf + dSf2_dQf*d2Qf_dVrf_dVit;

    d2St2_dVrf_dVrf = 2*dPt_dVrf*dPt_dVrf + dSt2_dPt*d2Pt_dVrf_dVrf +  2*dQt_dVrf*dQt_dVrf + dSt2_dQt*d2Qt_dVrf_dVrf;
    d2St2_dVrf_dVif = 2*dPt_dVif*dPt_dVrf + dSt2_dPt*d2Pt_dVrf_dVif +  2*dQt_dVif*dQt_dVrf + dSt2_dQt*d2Qt_dVrf_dVif;
    d2St2_dVrf_dVrt = 2*dPt_dVrt*dPt_dVrf + dSt2_dPt*d2Pt_dVrf_dVrt +  2*dQt_dVrt*dQt_dVrf + dSt2_dQt*d2Qt_dVrf_dVrt;
    d2St2_dVrf_dVit = 2*dPt_dVit*dPt_dVrf + dSt2_dPt*d2Pt_dVrf_dVit +  2*dQt_dVit*dQt_dVrf + dSt2_dQt*d2Qt_dVrf_dVit;

    val[0] = val[1] = val[2] = val[3] = 0.0;
    col[0] = xlocglobf; 
    col[1] = xlocglobf+1; 
    col[2] = xlocglobt;
    col[3] = xlocglobt+1;
    row[0] = xlocglobf;
    val[0] = lambda[gloc]*d2Sf2_dVrf_dVrf + lambda[gloc+1]*d2St2_dVrf_dVrf;
    val[1] = lambda[gloc]*d2Sf2_dVrf_dVif + lambda[gloc+1]*d2St2_dVrf_dVif;
    val[2] = lambda[gloc]*d2Sf2_dVrf_dVrt + lambda[gloc+1]*d2St2_dVrf_dVrt;
    val[3] = lambda[gloc]*d2Sf2_dVrf_dVit + lambda[gloc+1]*d2St2_dVrf_dVit;

    ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

    PetscScalar d2Sf2_dVif_dVrf,d2Sf2_dVif_dVif,d2Sf2_dVif_dVrt,d2Sf2_dVif_dVit;
    PetscScalar d2St2_dVif_dVrf,d2St2_dVif_dVif,d2St2_dVif_dVrt,d2St2_dVif_dVit;

    d2Sf2_dVif_dVrf = 2*dPf_dVrf*dPf_dVif + dSf2_dPf*d2Pf_dVif_dVrf +  2*dQf_dVrf*dQf_dVif + dSf2_dQf*d2Qf_dVif_dVrf;
    d2Sf2_dVif_dVif = 2*dPf_dVif*dPf_dVif + dSf2_dPf*d2Pf_dVif_dVif +  2*dQf_dVif*dQf_dVif + dSf2_dQf*d2Qf_dVif_dVif;
    d2Sf2_dVif_dVrt = 2*dPf_dVrt*dPf_dVif + dSf2_dPf*d2Pf_dVif_dVrt +  2*dQf_dVrt*dQf_dVif + dSf2_dQf*d2Qf_dVif_dVrt;
    d2Sf2_dVif_dVit = 2*dPf_dVit*dPf_dVif + dSf2_dPf*d2Pf_dVif_dVit +  2*dQf_dVit*dQf_dVif + dSf2_dQf*d2Qf_dVif_dVit;

    d2St2_dVif_dVrf = 2*dPt_dVrf*dPt_dVif + dSt2_dPt*d2Pt_dVif_dVrf +  2*dQt_dVrf*dQt_dVif + dSt2_dQt*d2Qt_dVif_dVrf;
    d2St2_dVif_dVif = 2*dPt_dVif*dPt_dVif + dSt2_dPt*d2Pt_dVif_dVif +  2*dQt_dVif*dQt_dVif + dSt2_dQt*d2Qt_dVif_dVif;
    d2St2_dVif_dVrt = 2*dPt_dVrt*dPt_dVif + dSt2_dPt*d2Pt_dVif_dVrt +  2*dQt_dVrt*dQt_dVif + dSt2_dQt*d2Qt_dVif_dVrt;
    d2St2_dVif_dVit = 2*dPt_dVit*dPt_dVif + dSt2_dPt*d2Pt_dVif_dVit +  2*dQt_dVit*dQt_dVif + dSt2_dQt*d2Qt_dVif_dVit;

    val[0] = val[1] = val[2] = val[3] = 0.0;
    col[0] = xlocglobf; 
    col[1] = xlocglobf+1; 
    col[2] = xlocglobt;
    col[3] = xlocglobt+1;
    row[0] = xlocglobf+1;
    val[0] = lambda[gloc]*d2Sf2_dVif_dVrf + lambda[gloc+1]*d2St2_dVif_dVrf;
    val[1] = lambda[gloc]*d2Sf2_dVif_dVif + lambda[gloc+1]*d2St2_dVif_dVif;
    val[2] = lambda[gloc]*d2Sf2_dVif_dVrt + lambda[gloc+1]*d2St2_dVif_dVrt;
    val[3] = lambda[gloc]*d2Sf2_dVif_dVit + lambda[gloc+1]*d2St2_dVif_dVit;
    ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

    PetscScalar d2Sf2_dVrt_dVrf,d2Sf2_dVrt_dVif,d2Sf2_dVrt_dVrt,d2Sf2_dVrt_dVit;
    PetscScalar d2St2_dVrt_dVrf,d2St2_dVrt_dVif,d2St2_dVrt_dVrt,d2St2_dVrt_dVit;

    d2Sf2_dVrt_dVrf = 2*dPf_dVrf*dPf_dVrt + dSf2_dPf*d2Pf_dVrt_dVrf +  2*dQf_dVrt*dQf_dVrf + dSf2_dQf*d2Qf_dVrt_dVrf;
    d2Sf2_dVrt_dVif = 2*dPf_dVif*dPf_dVrt + dSf2_dPf*d2Pf_dVrt_dVif +  2*dQf_dVrt*dQf_dVif + dSf2_dQf*d2Qf_dVrt_dVif;
    d2Sf2_dVrt_dVrt = 2*dPf_dVrt*dPf_dVrt + dSf2_dPf*d2Pf_dVrt_dVrt +  2*dQf_dVrt*dQf_dVrt + dSf2_dQf*d2Qf_dVrt_dVrt;
    d2Sf2_dVrt_dVit = 2*dPf_dVit*dPf_dVrt + dSf2_dPf*d2Pf_dVrt_dVit +  2*dQf_dVrt*dQf_dVit + dSf2_dQf*d2Qf_dVrt_dVit;

    d2St2_dVrt_dVrf = 2*dPt_dVrf*dPt_dVrt + dSt2_dPt*d2Pt_dVrt_dVrf +  2*dQt_dVrf*dQt_dVrt + dSt2_dQt*d2Qt_dVrt_dVrf;
    d2St2_dVrt_dVif = 2*dPt_dVif*dPt_dVrt + dSt2_dPt*d2Pt_dVrt_dVif +  2*dQt_dVif*dQt_dVrt + dSt2_dQt*d2Qt_dVrt_dVif;
    d2St2_dVrt_dVrt = 2*dPt_dVrt*dPt_dVrt + dSt2_dPt*d2Pt_dVrt_dVrt +  2*dQt_dVrt*dQt_dVrt + dSt2_dQt*d2Qt_dVrt_dVrt;
    d2St2_dVrt_dVit = 2*dPt_dVit*dPt_dVrt + dSt2_dPt*d2Pt_dVrt_dVit +  2*dQt_dVit*dQt_dVrt + dSt2_dQt*d2Qt_dVrt_dVit;

    val[0] = val[1] = val[2] = val[3] = 0.0;
    col[0] = xlocglobf; 
    col[1] = xlocglobf+1; 
    col[2] = xlocglobt;
    col[3] = xlocglobt+1;
    row[0] = xlocglobt;

    val[0] = lambda[gloc]*d2Sf2_dVrt_dVrf + lambda[gloc+1]*d2St2_dVrt_dVrf;
    val[1] = lambda[gloc]*d2Sf2_dVrt_dVif + lambda[gloc+1]*d2St2_dVrt_dVif;
    val[2] = lambda[gloc]*d2Sf2_dVrt_dVrt + lambda[gloc+1]*d2St2_dVrt_dVrt;
    val[3] = lambda[gloc]*d2Sf2_dVrt_dVit + lambda[gloc+1]*d2St2_dVrt_dVit;
    
    ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

    PetscScalar d2Sf2_dVit_dVrf,d2Sf2_dVit_dVif,d2Sf2_dVit_dVrt,d2Sf2_dVit_dVit;
    PetscScalar d2St2_dVit_dVrf,d2St2_dVit_dVif,d2St2_dVit_dVrt,d2St2_dVit_dVit;

    d2Sf2_dVit_dVrf = 2*dPf_dVrf*dPf_dVit + dSf2_dPf*d2Pf_dVit_dVrf +  2*dQf_dVrf*dQf_dVit + dSf2_dQf*d2Qf_dVit_dVrf;
    d2Sf2_dVit_dVif = 2*dPf_dVif*dPf_dVit + dSf2_dPf*d2Pf_dVit_dVif +  2*dQf_dVif*dQf_dVit + dSf2_dQf*d2Qf_dVit_dVif;
    d2Sf2_dVit_dVrt = 2*dPf_dVrt*dPf_dVit + dSf2_dPf*d2Pf_dVit_dVrt +  2*dQf_dVrt*dQf_dVit + dSf2_dQf*d2Qf_dVit_dVrt;
    d2Sf2_dVit_dVit = 2*dPf_dVit*dPf_dVit + dSf2_dPf*d2Pf_dVit_dVit +  2*dQf_dVit*dQf_dVit + dSf2_dQf*d2Qf_dVit_dVit;

    d2St2_dVit_dVrf = 2*dPt_dVrf*dPt_dVit + dSt2_dPt*d2Pt_dVit_dVrf +  2*dQt_dVrf*dQt_dVit + dSt2_dQt*d2Qt_dVit_dVrf;
    d2St2_dVit_dVif = 2*dPt_dVif*dPt_dVit + dSt2_dPt*d2Pt_dVit_dVif +  2*dQt_dVif*dQt_dVit + dSt2_dQt*d2Qt_dVit_dVif;
    d2St2_dVit_dVrt = 2*dPt_dVrt*dPt_dVit + dSt2_dPt*d2Pt_dVit_dVrt +  2*dQt_dVrt*dQt_dVit + dSt2_dQt*d2Qt_dVit_dVrt;
    d2St2_dVit_dVit = 2*dPt_dVit*dPt_dVit + dSt2_dPt*d2Pt_dVit_dVit +  2*dQt_dVit*dQt_dVit + dSt2_dQt*d2Qt_dVit_dVit;

    val[0] = val[1] = val[2] = val[3] = 0.0;
    col[0] = xlocglobf; 
    col[1] = xlocglobf+1; 
    col[2] = xlocglobt;
    col[3] = xlocglobt+1;
    row[0] = xlocglobt+1;

    val[0] = lambda[gloc]*d2Sf2_dVit_dVrf + lambda[gloc+1]*d2St2_dVit_dVrf;
    val[1] = lambda[gloc]*d2Sf2_dVit_dVif + lambda[gloc+1]*d2St2_dVit_dVif;
    val[2] = lambda[gloc]*d2Sf2_dVit_dVrt + lambda[gloc+1]*d2St2_dVit_dVrt;
    val[3] = lambda[gloc]*d2Sf2_dVit_dVit + lambda[gloc+1]*d2St2_dVit_dVit;
    
    ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

    gloc +=  2;
  }

  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Lambda,&lambda);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeObjectiveHessian - Computes the Hessian for the objective function part
  
  Input Parameters:
+ opflow - the OPFLOW object
- X        - solution vector X

  Output Parameters:
. H - the Hessian part for the objective function

*/
PetscErrorCode OPFLOWComputeObjectiveHessian_PBCAR(OPFLOW opflow,Vec X,Mat H) 
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscInt       i,k;
  PSBUS          bus;
  PSGEN          gen;
  PetscInt       xloc,xlocglob;
  const PetscScalar *x;
  PetscInt       row[2],col[2];
  PetscScalar    val[2];
  PetscScalar    obj_factor = opflow->obj_factor;
  PetscBool      isghost;
  Vec            localX;


  PetscFunctionBegin;

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);

  // for the part of objective
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);

    if(isghost) continue;

    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bus,&xlocglob);CHKERRQ(ierr);
   
    for(k=0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      xlocglob = xlocglob+2;
      row[0] = xlocglob;
      col[0] = xlocglob;
      val[0] = obj_factor*2.0*gen->cost_alpha*ps->MVAbase*ps->MVAbase;
      ierr = MatSetValues(H,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
    }

    if(opflow->include_loadloss_variables) {
      PSLOAD load;
      for(k=0; k < bus->nload; k++) {
	xlocglob = xlocglob+2;
	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
	if(!load->status) continue;
	row[0] = xlocglob;
	col[0] = xlocglob;
	val[0] = obj_factor*2.0*opflow->loadloss_penalty*ps->MVAbase*ps->MVAbase;
	ierr = MatSetValues(H,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
	
	row[0] = xlocglob+1;
	col[0] = xlocglob+1;
	val[0] = obj_factor*2.0*opflow->loadloss_penalty*ps->MVAbase*ps->MVAbase;
	ierr = MatSetValues(H,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);  
      }
    }

    if(opflow->include_powerimbalance_variables) {
      xlocglob = xlocglob+2;
      row[0] = xlocglob;
      col[0] = xlocglob;
      val[0] = obj_factor*2.0*opflow->powerimbalance_penalty*ps->MVAbase*ps->MVAbase;
      ierr = MatSetValues(H,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
	
      row[0] = xlocglob+1;
      col[0] = xlocglob+1;
      val[0] = obj_factor*2.0*opflow->powerimbalance_penalty*ps->MVAbase*ps->MVAbase;
      ierr = MatSetValues(H,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
    }
  }

  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeHessian_PBCAR(OPFLOW opflow,Vec X,Vec Lambdae,Vec Lambdai,Mat H)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = MatZeroEntries(H);CHKERRQ(ierr);

  /* Objective function Hessian */
  ierr = OPFLOWComputeObjectiveHessian_PBCAR(opflow,X,H);CHKERRQ(ierr);

  /* Equality constraints Hessian */
  ierr = OPFLOWComputeEqualityConstraintsHessian_PBCAR(opflow,X,Lambdae,H);CHKERRQ(ierr);
  
  /* Inequality constraints Hessian */
  if(opflow->nconineq) {
    ierr = OPFLOWComputeInequalityConstraintsHessian_PBCAR(opflow,X,Lambdai,H);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWFormulationCreate_PBCAR(OPFLOW opflow)
{
  PBCAR pbcar;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  
  ierr = PetscCalloc1(1,&pbcar);CHKERRQ(ierr);

  opflow->formulation = pbcar;

  /* Inherit Ops */
  opflow->formops.destroy = OPFLOWFormulationDestroy_PBCAR;
  opflow->formops.setnumvariables = OPFLOWFormulationSetNumVariables_PBCAR;
  opflow->formops.setnumconstraints = OPFLOWFormulationSetNumConstraints_PBCAR;
  opflow->formops.setvariablebounds = OPFLOWSetVariableBounds_PBCAR;
  opflow->formops.setconstraintbounds = OPFLOWSetConstraintBounds_PBCAR;
  opflow->formops.setvariableandconstraintbounds = OPFLOWSetVariableandConstraintBounds_PBCAR;
  opflow->formops.setinitialguess = OPFLOWSetInitialGuess_PBCAR;
  opflow->formops.computeequalityconstraints = OPFLOWComputeEqualityConstraints_PBCAR;
  opflow->formops.computeinequalityconstraints = OPFLOWComputeInequalityConstraints_PBCAR;
  opflow->formops.computeequalityconstraintjacobian = OPFLOWComputeEqualityConstraintJacobian_PBCAR;
  opflow->formops.computeinequalityconstraintjacobian = OPFLOWComputeInequalityConstraintJacobian_PBCAR;
  opflow->formops.computehessian = OPFLOWComputeHessian_PBCAR;
  opflow->formops.computeobjandgradient = OPFLOWComputeObjandGradient_PBCAR;
  opflow->formops.computeobjective = OPFLOWComputeObjective_PBCAR;
  opflow->formops.computegradient  = OPFLOWComputeGradient_PBCAR;
  
  PetscFunctionReturn(0);
}
