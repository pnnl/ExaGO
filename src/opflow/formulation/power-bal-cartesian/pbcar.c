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

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(Xu,&xu);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    PetscInt k;

    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);

    /* Bounds on real and imaginary part of voltages */
    xl[loc] = PETSC_NINFINITY; xu[loc] = PETSC_INFINITY;
    xl[loc+1] = PETSC_NINFINITY; xu[loc+1] = PETSC_INFINITY;

    if(bus->ide == ISOLATED_BUS) {
      xl[loc] = xu[loc] = bus->vm*PetscCosScalar(bus->va*PETSC_PI/180.0);
      xl[loc+1] = xu[loc+1] = bus->vm*PetscSinScalar(bus->va*PETSC_PI/180.0);
    }

    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;

      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      loc = loc+2;
      if(!gen->status) xl[loc] = xu[loc] = xl[loc+1] = xu[loc+1] = 0.0;
      else {
	xl[loc] = gen->pb; /* PGmin */
	xu[loc] = gen->pt; /* PGmax */
	xl[loc+1] = gen->qb; /* QGmin */
	xu[loc+1] = gen->qt; /* QGmax */
	/* pb, pt, qb, qt are converted in p.u. in ps.c */
      }
    }
  }
  
  ierr = VecRestoreArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Xu,&xu);CHKERRQ(ierr);

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

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(Gu,&gu);CHKERRQ(ierr);

  /* Equality constraint bounds for balance equations and ref. bus angle */
  for(i=0; i < ps->nbus; i++) {

    bus = &ps->bus[i];

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
    
    gl[gloc] = bus->Vmin*bus->Vmin;
    gu[gloc] = bus->Vmax*bus->Vmax;
    if(bus->ide == ISOLATED_BUS) {
      gl[gloc] = PETSC_NINFINITY;
      gu[gloc] = PETSC_INFINITY;
    }
    gloc++;
  }
  

    for(i=0; i < ps->Nbranch; i++) {
      line = &ps->line[i];
      
      /* Line flow inequality constraints */
      if(!line->status) gl[gloc] = gu[gloc] = gl[gloc+1] = gu[gloc+1] = 0.0;
      else {
	gl[gloc] = gl[gloc+1] = 0.0; 
	gu[gloc] = gu[gloc+1] = (line->rateA/ps->MVAbase)*(line->rateA/ps->MVAbase);
      }    
      gloc += 2;
    }
  }

  ierr = VecRestoreArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gu,&gu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetVariableandConstraintBounds_PBCAR(OPFLOW opflow,Vec Xl,Vec Xu, Vec Gl, Vec Gu)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetInitialGuess_PBCAR(OPFLOW opflow,Vec X)
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

    /* Initial guess for real and imaginary part of voltages */
    x[loc]   = 1.0;
    x[loc+1] = 0.0;
    
    if(bus->ide == ISOLATED_BUS) {
      x[loc]   = bus->vm*PetscCosScalar(bus->va*PETSC_PI/180.0);
      x[loc+1] = bus->vm*PetscCosScalar(bus->va*PETSC_PI/180.0);
    }

    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;

      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      loc = loc+2;

      x[loc]   = (xl[loc] + xu[loc])/2.0;   /* Initial guess for Pg */
      x[loc+1] = (xl[loc+1] + xu[loc+1])/2.0; /* Initial guess for Qg */
    }
  }

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(opflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(opflow->Xu,&xu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeEqualityConstraints_PBCAR(OPFLOW opflow,Vec X,Vec Ge)
{
  PetscErrorCode ierr;
  PetscInt       i,k,nconnlines;
  PetscInt       gloc=0,row[3];
  PetscInt       xloc,xlocf,xloct;
  PetscScalar    val[3];
  PetscScalar    Pg,Qg;
  PetscScalar    Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  PetscScalar    Vrf,Vif,Vrt,Vit;
  PetscScalar    Pf,Qf,Pt,Qt;
  PetscScalar    Vr,Vi,Vm;
  PS             ps=opflow->ps;
  PSLOAD         load;
  PSLINE         line;
  PSBUS          bus,busf,bust;
  const PSBUS    *connbuses;
  const PSLINE   *connlines;
  const PetscScalar *x;
  PetscInt       rstart;

  PetscFunctionBegin;
  ierr = VecSet(Ge,0.0);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(Ge,&rstart,NULL);CHKERRQ(ierr);

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);

  for (i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    row[0] = rstart + gloc; row[1] = row[0] + 1;

    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);

    /* Real and imaginary voltages for the bus */
    Vr = x[xloc];
    Vi = x[xloc+1];
    Vm = PetscSqrtScalar(Vr*Vr + Vi*Vi);
    
    if (bus->ide == ISOLATED_BUS) {
      val[0] = Vr - bus->vm*PetscCosScalar(bus->va*PETSC_PI/180.0);
      val[1] = Vi - bus->vm*PetscCosScalar(bus->va*PETSC_PI/180.0);
  
      ierr = VecSetValues(Ge,2,row,val,ADD_VALUES);CHKERRQ(ierr);
      gloc += 2;
      continue;
    }

    /* Shunt injections */
    val[0] = Vm*Vm*bus->gl;
    val[1] = -Vm*Vm*bus->bl;
    ierr = VecSetValues(Ge,2,row,val,ADD_VALUES);CHKERRQ(ierr);

    for (k=0; k < bus->ngen; k++) {
      xloc += 2;
      Pg = x[xloc];
      Qg = x[xloc+1];
      
      val[0] = -Pg;
      val[1] = -Qg;
      ierr = VecSetValues(Ge,2,row,val,ADD_VALUES);CHKERRQ(ierr);      
    }

    for (k=0; k < bus->nload; k++) {
      ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
      val[0] = load->pl;
      val[1] = load->ql;
      ierr = VecSetValues(Ge,2,row,val,ADD_VALUES);CHKERRQ(ierr);
    }

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

        val[0] = Pf;
        val[1] = Qf;
        ierr = VecSetValues(Ge,2,row,val,ADD_VALUES);CHKERRQ(ierr);
      } else {
	Pt =  Gtt*(Vrt*Vrt + Vit*Vit) + Vrt*(Gtf*Vrf - Btf*Vif) + Vit*(Btf*Vrf + Gtf*Vif);
	Qt = -Btt*(Vrt*Vrt + Vit*Vit) + Vit*(Gtf*Vrf - Btf*Vif) - Vrt*(Btf*Vrf + Gtf*Vif);

        val[0] = Pt;
        val[1] = Qt;
        ierr = VecSetValues(Ge,2,row,val,ADD_VALUES);CHKERRQ(ierr);
      }
    }
    gloc += 2;
    if(bus->ide == REF_BUS) {
      /* Equality constraint for angle */
      row[2] = row[1] + 1;
      val[2] = Vi - Vr*PetscTanScalar(bus->va*PETSC_PI/180.0); 
      ierr = VecSetValues(Ge,1,&row[2],&val[2],ADD_VALUES);CHKERRQ(ierr);
      gloc += 1;
    }
  }
  ierr = VecAssemblyBegin(Ge);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(Ge);CHKERRQ(ierr);

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeEqualityConstraintJacobian_PBCAR(OPFLOW opflow,Vec X,Mat Je)
{
  PetscErrorCode ierr;
  PetscInt       i,k,row[3],col[4],genctr,gloc=0;
  PetscInt       nconnlines,locglob,loc,locglobf,locglobt,locf,loct;
  PetscScalar    Vr,Vi,Vm,val[8],Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  PetscScalar    Vrf,Vrt,Vif,Vit;
  PS             ps=opflow->ps;
  PSBUS          bus;
  PSLINE         line;
  PSBUS          busf,bust;
  const PSLINE   *connlines;
  const PSBUS    *connbuses;
  const PetscScalar *xarr;
  PetscInt       rstart;

  PetscFunctionBegin;
  ierr = MatZeroEntries(Je);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(Je,&rstart,NULL);CHKERRQ(ierr);

  ierr = VecGetArrayRead(X,&xarr);CHKERRQ(ierr);
  for (i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    row[0] = rstart + gloc; row[1] = row[0] + 1;

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bus,&locglob);CHKERRQ(ierr);

    Vr = xarr[loc];
    Vi = xarr[loc+1];
    Vm = PetscSqrtScalar(Vr*Vr + Vi*Vi);

    col[0] = locglob; col[1] = locglob+1;
    /* Isolated and reference bus */
    if(bus->ide == ISOLATED_BUS) {
      val[0] = val[3] = 1.0;
      val[1] = val[2] = 0.0;
      ierr = MatSetValues(Je,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
      gloc += 2;
      continue;
    }
    /* Partial derivative for shunt contribution */
    val[0] = 2*Vr*bus->gl;  val[1] = 2*Vi*bus->gl;
    val[2] = -2*Vr*bus->bl; val[3]= -2*Vi*bus->bl; 
    ierr = MatSetValues(Je,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
    
    genctr = 0;
    for (k=0; k < bus->ngen; k++) {
      val[0] = -1;
      col[0] = locglob + 2 + genctr;
      ierr = MatSetValues(Je,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
      
      col[0] = locglob + 2 + genctr+1;
      ierr = MatSetValues(Je,1,row+1,1,col,val,ADD_VALUES);CHKERRQ(ierr);
      genctr += 2;
      
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
	val[0] = 2*Gff*Vrf + (Gft*Vrt - Bft*Vit);
	val[1] = 2*Gff*Vif + (Bft*Vrt + Gft*Vit);
	val[2] = Vrf*Gft + Vif*Bft;
	val[3] = -Vrf*Bft + Vif*Gft;

	val[4] = -2*Bff*Vrf - (Bft*Vrt + Gft*Vit);
	val[5] = -2*Bff*Vif + (Gft*Vrt - Bft*Vit);
	val[6] =  Vif*Gft - Vrf*Bft;
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
    gloc += 2;
    if(bus->ide == REF_BUS) {
      /* Partial of equality constraint for ref. bus angle */
      row[2] = row[1] + 1;
      col[0] = locglob; col[1] = locglob + 1;
      val[0] = -PetscTanScalar(bus->va*PETSC_PI/180.0);
      val[1] = 1.0;
      ierr = MatSetValues(Je,1,row+2,2,col,val,ADD_VALUES);CHKERRQ(ierr);
      gloc += 1;
    }
  }
  ierr = VecRestoreArrayRead(X,&xarr);CHKERRQ(ierr);

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

  PetscFunctionBegin;
  ierr = VecSet(Gi,0.0);CHKERRQ(ierr);

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Gi,&g);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);

    Vr = x[xloc];
    Vi = x[xloc+1];
    
    g[gloc] = Vr*Vr + Vi*Vi;
    gloc++;
  }

  for(i=0; i<ps->nbranch; i++) {
    line = &ps->line[i];

    if(!line->status) {
      gloc += 2;
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
  
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gi,&g);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeInequalityConstraintJacobian_PBCAR(OPFLOW opflow,Vec X,Mat Ji)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscInt       row[2],col[4];
  PetscInt       rstart,rend;
  PetscInt       gloc=0,xloc,xlocf,xloct,xlocglob;
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

  PetscFunctionBegin;
  ierr = MatZeroEntries(Ji);CHKERRQ(ierr);

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(Ji,&rstart,&rend);CHKERRQ(ierr);
  ierr = MPI_Barrier(comm);CHKERRQ(ierr);

  gloc = rstart;
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

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
    if(!line->status) {
      gloc += 2;
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

    Vrf  = x[xlocf];
    Vif     = x[xlocf+1];
    Vrt  = x[xloct];
    Vit     = x[xloct+1];

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
    col[0] = xlocf; col[1] = xlocf+1; col[2] = xloct; col[3] = xloct+1;
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
    col[0] = xloct; col[1] = xloct+1; col[2] = xlocf; col[3] = xlocf+1;
    val[0]   = dSt2_dVrt;
    val[1] = dSt2_dVit;
    val[2] = dSt2_dVrf;
    val[3] = dSt2_dVif;
    ierr = MatSetValues(Ji,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

    gloc += 2;
  }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(Ji,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Ji,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeConstraints_PBCAR(OPFLOW opflow,Vec X,Vec G)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeObjandGradient_PBCAR(OPFLOW opflow,Vec X,PetscScalar *obj,Vec Grad)
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
  PetscInt       loc;

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);

  *obj = 0.0;
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    
    PetscInt k;
    PSGEN    gen;
    PetscScalar Pg;
    for(k=0; k < bus->ngen; k++) {
      loc = loc+2;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      Pg = x[loc]*ps->MVAbase;
      *obj += gen->cost_alpha*Pg*Pg + gen->cost_beta*Pg + gen->cost_gamma;
    }
  }
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);

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
  PetscInt       loc;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(grad,&df);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    
    PetscInt k;
    PSGEN    gen;
    PetscScalar Pg;
    for(k=0; k < bus->ngen; k++) {
      loc = loc+2;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      Pg = x[loc]*ps->MVAbase;
      df[loc] = ps->MVAbase*(2*gen->cost_alpha*Pg + gen->cost_beta);
    }
  }
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(grad,&df);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWFormulationSetNumVariables_PBCAR(OPFLOW opflow,PetscInt *busnvar,PetscInt *branchnvar,PetscInt *nx)
{
  PetscInt i,ngen;
  PS       ps=opflow->ps;
  PSBUS    bus;
  PetscErrorCode ierr;

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
    ierr = PSBUSGetNGen(bus,&ngen);CHKERRQ(ierr);
    /* Number of variables = 2 + 2*ngen (2 variables for voltage + Pg, Qg for each gen) */
    busnvar[i] = 2 + 2*ngen;
    *nx += busnvar[i];
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWFormulationSetNumConstraints_PBCAR(OPFLOW opflow,PetscInt *nconeq,PetscInt *nconineq)
{
  PS  ps = opflow->ps;

  PetscFunctionBegin;
  *nconeq = 2*ps->nbus + ps->nref; /* Power balance constraint at each bus (2) + angle constraint at ref. bus */
  *nconineq = 2*ps->nbranch + ps->nbus; /* Line flow constraints + voltage magnitude constraints */

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
  PetscScalar    Vr,Vi;

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(Lambda,&lambda);CHKERRQ(ierr);

  // For equality constraints (power flow) */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    
    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bus,&xlocglob);CHKERRQ(ierr);
    
    Vr = x[xloc];
    Vi = x[xloc+1];
    
    /* Shunt */
    row[0] = xlocglob; row[1] = xlocglob+1;
    col[0] = xlocglob; col[1] = xlocglob+1;
    val[0] = lambda[gloc]*2*bus->gl + lambda[gloc+1]*(-2*bus->bl);
    val[3] = val[0];
    val[1] = val[2] = 0.0;
    ierr = MatSetValues(H,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
    
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
     
      PetscScalar Vrf,Vif,Vrt,Vit;
      
      Vrf  = x[xlocf];
      Vif  = x[xlocf+1];
      Vrt  = x[xloct];
      Vit  = x[xloct+1];
    
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

	val[0]  = lambda[gloc]*dPf_dVrf_dVrf + lambda[gloc+1]*dQf_dVrf_dVrf;
	val[1]  = lambda[gloc]*dPf_dVif_dVrf + lambda[gloc+1]*dQf_dVif_dVrf;
	val[2]  = lambda[gloc]*dPf_dVrt_dVrf + lambda[gloc+1]*dQf_dVrt_dVrf;
	val[3]  = lambda[gloc]*dPf_dVit_dVrf + lambda[gloc+1]*dQf_dVit_dVrf;

	ierr = MatSetValues(H,4,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	col[0] 	= xlocfglob+1;

	val[0]  = lambda[gloc]*dPf_dVrf_dVif + lambda[gloc+1]*dQf_dVrf_dVif;
	val[1]  = lambda[gloc]*dPf_dVif_dVif + lambda[gloc+1]*dQf_dVif_dVif;
	val[2]  = lambda[gloc]*dPf_dVrt_dVif + lambda[gloc+1]*dQf_dVrt_dVif;
	val[3]  = lambda[gloc]*dPf_dVit_dVif + lambda[gloc+1]*dQf_dVit_dVif;

	ierr = MatSetValues(H,4,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	col[0] 	= xloctglob;

	val[0]  = lambda[gloc]*dPf_dVrf_dVrt + lambda[gloc+1]*dQf_dVrf_dVrt;
	val[1]  = lambda[gloc]*dPf_dVif_dVrt + lambda[gloc+1]*dQf_dVif_dVrt;
	val[2]  = lambda[gloc]*dPf_dVrt_dVrt + lambda[gloc+1]*dQf_dVrt_dVrt;
	val[3]  = lambda[gloc]*dPf_dVit_dVrt + lambda[gloc+1]*dQf_dVit_dVrt;
	ierr = MatSetValues(H,4,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	col[0] 	= xloctglob+1;

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

	val[0]  = lambda[gloc]*dPt_dVrt_dVrt + lambda[gloc+1]*dQt_dVrt_dVrt;
	val[1]  = lambda[gloc]*dPt_dVit_dVrt + lambda[gloc+1]*dQt_dVit_dVrt;
	val[2]  = lambda[gloc]*dPt_dVrf_dVrt + lambda[gloc+1]*dQt_dVrf_dVrt;
	val[3]  = lambda[gloc]*dPt_dVif_dVrt + lambda[gloc+1]*dQt_dVif_dVrt;
	ierr = MatSetValues(H,4,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	col[0] 	= xloctglob+1;

	val[0]  = lambda[gloc]*dPt_dVrt_dVit + lambda[gloc+1]*dQt_dVrt_dVit;
	val[1]  = lambda[gloc]*dPt_dVit_dVit + lambda[gloc+1]*dQt_dVit_dVit;
	val[2]  = lambda[gloc]*dPt_dVrf_dVit + lambda[gloc+1]*dQt_dVrf_dVit;
	val[3]  = lambda[gloc]*dPt_dVif_dVit + lambda[gloc+1]*dQt_dVif_dVit;
	ierr = MatSetValues(H,4,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	col[0] 	= xlocfglob;

	val[0]  = lambda[gloc]*dPt_dVrt_dVrf + lambda[gloc+1]*dQt_dVrt_dVrf;
	val[1]  = lambda[gloc]*dPt_dVit_dVrf + lambda[gloc+1]*dQt_dVit_dVrf;
	val[2]  = lambda[gloc]*dPt_dVrf_dVrf + lambda[gloc+1]*dQt_dVrf_dVrf;
	val[3]  = lambda[gloc]*dPt_dVif_dVrf + lambda[gloc+1]*dQt_dVif_dVrf;
	ierr = MatSetValues(H,4,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	col[0] 	= xlocfglob+1;

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

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Lambda,&lambda);CHKERRQ(ierr);
  
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
  PetscInt       xlocf,xloct;
  PSBUS          busf,bust;
  const PetscScalar *x;
  const PetscScalar *lambda;
  PetscInt       gloc;
  PetscInt       row[12],col[12];
  PetscScalar    val[12];

  PetscFunctionBegin;

  gloc = opflow->nconeq; /* offset for the inequality constraints in the Lambda vector */
      
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(Lambda,&lambda);CHKERRQ(ierr);

  // for the part of line constraints
  for(i=0; i < ps->Nbranch; i++) {
    line = &ps->line[i];

    if(!line->status) {
      gloc += 2;
      continue;
    }

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

    PetscScalar Vmf,Vmt,thetaf,thetat,thetaft,thetatf;
      
    thetaf  = x[xlocf];
    Vmf     = x[xlocf+1];
    thetat  = x[xloct];
    Vmt     = x[xloct+1];
    thetaft = thetaf - thetat;
    thetatf = thetat - thetaf;

    // Sf2 and St2 are the constraints	
    PetscScalar Pf,Qf,Pt,Qt,Sf2,St2;

    Pf =  Gff*Vmf*Vmf + Vmf*Vmt*( Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
    Qf = -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
	
    Pt =  Gtt*Vmt*Vmt + Vmt*Vmf*( Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
    Qt = -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));

    Sf2 = Pf*Pf + Qf*Qf;
    St2 = Pt*Pt + Qt*Qt;

    PetscScalar dSf2_dPf, dSf2_dQf, dSt2_dPt, dSt2_dQt;

    dSf2_dPf = 2.*Pf;
    dSf2_dQf = 2.*Qf;
    dSt2_dPt = 2.*Pt;
    dSt2_dQt = 2.*Qt;

    PetscScalar dSf2_dPf_dPf, dSt2_dPt_dPt;
    PetscScalar dSf2_dQf_dQf, dSt2_dQt_dQt;

    dSf2_dPf_dPf = 2.;
    dSf2_dQf_dQf = 2.;	
    dSt2_dPt_dPt = 2.;
    dSt2_dQt_dQt = 2.;
	

    PetscScalar dPf_dthetaf,dPf_dVmf,dPf_dthetat,dPf_dVmt;
    PetscScalar dQf_dthetaf,dQf_dVmf,dQf_dthetat,dQf_dVmt;
    PetscScalar dPt_dthetaf,dPt_dVmf,dPt_dthetat,dPt_dVmt;
    PetscScalar dQt_dthetaf,dQt_dVmf,dQt_dthetat,dQt_dVmt;

    dPf_dthetaf = 			Vmf*Vmt*(-Gft*PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft));
    dPf_dVmf    = 2.*Gff*Vmf + 	Vmt*( Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
    dPf_dthetat = 			Vmf*Vmt*( Gft*PetscSinScalar(thetaft) - Bft*PetscCosScalar(thetaft));
    dPf_dVmt    = 				Vmf*( Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));

    dQf_dthetaf = 			Vmf*Vmt*( Bft*PetscSinScalar(thetaft) + Gft*PetscCosScalar(thetaft));
    dQf_dVmf    = -2.*Bff*Vmf + 	Vmt*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
    dQf_dthetat = 			Vmf*Vmt*(-Bft*PetscSinScalar(thetaft) - Gft*PetscCosScalar(thetaft));
    dQf_dVmt    = 				Vmf*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));

    dPt_dthetat = 			Vmt*Vmf*(-Gtf*PetscSinScalar(thetatf) + Btf*PetscCosScalar(thetatf));
    dPt_dVmt    = 2.*Gtt*Vmt + 	Vmf*( Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
    dPt_dthetaf = 			Vmt*Vmf*( Gtf*PetscSinScalar(thetatf) - Btf*PetscCosScalar(thetatf));
    dPt_dVmf    = 				Vmt*( Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));

    dQt_dthetat = 			Vmt*Vmf*( Btf*PetscSinScalar(thetatf) + Gtf*PetscCosScalar(thetatf));
    dQt_dVmt    = -2.*Btt*Vmt + 	Vmf*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
    dQt_dthetaf = 			Vmt*Vmf*(-Btf*PetscSinScalar(thetatf) - Gtf*PetscCosScalar(thetatf));
    dQt_dVmf    = 				Vmt*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));

    PetscScalar dPf_dthetaf_dthetaf,dPf_dVmf_dthetaf,dPf_dthetat_dthetaf,dPf_dVmt_dthetaf;
    PetscScalar dPf_dVmf_dVmf,dPf_dthetat_dVmf,dPf_dVmt_dVmf;
    PetscScalar dPf_dthetat_dthetat,dPf_dVmt_dthetat;
    PetscScalar dPf_dVmt_dVmt;

    dPf_dthetaf_dthetaf = 			Vmf*Vmt*(-Gft*PetscCosScalar(thetaft) - Bft*PetscSinScalar(thetaft));
    dPf_dVmf_dthetaf    =				Vmt*(-Gft*PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft));
    dPf_dthetat_dthetaf = 			Vmf*Vmt*( Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
	dPf_dVmt_dthetaf    = 				Vmf*(-Gft*PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft));

    dPf_dVmf_dVmf    	= 2.*Gff;
    dPf_dthetat_dVmf 	= 				Vmt*( Gft*PetscSinScalar(thetaft) - Bft*PetscCosScalar(thetaft));
    dPf_dVmt_dVmf    	= 					( Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));

    dPf_dthetat_dthetat = 			Vmf*Vmt*(-Gft*PetscCosScalar(thetaft) - Bft*PetscSinScalar(thetaft));
    dPf_dVmt_dthetat    = 				Vmf*( Gft*PetscSinScalar(thetaft) - Bft*PetscCosScalar(thetaft));

    dPf_dVmt_dVmt    	= 				0.;

    PetscScalar dQf_dthetaf_dthetaf,dQf_dVmf_dthetaf,dQf_dthetat_dthetaf,dQf_dVmt_dthetaf;
    PetscScalar dQf_dVmf_dVmf,dQf_dthetat_dVmf,dQf_dVmt_dVmf;
    PetscScalar dQf_dthetat_dthetat,dQf_dVmt_dthetat;
    PetscScalar dQf_dVmt_dVmt;

    dQf_dthetaf_dthetaf = 			Vmf*Vmt*( Bft*PetscCosScalar(thetaft) - Gft*PetscSinScalar(thetaft));
    dQf_dVmf_dthetaf    = 				Vmt*( Bft*PetscSinScalar(thetaft) + Gft*PetscCosScalar(thetaft));
    dQf_dthetat_dthetaf = 			Vmf*Vmt*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
    dQf_dVmt_dthetaf    = 				Vmf*( Bft*PetscSinScalar(thetaft) + Gft*PetscCosScalar(thetaft));

    dQf_dVmf_dVmf    	= -2.*Bff;
    dQf_dthetat_dVmf 	= 				Vmt*(-Bft*PetscSinScalar(thetaft) - Gft*PetscCosScalar(thetaft));
    dQf_dVmt_dVmf    	= 					(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));

    dQf_dthetat_dthetat = 			Vmf*Vmt*( Bft*PetscCosScalar(thetaft) - Gft*PetscSinScalar(thetaft));
    dQf_dVmt_dthetat    = 				Vmf*(-Bft*PetscSinScalar(thetaft) - Gft*PetscCosScalar(thetaft));

    dQf_dVmt_dVmt    	= 				0.;
	
    PetscScalar dPt_dthetat_dthetat,dPt_dVmt_dthetat,dPt_dthetaf_dthetat,dPt_dVmf_dthetat;
    PetscScalar dPt_dVmt_dVmt,dPt_dthetaf_dVmt,dPt_dVmf_dVmt;
    PetscScalar dPt_dthetaf_dthetaf,dPt_dVmf_dthetaf;
    PetscScalar dPt_dVmf_dVmf;
	
    dPt_dthetat_dthetat = 			Vmt*Vmf*(-Gtf*PetscCosScalar(thetatf) - Btf*PetscSinScalar(thetatf));
    dPt_dVmt_dthetat    = 				Vmf*(-Gtf*PetscSinScalar(thetatf) + Btf*PetscCosScalar(thetatf));
    dPt_dthetaf_dthetat = 			Vmt*Vmf*( Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
    dPt_dVmf_dthetat    = 				Vmt*(-Gtf*PetscSinScalar(thetatf) + Btf*PetscCosScalar(thetatf));

    dPt_dVmt_dVmt   	= 2.*Gtt;
    dPt_dthetaf_dVmt 	= 				Vmf*( Gtf*PetscSinScalar(thetatf) - Btf*PetscCosScalar(thetatf));
    dPt_dVmf_dVmt    	= 					( Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));

    dPt_dthetaf_dthetaf = 			Vmt*Vmf*(-Gtf*PetscCosScalar(thetatf) - Btf*PetscSinScalar(thetatf));
    dPt_dVmf_dthetaf    = 				Vmt*( Gtf*PetscSinScalar(thetatf) - Btf*PetscCosScalar(thetatf));

    dPt_dVmf_dVmf    	= 				0.;

    PetscScalar dQt_dthetat_dthetat,dQt_dVmt_dthetat,dQt_dthetaf_dthetat,dQt_dVmf_dthetat;
    PetscScalar dQt_dVmt_dVmt,dQt_dthetaf_dVmt,dQt_dVmf_dVmt;
    PetscScalar dQt_dthetaf_dthetaf,dQt_dVmf_dthetaf;
    PetscScalar dQt_dVmf_dVmf;	

    dQt_dthetat_dthetat = 			Vmt*Vmf*( Btf*PetscCosScalar(thetatf) - Gtf*PetscSinScalar(thetatf));
    dQt_dVmt_dthetat    = 				Vmf*( Btf*PetscSinScalar(thetatf) + Gtf*PetscCosScalar(thetatf));
    dQt_dthetaf_dthetat = 			Vmt*Vmf*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
    dQt_dVmf_dthetat    = 				Vmt*( Btf*PetscSinScalar(thetatf) + Gtf*PetscCosScalar(thetatf));

    dQt_dVmt_dVmt    	= -2.*Btt;
    dQt_dthetaf_dVmt 	= 				Vmf*(-Btf*PetscSinScalar(thetatf) - Gtf*PetscCosScalar(thetatf));
    dQt_dVmf_dVmt    	= 					(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));

    dQt_dthetaf_dthetaf = 			Vmt*Vmf*( Btf*PetscCosScalar(thetatf) - Gtf*PetscSinScalar(thetatf));
    dQt_dVmf_dthetaf    = 				Vmt*(-Btf*PetscSinScalar(thetatf) - Gtf*PetscCosScalar(thetatf));

    dQt_dVmf_dVmf    	= 				0.;


    PetscScalar dSf2_dPf_dthetaf, dSf2_dQf_dthetaf, dSt2_dPt_dthetaf, dSt2_dQt_dthetaf;
    dSf2_dPf_dthetaf = dSf2_dPf_dPf*dPf_dthetaf;
    dSf2_dQf_dthetaf = dSf2_dQf_dQf*dQf_dthetaf;
    dSt2_dPt_dthetaf = dSt2_dPt_dPt*dPt_dthetaf;
    dSt2_dQt_dthetaf = dSt2_dQt_dQt*dQt_dthetaf;

    PetscScalar dSf2_dPf_dVmf, dSf2_dQf_dVmf, dSt2_dPt_dVmf, dSt2_dQt_dVmf;
    dSf2_dPf_dVmf = dSf2_dPf_dPf*dPf_dVmf;
    dSf2_dQf_dVmf = dSf2_dQf_dQf*dQf_dVmf;
    dSt2_dPt_dVmf = dSt2_dPt_dPt*dPt_dVmf;
    dSt2_dQt_dVmf = dSt2_dQt_dQt*dQt_dVmf;

    PetscScalar dSf2_dPf_dthetat, dSf2_dQf_dthetat, dSt2_dPt_dthetat, dSt2_dQt_dthetat;
    dSf2_dPf_dthetat = dSf2_dPf_dPf*dPf_dthetat;
    dSf2_dQf_dthetat = dSf2_dQf_dQf*dQf_dthetat;
    dSt2_dPt_dthetat = dSt2_dPt_dPt*dPt_dthetat;
    dSt2_dQt_dthetat = dSt2_dQt_dQt*dQt_dthetat;

    PetscScalar dSf2_dPf_dVmt, dSf2_dQf_dVmt, dSt2_dPt_dVmt, dSt2_dQt_dVmt;
    dSf2_dPf_dVmt = dSf2_dPf_dPf*dPf_dVmt;
    dSf2_dQf_dVmt = dSf2_dQf_dQf*dQf_dVmt;
    dSt2_dPt_dVmt = dSt2_dPt_dPt*dPt_dVmt;
    dSt2_dQt_dVmt = dSt2_dQt_dQt*dQt_dVmt;


    PetscScalar dSf2_dthetaf_dthetaf,dSf2_dVmf_dthetaf,dSf2_dthetat_dthetaf,dSf2_dVmt_dthetaf;
    dSf2_dthetaf_dthetaf =  	dSf2_dPf_dthetaf*dPf_dthetaf 	+ dSf2_dPf*dPf_dthetaf_dthetaf
						  	+ 	dSf2_dQf_dthetaf*dQf_dthetaf 	+ dSf2_dQf*dQf_dthetaf_dthetaf;
    dSf2_dVmf_dthetaf    =  	dSf2_dPf_dthetaf*dPf_dVmf 		+ dSf2_dPf*dPf_dVmf_dthetaf    
						  	+ 	dSf2_dQf_dthetaf*dQf_dVmf 		+ dSf2_dQf*dQf_dVmf_dthetaf;	
    dSf2_dthetat_dthetaf =  	dSf2_dPf_dthetaf*dPf_dthetat 	+ dSf2_dPf*dPf_dthetat_dthetaf 
						  	+ 	dSf2_dQf_dthetaf*dQf_dthetat 	+ dSf2_dQf*dQf_dthetat_dthetaf;
    dSf2_dVmt_dthetaf    = 		dSf2_dPf_dthetaf*dPf_dVmt 		+ dSf2_dPf*dPf_dVmt_dthetaf    
							+ 	dSf2_dQf_dthetaf*dQf_dVmt    	+ dSf2_dQf*dQf_dVmt_dthetaf;

    PetscScalar dSf2_dVmf_dVmf,dSf2_dthetat_dVmf,dSf2_dVmt_dVmf;
    dSf2_dVmf_dVmf    	= 		dSf2_dPf_dVmf*dPf_dVmf    + dSf2_dPf*dPf_dVmf_dVmf    
							+   dSf2_dQf_dVmf*dQf_dVmf	  + dSf2_dQf*dQf_dVmf_dVmf;
    dSf2_dthetat_dVmf 	= 		dSf2_dPf_dVmf*dPf_dthetat + dSf2_dPf*dPf_dthetat_dVmf
							+	dSf2_dQf_dVmf*dQf_dthetat + dSf2_dQf*dQf_dthetat_dVmf;	
    dSf2_dVmt_dVmf    	= 		dSf2_dPf_dVmf*dPf_dVmt    + dSf2_dPf*dPf_dVmt_dVmf
							+	dSf2_dQf_dVmf*dQf_dVmt    + dSf2_dQf*dQf_dVmt_dVmf;

    PetscScalar dSf2_dthetat_dthetat,dSf2_dVmt_dthetat;
    dSf2_dthetat_dthetat = 		dSf2_dPf_dthetat*dPf_dthetat + dSf2_dPf*dPf_dthetat_dthetat
							+	dSf2_dQf_dthetat*dQf_dthetat + dSf2_dQf*dQf_dthetat_dthetat;	
    dSf2_dVmt_dthetat    = 		dSf2_dPf_dthetat*dPf_dVmt    + dSf2_dPf*dPf_dVmt_dthetat
							+	dSf2_dQf_dthetat*dQf_dVmt    + dSf2_dQf*dQf_dVmt_dthetat;
	
    PetscScalar dSf2_dVmt_dVmt;
    dSf2_dVmt_dVmt   	 = 		dSf2_dPf_dVmt*dPf_dVmt    + dSf2_dPf*dPf_dVmt_dVmt
							+	dSf2_dQf_dVmt*dQf_dVmt    + dSf2_dQf*dQf_dVmt_dVmt;

    PetscScalar dSt2_dthetaf_dthetaf,dSt2_dVmf_dthetaf,dSt2_dthetat_dthetaf,dSt2_dVmt_dthetaf;
    dSt2_dthetaf_dthetaf =  	dSt2_dPt_dthetaf*dPt_dthetaf 	+ dSt2_dPt*dPt_dthetaf_dthetaf
						  	+ 	dSt2_dQt_dthetaf*dQt_dthetaf 	+ dSt2_dQt*dQt_dthetaf_dthetaf;
    dSt2_dVmf_dthetaf    =  	dSt2_dPt_dthetaf*dPt_dVmf 		+ dSt2_dPt*dPt_dVmf_dthetaf    
						  	+ 	dSt2_dQt_dthetaf*dQt_dVmf 		+ dSt2_dQt*dQt_dVmf_dthetaf;	
    dSt2_dthetat_dthetaf =  	dSt2_dPt_dthetaf*dPt_dthetat 	+ dSt2_dPt*dPt_dthetaf_dthetat 
						  	+ 	dSt2_dQt_dthetaf*dQt_dthetat 	+ dSt2_dQt*dQt_dthetaf_dthetat;
    dSt2_dVmt_dthetaf    = 		dSt2_dPt_dthetaf*dPt_dVmt 		+ dSt2_dPt*dPt_dthetaf_dVmt    
							+ 	dSt2_dQt_dthetaf*dQt_dVmt    	+ dSt2_dQt*dQt_dthetaf_dVmt;

    PetscScalar dSt2_dVmf_dVmf,dSt2_dthetat_dVmf,dSt2_dVmt_dVmf;
    dSt2_dVmf_dVmf    	= 		dSt2_dPt_dVmf*dPt_dVmf    + dSt2_dPt*dPt_dVmf_dVmf 
							+	dSt2_dQt_dVmf*dQt_dVmf    + dSt2_dQt*dQt_dVmf_dVmf;
    dSt2_dthetat_dVmf 	= 		dSt2_dPt_dVmf*dPt_dthetat + dSt2_dPt*dPt_dVmf_dthetat
							+	dSt2_dQt_dVmf*dQt_dthetat + dSt2_dQt*dQt_dVmf_dthetat;	
    dSt2_dVmt_dVmf    	= 		dSt2_dPt_dVmf*dPt_dVmt    + dSt2_dPt*dPt_dVmf_dVmt
							+	dSt2_dQt_dVmf*dQt_dVmt    + dSt2_dQt*dQt_dVmf_dVmt;

    PetscScalar dSt2_dthetat_dthetat,dSt2_dVmt_dthetat;
    dSt2_dthetat_dthetat = 		dSt2_dPt_dthetat*dPt_dthetat + dSt2_dPt*dPt_dthetat_dthetat
							+	dSt2_dQt_dthetat*dQt_dthetat + dSt2_dQt*dQt_dthetat_dthetat;	
    dSt2_dVmt_dthetat    = 		dSt2_dPt_dthetat*dPt_dVmt    + dSt2_dPt*dPt_dVmt_dthetat
							+	dSt2_dQt_dthetat*dQt_dVmt    + dSt2_dQt*dQt_dVmt_dthetat;
	
    PetscScalar dSt2_dVmt_dVmt;
    dSt2_dVmt_dVmt   	 = 		dSt2_dPt_dVmt*dPt_dVmt    + dSt2_dPt*dPt_dVmt_dVmt
							+	dSt2_dQt_dVmt*dQt_dVmt    + dSt2_dQt*dQt_dVmt_dVmt;

    row[0] 	= xlocf; 
    row[1] 	= xlocf+1; 
    col[0] 	= xlocf;
    val[0]   = lambda[gloc]*dSf2_dthetaf_dthetaf 	+ lambda[gloc+1]*dSt2_dthetaf_dthetaf;
    val[1] = lambda[gloc]*dSf2_dVmf_dthetaf 		+ lambda[gloc+1]*dSt2_dVmf_dthetaf;
    ierr = MatSetValues(H,2,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

    row[0] 	= xlocf+1; 	col[0] 	= xlocf+1;
    val[0] = lambda[gloc]*dSf2_dVmf_dVmf 		+ lambda[gloc+1]*dSt2_dVmf_dVmf;
    ierr = MatSetValues(H,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

    row[0] 	= xloct; 	
    col[0] 	= xlocf;
    col[1] 	= xlocf+1;
    col[2] 	= xloct;
    val[0] = lambda[gloc]*dSf2_dthetat_dthetaf 	+ lambda[gloc+1]*dSt2_dthetat_dthetaf;
    val[1] = lambda[gloc]*dSf2_dthetat_dVmf		+ lambda[gloc+1]*dSt2_dthetat_dVmf;
    val[2] = lambda[gloc]*dSf2_dthetat_dthetat 	+ lambda[gloc+1]*dSt2_dthetat_dthetat;
    ierr = MatSetValues(H,1,row,3,col,val,ADD_VALUES);CHKERRQ(ierr);


    row[0] 	= xloct+1; 	
    col[0] 	= xlocf;
    col[1] 	= xlocf+1;
    col[2] 	= xloct;
    col[3] 	= xloct+1;
    val[0] = lambda[gloc]*dSf2_dVmt_dthetaf 		+ lambda[gloc+1]*dSt2_dVmt_dthetaf;
    val[1] = lambda[gloc]*dSf2_dVmt_dVmf 		+ lambda[gloc+1]*dSt2_dVmt_dVmf;
    val[2] = lambda[gloc]*dSf2_dVmt_dthetat 		+ lambda[gloc+1]*dSt2_dVmt_dthetat;
    val[3] = lambda[gloc]*dSf2_dVmt_dVmt 		+ lambda[gloc+1]*dSt2_dVmt_dVmt;
    ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
    
    gloc +=  2;
  }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Lambda,&lambda);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeObjectiveHessian - Computes the Hessian for the objective function part
  
  Input Parameters:
+ opflow - the OPFLOW object
- X        - solution vecto X

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
  PetscInt       xloc;
  const PetscScalar *x;
  PetscInt       row[2],col[2];
  PetscScalar    val[2];
  PetscScalar    obj_factor = opflow->obj_factor;


  PetscFunctionBegin;

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);

  // for the part of objective
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);
   
    for(k=0; k < bus->ngen; k++) {
      xloc = xloc+2;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      row[0] = xloc;
      col[0] = xloc;
      val[0] = obj_factor*2.0*gen->cost_alpha*ps->MVAbase*ps->MVAbase;
      ierr = MatSetValues(H,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
    }
  }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeHessian_PBCAR(OPFLOW opflow,Vec X,Vec Lambda,Mat H)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = MatZeroEntries(H);CHKERRQ(ierr);

  /* Objective function Hessian */
  ierr = OPFLOWComputeObjectiveHessian_PBCAR(opflow,X,H);CHKERRQ(ierr);

  /* Equality constraints Hessian */
  ierr = OPFLOWComputeEqualityConstraintsHessian_PBCAR(opflow,X,Lambda,H);CHKERRQ(ierr);
  
  /* Inequality constraints Hessian */
  if(opflow->nconineq) {
    ierr = OPFLOWComputeInequalityConstraintsHessian_PBCAR(opflow,X,Lambda,H);CHKERRQ(ierr);
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
