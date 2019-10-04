#include <private/opflowimpl.h>
#include "pbpol.h"

PetscErrorCode OPFLOWFormulationDestroy_PBPOL(OPFLOW opflow)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscFree(opflow->formulation);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetVariableBounds_PBPOL(OPFLOW opflow,Vec Xl,Vec Xu)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetConstraintBounds_PBPOL(OPFLOW opflow,Vec Gl,Vec Gu)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetVariableandConstraintBounds_PBPOL(OPFLOW opflow,Vec Xl,Vec Xu, Vec Gl, Vec Gu)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscScalar    *xl,*xu,*gl,*gu;
  PetscInt       i;
  PSLINE         line;
  PSBUS          bus;
  PetscInt       loc,gloc=0;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(Xu,&xu);CHKERRQ(ierr);
  ierr = VecGetArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(Gu,&gu);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    PetscInt k;

    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);

    /* Bounds on voltage angles and bounds on real power mismatch equality constraints */
    xl[loc] = PETSC_NINFINITY; xu[loc] = PETSC_INFINITY;
    gl[gloc] = 0.0;   gu[gloc] = 0.0;

    /* Bounds on voltage magnitudes and bounds on reactive power mismatch equality constraints */
    xl[loc+1] = bus->Vmin; xu[loc+1] = bus->Vmax;
    gl[gloc+1] = 0.0;       gu[gloc+1] = 0.0;

    if(bus->ide == REF_BUS || bus->ide == ISOLATED_BUS) xl[loc] = xu[loc] = bus->va*PETSC_PI/180.0;
    if(bus->ide == ISOLATED_BUS) xl[loc+1] = xu[loc+1] = bus->vm;
    
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
    gloc += 2;
  }
  
  if(opflow->nconineq) {
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

  ierr = VecRestoreArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Xu,&xu);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gu,&gu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetInitialGuess_PBPOL(OPFLOW opflow,Vec X)
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

    /* Initial guess for voltage angles and bounds on voltage magnitudes */
    x[loc]   = (xl[loc] + xu[loc])/2.0;
    x[loc+1] = (xl[loc+1] + xu[loc+1])/2.0;

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

PetscErrorCode OPFLOWComputeEqualityConstraints_PBPOL(OPFLOW opflow,Vec X,Vec Ge)
{
  PetscErrorCode ierr;
  PetscInt       i,k,nconnlines;
  PetscInt       gloc,row[2];
  PetscInt       xloc,xlocf,xloct;
  PetscScalar    val[2];
  PetscScalar    Pg,Qg;
  PetscScalar    Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  PetscScalar    Vmf,Vmt,thetaf,thetat,thetaft,thetatf;
  PetscScalar    Pf,Qf,Pt,Qt;
  PetscScalar    theta,Vm;
  PS             ps=opflow->ps;
  PSLOAD         load;
  PSLINE         line;
  PSBUS          bus,busf,bust;
  const PSBUS    *connbuses;
  const PSLINE   *connlines;
  const PetscScalar *x;

  PetscFunctionBegin;
  ierr = VecSet(Ge,0.0);CHKERRQ(ierr);

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);

  for (i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = DMNetworkGetVertexLocalToGlobalOrdering(ps->networkdm,i,&gloc);CHKERRQ(ierr);
    gloc *= 2;
    row[0] = gloc; row[1] = row[0] + 1;

    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);
    theta = x[xloc];
    Vm    = x[xloc+1];
    
    if (bus->ide == ISOLATED_BUS) {
      row[0]  = gloc; row[1] = row[0] + 1;
      val[0] = theta - bus->va*PETSC_PI/180.0;
      val[1] = Vm - bus->vm;
      ierr = VecSetValues(Ge,2,row,val,ADD_VALUES);CHKERRQ(ierr);
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

      thetaf  = x[xlocf];
      Vmf     = x[xlocf+1];
      thetat  = x[xloct];
      Vmt     = x[xloct+1];
      thetaft = thetaf - thetat;
      thetatf = thetat - thetaf;

      if (bus == busf) {
      	Pf = Gff*Vmf*Vmf  + Vmf*Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
      	Qf = -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));

        val[0] = Pf;
        val[1] = Qf;
        ierr = VecSetValues(Ge,2,row,val,ADD_VALUES);CHKERRQ(ierr);
      } else {
      	Pt = Gtt*Vmt*Vmt  + Vmt*Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf));
      	Qt = -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));

        val[0] = Pt;
        val[1] = Qt;
        ierr = VecSetValues(Ge,2,row,val,ADD_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = VecAssemblyBegin(Ge);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(Ge);CHKERRQ(ierr);

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);


  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeEqualityConstraintJacobian_PBPOL(OPFLOW opflow,Vec X,Mat Je)
{
  PetscErrorCode ierr;
  PetscInt       i,k,row[2],col[4],genctr,gidx;
  PetscInt       nconnlines,locglob,loc,locglobf,locglobt,locf,loct;
  PetscScalar    Vm,val[8],Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  PetscScalar    thetaf,thetat,Vmf,Vmt,thetaft,thetatf;
  PS             ps=opflow->ps;
  PSBUS          bus;
  PSLINE         line;
  PSBUS          busf,bust;
  DM             networkdm=ps->networkdm;
  const PSLINE   *connlines;
  const PSBUS    *connbuses;
  const PetscScalar *xarr;

  PetscFunctionBegin;
  ierr = MatZeroEntries(Je);CHKERRQ(ierr);

  ierr = VecGetArrayRead(X,&xarr);CHKERRQ(ierr);
  for (i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = DMNetworkGetVertexLocalToGlobalOrdering(networkdm,i,&gidx);CHKERRQ(ierr);
    gidx *= 2;
    row[0] = gidx; row[1] = row[0] + 1;

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bus,&locglob);CHKERRQ(ierr);

    Vm = xarr[loc+1];

    col[0] = locglob; col[1] = locglob+1;
    /* Isolated and reference bus */
    if(bus->ide == ISOLATED_BUS) {
      val[0] = val[3] = 1.0;
      val[1] = val[2] = 0.0;
      ierr = MatSetValues(Je,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
      continue;
    }
    /* Shunt injections */
    val[0] = 0.0; val[1] = 2*Vm*bus->gl;
    val[2] = 0.0; val[3]= -2*Vm*bus->bl; /* Partial derivative for shunt contribution */
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

      thetaf = xarr[locf];  Vmf = xarr[locf+1];
      thetat = xarr[loct];  Vmt = xarr[loct+1];
      thetaft = thetaf - thetat;
      thetatf = thetat - thetaf;

      if (bus == busf) {
      	col[0] = locglobf; col[1] = locglobf+1; col[2] = locglobt; col[3] = locglobt+1;
      	val[0] = Vmf*Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft));
      	val[1] = 2*Gff*Vmf + Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
      	val[2] = Vmf*Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft));
      	val[3] = Vmf*(Gft*cos(thetaft) + Bft*sin(thetaft));

        val[4] = Vmf*Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft));
        val[5] = -2*Bff*Vmf + Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
        val[6] = Vmf*Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft));
        val[7] = Vmf*(-Bft*cos(thetaft) + Gft*sin(thetaft));
        ierr = MatSetValues(Je,2,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
      } else {
      	col[0] = locglobt; col[1] = locglobt+1; col[2] = locglobf; col[3] = locglobf+1;
      	val[0] = Vmt*Vmf*(-Gtf*sin(thetatf) + Btf*cos(thetatf));
      	val[1] = 2*Gtt*Vmt + Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf));
      	val[2] = Vmt*Vmf*(Gtf*sin(thetatf) - Btf*cos(thetatf));
      	val[3] = Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf));

        val[4] = Vmt*Vmf*(Btf*sin(thetatf) + Gtf*cos(thetatf));
        val[5] = -2*Btt*Vmt + Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
        val[6] = Vmt*Vmf*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
        val[7] = Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
        ierr = MatSetValues(Je,2,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = VecRestoreArrayRead(X,&xarr);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(Je,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Je,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeInequalityConstraints_PBPOL(OPFLOW opflow,Vec X,Vec Gi)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscInt       gloc=0;
  PetscInt       xlocf,xloct;
  PetscScalar    *g;
  PetscScalar    Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  PetscScalar    Vmf,Vmt,thetaf,thetat,thetaft,thetatf;
  PetscScalar    Pf,Qf,Pt,Qt,Sf2,St2;
  PS             ps=opflow->ps;
  PSLINE         line;
  PSBUS          busf,bust;
  const PSBUS    *connbuses;
  const PetscScalar *x;

  PetscFunctionBegin;
  ierr = VecSet(Gi,0.0);CHKERRQ(ierr);

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Gi,&g);CHKERRQ(ierr);

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

PetscErrorCode OPFLOWComputeInequalityConstraintJacobian_PBPOL(OPFLOW opflow,Vec X,Mat Ji)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscInt       row[2],col[4];
  PetscInt       rstart,rend;
  PetscInt       gloc=0,xlocf,xloct;
  PetscScalar    val[4];
  PetscScalar    Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  PetscScalar    Vmf,Vmt,thetaf,thetat,thetaft,thetatf;
  PetscScalar    Pf,Qf,Pt,Qt;
  PetscScalar    dSf2_dPf, dSf2_dQf, dSt2_dPt, dSt2_dQt;
  PetscScalar    dPf_dthetaf,dPf_dVmf,dPf_dthetat,dPf_dVmt;
  PetscScalar    dQf_dthetaf,dQf_dVmf,dQf_dthetat,dQf_dVmt;
  PetscScalar    dPt_dthetaf,dPt_dVmf,dPt_dthetat,dPt_dVmt;
  PetscScalar    dQt_dthetaf,dQt_dVmf,dQt_dthetat,dQt_dVmt;
  PetscScalar    dSf2_dthetaf,dSf2_dVmf,dSf2_dthetat,dSf2_dVmt;
  PetscScalar    dSt2_dthetaf,dSt2_dVmf,dSt2_dthetat,dSt2_dVmt;
  PS             ps=opflow->ps;
  MPI_Comm       comm=opflow->comm->type;
  PSLINE         line;
  PSBUS          busf,bust;
  const PSBUS    *connbuses;
  const PetscScalar *x;

  PetscFunctionBegin;
  ierr = MatZeroEntries(Ji);CHKERRQ(ierr);

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(Ji,&rstart,&rend);CHKERRQ(ierr);
  ierr = MPI_Barrier(comm);CHKERRQ(ierr);

  gloc = rstart;
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

    dSf2_dPf = 2*Pf;
    dSf2_dQf = 2*Qf;
    dSt2_dPt = 2*Pt;
    dSt2_dQt = 2*Qt;

    dPf_dthetaf = Vmf*Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft));
    dPf_dVmf    = 2*Gff*Vmf + Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
    dPf_dthetat = Vmf*Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft));
    dPf_dVmt    = Vmf*(Gft*cos(thetaft) + Bft*sin(thetaft));

    dQf_dthetaf = Vmf*Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft));
    dQf_dVmf    = -2*Bff*Vmf + Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
    dQf_dthetat = Vmf*Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft));
    dQf_dVmt    = Vmf*(-Bft*cos(thetaft) + Gft*sin(thetaft));

    dPt_dthetat = Vmt*Vmf*(-Gtf*sin(thetatf) + Btf*cos(thetatf));
    dPt_dVmt    = 2*Gtt*Vmt + Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf));
    dPt_dthetaf = Vmt*Vmf*(Gtf*sin(thetatf) - Btf*cos(thetatf));
    dPt_dVmf    = Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf));

    dQt_dthetat = Vmt*Vmf*(Btf*sin(thetatf) + Gtf*cos(thetatf));
    dQt_dVmt    = -2*Btt*Vmt + Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
    dQt_dthetaf = Vmt*Vmf*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
    dQt_dVmf    = Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf));

    dSf2_dthetaf = dSf2_dPf*dPf_dthetaf + dSf2_dQf*dQf_dthetaf;
    dSf2_dthetat = dSf2_dPf*dPf_dthetat + dSf2_dQf*dQf_dthetat;
    dSf2_dVmf    = dSf2_dPf*dPf_dVmf    + dSf2_dQf*dQf_dVmf;
    dSf2_dVmt    = dSf2_dPf*dPf_dVmt    + dSf2_dQf*dQf_dVmt;

    row[0] = gloc;
    col[0] = xlocf; col[1] = xlocf+1; col[2] = xloct; col[3] = xloct+1;
    val[0] = dSf2_dthetaf;
    val[1] = dSf2_dVmf;
    val[2] = dSf2_dthetat;
    val[3] = dSf2_dVmt;
    ierr = MatSetValues(Ji,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

    dSt2_dthetaf = dSt2_dPt*dPt_dthetaf + dSt2_dQt*dQt_dthetaf;
    dSt2_dthetat = dSt2_dPt*dPt_dthetat + dSt2_dQt*dQt_dthetat;
    dSt2_dVmf    = dSt2_dPt*dPt_dVmf    + dSt2_dQt*dQt_dVmf;
    dSt2_dVmt    = dSt2_dPt*dPt_dVmt    + dSt2_dQt*dQt_dVmt;

    row[0] = gloc+1;
    col[0] = xloct; col[1] = xloct+1; col[2] = xlocf; col[3] = xlocf+1;
    val[0]   = dSt2_dthetat;
    val[1] = dSt2_dVmt;
    val[2] = dSt2_dthetaf;
    val[3] = dSt2_dVmf;
    ierr = MatSetValues(Ji,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

    gloc += 2;
  }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(Ji,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Ji,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeConstraints_PBPOL(OPFLOW opflow,Vec X,Vec G)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeObjandGradient_PBPOL(OPFLOW opflow,Vec X,PetscScalar *obj,Vec Grad)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeObjective_PBPOL(OPFLOW opflow,Vec X,PetscScalar *obj)
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

PetscErrorCode OPFLOWComputeGradient_PBPOL(OPFLOW opflow,Vec X,Vec grad)
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

PetscErrorCode OPFLOWFormulationSetNumVariables_PBPOL(OPFLOW opflow,PetscInt *busnvar,PetscInt *branchnvar,PetscInt *nx)
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
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetNGen(bus,&ngen);CHKERRQ(ierr);
    /* Number of variables = 2 + 2*ngen (2 variables for voltage + Pg, Qg for each gen) */
    busnvar[i] = 2 + 2*ngen;
    *nx += busnvar[i];
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWFormulationSetNumConstraints_PBPOL(OPFLOW opflow,PetscInt *nconeq,PetscInt *nconineq)
{
  PS  ps = opflow->ps;

  PetscFunctionBegin;
  *nconeq = 2*ps->nbus;
  *nconineq = 2*ps->nbranch;

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
PetscErrorCode OPFLOWComputeEqualityConstraintsHessian_PBPOL(OPFLOW opflow,Vec X,Vec Lambda,Mat H) 
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscInt       i,k;
  PSLINE         line;
  PSBUS          bus;
  PetscInt       nconnlines;
  const PSLINE   *connlines;
  const PSBUS    *connbuses;
  PetscInt       xloc,xlocf,xloct;
  PSBUS          busf,bust;
  const PetscScalar *x;
  const PetscScalar *lambda;
  PetscInt       gloc=0;
  PetscInt       row[12],col[12];
  PetscScalar    val[12];

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(Lambda,&lambda);CHKERRQ(ierr);

  // For equality constraints (power flow) */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    
    PetscScalar theta,Vm;
    
    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);
    
    theta = x[xloc];
    Vm    = x[xloc+1];
    
    row[0] = xloc + 1; col[0] = xloc + 1;
    val[0] = lambda[gloc]*2*bus->gl + lambda[gloc+1]*(-2*bus->bl);
    ierr = MatSetValues(H,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
    
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
     
      PetscScalar Vmf,Vmt,thetaf,thetat,thetaft,thetatf;
      
      thetaf  = x[xlocf];
      Vmf     = x[xlocf+1];
      thetat  = x[xloct];
      Vmt     = x[xloct+1];
      thetaft = thetaf - thetat;
      thetatf = thetat - thetaf;
    
      if(bus == busf) {
	 
	PetscScalar dPf_dthetaf_dthetaf,dPf_dVmf_dthetaf,dPf_dthetat_dthetaf,dPf_dVmt_dthetaf;
	PetscScalar dPf_dVmf_dVmf,dPf_dthetat_dVmf,dPf_dVmt_dVmf;
	PetscScalar dPf_dthetat_dthetat,dPf_dVmt_dthetat;
	PetscScalar dPf_dVmt_dVmt;
	
	dPf_dthetaf_dthetaf = Vmf*Vmt*(-Gft*PetscCosScalar(thetaft) - Bft*PetscSinScalar(thetaft));
	dPf_dVmf_dthetaf    = Vmt*(-Gft*PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft));
	dPf_dthetat_dthetaf = Vmf*Vmt*( Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
	dPf_dVmt_dthetaf    = Vmf*(-Gft*PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft));
	
	dPf_dVmf_dVmf	    = 2*Gff;
	dPf_dthetat_dVmf    = Vmt*( Gft*PetscSinScalar(thetaft) - Bft*PetscCosScalar(thetaft));
	dPf_dVmt_dVmf	    =	  ( Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
	
	dPf_dthetat_dthetat = Vmf*Vmt*(-Gft*PetscCosScalar(thetaft) - Bft*PetscSinScalar(thetaft));
	dPf_dVmt_dthetat    = Vmf*( Gft*PetscSinScalar(thetaft) - Bft*PetscCosScalar(thetaft));
	
	dPf_dVmt_dVmt	    = 0.;
	
	PetscScalar dQf_dthetaf_dthetaf,dQf_dVmf_dthetaf,dQf_dthetat_dthetaf,dQf_dVmt_dthetaf;
	PetscScalar dQf_dVmf_dVmf,dQf_dthetat_dVmf,dQf_dVmt_dVmf;
	PetscScalar dQf_dthetat_dthetat,dQf_dVmt_dthetat;
	PetscScalar dQf_dVmt_dVmt;
	
	dQf_dthetaf_dthetaf     = Vmf*Vmt*( Bft*PetscCosScalar(thetaft) - Gft*PetscSinScalar(thetaft));
	dQf_dVmf_dthetaf	=     Vmt*( Bft*PetscSinScalar(thetaft) + Gft*PetscCosScalar(thetaft));
	dQf_dthetat_dthetaf     = Vmf*Vmt*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
	dQf_dVmt_dthetaf	=     Vmf*( Bft*PetscSinScalar(thetaft) + Gft*PetscCosScalar(thetaft));
	
	dQf_dVmf_dVmf		= -2*Bff;
	dQf_dthetat_dVmf	=     Vmt*(-Bft*PetscSinScalar(thetaft) - Gft*PetscCosScalar(thetaft));
	dQf_dVmt_dVmf		=	  (-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
	
	dQf_dthetat_dthetat     = Vmf*Vmt*( Bft*PetscCosScalar(thetaft) - Gft*PetscSinScalar(thetaft));
	dQf_dVmt_dthetat	=				Vmf*(-Bft*PetscSinScalar(thetaft) - Gft*PetscCosScalar(thetaft));
	
	dQf_dVmt_dVmt		=				0.;

	row[0] 	= xlocf; 	
	row[1] 	= xlocf+1; 	col[0] 	= xlocf;
	val[0]   = lambda[gloc]*dPf_dthetaf_dthetaf 	+ lambda[gloc+1]*dQf_dthetaf_dthetaf;
	val[1] = lambda[gloc]*dPf_dVmf_dthetaf 		+ lambda[gloc+1]*dQf_dVmf_dthetaf;
	ierr = MatSetValues(H,2,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	row[0] 	= xlocf+1; 	col[0] 	= xlocf+1;
	val[0] = lambda[gloc]*dPf_dVmf_dVmf 			+ lambda[gloc+1]*dQf_dVmf_dVmf;
	ierr = MatSetValues(H,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);


	row[0] 	= xloct; 	
	col[0] 	= xlocf;
	col[1] 	= xlocf+1;
	col[2] 	= xloct;
	val[0] = lambda[gloc]*dPf_dthetat_dthetaf 	+ lambda[gloc+1]*dQf_dthetat_dthetaf;
	val[1] = lambda[gloc]*dPf_dthetat_dVmf		+ lambda[gloc+1]*dQf_dthetat_dVmf;
	val[2] = lambda[gloc]*dPf_dthetat_dthetat 	+ lambda[gloc+1]*dQf_dthetat_dthetat;
	ierr = MatSetValues(H,1,row,3,col,val,ADD_VALUES);CHKERRQ(ierr);

	row[0] 	= xloct+1;
 	col[0] 	= xlocf;
	col[1] 	= xlocf+1;
	col[2] 	= xloct;
	val[0] = lambda[gloc]*dPf_dVmt_dthetaf 		+ lambda[gloc+1]*dQf_dVmt_dthetaf;
	val[1] = lambda[gloc]*dPf_dVmt_dVmf 		+ lambda[gloc+1]*dQf_dVmt_dVmf;
	val[2] = lambda[gloc]*dPf_dVmt_dthetat 		+ lambda[gloc+1]*dQf_dVmt_dthetat;
	ierr = MatSetValues(H,1,row,3,col,val,ADD_VALUES);CHKERRQ(ierr);

      } else {
	
	PetscScalar dPt_dthetat_dthetat,dPt_dVmt_dthetat,dPt_dthetaf_dthetat,dPt_dVmf_dthetat;
	PetscScalar dPt_dVmt_dVmt,dPt_dthetaf_dVmt,dPt_dVmf_dVmt;
	PetscScalar dPt_dthetaf_dthetaf,dPt_dVmf_dthetaf;
	PetscScalar dPt_dVmf_dVmf;
	
	dPt_dthetat_dthetat = 		  Vmt*Vmf*(-Gtf*PetscCosScalar(thetatf) - Btf*PetscSinScalar(thetatf));
	dPt_dVmt_dthetat	= 			  Vmf*(-Gtf*PetscSinScalar(thetatf) + Btf*PetscCosScalar(thetatf));
	dPt_dthetaf_dthetat = 		  Vmt*Vmf*( Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
	dPt_dVmf_dthetat	= 			  Vmt*(-Gtf*PetscSinScalar(thetatf) + Btf*PetscCosScalar(thetatf));
	
	dPt_dVmt_dVmt 	  	= 2*Gtt;
	dPt_dthetaf_dVmt	= 			  Vmf*( Gtf*PetscSinScalar(thetatf) - Btf*PetscCosScalar(thetatf));
	dPt_dVmf_dVmt 	  	= 				  ( Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
	
	dPt_dthetaf_dthetaf = 		  Vmt*Vmf*(-Gtf*PetscCosScalar(thetatf) - Btf*PetscSinScalar(thetatf));
	dPt_dVmf_dthetaf	= 			  Vmt*( Gtf*PetscSinScalar(thetatf) - Btf*PetscCosScalar(thetatf));
	
	dPt_dVmf_dVmf 	  	= 			  0.;
	
	PetscScalar dQt_dthetat_dthetat,dQt_dVmt_dthetat,dQt_dthetaf_dthetat,dQt_dVmf_dthetat;
	PetscScalar dQt_dVmt_dVmt,dQt_dthetaf_dVmt,dQt_dVmf_dVmt;
	PetscScalar dQt_dthetaf_dthetaf,dQt_dVmf_dthetaf;
	PetscScalar dQt_dVmf_dVmf;  
	
	dQt_dthetat_dthetat = 		  Vmt*Vmf*( Btf*PetscCosScalar(thetatf) - Gtf*PetscSinScalar(thetatf));
	dQt_dVmt_dthetat	= 			  Vmf*( Btf*PetscSinScalar(thetatf) + Gtf*PetscCosScalar(thetatf));
	dQt_dthetaf_dthetat = 		  Vmt*Vmf*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
	dQt_dVmf_dthetat	= 			  Vmt*( Btf*PetscSinScalar(thetatf) + Gtf*PetscCosScalar(thetatf));
	
	dQt_dVmt_dVmt 	  	= -2*Btt;
	dQt_dthetaf_dVmt	= 			  Vmf*(-Btf*PetscSinScalar(thetatf) - Gtf*PetscCosScalar(thetatf));
	dQt_dVmf_dVmt 	  	= 				  (-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
	
	dQt_dthetaf_dthetaf = 		  Vmt*Vmf*( Btf*PetscCosScalar(thetatf) - Gtf*PetscSinScalar(thetatf));
	dQt_dVmf_dthetaf	= 			  Vmt*(-Btf*PetscSinScalar(thetatf) - Gtf*PetscCosScalar(thetatf));
	
	dQt_dVmf_dVmf 	  	= 			  0.;
	
	row[0]	= xlocf;
	row[1]  = xlocf+1;  col[0]  = xlocf;
	val[0]  = lambda[gloc]*dPt_dthetaf_dthetaf 	+ lambda[gloc+1]*dQt_dthetaf_dthetaf;
	val[1]  = lambda[gloc]*dPt_dVmf_dthetaf 	+ lambda[gloc+1]*dQt_dVmf_dthetaf;
	ierr = MatSetValues(H,2,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	row[0]  = xloct;
	col[0]  = xlocf;
	col[1]  = xlocf+1;
	col[2]  = xloct;
	val[0] = lambda[gloc]*dPt_dthetaf_dthetat 	+ lambda[gloc+1]*dQt_dthetaf_dthetat;
	val[1] = lambda[gloc]*dPt_dVmf_dthetat 		+ lambda[gloc+1]*dQt_dVmf_dthetat;
	val[2] = lambda[gloc]*dPt_dthetat_dthetat 	+ lambda[gloc+1]*dQt_dthetat_dthetat;
	ierr = MatSetValues(H,1,row,3,col,val,ADD_VALUES);CHKERRQ(ierr);

	row[0] = xloct + 1;
	col[0]  = xlocf;
	col[1]  = xlocf+1;
	col[2]  = xloct;
	col[3]  = xloct+1;
	val[0] = lambda[gloc]*dPt_dthetaf_dVmt 		+ lambda[gloc+1]*dQt_dthetaf_dVmt;
	val[1] = lambda[gloc]*dPt_dVmf_dVmt 			+ lambda[gloc+1]*dQt_dVmf_dVmt;
	val[2] = lambda[gloc]*dPt_dVmt_dthetat 		+ lambda[gloc+1]*dQt_dVmt_dthetat;
	val[3] = lambda[gloc]*dPt_dVmt_dVmt		 	+ lambda[gloc+1]*dQt_dVmt_dVmt;
	ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

      }
    }
    gloc += 2;
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
PetscErrorCode OPFLOWComputeInequalityConstraintsHessian_PBPOL(OPFLOW opflow, Vec X, Vec Lambda,Mat H)
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
PetscErrorCode OPFLOWComputeObjectiveHessian_PBPOL(OPFLOW opflow,Vec X,Mat H) 
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

PetscErrorCode OPFLOWComputeHessian_PBPOL(OPFLOW opflow,Vec X,Vec Lambda,Mat H)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = MatZeroEntries(H);CHKERRQ(ierr);

  /* Objective function Hessian */
  ierr = OPFLOWComputeObjectiveHessian_PBPOL(opflow,X,H);CHKERRQ(ierr);

  /* Equality constraints Hessian */
  ierr = OPFLOWComputeEqualityConstraintsHessian_PBPOL(opflow,X,Lambda,H);CHKERRQ(ierr);
  
  /* Inequality constraints Hessian */
  if(opflow->nconineq) {
    ierr = OPFLOWComputeInequalityConstraintsHessian_PBPOL(opflow,X,Lambda,H);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


PetscErrorCode OPFLOWFormulationCreate_PBPOL(OPFLOW opflow)
{
  PBPOL pbpol;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  
  ierr = PetscCalloc1(1,&pbpol);CHKERRQ(ierr);

  opflow->formulation = pbpol;

  /* Inherit Ops */
  opflow->formops.destroy = OPFLOWFormulationDestroy_PBPOL;
  opflow->formops.setnumvariables = OPFLOWFormulationSetNumVariables_PBPOL;
  opflow->formops.setnumconstraints = OPFLOWFormulationSetNumConstraints_PBPOL;
  opflow->formops.setvariablebounds = OPFLOWSetVariableBounds_PBPOL;
  opflow->formops.setconstraintbounds = OPFLOWSetConstraintBounds_PBPOL;
  opflow->formops.setvariableandconstraintbounds = OPFLOWSetVariableandConstraintBounds_PBPOL;
  opflow->formops.setinitialguess = OPFLOWSetInitialGuess_PBPOL;
  opflow->formops.computeequalityconstraints = OPFLOWComputeEqualityConstraints_PBPOL;
  opflow->formops.computeinequalityconstraints = OPFLOWComputeInequalityConstraints_PBPOL;
  opflow->formops.computeequalityconstraintjacobian = OPFLOWComputeEqualityConstraintJacobian_PBPOL;
  opflow->formops.computeinequalityconstraintjacobian = OPFLOWComputeInequalityConstraintJacobian_PBPOL;
  opflow->formops.computehessian = OPFLOWComputeHessian_PBPOL;
  opflow->formops.computeobjandgradient = OPFLOWComputeObjandGradient_PBPOL;
  opflow->formops.computeobjective = OPFLOWComputeObjective_PBPOL;
  opflow->formops.computegradient  = OPFLOWComputeGradient_PBPOL;
  
  PetscFunctionReturn(0);
}
