

#include <private/opflowimpl.h>
#include "pbpol.h"

PetscErrorCode OPFLOWModelDestroy_PBPOL(OPFLOW opflow)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscFree(opflow->model);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetVariableBounds_PBPOL(OPFLOW opflow,Vec Xl,Vec Xu)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscScalar    *xl,*xu;
  PetscInt       i;
  PSBUS          bus;
  PSGEN          gen;
  PetscInt       loc;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(Xu,&xu);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    PetscInt k;

    bus = &ps->bus[i];

    /* Bounds on bus variables */
    loc = bus->startxVloc;

    xl[loc] = PETSC_NINFINITY; xu[loc] = PETSC_INFINITY;

    if(!opflow->genbusVmfixed) {
      xl[loc+1] = bus->Vmin; xu[loc+1] = bus->Vmax;
    } else {
      if(bus->ide == REF_BUS || bus->ide == PV_BUS) {
	/* Hold bus voltage at generator set point voltage */
	xl[loc+1] = xu[loc+1] = bus->vm; /* Note: For pv buses bus->vm is already set to generator set point voltage when the data is read */
      } else {
	xl[loc+1] = bus->Vmin; xu[loc+1] = bus->Vmax; /* PQ buses */
      }
    }

    if(bus->ide == REF_BUS || bus->ide == ISOLATED_BUS) xl[loc] = xu[loc] = bus->va*PETSC_PI/180.0;
    if(bus->ide == ISOLATED_BUS) xl[loc+1] = xu[loc+1] = bus->vm;

    if(opflow->include_powerimbalance_variables) {
      loc = bus->startxpimbloc;

      xl[loc] = xl[loc+1] = PETSC_NINFINITY;
      xu[loc] = xu[loc+1] = PETSC_INFINITY;
      loc += bus->nxpimb;
    }

    /* Bounds on generator variables */

    for(k=0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;

      loc = gen->startxpowloc;

      xl[loc] = gen->pb; /* PGmin */
      xu[loc] = gen->pt; /* PGmax */
      xl[loc+1] = gen->qb; /* QGmin */
      xu[loc+1] = gen->qt; /* QGmax */
      /* pb, pt, qb, qt are converted in p.u. in ps.c */

      if(opflow->has_gensetpoint) {
	loc = gen->startxpdevloc;

	/* Lower and upper bounds for ramp up and down set to the
	   available headroom. This should be modified later on
	*/
	xl[loc] = xl[loc+1] = 0.0;
	xu[loc]   = gen->pt - gen->pgs;
	xu[loc+1] = gen->pgs - gen->pb; 
      }
    }

    if(opflow->include_loadloss_variables) {
      for(k=0; k < bus->nload; k++) {
	PSLOAD load;
	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
	loc = load->startxloadlossloc;
	if(!load->status) xl[loc] = xu[loc] = xl[loc+1] = xu[loc+1] = 0.0;
	else {
	  xl[loc] = PetscMin(0.0,load->pl);
	  xu[loc] = PetscMax(0.0,load->pl);
	  xl[loc+1] = PetscMin(0.0,load->ql);
	  xu[loc+1] = PetscMax(0.0,load->ql);
	}
      }
    }
  }

  ierr = VecRestoreArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Xu,&xu);CHKERRQ(ierr);

  PetscFunctionReturn(0);

}

PetscErrorCode OPFLOWSetConstraintBounds_PBPOL(OPFLOW opflow,Vec Gl,Vec Gu)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscScalar    *gl,*gu;
  PetscInt       i,k,ngen;
  PSBUS          bus;
  PSLINE         line;
  PSGEN          gen;
  PetscInt       gloc;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(Gu,&gu);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    gloc = bus->starteqloc;
    /* Equality constraint bounds for real and reactive power mismatch at the bus */
    gl[gloc] = 0.0;   gu[gloc] = 0.0;
    gl[gloc+1] = 0.0; gu[gloc+1] = 0.0;

    if(opflow->has_gensetpoint) {
      ierr = PSBUSGetNGen(bus,&ngen);CHKERRQ(ierr);
      for(k=0; k < ngen; k++) {
	ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	if(!gen->status) continue;
	gloc = gen->starteqloc;
	gl[gloc] = gu[gloc] = 0.0;
      }
    }
  }
  
  if(!opflow->ignore_lineflow_constraints) {
    for(i=0; i < ps->nline; i++) {
      line = &ps->line[i];
      if(!line->status || line->rateA > 1e5) continue;
      
      gloc = opflow->nconeq + line->startineqloc;
      /* Line flow inequality constraints */
      gl[gloc] = gl[gloc+1] = 0.0; 
      gu[gloc] = gu[gloc+1] = (line->rateA/ps->MVAbase)*(line->rateA/ps->MVAbase);    
    }
  }

  ierr = VecRestoreArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gu,&gu);CHKERRQ(ierr);

  PetscFunctionReturn(0);

}

PetscErrorCode OPFLOWSetVariableandConstraintBounds_PBPOL(OPFLOW opflow,Vec Xl,Vec Xu, Vec Gl, Vec Gu)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = OPFLOWSetVariableBounds_PBPOL(opflow,Xl,Xu);CHKERRQ(ierr);
  ierr = OPFLOWSetConstraintBounds_PBPOL(opflow,Gl,Gu);CHKERRQ(ierr);

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

    loc = bus->startxVloc;

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

    if(opflow->include_powerimbalance_variables) {
      loc = bus->startxpimbloc;
      x[loc] = x[loc+1] = 0.0;
    }

    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;

      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;

      loc = gen->startxpowloc;

      if(opflow->initializationtype == OPFLOWINIT_MIDPOINT || opflow->initializationtype == OPFLOWINIT_FLATSTART) {
	x[loc]   = 0.5*(xl[loc] + xu[loc]);
	x[loc+1] = 0.5*(xl[loc+1] + xu[loc+1]);
      } else if(opflow->initializationtype == OPFLOWINIT_FROMFILE || opflow->initializationtype == OPFLOWINIT_ACPF) {
	x[loc] = PetscMax(gen->pb,PetscMin(gen->pg,gen->pt));
	x[loc+1] = PetscMax(gen->qb,PetscMin(gen->qg,gen->qt));
      }
      loc += gen->nxpow;

      if(opflow->has_gensetpoint) {
	loc = gen->startxpdevloc;
	x[loc] = x[loc+1] = 0.0;
      }
    }

    if(opflow->include_loadloss_variables) {
      for(k=0; k < bus->nload; k++) {
	PSLOAD load;
	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
	if(!load->status) continue;
	loc = load->startxloadlossloc;
	/* Initial value for real and reactive power load loss */
	x[loc] = 0.0;
	x[loc+1] = 0.0;
	loc += load->nxloadloss;
      }
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
  PetscScalar    Pg,Qg,Pd,Qd;
  PetscScalar    Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  PetscScalar    Vmf,Vmt,thetaf,thetat,thetaft,thetatf;
  PetscScalar    Pf,Qf,Pt,Qt;
  PetscScalar    theta,Vm;
  PS             ps=opflow->ps;
  PSLOAD         load;
  PSLINE         line;
  PSBUS          bus,busf,bust;
  PSGEN          gen;
  const PSBUS    *connbuses;
  const PSLINE   *connlines;
  const PetscScalar *x;
  double flps=0.0;

  PetscFunctionBegin;
  ierr = VecSet(Ge,0.0);CHKERRQ(ierr);

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);

  for (i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    gloc = bus->starteqloc;

    row[0] = gloc; row[1] = row[0] + 1;

    xloc = bus->startxVloc;

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

    flps += 5.0;

    /* Power imbalance addition */
    if(opflow->include_powerimbalance_variables) {
      PetscScalar Pimb,Qimb;
      xloc = bus->startxpimbloc;
      Pimb = x[xloc];
      Qimb = x[xloc+1];
      val[0] = Pimb;
      val[1] = Qimb;
      ierr = VecSetValues(Ge,2,row,val,ADD_VALUES);CHKERRQ(ierr);

      flps += 2.0;
    }

    for (k=0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;

      xloc = gen->startxpowloc;

      Pg = x[xloc];
      Qg = x[xloc+1];

      val[0] = -Pg;
      val[1] = -Qg;
      ierr = VecSetValues(Ge,2,row,val,ADD_VALUES);CHKERRQ(ierr);

      flps += 2.0;
    }

    for (k=0; k < bus->nload; k++) {
      ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
      xloc = load->startxloadlossloc;
      if(opflow->include_loadloss_variables) {
	Pd = load->pl - x[xloc];
	Qd = load->ql - x[xloc+1];
	xloc += load->nxloadloss;
      } else {
	Pd = load->pl;
	Qd = load->ql;
      }

      val[0] = Pd;
      val[1] = Qd;
      ierr = VecSetValues(Ge,2,row,val,ADD_VALUES);CHKERRQ(ierr);
      flps += 2.0;
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

      xlocf = busf->startxVloc;
      xloct = bust->startxVloc;

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

	flps += 78.0;
      } else {
      	Pt = Gtt*Vmt*Vmt  + Vmt*Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf));
      	Qt = -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));

        val[0] = Pt;
        val[1] = Qt;
        ierr = VecSetValues(Ge,2,row,val,ADD_VALUES);CHKERRQ(ierr);
	flps += 78.0;
      }
    }

    if(opflow->has_gensetpoint) {
      for (k=0; k < bus->ngen; k++) {
	PetscScalar delPgup,delPgdown; /* Ramp up and down */
	ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	if(!gen->status) continue;
	
	Pg = x[gen->startxpowloc];

	delPgup = x[gen->startxpdevloc];
	delPgdown = x[gen->startxpdevloc + 1];

	gloc = gen->starteqloc;

	row[0] = gloc;
	val[0] = gen->pgs + delPgup - delPgdown - Pg;

	ierr = VecSetValues(Ge,1,row,val,ADD_VALUES);CHKERRQ(ierr);

	flps += 3.0;
      }
    }

  }
  ierr = VecAssemblyBegin(Ge);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(Ge);CHKERRQ(ierr);

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);

  PetscLogFlops(flps);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeEqualityConstraintJacobian_PBPOL(OPFLOW opflow,Vec X,Mat Je)
{
  PetscErrorCode ierr;
  PetscInt       i,k,row[2],col[4],gidx;
  PetscInt       nconnlines,locglob,loc,locglobf,locglobt,locf,loct;
  PetscScalar    Vm,val[8],Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  PetscScalar    thetaf,thetat,Vmf,Vmt,thetaft,thetatf;
  PS             ps=opflow->ps;
  PSBUS          bus;
  PSLINE         line;
  PSBUS          busf,bust;
  PSGEN          gen;
  PSLOAD         load;
  DM             networkdm=ps->networkdm;
  const PSLINE   *connlines;
  const PSBUS    *connbuses;
  const PetscScalar *xarr;

  PetscFunctionBegin;
  ierr = MatZeroEntries(Je);CHKERRQ(ierr);

  ierr = VecGetArrayRead(X,&xarr);CHKERRQ(ierr);

  for (i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    gidx = bus->starteqloc;

    row[0] = gidx; row[1] = row[0] + 1;

    loc = bus->startxVloc;
    locglob = bus->startxVlocglob;

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
    
    /* Power imbalance Jacobian terms */
    if(opflow->include_powerimbalance_variables) {
      locglob = bus->startxpimblocglob;

      val[0] = 1;
      col[0] = locglob;
      ierr = MatSetValues(Je,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

      col[0] = locglob + 1;
      ierr = MatSetValues(Je,1,row+1,1,col,val,ADD_VALUES);CHKERRQ(ierr);

    }

    for (k=0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;

      locglob = gen->startxpowlocglob;

      val[0] = -1;
      col[0] = locglob;
      ierr = MatSetValues(Je,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
      
      col[0] = locglob + 1;
      ierr = MatSetValues(Je,1,row+1,1,col,val,ADD_VALUES);CHKERRQ(ierr);
    }

    if(opflow->include_loadloss_variables) {
      for(k=0; k < bus->nload; k++) {
	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);

	locglob = load->startxloadlosslocglob;

	val[0] = -1;
	col[0] = locglob;
	ierr = MatSetValues(Je,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	col[0] = locglob + 1;
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

      locf = busf->startxVloc;
      loct = bust->startxVloc;

      locglobf = busf->startxVlocglob;
      locglobt = bust->startxVlocglob;

      thetaf = xarr[locf];  Vmf = xarr[locf+1];
      thetat = xarr[loct];  Vmt = xarr[loct+1];
      thetaft = thetaf - thetat;
      thetatf = thetat - thetaf;

      if (bus == busf) {
      	col[0] = locglobf; col[1] = locglobf+1; col[2] = locglobt; col[3] = locglobt+1;
	/* dPf_dthetaf */
      	val[0] = Vmf*Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft));
	/*dPf_dVmf */
      	val[1] = 2*Gff*Vmf + Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
	/*dPf_dthetat */
      	val[2] = Vmf*Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft));
	/* dPf_dVmt */
      	val[3] = Vmf*(Gft*cos(thetaft) + Bft*sin(thetaft));

	/* dQf_dthetaf */
        val[4] = Vmf*Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft));
	/* dQf_dVmf */
        val[5] = -2*Bff*Vmf + Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
	/* dQf_dthetat */
        val[6] = Vmf*Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft));
	/* dQf_dVmt */
        val[7] = Vmf*(-Bft*cos(thetaft) + Gft*sin(thetaft));
        ierr = MatSetValues(Je,2,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
      } else {
      	col[0] = locglobt; col[1] = locglobt+1; col[2] = locglobf; col[3] = locglobf+1;
	/* dPt_dthetat */
      	val[0] = Vmt*Vmf*(-Gtf*sin(thetatf) + Btf*cos(thetatf));
	/* dPt_dVmt */
      	val[1] = 2*Gtt*Vmt + Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf));
	/* dPt_dthetaf */
      	val[2] = Vmt*Vmf*(Gtf*sin(thetatf) - Btf*cos(thetatf));
	/* dPt_dVmf */
      	val[3] = Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf));

	/* dQt_dthetat */
        val[4] = Vmt*Vmf*(Btf*sin(thetatf) + Gtf*cos(thetatf));
	/* dQt_dVmt */
        val[5] = -2*Btt*Vmt + Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
	/* dQt_dthetaf */
        val[6] = Vmt*Vmf*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
	/* dQt_dVmf */
        val[7] = Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
        ierr = MatSetValues(Je,2,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
      }
    }

    if(opflow->has_gensetpoint) {
      ierr = PSBUSGetVariableGlobalLocation(bus,&locglob);CHKERRQ(ierr);

      for (k=0; k < bus->ngen; k++) {
	ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	if(!gen->status) continue;

	locglob = gen->startxpowlocglob;

	gidx = gen->starteqloc;
	row[0] = gidx;

	col[0] = locglob;
	val[0] = -1.0;

	locglob = gen->startxpdevlocglob;

	col[1] = locglob;
	col[2] = locglob + 1;
	val[1] = 1.0;
	val[2] = -1.0;

	ierr = MatSetValues(Je,1,row,3,col,val,ADD_VALUES);CHKERRQ(ierr);
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
  PetscInt       gloc;
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
  double         flps=0.0;

  PetscFunctionBegin;
  ierr = VecSet(Gi,0.0);CHKERRQ(ierr);

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Gi,&g);CHKERRQ(ierr);

  if(!opflow->ignore_lineflow_constraints) {
    for(i=0; i<ps->nline; i++) {
      line = &ps->line[i];
      if(!line->status || line->rateA > 1e5) continue;

      gloc = line->startineqloc;

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
      
      xlocf = busf->startxVloc;
      xloct = bust->startxVloc;

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
      
      flps += 158.0;
    }
  }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gi,&g);CHKERRQ(ierr);

  PetscLogFlops(flps);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeInequalityConstraintJacobian_PBPOL(OPFLOW opflow,Vec X,Mat Ji)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscInt       row[2],col[4];
  PetscInt       rstart,rend;
  PetscInt       gloc,xlocf,xloct;
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

  if(!opflow->ignore_lineflow_constraints) {
    for (i=0; i < ps->nline; i++) {
      line = &ps->line[i];
      if(!line->status || line->rateA > 1e5) continue;

      gloc = line->startineqloc;

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
      
      xlocf = busf->startxVloc;
      xloct = bust->startxVloc;

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
      
    }
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

PetscErrorCode OPFLOWComputeObjective_PBPOL(OPFLOW opflow,Vec X,PetscScalar *obj)
{
  PetscErrorCode ierr;
  const PetscScalar *x;
  PS             ps=opflow->ps;
  PetscInt       i;
  PSBUS          bus;
  PSGEN          gen;
  PetscInt       loc;
  PetscInt       k;
  PetscScalar    Pg;
  double         flps=0.0;

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);

  *obj = 0.0;
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    
    if(opflow->include_powerimbalance_variables) {
      loc = bus->startxpimbloc;
      PetscScalar Pimb,Qimb;
      Pimb = x[loc];
      Qimb = x[loc+1];
      *obj += opflow->powerimbalance_penalty*ps->MVAbase*ps->MVAbase*(Pimb*Pimb + Qimb*Qimb);
      flps += 7.0;
    }

    for(k=0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;

      if(opflow->objectivetype == OPFLOWOBJ_MIN_GEN_COST) {
	loc = gen->startxpowloc;
	Pg = x[loc]*ps->MVAbase;
	*obj += gen->cost_alpha*Pg*Pg + gen->cost_beta*Pg + gen->cost_gamma;
	flps += 7.0;
      } else if(opflow->objectivetype == OPFLOWOBJ_MIN_GENSETPOINT_DEVIATION) {
	loc = gen->startxpdevloc;
	PetscScalar delPgup,delPgdown;
	delPgup   = x[loc];
	delPgdown = x[loc+1];

	*obj += delPgup*delPgup + delPgdown*delPgdown;
	flps += 3.0;
      }
    }

    if(opflow->include_loadloss_variables) {
      PSLOAD load;
      PetscScalar Pdloss,Qdloss;
      for(k=0; k < bus->nload; k++) {
	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
	if(!load->status) continue;
	loc = load->startxloadlossloc;
	Pdloss = x[loc];
	Qdloss = x[loc+1];
	*obj += opflow->loadloss_penalty*ps->MVAbase*ps->MVAbase*(Pdloss*Pdloss + Qdloss*Qdloss);
	flps += 7.0;
      }
    }
  }
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);

  PetscLogFlops(flps);
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
  PSGEN          gen;
  PSLOAD         load;
  PetscInt       loc;
  PetscInt       k;
  PetscScalar    Pg;
  double         flps=0.0;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(grad,&df);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    if(opflow->include_powerimbalance_variables) {
      PetscScalar Pimb,Qimb;
      loc = bus->startxpimbloc;
      Pimb = x[loc];
      Qimb = x[loc+1];
      df[loc] = opflow->powerimbalance_penalty*ps->MVAbase*ps->MVAbase*2*Pimb;
      df[loc+1] = opflow->powerimbalance_penalty*ps->MVAbase*ps->MVAbase*2*Qimb;
      flps += 8.0;
    }

    for(k=0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;

      if(opflow->objectivetype == OPFLOWOBJ_MIN_GEN_COST) {
	loc = gen->startxpowloc;
	Pg = x[loc]*ps->MVAbase;
	df[loc] = ps->MVAbase*(2*gen->cost_alpha*Pg + gen->cost_beta);
	flps += 5.0;
      } else if(opflow->objectivetype == OPFLOWOBJ_MIN_GENSETPOINT_DEVIATION) {
	PetscScalar delPgup,delPgdown;
	loc = gen->startxpdevloc;
	delPgup   = x[loc];
	delPgdown = x[loc+1];

	df[loc]     = 2*delPgup;
	df[loc + 1] = 2*delPgdown; 
	flps += 2.0;
      }
    }
    
    if(opflow->include_loadloss_variables) {
      PSLOAD load;
      PetscScalar Pdloss,Qdloss;
      for(k=0; k < bus->nload; k++) {
	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
	if(!load->status) continue;
	loc = load->startxloadlossloc;
	Pdloss = x[loc];
	Qdloss = x[loc+1];
	df[loc]   = opflow->loadloss_penalty*ps->MVAbase*ps->MVAbase*2*Pdloss;
	df[loc+1] = opflow->loadloss_penalty*ps->MVAbase*ps->MVAbase*2*Qdloss; 
	flps += 8.0;
      }
    }
  }
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(grad,&df);CHKERRQ(ierr);

  PetscLogFlops(flps);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeObjandGradient_PBPOL(OPFLOW opflow,Vec X,PetscScalar *obj,Vec Grad)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = OPFLOWComputeObjective_PBPOL(opflow,X,obj);CHKERRQ(ierr);
  ierr = OPFLOWComputeGradient_PBPOL(opflow,X,Grad);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* Note: This function also sets the number of variables for each component, i.e., bus, gen, load, and line */
PetscErrorCode OPFLOWModelSetNumVariables_PBPOL(OPFLOW opflow,PetscInt *busnvar,PetscInt *branchnvar,PetscInt *nx)
{
  PetscInt i,ngen,nload,k;
  PS       ps=opflow->ps;
  PSBUS    bus;
  PSGEN    gen;
  PSLOAD   load;
  PSLINE   line;
  PetscErrorCode ierr;
  PetscBool isghost;

  PetscFunctionBegin;
  
  *nx = 0;
  /* No variables for the branches */
  for(i=0; i < ps->nline; i++) {
    line = &ps->line[i];
    branchnvar[i] = line->nx = 0;
    *nx += branchnvar[i];
  }

  /* Variables for the buses */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);
    //    if(isghost) continue;
    busnvar[i] = 2; /* 2 variables for the bus voltages */
    bus->nxV = busnvar[i];

    if(opflow->include_powerimbalance_variables) {
      busnvar[i] += 2;
      bus->nxpimb = 2;
    }

    bus->nxshunt = 0;

    bus->nx = bus->nxV + bus->nxshunt + bus->nxpimb;

    ierr = PSBUSGetNGen(bus,&ngen);CHKERRQ(ierr);
    for(k=0; k < ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      busnvar[i] += 2; /* (2 variables for voltage + Pg, Qg for each gen) */
      gen->nxpow = 2;
      if(opflow->has_gensetpoint) {
	gen->nxpdev = 2;
	busnvar[i] += 2;
      }
      gen->nx = gen->nxpow + gen->nxpdev;
    }

    if(opflow->include_loadloss_variables) {
      for(k=0; k < bus->nload; k++) {
	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
	if(!load->status) continue;
	/* Load loss variables..Real and imaginary part of the load loss */
	busnvar[i] += 2;
	load->nxloadloss = 2;
	load->nx = load->nxloadloss;
      }
    }

    if(!isghost) *nx += busnvar[i];
  }

  PetscFunctionReturn(0);
}

/* Note: This function also sets the number of constraints for each component, i.e., bus, gen, load, and line */
PetscErrorCode OPFLOWModelSetNumConstraints_PBPOL(OPFLOW opflow,PetscInt *branchnconeq,PetscInt *busnconeq,PetscInt *nconeq,PetscInt *nconineq)
{
  PetscInt i,k,ngen;
  PS       ps=opflow->ps;
  PSBUS    bus;
  PSLINE   line;
  PSGEN    gen;
  PetscErrorCode ierr;
  PetscBool isghost;

  PetscFunctionBegin;
  *nconeq = *nconineq = 0;

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);
    if(isghost) continue;
    *nconeq += 2;
    bus->nconeq = 2;
    bus->nconineq = 0;
    ierr = PSBUSGetNGen(bus,&ngen);CHKERRQ(ierr);
    for(k=0; k < ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      if(opflow->has_gensetpoint) {
	*nconeq += 1;
	gen->nconeq = 1;
      }
    }
  }

  if(!opflow->ignore_lineflow_constraints) {
    for(i=0; i < ps->nline; i++) {
      line = &ps->line[i];
      if(line->status && line->rateA < 1e5) {
	*nconineq += 2; /* Line flow constraints */
	line->nconineq = 2;
	line->nconeq = 0;
      }
    }
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
  PetscInt       gloc;
  PetscInt       row[16],col[16];
  PetscScalar    val[16];
  PSGEN          gen;
  PetscInt       ngen;

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(Lambda,&lambda);CHKERRQ(ierr);

  // For equality constraints (power flow) */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    gloc = bus->starteqloc;

    xloc = bus->startxVloc;
    
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

      xlocf = busf->startxVloc;
      xloct = bust->startxVloc;

      PetscScalar Vmf,Vmt,thetaf,thetat,thetaft,thetatf;
      
      thetaf  = x[xlocf];
      Vmf     = x[xlocf+1];
      thetat  = x[xloct];
      Vmt     = x[xloct+1];
      thetaft = thetaf - thetat;
      thetatf = thetat - thetaf;
    
      if(bus == busf) {
	
	PetscScalar dPf_dthetaf_dthetaf,dPf_dthetaf_dVmf,dPf_dthetaf_dthetat,dPf_dthetaf_dVmt;
	PetscScalar dPf_dVmf_dthetaf,   dPf_dVmf_dVmf,   dPf_dVmf_dthetat,   dPf_dVmf_dVmt;
	PetscScalar dPf_dthetat_dthetaf,dPf_dthetat_dVmf,dPf_dthetat_dthetat,dPf_dthetat_dVmt;
	PetscScalar dPf_dVmt_dthetaf,   dPf_dVmt_dVmf,   dPf_dVmt_dthetat,   dPf_dVmt_dVmt;

	/* dPf_dthetaf = Vmf*Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft)); */
	dPf_dthetaf_dthetaf = -Vmf*Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
	dPf_dthetaf_dVmf    =     Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft));
	dPf_dthetaf_dthetat =  Vmf*Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
	dPf_dthetaf_dVmt    =     Vmf*(-Gft*sin(thetaft) + Bft*cos(thetaft));

	/* dPf_Vmf  = 2*Gff*Vmf + Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
	dPf_dVmf_dthetaf    =  Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft)); 
	dPf_dVmf_dVmf       =  2*Gff;
	dPf_dVmf_dthetat    =  Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft));
	dPf_dVmf_dVmt       =      (Gft*cos(thetaft) + Bft*sin(thetaft));

      	/* dPf_dthetat = Vmf*Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft)); */
	dPf_dthetat_dthetaf = Vmf*Vmt*(Gft*cos(thetaft)  + Bft*sin(thetaft));
	dPf_dthetat_dVmf    =     Vmt*(Gft*sin(thetaft)  - Bft*cos(thetaft));
	dPf_dthetat_dthetat = Vmf*Vmt*(-Gft*cos(thetaft) - Bft*sin(thetaft));
	dPf_dthetat_dVmt    =     Vmf*(Gft*sin(thetaft)  - Bft*cos(thetaft));

	/* dPf_dVmt = Vmf*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
	dPf_dVmt_dthetaf    = Vmf*(-Gft*sin(thetaft) + Bft*cos(thetaft)); 
	dPf_dVmt_dVmf       =      (Gft*cos(thetaft) + Bft*sin(thetaft));
	dPf_dVmt_dthetat    = Vmf*(Gft*sin(thetaft) - Bft*cos(thetaft));
	dPf_dVmt_dVmt       = 0.0;
	
	PetscScalar dQf_dthetaf_dthetaf,dQf_dthetaf_dVmf,dQf_dthetaf_dthetat,dQf_dthetaf_dVmt;
	PetscScalar dQf_dVmf_dthetaf,   dQf_dVmf_dVmf,   dQf_dVmf_dthetat,   dQf_dVmf_dVmt;
	PetscScalar dQf_dthetat_dthetaf,dQf_dthetat_dVmf,dQf_dthetat_dthetat,dQf_dthetat_dVmt;
	PetscScalar dQf_dVmt_dthetaf,   dQf_dVmt_dVmf,   dQf_dVmt_dthetat,   dQf_dVmt_dVmt;

	/* dQf_dthetaf = Vmf*Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft)); */
	dQf_dthetaf_dthetaf = Vmf*Vmt*(Bft*cos(thetaft)  - Gft*sin(thetaft));
	dQf_dthetaf_dVmf    =     Vmt*(Bft*sin(thetaft)  + Gft*cos(thetaft));
	dQf_dthetaf_dthetat = Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
	dQf_dthetaf_dVmt    =     Vmf*(Bft*sin(thetaft) + Gft*cos(thetaft));

	/* dQf_dVmf = -2*Bff*Vmf + Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft)); */
	dQf_dVmf_dthetaf    =  Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft)); 
	dQf_dVmf_dVmf       = -2*Bff;
	dQf_dVmf_dthetat    =  Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft));
	dQf_dVmf_dVmt       =      (-Bft*cos(thetaft) + Gft*sin(thetaft));

	/* dQf_dthetat = Vmf*Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft)); */
	dQf_dthetat_dthetaf = Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
	dQf_dthetat_dVmf    =     Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft));
	dQf_dthetat_dthetat = Vmf*Vmt*( Bft*cos(thetaft) - Gft*sin(thetaft));
	dQf_dthetat_dVmt    =     Vmf*(-Bft*sin(thetaft) - Gft*cos(thetaft));

	/* dQf_dVmt = Vmf*(-Bft*cos(thetaft) + Gft*sin(thetaft)); */
	dQf_dVmt_dthetaf    = Vmf*(Bft*sin(thetaft) + Gft*cos(thetaft)); 
	dQf_dVmt_dVmf       =    (-Bft*cos(thetaft) + Gft*sin(thetaft));
	dQf_dVmt_dthetat    = Vmf*(-Bft*sin(thetaft) - Gft*cos(thetaft));
	dQf_dVmt_dVmt       = 0.0;

	row[0] = xlocf; row[1] = xlocf+1;
	col[0] = xlocf; col[1] = xlocf+1; col[2] = xloct; col[3] = xloct+1;

	val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.0;

	val[0] = lambda[gloc]*dPf_dthetaf_dthetaf + lambda[gloc+1]*dQf_dthetaf_dthetaf;
	val[1] = lambda[gloc]*dPf_dthetaf_dVmf    + lambda[gloc+1]*dQf_dthetaf_dVmf;
	val[2] = lambda[gloc]*dPf_dthetaf_dthetat + lambda[gloc+1]*dQf_dthetaf_dthetat;
	val[3] = lambda[gloc]*dPf_dthetaf_dVmt    + lambda[gloc+1]*dQf_dthetaf_dVmt;

	val[4] = lambda[gloc]*dPf_dVmf_dthetaf + lambda[gloc+1]*dQf_dVmf_dthetaf;
	val[5] = lambda[gloc]*dPf_dVmf_dVmf    + lambda[gloc+1]*dQf_dVmf_dVmf;
	val[6] = lambda[gloc]*dPf_dVmf_dthetat + lambda[gloc+1]*dQf_dVmf_dthetat;
	val[7] = lambda[gloc]*dPf_dVmf_dVmt    + lambda[gloc+1]*dQf_dVmf_dVmt;

	ierr = MatSetValues(H,2,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

	row[0] = xloct; row[1] = xloct+1;
	col[0] = xlocf; col[1] = xlocf+1; col[2] = xloct; col[3] = xloct+1;

	val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.0;

	val[0] = lambda[gloc]*dPf_dthetat_dthetaf + lambda[gloc+1]*dQf_dthetat_dthetaf;
	val[1] = lambda[gloc]*dPf_dthetat_dVmf    + lambda[gloc+1]*dQf_dthetat_dVmf;
	val[2] = lambda[gloc]*dPf_dthetat_dthetat + lambda[gloc+1]*dQf_dthetat_dthetat;
	val[3] = lambda[gloc]*dPf_dthetat_dVmt    + lambda[gloc+1]*dQf_dthetat_dVmt;

	val[4] = lambda[gloc]*dPf_dVmt_dthetaf + lambda[gloc+1]*dQf_dVmt_dthetaf;
	val[5] = lambda[gloc]*dPf_dVmt_dVmf    + lambda[gloc+1]*dQf_dVmt_dVmf;
	val[6] = lambda[gloc]*dPf_dVmt_dthetat + lambda[gloc+1]*dQf_dVmt_dthetat;
	val[7] = lambda[gloc]*dPf_dVmt_dVmt    + lambda[gloc+1]*dQf_dVmt_dVmt;

	ierr = MatSetValues(H,2,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

      } else {
	
	PetscScalar dPt_dthetat_dthetat,dPt_dthetat_dVmt,dPt_dthetat_dthetaf,dPt_dthetat_dVmf;
	PetscScalar dPt_dVmt_dthetat,   dPt_dVmt_dVmt,   dPt_dVmt_dthetaf,   dPt_dVmt_dVmf;
	PetscScalar dPt_dthetaf_dthetat,dPt_dthetaf_dVmt,dPt_dthetaf_dthetaf,dPt_dthetaf_dVmf;
	PetscScalar dPt_dVmf_dthetat,   dPt_dVmf_dVmt,   dPt_dVmf_dthetaf,   dPt_dVmf_dVmf;

	/* dPt_dthetat = Vmf*Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf)); */
	dPt_dthetat_dthetat = Vmf*Vmt*(-Gtf*cos(thetatf) - Btf*sin(thetatf));
	dPt_dthetat_dVmt    =     Vmf*(-Gtf*sin(thetatf) + Btf*cos(thetatf));
	dPt_dthetat_dthetaf =  Vmf*Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf));
	dPt_dthetat_dVmf    =     Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf));

	/* dPt_Vmt  = 2*Gtt*Vmt + Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
	dPt_dVmt_dthetat    =  Vmf*(-Gtf*sin(thetatf) + Bft*cos(thetatf)); 
	dPt_dVmt_dVmt       =  2*Gtt;
	dPt_dVmt_dthetaf    =  Vmf*(Gtf*sin(thetatf) - Btf*cos(thetatf));
	dPt_dVmt_dVmf       =      (Gtf*cos(thetatf) + Btf*sin(thetatf));

      	/* dPt_dthetaf = Vmf*Vmt*(Gtf*sin(thetatf) - Btf*cos(thetatf)); */
	dPt_dthetaf_dthetat = Vmf*Vmt*(Gtf*cos(thetatf)  + Btf*sin(thetatf));
	dPt_dthetaf_dVmt    =     Vmf*(Gtf*sin(thetatf)  - Btf*cos(thetatf));
	dPt_dthetaf_dthetaf = Vmf*Vmt*(-Gtf*cos(thetatf) - Btf*sin(thetatf));
	dPt_dthetaf_dVmf    =     Vmt*(Gtf*sin(thetatf)  - Btf*cos(thetatf));

	/* dPt_dVmf = Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
	dPt_dVmf_dthetat    = Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf)); 
	dPt_dVmf_dVmt       =     (Gtf*cos(thetatf) + Btf*sin(thetatf));
	dPt_dVmf_dthetaf    = Vmt*(Gtf*sin(thetatf) - Btf*cos(thetatf));
	dPt_dVmf_dVmf       = 0.0;
	
	PetscScalar dQt_dthetaf_dthetaf,dQt_dthetaf_dVmf,dQt_dthetaf_dthetat,dQt_dthetaf_dVmt;
	PetscScalar dQt_dVmf_dthetaf,   dQt_dVmf_dVmf,   dQt_dVmf_dthetat,   dQt_dVmf_dVmt;
	PetscScalar dQt_dthetat_dthetaf,dQt_dthetat_dVmf,dQt_dthetat_dthetat,dQt_dthetat_dVmt;
	PetscScalar dQt_dVmt_dthetaf,   dQt_dVmt_dVmf,   dQt_dVmt_dthetat,   dQt_dVmt_dVmt;

	/* dQt_dthetat = Vmf*Vmt*(Btf*sin(thetatf) + Gtf*cos(thetatf)); */
	dQt_dthetat_dthetat = Vmf*Vmt*(Btf*cos(thetatf)  - Gtf*sin(thetatf));
	dQt_dthetat_dVmt    =     Vmf*(Btf*sin(thetatf)  + Gtf*cos(thetatf));
	dQt_dthetat_dthetaf = Vmf*Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
	dQt_dthetat_dVmf    =     Vmf*(Btf*sin(thetatf) + Gtf*cos(thetatf));

	/* dQt_dVmt = -2*Btt*Vmt + Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf)); */
	dQt_dVmt_dthetat    =  Vmf*(Btf*sin(thetatf) + Gtf*cos(thetatf)); 
	dQt_dVmt_dVmt       = -2*Btt;
	dQt_dVmt_dthetaf    =  Vmf*(-Btf*sin(thetatf) + Gtf*cos(thetatf));
	dQt_dVmt_dVmf       =      (-Btf*cos(thetatf) + Gtf*sin(thetatf));

	/* dQt_dthetaf = Vmf*Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf)); */
	dQt_dthetaf_dthetat = Vmf*Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
	dQt_dthetaf_dVmt    =     Vmf*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
	dQt_dthetaf_dthetaf = Vmf*Vmt*( Btf*cos(thetatf) - Gtf*sin(thetatf));
	dQt_dthetaf_dVmf    =     Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf));

	/* dQt_dVmf = Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf)); */
	dQt_dVmf_dthetat    = Vmt*(Btf*sin(thetatf) + Gtf*cos(thetatf)); 
	dQt_dVmf_dVmt       =    (-Btf*cos(thetatf) + Gtf*sin(thetatf));
	dQt_dVmf_dthetaf    = Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
	dQt_dVmf_dVmf       = 0.0;

	row[0] = xloct; row[1] = xloct+1;
	col[0] = xloct; col[1] = xloct+1; col[2] = xlocf; col[3] = xlocf+1;

	val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.0;

	val[0] = lambda[gloc]*dPt_dthetat_dthetat + lambda[gloc+1]*dQt_dthetat_dthetat;
	val[1] = lambda[gloc]*dPt_dthetat_dVmt    + lambda[gloc+1]*dQt_dthetat_dVmt;
	val[2] = lambda[gloc]*dPt_dthetat_dthetaf + lambda[gloc+1]*dQt_dthetat_dthetaf;
	val[3] = lambda[gloc]*dPt_dthetat_dVmf    + lambda[gloc+1]*dQt_dthetat_dVmf;

	val[4] = lambda[gloc]*dPt_dVmt_dthetat + lambda[gloc+1]*dQt_dVmt_dthetat;
	val[5] = lambda[gloc]*dPt_dVmt_dVmt    + lambda[gloc+1]*dQt_dVmt_dVmt;
	val[6] = lambda[gloc]*dPt_dVmt_dthetaf + lambda[gloc+1]*dQt_dVmt_dthetaf;
	val[7] = lambda[gloc]*dPt_dVmt_dVmf    + lambda[gloc+1]*dQt_dVmt_dVmf;

	ierr = MatSetValues(H,2,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

	row[0] = xlocf; row[1] = xlocf+1;
	col[0] = xloct; col[1] = xloct+1; col[2] = xlocf; col[3] = xlocf+1;

	val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.0;

	val[0] = lambda[gloc]*dPt_dthetaf_dthetat + lambda[gloc+1]*dQt_dthetaf_dthetat;
	val[1] = lambda[gloc]*dPt_dthetaf_dVmt    + lambda[gloc+1]*dQt_dthetaf_dVmt;
	val[2] = lambda[gloc]*dPt_dthetaf_dthetaf + lambda[gloc+1]*dQt_dthetaf_dthetaf;
	val[3] = lambda[gloc]*dPt_dthetaf_dVmf    + lambda[gloc+1]*dQt_dthetaf_dVmf;

	val[4] = lambda[gloc]*dPt_dVmf_dthetat + lambda[gloc+1]*dQt_dVmf_dthetat;
	val[5] = lambda[gloc]*dPt_dVmf_dVmt    + lambda[gloc+1]*dQt_dVmf_dVmt;
	val[6] = lambda[gloc]*dPt_dVmf_dthetaf + lambda[gloc+1]*dQt_dVmf_dthetaf;
	val[7] = lambda[gloc]*dPt_dVmf_dVmf    + lambda[gloc+1]*dQt_dVmf_dVmf;

	ierr = MatSetValues(H,2,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
     }
    }
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
  PetscInt       xlocf,xloct,xlocglobf,xlocglobt;
  PSBUS          busf,bust;
  const PetscScalar *x;
  const PetscScalar *lambda;
  PetscInt       gloc;
  PetscInt       row[12],col[12];
  PetscScalar    val[12];

  PetscFunctionBegin;
      
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(Lambda,&lambda);CHKERRQ(ierr);

  // for the part of line constraints
  if(!opflow->ignore_lineflow_constraints) {
    for(i=0; i < ps->nline; i++) {
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

      gloc = line->startineqloc;

      xlocf = busf->startxVloc;
      xloct = bust->startxVloc;
      xlocglobf = busf->startxVlocglob;
      xlocglobt = bust->startxVlocglob;

      PetscScalar Vmf,Vmt,thetaf,thetat,thetaft,thetatf;
      
      thetaf  = x[xlocf];
      Vmf     = x[xlocf+1];
      thetat  = x[xloct];
      Vmt     = x[xloct+1];
      thetaft = thetaf - thetat;
      thetatf = thetat - thetaf;
      
      // Sf2 and St2 are the constraints	
      PetscScalar Pf,Qf,Pt,Qt;
      
      Pf =  Gff*Vmf*Vmf + Vmf*Vmt*( Gft*cos(thetaft) + Bft*sin(thetaft));
      Qf = -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
      
      Pt =  Gtt*Vmt*Vmt + Vmt*Vmf*( Gtf*cos(thetatf) + Btf*sin(thetatf));
      Qt = -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
      
      PetscScalar dSf2_dPf, dSf2_dQf, dSt2_dPt, dSt2_dQt;
      
      dSf2_dPf = 2.*Pf;
      dSf2_dQf = 2.*Qf;
      dSt2_dPt = 2.*Pt;
      dSt2_dQt = 2.*Qt;
      
      PetscScalar dPf_dthetaf,dPf_dVmf,dPf_dthetat,dPf_dVmt;
      PetscScalar dQf_dthetaf,dQf_dVmf,dQf_dthetat,dQf_dVmt;
      PetscScalar dPt_dthetaf,dPt_dVmf,dPt_dthetat,dPt_dVmt;
      PetscScalar dQt_dthetaf,dQt_dVmf,dQt_dthetat,dQt_dVmt;
      
      dPf_dthetaf = 			Vmf*Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft));
      dPf_dVmf    = 2.*Gff*Vmf + 	Vmt*( Gft*cos(thetaft) + Bft*sin(thetaft));
      dPf_dthetat = 			Vmf*Vmt*( Gft*sin(thetaft) - Bft*cos(thetaft));
      dPf_dVmt    = 				Vmf*( Gft*cos(thetaft) + Bft*sin(thetaft));
      
      dQf_dthetaf = 			Vmf*Vmt*( Bft*sin(thetaft) + Gft*cos(thetaft));
      dQf_dVmf    = -2.*Bff*Vmf + 	Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
      dQf_dthetat = 			Vmf*Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft));
      dQf_dVmt    = 				Vmf*(-Bft*cos(thetaft) + Gft*sin(thetaft));
      
      dPt_dthetat = 			Vmt*Vmf*(-Gtf*sin(thetatf) + Btf*cos(thetatf));
      dPt_dVmt    = 2.*Gtt*Vmt + 	Vmf*( Gtf*cos(thetatf) + Btf*sin(thetatf));
      dPt_dthetaf = 			Vmt*Vmf*( Gtf*sin(thetatf) - Btf*cos(thetatf));
      dPt_dVmf    = 				Vmt*( Gtf*cos(thetatf) + Btf*sin(thetatf));
      
      dQt_dthetat = 			Vmt*Vmf*( Btf*sin(thetatf) + Gtf*cos(thetatf));
      dQt_dVmt    = -2.*Btt*Vmt + 	Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
      dQt_dthetaf = 			Vmt*Vmf*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
      dQt_dVmf    = 				Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
      
      PetscScalar d2Pf_dthetaf_dthetaf,d2Pf_dthetaf_dVmf,d2Pf_dthetaf_dthetat,d2Pf_dthetaf_dVmt;
      PetscScalar d2Pf_dVmf_dthetaf,   d2Pf_dVmf_dVmf,   d2Pf_dVmf_dthetat,   d2Pf_dVmf_dVmt;
      PetscScalar d2Pf_dthetat_dthetaf,d2Pf_dthetat_dVmf,d2Pf_dthetat_dthetat,d2Pf_dthetat_dVmt;
      PetscScalar d2Pf_dVmt_dthetaf,   d2Pf_dVmt_dVmf,   d2Pf_dVmt_dthetat,   d2Pf_dVmt_dVmt;
      
      /* dPf_dthetaf = Vmf*Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft)); */
      d2Pf_dthetaf_dthetaf = -Vmf*Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
      d2Pf_dthetaf_dVmf    =     Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft));
      d2Pf_dthetaf_dthetat =  Vmf*Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
      d2Pf_dthetaf_dVmt    =     Vmf*(-Gft*sin(thetaft) + Bft*cos(thetaft));
      
      /* dPf_Vmf  = 2*Gff*Vmf + Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
      d2Pf_dVmf_dthetaf    =  Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft)); 
      d2Pf_dVmf_dVmf       =  2*Gff;
      d2Pf_dVmf_dthetat    =  Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft));
      d2Pf_dVmf_dVmt       =      (Gft*cos(thetaft) + Bft*sin(thetaft));
      
      /* dPf_dthetat = Vmf*Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft)); */
      d2Pf_dthetat_dthetaf = Vmf*Vmt*(Gft*cos(thetaft)  + Bft*sin(thetaft));
      d2Pf_dthetat_dVmf    =     Vmt*(Gft*sin(thetaft)  - Bft*cos(thetaft));
      d2Pf_dthetat_dthetat = Vmf*Vmt*(-Gft*cos(thetaft) - Bft*sin(thetaft));
      d2Pf_dthetat_dVmt    =     Vmf*(Gft*sin(thetaft)  - Bft*cos(thetaft));
      
      /* dPf_dVmt = Vmf*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
      d2Pf_dVmt_dthetaf    = Vmf*(-Gft*sin(thetaft) + Bft*cos(thetaft)); 
      d2Pf_dVmt_dVmf       =      (Gft*cos(thetaft) + Bft*sin(thetaft));
      d2Pf_dVmt_dthetat    = Vmf*(Gft*sin(thetaft) - Bft*cos(thetaft));
      d2Pf_dVmt_dVmt       = 0.0;
      
      PetscScalar d2Qf_dthetaf_dthetaf,d2Qf_dthetaf_dVmf,d2Qf_dthetaf_dthetat,d2Qf_dthetaf_dVmt;
      PetscScalar d2Qf_dVmf_dthetaf,   d2Qf_dVmf_dVmf,   d2Qf_dVmf_dthetat,   d2Qf_dVmf_dVmt;
      PetscScalar d2Qf_dthetat_dthetaf,d2Qf_dthetat_dVmf,d2Qf_dthetat_dthetat,d2Qf_dthetat_dVmt;
      PetscScalar d2Qf_dVmt_dthetaf,   d2Qf_dVmt_dVmf,   d2Qf_dVmt_dthetat,   d2Qf_dVmt_dVmt;
      
      /* dQf_dthetaf = Vmf*Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft)); */
      d2Qf_dthetaf_dthetaf = Vmf*Vmt*(Bft*cos(thetaft)  - Gft*sin(thetaft));
      d2Qf_dthetaf_dVmf    =     Vmt*(Bft*sin(thetaft)  + Gft*cos(thetaft));
      d2Qf_dthetaf_dthetat = Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
      d2Qf_dthetaf_dVmt    =     Vmf*(Bft*sin(thetaft) + Gft*cos(thetaft));
      
      /* dQf_dVmf = -2*Bff*Vmf + Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft)); */
      d2Qf_dVmf_dthetaf    =  Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft)); 
      d2Qf_dVmf_dVmf       = -2*Bff;
      d2Qf_dVmf_dthetat    =  Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft));
      d2Qf_dVmf_dVmt       =      (-Bft*cos(thetaft) + Gft*sin(thetaft));
      
      /* dQf_dthetat = Vmf*Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft)); */
      d2Qf_dthetat_dthetaf = Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
      d2Qf_dthetat_dVmf    =     Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft));
      d2Qf_dthetat_dthetat = Vmf*Vmt*( Bft*cos(thetaft) - Gft*sin(thetaft));
      d2Qf_dthetat_dVmt    =     Vmf*(-Bft*sin(thetaft) - Gft*cos(thetaft));
      
      /* dQf_dVmt = Vmf*(-Bft*cos(thetaft) + Gft*sin(thetaft)); */
      d2Qf_dVmt_dthetaf    = Vmf*(Bft*sin(thetaft) + Gft*cos(thetaft)); 
      d2Qf_dVmt_dVmf       =    (-Bft*cos(thetaft) + Gft*sin(thetaft));
      d2Qf_dVmt_dthetat    = Vmf*(-Bft*sin(thetaft) - Gft*cos(thetaft));
      d2Qf_dVmt_dVmt       = 0.0;
      
      PetscScalar d2Pt_dthetat_dthetat,d2Pt_dthetat_dVmt,d2Pt_dthetat_dthetaf,d2Pt_dthetat_dVmf;
      PetscScalar d2Pt_dVmt_dthetat,   d2Pt_dVmt_dVmt,   d2Pt_dVmt_dthetaf,   d2Pt_dVmt_dVmf;
      PetscScalar d2Pt_dthetaf_dthetat,d2Pt_dthetaf_dVmt,d2Pt_dthetaf_dthetaf,d2Pt_dthetaf_dVmf;
      PetscScalar d2Pt_dVmf_dthetat,   d2Pt_dVmf_dVmt,   d2Pt_dVmf_dthetaf,   d2Pt_dVmf_dVmf;
      
      /* dPt_dthetat = Vmf*Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf)); */
      d2Pt_dthetat_dthetat = Vmf*Vmt*(-Gtf*cos(thetatf) - Btf*sin(thetatf));
      d2Pt_dthetat_dVmt    =     Vmf*(-Gtf*sin(thetatf) + Btf*cos(thetatf));
      d2Pt_dthetat_dthetaf =  Vmf*Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf));
      d2Pt_dthetat_dVmf    =     Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf));
      
      /* dPt_Vmt  = 2*Gtt*Vmt + Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
      d2Pt_dVmt_dthetat    =  Vmf*(-Gtf*sin(thetatf) + Bft*cos(thetatf)); 
      d2Pt_dVmt_dVmt       =  2*Gtt;
      d2Pt_dVmt_dthetaf    =  Vmf*(Gtf*sin(thetatf) - Btf*cos(thetatf));
      d2Pt_dVmt_dVmf       =      (Gtf*cos(thetatf) + Btf*sin(thetatf));
      
      /* dPt_dthetaf = Vmf*Vmt*(Gtf*sin(thetatf) - Btf*cos(thetatf)); */
      d2Pt_dthetaf_dthetat = Vmf*Vmt*(Gtf*cos(thetatf)  + Btf*sin(thetatf));
      d2Pt_dthetaf_dVmt    =     Vmf*(Gtf*sin(thetatf)  - Btf*cos(thetatf));
      d2Pt_dthetaf_dthetaf = Vmf*Vmt*(-Gtf*cos(thetatf) - Btf*sin(thetatf));
      d2Pt_dthetaf_dVmf    =     Vmt*(Gtf*sin(thetatf)  - Btf*cos(thetatf));
      
      /* dPt_dVmf = Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
      d2Pt_dVmf_dthetat    = Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf)); 
      d2Pt_dVmf_dVmt       =     (Gtf*cos(thetatf) + Btf*sin(thetatf));
      d2Pt_dVmf_dthetaf    = Vmt*(Gtf*sin(thetatf) - Btf*cos(thetatf));
      d2Pt_dVmf_dVmf       = 0.0;
      
      PetscScalar d2Qt_dthetaf_dthetaf,d2Qt_dthetaf_dVmf,d2Qt_dthetaf_dthetat,d2Qt_dthetaf_dVmt;
      PetscScalar d2Qt_dVmf_dthetaf,   d2Qt_dVmf_dVmf,   d2Qt_dVmf_dthetat,   d2Qt_dVmf_dVmt;
      PetscScalar d2Qt_dthetat_dthetaf,d2Qt_dthetat_dVmf,d2Qt_dthetat_dthetat,d2Qt_dthetat_dVmt;
      PetscScalar d2Qt_dVmt_dthetaf,   d2Qt_dVmt_dVmf,   d2Qt_dVmt_dthetat,   d2Qt_dVmt_dVmt;
      
      /* dQt_dthetat = Vmf*Vmt*(Btf*sin(thetatf) + Gtf*cos(thetatf)); */
      d2Qt_dthetat_dthetat = Vmf*Vmt*(Btf*cos(thetatf)  - Gtf*sin(thetatf));
      d2Qt_dthetat_dVmt    =     Vmf*(Btf*sin(thetatf)  + Gtf*cos(thetatf));
      d2Qt_dthetat_dthetaf = Vmf*Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
      d2Qt_dthetat_dVmf    =     Vmf*(Btf*sin(thetatf) + Gtf*cos(thetatf));
      
      /* dQt_dVmt = -2*Btt*Vmt + Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf)); */
      d2Qt_dVmt_dthetat    =  Vmf*(Btf*sin(thetatf) + Gtf*cos(thetatf)); 
      d2Qt_dVmt_dVmt       = -2*Btt;
      d2Qt_dVmt_dthetaf    =  Vmf*(-Btf*sin(thetatf) + Gtf*cos(thetatf));
      d2Qt_dVmt_dVmf       =      (-Btf*cos(thetatf) + Gtf*sin(thetatf));
      
      /* dQt_dthetaf = Vmf*Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf)); */
      d2Qt_dthetaf_dthetat = Vmf*Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
      d2Qt_dthetaf_dVmt    =     Vmf*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
      d2Qt_dthetaf_dthetaf = Vmf*Vmt*( Btf*cos(thetatf) - Gtf*sin(thetatf));
      d2Qt_dthetaf_dVmf    =     Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
      
      /* dQt_dVmf = Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf)); */
      d2Qt_dVmf_dthetat    = Vmt*(Btf*sin(thetatf) + Gtf*cos(thetatf)); 
      d2Qt_dVmf_dVmt       =    (-Btf*cos(thetatf) + Gtf*sin(thetatf));
      d2Qt_dVmf_dthetaf    = Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
      d2Qt_dVmf_dVmf       = 0.0;
      
      PetscScalar d2Sf2_dthetaf_dthetaf=0.0,d2Sf2_dthetaf_dVmf=0.0,d2Sf2_dthetaf_dthetat=0.0,d2Sf2_dthetaf_dVmt=0.0;
      PetscScalar d2St2_dthetaf_dthetaf=0.0,d2St2_dthetaf_dVmf=0.0,d2St2_dthetaf_dthetat=0.0,d2St2_dthetaf_dVmt=0.0;
      
      d2Sf2_dthetaf_dthetaf = 2*dPf_dthetaf*dPf_dthetaf + dSf2_dPf*d2Pf_dthetaf_dthetaf +  2*dQf_dthetaf*dQf_dthetaf + dSf2_dQf*d2Qf_dthetaf_dthetaf;
      d2Sf2_dthetaf_dVmf = 2*dPf_dVmf*dPf_dthetaf + dSf2_dPf*d2Pf_dthetaf_dVmf +  2*dQf_dVmf*dQf_dthetaf + dSf2_dQf*d2Qf_dthetaf_dVmf;
      d2Sf2_dthetaf_dthetat = 2*dPf_dthetat*dPf_dthetaf + dSf2_dPf*d2Pf_dthetaf_dthetat +  2*dQf_dthetat*dQf_dthetaf + dSf2_dQf*d2Qf_dthetaf_dthetat;
      d2Sf2_dthetaf_dVmt = 2*dPf_dVmt*dPf_dthetaf + dSf2_dPf*d2Pf_dthetaf_dVmt +  2*dQf_dVmt*dQf_dthetaf + dSf2_dQf*d2Qf_dthetaf_dVmt;
      
      d2St2_dthetaf_dthetaf = 2*dPt_dthetaf*dPt_dthetaf + dSt2_dPt*d2Pt_dthetaf_dthetaf +  2*dQt_dthetaf*dQt_dthetaf + dSt2_dQt*d2Qt_dthetaf_dthetaf;
      d2St2_dthetaf_dVmf = 2*dPt_dVmf*dPt_dthetaf + dSt2_dPt*d2Pt_dthetaf_dVmf +  2*dQt_dVmf*dQt_dthetaf + dSt2_dQt*d2Qt_dthetaf_dVmf;
      d2St2_dthetaf_dthetat = 2*dPt_dthetat*dPt_dthetaf + dSt2_dPt*d2Pt_dthetaf_dthetat +  2*dQt_dthetat*dQt_dthetaf + dSt2_dQt*d2Qt_dthetaf_dthetat;
      d2St2_dthetaf_dVmt = 2*dPt_dVmt*dPt_dthetaf + dSt2_dPt*d2Pt_dthetaf_dVmt +  2*dQt_dVmt*dQt_dthetaf + dSt2_dQt*d2Qt_dthetaf_dVmt;
      
      val[0] = val[1] = val[2] = val[3] = 0.0;
      col[0] = xlocglobf; 
      col[1] = xlocglobf+1; 
      col[2] = xlocglobt;
      col[3] = xlocglobt+1;
      row[0] = xlocglobf;
      val[0] = lambda[gloc]*d2Sf2_dthetaf_dthetaf + lambda[gloc+1]*d2St2_dthetaf_dthetaf;
      val[1] = lambda[gloc]*d2Sf2_dthetaf_dVmf + lambda[gloc+1]*d2St2_dthetaf_dVmf;
      val[2] = lambda[gloc]*d2Sf2_dthetaf_dthetat + lambda[gloc+1]*d2St2_dthetaf_dthetat;
      val[3] = lambda[gloc]*d2Sf2_dthetaf_dVmt + lambda[gloc+1]*d2St2_dthetaf_dVmt;
      
      ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
      
      PetscScalar d2Sf2_dVmf_dthetaf,d2Sf2_dVmf_dVmf,d2Sf2_dVmf_dthetat,d2Sf2_dVmf_dVmt;
      PetscScalar d2St2_dVmf_dthetaf,d2St2_dVmf_dVmf,d2St2_dVmf_dthetat,d2St2_dVmf_dVmt;
      
      d2Sf2_dVmf_dthetaf = 2*dPf_dthetaf*dPf_dVmf + dSf2_dPf*d2Pf_dVmf_dthetaf +  2*dQf_dthetaf*dQf_dVmf + dSf2_dQf*d2Qf_dVmf_dthetaf;
      d2Sf2_dVmf_dVmf = 2*dPf_dVmf*dPf_dVmf + dSf2_dPf*d2Pf_dVmf_dVmf +  2*dQf_dVmf*dQf_dVmf + dSf2_dQf*d2Qf_dVmf_dVmf;
      d2Sf2_dVmf_dthetat = 2*dPf_dthetat*dPf_dVmf + dSf2_dPf*d2Pf_dVmf_dthetat +  2*dQf_dthetat*dQf_dVmf + dSf2_dQf*d2Qf_dVmf_dthetat;
      d2Sf2_dVmf_dVmt = 2*dPf_dVmt*dPf_dVmf + dSf2_dPf*d2Pf_dVmf_dVmt +  2*dQf_dVmt*dQf_dVmf + dSf2_dQf*d2Qf_dVmf_dVmt;
      
      d2St2_dVmf_dthetaf = 2*dPt_dthetaf*dPt_dVmf + dSt2_dPt*d2Pt_dVmf_dthetaf +  2*dQt_dthetaf*dQt_dVmf + dSt2_dQt*d2Qt_dVmf_dthetaf;
      d2St2_dVmf_dVmf = 2*dPt_dVmf*dPt_dVmf + dSt2_dPt*d2Pt_dVmf_dVmf +  2*dQt_dVmf*dQt_dVmf + dSt2_dQt*d2Qt_dVmf_dVmf;
      d2St2_dVmf_dthetat = 2*dPt_dthetat*dPt_dVmf + dSt2_dPt*d2Pt_dVmf_dthetat +  2*dQt_dthetat*dQt_dVmf + dSt2_dQt*d2Qt_dVmf_dthetat;
      d2St2_dVmf_dVmt = 2*dPt_dVmt*dPt_dVmf + dSt2_dPt*d2Pt_dVmf_dVmt +  2*dQt_dVmt*dQt_dVmf + dSt2_dQt*d2Qt_dVmf_dVmt;
      
      val[0] = val[1] = val[2] = val[3] = 0.0;
      col[0] = xlocglobf; 
      col[1] = xlocglobf+1; 
      col[2] = xlocglobt;
      col[3] = xlocglobt+1;
      row[0] = xlocglobf+1;
      val[0] = lambda[gloc]*d2Sf2_dVmf_dthetaf + lambda[gloc+1]*d2St2_dVmf_dthetaf;
      val[1] = lambda[gloc]*d2Sf2_dVmf_dVmf + lambda[gloc+1]*d2St2_dVmf_dVmf;
      val[2] = lambda[gloc]*d2Sf2_dVmf_dthetat + lambda[gloc+1]*d2St2_dVmf_dthetat;
      val[3] = lambda[gloc]*d2Sf2_dVmf_dVmt + lambda[gloc+1]*d2St2_dVmf_dVmt;
      ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
      
      PetscScalar d2Sf2_dthetat_dthetaf,d2Sf2_dthetat_dVmf,d2Sf2_dthetat_dthetat,d2Sf2_dthetat_dVmt;
      PetscScalar d2St2_dthetat_dthetaf,d2St2_dthetat_dVmf,d2St2_dthetat_dthetat,d2St2_dthetat_dVmt;
      
      d2Sf2_dthetat_dthetaf = 2*dPf_dthetaf*dPf_dthetat + dSf2_dPf*d2Pf_dthetat_dthetaf +  2*dQf_dthetat*dQf_dthetaf + dSf2_dQf*d2Qf_dthetat_dthetaf;
      d2Sf2_dthetat_dVmf = 2*dPf_dVmf*dPf_dthetat + dSf2_dPf*d2Pf_dthetat_dVmf +  2*dQf_dthetat*dQf_dVmf + dSf2_dQf*d2Qf_dthetat_dVmf;
      d2Sf2_dthetat_dthetat = 2*dPf_dthetat*dPf_dthetat + dSf2_dPf*d2Pf_dthetat_dthetat +  2*dQf_dthetat*dQf_dthetat + dSf2_dQf*d2Qf_dthetat_dthetat;
      d2Sf2_dthetat_dVmt = 2*dPf_dVmt*dPf_dthetat + dSf2_dPf*d2Pf_dthetat_dVmt +  2*dQf_dthetat*dQf_dVmt + dSf2_dQf*d2Qf_dthetat_dVmt;
      
      d2St2_dthetat_dthetaf = 2*dPt_dthetaf*dPt_dthetat + dSt2_dPt*d2Pt_dthetat_dthetaf +  2*dQt_dthetaf*dQt_dthetat + dSt2_dQt*d2Qt_dthetat_dthetaf;
      d2St2_dthetat_dVmf = 2*dPt_dVmf*dPt_dthetat + dSt2_dPt*d2Pt_dthetat_dVmf +  2*dQt_dVmf*dQt_dthetat + dSt2_dQt*d2Qt_dthetat_dVmf;
      d2St2_dthetat_dthetat = 2*dPt_dthetat*dPt_dthetat + dSt2_dPt*d2Pt_dthetat_dthetat +  2*dQt_dthetat*dQt_dthetat + dSt2_dQt*d2Qt_dthetat_dthetat;
      d2St2_dthetat_dVmt = 2*dPt_dVmt*dPt_dthetat + dSt2_dPt*d2Pt_dthetat_dVmt +  2*dQt_dVmt*dQt_dthetat + dSt2_dQt*d2Qt_dthetat_dVmt;
      
      val[0] = val[1] = val[2] = val[3] = 0.0;
      col[0] = xlocglobf; 
      col[1] = xlocglobf+1; 
      col[2] = xlocglobt;
      col[3] = xlocglobt+1;
      row[0] = xlocglobt;
      
      val[0] = lambda[gloc]*d2Sf2_dthetat_dthetaf + lambda[gloc+1]*d2St2_dthetat_dthetaf;
      val[1] = lambda[gloc]*d2Sf2_dthetat_dVmf + lambda[gloc+1]*d2St2_dthetat_dVmf;
      val[2] = lambda[gloc]*d2Sf2_dthetat_dthetat + lambda[gloc+1]*d2St2_dthetat_dthetat;
      val[3] = lambda[gloc]*d2Sf2_dthetat_dVmt + lambda[gloc+1]*d2St2_dthetat_dVmt;
      
      ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
      
      PetscScalar d2Sf2_dVmt_dthetaf,d2Sf2_dVmt_dVmf,d2Sf2_dVmt_dthetat,d2Sf2_dVmt_dVmt;
      PetscScalar d2St2_dVmt_dthetaf,d2St2_dVmt_dVmf,d2St2_dVmt_dthetat,d2St2_dVmt_dVmt;
      
      d2Sf2_dVmt_dthetaf = 2*dPf_dthetaf*dPf_dVmt + dSf2_dPf*d2Pf_dVmt_dthetaf +  2*dQf_dthetaf*dQf_dVmt + dSf2_dQf*d2Qf_dVmt_dthetaf;
      d2Sf2_dVmt_dVmf = 2*dPf_dVmf*dPf_dVmt + dSf2_dPf*d2Pf_dVmt_dVmf +  2*dQf_dVmf*dQf_dVmt + dSf2_dQf*d2Qf_dVmt_dVmf;
      d2Sf2_dVmt_dthetat = 2*dPf_dthetat*dPf_dVmt + dSf2_dPf*d2Pf_dVmt_dthetat +  2*dQf_dthetat*dQf_dVmt + dSf2_dQf*d2Qf_dVmt_dthetat;
      d2Sf2_dVmt_dVmt = 2*dPf_dVmt*dPf_dVmt + dSf2_dPf*d2Pf_dVmt_dVmt +  2*dQf_dVmt*dQf_dVmt + dSf2_dQf*d2Qf_dVmt_dVmt;
      
      d2St2_dVmt_dthetaf = 2*dPt_dthetaf*dPt_dVmt + dSt2_dPt*d2Pt_dVmt_dthetaf +  2*dQt_dthetaf*dQt_dVmt + dSt2_dQt*d2Qt_dVmt_dthetaf;
      d2St2_dVmt_dVmf = 2*dPt_dVmf*dPt_dVmt + dSt2_dPt*d2Pt_dVmt_dVmf +  2*dQt_dVmf*dQt_dVmt + dSt2_dQt*d2Qt_dVmt_dVmf;
      d2St2_dVmt_dthetat = 2*dPt_dthetat*dPt_dVmt + dSt2_dPt*d2Pt_dVmt_dthetat +  2*dQt_dthetat*dQt_dVmt + dSt2_dQt*d2Qt_dVmt_dthetat;
      d2St2_dVmt_dVmt = 2*dPt_dVmt*dPt_dVmt + dSt2_dPt*d2Pt_dVmt_dVmt +  2*dQt_dVmt*dQt_dVmt + dSt2_dQt*d2Qt_dVmt_dVmt;
      
      val[0] = val[1] = val[2] = val[3] = 0.0;
      col[0] = xlocglobf; 
      col[1] = xlocglobf+1; 
      col[2] = xlocglobt;
      col[3] = xlocglobt+1;
      row[0] = xlocglobt+1;
      
      val[0] = lambda[gloc]*d2Sf2_dVmt_dthetaf + lambda[gloc+1]*d2St2_dVmt_dthetaf;
      val[1] = lambda[gloc]*d2Sf2_dVmt_dVmf + lambda[gloc+1]*d2St2_dVmt_dVmf;
      val[2] = lambda[gloc]*d2Sf2_dVmt_dthetat + lambda[gloc+1]*d2St2_dVmt_dthetat;
      val[3] = lambda[gloc]*d2Sf2_dVmt_dVmt + lambda[gloc+1]*d2St2_dVmt_dVmt;
      
      ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
    }
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
  PetscInt       xlocglob;
  const PetscScalar *x;
  PetscInt       row[2],col[2];
  PetscScalar    val[2];
  PetscScalar    obj_factor = opflow->obj_factor;


  PetscFunctionBegin;

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);

  // for the part of objective
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableGlobalLocation(bus,&xlocglob);CHKERRQ(ierr);

    if(opflow->include_powerimbalance_variables) {
      xlocglob = bus->startxpimblocglob;

      row[0] = xlocglob;
      col[0] = xlocglob;
      val[0] = obj_factor*2.0*opflow->powerimbalance_penalty*ps->MVAbase*ps->MVAbase;
      ierr = MatSetValues(H,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
	
      row[0] = xlocglob+1;
      col[0] = xlocglob+1;
      val[0] = obj_factor*2.0*opflow->powerimbalance_penalty*ps->MVAbase*ps->MVAbase;
      ierr = MatSetValues(H,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

    }

    for(k=0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;

      if(opflow->objectivetype == OPFLOWOBJ_MIN_GEN_COST) {
	xlocglob = gen->startxpowlocglob;

	row[0] = xlocglob;
	col[0] = xlocglob;

	val[0] = obj_factor*2.0*gen->cost_alpha*ps->MVAbase*ps->MVAbase;
	ierr = MatSetValues(H,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
      } else if(opflow->objectivetype == OPFLOWOBJ_MIN_GENSETPOINT_DEVIATION) {
	xlocglob = gen->startxpdevlocglob;
	row[0] = xlocglob;
	col[0] = xlocglob;
	val[0] = obj_factor*2.0;
	ierr = MatSetValues(H,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	row[0] = xlocglob + 1;
	col[0] = xlocglob + 1;
	val[0] = obj_factor*2.0;
	ierr = MatSetValues(H,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
      }
    }
    
    if(opflow->include_loadloss_variables) {
      PSLOAD load;
      for(k=0; k < bus->nload; k++) {
	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
	if(!load->status) continue;
	xlocglob = load->startxloadlosslocglob;
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
  }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeHessian_PBPOL(OPFLOW opflow,Vec X,Vec Lambdae,Vec Lambdai,Mat H)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = MatZeroEntries(H);CHKERRQ(ierr);

  /* Objective function Hessian */
  ierr = OPFLOWComputeObjectiveHessian_PBPOL(opflow,X,H);CHKERRQ(ierr);

  /* Equality constraints Hessian */
  ierr = OPFLOWComputeEqualityConstraintsHessian_PBPOL(opflow,X,Lambdae,H);CHKERRQ(ierr);
  
  /* Inequality constraints Hessian */
  if(opflow->nconineq) {
    ierr = OPFLOWComputeInequalityConstraintsHessian_PBPOL(opflow,X,Lambdai,H);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolutionToPS_PBPOL(OPFLOW opflow)
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

    loc = bus->startxVloc;

    bus->va = x[loc];
    bus->vm = x[loc+1];

    gloc = bus->starteqloc;
    bus->mult_pmis = lambdae[gloc];
    bus->mult_qmis = lambdae[gloc+1];

    if(opflow->include_powerimbalance_variables) {
      loc = bus->startxpimbloc;
      bus->pimb = x[loc];
      bus->qimb = x[loc+1];
    }


    for(k=0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) {
	gen->pg = gen->qg = 0.0;
	continue;
      }
      loc = gen->startxpowloc;

      gen->pg = x[loc];
      gen->qg = x[loc+1];

      if(opflow->has_gensetpoint) {
	gloc += gen->nconeq;
      }
    }

    if(opflow->include_loadloss_variables) {
      for(k=0; k < bus->nload; k++) {
	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
	loc = load->startxloadlossloc;
	load->pl = load->pl - x[loc];
	load->ql = load->ql - x[loc+1];
      }
    }
  }

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
      
      xlocf = busf->startxVloc;
      xloct = bust->startxVloc;

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
	gloc = opflow->nconeq + line->startineqloc;
	line->mult_sf = lambdai[gloc];
	line->mult_st = lambdai[gloc+1];
      }
    }
  }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Lambda,&lambda);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelSetUp_PBPOL(OPFLOW opflow)
{
  PetscErrorCode ierr;
  PS             ps=(PS)opflow->ps;
  PSBUS          bus;
  PSGEN          gen;
  PSLOAD         load;
  PSLINE         line;
  PetscInt       i,k;
  PetscInt       loc,locglob;
  PetscInt       eqloc=0,ineqloc=0;

  PetscFunctionBegin;

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bus,&locglob);CHKERRQ(ierr);

    bus->startxVloc = loc;
    bus->startxVlocglob = locglob;

    bus->startxshuntloc = bus->startxVloc + bus->nxV;
    bus->startxshuntlocglob = bus->startxVlocglob + bus->nxV;

    bus->startxpimbloc = bus->startxshuntloc + bus->nxshunt;
    bus->startxpimblocglob = bus->startxshuntlocglob + bus->nxshunt;

    bus->nx = bus->nxV + bus->nxshunt + bus->nxpimb;

    loc     += bus->nx;
    locglob += bus->nx;

    /* Equality constraints */
    bus->starteqloc = eqloc;
    eqloc += bus->nconeq;
 
    for(k=0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;

      gen->startxpowloc     = loc;
      gen->startxpowlocglob = locglob;
  
      gen->startxpdevloc     = gen->startxpowloc + gen->nxpow;
      gen->startxpdevlocglob = gen->startxpowlocglob + gen->nxpow;

      if(opflow->has_gensetpoint) {
	gen->starteqloc = eqloc;
	eqloc += gen->nconeq;
      }

      gen->nx = gen->nxpow + gen->nxpdev;

      loc += gen->nx;
      locglob += gen->nx;
    }

    for(k=0; k < bus->nload; k++) {
      ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
      if(!load->status) continue;
      if(opflow->include_loadloss_variables) {
	load->startxloadlossloc     = loc;
	load->startxloadlosslocglob = locglob;

	load->nx = load->nxloadloss;

	loc     += load->nx;
	locglob += load->nx;
      }
    }
  }

  for(i=0; i < ps->nline; i++) {
    line = &ps->line[i];
    if(line->status && line->rateA < 1e5) {
      line->startineqloc = ineqloc;
      ineqloc += line->nconineq;
    }
  }

     
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelCreate_PBPOL(OPFLOW opflow)
{
  PBPOL pbpol;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  
  ierr = PetscCalloc1(1,&pbpol);CHKERRQ(ierr);

  opflow->model = pbpol;

  /* Inherit Ops */
  opflow->modelops.destroy = OPFLOWModelDestroy_PBPOL;
  opflow->modelops.setnumvariables = OPFLOWModelSetNumVariables_PBPOL;
  opflow->modelops.setnumconstraints = OPFLOWModelSetNumConstraints_PBPOL;
  opflow->modelops.setvariablebounds = OPFLOWSetVariableBounds_PBPOL;
  opflow->modelops.setconstraintbounds = OPFLOWSetConstraintBounds_PBPOL;
  opflow->modelops.setvariableandconstraintbounds = OPFLOWSetVariableandConstraintBounds_PBPOL;
  opflow->modelops.setinitialguess = OPFLOWSetInitialGuess_PBPOL;
  opflow->modelops.setup = OPFLOWModelSetUp_PBPOL;
  opflow->modelops.computeequalityconstraints = OPFLOWComputeEqualityConstraints_PBPOL;
  opflow->modelops.computeinequalityconstraints = OPFLOWComputeInequalityConstraints_PBPOL;
  opflow->modelops.computeequalityconstraintjacobian = OPFLOWComputeEqualityConstraintJacobian_PBPOL;
  opflow->modelops.computeinequalityconstraintjacobian = OPFLOWComputeInequalityConstraintJacobian_PBPOL;
  opflow->modelops.computehessian = OPFLOWComputeHessian_PBPOL;
  opflow->modelops.computeobjandgradient = OPFLOWComputeObjandGradient_PBPOL;
  opflow->modelops.computeobjective = OPFLOWComputeObjective_PBPOL;
  opflow->modelops.computegradient  = OPFLOWComputeGradient_PBPOL;
  opflow->modelops.solutiontops     = OPFLOWSolutionToPS_PBPOL;
  
  PetscFunctionReturn(0);
}
