#include <private/psimpl.h>
#include <private/opflowimpl.h>
#include <private/scopflowpipsimpl.h>
#include <math.h>
#include <../src/mat/impls/aij/seq/aij.h>

/* 
   All functions related to SCOPFLOW constraints
   are in this file
*/

/*
  OPFLOWEqualityConstraints - Evalulates the equality constraints for the optimal power flow

  Input Parameters:
+ opflow - the opflow application object
- X   - the current iterate


  Output Parameters:
. Ge  - vector of equality constraints
*/
PetscErrorCode OPFLOWEqualityConstraints(OPFLOW opflow,Vec X,Vec Ge)
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

    if (!bus->isghost) {
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

/*
  OPFLOWInequalityConstraints - Evalulates the inequality constraints for the optimal power flow

  Input Parameters:
+ opflow - the opflow application object
- X   - the current iterate

  Output Parameters:
. Gi  - vector of inequality constraints
*/
PetscErrorCode OPFLOWInequalityConstraints(OPFLOW opflow,Vec X,Vec Gi)
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


/* Coupling constraints */
PetscErrorCode SCOPFLOWAddCouplingConstraints(OPFLOW opflow,Vec X0,Vec X1,Vec Gi)
{
  PetscErrorCode ierr;
  PetscScalar    *x0,*x1,*gi;
  PS             ps = opflow->ps; /* PS for the scenario */
  PetscInt       gloc; /* starting location for coupling constraints in G vector */
  PetscInt       xloc;
  PetscInt       i;
  PSBUS          bus;

  PetscFunctionBegin;

  gloc = opflow->Nconineq;

  ierr = VecGetArray(X0,&x0);CHKERRQ(ierr);
  ierr = VecGetArray(X1,&x1);CHKERRQ(ierr);
  ierr = VecGetArray(Gi,&gi);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    PetscInt k;

    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);

    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;
      xloc += 2;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      gi[gloc] = x1[xloc] - x0[xloc];

      gloc += 1;
    }
  }

  ierr = VecRestoreArray(X0,&x0);CHKERRQ(ierr);
  ierr = VecRestoreArray(X1,&x1);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gi,&gi);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

/*
  OPFLOWInequalityConstraintsJacobian - Jacobian evaluation routine for the inequality constraints

  Input Parameters:
+ opflow - the opflow solver object
. X   - the current iterate
- scenario - scenario number


  Output Parameters:
. Ji - Jacobian for inequality constraints

*/
PetscErrorCode OPFLOWInequalityConstraintsJacobian(OPFLOW opflow, Vec X, Mat Ji,PetscInt scenario)
{
  PetscErrorCode ierr;
  PetscInt       ctr=0,i;
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
    if(!line->status) gloc += 2;

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
    val[ctr]   = dSt2_dthetat;
    val[ctr+1] = dSt2_dVmt;
    val[ctr+2] = dSt2_dthetaf;
    val[ctr+3] = dSt2_dVmf;
    ierr = MatSetValues(Ji,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

    gloc += 2;
  }

  if(scenario != 0) {
    PSBUS bus;
    PetscInt k,xloc;
    for(i=0; i < ps->nbus; i++) {
      bus = &ps->bus[i];

      ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);
      for(k=0; k < bus->ngen; k++) {
	xloc += 2;
	val[0] = 1;
	row[0] = gloc;
	col[0] = xloc;
	ierr = MatSetValues(Ji,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
	gloc += 1;
      }
    }
  }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(Ji,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Ji,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWEqualityConstraintsJacobian - Sets the nonzero values for the
              equality constraints Jacobian

  Input Parameters:
+ OPFLOW - opflow solver object
- X   - the current iterate

  Output Parameters:
. Je - Jacobian of equality constraints
*/
PetscErrorCode OPFLOWEqualityConstraintsJacobian(OPFLOW opflow, Vec X,Mat Je)
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

    if (!bus->isghost) {
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



/*
  OPFLOWGetJacobianNonzeros - Gets the number of nonzeros in the constraint jacobian matrix

  Input Paramereters:
. opflow - the optimal power flow application object

  Output Parameters:
+ e_nnz  - number of nonzeros for the equality constraint
- i_nnz - number of nonzeros for the inequality constraint

*/
PetscErrorCode OPFLOWGetJacobianNonzeros(OPFLOW opflow,PetscInt *e_nnz,PetscInt *i_nnz)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscInt       i,k;
  PSBUS          bus;
  PSLINE         line;
  PetscInt       nconnlines;
  const PSLINE   *connlines;

  PetscFunctionBegin;

  *e_nnz = *i_nnz = 0;

  for(i=0; i < ps->Nbranch; i++) *i_nnz += 8;

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    if(bus->ide == ISOLATED_BUS) {
      *e_nnz += 2;
      continue;
    }	
    *e_nnz += 2 + 2*bus->ngen;
    ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);
    for(k=0; k < nconnlines;k++) {
      line = connlines[k];
      if(!line->status) continue;
      *e_nnz += 8;
    }
  }

  PetscFunctionReturn(0);
}


int str_eval_g(double* x0, double* x1, double* eq_g, double* inq_g,
	       CallBackDataPtr cbd) 
{
  PetscErrorCode ierr;
  int row = cbd->row_node_id;
  SCOPFLOW scopflow=(SCOPFLOW)cbd->prob;
  OPFLOW   opflow=scopflow->opflows[row];
  double   *x;

  if(row == 0) x = x0;
  else x = x1;

  /* Equality constraints */
  ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
  ierr = VecPlaceArray(opflow->Ge,eq_g);CHKERRQ(ierr);
  ierr = OPFLOWEqualityConstraints(opflow,opflow->X,opflow->Ge);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->Ge);CHKERRQ(ierr);

  /* Inequality Constraints */
  ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
  ierr = VecPlaceArray(opflow->Gi,inq_g);CHKERRQ(ierr);
  ierr = OPFLOWInequalityConstraints(opflow,opflow->X,opflow->Gi);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->Gi);CHKERRQ(ierr);

  if(row != 0) {
    Vec  X0=scopflow->opflows[0]->X;
    Vec  X1=opflow->X;
    ierr = VecPlaceArray(X1,x1);CHKERRQ(ierr);
    ierr = VecPlaceArray(X0,x0);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Gi,inq_g);CHKERRQ(ierr);
    /* Add coupling constraints */
    ierr = SCOPFLOWAddCouplingConstraints(opflow,X0,X1,opflow->Gi);CHKERRQ(ierr);
    
    ierr = VecResetArray(X1);CHKERRQ(ierr);
    ierr = VecResetArray(X0);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Gi);CHKERRQ(ierr);
  }
  return 1;
}


int str_eval_jac_g(double* x0, double* x1, int* e_nz, double* e_elts,
		int* e_rowidx, int* e_colptr, int* i_nz, double* i_elts, int* i_rowidx,
		int* i_colptr, CallBackDataPtr cbd) {
  int row = cbd->row_node_id;
  int col = cbd->col_node_id;
  PetscErrorCode ierr;
  SCOPFLOW scopflow=(SCOPFLOW)cbd->prob;
  /*  PetscInt rank=scopflow->comm->rank; */
  OPFLOW   opflow;
  Mat_SeqAIJ *aij;
  PetscInt    nrow,ncol;

  if(e_elts==NULL && i_elts == NULL) {
    /* Number of non-zeros in equality and inequality constraint Jacobian */
    opflow = scopflow->opflows[row];
    if(row == col) {
      ierr = VecPlaceArray(opflow->X,x1);CHKERRQ(ierr);
      /* Equality constraints Jacobian */
      ierr = OPFLOWEqualityConstraintsJacobian(opflow,opflow->X,opflow->Jac_Ge);CHKERRQ(ierr);
      ierr = MatTranspose(opflow->Jac_Ge,MAT_INITIAL_MATRIX,&opflow->Jac_GeT);CHKERRQ(ierr);
      aij = (Mat_SeqAIJ*)opflow->Jac_GeT->data;
      *e_nz = aij->nz;

      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);

      ierr = VecPlaceArray(opflow->X,x1);CHKERRQ(ierr);

      /* Inequality Constraints Jacobian */
      ierr = OPFLOWInequalityConstraintsJacobian(opflow,opflow->X,opflow->Jac_Gi,row);CHKERRQ(ierr);
      ierr = MatTranspose(opflow->Jac_Gi,MAT_INITIAL_MATRIX,&opflow->Jac_GiT);CHKERRQ(ierr);
      aij = (Mat_SeqAIJ*)opflow->Jac_GiT->data;
      *i_nz = aij->nz;

      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);

    } else {
      if(col == 0) {
	*e_nz = 0;
	aij = (Mat_SeqAIJ*)scopflow->JcoupT->data;
	*i_nz = aij->nz;

      }
    }
  } else {
    opflow = scopflow->opflows[row];
    if(row == col) {
      ierr = VecPlaceArray(opflow->X,x0);CHKERRQ(ierr);
      /* Equality constraints Jacobian */
      ierr = OPFLOWEqualityConstraintsJacobian(opflow,opflow->X,opflow->Jac_Ge);CHKERRQ(ierr);
      ierr = MatTranspose(opflow->Jac_Ge,MAT_REUSE_MATRIX,&opflow->Jac_GeT);CHKERRQ(ierr);
      ierr = MatGetSize(opflow->Jac_GeT,&nrow,&ncol);CHKERRQ(ierr);
      aij = (Mat_SeqAIJ*)opflow->Jac_GeT->data;
      ierr = PetscMemcpy(e_rowidx,aij->j,aij->nz*sizeof(PetscInt));CHKERRQ(ierr);
      ierr = PetscMemcpy(e_colptr,aij->i,(nrow+1)*sizeof(PetscInt));CHKERRQ(ierr);
      ierr = PetscMemcpy(e_elts,aij->a,aij->nz*sizeof(PetscScalar));CHKERRQ(ierr);

      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);

      ierr = VecPlaceArray(opflow->X,x0);CHKERRQ(ierr);
      /* Inequality Constraints Jacobian */
      ierr = OPFLOWInequalityConstraintsJacobian(opflow,opflow->X,opflow->Jac_Gi,row);CHKERRQ(ierr);
      ierr = MatTranspose(opflow->Jac_Gi,MAT_REUSE_MATRIX,&opflow->Jac_GiT);CHKERRQ(ierr);
      aij = (Mat_SeqAIJ*)opflow->Jac_GiT->data;
      ierr = PetscMemcpy(i_rowidx,aij->j,aij->nz*sizeof(PetscInt));CHKERRQ(ierr);
      ierr = PetscMemcpy(i_colptr,aij->i,(nrow+1)*sizeof(PetscInt));CHKERRQ(ierr);
      ierr = PetscMemcpy(i_elts,aij->a,aij->nz*sizeof(PetscScalar));CHKERRQ(ierr);

      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);


    } else {
      if(col == 0) {
	ierr = MatGetSize(scopflow->JcoupT,&nrow,&ncol);CHKERRQ(ierr);
	aij = (Mat_SeqAIJ*)scopflow->JcoupT->data;
	ierr = PetscMemcpy(i_rowidx,aij->j,aij->nz*sizeof(PetscInt));CHKERRQ(ierr);
	ierr = PetscMemcpy(i_colptr,aij->i,(nrow+1)*sizeof(PetscInt));CHKERRQ(ierr);
	ierr = PetscMemcpy(i_elts,aij->a,aij->nz*sizeof(PetscScalar));CHKERRQ(ierr);
      }
    }
  }

  return 1;
}


