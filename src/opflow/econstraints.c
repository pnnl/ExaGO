#include <private/opflowimpl.h>

/****************************************************
  Designed to be a side file dedicated to
  the functions dealing with the equality constraints
  Order to file:
    OPFLOWEqualityConstraintsFunction
    OPFLOWCreateEqualityConstraintsJacobian
    OPFLOWEqualityConstraintsJacobianFunction
*****************************************************/

/*
  OPFLOWEqualityConstraintsFunction - Evalulates the equality constraints for the optimal power flow

  Input Parameters:
+ Tao - the Tao nlp solver object
. X   - the current iterate
- ctx - application data set with TaoSetEqualityConstraintsRoutine

  Output Parameters:
. Ge  - vector of equality constraints
*/
PetscErrorCode OPFLOWEqualityConstraintsFunction(Tao nlp,Vec X,Vec Ge,void* ctx)
{
  PetscErrorCode ierr;
  PetscInt       i,k,nconnlines;
  PetscInt       gloc,row[2];
  PetscInt       xloc,xlocf,xloct;
  PetscMPIInt    rank,size;
  PetscScalar    val[2];
  PetscScalar    Pg,Qg;
  PetscScalar    Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  PetscScalar    Vmf,Vmt,thetaf,thetat,thetaft,thetatf;
  PetscScalar    Pf,Qf,Pt,Qt;
  PetscScalar    theta,Vm;
  OPFLOW         opflow=(OPFLOW)ctx;
  PS             ps=opflow->ps;
  MPI_Comm       comm=opflow->comm->type;
  Vec            localX = opflow->localX;
  PSLOAD         load;
  PSLINE         line;
  PSBUS          bus,busf,bust;
  const PSBUS    *connbuses;
  const PSLINE   *connlines;
  const PetscScalar *x;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
  ierr = VecSet(Ge,0.0);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);

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

  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWCreateEqualityConstraintsJacobian - Returns a distributed matrix of appropriate size that can be used as the Jacobian for the equality constraints

  Input Paramereters:
. opflow - the optimal power flow application object

  Output Parameters:
. mat - the jacobian of the equality constraints

  Notes:
  OPFLOWSetUp() must be called before calling this routine.
*/
PetscErrorCode OPFLOWCreateEqualityConstraintsJacobian(OPFLOW opflow,Mat *mat)
{
  PetscErrorCode ierr;
  PetscInt       nconeq=opflow->nconeq; /* Local number of equality constraints */
  PetscInt       nvar=opflow->nvar;     /* Local number of variables */
  PetscInt       i,k,nconnlines;
  PetscInt       loc,locf,loct,gidx;
  PetscInt       row[2],col[2];
  PetscInt       *dnnz,*onnz;
  PetscMPIInt    rank,size;
  PetscScalar    val[4];
  PS             ps=opflow->ps;
  MPI_Comm       comm=opflow->comm->type;
  Mat            jac;
  PSLINE         line;
  PSBUS          bus,busf,bust;
  DM             networkdm=ps->networkdm;
  const PSBUS    *connbuses;
  const PSLINE   *connlines;
  Vec            Vdnz,Vonz;
  const PetscScalar *varr;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);

  ierr = MatCreate(opflow->comm->type,&jac);CHKERRQ(ierr);
  ierr = MatSetSizes(jac,nconeq,nvar,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = MatSetType(jac,MATAIJ);CHKERRQ(ierr);
#if 1
  /* Set up preallocation */
  ierr = PetscCalloc2(nconeq,&dnnz,nconeq,&onnz);CHKERRQ(ierr);

  ierr = VecCreate(comm,&Vdnz);CHKERRQ(ierr);
  ierr = VecSetSizes(Vdnz,nconeq,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vdnz);CHKERRQ(ierr);
  ierr = VecSet(Vdnz,0.0);CHKERRQ(ierr);

  ierr = VecCreate(comm,&Vonz);CHKERRQ(ierr);
  ierr = VecSetSizes(Vonz,nconeq,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vonz);CHKERRQ(ierr);
  ierr = VecSet(Vonz,0.0);CHKERRQ(ierr);


  for (i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = DMNetworkGetVertexLocalToGlobalOrdering(networkdm,i,&gidx);CHKERRQ(ierr);
    gidx *= 2;
    row[0] = gidx; row[1] = row[0] + 1;

    if (!bus->isghost) {
      val[0] = val[1] = (PetscScalar)(2 + bus->ngen);
      ierr = VecSetValues(Vdnz,2,row,val,ADD_VALUES);CHKERRQ(ierr);
    }

    ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);
    for (k=0; k < nconnlines; k++) {
      line = connlines[k];
      if (!line->status) continue;

      /* Get the connected buses to this line */
      ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];
      val[0] = 2.0; val[1] = 2.0;

      if (!bus->isghost) {
        if (bus == busf) {
          if (!bust->isghost) {
            ierr = VecSetValues(Vdnz,2,row,val,ADD_VALUES);CHKERRQ(ierr);
          } else {
            ierr = VecSetValues(Vonz,2,row,val,ADD_VALUES);CHKERRQ(ierr);
          }
        } else {
          if (!busf->isghost) {
            ierr = VecSetValues(Vdnz,2,row,val,ADD_VALUES);CHKERRQ(ierr);
          } else {
            ierr = VecSetValues(Vonz,2,row,val,ADD_VALUES);CHKERRQ(ierr);
          }
        }
      } else { /* bus->isghost */
        ierr = VecSetValues(Vonz,2,row,val,ADD_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = VecAssemblyBegin(Vdnz);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(Vdnz);CHKERRQ(ierr);
  //ierr = VecView(Vdnz,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  ierr = VecAssemblyBegin(Vonz);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(Vonz);CHKERRQ(ierr);
  //ierr = VecView(Vonz,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  ierr = VecGetArrayRead(Vdnz,&varr);CHKERRQ(ierr);
  for (i=0; i<nconeq; i++) dnnz[i] = (PetscInt)varr[i];
  ierr = VecRestoreArrayRead(Vdnz,&varr);CHKERRQ(ierr);
  ierr = VecDestroy(&Vdnz);CHKERRQ(ierr);

  ierr = VecGetArrayRead(Vonz,&varr);CHKERRQ(ierr);
  for (i=0; i<nconeq; i++) onnz[i] = (PetscInt)varr[i];
  ierr = VecRestoreArrayRead(Vonz,&varr);CHKERRQ(ierr);
  ierr = VecDestroy(&Vonz);CHKERRQ(ierr);

  ierr = MatSeqAIJSetPreallocation(jac,0,dnnz);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(jac,0,dnnz,0,onnz);CHKERRQ(ierr);
  ierr = PetscFree2(dnnz,onnz);CHKERRQ(ierr);
  CHKMEMQ;
  //ierr = MPI_Barrier(comm);CHKERRQ(ierr);

  /* TODO: 'mpiexec -n 5 ./OPFLOW -petscpartitioner_type simple' throws an error;
     Thus enable 'NEW_NONZERO_ALLOCATION' below until the bug is fixed. */
  ierr = MatSetOption(jac,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);CHKERRQ(ierr);
#endif

  for (i=0; i< 4; i++) val[i] = 0.0;
  for (i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = DMNetworkGetVertexLocalToGlobalOrdering(networkdm,i,&gidx);CHKERRQ(ierr);
    gidx = 2*gidx;
    row[0] = gidx; row[1] = row[0] + 1;

    if (!bus->isghost) {
      ierr = PSBUSGetVariableGlobalLocation(bus,&loc);CHKERRQ(ierr);
      col[0] = loc;  col[1] = loc + 1;
      ierr = MatSetValues(jac,2,row,2,col,val,INSERT_VALUES);CHKERRQ(ierr);

      for (k=0; k < bus->ngen; k++) { /* one entry per row! */
        /* Skip over OFF generators OR use status*x[i] ? */
        loc += 2;
        col[0] = loc;
        ierr = MatSetValues(jac,1,row,1,col,val,INSERT_VALUES);CHKERRQ(ierr);
        col[1] = loc+1;
        ierr = MatSetValues(jac,1,row+1,1,col+1,val,INSERT_VALUES);CHKERRQ(ierr);
      }
    }

    ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);
    for (k=0; k < nconnlines; k++) {
      line = connlines[k];
      ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];

      ierr = PSBUSGetVariableGlobalLocation(busf,&locf);CHKERRQ(ierr);
      ierr = PSBUSGetVariableGlobalLocation(bust,&loct);CHKERRQ(ierr);
      if (bus == busf) {
        col[0] = loct;
      } else col[0] = locf;
      col[1] = col[0] + 1;
      ierr = MatSetValues(jac,2,row,2,col,val,INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  CHKMEMQ;

  *mat = jac;
  //ierr = PetscPrintf(comm,"Je structure:\n");CHKERRQ(ierr);
  //ierr = MatView(jac,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*
  OPFLOWEqualityConstraintsJacobianFunction - Sets the nonzero values for the
              equality constraints Jacobian

  Input Parameters:
+ nlp - Tao nlp solver object
. X   - the current iterate
- ctx - application data set with OPFLOWEqualityConstraintsJacobianRoutine


  Output Parameters:
. Je - Jacobian of equality constraints
*/
PetscErrorCode OPFLOWEqualityConstraintsJacobianFunction(Tao nlp, Vec X,Mat Je, Mat Je_pre, void* ctx)
{
  PetscErrorCode ierr;
  PetscInt       i,k,row[2],col[4],genctr,gidx;
  PetscInt       nconnlines,locglob,loc,locglobf,locglobt,locf,loct;
  PetscScalar    Vm,val[8],Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  PetscScalar    thetaf,thetat,Vmf,Vmt,thetaft,thetatf;
  OPFLOW         opflow=(OPFLOW)ctx;
  Vec            localX=opflow->localX;
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

  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);
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
  ierr = VecRestoreArrayRead(localX,&xarr);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(Je,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Je,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
