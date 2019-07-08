#include <private/psimpl.h>
#include <private/opflowimpl.h>

/*********************************************
  The optimization problem for Tao should be in
  the following form

  min f(x)
  s.t.
    g(x) = 0
    h(x) >= 0
    xl <= x <= xu

************************************/

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
  OPFLOW         opflow=(OPFLOW)ctx;
  PS             ps=opflow->ps;
  Vec            localX;
  const PetscScalar *xarr;
  PetscBool      ghostbus;
  PetscInt       i,k;
  PetscInt       rowctr=0;
  PetscInt       rstart;

  PetscFunctionBegin;
  ierr = MatZeroEntries(Je);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(Je,&rstart,NULL);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    PetscScalar Vm;
    PetscInt    loc;
    PetscInt    locglob;
    PSBUS       bus;
    PetscInt    row[2],col[4];
    PetscScalar val[4];
    PetscInt    genctr;

    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bus,&locglob);CHKERRQ(ierr);
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);

    Vm = xarr[loc+1];

    if(!ghostbus) {
      row[0] = rstart + rowctr;
      row[1] = rstart + rowctr + 1;

      col[0] = locglob; col[1] = locglob+1;
      /* Isolated and reference bus */
      if(bus->ide == ISOLATED_BUS) {
	val[0] = val[3] = 1.0;
	val[1] = val[2] = 0.0;
	ierr = MatSetValues(Je,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
	rowctr += 2;
	continue;
      }

      /* Shunt injections */
      val[0] = 0.0; val[1] = 2*Vm*bus->gl;
      val[2] = 0.0; val[3]= -2*Vm*bus->bl; /* Partial derivative for shunt contribution */
      ierr = MatSetValues(Je,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);

      genctr = 0;
      for(k=0; k < bus->ngen; k++) {
	col[0] = locglob + 2 + genctr;
	val[0] = -1;
	ierr = MatSetValues(Je,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	col[0] = locglob + 2 + genctr+1;
	val[0] = -1;
	ierr = MatSetValues(Je,1,row+1,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	genctr += 2;
      }
    }

    /* Partial derivatives of network equations */
    PetscInt nconnlines;
    const PSLINE *connlines;
    PSLINE line;

    /* Get the lines supporting the bus */
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

      const PSBUS *connbuses;
      PSBUS busf,bust;
      PetscInt locglobf,locglobt,locf,loct;
      PetscScalar thetaf,thetat,Vmf,Vmt,thetaft,thetatf;

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

      if(bus == busf) {
	row[0] = rstart + rowctr;
	col[0] = locglobf; col[1] = locglobf+1; col[2] = locglobt; col[3] = locglobt+1;
	val[0] = Vmf*Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft));
	val[1] = 2*Gff*Vmf + Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
	val[2] = Vmf*Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft));
	val[3] = Vmf*(Gft*cos(thetaft) + Bft*sin(thetaft));
	ierr = MatSetValues(Je,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

	row[0] = rstart + rowctr + 1;
	val[0] = Vmf*Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft));
	val[1] = -2*Bff*Vmf + Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
	val[2] = Vmf*Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft));
	val[3] = Vmf*(-Bft*cos(thetaft) + Gft*sin(thetaft));
	ierr = MatSetValues(Je,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
      } else {
	row[0] = rstart + rowctr;
	col[0] = locglobt; col[1] = locglobt+1; col[2] = locglobf; col[3] = locglobf+1;
	val[0] = Vmt*Vmf*(-Gtf*sin(thetatf) + Btf*cos(thetatf));
	val[1] = 2*Gtt*Vmt + Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf));
	val[2] = Vmt*Vmf*(Gtf*sin(thetatf) - Btf*cos(thetatf));
	val[3] = Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf));
	ierr = MatSetValues(Je,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

	row[0] = rstart + rowctr + 1;
	val[0] = Vmt*Vmf*(Btf*sin(thetatf) + Gtf*cos(thetatf));
	val[1] = -2*Btt*Vmt + Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
	val[2] = Vmt*Vmf*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
	val[3] = Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
	ierr = MatSetValues(Je,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
      }
    }
    rowctr += 2;
  }

  ierr = VecRestoreArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(Je,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Je,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  //  ierr = MatView(Je,0);
  //  exit(1);
  PetscFunctionReturn(0);
}

/*
  OPFLOWCreate - Creates an optimal power flow application object

  Input Parameters
. mpicomm - The MPI communicator

  Output Parameters
. opflowout - The optimal power flow application object
*/
PetscErrorCode OPFLOWCreate(MPI_Comm mpicomm, OPFLOW *opflowout)
{
  PetscErrorCode ierr;
  OPFLOW         opflow;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&opflow);CHKERRQ(ierr);

  ierr = COMMCreate(mpicomm,&opflow->comm);CHKERRQ(ierr);

  //if(opflow->comm->size > 1) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Optimal Power Flow not supported in parallel");

  ierr = PSCreate(mpicomm,&opflow->ps);CHKERRQ(ierr);

  /* Set the application with the PS object */
  ierr = PSSetApplication(opflow->ps,APP_ACOPF);CHKERRQ(ierr);

  opflow->Nconeq    = -1;
  opflow->Nconineq  = -1;
  opflow->Ncon      = -1;
  opflow->nlp       = NULL;

  /* Create the optimization solver */
  ierr = TaoCreate(mpicomm,&opflow->nlp);CHKERRQ(ierr);
  ierr = TaoSetType(opflow->nlp,TAOIPM);CHKERRQ(ierr);
  /* Set the prefix for tao.. All TAO runtime options will need to have the prefix "-opflow_" */
  ierr = TaoSetOptionsPrefix(opflow->nlp,"opflow_");CHKERRQ(ierr);

  opflow->setupcalled = PETSC_FALSE;

  *opflowout = opflow;
  PetscFunctionReturn(0);
}

/*
  OPFLOWDestroy - Destroys the optimal power flow application object

  Input Parameter
. opflow - The OPFLOW object to destroy
*/
PetscErrorCode OPFLOWDestroy(OPFLOW *opflow)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = COMMDestroy(&(*opflow)->comm);CHKERRQ(ierr);

  /* Solution vector */
  ierr = VecDestroy(&(*opflow)->X);CHKERRQ(ierr);

  /* Lower and upper bounds on X */
  ierr = VecDestroy(&(*opflow)->Xl);CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Xu);CHKERRQ(ierr);

  /* Constraints vector */
  ierr = VecDestroy(&(*opflow)->Ge);CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Gi);CHKERRQ(ierr);

  /* Jacobian of constraints */
  ierr = MatDestroy(&(*opflow)->Jac_Ge);CHKERRQ(ierr);
  ierr = MatDestroy(&(*opflow)->Jac_Gi);CHKERRQ(ierr);

  ierr = PSDestroy(&(*opflow)->ps);CHKERRQ(ierr);

  ierr = TaoDestroy(&(*opflow)->nlp);CHKERRQ(ierr);

  ierr = PetscFree(*opflow);CHKERRQ(ierr);
  *opflow = 0;

  PetscFunctionReturn(0);
}

/*
  OPFLOWReadMatPowerData - Reads the network data given in MATPOWER data format

  Input Parameter
+  opflow - The OPFLOW object
-  netfile - The name of the network file

*/
PetscErrorCode OPFLOWReadMatPowerData(OPFLOW opflow,const char netfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* Read MatPower data file and populate the PS data structure */
  ierr = PSReadMatPowerData(opflow->ps,netfile);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWCreateGlobalVector - Returns a global vector of the appropriate size
  and distribution conforming to the distribution of the PS object.

  Input Paramereters:
. opflow - the optimal power flow application object

  Output Parameters:
. vec - the global vector

  Notes:
  OPFLOWSetUp() must be called before calling this routine.

  If a vector of size X is needed by the OPFLOW application then this routine can be called.
*/
PetscErrorCode OPFLOWCreateGlobalVector(OPFLOW opflow,Vec *vec)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if(!opflow->setupcalled) SETERRQ(opflow->comm->type,0,"OPFLOWSetUp() must be called before calling PFLOWCreateGlobalVector");
  ierr = PSCreateGlobalVector(opflow->ps,vec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWCreateInequalityConstraintsJacobian - Returns a distributed matrix of appropriate size that can be used as the Jacobian for the inequality constraints


  Input Paramereters:
. opflow - the optimal power flow application object

  Output Parameters:
. mat - the jacobian of the inequality constraints

  Notes:
  OPFLOWSetUp() must be called before calling this routine.
*/
PetscErrorCode OPFLOWCreateInequalityConstraintsJacobian(OPFLOW opflow,Mat *mat)
{
  PetscErrorCode ierr;
  Mat            jac;
  PetscInt       nconineq=opflow->nconineq,Nconineq=opflow->Nconineq; /* Number of constraints */
  PetscInt       nvar=opflow->nvar,Nvar=opflow->Nvar; /* Number of variables */
  PS             ps=opflow->ps;
  PetscInt       i;
  PSLINE         line;
  PetscInt       mloc=0;
  const PSBUS    *connbuses;
  PSBUS          busf,bust;
  PetscInt       locf,loct;
  PetscInt       *nnz;
  PetscInt       row[2],col[4];
  PetscScalar    val[8];
  MPI_Comm       comm=opflow->comm->type;
  PetscMPIInt    rank;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  printf("[%d] Ji: Nconineq %d %d, Nvar %d %d\n",rank,nconineq,Nconineq,nvar,Nvar);

  ierr = MatCreate(opflow->comm->type,&jac);CHKERRQ(ierr);
  ierr = MatSetSizes(jac,nconineq,nvar,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = MatSetType(jac,MATAIJ);CHKERRQ(ierr);
  ierr = MatSetUp(jac);CHKERRQ(ierr);
#if 0
  /* Set up preallocation */
  ierr = PetscCalloc1(Nconineq,&nnz);CHKERRQ(ierr);

  for(i=0; i < ps->Nbranch; i++) {
    line = &ps->line[i];

    nnz[mloc] += 4;
    nnz[mloc+1] += 4;
    nnz[mloc+2] += 4;
    nnz[mloc+3] += 4;

    mloc += 4;
  }

  ierr = MatSeqAIJSetPreallocation(jac,0,nnz);CHKERRQ(ierr);
  ierr = PetscFree(nnz);CHKERRQ(ierr);
#endif
  PetscInt rstart,rend;
  ierr = MatGetOwnershipRange(jac,&rstart,&rend);CHKERRQ(ierr);
  //printf("[%d] Ji rstart/end: %d %d;\n",rank,rstart,rend);
  ierr = MPI_Barrier(comm);CHKERRQ(ierr);

  for (i=0; i<8; i++) val[i] = 0.0;

  /* Set up nonzero structure. 0 is inserted in the non-zero location */
  mloc = rstart;
  for(i=0; i < ps->nbranch; i++) {
    line = &ps->line[i];

    ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
    busf = connbuses[0];
    bust = connbuses[1];

    ierr = PSBUSGetVariableGlobalLocation(busf,&locf);CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bust,&loct);CHKERRQ(ierr);

    /* From bus Sf lower and upper bound*/
    row[0] = mloc, row[1] = mloc+1;
    col[0] = locf; col[1] = locf+1; col[2] = loct; col[3] = loct+1;
    ierr = MatSetValues(jac,2,row,4,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /* To bus St lower and upper bound */
    row[0] = mloc+2, row[1] = mloc+3;
    col[0] = loct; col[1] = loct+1; col[2] = locf; col[3] = locf+1;
    ierr = MatSetValues(jac,2,row,4,col,val,INSERT_VALUES);CHKERRQ(ierr);

    mloc += 4;
  }

  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  *mat = jac;
  //ierr = PetscPrintf(comm,"Ji structure:\n");CHKERRQ(ierr);
  //ierr = MatView(jac,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWInequalityConstraintsJacobianFunction - Jacobian evaluation routine for the inequality constraints

  Input Parameters:
+ nlp - the TAO nlp solver object
. X   - the current iterate
- ctx - the application data set with TaoInequalityConstraintsJacobianRoutine

  Output Parameters:
+ Ji - Jacobian for inequality constraints
- Ji_pre - Preconditioner (same as Jacobian, no change needed) for
           inequality constraints

*/
PetscErrorCode OPFLOWInequalityConstraintsJacobianFunction(Tao nlp, Vec X, Mat Ji, Mat Ji_pre, void* ctx)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=(OPFLOW)ctx;
  PS             ps=opflow->ps;
  PetscInt       ctr=0;
  PetscInt       i;
  PSLINE         line;
  const PSBUS    *connbuses;
  PetscInt       gloc=0,xlocf,xloct;
  PSBUS          busf,bust;
  PetscInt       row[2],col[4];
  PetscScalar    val[4];
  const PetscScalar *x;

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);

  for(i=0; i < ps->Nbranch; i++) {
    line = &ps->line[i];

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

    PetscScalar Pf,Qf,Pt,Qt;

    Pf = Gff*Vmf*Vmf  + Vmf*Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
    Qf = -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));

    Pt = Gtt*Vmt*Vmt  + Vmt*Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf));
    Qt = -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));

    PetscScalar dSf2_dPf, dSf2_dQf, dSt2_dPt, dSt2_dQt;

    dSf2_dPf = 2*Pf;
    dSf2_dQf = 2*Qf;
    dSt2_dPt = 2*Pt;
    dSt2_dQt = 2*Qt;

    PetscScalar dPf_dthetaf,dPf_dVmf,dPf_dthetat,dPf_dVmt;
    PetscScalar dQf_dthetaf,dQf_dVmf,dQf_dthetat,dQf_dVmt;
    PetscScalar dPt_dthetaf,dPt_dVmf,dPt_dthetat,dPt_dVmt;
    PetscScalar dQt_dthetaf,dQt_dVmf,dQt_dthetat,dQt_dVmt;

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

    PetscScalar dSf2_dthetaf,dSf2_dVmf,dSf2_dthetat,dSf2_dVmt;

    dSf2_dthetaf = dSf2_dPf*dPf_dthetaf + dSf2_dQf*dQf_dthetaf;
    dSf2_dthetat = dSf2_dPf*dPf_dthetat + dSf2_dQf*dQf_dthetat;
    dSf2_dVmf    = dSf2_dPf*dPf_dVmf    + dSf2_dQf*dQf_dVmf;
    dSf2_dVmt    = dSf2_dPf*dPf_dVmt    + dSf2_dQf*dQf_dVmt;

    /* g[gloc] */
    row[0] = gloc;
    col[0] = xlocf; col[1] = xlocf+1; col[2] = xloct; col[3] = xloct+1;
    val[0] = dSf2_dthetaf;
    val[1] = dSf2_dVmf;
    val[2] = dSf2_dthetat;
    val[3] = dSf2_dVmt;

    ierr = MatSetValues(Ji,1,row,4,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /* g[gloc+1] */
    row[0] = gloc+1;
    col[0] = xlocf; col[1] = xlocf+1; col[2] = xloct; col[3] = xloct+1;
    val[0] = -dSf2_dthetaf;
    val[1] = -dSf2_dVmf;
    val[2] = -dSf2_dthetat;
    val[3] = -dSf2_dVmt;

    ierr = MatSetValues(Ji,1,row,4,col,val,INSERT_VALUES);CHKERRQ(ierr);

    PetscScalar dSt2_dthetaf,dSt2_dVmf,dSt2_dthetat,dSt2_dVmt;

    dSt2_dthetaf = dSt2_dPt*dPt_dthetaf + dSt2_dQt*dQt_dthetaf;
    dSt2_dthetat = dSt2_dPt*dPt_dthetat + dSt2_dQt*dQt_dthetat;
    dSt2_dVmf    = dSt2_dPt*dPt_dVmf    + dSt2_dQt*dQt_dVmf;
    dSt2_dVmt    = dSt2_dPt*dPt_dVmt    + dSt2_dQt*dQt_dVmt;

    /* g[gloc+2] */
    row[0] = gloc+2;
    col[0] = xloct; col[1] = xloct+1; col[2] = xlocf; col[3] = xlocf+1;

    val[ctr]   = dSt2_dthetat;
    val[ctr+1] = dSt2_dVmt;
    val[ctr+2] = dSt2_dthetaf;
    val[ctr+3] = dSt2_dVmf;

    ierr = MatSetValues(Ji,1,row,4,col,val,INSERT_VALUES);CHKERRQ(ierr);

    /* g[gloc+3] */
    row[0] = gloc+3;
    col[0] = xloct; col[1] = xloct+1; col[2] = xlocf; col[3] = xlocf+1;

    val[ctr]   = -dSt2_dthetat;
    val[ctr+1] = -dSt2_dVmt;
    val[ctr+2] = -dSt2_dthetaf;
    val[ctr+3] = -dSt2_dVmf;

    ierr = MatSetValues(Ji,1,row,4,col,val,INSERT_VALUES);CHKERRQ(ierr);

    gloc += 4;
  }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(Ji,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Ji,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

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
  Mat            jac;
  PetscInt       nconeq=opflow->nconeq,Nconeq=opflow->Nconeq; /* Local and global number of equality constraints */
  PetscInt       nvar=opflow->nvar,Nvar=opflow->Nvar; /* Number of variables */
  PS             ps=opflow->ps;
  PetscInt       i,k;
  PSLINE         line;
  PSBUS          bus;
  PetscInt       mloc;
  const PSBUS    *connbuses;
  PetscInt       nconnlines;
  const PSLINE   *connlines;
  PSBUS          busf,bust;
  PetscInt       loc,locf,loct;
  PetscInt       *nnz;
  PetscInt       row[2],col[4];
  PetscScalar    val[4];
  MPI_Comm       comm=opflow->comm->type;
  PetscMPIInt    rank;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  printf("[%d] Je: nconeq %d, Nconeq %d,Nvar %d\n",rank,nconeq,Nconeq,Nvar);
  ierr = MatCreate(opflow->comm->type,&jac);CHKERRQ(ierr);
  ierr = MatSetSizes(jac,nconeq,nvar,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = MatSetType(jac,MATAIJ);CHKERRQ(ierr);
  ierr = MatSetUp(jac);CHKERRQ(ierr);
#if 0
  /* Set up preallocation */
  ierr = PetscCalloc1(nconeq,&nnz);CHKERRQ(ierr);

  mloc = 0;
  for(i=0; i < ps->nbus; i++) {
    if (bus->isghost) continue;
    bus = &ps->bus[i];
    ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);

    nnz[mloc]   += 2*nconnlines + 2 + 2*bus->ngen;
    nnz[mloc+1] += 2*nconnlines + 2 + 2*bus->ngen;
    mloc += 2;
  }

  ierr = MatSeqAIJSetPreallocation(jac,0,nnz);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(jac,0,nnz,0,nnz);CHKERRQ(ierr);
  ierr = PetscFree(nnz);CHKERRQ(ierr);
#endif
  PetscInt rstart,rend;
  ierr = MatGetOwnershipRange(jac,&rstart,&rend);CHKERRQ(ierr);
  printf("[%d] Je rstart/end: %d %d;\n",rank,rstart,rend);

  mloc = rstart; /* Should use MatGetOwnershipRange here for parallel */
  val[0] = val[1] = val[2] = val[3] = 0.0;
  for(i=0; i < ps->nbus; i++) {
    //if (bus->isghost) continue; /* skip ghost buses */
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableGlobalLocation(bus,&loc);CHKERRQ(ierr);
    row[0] = mloc;
    row[1] = mloc+1;
    col[0] = loc; col[1] = loc+1;
    ierr = MatSetValues(jac,2,row,2,col,val,INSERT_VALUES);CHKERRQ(ierr);

    for(k=0; k < bus->ngen; k++) {
      /* Skip over OFF generators OR use status*x[i] ? */
      loc += 2;
      col[0] = loc;
      ierr = MatSetValues(jac,1,row,1,col,val,INSERT_VALUES);CHKERRQ(ierr);
      col[1] = loc+1;
      ierr = MatSetValues(jac,1,row+1,1,col+1,val,INSERT_VALUES);CHKERRQ(ierr);
    }

    ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);
    for(k=0; k < nconnlines; k++) {
      line = connlines[k];
      ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];

      ierr = PSBUSGetVariableGlobalLocation(busf,&locf);CHKERRQ(ierr);
      ierr = PSBUSGetVariableGlobalLocation(bust,&loct);CHKERRQ(ierr);

      if(bus == busf) {
	row[0] = mloc;
	col[0] = locf; col[1] = locf+1; col[2] = loct; col[3] = loct+1;
	ierr = MatSetValues(jac,1,row,4,col,val,INSERT_VALUES);CHKERRQ(ierr);
	row[0] = mloc+1;
	col[0] = locf; col[1] = locf+1; col[2] = loct; col[3] = loct+1;
	ierr = MatSetValues(jac,1,row,4,col,val,INSERT_VALUES);CHKERRQ(ierr);
      } else {
	row[0] = mloc;
	col[0] = loct; col[1] = loct+1; col[2] = locf; col[3] = locf+1;
	ierr = MatSetValues(jac,1,row,4,col,val,INSERT_VALUES);CHKERRQ(ierr);
	row[0] = mloc+1;
	col[0] = loct; col[1] = loct+1; col[2] = locf; col[3] = locf+1;
	ierr = MatSetValues(jac,1,row,4,col,val,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
    mloc += 2;
  }

  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  *mat = jac;
  ierr = PetscPrintf(comm,"Je structure:\n");CHKERRQ(ierr);
  ierr = MatView(jac,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWSetVariableBounds - Sets the bounds on variables

  Input Parameters:
+ opflow - the OPFLOW object
. Xl     - vector of lower bound on variables
- Xu     - vector of upper bound on variables

  This routine inserts the bounds on variables in the Xl and Xu vectors
*/
PetscErrorCode OPFLOWSetVariableBounds(OPFLOW opflow, Vec Xl, Vec Xu)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscScalar    *xl,*xu;
  PetscInt       i,k;
  PSBUS          bus;
  PetscInt       loc;
  DM             networkdm=ps->networkdm;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(Xu,&xu);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    if (bus->isghost) continue;

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);

    /* Bounds on voltage angles */
    xl[loc] = -PETSC_PI; xu[loc] = PETSC_PI;

    /* Bounds on voltage magnitudes and bounds on reactive power mismatch equality constraints */
    xl[loc+1] = bus->Vmin; xu[loc+1] = bus->Vmax;

    if(bus->ide == REF_BUS || bus->ide == ISOLATED_BUS) xl[loc] = xu[loc] = bus->va*PETSC_PI/180.0;
    if(bus->ide == ISOLATED_BUS) xl[loc+1] = xu[loc+1] = bus->vm;

    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      loc = loc+2;
      if(!gen->status) xl[loc] = xu[loc] = xl[loc+1] = xu[loc+1] = 0.0;
      else {
	xl[loc] = gen->pb;
	xu[loc] = gen->pt;
	xl[loc+1] = gen->qb;
	xu[loc+1] = gen->qt;
      }
    }
  }

  ierr = VecRestoreArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Xu,&xu);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Xl:\n");
  ierr = VecView(Xl,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Xu:\n");
  ierr = VecView(Xu,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*
  OPFLOWSetInitialGuess - Sets the initial guess for the optimization

  Input Parameters:
. opflow - the OPFLOW object

  Output Parameters:
+ X     - initial guess

  Notes:
   Sets X[i] = (Xl[i] + Xu[i])/2
*/
PetscErrorCode OPFLOWSetInitialGuess(OPFLOW opflow, Vec X)
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

    /* Guess for voltage angles and bounds on voltage magnitudes */
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

/*
  OPFLOWObjectiveandGradientFunction - The objective and gradient function for the optimal power flow

  Input Parameters:
+ nlp - the TAO object
. X      - the current iterate

  Output Parameters:
+ obj - the objective function value (scalar)
- grad - the gradient vector
*/
PetscErrorCode OPFLOWObjectiveandGradientFunction(Tao nlp,Vec X, PetscScalar* obj,Vec grad,void* ctx)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=(OPFLOW)ctx;
  PetscScalar    *df,Pg;
  PS             ps=opflow->ps;
  PetscInt       i,k;
  PSBUS          bus;
  PetscInt       loc;
  Vec            localX,localgrad;
  PSGEN          gen;
  const PetscScalar *x;

  PetscFunctionBegin;
  *obj = 0.0;

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ps->networkdm,&localgrad);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = VecSet(localgrad,0.0);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = VecGetArray(localgrad,&df);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);

    for(k=0; k < bus->ngen; k++) {
      loc = loc+2;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      Pg = x[loc]*ps->MVAbase;
      *obj += gen->cost_alpha*Pg*Pg + gen->cost_beta*Pg + gen->cost_gamma;
      df[loc] = ps->MVAbase*(2*gen->cost_alpha*Pg + gen->cost_beta);
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
  OPFLOW         opflow=(OPFLOW)ctx;
  PetscScalar    *g;
  PS             ps=opflow->ps;
  PetscInt       i,k;
  PSLINE         line;
  PSBUS          bus;
  PetscInt       gloc=0;
  const PSBUS    *connbuses;
  PetscInt       nconnlines;
  const PSLINE   *connlines;
  PSBUS          busf,bust;
  PetscInt       xloc,xlocf,xloct;
  Vec            localX;//,localGe;
  PSLOAD         load;
  PetscScalar    theta,Vm;
  const PetscScalar *x;

  PetscFunctionBegin;
  ierr = VecSet(Ge,0.0);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  //ierr = DMGetLocalVector(ps->networkdm,&localGe);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  //ierr = VecSet(Ge,0.0);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Ge,&g);CHKERRQ(ierr);
  /* Need to use DMLocaltoGlobal stuff in parallel */

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);
    theta = x[xloc];
    Vm    = x[xloc+1];

    if(bus->ide == ISOLATED_BUS) {
      g[gloc] =  theta - bus->va*PETSC_PI/180.0;
      g[gloc+1] = Vm - bus->vm;
      continue;
    }

    /* Shunt injections */
    g[gloc]   += Vm*Vm*bus->gl;
    g[gloc+1] += -Vm*Vm*bus->bl;

    for(k=0; k < bus->ngen; k++) {
      xloc += 2;
      PetscScalar Pg,Qg;
      Pg = x[xloc];
      Qg = x[xloc+1];

      g[gloc]   -= Pg;
      g[gloc+1] -= Qg;
    }

    for(k=0; k < bus->nload; k++) {
      ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);

      g[gloc]   += load->pl;
      g[gloc+1] += load->ql;
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

      PetscScalar Vmf,Vmt,thetaf,thetat,thetaft,thetatf;

      thetaf  = x[xlocf];
      Vmf     = x[xlocf+1];
      thetat  = x[xloct];
      Vmt     = x[xloct+1];
      thetaft = thetaf - thetat;
      thetatf = thetat - thetaf;

      PetscScalar Pf,Qf,Pt,Qt;

      if(bus == busf) {
	Pf = Gff*Vmf*Vmf  + Vmf*Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
	Qf = -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));

	g[gloc]   += Pf;
	g[gloc+1] += Qf;
      } else {
	Pt = Gtt*Vmt*Vmt  + Vmt*Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf));
	Qt = -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));

	g[gloc]   += Pt;
	g[gloc+1] += Qt;
      }
    }
    gloc += 2;

  }

  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Ge,&g);CHKERRQ(ierr);
  printf("localGe:\n");
  //VecView(localGe,0);

  //ierr = DMLocalToGlobalBegin(ps->networkdm,localGe,INSERT_VALUES,Ge);CHKERRQ(ierr);
  //ierr = DMLocalToGlobalEnd(ps->networkdm,localGe,INSERT_VALUES,Ge);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  //ierr = DMRestoreLocalVector(ps->networkdm,&localGe);CHKERRQ(ierr);
  
  ierr = VecView(Ge,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

/*
  OPFLOWInequalityConstraintsFunction - Evalulates the inequality constraints for the optimal power flow

  Input Parameters:
+ nlp - the TAO nlp solver object
. X   - the current iterate
- ctx - Application data set with TaoSetInequalityConstraintsRoutine

  Output Parameters:
. Gi  - vector of inequality constraints
*/
PetscErrorCode OPFLOWInequalityConstraintsFunction(Tao nlp,Vec X,Vec Gi,void* ctx)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=(OPFLOW)ctx;
  const PetscScalar *x;
  PetscScalar    *g;
  PS             ps=opflow->ps;
  PetscInt       i;
  PSLINE         line;
  PetscInt       gloc=0;
  const PSBUS    *connbuses;
  PSBUS          busf,bust;
  PetscInt       xlocf,xloct;


  PetscFunctionBegin;
  ierr = VecSet(Gi,0.0);CHKERRQ(ierr);

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Gi,&g);CHKERRQ(ierr);

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

    PetscScalar Pf,Qf,Pt,Qt,Sf2,St2;

    Pf = Gff*Vmf*Vmf  + Vmf*Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
    Qf = -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));

    Pt = Gtt*Vmt*Vmt  + Vmt*Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf));
    Qt = -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));

    Sf2 = Pf*Pf + Qf*Qf;
    St2 = Pt*Pt + Qt*Qt;

    PetscScalar Slim2;
    Slim2 = (line->rateA/ps->MVAbase)*(line->rateA/ps->MVAbase);

    /* TAO requires inequalities of the form
               h(x) >= 0
       The line flow inequality constraint is
           0 <= S <= Slim
       Converting it to TAO form, we get two constraints
          S >= 0 and Slim - S >= 0
    */
    g[gloc]   = Sf2;
    g[gloc+1] = -Sf2 + Slim2;
    g[gloc+2] = St2;
    g[gloc+3] = -St2 + Slim2;

    gloc = gloc + 4;
  }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gi,&g);CHKERRQ(ierr);
  //  ierr = VecView(Gi,0);CHKERRQ(ierr);
  //  exit(1);

  PetscFunctionReturn(0);
}


/*
  OPFLOWSolve - Solves the AC optimal power flow

  Input Parameters:
. opflow - the optimal power flow application object
*/
PetscErrorCode OPFLOWSolve(OPFLOW opflow)
{
  PetscErrorCode ierr;
  TaoConvergedReason reason;

  PetscFunctionBegin;
  if(!opflow->setupcalled) {
    ierr = OPFLOWSetUp(opflow);CHKERRQ(ierr);
  }

  /* Set variable bounds */
  ierr = OPFLOWSetVariableBounds(opflow,opflow->Xl,opflow->Xu);CHKERRQ(ierr);
  ierr = TaoSetVariableBounds(opflow->nlp,opflow->Xl,opflow->Xu);CHKERRQ(ierr);

  /* Set Initial Guess */
  ierr = OPFLOWSetInitialGuess(opflow,opflow->X);CHKERRQ(ierr);

  ierr = TaoSolve(opflow->nlp);CHKERRQ(ierr);
  ierr = TaoGetConvergedReason(opflow->nlp,&reason);CHKERRQ(ierr);
  opflow->converged = reason< 0 ? PETSC_FALSE:PETSC_TRUE;

  ierr = VecView(opflow->X,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWSetUp - Sets up an optimal power flow application object

  Input Parameters:
. opflow - the OPFLOW object

  Notes:
  This routine sets up the OPFLOW object and the underlying PS object. It
  also distributes the PS object when used in parallel.
*/
PetscErrorCode OPFLOWSetUp(OPFLOW opflow)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscMPIInt    rank;
  PetscInt       i,nbus=0;
  PSBUS          bus;

  PetscFunctionBegin;
  /* Set up PS object */
  ierr = PSSetUp(ps);CHKERRQ(ierr);

  /* Get local nconeq */
  ierr = MPI_Comm_rank(ps->comm->type,&rank);CHKERRQ(ierr);
  for (i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    if (!bus->isghost) nbus++;
  }
  printf("[%d] nbus/Nbus %d (%d), %d; nbranch/Nbranch %d %d\n",rank,ps->nbus,nbus,ps->Nbus,ps->nbranch,ps->Nbranch);
  ierr = MPI_Barrier(ps->comm->type);CHKERRQ(ierr);

  opflow->nconeq   = 2*nbus;
  opflow->Nconeq   = 2*ps->Nbus;
  opflow->nconineq = 2*2*ps->nbranch;
  opflow->Nconineq = 2*2*ps->Nbranch; /* 0 <= Sf2 <= Smax2, 0 <= St2 <= Smax2 */
  printf("[%d] nconeq/Nconeq %d, %d; nconineq/Nconineq %d, %d\n",rank,opflow->nconeq,opflow->Nconeq,opflow->nconineq,opflow->Nconineq);

  /* Create the solution vector */
  ierr = PSCreateGlobalVector(opflow->ps,&opflow->X);CHKERRQ(ierr);
  ierr = VecView(opflow->X,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  Vec vtmp;
  ierr = VecDuplicate(opflow->X,&vtmp);CHKERRQ(ierr);
  ierr = VecView(vtmp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = VecDestroy(&vtmp);CHKERRQ(ierr);

  ierr = TaoSetInitialVector(opflow->nlp,opflow->X);CHKERRQ(ierr);

  /* Get the size of the solution vector */
  ierr = VecGetLocalSize(opflow->X,&opflow->nvar);CHKERRQ(ierr);
  ierr = VecGetSize(opflow->X,&opflow->Nvar);CHKERRQ(ierr);

  /* Create the vector for upper and lower bounds on X */
  ierr = VecDuplicate(opflow->X,&opflow->Xl);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->Xu);CHKERRQ(ierr);

  /* Create the equality constraint vector */
  ierr = VecCreate(opflow->comm->type,&opflow->Ge);CHKERRQ(ierr);
  //ierr = VecSetSizes(opflow->Ge,PETSC_DETERMINE,opflow->Nconeq);CHKERRQ(ierr);
  ierr = VecSetSizes(opflow->Ge,opflow->nconeq,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(opflow->Ge);CHKERRQ(ierr);

  /* Create the inequality constraint vector */
  ierr = VecCreate(opflow->comm->type,&opflow->Gi);CHKERRQ(ierr);
  //ierr = VecSetSizes(opflow->Gi,opflow->Nconineq,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetSizes(opflow->Gi,opflow->nconineq,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(opflow->Gi);CHKERRQ(ierr);

  /* Create the equality constraint Jacobian */
  ierr = OPFLOWCreateEqualityConstraintsJacobian(opflow,&opflow->Jac_Ge);CHKERRQ(ierr);

  /* Create the inequality constraint Jacobian */
  ierr = OPFLOWCreateInequalityConstraintsJacobian(opflow,&opflow->Jac_Gi);CHKERRQ(ierr);

  /* Set the different callbacks Tao requires */
  /* Objective and Gradient */
  ierr = TaoSetObjectiveAndGradientRoutine(opflow->nlp,OPFLOWObjectiveandGradientFunction,(void*)opflow);CHKERRQ(ierr);

  /* Equality Constraints */
  ierr = TaoSetEqualityConstraintsRoutine(opflow->nlp,opflow->Ge,OPFLOWEqualityConstraintsFunction,(void*)opflow);CHKERRQ(ierr);

  /* Inequality Constraints */
  ierr = TaoSetInequalityConstraintsRoutine(opflow->nlp,opflow->Gi,OPFLOWInequalityConstraintsFunction,(void*)opflow);CHKERRQ(ierr);

  /* Equality Jacobian */
  ierr = TaoSetJacobianEqualityRoutine(opflow->nlp,opflow->Jac_Ge,opflow->Jac_Ge,OPFLOWEqualityConstraintsJacobianFunction,(void*)opflow);CHKERRQ(ierr);

  /* Inequality Jacobian */
  ierr = TaoSetJacobianInequalityRoutine(opflow->nlp,opflow->Jac_Gi,opflow->Jac_Gi,OPFLOWInequalityConstraintsJacobianFunction,(void*)opflow);CHKERRQ(ierr);
  ierr = TaoSetFromOptions(opflow->nlp);CHKERRQ(ierr);

  opflow->setupcalled = PETSC_TRUE;
  PetscFunctionReturn(0);
}
