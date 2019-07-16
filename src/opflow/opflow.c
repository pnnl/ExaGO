#include <private/psimpl.h>
#include <private/opflowimpl.h>
#include <petsc/private/matimpl.h>

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
  PetscInt       i,k,row[2],col[4],genctr,gidx;
  PetscInt       nconnlines,locglob,loc,locglobf,locglobt,locf,loct;
  PetscScalar    Vm,val[8],Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  PetscScalar    thetaf,thetat,Vmf,Vmt,thetaft,thetatf;
  OPFLOW         opflow=(OPFLOW)ctx;
  Vec            localX;
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

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
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
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(Je,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Je,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
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
  PetscInt       nconineq=opflow->nconineq; /* Number of local constraints */
  PetscInt       nvar=opflow->nvar; /* Number of local variables */
  PetscInt       i,mloc=0,locf,loct;
  PetscInt       rstart,rend,row[2],col[4];
  PetscInt       *dnnz,*onnz;
  PetscScalar    val[8];
  MPI_Comm       comm=opflow->comm->type;
  Mat            jac;
  PS             ps=opflow->ps;
  PSLINE         line;
  PSBUS          busf,bust;
  const PSBUS    *connbuses;

  PetscFunctionBegin;
  ierr = MatCreate(comm,&jac);CHKERRQ(ierr);
  ierr = MatSetSizes(jac,nconineq,nvar,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = MatSetType(jac,MATAIJ);CHKERRQ(ierr);

  /* Set up preallocation */
  ierr = PetscCalloc2(nconineq,&dnnz,nconineq,&onnz);CHKERRQ(ierr);
  mloc = 0;
  for (i=0; i < ps->nbranch; i++) {
    line = &ps->line[i];
    ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
    busf = connbuses[0];
    bust = connbuses[1];

    if (busf->isghost && bust->isghost) {
      onnz[mloc]   += 4;
      onnz[mloc+1] += 4;
      onnz[mloc+2] += 4;
      onnz[mloc+3] += 4;
    } else if (!busf->isghost && !bust->isghost) {
      dnnz[mloc]   += 4;
      dnnz[mloc+1] += 4;
      dnnz[mloc+2] += 4;
      dnnz[mloc+3] += 4;
    } else {
      onnz[mloc]   += 2;
      onnz[mloc+1] += 2;
      onnz[mloc+2] += 2;
      onnz[mloc+3] += 2;
      dnnz[mloc]   += 2;
      dnnz[mloc+1] += 2;
      dnnz[mloc+2] += 2;
      dnnz[mloc+3] += 2;
    }
    mloc += 4;
  }

  ierr = MatSeqAIJSetPreallocation(jac,0,dnnz);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(jac,0,dnnz,0,onnz);CHKERRQ(ierr);
  ierr = PetscFree2(dnnz,onnz);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(jac,&rstart,&rend);CHKERRQ(ierr);
  for (i=0; i<8; i++) val[i] = 0.0;

  /* Set up nonzero structure. 0 is inserted in the non-zero location */
  mloc = rstart;
  for (i=0; i < ps->nbranch; i++) {
    line = &ps->line[i];
    ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
    busf = connbuses[0];
    bust = connbuses[1];

    ierr = PSBUSGetVariableGlobalLocation(busf,&locf);CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bust,&loct);CHKERRQ(ierr);

    /* From bus Sf lower and upper bound */
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
  ierr = MatSetOption(jac,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);

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
  PetscInt       ctr=0,i;
  PetscInt       row[2],col[4];
  PetscInt       rstart,rend;
  PetscInt       gloc=0,xlocf,xloct,glocf,gloct;
  PetscMPIInt    rank,size;
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
  OPFLOW         opflow=(OPFLOW)ctx;
  PS             ps=opflow->ps;
  MPI_Comm       comm=opflow->comm->type;
  Vec            localX;
  PSLINE         line;
  PSBUS          busf,bust;
  const PSBUS    *connbuses;
  const PetscScalar *x;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(Ji,&rstart,&rend);CHKERRQ(ierr);
  ierr = MPI_Barrier(comm);CHKERRQ(ierr);

  gloc = rstart;
  for (i=0; i < ps->nbranch; i++) {
    line = &ps->line[i];

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
    ierr = PSBUSGetVariableGlobalLocation(busf,&glocf);CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bust,&gloct);CHKERRQ(ierr);

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
    col[0] = glocf; col[1] = glocf+1; col[2] = gloct; col[3] = gloct+1;
    val[0] = dSf2_dthetaf;
    val[1] = dSf2_dVmf;
    val[2] = dSf2_dthetat;
    val[3] = dSf2_dVmt;
    ierr = MatSetValues(Ji,1,row,4,col,val,INSERT_VALUES);CHKERRQ(ierr);

    row[0] = gloc+1;
    col[0] = glocf; col[1] = glocf+1; col[2] = gloct; col[3] = gloct+1;
    val[0] = -dSf2_dthetaf;
    val[1] = -dSf2_dVmf;
    val[2] = -dSf2_dthetat;
    val[3] = -dSf2_dVmt;
    ierr = MatSetValues(Ji,1,row,4,col,val,INSERT_VALUES);CHKERRQ(ierr);

    dSt2_dthetaf = dSt2_dPt*dPt_dthetaf + dSt2_dQt*dQt_dthetaf;
    dSt2_dthetat = dSt2_dPt*dPt_dthetat + dSt2_dQt*dQt_dthetat;
    dSt2_dVmf    = dSt2_dPt*dPt_dVmf    + dSt2_dQt*dQt_dVmf;
    dSt2_dVmt    = dSt2_dPt*dPt_dVmt    + dSt2_dQt*dQt_dVmt;

    row[0] = gloc+2;
    col[0] = gloct; col[1] = gloct+1; col[2] = glocf; col[3] = glocf+1;
    val[ctr]   = dSt2_dthetat;
    val[ctr+1] = dSt2_dVmt;
    val[ctr+2] = dSt2_dthetaf;
    val[ctr+3] = dSt2_dVmf;
    ierr = MatSetValues(Ji,1,row,4,col,val,INSERT_VALUES);CHKERRQ(ierr);

    row[0] = gloc+3;
    col[0] = gloct; col[1] = gloct+1; col[2] = glocf; col[3] = glocf+1;
    val[ctr]   = -dSt2_dthetat;
    val[ctr+1] = -dSt2_dVmt;
    val[ctr+2] = -dSt2_dthetaf;
    val[ctr+3] = -dSt2_dVmf;
    ierr = MatSetValues(Ji,1,row,4,col,val,INSERT_VALUES);CHKERRQ(ierr);

    gloc += 4;
  }

  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

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
  //ierr = MatSetUp(jac);CHKERRQ(ierr);
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
    } else {
      val[0] = val[1] = 2.0;
      ierr = VecSetValues(Vonz,2,row,val,ADD_VALUES);CHKERRQ(ierr);
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

  /* TODO: 'mpiexec -n 5 ./OPFLOW -petscpartitioner_type simple' throws an error; 
     Thus enable 'NEW_NONZERO_ALLOCATION' below until the bug is fixed. */
  //ierr = MatSetOption(jac,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);CHKERRQ(ierr);
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

  *mat = jac;
  //ierr = PetscPrintf(comm,"Je structure:\n");CHKERRQ(ierr);
  //ierr = MatView(jac,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
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
  PetscScalar    vals[2];
  PetscInt       i,k,row[2];
  PSBUS          bus;
  PetscInt       loc;

  PetscFunctionBegin;
  /* Get array pointers */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    if (bus->isghost) continue;
    ierr = PSBUSGetVariableGlobalLocation(bus,&loc);CHKERRQ(ierr);

    /* Bounds on voltage angles */
    row[0] = loc; row[1] = loc+1;
    vals[0] = -PETSC_PI; vals[1] = PETSC_PI;
    ierr = VecSetValues(Xl,1,row,vals,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValues(Xu,1,row,vals+1,INSERT_VALUES);CHKERRQ(ierr);

    /* Bounds on voltage magnitudes and bounds on reactive power mismatch equality constraints */
    vals[0] = bus->Vmin; vals[1] = bus->Vmax;
    ierr = VecSetValues(Xl,1,row+1,vals,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValues(Xu,1,row+1,vals+1,INSERT_VALUES);CHKERRQ(ierr);

    if(bus->ide == REF_BUS || bus->ide == ISOLATED_BUS){
      vals[0] = bus->va*PETSC_PI/180.0;
      ierr = VecSetValues(Xl,1,row,vals,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecSetValues(Xu,1,row,vals,INSERT_VALUES);CHKERRQ(ierr);
    }
    if(bus->ide == ISOLATED_BUS){
      vals[0] = bus->vm;
      ierr = VecSetValues(Xl,1,row+1,vals,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecSetValues(Xu,1,row+1,vals,INSERT_VALUES);CHKERRQ(ierr);
    }

    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      row[0] = row[0]+2; row[1] = row[1]+2;
      if(!gen->status){
        //xl[loc] = xu[loc] = xl[loc+1] = xu[loc+1] = 0.0;
        vals[0] =0.0; vals[1] = 0.0;
        ierr = VecSetValues(Xl,2,row,vals,INSERT_VALUES);CHKERRQ(ierr);
        ierr = VecSetValues(Xu,2,row,vals,INSERT_VALUES);CHKERRQ(ierr);
      }
      else {
	//xl[loc] = gen->pb;  xu[loc] = gen->pt;
   vals[0] = gen->pb; vals[1] = gen->pt;
   ierr = VecSetValues(Xl,1,row,vals,INSERT_VALUES);CHKERRQ(ierr);
   ierr = VecSetValues(Xu,1,row,vals+1,INSERT_VALUES);CHKERRQ(ierr);
	//xl[loc+1] = gen->qb; xu[loc+1] = gen->qt;
   vals[0] = gen->qb; vals[1] = gen->qt;
   ierr = VecSetValues(Xl,1,row+1,vals,INSERT_VALUES);CHKERRQ(ierr);
   ierr = VecSetValues(Xu,1,row+1,vals+1,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }
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
  PetscInt       i,row[2];
  PetscInt       loc,gloc,k;
  PetscMPIInt    rank;
  PetscScalar    vals[2];
  PS             ps=opflow->ps;
  MPI_Comm       comm=opflow->comm->type;
  Vec            localXl,localXu;
  PSBUS          bus;
  PSGEN          gen;
  const PetscScalar *xl,*xu;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = DMGetLocalVector(ps->networkdm,&localXl);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ps->networkdm,&localXu);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,opflow->Xl,INSERT_VALUES,localXl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,opflow->Xu,INSERT_VALUES,localXu);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,opflow->Xl,INSERT_VALUES,localXl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,opflow->Xu,INSERT_VALUES,localXu);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localXl,&xl);CHKERRQ(ierr);
  ierr = VecGetArrayRead(localXu,&xu);CHKERRQ(ierr);

  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {

    bus = &ps->bus[i];
    if (bus->isghost)continue;

    ierr = PSBUSGetVariableGlobalLocation(bus,&gloc);CHKERRQ(ierr);
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    /* Guess for voltage angles and bounds on voltage magnitudes */

    row[0] = gloc; row[1] = gloc+1;
    vals[0] = (xl[loc] + xu[loc])/2.0;
    vals[1] = (xl[loc+1] + xu[loc+1])/2.0;
    ierr = VecSetValues(X,2,row,vals,INSERT_VALUES);CHKERRQ(ierr);

    for(k=0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      loc = loc+2; row[0] = row[0]+2; row[1] = row[1]+2;
      vals[0] = (xl[loc] + xu[loc])/2.0;
      vals[1] = (xl[loc+1] + xu[loc+1])/2.0;
      ierr = VecSetValues(X,2,row,vals,INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  ierr = VecRestoreArrayRead(localXl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(localXu,&xu);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localXl);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localXu);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(X);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(X);CHKERRQ(ierr);
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
  PetscInt       i,k;
  PetscInt       loc;
  PetscScalar    *df,Pg,sobj;
  OPFLOW         opflow=(OPFLOW)ctx;
  PS             ps=opflow->ps;
  Vec            localX,localgrad;
  PSBUS          bus;
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
      sobj += gen->cost_alpha*Pg*Pg + gen->cost_beta*Pg + gen->cost_gamma;
      df[loc] = ps->MVAbase*(2*gen->cost_alpha*Pg + gen->cost_beta);
    }
  }
  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(localgrad,&df);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(ps->networkdm,localgrad,INSERT_VALUES,grad);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(ps->networkdm,localgrad,INSERT_VALUES,grad);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localgrad);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&sobj,obj,1,MPI_DOUBLE,MPI_SUM,opflow->comm->type);CHKERRQ(ierr);

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
  Vec            localX;
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

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
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

  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
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
  PetscInt       i;
  PetscInt       gloc=0;
  PetscInt       xlocf,xloct;
  PetscScalar    *g;
  PetscScalar    Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  PetscScalar    Vmf,Vmt,thetaf,thetat,thetaft,thetatf;
  PetscScalar    Pf,Qf,Pt,Qt,Sf2,St2;
  PetscScalar    Slim2;
  OPFLOW         opflow=(OPFLOW)ctx;
  PS             ps=opflow->ps;
  Vec            localX;
  PSLINE         line;
  PSBUS          busf,bust;
  const PSBUS    *connbuses;
  const PetscScalar *x;

  PetscFunctionBegin;
  ierr = VecSet(Gi,0.0);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Gi,&g);CHKERRQ(ierr);

  for(i=0; i < ps->nbranch; i++) {
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

  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gi,&g);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
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

  ierr = VecView(opflow->X,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
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
  PetscInt       i,nbus=0;
  PetscMPIInt    rank;
  PS             ps=opflow->ps;
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
  ierr = PetscSynchronizedPrintf(ps->comm->type,"[%d] nbus/Nbus %d (%d !ghost), %d; nbranch/Nbranch %d %d\n",rank,ps->nbus,nbus,ps->Nbus,ps->nbranch,ps->Nbranch);CHKERRQ(ierr);
  ierr = PetscSynchronizedFlush(ps->comm->type,PETSC_STDOUT);CHKERRQ(ierr);
  ierr = MPI_Barrier(ps->comm->type);CHKERRQ(ierr);

  opflow->nconeq   = 2*nbus;
  opflow->Nconeq   = 2*ps->Nbus;
  opflow->nconineq = 2*2*ps->nbranch;
  opflow->Nconineq = 2*2*ps->Nbranch; /* 0 <= Sf2 <= Smax2, 0 <= St2 <= Smax2 */
  printf("[%d] nconeq/Nconeq %d, %d; nconineq/Nconineq %d, %d\n",rank,opflow->nconeq,opflow->Nconeq,opflow->nconineq,opflow->Nconineq);

  /* Create the solution vector */
  ierr = PSCreateGlobalVector(opflow->ps,&opflow->X);CHKERRQ(ierr);

  ierr = TaoSetInitialVector(opflow->nlp,opflow->X);CHKERRQ(ierr);

  /* Get the size of the solution vector */
  ierr = VecGetLocalSize(opflow->X,&opflow->nvar);CHKERRQ(ierr);
  ierr = VecGetSize(opflow->X,&opflow->Nvar);CHKERRQ(ierr);

  /* Create the vector for upper and lower bounds on X */
  ierr = VecDuplicate(opflow->X,&opflow->Xl);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->Xu);CHKERRQ(ierr);

  /* Create the equality constraint vector */
  ierr = VecCreate(opflow->comm->type,&opflow->Ge);CHKERRQ(ierr);
  ierr = VecSetSizes(opflow->Ge,opflow->nconeq,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(opflow->Ge);CHKERRQ(ierr);
  ierr = DMNetworkSetVertexLocalToGlobalOrdering(ps->networkdm);CHKERRQ(ierr);

  /* Create the inequality constraint vector */
  ierr = VecCreate(opflow->comm->type,&opflow->Gi);CHKERRQ(ierr);
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
