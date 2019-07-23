#include <private/opflowimpl.h>

/****************************************************
  Designed to be a side file dedicated to
  the functions dealing with the inequality constraints
  Order to file:
    OPFLOWInequalityConstraintsFunction
    OPFLOWCreateInequalityConstraintsJacobian
    OPFLOWInequalityConstraintsJacobianFunction
*****************************************************/

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
  Vec            localX = opflow->localX;
  PSLINE         line;
  PSBUS          busf,bust;
  const PSBUS    *connbuses;
  const PetscScalar *x;

  PetscFunctionBegin;
  ierr = VecSet(Gi,0.0);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);
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
  Vec            localX=opflow->localX;
  PSLINE         line;
  PSBUS          busf,bust;
  const PSBUS    *connbuses;
  const PetscScalar *x;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);

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
  ierr = MatAssemblyBegin(Ji,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Ji,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
