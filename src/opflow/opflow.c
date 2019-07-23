#include <private/opflowimpl.h>

/*********************************************
  The optimization problem for Tao should be in
  the following form

  min f(x)
  s.t.
    g(x)  = 0
    h(x) >= 0
    xl <= x <= xu

************************************/

#include <../src/tao/constrained/impls/ipm/ipm.h> /*I "ipm.h" I*/

/*
  f(x) = SUM(a*(Pg*MVAbase)^2 + b*(Pg*MVAbase) + c); for each xi
  H(x) = f_xx = constant diagonal matrix: 2*a*MVAbase*MVAbase;
 */
PetscErrorCode OPFLOWCreateObjectiveHessian(OPFLOW opflow)
{
  PetscErrorCode ierr;
  PetscInt       i,k,n,loc;
  PetscScalar    val;
  PS             ps=opflow->ps;
  PSBUS          bus;
  PSGEN          gen;
  Tao            tao=opflow->nlp;
  Mat            H=tao->hessian;

  PetscFunctionBegin;
  ierr = VecGetLocalSize(opflow->X,&n);CHKERRQ(ierr);
  ierr = MatSetSizes(H,n,n,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = MatSetType(H,MATAIJ);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(H,1,NULL);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(H,1,NULL,0,NULL);CHKERRQ(ierr);

  for (i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableGlobalLocation(bus,&loc);CHKERRQ(ierr);
    for (k=0; k < bus->ngen; k++) {
      loc = loc + 2;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if (!gen->status) continue;
      val = 2*ps->MVAbase*ps->MVAbase*gen->cost_alpha;
      ierr = MatSetValue(H,loc,loc,val,INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  ierr = MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  //ierr = MatView(H,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWHessian(Tao tao, Vec X, Mat H, Mat H_pre, void* ctx)
{
  PetscFunctionBegin;
  H     = tao->hessian;
  H_pre = tao->hessian;
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
  ierr = VecDestroy(&(*opflow)->localX);CHKERRQ(ierr);

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
  PSGEN          gen;
  PetscInt       loc;

  PetscFunctionBegin;
  /* Get array pointers */
  for (i=0; i < ps->nbus; i++) {
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

    if (bus->ide == REF_BUS || bus->ide == ISOLATED_BUS){
      vals[0] = bus->va*PETSC_PI/180.0;
      ierr = VecSetValues(Xl,1,row,vals,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecSetValues(Xu,1,row,vals,INSERT_VALUES);CHKERRQ(ierr);
    }
    if (bus->ide == ISOLATED_BUS){
      vals[0] = bus->vm;
      ierr = VecSetValues(Xl,1,row+1,vals,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecSetValues(Xu,1,row+1,vals,INSERT_VALUES);CHKERRQ(ierr);
    }

    for (k=0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      row[0] = row[0]+2; row[1] = row[1]+2;
      if (!gen->status){
        vals[0] =0.0; vals[1] = 0.0;
        ierr = VecSetValues(Xl,2,row,vals,INSERT_VALUES);CHKERRQ(ierr);
        ierr = VecSetValues(Xu,2,row,vals,INSERT_VALUES);CHKERRQ(ierr);
      } else {
        vals[0] = gen->pb; vals[1] = gen->pt;
        ierr = VecSetValues(Xl,1,row,vals,INSERT_VALUES);CHKERRQ(ierr);
        ierr = VecSetValues(Xu,1,row,vals+1,INSERT_VALUES);CHKERRQ(ierr);

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

  PetscFunctionBegin;
  ierr = VecAXPBYPCZ(X,0.5,0.5,0.0,opflow->Xl,opflow->Xu);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWObjectiveandGradientFunction - The objective and gradient function for the optimal power flow
    f(x) = SUM(a*(Pg*MVAbase)^2 + b*(Pg*MVAbase) + c); for each xi
    grad = f_x st. grad[i] = 2*a*(Pg*MVAbase)*MVAbase + b*MVAbase

  Input Parameters:
+ nlp - the TAO object
. X      - the current iterate

  Output Parameters:
+ obj - the objective function value (scalar)
- grad - the gradient vector
*/
PetscErrorCode OPFLOWObjectiveandGradientFunction(Tao nlp,Vec X,PetscScalar* obj,Vec grad,void* ctx)
{
  PetscErrorCode ierr;
  PetscInt       i,k;
  PetscInt       loc;
  PetscScalar    *df,Pg,sobj;
  OPFLOW         opflow=(OPFLOW)ctx;
  PS             ps=opflow->ps;
  Vec            localgrad;
  Vec            localX=opflow->localX;
  PSBUS          bus;
  PSGEN          gen;
  const PetscScalar *x;

  PetscFunctionBegin;
  *obj = 0.0;

  ierr = VecDuplicate(localX,&localgrad);CHKERRQ(ierr);

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

  ierr = VecDestroy(&localgrad);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&sobj,obj,1,MPI_DOUBLE,MPI_SUM,opflow->comm->type);CHKERRQ(ierr);
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
#if defined(PETSC_USE_LOG)
  PetscLogStage stages[3];
#endif

  PetscFunctionBegin;
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  ierr = PetscLogStageRegister(" OPFLOWSetUp",&stages[0]);CHKERRQ(ierr);
  ierr = PetscLogStageRegister(" TaoSetUp",&stages[1]);CHKERRQ(ierr);
  ierr = PetscLogStageRegister(" TaoSolve",&stages[2]);CHKERRQ(ierr);

  if (!opflow->setupcalled) {
    ierr = PetscLogStagePush(stages[0]);CHKERRQ(ierr);
    ierr = OPFLOWSetUp(opflow);CHKERRQ(ierr);
    ierr = PetscLogStagePop();CHKERRQ(ierr);
  }

  ierr = PetscLogStagePush(stages[1]);CHKERRQ(ierr);

  /* Set variable bounds */
  ierr = OPFLOWSetVariableBounds(opflow,opflow->Xl,opflow->Xu);CHKERRQ(ierr);
  ierr = TaoSetVariableBounds(opflow->nlp,opflow->Xl,opflow->Xu);CHKERRQ(ierr);

  /* Set Initial Guess */
  ierr = OPFLOWSetInitialGuess(opflow,opflow->X);CHKERRQ(ierr);

  ierr = TaoSetUp(opflow->nlp);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  ierr = PetscLogStagePush(stages[2]);CHKERRQ(ierr);
  ierr = TaoSolve(opflow->nlp);CHKERRQ(ierr);
  ierr = TaoGetConvergedReason(opflow->nlp,&reason);CHKERRQ(ierr);
  opflow->converged = reason< 0 ? PETSC_FALSE:PETSC_TRUE;

  PetscBool solu_view=PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL, "-solu_view", &solu_view,NULL);CHKERRQ(ierr);
  if (solu_view) {
    ierr = VecView(opflow->X,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
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
  ierr = PetscSynchronizedPrintf(ps->comm->type,"[%d] nbus %d (%d !ghost), Nbus %d; nbranch %d, Nbranch %d\n",rank,ps->nbus,nbus,ps->Nbus,ps->nbranch,ps->Nbranch);CHKERRQ(ierr);
  ierr = PetscSynchronizedFlush(ps->comm->type,PETSC_STDOUT);CHKERRQ(ierr);

  opflow->nconeq   = 2*nbus;
  opflow->Nconeq   = 2*ps->Nbus;
  opflow->nconineq = 2*2*ps->nbranch;
  opflow->Nconineq = 2*2*ps->Nbranch; /* 0 <= Sf2 <= Smax2, 0 <= St2 <= Smax2 */
  ierr = PetscSynchronizedPrintf(ps->comm->type,"[%d] nconeq/Nconeq %d, %d; nconineq/Nconineq %d, %d\n",rank,opflow->nconeq,opflow->Nconeq,opflow->nconineq,opflow->Nconineq);CHKERRQ(ierr);
  ierr = PetscSynchronizedFlush(ps->comm->type,PETSC_STDOUT);CHKERRQ(ierr);

  /* Create the solution vector */
  ierr = PSCreateGlobalVector(opflow->ps,&opflow->X);CHKERRQ(ierr);
  ierr = DMCreateLocalVector(ps->networkdm,&opflow->localX);CHKERRQ(ierr);

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

  /* End of Second Stage Start of Third */
  // ierr = PetscLogStagePop();CHKERRQ(ierr);
  // ierr = PetscLogStagePush(stages[2]);CHKERRQ(ierr);

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

  /* Create Hessian for objective function and set Hessian routine */
  ierr = OPFLOWCreateObjectiveHessian(opflow);CHKERRQ(ierr);
  ierr = TaoSetHessianRoutine(opflow->nlp,NULL,NULL,OPFLOWHessian,NULL);CHKERRQ(ierr);

  opflow->setupcalled = PETSC_TRUE;
  PetscFunctionReturn(0);
}
