#include <private/opflowimpl.h>
#include <private/tcopflowimpl.h>

/*
  TCOPFLOWSetLoadProfiles - Sets the data files for time-varying load profiles

  Input Parameter
+  tcopflow - The TCOPFLOW object
.  ploadproffile - The name of the file with real power load variationt
-  qloadproffile - The name of the file with reactive power load variationt
*/
PetscErrorCode TCOPFLOWSetLoadProfiles(TCOPFLOW tcopflow,const char ploadprofile[],const char qloadprofile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  
  if(ploadprofile) {
    ierr = PetscMemcpy(tcopflow->ploadprofile,ploadprofile,100*sizeof(char));CHKERRQ(ierr);
    tcopflow->ploadprofileset = PETSC_TRUE;
  }

  if(qloadprofile) {
    ierr = PetscMemcpy(tcopflow->qloadprofile,qloadprofile,100*sizeof(char));CHKERRQ(ierr);
    tcopflow->qloadprofileset = PETSC_TRUE;
  }

  PetscFunctionReturn(0);
}

/*
  TCOPFLOWWindGenProfiles - Sets the data file for wind generation profiles

  Input Parameter
+  tcopflow - The TCOPFLOW object
-  windgenproffile - The name of the file with wind generation profile
*/
PetscErrorCode TCOPFLOWSetWindGenProfiles(TCOPFLOW tcopflow,const char windgenprofile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  
  if(windgenprofile) {
    ierr = PetscMemcpy(tcopflow->windgenprofile,windgenprofile,100*sizeof(char));CHKERRQ(ierr);
    tcopflow->windgenprofileset = PETSC_TRUE;
  }

  PetscFunctionReturn(0);
}


/*
  TCOPFLOWCreate - Creates a multi-period optimal power flow application object

  Input Parameters
. mpicomm - The MPI communicator

  Output Parameters
. tcopflowout - The multi-period optimal power flow application object
*/
PetscErrorCode TCOPFLOWCreate(MPI_Comm mpicomm, TCOPFLOW *tcopflowout)
{
  PetscErrorCode ierr;
  TCOPFLOW         tcopflow;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&tcopflow);CHKERRQ(ierr);

  ierr = COMMCreate(mpicomm,&tcopflow->comm);CHKERRQ(ierr);

  tcopflow->Nconeq   = tcopflow->nconeq  = 0;
  tcopflow->Nconineq = tcopflow->nconineq = 0;
  tcopflow->Ncon     = tcopflow->ncon     = 0;
  tcopflow->Nx       = tcopflow->nx       = 0;
  tcopflow->Gi       = NULL;
  tcopflow->Lambdai  = NULL;
  
  tcopflow->obj_factor = 1.0;
  tcopflow->obj = 0.0;

  tcopflow->duration = 1.0;
  tcopflow->dT       = 60.0;
  tcopflow->Nt       = 1;
  tcopflow->ploadprofileset = PETSC_FALSE;
  tcopflow->qloadprofileset = PETSC_FALSE;
  tcopflow->windgenprofileset = PETSC_FALSE;

  tcopflow->solver   = NULL;

  tcopflow->nsolversregistered = 0;
  tcopflow->TCOPFLOWSolverRegisterAllCalled = PETSC_FALSE;

  /* Register all solvers */
  ierr = TCOPFLOWSolverRegisterAll(tcopflow);

  tcopflow->solutiontops = PETSC_FALSE;
  /* Run-time optiont */
  tcopflow->iscoupling = PETSC_FALSE;

  tcopflow->setupcalled = PETSC_FALSE;
  *tcopflowout = tcopflow;

  ierr = PetscPrintf(tcopflow->comm->type,"TCOPFLOW: Application created\n");
  PetscFunctionReturn(0);
}

/*
  TCOPFLOWDestroy - Destroys the multi-period optimal power flow application object

  Input Parameter
. tcopflow - The TCOPFLOW object to destroy
*/
PetscErrorCode TCOPFLOWDestroy(TCOPFLOW *tcopflow)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  ierr = COMMDestroy(&(*tcopflow)->comm);CHKERRQ(ierr);

  /* Solution vector */
  ierr = VecDestroy(&(*tcopflow)->X);CHKERRQ(ierr);
  ierr = VecDestroy(&(*tcopflow)->localX);CHKERRQ(ierr);
  ierr = VecDestroy(&(*tcopflow)->gradobj);CHKERRQ(ierr);

  /* Lower and upper bounds on X */
  ierr = VecDestroy(&(*tcopflow)->Xl);CHKERRQ(ierr);
  ierr = VecDestroy(&(*tcopflow)->Xu);CHKERRQ(ierr);

  /* Conttraints vector */
  ierr = VecDestroy(&(*tcopflow)->G);CHKERRQ(ierr);
  ierr = VecDestroy(&(*tcopflow)->Ge);CHKERRQ(ierr);
  ierr = VecDestroy(&(*tcopflow)->Gelocal);CHKERRQ(ierr);
  ierr = VecDestroy(&(*tcopflow)->Lambdae);CHKERRQ(ierr);
  ierr = VecDestroy(&(*tcopflow)->Lambdaelocal);CHKERRQ(ierr);
  if((*tcopflow)->Nconineq) {
    ierr = VecDestroy(&(*tcopflow)->Gi);CHKERRQ(ierr);
    ierr = VecDestroy(&(*tcopflow)->Lambdai);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&(*tcopflow)->Gl);CHKERRQ(ierr);
  ierr = VecDestroy(&(*tcopflow)->Gu);CHKERRQ(ierr);
  ierr = VecDestroy(&(*tcopflow)->Lambda);CHKERRQ(ierr);

  /* Jacobian of conttraints */
  ierr = MatDestroy(&(*tcopflow)->Jac);CHKERRQ(ierr);
  ierr = MatDestroy(&(*tcopflow)->Jac_Ge);CHKERRQ(ierr);
  ierr = MatDestroy(&(*tcopflow)->Jac_Gi);CHKERRQ(ierr);

  ierr = MatDestroy(&(*tcopflow)->Hes);CHKERRQ(ierr);

  if((*tcopflow)->solverops.destroy) {
    ierr = ((*tcopflow)->solverops.destroy)(*tcopflow);
  }

  /* Destroy OPFLOW objects */
  for(i=0; i < (*tcopflow)->Nt; i++) {
    ierr = OPFLOWDestroy(&(*tcopflow)->opflows[i]);CHKERRQ(ierr);
  }

  ierr = PetscFree((*tcopflow)->nconineqcoup);CHKERRQ(ierr);
  ierr = PetscFree((*tcopflow)->opflows);CHKERRQ(ierr);
  ierr = PetscFree(*tcopflow);CHKERRQ(ierr);
  //  *tcopflow = 0;
  PetscFunctionReturn(0);
}

/*
  TCOPFLOWSetSolver - Sets the solver for TCOPFLOW

  Input Parameters:
+ tcopflow - tcopflow application object
- solvername - name of the solver
*/
PetscErrorCode TCOPFLOWSetSolver(TCOPFLOW tcopflow,const char* solvername)
{
  PetscErrorCode ierr,(*r)(TCOPFLOW)=NULL;
  PetscInt       i;
  PetscFunctionBegin;
  PetscBool match;
  for(i=0;i < tcopflow->nsolversregistered;i++) {
    ierr = PetscStrcmp(tcopflow->TCOPFLOWSolverList[i].name,solvername,&match);CHKERRQ(ierr);
    if(match) {
      r = tcopflow->TCOPFLOWSolverList[i].create;
      break;
    }
  }

  if(!r) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown type for TCOPFLOW Solver %s",solvername);

  /* Initialize (Null) the function pointers */
  tcopflow->solverops.destroy = 0;
  tcopflow->solverops.solve   = 0;
  tcopflow->solverops.setup   = 0;

  ierr = PetscStrcpy(tcopflow->solvername,solvername);CHKERRQ(ierr);
  /* Call the underlying implementation constructor */
  ierr = (*r)(tcopflow);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  TCOPFLOWSetNetworkData - Sets and reads the network data

  Input Parameter
+  tcopflow - The TCOPFLOW object
-  netfile - The name of the network file

  Notes: The input data must be in MATPOWER data format
*/
PetscErrorCode TCOPFLOWSetNetworkData(TCOPFLOW tcopflow,const char netfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscMemcpy(tcopflow->netfile,netfile,100*sizeof(char));CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  TCOPFLOWSetUp - Sets up an multi-period optimal power flow application object

  Input Parameters:
. tcopflow - the TCOPFLOW object

  Notes:
  This routine sets up the TCOPFLOW object and the underlying PS object. It
  also distributes the PS object when used in parallel.
*/
PetscErrorCode TCOPFLOWSetUp(TCOPFLOW tcopflow)
{
  PetscErrorCode ierr;
  PetscBool      solverset;
  char           modelname[32]="POWER_BALANCE_POLAR";
  char           solvername[32]="IPOPT";
  PetscInt       i,j;
  PS             ps;
  PetscBool      flg;

  PetscFunctionBegin;

  ierr = PetscOptionsBegin(tcopflow->comm->type,NULL,"TCOPFLOW options",NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-tcopflow_model","TCOPFLOW model type","",modelname,modelname,32,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-tcopflow_solver","TCOPFLOW solver type","",solvername,solvername,32,&solverset);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-tcopflow_iscoupling","Include coupling between first stage and second stage","",tcopflow->iscoupling,&tcopflow->iscoupling,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-tcopflow_dT","Length of time-step (minutes)","",tcopflow->dT,&tcopflow->dT,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-tcopflow_duration","Time horizon (hours)","",tcopflow->duration,&tcopflow->duration,NULL);CHKERRQ(ierr);
  PetscOptionsEnd();

  tcopflow->Nt = round(tcopflow->duration*60.0/tcopflow->dT);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"TCOPFLOW: Duration = %lf hours, timestep = %lf minutes, number of time-steps = %d\n",tcopflow->duration,tcopflow->dT,tcopflow->Nt);CHKERRQ(ierr);

  /* Set solver */
  if(solverset) {
    if(tcopflow->solver) ierr = (*tcopflow->solverops.destroy)(tcopflow);
    ierr = TCOPFLOWSetSolver(tcopflow,solvername);CHKERRQ(ierr);
    ierr = PetscPrintf(tcopflow->comm->type,"TCOPFLOW: Using %s solver\n",solvername);CHKERRQ(ierr);
  } else {
    if(!tcopflow->solver) {
      ierr = TCOPFLOWSetSolver(tcopflow,TCOPFLOWSOLVER_IPOPT);CHKERRQ(ierr);
      ierr = PetscPrintf(tcopflow->comm->type,"TCOPFLOW: Using %s solver\n",TCOPFLOWSOLVER_IPOPT);CHKERRQ(ierr); 
    }
  }

  /* Create OPFLOW objects */
  ierr = PetscCalloc1(tcopflow->Nt,&tcopflow->opflows);CHKERRQ(ierr);
  for(i=0; i < tcopflow->Nt; i++) {
    ierr = OPFLOWCreate(PETSC_COMM_SELF,&tcopflow->opflows[i]);CHKERRQ(ierr);
    ierr = OPFLOWSetModel(tcopflow->opflows[i],modelname);CHKERRQ(ierr);

    /* Read network data file */
    ierr = OPFLOWReadMatPowerData(tcopflow->opflows[i],tcopflow->netfile);CHKERRQ(ierr);
    /* Set up the PS object for opflow */
    ps = tcopflow->opflows[i]->ps;
    ierr = PSSetUp(ps);CHKERRQ(ierr);

    /* Set up OPFLOW object */
    ierr = OPFLOWSetUp(tcopflow->opflows[i]);CHKERRQ(ierr);
  }

  /* Read Pload profiles */
  if(tcopflow->ploadprofileset) {
    ierr = TCOPFLOWReadPloadProfile(tcopflow,tcopflow->ploadprofile);
  }

  /* Read Qload profiles */
  if(tcopflow->qloadprofileset) {
    ierr = TCOPFLOWReadQloadProfile(tcopflow,tcopflow->qloadprofile);
  }

  /* Read wind generation profiles */
  if(tcopflow->windgenprofileset) {
    ierr = TCOPFLOWReadWindGenProfile(tcopflow,tcopflow->windgenprofile);
  }

  
  ierr = PetscCalloc1(tcopflow->Nt,&tcopflow->nconineqcoup);CHKERRQ(ierr);
  ierr = (*tcopflow->solverops.setup)(tcopflow);CHKERRQ(ierr);
  ierr = PetscPrintf(tcopflow->comm->type,"TCOPFLOW: Setup completed\n");CHKERRQ(ierr);
  
  tcopflow->setupcalled = PETSC_TRUE;
  PetscFunctionReturn(0);
}


/*
  TCOPFLOWSolve - Solves the AC multi-period optimal power flow

  Input Parameters:
. tcopflow - the multi-period optimal power flow application object
*/
PetscErrorCode TCOPFLOWSolve(TCOPFLOW tcopflow)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if(!tcopflow->setupcalled) {
    ierr = TCOPFLOWSetUp(tcopflow);
  }

  /* Solve */
  ierr = (*tcopflow->solverops.solve)(tcopflow);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  TCOPFLOWGetObjective - Returns the objective function value

  Input Parameters:
+ TCOPFLOW - the TCOPFLOW object
- obj    - the objective function value

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode TCOPFLOWGetObjective(TCOPFLOW tcopflow,PetscReal *obj)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*tcopflow->solverops.getobjective)(tcopflow,obj);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  TCOPFLOWGetSolution - Returns the TCOPFLOW solution for a given time-step

  Input Parameters:
+ TCOPFLOW - the TCOPFLOW object
. tnum     - time-step (0 is the first time-step)
- X        - the tcopflow solution for the given time-step

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode TCOPFLOWGetSolution(TCOPFLOW tcopflow,PetscInt tnum,Vec *X)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*tcopflow->solverops.getsolution)(tcopflow,tnum,X);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  TCOPFLOWGetConstraints - Returns the TCOPFLOW constraints for a given time-step

  Input Parameters:
+ TCOPFLOW - the TCOPFLOW object
. tnum     - time-step (0 for the first time-step)
- G        - the tcopflow constraints

  Notes: Should be called after the optimization finishes.
         Equality constraints first followed by inequality constraints
*/
PetscErrorCode TCOPFLOWGetConstraints(TCOPFLOW tcopflow,PetscInt tnum,Vec *G)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*tcopflow->solverops.getconstraints)(tcopflow,tnum,G);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  TCOPFLOWGetConstraintMultipliers - Returns the TCOPFLOW constraint multipliers for a given time-step

  Input Parameters:
+ TCOPFLOW - the TCOPFLOW object
. tnum     - time-step (0 for the first time-step)
- G    - the tcopflow constraint lagrange multipliers

  Notes: Should be called after the optimization finishes.
    Equality constraint multipliers first followed by inequality constraint multipliers
*/
PetscErrorCode TCOPFLOWGetConstraintMultipliers(TCOPFLOW tcopflow,PetscInt tnum,Vec *Lambda)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*tcopflow->solverops.getconstraintmultipliers)(tcopflow,tnum,Lambda);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  TCOPFLOWGetConvergenceStatus - Did TCOPFLOW converge?

  Input Parameters:
+ TCOPFLOW - the TCOPFLOW object
- status - PETSC_TRUE if converged, PETSC_FALSE otherwise

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode TCOPFLOWGetConvergenceStatus(TCOPFLOW tcopflow,PetscBool *status)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*tcopflow->solverops.getconvergencestatus)(tcopflow,status);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*
  TCOPFLOWSolutionToPS - Updates the PS struct from TCOPFLOW solution

  Input Parameters:
. tcopflow - the TCOPFLOW object

  Notes: Updates the different fields in the PS struct from the TCOPFLOW solution
*/
PetscErrorCode TCOPFLOWSolutionToPS(TCOPFLOW tcopflow)
{
  PetscErrorCode     ierr;
  PetscInt           i;
  OPFLOW             opflow;

  PetscFunctionBegin;

  ierr = (*opflow->modelops.solutiontops)(opflow);

  tcopflow->solutiontops = PETSC_TRUE;
  PetscFunctionReturn(0);
}
