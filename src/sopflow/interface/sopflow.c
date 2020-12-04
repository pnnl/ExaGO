#include <private/opflowimpl.h>
#include <private/sopflowimpl.h>

/*
  SOPFLOWCreate - Creates a stochastic optimal power flow application object

  Input Parameters
. mpicomm - The MPI communicator

  Output Parameters
. sopflowout - The stochastic optimal power flow application object
*/
PetscErrorCode SOPFLOWCreate(MPI_Comm mpicomm, SOPFLOW *sopflowout)
{
  PetscErrorCode ierr;
  SOPFLOW         sopflow;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&sopflow);CHKERRQ(ierr);

  ierr = COMMCreate(mpicomm,&sopflow->comm);CHKERRQ(ierr);

  sopflow->Nconeq   = sopflow->nconeq   = 0;
  sopflow->Nconineq = sopflow->nconineq = 0;
  sopflow->Ncon     = sopflow->ncon     = 0;
  sopflow->Nx       = sopflow->nx       = 0;
  sopflow->Ns       = sopflow->ns       = 0;
  sopflow->Gi       = NULL;
  sopflow->Lambdai  = NULL;
  
  sopflow->obj_factor = 1.0;
  sopflow->obj = 0.0;


  sopflow->makeup_power_source = 1;

  sopflow->solver   = NULL;
  sopflow->model    = NULL;

  sopflow->nmodelsregistered = 0;
  sopflow->SOPFLOWModelRegisterAllCalled = PETSC_FALSE;

  /* Register all models */
  ierr = SOPFLOWModelRegisterAll(sopflow);

  sopflow->nsolversregistered = 0;
  sopflow->SOPFLOWSolverRegisterAllCalled = PETSC_FALSE;

  /* Register all solvers */
  ierr = SOPFLOWSolverRegisterAll(sopflow);

  /* Run-time options */
  sopflow->iscoupling = PETSC_FALSE;

  sopflow->scenfileset = PETSC_FALSE;
  sopflow->scenunctype = NONE;
  sopflow->setupcalled = PETSC_FALSE;
  *sopflowout = sopflow;

  ierr = PetscPrintf(sopflow->comm->type,"SOPFLOW: Application created\n");
  PetscFunctionReturn(0);
}

/*
  SOPFLOWDestroy - Destroys the stochastic optimal power flow application object

  Input Parameter
. sopflow - The SOPFLOW object to destroy
*/
PetscErrorCode SOPFLOWDestroy(SOPFLOW *sopflow)
{
  PetscErrorCode ierr;
  PetscInt       s;

  PetscFunctionBegin;
  ierr = COMMDestroy(&(*sopflow)->comm);CHKERRQ(ierr);

  /* Solution vector */
  ierr = VecDestroy(&(*sopflow)->X);CHKERRQ(ierr);
  ierr = VecDestroy(&(*sopflow)->gradobj);CHKERRQ(ierr);

  /* Lower and upper bounds on X */
  ierr = VecDestroy(&(*sopflow)->Xl);CHKERRQ(ierr);
  ierr = VecDestroy(&(*sopflow)->Xu);CHKERRQ(ierr);

  /* Constraints vector */
  ierr = VecDestroy(&(*sopflow)->G);CHKERRQ(ierr);

  ierr = VecDestroy(&(*sopflow)->Gl);CHKERRQ(ierr);
  ierr = VecDestroy(&(*sopflow)->Gu);CHKERRQ(ierr);
  ierr = VecDestroy(&(*sopflow)->Lambda);CHKERRQ(ierr);

  /* Jacobian of constraints */
  ierr = MatDestroy(&(*sopflow)->Jac);CHKERRQ(ierr);
  ierr = MatDestroy(&(*sopflow)->Hes);CHKERRQ(ierr);

  if((*sopflow)->solverops.destroy) {
    ierr = ((*sopflow)->solverops.destroy)(*sopflow);
  }

  if((*sopflow)->modelops.destroy) {
    ierr = ((*sopflow)->modelops.destroy)(*sopflow);
  }

  /* Destroy OPFLOW objects */
  for(s=0; s < (*sopflow)->ns; s++) {
    ierr = OPFLOWDestroy(&(*sopflow)->opflows[s]);CHKERRQ(ierr);
  }

  ierr = PetscFree((*sopflow)->xstarti);CHKERRQ(ierr);
  ierr = PetscFree((*sopflow)->gstarti);CHKERRQ(ierr);
  ierr = PetscFree((*sopflow)->nxi);CHKERRQ(ierr);
  ierr = PetscFree((*sopflow)->ngi);CHKERRQ(ierr);
  ierr = PetscFree((*sopflow)->nconineqcoup);CHKERRQ(ierr);
  ierr = PetscFree((*sopflow)->opflows);CHKERRQ(ierr);
  ierr = PetscFree(*sopflow);CHKERRQ(ierr);
  //  *sopflow = 0;
  PetscFunctionReturn(0);
}

/*
  SOPFLOWSetModel - Sets the model for SOPFLOW

  Input Parameters:
+ sopflow - opflow application object
- modelname - name of the model
*/
PetscErrorCode SOPFLOWSetModel(SOPFLOW sopflow,const char* modelname)
{
  PetscErrorCode ierr,(*r)(SOPFLOW)=NULL;
  PetscInt       i;
  PetscFunctionBegin;
  PetscBool match;
  for(i=0;i < sopflow->nmodelsregistered;i++) {
    ierr = PetscStrcmp(sopflow->SOPFLOWModelList[i].name,modelname,&match);CHKERRQ(ierr);
    if(match) {
      r = sopflow->SOPFLOWModelList[i].create;
      break;
    }
  }

  if(!r) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown type for SOPFLOW Model %s",modelname);

  /* Null the function pointers */
  sopflow->modelops.destroy                        = 0;
  sopflow->modelops.setup                          = 0;
  sopflow->modelops.setnumvariablesandconstraints  = 0;
  sopflow->modelops.setvariablebounds              = 0;
  sopflow->modelops.setconstraintbounds            = 0;
  sopflow->modelops.setvariableandconstraintbounds = 0;
  sopflow->modelops.setinitialguess                = 0;
  sopflow->modelops.computeconstraints             = 0;
  sopflow->modelops.computejacobian                = 0;
  sopflow->modelops.computehessian                 = 0;
  sopflow->modelops.computeobjandgradient          = 0;
  sopflow->modelops.computeobjective               = 0;
  sopflow->modelops.computegradient                = 0;

  ierr = PetscStrcpy(sopflow->modelname,modelname);CHKERRQ(ierr);
  /* Call the underlying implementation constructor */
  ierr = (*r)(sopflow);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  SOPFLOWSetSolver - Sets the solver for SOPFLOW

  Input Parameters:
+ sopflow - sopflow application object
- solvername - name of the solver
*/
PetscErrorCode SOPFLOWSetSolver(SOPFLOW sopflow,const char* solvername)
{
  PetscErrorCode ierr,(*r)(SOPFLOW)=NULL;
  PetscInt       i;
  PetscFunctionBegin;
  PetscBool match;
  for(i=0;i < sopflow->nsolversregistered;i++) {
    ierr = PetscStrcmp(sopflow->SOPFLOWSolverList[i].name,solvername,&match);CHKERRQ(ierr);
    if(match) {
      r = sopflow->SOPFLOWSolverList[i].create;
      break;
    }
  }

  if(!r) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown type for SOPFLOW Solver %s",solvername);

  /* Initialize (Null) the function pointers */
  sopflow->solverops.destroy = 0;
  sopflow->solverops.solve   = 0;
  sopflow->solverops.setup   = 0;

  ierr = PetscStrcpy(sopflow->solvername,solvername);CHKERRQ(ierr);
  /* Call the underlying implementation constructor */
  ierr = (*r)(sopflow);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  SOPFLOWSetNetworkData - Sets and reads the network data

  Input Parameter
+  sopflow - The SOPFLOW object
-  netfile - The name of the network file

  Notes: The input data must be in MATPOWER data format
*/
PetscErrorCode SOPFLOWSetNetworkData(SOPFLOW sopflow,const char netfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscMemcpy(sopflow->netfile,netfile,100*sizeof(char));CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

extern PetscErrorCode SOPFLOWGetNumScenarios(SOPFLOW,ScenarioFileInputFormat,const char scenfile[],PetscInt*);

/*
  SOPFLOWSetUp - Sets up an stochastic optimal power flow application object

  Input Parameters:
. sopflow - the SOPFLOW object

  Notes:
  This routine sets up the SOPFLOW object and the underlying PS object. It
  also distributes the PS object when used in parallel.
*/
PetscErrorCode SOPFLOWSetUp(SOPFLOW sopflow)
{
  PetscErrorCode ierr;
  PetscBool      sopflowsolverset;
  char           opflowmodelname[32]="POWER_BALANCE_POLAR";
  char           sopflowsolvername[32]="IPOPT";
  PetscInt       c,i,j;
  char           sopflowmodelname[32]="GENRAMP";
  PS             ps;
  OPFLOW         opflow;

  PetscFunctionBegin;

  ierr = PetscOptionsBegin(sopflow->comm->type,NULL,"SOPFLOW options",NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-sopflow_model","SOPFLOW model type","",sopflowmodelname,sopflowmodelname,32,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-sopflow_solver","SOPFLOW solver type","",sopflowsolvername,sopflowsolvername,32,&sopflowsolverset);CHKERRQ(ierr);

  ierr = PetscOptionsBool("-sopflow_iscoupling","Include coupling between first stage and second stage","",sopflow->iscoupling,&sopflow->iscoupling,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-sopflow_Ns","Number of second stage scenarios","",sopflow->Ns,&sopflow->Ns,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-sopflow_mode","Operation mode:Preventive (0) or Corrective (1)","",sopflow->mode,&sopflow->mode,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-sopflow_makeup_power_source","Make-up power source for scenarios","",sopflow->makeup_power_source,&sopflow->makeup_power_source,NULL);CHKERRQ(ierr);
  
  PetscOptionsEnd();

  if(sopflow->Ns >= 0) sopflow->Ns += 1;
  if(sopflow->scenfileset) {
    if(sopflow->Ns == -1) { 
      ierr = SOPFLOWGetNumScenarios(sopflow,sopflow->scenfileformat,sopflow->scenfile,&sopflow->Ns);CHKERRQ(ierr);
      sopflow->Ns += 1;
    }
  } else {
    if(sopflow->Ns == -1) sopflow->Ns = 1;
  }

  int q = sopflow->Ns/sopflow->comm->size;
  int d = sopflow->Ns%sopflow->comm->size;
  if(d) {
    sopflow->ns = q + ((sopflow->comm->rank < d)?1:0); 
  } else {
    sopflow->ns = q;
  }
  ierr = MPI_Scan(&sopflow->ns,&sopflow->send,1,MPIU_INT,MPI_SUM,sopflow->comm->type);CHKERRQ(ierr);
  sopflow->sstart = sopflow->send - sopflow->ns;

  ierr = PetscPrintf(sopflow->comm->type,"SOPFLOW running with %d scenarios (base case + %d scenarios)\n",sopflow->Ns,sopflow->Ns-1);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF,"Rank %d has %d scenarios, range [%d -- %d]\n",sopflow->comm->rank,sopflow->ns,sopflow->sstart,sopflow->send);CHKERRQ(ierr);

  /* Set Model */
  ierr = SOPFLOWSetModel(sopflow,sopflowmodelname);CHKERRQ(ierr);

  /* Set solver */
  if(sopflowsolverset) {
    if(sopflow->solver) ierr = (*sopflow->solverops.destroy)(sopflow);
    ierr = SOPFLOWSetSolver(sopflow,sopflowsolvername);CHKERRQ(ierr);
    ierr = PetscPrintf(sopflow->comm->type,"SOPFLOW: Using %s solver\n",sopflowsolvername);CHKERRQ(ierr);
  } else {
    if(!sopflow->solver) {
      ierr = SOPFLOWSetSolver(sopflow,SOPFLOWSOLVER_IPOPT);CHKERRQ(ierr);
      ierr = PetscPrintf(sopflow->comm->type,"SOPFLOW: Using %s solver\n",SOPFLOWSOLVER_IPOPT);CHKERRQ(ierr); 
    }
  }

  /* Create OPFLOW objects */
  ierr = PetscCalloc1(sopflow->ns,&sopflow->opflows);CHKERRQ(ierr);
  for(c=0; c < sopflow->ns; c++) {
    ierr = OPFLOWCreate(PETSC_COMM_SELF,&sopflow->opflows[c]);CHKERRQ(ierr);
    ierr = OPFLOWSetModel(sopflow->opflows[c],OPFLOWMODEL_PBPOL);CHKERRQ(ierr);
    //    ierr = OPFLOWSetSolver(sopflow->opflows[c],opflowsolvername);CHKERRQ(ierr);

    ierr = OPFLOWReadMatPowerData(sopflow->opflows[c],sopflow->netfile);CHKERRQ(ierr);
    /* Set up the PS object for opflow */
    ps = sopflow->opflows[c]->ps;
    ierr = PSSetUp(ps);CHKERRQ(ierr);

    /* Set up OPFLOW object */
    //    ierr = OPFLOWGenbusVoltageFixed(sopflow->opflows[c],PETSC_TRUE);CHKERRQ(ierr);
    //    if(sopflow->sstart+c > 0) sopflow->opflows[c]->obj_gencost = PETSC_FALSE; /* No gen. cost minimization for second stage */
    ierr = OPFLOWSetUp(sopflow->opflows[c]);CHKERRQ(ierr);
  }
  
  if(sopflow->scenfileset) {
    ierr = SOPFLOWReadScenarioData(sopflow,sopflow->scenfileformat,sopflow->scenfile);CHKERRQ(ierr);
  }

  ierr = PetscCalloc1(sopflow->ns,&sopflow->nconineqcoup);CHKERRQ(ierr);
  ierr = PetscCalloc1(sopflow->ns,&sopflow->nxi);CHKERRQ(ierr);
  ierr = PetscCalloc1(sopflow->ns,&sopflow->ngi);CHKERRQ(ierr);
  ierr = PetscCalloc1(sopflow->ns,&sopflow->xstarti);CHKERRQ(ierr);
  ierr = PetscCalloc1(sopflow->ns,&sopflow->gstarti);CHKERRQ(ierr);

  /* Set number of variables and constraints */
  ierr = (*sopflow->modelops.setnumvariablesandconstraints)(sopflow,sopflow->nxi,sopflow->ngi,sopflow->nconineqcoup);

  sopflow->nx = sopflow->nxi[0];
  sopflow->ncon = sopflow->ngi[0];
  sopflow->nconcoup = sopflow->nconineqcoup[0];
  opflow = sopflow->opflows[0];
  sopflow->nconeq = opflow->nconeq;
  sopflow->nconineq = opflow->nconineq;

  for(i=1; i < sopflow->ns; i++) {
    sopflow->xstarti[i] = sopflow->xstarti[i-1] + sopflow->nxi[i-1];
    sopflow->gstarti[i] = sopflow->gstarti[i-1] + sopflow->ngi[i-1];
    sopflow->nx += sopflow->nxi[i];
    sopflow->ncon += sopflow->ngi[i];
    sopflow->nconcoup += sopflow->nconineqcoup[i];
    opflow = sopflow->opflows[i];
    sopflow->nconeq   += opflow->nconeq;
    sopflow->nconineq += opflow->nconineq;
  }

  /* Create vector X */
  ierr = VecCreate(sopflow->comm->type,&sopflow->X);CHKERRQ(ierr);
  ierr = VecSetSizes(sopflow->X,sopflow->nx,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(sopflow->X);CHKERRQ(ierr);
  ierr = VecGetSize(sopflow->X,&sopflow->Nx);CHKERRQ(ierr);

  ierr = VecDuplicate(sopflow->X,&sopflow->Xl);CHKERRQ(ierr);
  ierr = VecDuplicate(sopflow->X,&sopflow->Xu);CHKERRQ(ierr);
  ierr = VecDuplicate(sopflow->X,&sopflow->gradobj);CHKERRQ(ierr);

  /* vector for constraints */
  ierr = VecCreate(sopflow->comm->type,&sopflow->G);CHKERRQ(ierr);
  ierr = VecSetSizes(sopflow->G,sopflow->ncon,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(sopflow->G);CHKERRQ(ierr);
  ierr = VecGetSize(sopflow->G,&sopflow->Ncon);CHKERRQ(ierr);

  /* Constraint bounds vectors  */
  ierr = VecDuplicate(sopflow->G,&sopflow->Gl);CHKERRQ(ierr);
  ierr = VecDuplicate(sopflow->G,&sopflow->Gu);CHKERRQ(ierr);

  /* Constraint Jacobian */
  ierr = MatCreate(sopflow->comm->type,&sopflow->Jac);CHKERRQ(ierr);
  ierr = MatSetSizes(sopflow->Jac,sopflow->ncon,sopflow->nx,sopflow->Ncon,sopflow->Nx);CHKERRQ(ierr);
  ierr = MatSetUp(sopflow->Jac);CHKERRQ(ierr);
  ierr = MatSetFromOptions(sopflow->Jac);CHKERRQ(ierr);

  /* Hessian */
  ierr = MatCreate(sopflow->comm->type,&sopflow->Hes);CHKERRQ(ierr);
  ierr = MatSetSizes(sopflow->Hes,sopflow->nx,sopflow->nx,sopflow->Nx,sopflow->Nx);CHKERRQ(ierr);
  ierr = MatSetUp(sopflow->Hes);CHKERRQ(ierr);
  ierr = MatSetFromOptions(sopflow->Hes);CHKERRQ(ierr);

  /* Lagrangian multipliers */
  ierr = VecDuplicate(sopflow->G,&sopflow->Lambda);CHKERRQ(ierr);

  ierr = (*sopflow->solverops.setup)(sopflow);CHKERRQ(ierr);
  ierr = PetscPrintf(sopflow->comm->type,"SOPFLOW: Setup completed\n");CHKERRQ(ierr);
  
  sopflow->setupcalled = PETSC_TRUE;
  PetscFunctionReturn(0);
}


/*
  SOPFLOWSolve - Solves the AC stochastic optimal power flow

  Input Parameters:
. sopflow - the stochastic optimal power flow application object
*/
PetscErrorCode SOPFLOWSolve(SOPFLOW sopflow)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if(!sopflow->setupcalled) {
    ierr = SOPFLOWSetUp(sopflow);
  }

  /* Set bounds on variables */
  if(sopflow->modelops.setvariablebounds) {
    ierr = (*sopflow->modelops.setvariablebounds)(sopflow,sopflow->Xl,sopflow->Xu);CHKERRQ(ierr);
  }
  
  /* Set bounds on constraints */
  if(sopflow->modelops.setconstraintbounds) {
    ierr = (*sopflow->modelops.setconstraintbounds)(sopflow,sopflow->Gl,sopflow->Gu);CHKERRQ(ierr);
  }
  
  /* Set initial guess */
  if(sopflow->modelops.setinitialguess) {
    ierr = (*sopflow->modelops.setinitialguess)(sopflow,sopflow->X);CHKERRQ(ierr);
  }

  ierr = VecSet(sopflow->Lambda,1.0);CHKERRQ(ierr);

  /* Solve */
  ierr = (*sopflow->solverops.solve)(sopflow);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  SOPFLOWGetObjective - Returns the objective function value

  Input Parameters:
+ SOPFLOW - the SOPFLOW object
- obj    - the objective function value

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode SOPFLOWGetObjective(SOPFLOW sopflow,PetscReal *obj)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*sopflow->solverops.getobjective)(sopflow,obj);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SOPFLOWGetBaseSolution - Returns the SOPFLOW solution for a given scenario

  Input Parameters:
+ SOPFLOW - the SOPFLOW object
. scennum  - Scenario number (0 for base)
- X        - the sopflow solution for the given scenario

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode SOPFLOWGetSolution(SOPFLOW sopflow,PetscInt scennum,Vec *X)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*sopflow->solverops.getsolution)(sopflow,scennum,X);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SOPFLOWGetConstraints - Returns the SOPFLOW constraints for a given scenario

  Input Parameters:
+ SOPFLOW - the SOPFLOW object
. scennum  - Scenario number (0 for base)
- G    - the sopflow constraints

  Notes: Should be called after the optimization finishes.
         Equality constraints first followed by inequality constraints
*/
PetscErrorCode SOPFLOWGetConstraints(SOPFLOW sopflow,PetscInt scennum,Vec *G)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*sopflow->solverops.getconstraints)(sopflow,scennum,G);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SOPFLOWGetConstraintMultipliers - Returns the SOPFLOW constraint multipliers for a given scenario

  Input Parameters:
+ SOPFLOW - the SOPFLOW object
. scennum  - Scenario number (0 for base)
- G    - the sopflow constraint lagrange multipliers

  Notes: Should be called after the optimization finishes.
    Equality constraint multipliers first followed by inequality constraint multipliers
*/
PetscErrorCode SOPFLOWGetConstraintMultipliers(SOPFLOW sopflow,PetscInt scennum,Vec *Lambda)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*sopflow->solverops.getconstraintmultipliers)(sopflow,scennum,Lambda);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SOPFLOWGetConvergenceStatus - Did SOPFLOW converge?

  Input Parameters:
+ SOPFLOW - the SOPFLOW object
- status - PETSC_TRUE if converged, PETSC_FALSE otherwise

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode SOPFLOWGetConvergenceStatus(SOPFLOW sopflow,PetscBool *status)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*sopflow->solverops.getconvergencestatus)(sopflow,status);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
