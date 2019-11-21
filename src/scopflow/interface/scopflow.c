#include <private/scopflowimpl.h>
#include <petsc/private/dmnetworkimpl.h>

/*
  SCOPFLOWCreate - Creates an security constrained optimal power flow application object

  Input Parameters
. mpicomm - The MPI communicator

  Output Parameters
. scopflowout - The security constrained optimal power flow application object
*/
PetscErrorCode SCOPFLOWCreate(MPI_Comm mpicomm, SCOPFLOW *scopflowout)
{
  PetscErrorCode ierr;
  SCOPFLOW         scopflow;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&scopflow);CHKERRQ(ierr);

  ierr = COMMCreate(mpicomm,&scopflow->comm);CHKERRQ(ierr);

  scopflow->Nconeq   = scopflow->nconeq  = -1;
  scopflow->Nconineq = scopflow->nconineq = -1;
  scopflow->Ncon     = scopflow->ncon     = -1;
  scopflow->Nx       = scopflow->nx       = -1;
  scopflow->Gi       = NULL;
  scopflow->Lambdai  = NULL;
  
  scopflow->obj_factor = 1.0;
  scopflow->obj = 0.0;

  scopflow->solver   = NULL;

  scopflow->nsolversregistered = 0;
  scopflow->SCOPFLOWSolverRegisterAllCalled = PETSC_FALSE;

  /* Register all solvers */
  ierr = SCOPFLOWSolverRegisterAll(scopflow);

  /* Run-time options */
  scopflow->iscoupling = PETSC_FALSE;
  scopflow->first_stage_gen_cost_only = PETSC_FALSE;
  scopflow->ignore_line_flow_constraints = PETSC_FALSE;
  scopflow->replicate_basecase = PETSC_FALSE;

  scopflow->ctgcfileset = PETSC_FALSE;
  scopflow->setupcalled = PETSC_FALSE;
  *scopflowout = scopflow;

  ierr = PetscPrintf(scopflow->comm->type,"SCOPFLOW: Application created\n");
  PetscFunctionReturn(0);
}

/*
  SCOPFLOWDestroy - Destroys the security constrained optimal power flow application object

  Input Parameter
. scopflow - The SCOPFLOW object to destroy
*/
PetscErrorCode SCOPFLOWDestroy(SCOPFLOW *scopflow)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = COMMDestroy(&(*scopflow)->comm);CHKERRQ(ierr);

  /* Solution vector */
  ierr = VecDestroy(&(*scopflow)->X);CHKERRQ(ierr);
  ierr = VecDestroy(&(*scopflow)->localX);CHKERRQ(ierr);
  ierr = VecDestroy(&(*scopflow)->gradobj);CHKERRQ(ierr);

  /* Lower and upper bounds on X */
  ierr = VecDestroy(&(*scopflow)->Xl);CHKERRQ(ierr);
  ierr = VecDestroy(&(*scopflow)->Xu);CHKERRQ(ierr);

  /* Constraints vector */
  ierr = VecDestroy(&(*scopflow)->G);CHKERRQ(ierr);
  ierr = VecDestroy(&(*scopflow)->Ge);CHKERRQ(ierr);
  ierr = VecDestroy(&(*scopflow)->Gelocal);CHKERRQ(ierr);
  ierr = VecDestroy(&(*scopflow)->Lambdae);CHKERRQ(ierr);
  ierr = VecDestroy(&(*scopflow)->Lambdaelocal);CHKERRQ(ierr);
  if((*scopflow)->nconineq) {
    ierr = VecDestroy(&(*scopflow)->Gi);CHKERRQ(ierr);
    ierr = VecDestroy(&(*scopflow)->Lambdai);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&(*scopflow)->Gl);CHKERRQ(ierr);
  ierr = VecDestroy(&(*scopflow)->Gu);CHKERRQ(ierr);
  ierr = VecDestroy(&(*scopflow)->Lambda);CHKERRQ(ierr);

  /* Jacobian of constraints */
  ierr = MatDestroy(&(*scopflow)->Jac);CHKERRQ(ierr);
  ierr = MatDestroy(&(*scopflow)->Jac_Ge);CHKERRQ(ierr);
  ierr = MatDestroy(&(*scopflow)->Jac_Gi);CHKERRQ(ierr);

  ierr = MatDestroy(&(*scopflow)->Hes);CHKERRQ(ierr);

  if((*scopflow)->solverops.destroy) {
    ierr = ((*scopflow)->solverops.destroy)(*scopflow);
  }

  ierr = PetscFree(*scopflow);CHKERRQ(ierr);
  *scopflow = 0;
  PetscFunctionReturn(0);
}

/*
  SCOPFLOWSetSolver - Sets the solver for SCOPFLOW

  Input Parameters:
+ scopflow - scopflow application object
- solvername - name of the solver
*/
PetscErrorCode SCOPFLOWSetSolver(SCOPFLOW scopflow,const char* solvername)
{
  PetscErrorCode ierr,(*r)(SCOPFLOW)=NULL;
  PetscInt       i;
  PetscFunctionBegin;
  PetscBool match;
  for(i=0;i < scopflow->nsolversregistered;i++) {
    ierr = PetscStrcmp(scopflow->SCOPFLOWSolverList[i].name,solvername,&match);CHKERRQ(ierr);
    if(match) {
      r = scopflow->SCOPFLOWSolverList[i].create;
      break;
    }
  }

  if(!r) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown type for SCOPFLOW Solver %s",solvername);

  /* Initialize (Null) the function pointers */
  scopflow->solverops.destroy = 0;
  scopflow->solverops.solve   = 0;
  scopflow->solverops.setup   = 0;

  /* Call the underlying implementation constructor */
  ierr = (*r)(scopflow);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  SCOPFLOWSetNetworkData - Sets and reads the network data

  Input Parameter
+  scopflow - The SCOPFLOW object
-  netfile - The name of the network file

  Notes: The input data must be in MATPOWER data format
*/
PetscErrorCode SCOPFLOWSetNetworkData(SCOPFLOW scopflow,const char netfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscMemcpy(scopflow->netfile,netfile,100*sizeof(char));CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  SCOPFLOWSetContingencyData - Sets the contingency data

  Input Parameter
+  scopflow - The SCOPFLOW object
-  ctgcfile - The name of the contingency list file

*/
PetscErrorCode SCOPFLOWSetContingencyData(SCOPFLOW scopflow,const char ctgcfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscMemcpy(scopflow->ctgcfile,ctgcfile,100*sizeof(char));CHKERRQ(ierr);

  scopflow->ctgcfileset = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/*
  SCOPFLOWSetUp - Sets up an security constrained optimal power flow application object

  Input Parameters:
. scopflow - the SCOPFLOW object

  Notes:
  This routine sets up the SCOPFLOW object and the underlying PS object. It
  also distributes the PS object when used in parallel.
*/
PetscErrorCode SCOPFLOWSetUp(SCOPFLOW scopflow)
{
  PetscErrorCode ierr;
  PetscBool      formulationset=PETSC_FALSE;
  PetscBool      solverset=PETSC_FALSE;
  char           solvername[32];

  PetscFunctionBegin;

  ierr = PetscOptionsGetString(NULL,NULL,"-scopflow_formulation",formulationname,32,&formulationset);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL,NULL,"-scopflow_solver",solvername,32,&solverset);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-scopflow_iscoupling",&scopflow->iscoupling,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-scopflow_first_stage_gen_cost_only",&scopflow->first_stage_gen_cost_only,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-scopflow_ignore_line_flow_constraints",&scopflow->ignore_line_flow_constraints,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-scopflow_Ns",&scopflow->Ns,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-scopflow_replicate_basecase",&scopflow->replicate_basecase,NULL);CHKERRQ(ierr);

  /* Set solver */
  if(solverset) {
    if(scopflow->solver) ierr = (*scopflow->solverops.destroy)(scopflow);
    ierr = SCOPFLOWSetSolver(scopflow,solvername);CHKERRQ(ierr);
    ierr = PetscPrintf(scopflow->comm->type,"SCOPFLOW: Using %s solver\n",solvername);CHKERRQ(ierr);
  } else {
    if(!scopflow->solver) {
      ierr = SCOPFLOWSetSolver(scopflow,SCOPFLOWSOLVER_IPOPT);CHKERRQ(ierr);
      ierr = PetscPrintf(scopflow->comm->type,"SCOPFLOW: Using %s solver\n",SCOPFLOWSOLVER_IPOPT);CHKERRQ(ierr); 
    }
  }

  ierr = (*scopflow->solverops.setup)(scopflow);CHKERRQ(ierr);
  ierr = PetscPrintf(scopflow->comm->type,"SCOPFLOW: Setup completed\n");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


/*
  SCOPFLOWSolve - Solves the AC security constrained optimal power flow

  Input Parameters:
. scopflow - the security constrained optimal power flow application object
*/
PetscErrorCode SCOPFLOWSolve(SCOPFLOW scopflow)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if(!scopflow->setupcalled) {
    ierr = SCOPFLOWSetUp(scopflow);
  }

  /* Solve */
  ierr = (*scopflow->solverops.solve)(scopflow);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

