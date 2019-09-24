#include <private/opflowimpl.h>

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

  opflow->Nconeq   = opflow->nconeq  = -1;
  opflow->Nconineq = opflow->nconineq = -1;
  opflow->Ncon     = opflow->ncon     = -1;
  opflow->Nx       = opflow->nx       = -1;

  opflow->solver   = NULL;
  opflow->formulation = NULL;

  opflow->nformulationsregistered = opflow->nsolversregistered = 0;
  opflow->OPFLOWFormulationRegisterAllCalled = opflow->OPFLOWSolverRegisterAllCalled = PETSC_FALSE;

  /* Register all formulations */
  ierr = OPFLOWFormulationRegisterAll(opflow);

  /* Register all solvers */
  ierr = OPFLOWSolverRegisterAll(opflow);

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

  if((*opflow)->solverops.destroy) {
    ierr = ((*opflow)->solverops.destroy)(*opflow);
  }

  if((*opflow)->formops.destroy) {
    ierr = ((*opflow)->formops.destroy)(*opflow);
  }


  ierr = PetscFree(*opflow);CHKERRQ(ierr);
  *opflow = 0;
  PetscFunctionReturn(0);
}

/*
  OPFLOWSetSolver - Sets the solver for OPFLOW

  Input Parameters:
+ opflow - opflow application object
- solvername - name of the solver
*/
PetscErrorCode OPFLOWSetSolver(OPFLOW opflow,const char* solvername)
{
  PetscErrorCode ierr,(*r)(OPFLOW)=NULL;
  PetscInt       i;
  PetscFunctionBegin;
  PetscBool match;
  for(i=0;i < opflow->nsolversregistered;i++) {
    ierr = PetscStrcmp(opflow->OPFLOWSolverList[i].name,solvername,&match);CHKERRQ(ierr);
    if(match) {
      r = opflow->OPFLOWSolverList[i].create;
      break;
    }
  }

  if(!r) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown type for OPFLOW Solver %s",solvername);

  /* Initialize (Null) the function pointers */
  opflow->solverops.destroy = 0;
  opflow->solverops.solve   = 0;
  opflow->solverops.setup   = 0;

  /* Call the underlying implementation constructor */
  ierr = (*r)(opflow);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  OPFLOWSetFormulation - Sets the formulation for OPFLOW

  Input Parameters:
+ opflow - opflow application object
- solvername - name of the formulation
*/
PetscErrorCode OPFLOWSetFormulation(OPFLOW opflow,const char* formulationname)
{
  PetscErrorCode ierr,(*r)(OPFLOW)=NULL;
  PetscInt       i;
  PetscFunctionBegin;
  PetscBool match;
  for(i=0;i < opflow->nformulationsregistered;i++) {
    ierr = PetscStrcmp(opflow->OPFLOWFormulationList[i].name,formulationname,&match);CHKERRQ(ierr);
    if(match) {
      r = opflow->OPFLOWFormulationList[i].create;
      break;
    }
  }

  if(!r) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown type for OPFLOW Formulation %s",formulationname);

  /* Null the function pointers */
  opflow->formops.destroy                        = 0;
  opflow->formops.setvariablebounds              = 0;
  opflow->formops.setconstraintbounds            = 0;
  opflow->formops.setvariableandconstraintbounds = 0;
  opflow->formops.setinitialguess                = 0;
  opflow->formops.computeequalityconstraints     = 0;
  opflow->formops.computeinequalityconstraints   = 0;
  opflow->formops.computeconstraints             = 0;
  opflow->formops.computeobjandgradient          = 0;
  opflow->formops.computeobjective               = 0;
  opflow->formops.computegradient                = 0;
  opflow->formops.computejacobian                = 0;
  opflow->formops.computeequalityconstraintjacobian = 0;
  opflow->formops.computeinequalityconstraintjacobian = 0;
  opflow->formops.computehessian                 = 0;

  /* Call the underlying implementation constructor */
  ierr = (*r)(opflow);CHKERRQ(ierr);

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
  PetscBool      formulationset=PETSC_FALSE;
  PetscBool      solverset=PETSC_FALSE;
  char           formulationname[32],solvername[32];

  PetscFunctionBegin;

  ierr = PetscOptionsGetString(NULL,NULL,"-opflow_formulation",formulationname,32,&formulationset);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL,NULL,"-opflow_solver",solvername,32,&solverset);CHKERRQ(ierr);

  /* Set formulation */
  if(formulationset) {
    if(opflow->formulation) ierr = (*opflow->formops.destroy)(opflow);
    ierr = OPFLOWSetFormulation(opflow,formulationname);CHKERRQ(ierr);
  } else {
    if(!opflow->formulation) {
      ierr = OPFLOWSetFormulation(opflow,OPFLOWFORMULATION_PBPOL);CHKERRQ(ierr);
    }
  }

  /* Set solver */
  if(solverset) {
    if(opflow->solver) ierr = (*opflow->solverops.destroy)(opflow);
    ierr = OPFLOWSetSolver(opflow,solvername);CHKERRQ(ierr);
  } else {
    if(!opflow->solver) {
      ierr = OPFLOWSetSolver(opflow,OPFLOWSOLVER_IPOPT);CHKERRQ(ierr); 
    }
  }

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

  PetscFunctionBegin;

  if(!opflow->setupcalled) {
    ierr = OPFLOWSetUp(opflow);
  }

  PetscFunctionReturn(0);
}

