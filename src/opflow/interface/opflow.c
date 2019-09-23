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

  /*  if((*opflow)->solverops.destroy) {
    ierr = ((*opflow)->solverops.destroy)(*opflow);
  }
  */

  //  if((*opflow)->formops.destroy) {
  //    ierr = (


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

