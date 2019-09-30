#include <private/opflowimpl.h>
#include <petsc/private/dmnetworkimpl.h>

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
  ierr = PSSetApplication(opflow->ps,opflow,APP_ACOPF);CHKERRQ(ierr);

  opflow->Nconeq   = opflow->nconeq  = -1;
  opflow->Nconineq = opflow->nconineq = -1;
  opflow->Ncon     = opflow->ncon     = -1;
  opflow->Nx       = opflow->nx       = -1;
  
  opflow->obj_factor = 1.0;
  opflow->obj = 0.0;

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
  ierr = VecDestroy(&(*opflow)->gradobj);CHKERRQ(ierr);

  /* Lower and upper bounds on X */
  ierr = VecDestroy(&(*opflow)->Xl);CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Xu);CHKERRQ(ierr);

  /* Constraints vector */
  ierr = VecDestroy(&(*opflow)->G);CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Ge);CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Gi);CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Gl);CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Gu);CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Lambda);CHKERRQ(ierr);

  /* Jacobian of constraints */
  ierr = MatDestroy(&(*opflow)->Jac);CHKERRQ(ierr);
  ierr = MatDestroy(&(*opflow)->Jac_Ge);CHKERRQ(ierr);
  ierr = MatDestroy(&(*opflow)->Jac_Gi);CHKERRQ(ierr);

  ierr = MatDestroy(&(*opflow)->Hes);CHKERRQ(ierr);
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
  opflow->formops.setnumvariables                = 0;
  opflow->formops.setnumconstraints              = 0;
  opflow->formops.setvariablebounds              = 0;
  opflow->formops.setconstraintbounds            = 0;
  opflow->formops.setvariableandconstraintbounds = 0;
  opflow->formops.setinitialguess                = 0;
  opflow->formops.computeequalityconstraints     = 0;
  opflow->formops.computeinequalityconstraints   = 0;
  opflow->formops.computeconstraints             = 0;
  opflow->formops.computeequalityconstraintjacobian = 0;
  opflow->formops.computeinequalityconstraintjacobian = 0;
  opflow->formops.computehessian                 = 0;
  opflow->formops.computeobjandgradient          = 0;
  opflow->formops.computeobjective               = 0;
  opflow->formops.computegradient                = 0;
  opflow->formops.computejacobian                = 0;

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
   OPFLOWSetNumConstraints - Sets the number of constraints for the OPFLOW problem

   Input Parameters:
.  OPFLOW - the opflow application object

   Output Parameters:
+  nconeq   - number of equality constraints
-  nconineq - number of inequality constraints

*/
PetscErrorCode OPFLOWSetNumConstraints(OPFLOW opflow,PetscInt *nconeq,PetscInt *nconineq)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = (*opflow->formops.setnumconstraints)(opflow,nconeq,nconineq);
  PetscFunctionReturn(0);
}
	  
/* 
   OPFLOWSetNumVariables - Sets the number of variables for the OPFLOW problem

   Input Parameters:
. opflow - OPFLOW application object

   Output Parameters:
+ busnvararray    - array of number of variables at each bus
. branchnvararray - array of number of variables at each branch 
- nx              - total number of variables
*/
PetscErrorCode OPFLOWSetNumVariables(OPFLOW opflow,PetscInt *busnvararray,PetscInt *branchnvararray,PetscInt *nx)
{
  PetscErrorCode ierr;
  PetscSection   varsection,oldvarsection;
  PetscInt       i,vStart,vEnd,eStart,eEnd;
  DM             networkdm=opflow->ps->networkdm;
  PS             ps=opflow->ps;
  DM             plexdm;
  PetscSection   globalvarsection;

  PetscFunctionBegin;

  ierr = PetscSectionCreate(opflow->comm->type,&varsection);CHKERRQ(ierr);

  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
  ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);

  ierr = PetscSectionSetChart(varsection,eStart,vEnd);CHKERRQ(ierr);

  ierr = (*opflow->formops.setnumvariables)(opflow,busnvararray,branchnvararray,nx);

  for(i=0; i < ps->nbranch; i++) {
    ierr = PetscSectionSetDof(varsection,eStart+i,branchnvararray[i]);
  }

  for(i=0; i < ps->nbus; i++) {
    ierr = PetscSectionSetDof(varsection,vStart+i,busnvararray[i]);CHKERRQ(ierr);
  }
  ierr = PetscSectionSetUp(varsection);CHKERRQ(ierr);

  ierr = DMNetworkGetPlex(networkdm,&plexdm);CHKERRQ(ierr);
  ierr = DMGetSection(plexdm,&oldvarsection);CHKERRQ(ierr);
  ierr = PetscSectionDestroy(&oldvarsection);CHKERRQ(ierr);
  /* This is hack to null the networkdm->DofSection pointer otherwise
     DMDestroy_Network throws an error. The issue is that though the
     section is destroyed, there is a dangling pointer for networkdm->DofSection.
     This hack needs to be fixed correctly via changes to DMNetwork
  */
  DM_Network *networkdmdata = (DM_Network*)(networkdm->data);
  networkdmdata->DofSection = 0;
  

  /* Set the section with number of variables */
  ierr = DMSetSection(plexdm,varsection);CHKERRQ(ierr);
  ierr = DMGetGlobalSection(plexdm,&globalvarsection);CHKERRQ(ierr);
  networkdmdata->DofSection = varsection;

  /* Update starting locations for variables at each bus */
  for(i=0; i < ps->nbus; i++) {
    ierr = DMNetworkGetVariableOffset(networkdm,vStart+i,&ps->bus[i].startloc);CHKERRQ(ierr);
    ierr = DMNetworkGetVariableGlobalOffset(networkdm,vStart+i,&ps->bus[i].startlocglob);CHKERRQ(ierr);
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
  PS             ps=opflow->ps;
  PetscBool      formulationset=PETSC_FALSE;
  PetscBool      solverset=PETSC_FALSE;
  char           formulationname[32],solvername[32];
  PetscInt       *busnvararray,*branchnvararray;

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

  /* Set up underlying PS object */
  ierr = PSSetUp(ps);CHKERRQ(ierr);

  ierr = PetscCalloc1(ps->nbus,&busnvararray);CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->nbranch,&branchnvararray);CHKERRQ(ierr);

  /* Set up number of variables for branches and buses */
  ierr = OPFLOWSetNumVariables(opflow,busnvararray,branchnvararray,&opflow->nx);CHKERRQ(ierr);

  /* Set number of constraints */
  ierr = OPFLOWSetNumConstraints(opflow,&opflow->nconeq,&opflow->nconineq);CHKERRQ(ierr);
  opflow->ncon = opflow->nconeq + opflow->nconineq;

  ierr = PetscFree(busnvararray);CHKERRQ(ierr); 
  ierr = PetscFree(branchnvararray);CHKERRQ(ierr);

  /* Set vertex local to global ordering */
  ierr = DMNetworkSetVertexLocalToGlobalOrdering(opflow->ps->networkdm);CHKERRQ(ierr);

  /* Create solution vector and upper/lower bounds */
  ierr = PSCreateGlobalVector(opflow->ps,&opflow->X);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->gradobj);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->Xl);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->Xu);CHKERRQ(ierr);

  /* Vectors for constraints */
  ierr = VecCreate(ps->comm->type,&opflow->G);CHKERRQ(ierr);
  ierr = VecSetSizes(opflow->G,PETSC_DECIDE,opflow->ncon);CHKERRQ(ierr);
  ierr = VecSetFromOptions(opflow->G);CHKERRQ(ierr);
  
  ierr = VecDuplicate(opflow->G,&opflow->Gl);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->G,&opflow->Gu);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->G,&opflow->Lambda);CHKERRQ(ierr);

  ierr = VecCreate(ps->comm->type,&opflow->Ge);CHKERRQ(ierr);
  ierr = VecSetSizes(opflow->Ge,PETSC_DECIDE,opflow->nconeq);CHKERRQ(ierr);
  ierr = VecSetFromOptions(opflow->Ge);CHKERRQ(ierr);

  ierr = VecCreate(ps->comm->type,&opflow->Gi);CHKERRQ(ierr);
  ierr = VecSetSizes(opflow->Gi,PETSC_DECIDE,opflow->ncon);CHKERRQ(ierr);
  ierr = VecSetFromOptions(opflow->Gi);CHKERRQ(ierr);

  /* Create equality and inequality constraint Jacobian matrices */
  ierr = MatCreate(opflow->comm->type,&opflow->Jac_Ge);CHKERRQ(ierr);
  ierr = MatSetSizes(opflow->Jac_Ge,opflow->nconeq,opflow->nx,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = MatSetUp(opflow->Jac_Ge);CHKERRQ(ierr);
  ierr = MatSetFromOptions(opflow->Jac_Ge);CHKERRQ(ierr);

  if(opflow->nconineq) {
    ierr = MatCreate(opflow->comm->type,&opflow->Jac_Gi);CHKERRQ(ierr);
    ierr = MatSetSizes(opflow->Jac_Gi,opflow->nconineq,opflow->nx,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
    ierr = MatSetUp(opflow->Jac_Gi);CHKERRQ(ierr);
    ierr = MatSetFromOptions(opflow->Jac_Gi);CHKERRQ(ierr);
  }

  /* Create Hessian */
  ierr = PSCreateMatrix(opflow->ps,&opflow->Hes);CHKERRQ(ierr);

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
  /* Set bounds on variables */
  if(opflow->formops.setvariablebounds) {
    ierr = (*opflow->formops.setvariablebounds)(opflow,opflow->Xl,opflow->Xu);CHKERRQ(ierr);
  }

  /* Set bounds on constraints */
  if(opflow->formops.setconstraintbounds) {
    ierr = (*opflow->formops.setconstraintbounds)(opflow,opflow->Gl,opflow->Gu);CHKERRQ(ierr);
  }

  /* Set bounds on variables and constraints */
  if(opflow->formops.setvariableandconstraintbounds) {
    ierr = (*opflow->formops.setvariableandconstraintbounds)(opflow,opflow->Xl,opflow->Xu,opflow->Gl,opflow->Gu);CHKERRQ(ierr);
  }

  /* Set initial guess */
  if(opflow->formops.setinitialguess) {
    ierr = (*opflow->formops.setinitialguess)(opflow,opflow->X);CHKERRQ(ierr);
    ierr = VecSet(opflow->Lambda,1.0);CHKERRQ(ierr);
  }

  /* Solve */
  ierr = (*opflow->solverops.solve)(opflow);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

