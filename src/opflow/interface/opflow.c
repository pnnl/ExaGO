#include <exago_config.h>
#include <private/pflowimpl.h>
#include <private/opflowimpl.h>
#include <petsc/private/dmnetworkimpl.h>

const char *const OPFLOWInitializationTypes[] = {"MIDPOINT","FROMFILE","ACPF","FLATSTART","OPFLOWInitializationType","OPFLOWINIT_",NULL};

const char *const OPFLOWObjectiveTypes[] = {"MIN_GEN_COST","MIN_GENSETPOINT_DEVIATION","OPFLOWObjectiveType","",NULL};

const char *const OPFLOWGenBusVoltageTypes[] = {"VARIABLE_WITHIN_BOUNDS","FIXED_WITHIN_QBOUNDS","FIXED_AT_SETPOINT","OPFLOWGenBusVoltageType","",NULL};

void swap_dm(DM *dm1, DM *dm2)
{
  DM temp = *dm1;
  *dm1 = *dm2;
  *dm2 = temp;
}

/* 
   OPFLOWSetGenBusVoltageType - Sets the voltage control mode for generator buses
   
   Input Parameters:
+  opflow - the opflow object
-  vtype  - voltage control type

   Command-line option: -opflow_genbusvoltage
   Notes: Must be called before OPFLOWSetUp()
*/
PetscErrorCode OPFLOWSetGenBusVoltageType(OPFLOW opflow,OPFLOWGenBusVoltageType vtype)
{
  PetscFunctionBegin;
  opflow->genbusvoltagetype = vtype;
  PetscFunctionReturn(0);
}

/*
  OPFLOWHasGenSetPoint - Use gen. set point in the OPFLOW formulation

  Input Parameters:
+ opflow - OPFLOW object
- hassetpoint - Use set-point?

  Notes: Should be called before OPFLOWSetUp
*/
PetscErrorCode OPFLOWHasGenSetPoint(OPFLOW opflow,PetscBool hassetpoint)
{
  PetscFunctionBegin;
  opflow->has_gensetpoint = hassetpoint;
  PetscFunctionReturn(0);
}

/*
  OPFLOWUseAGC - Uses AGC to proportionally redispatch the generators

  Input Parameters:
+ opflow - OPFLOW object
- useagc - Use AGC?

  Notes: Should be called before OPFLOWSetUp
*/
PetscErrorCode OPFLOWUseAGC(OPFLOW opflow,PetscBool useagc)
{
  PetscFunctionBegin;
  opflow->use_agc = useagc;
  if(useagc) opflow->has_gensetpoint = PETSC_TRUE;
  PetscFunctionReturn(0);
}


/*
  OPFLOWGetSizes - Returns the size of solution vector, and constraints

  Input Parameters:
- opflow - the OPFLOW object

  Output Parameters:
+ nx - size of X vector
. nconeq - size of equality constraints vector
- nconineq - size of inequality constraints vector

  Notes:
  Should be called after OPFLOWSetUp()
*/
PetscErrorCode OPFLOWGetSizes(OPFLOW opflow,int *nx,int *nconeq,int *nconineq)
{
  PetscFunctionBegin;
  *nx = opflow->nx;
  *nconeq = opflow->nconeq;
  *nconineq = opflow->nconineq;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWGetVariableOrdering(OPFLOW opflow,int **ordering)
{
  PetscFunctionBegin;
  *ordering = opflow->idxn2sd_map;
  PetscFunctionReturn(0);
}
  

/*
  OPFLOWSetUpInitPflow - Sets up the power flow solver object used in obtaining initial conditions.

  Note: This power flow solver object shares the underlying power system object and all the data. It
  is created at the end of OPFLOWSetUp. So, it already has the distributed ps object. The only
  difference is that it uses a different PetscSection for storing the degrees of freedom that
  gets associated with the dmnetwork (and plex)
*/
PetscErrorCode OPFLOWSetUpInitPflow(OPFLOW opflow)
{
  PetscErrorCode ierr;
  PetscInt       vStart,vEnd,eStart,eEnd;
  DM             networkdm,plexdm;
  PetscInt       i;

  PetscFunctionBegin;

  networkdm = opflow->ps->networkdm;

  ierr = PFLOWCreate(opflow->comm->type,&opflow->initpflow);CHKERRQ(ierr);

  /* PFLOW creates a new PS object. Destroy it so that we can associate the
     PS from OPFLOW with initpflow
  */
  ierr = PSDestroy(&opflow->initpflow->ps);CHKERRQ(ierr);

  opflow->initpflow->ps = opflow->ps;
  /* Increase the reference count for opflow->ps */
  ierr = PSIncreaseReferenceCount(opflow->ps);CHKERRQ(ierr);

  ierr = PetscSectionCreate(opflow->comm->type,&opflow->initpflowpsection);CHKERRQ(ierr);

  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
  ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);

  ierr = PetscSectionSetChart(opflow->initpflowpsection,eStart,vEnd);CHKERRQ(ierr);

  for(i=vStart; i < vEnd; i++) {
    /* Two variables at each bus/vertex */
    ierr = PetscSectionSetDof(opflow->initpflowpsection,i,2);CHKERRQ(ierr);
  }
  ierr = PetscSectionSetUp(opflow->initpflowpsection);CHKERRQ(ierr);

  /* Clone DM to be used with initpflow */
  ierr = DMClone(networkdm,&opflow->initpflowdm);CHKERRQ(ierr);

  /* Set initpflowdm in opflow->ps->networkdm, the previous networkdm get stashed in opflow->initpflowdm */
  swap_dm(&opflow->ps->networkdm,&opflow->initpflowdm);
  networkdm = opflow->ps->networkdm;

  /* Get the plex dm */
  ierr = DMNetworkGetPlex(networkdm,&plexdm);CHKERRQ(ierr);

  /* Get default sections associated with this plex */
  ierr = DMGetSection(plexdm,&opflow->defaultsection);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)opflow->defaultsection);CHKERRQ(ierr);

  ierr = DMGetGlobalSection(plexdm,&opflow->defaultglobalsection);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)opflow->defaultglobalsection);CHKERRQ(ierr);

  /* Set the new section created for initial power flow */
  ierr = DMSetSection(plexdm,opflow->initpflowpsection);CHKERRQ(ierr);
  ierr = DMGetGlobalSection(plexdm,&opflow->initpflowpglobsection);CHKERRQ(ierr);

  /* Set up PFLOW object. Note pflow->ps will not be set up again as it has
     been already set up by opflow
  */
  ierr = PFLOWSetUp(opflow->initpflow);CHKERRQ(ierr);

  /*
     Update the ref. counts for init pflow sections so that they do not
     get destroyed when DMSetDefaultSection is called
  */
  // ierr = PetscObjectReference((PetscObject)opflow->initpflowpsection);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)opflow->initpflowpglobsection);CHKERRQ(ierr);
  
  /* Reset the sections */
  ierr = DMSetSection(plexdm,opflow->defaultsection);CHKERRQ(ierr);
  ierr = DMSetGlobalSection(plexdm,opflow->defaultglobalsection);CHKERRQ(ierr);

  /* Reset dm */
  swap_dm(&opflow->ps->networkdm,&opflow->initpflowdm);

  PetscFunctionReturn(0);
}

/*
  OPFLOWComputePrePflow - Computes the steady-state power flow for initializing optimal power flow

  Input Parameters:
. opflow - the OPFLOW object
*/
PetscErrorCode OPFLOWComputePrePflow(OPFLOW opflow,PetscBool *converged)
{
  PetscErrorCode ierr;
  DM             plexdm;
  PFLOW          initpflow=opflow->initpflow;
  PS             ps=opflow->ps;

  PetscFunctionBegin;
  /* Set initpflowdm for solving the power flow */
  swap_dm(&opflow->ps->networkdm,&opflow->initpflowdm);

  ierr = DMNetworkGetPlex(opflow->ps->networkdm,&plexdm);CHKERRQ(ierr);

  /* Increase the ref. counts for the default sections so that they do not get
     destroyed when DMSetDefaultXXX is called with the initpflowxxx sections
  */
  ierr = PetscObjectReference((PetscObject)opflow->defaultsection);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)opflow->defaultglobalsection);CHKERRQ(ierr);

  /* Set the new section created for initial power flow. */
  ierr = DMSetSection(plexdm,opflow->initpflowpsection);CHKERRQ(ierr);
  ierr = DMSetGlobalSection(plexdm,opflow->initpflowpglobsection);CHKERRQ(ierr);

  ierr = DMCreateSectionSF(plexdm,opflow->initpflowpsection,opflow->initpflowpglobsection);CHKERRQ(ierr);

  /* Reset the edge and bus starting locations of variables */
  ierr = PSSetEdgeandBusStartLoc(ps);CHKERRQ(ierr);

  /* Solve */
  ierr = PFLOWSolve(initpflow);CHKERRQ(ierr);
  ierr = PFLOWConverged(initpflow,converged);CHKERRQ(ierr);

  /* Update bus and gen structs in pflow->ps */
  ierr = PFLOWPostSolve(initpflow);CHKERRQ(ierr);

  /*
     Update the ref. counts for init pflow sections so that they do not
     get destroyed when DMSetSection is called
  */
  ierr = PetscObjectReference((PetscObject)opflow->initpflowpsection);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)opflow->initpflowpglobsection);CHKERRQ(ierr);

  /* Reset the sections */
  ierr = DMSetSection(plexdm,opflow->defaultsection);CHKERRQ(ierr);
  ierr = DMSetGlobalSection(plexdm,opflow->defaultglobalsection);CHKERRQ(ierr);

  /* Reset SF */
  ierr = DMCreateSectionSF(plexdm,opflow->defaultsection,opflow->defaultglobalsection);CHKERRQ(ierr);

  /* Reset the bus and edge starting locations for variables */
  ierr = PSSetEdgeandBusStartLoc(ps);CHKERRQ(ierr);

  /* Reset dm */
  swap_dm(&opflow->initpflowdm,&opflow->ps->networkdm);
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
  ierr = PSSetApplication(opflow->ps,opflow,APP_ACOPF);CHKERRQ(ierr);

  opflow->Nconeq   = opflow->nconeq  = -1;
  opflow->Nconineq = opflow->nconineq = -1;
  opflow->Ncon     = opflow->ncon     = -1;
  opflow->Nx       = opflow->nx       = -1;
  opflow->Gi       = NULL;
  opflow->Lambdai  = NULL;
  
  opflow->obj_factor = 1.0;
  opflow->obj = 0.0;
  opflow->obj_gencost = PETSC_TRUE; /* Generation cost minimization ON by default */
  opflow->solutiontops = PETSC_FALSE;
  opflow->tolerance = 1e-6;

  opflow->has_gensetpoint = PETSC_FALSE;
  opflow->use_agc = PETSC_FALSE;

  opflow->solver   = NULL;
  opflow->model = NULL;

  opflow->initializationtype = OPFLOWINIT_MIDPOINT;

  opflow->objectivetype = MIN_GEN_COST;

  opflow->genbusvoltagetype = FIXED_WITHIN_QBOUNDS;

  opflow->nmodelsregistered = opflow->nsolversregistered = 0;
  opflow->OPFLOWModelRegisterAllCalled = opflow->OPFLOWSolverRegisterAllCalled = PETSC_FALSE;

  /* Register all models */
  ierr = OPFLOWModelRegisterAll(opflow);

  /* Register all solvers */
  ierr = OPFLOWSolverRegisterAll(opflow);

  /* Run-time options */
  opflow->ignore_lineflow_constraints = PETSC_FALSE;
  opflow->include_loadloss_variables = PETSC_FALSE;
  opflow->include_powerimbalance_variables = PETSC_FALSE;
  opflow->loadloss_penalty = 1e1;
  opflow->powerimbalance_penalty = 1e2;
  opflow->spdnordering = PETSC_FALSE;
  opflow->setupcalled = PETSC_FALSE;

  *opflowout = opflow;

  //  ierr = PetscPrintf(opflow->comm->type,"OPFLOW: Application created\n");
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

  if((*opflow)->initializationtype == OPFLOWINIT_ACPF) {
    /* Destroy objects created for initial power flow */
    ierr = PFLOWDestroy(&(*opflow)->initpflow);CHKERRQ(ierr);
    
    ierr = DMDestroy(&(*opflow)->initpflowdm);CHKERRQ(ierr);

    PetscObjectDereference((PetscObject)((*opflow)->initpflowpsection));
    PetscObjectDereference((PetscObject)((*opflow)->initpflowpglobsection));
    ierr = PetscSectionDestroy(&(*opflow)->initpflowpsection);CHKERRQ(ierr);
    ierr = PetscSectionDestroy(&(*opflow)->initpflowpglobsection);CHKERRQ(ierr);
  }


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
  ierr = VecDestroy(&(*opflow)->Gelocal);CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Lambdae);CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Lambdaelocal);CHKERRQ(ierr);
  if((*opflow)->nconineq) {
    ierr = VecDestroy(&(*opflow)->Gi);CHKERRQ(ierr);
    ierr = VecDestroy(&(*opflow)->Lambdai);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&(*opflow)->Gl);CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Gu);CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Lambda);CHKERRQ(ierr);

  /* Index sets and vecscatter for equality constraints */
  ierr = ISDestroy(&(*opflow)->isconeqlocal);CHKERRQ(ierr);
  ierr = ISDestroy(&(*opflow)->isconeqglob);CHKERRQ(ierr);
  ierr = VecScatterDestroy(&(*opflow)->scattereqcon);CHKERRQ(ierr);

  /* Jacobian of constraints */
  ierr = MatDestroy(&(*opflow)->Jac);CHKERRQ(ierr);
  ierr = MatDestroy(&(*opflow)->Jac_Ge);CHKERRQ(ierr);
  ierr = MatDestroy(&(*opflow)->Jac_Gi);CHKERRQ(ierr);

  ierr = MatDestroy(&(*opflow)->Hes);CHKERRQ(ierr);
  ierr = PSDestroy(&(*opflow)->ps);CHKERRQ(ierr);

  if((*opflow)->modelops.destroy) {
    ierr = ((*opflow)->modelops.destroy)(*opflow);
  }

  if((*opflow)->solverops.destroy) {
    ierr = ((*opflow)->solverops.destroy)(*opflow);
  }

  ierr = PetscFree((*opflow)->idxn2sd_map);CHKERRQ(ierr);
  ierr = PetscFree((*opflow)->busnvararray);CHKERRQ(ierr); 
  ierr = PetscFree((*opflow)->branchnvararray);CHKERRQ(ierr);
  ierr = PetscFree((*opflow)->eqconglobloc);CHKERRQ(ierr);

  ierr = PetscFree(*opflow);CHKERRQ(ierr);
  *opflow = 0;
  PetscFunctionReturn(0);
}

/* 
   OPFLOWSetTolerance - Set the solver tolerance

  Input Parameters:
+ opflow - opflow application object
- tol    - solver tolerance
*/
PetscErrorCode OPFLOWSetTolerance(OPFLOW opflow,PetscReal tol)
{
  PetscFunctionBegin;
  opflow->tolerance = tol;
  PetscFunctionReturn(0);
}

/**
 * OPFLOWGetNumIterations - Get the solver iterations
 * Input Parameters:
 * opflow - opflow application object
 * tol    - pointer to solver iterations variable that will be set
 */
PetscErrorCode OPFLOWGetNumIterations(OPFLOW opflow,PetscInt *its)
{
  PetscFunctionBegin;
  *its = opflow->numits;
  PetscFunctionReturn(0);
}

/**
 * OPFLOWGetTolerance - Get the solver tolerance
 * Input Parameters:
 * opflow - opflow application object
 * tol    - pointer to solver tolerance variable that will be set
 */
PetscErrorCode OPFLOWGetTolerance(OPFLOW opflow,PetscReal *tol)
{
  PetscFunctionBegin;
  *tol = opflow->tolerance;
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

  ierr = PetscStrcpy(opflow->solvername,solvername);CHKERRQ(ierr);
  /* Call the underlying implementation constructor */
  ierr = (*r)(opflow);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  OPFLOWSetModel - Sets the model for OPFLOW

  Input Parameters:
+ opflow - opflow application object
- modelname - name of the model
*/
PetscErrorCode OPFLOWSetModel(OPFLOW opflow,const char* modelname)
{
  PetscErrorCode ierr,(*r)(OPFLOW)=NULL;
  PetscInt       i;
  PetscFunctionBegin;
  PetscBool match;
  for(i=0;i < opflow->nmodelsregistered;i++) {
    ierr = PetscStrcmp(opflow->OPFLOWModelList[i].name,modelname,&match);CHKERRQ(ierr);
    if(match) {
      r = opflow->OPFLOWModelList[i].create;
      break;
    }
  }

  if(!r) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown type for OPFLOW Model %s",modelname);

  /* Null the function pointers */
  opflow->modelops.destroy                        = 0;
  opflow->modelops.setup                          = 0;
  opflow->modelops.setnumvariables                = 0;
  opflow->modelops.setnumconstraints              = 0;
  opflow->modelops.setvariablebounds              = 0;
  opflow->modelops.setvariableboundsarray         = 0;
  opflow->modelops.setconstraintbounds            = 0;
  opflow->modelops.setconstraintboundsarray       = 0;
  opflow->modelops.setvariableandconstraintbounds = 0;
  opflow->modelops.setvariableandconstraintboundsarray = 0;
  opflow->modelops.setinitialguess                = 0;
  opflow->modelops.setinitialguessarray           = 0;
  opflow->modelops.computeequalityconstraints     = 0;
  opflow->modelops.computeequalityconstraintsarray = 0;
  opflow->modelops.computeinequalityconstraints   = 0;
  opflow->modelops.computeinequalityconstraintsarray = 0;
  opflow->modelops.computeconstraints             = 0;
  opflow->modelops.computeconstraintsarray        = 0;
  opflow->modelops.computeequalityconstraintjacobian = 0;
  opflow->modelops.computeinequalityconstraintjacobian = 0;
  opflow->modelops.computehessian                 = 0;
  opflow->modelops.computeobjandgradient          = 0;
  opflow->modelops.computeobjective               = 0;
  opflow->modelops.computeobjectivearray          = 0;
  opflow->modelops.computegradient                = 0;
  opflow->modelops.computegradientarray           = 0;
  opflow->modelops.computejacobian                = 0;
  opflow->modelops.solutiontops                   = 0;

  ierr = PetscStrcpy(opflow->modelname,modelname);CHKERRQ(ierr);
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
  //  ierr = PetscPrintf(opflow->comm->type,"OPFLOW: Finished reading network data file %s\n",netfile);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* 
   OPFLOWSetNumConstraints - Sets the number of constraints for the OPFLOW problem

   Input Parameters:
.  OPFLOW - the opflow application object

   Output Parameters:
+  busnconeq - number of equality constraints at each bus
.  branchnconeq - number of equality constraints at each branch
.  nconeq   -  local number of equality constraints
-  nconineq -  local number of inequality constraints

*/
PetscErrorCode OPFLOWSetNumConstraints(OPFLOW opflow,PetscInt *branchnconeq, PetscInt *busnconeq,PetscInt *nconeq,PetscInt *nconineq)
{
  PetscErrorCode ierr;
  PetscInt       i,vStart,vEnd,eStart,eEnd;
  DM             networkdm=opflow->ps->networkdm;
  PS             ps=opflow->ps;
  DM             plexdm;
  PetscSection   buseqconsection,buseqconglobsection;
  PetscSection   varsection,varglobsection;
  PetscInt       nconeqloc=0;

  PetscFunctionBegin;
  ierr = (*opflow->modelops.setnumconstraints)(opflow,branchnconeq,busnconeq,nconeq,nconineq);

  ierr = PetscSectionCreate(opflow->comm->type,&buseqconsection);CHKERRQ(ierr);

  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
  ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);

  ierr = PetscSectionSetChart(buseqconsection,eStart,vEnd);CHKERRQ(ierr);

  for(i=0; i < ps->nline; i++) {
    ierr = PetscSectionSetDof(buseqconsection,eStart+i,0);
  }

  for(i=0; i < ps->nbus; i++) {
    ierr = PetscSectionSetDof(buseqconsection,vStart+i,busnconeq[i]);CHKERRQ(ierr);
    nconeqloc += busnconeq[i];
  }
  ierr = PetscSectionSetUp(buseqconsection);CHKERRQ(ierr);

  /* What we want to do here is have the DM set up the global section
     for the equality constraints (busnconeqglobsection) to get the global
     indices for equality constraints. We will use it later when operating 
     with equality constraints. 
     To do so, we need to 
     i) Get the local and global variables section (varsection and globvarsection) from the DM.
        Increase the reference count for the sections so that they do not get destroyed when
        the new section is set in ii
     ii) Set the section for equality constraints (busnconeqsection) on the DM
     iii) Have the DM generate the Global section and save it to busnconeqglobsection
     iv) Copy over the global indices from busnconeqglobsection.
      v) Reset the variable sections on the DM and decrement the reference count otherwise
         there will be a memory leak.
  */
  /* (i) */
  ierr = DMNetworkGetPlex(networkdm,&plexdm);CHKERRQ(ierr);
  ierr = DMGetSection(plexdm,&varsection);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)varsection);CHKERRQ(ierr);
  ierr = DMGetGlobalSection(plexdm,&varglobsection);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)varglobsection);CHKERRQ(ierr);

  /*  (ii) */
  ierr = DMSetSection(plexdm,buseqconsection);CHKERRQ(ierr);

  /* (iii) */
  ierr = DMGetGlobalSection(plexdm,&buseqconglobsection);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)buseqconglobsection);CHKERRQ(ierr);

  /* (iv) Set the global indices */
  ierr = PetscCalloc1(ps->nbus,&opflow->eqconglobloc);CHKERRQ(ierr);
  for(i=0; i < ps->nbus; i++) {
    ierr = DMNetworkGetVariableGlobalOffset(networkdm,vStart+i,&opflow->eqconglobloc[i]);CHKERRQ(ierr);
  }

  /* (v) */
  ierr = DMSetSection(plexdm,varsection);CHKERRQ(ierr);
  ierr = DMSetGlobalSection(plexdm,varglobsection);CHKERRQ(ierr);
  ierr = PetscObjectDereference((PetscObject)varsection);CHKERRQ(ierr);
  ierr = PetscObjectDereference((PetscObject)varglobsection);CHKERRQ(ierr);

  ierr = PetscSectionDestroy(&buseqconsection);CHKERRQ(ierr);
  ierr = PetscSectionDestroy(&buseqconglobsection);CHKERRQ(ierr);

  /* Create VecScatter for scattering values from global eq. constraints vector to local eq. constraints vector. */
  PetscInt *eqconloc,*eqconglobloc,ctr=0,j;
  ierr = PetscCalloc1(nconeqloc,&eqconloc);CHKERRQ(ierr);
  ierr = PetscCalloc1(nconeqloc,&eqconglobloc);CHKERRQ(ierr);
  for(i=0; i < ps->nbus; i++) {
    for(j = 0; j < busnconeq[i]; j++) {
      eqconloc[ctr] = ctr;
      eqconglobloc[ctr] = opflow->eqconglobloc[i] + j;
      ctr++;
    }
  }

  ierr = ISCreateGeneral(opflow->comm->type,nconeqloc,eqconloc,PETSC_OWN_POINTER,&opflow->isconeqlocal);CHKERRQ(ierr);
  ierr = ISCreateGeneral(opflow->comm->type,nconeqloc,eqconglobloc,PETSC_OWN_POINTER,&opflow->isconeqglob);CHKERRQ(ierr);
  
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

  ierr = (*opflow->modelops.setnumvariables)(opflow,busnvararray,branchnvararray,nx);

  for(i=0; i < ps->nline; i++) {
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

  /* Update starting locations for variables at each line */
  for(i=0; i < ps->nline; i++) {
    ierr = DMNetworkGetVariableOffset(networkdm,eStart+i,&ps->line[i].startloc);CHKERRQ(ierr);
    ierr = DMNetworkGetVariableGlobalOffset(networkdm,eStart+i,&ps->line[i].startlocglob);CHKERRQ(ierr);
  }

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
  char           modelname[32],solvername[32];
  PetscBool      modelset=PETSC_FALSE,solverset=PETSC_FALSE;
  PS             ps=opflow->ps;
  PetscInt       *branchnconeq,*busnconeq;
  PetscInt       sendbuf[3],recvbuf[3],i;

  PetscFunctionBegin;

  /* Default model and solver */
  PetscStrcpy(modelname,"POWER_BALANCE_POLAR");
  PetscStrcpy(solvername,"IPOPT");

  /* Read run-time options */
  ierr =  PetscOptionsBegin(opflow->comm->type,NULL,"OPFLOW options",NULL);CHKERRQ(ierr);

  ierr = PetscOptionsString("-opflow_model","OPFLOW model type","",modelname,modelname,32,&modelset);CHKERRQ(ierr);

  ierr = PetscOptionsString("-opflow_solver","OPFLOW solver type","",solvername,solvername,32,&solverset);CHKERRQ(ierr);
  ierr = PetscOptionsEnum("-opflow_initialization","Type of OPFLOW initialization","",OPFLOWInitializationTypes,(PetscEnum)opflow->initializationtype,(PetscEnum*)&opflow->initializationtype,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnum("-opflow_objective","Type of OPFLOW objective","",OPFLOWObjectiveTypes,(PetscEnum)opflow->objectivetype,(PetscEnum*)&opflow->objectivetype,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnum("-opflow_genbusvoltage","Type of OPFLOW gen bus voltage control","",OPFLOWGenBusVoltageTypes,(PetscEnum)opflow->genbusvoltagetype,(PetscEnum*)&opflow->genbusvoltagetype,NULL);CHKERRQ(ierr);

  ierr = PetscOptionsBool("-opflow_has_gensetpoint","Use set-points for generator real power","",opflow->has_gensetpoint,&opflow->has_gensetpoint,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-opflow_use_agc","Use automatic generation control (AGC)","",opflow->use_agc,&opflow->use_agc,NULL);CHKERRQ(ierr);
  if(opflow->use_agc) opflow->has_gensetpoint = PETSC_TRUE;
  ierr = PetscOptionsReal("-opflow_tolerance","optimization tolerance","",opflow->tolerance,&opflow->tolerance,NULL);CHKERRQ(ierr);


  if(opflow->objectivetype == MIN_GEN_COST) {
    opflow->obj_gencost = PETSC_TRUE;
  } else if(opflow->objectivetype == MIN_GENSETPOINT_DEVIATION) {
    opflow->obj_gencost = PETSC_FALSE;
    opflow->has_gensetpoint = PETSC_TRUE;
  }

  ierr = PetscOptionsBool("-opflow_ignore_lineflow_constraints","Ignore line flow constraints?","",opflow->ignore_lineflow_constraints,&opflow->ignore_lineflow_constraints,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-opflow_include_loadloss_variables","Include load loss?","",opflow->include_loadloss_variables,&opflow->include_loadloss_variables,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-opflow_loadloss_penalty","Penalty for load loss","",opflow->loadloss_penalty,&opflow->loadloss_penalty,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-opflow_include_powerimbalance_variables","Allow power imbalance?","",opflow->include_powerimbalance_variables,&opflow->include_powerimbalance_variables,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-opflow_powerimbalance_penalty","Power imbalance penalty","",opflow->powerimbalance_penalty,&opflow->powerimbalance_penalty,NULL);CHKERRQ(ierr);
  PetscOptionsEnd();

  /* Set model */
  if(modelset) {
    if(opflow->model) ierr = (*opflow->modelops.destroy)(opflow);
    ierr = OPFLOWSetModel(opflow,modelname);CHKERRQ(ierr);
    //    ierr = PetscPrintf(opflow->comm->type,"OPFLOW: Using %s model\n",modelname);CHKERRQ(ierr);
  } else {
    if(!opflow->model) {
      ierr = OPFLOWSetModel(opflow,OPFLOWMODEL_PBPOL);CHKERRQ(ierr);
      //      ierr = PetscPrintf(opflow->comm->type,"OPFLOW: Using %s model\n",OPFLOWMODEL_PBCAR);CHKERRQ(ierr);
    }
  }

  /* Set solver */
  if(solverset) {
    if(opflow->solver) ierr = (*opflow->solverops.destroy)(opflow);
    ierr = OPFLOWSetSolver(opflow,solvername);CHKERRQ(ierr);
    //    ierr = PetscPrintf(opflow->comm->type,"OPFLOW: Using %s solver\n",solvername);CHKERRQ(ierr);
  } else {
    if(!opflow->solver) {
#if defined(EXAGO_ENABLE_IPOPT)
      ierr = OPFLOWSetSolver(opflow,OPFLOWSOLVER_IPOPT);CHKERRQ(ierr);
#else
      ierr = OPFLOWSetSolver(opflow,OPFLOWSOLVER_TAO);CHKERRQ(ierr);
#endif
      //      ierr = PetscPrintf(opflow->comm->type,"OPFLOW: Using %s solver\n",OPFLOWSOLVER_IPOPT);CHKERRQ(ierr); 
    }
  }


  /* Set up underlying PS object */
  ierr = PSSetUp(ps);CHKERRQ(ierr);

  ierr = PetscCalloc1(ps->nbus,&opflow->busnvararray);CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->nline,&opflow->branchnvararray);CHKERRQ(ierr);

  /* Set up number of variables for branches and buses */
  ierr = OPFLOWSetNumVariables(opflow,opflow->busnvararray,opflow->branchnvararray,&opflow->nx);CHKERRQ(ierr);

  ierr = PetscCalloc1(ps->nline,&branchnconeq);CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->nbus,&busnconeq);CHKERRQ(ierr);
  /* Set up number of equality and inequality constraints and 
     number of equality constraints at each bus */

  /* Set number of constraints */
  ierr = OPFLOWSetNumConstraints(opflow,branchnconeq,busnconeq,&opflow->nconeq,&opflow->nconineq);CHKERRQ(ierr);
  opflow->ncon = opflow->nconeq + opflow->nconineq;

  ierr = PetscFree(branchnconeq);CHKERRQ(ierr);
  ierr = PetscFree(busnconeq);CHKERRQ(ierr);

  sendbuf[0] = opflow->nx;
  sendbuf[1] = opflow->nconeq;
  sendbuf[2] = opflow->nconineq;
  ierr = MPI_Allreduce(sendbuf,recvbuf,3,MPIU_INT,MPI_SUM,opflow->comm->type);CHKERRQ(ierr);
  opflow->Nx = recvbuf[0];
  opflow->Nconeq = recvbuf[1];
  opflow->Nconineq = recvbuf[2];
  opflow->Ncon = opflow->Nconeq + opflow->Nconineq;
  //  ierr = PetscPrintf(PETSC_COMM_SELF,"OPFLOW: Rank %d: nx = %d nconeq = %d, nconineq = %d, ncon = %d\n",opflow->comm->rank,opflow->nx,opflow->nconeq,opflow->nconineq,opflow->ncon);CHKERRQ(ierr);

  /* Set vertex local to global ordering */
  ierr = DMNetworkSetVertexLocalToGlobalOrdering(opflow->ps->networkdm);CHKERRQ(ierr);

  /* Create solution vector and upper/lower bounds */
  ierr = PSCreateGlobalVector(opflow->ps,&opflow->X);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->gradobj);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->Xl);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->Xu);CHKERRQ(ierr);
  /* Default variable bounds are -inf <= x <= +inf */
  ierr = VecSet(opflow->Xl,PETSC_NINFINITY);CHKERRQ(ierr);
  ierr = VecSet(opflow->Xu,PETSC_INFINITY);CHKERRQ(ierr);

  /* Vectors for constraints */
  ierr = VecCreate(ps->comm->type,&opflow->G);CHKERRQ(ierr);
  ierr = VecSetSizes(opflow->G,opflow->ncon,opflow->Ncon);CHKERRQ(ierr);
  ierr = VecSetFromOptions(opflow->G);CHKERRQ(ierr);
  
  ierr = VecDuplicate(opflow->G,&opflow->Gl);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->G,&opflow->Gu);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->G,&opflow->Lambda);CHKERRQ(ierr);

  ierr = VecCreate(ps->comm->type,&opflow->Ge);CHKERRQ(ierr);
  ierr = VecSetSizes(opflow->Ge,opflow->nconeq,opflow->Nconeq);CHKERRQ(ierr);
  ierr = VecSetFromOptions(opflow->Ge);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->Ge,&opflow->Lambdae);CHKERRQ(ierr);

  PetscInt nconeqlocal;

  ierr = ISGetSize(opflow->isconeqlocal,&nconeqlocal);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_SELF,&opflow->Gelocal);CHKERRQ(ierr);
  ierr = VecSetSizes(opflow->Gelocal,nconeqlocal,nconeqlocal);CHKERRQ(ierr);
  ierr = VecSetFromOptions(opflow->Gelocal);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->Gelocal,&opflow->Lambdaelocal);CHKERRQ(ierr);

  ierr = VecScatterCreate(opflow->Ge,opflow->isconeqglob,opflow->Gelocal,opflow->isconeqlocal,&opflow->scattereqcon);CHKERRQ(ierr);

  if(opflow->Nconineq) {
    ierr = VecCreate(ps->comm->type,&opflow->Gi);CHKERRQ(ierr);
    ierr = VecSetSizes(opflow->Gi,opflow->nconineq,opflow->Nconineq);CHKERRQ(ierr);
    ierr = VecSetFromOptions(opflow->Gi);CHKERRQ(ierr);
    ierr = VecDuplicate(opflow->Gi,&opflow->Lambdai);CHKERRQ(ierr);
  }

  /* Create equality and inequality constraint Jacobian matrices */
  ierr = MatCreate(opflow->comm->type,&opflow->Jac_Ge);CHKERRQ(ierr);
  ierr = MatSetSizes(opflow->Jac_Ge,opflow->nconeq,opflow->nx,opflow->Nconeq,opflow->Nx);CHKERRQ(ierr);
  ierr = MatSetUp(opflow->Jac_Ge);CHKERRQ(ierr);
  ierr = MatSetFromOptions(opflow->Jac_Ge);CHKERRQ(ierr);

  if(opflow->Nconineq) {
    ierr = MatCreate(opflow->comm->type,&opflow->Jac_Gi);CHKERRQ(ierr);
    ierr = MatSetSizes(opflow->Jac_Gi,opflow->nconineq,opflow->nx,opflow->Nconineq,opflow->Nx);CHKERRQ(ierr);
    ierr = MatSetUp(opflow->Jac_Gi);CHKERRQ(ierr);
    ierr = MatSetFromOptions(opflow->Jac_Gi);CHKERRQ(ierr);
  }

  /* Create Hessian */
  ierr = PSCreateMatrix(opflow->ps,&opflow->Hes);CHKERRQ(ierr);


  /* Create natural to sparse dense ordering mapping (needed for some models)
     Here, we create the mapping array and fill up natural ordering. The model
     set up function should change the mapping to what it needs
  */
  ierr = PetscMalloc1(opflow->nx,&opflow->idxn2sd_map);CHKERRQ(ierr);
  for(i=0; i < opflow->nx; i++) opflow->idxn2sd_map[i] = i;

  /* Model set up */
  if(opflow->modelops.setup) {
    ierr = (*opflow->modelops.setup)(opflow);CHKERRQ(ierr);
  }

  /* Set bounds on variables */
  ierr = OPFLOWComputeVariableBounds(opflow,opflow->Xl,opflow->Xu);CHKERRQ(ierr);

  /* Set bounds on constraints */
  ierr = OPFLOWComputeConstraintBounds(opflow,opflow->Gl,opflow->Gu);CHKERRQ(ierr);

  /* Set initial guess */
  if(opflow->initializationtype == OPFLOWINIT_ACPF) {
    ierr = OPFLOWSetUpInitPflow(opflow);CHKERRQ(ierr);
  }
  ierr = OPFLOWSetInitialGuess(opflow,opflow->X);CHKERRQ(ierr);

  /* Initial guess for multipliers */
  ierr = VecSet(opflow->Lambda,1.0);CHKERRQ(ierr);
  ierr = VecSet(opflow->Lambdae,1.0);CHKERRQ(ierr);
  if(opflow->Nconineq) {
    ierr = VecSet(opflow->Lambdai,1.0);CHKERRQ(ierr);
  }

  /* Solver set up */
  if(opflow->solverops.setup) {
    ierr = (*opflow->solverops.setup)(opflow);CHKERRQ(ierr);
  }
  //  ierr = PetscPrintf(opflow->comm->type,"OPFLOW: Setup completed\n");CHKERRQ(ierr);


  /* Register events for logging */
  ierr = PetscLogEventRegister("OPFLOWObj",0,&opflow->objlogger);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OPFLOWGrad",0,&opflow->gradlogger);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OPFLOWEqCons",0,&opflow->eqconslogger);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OPFLOWIneqCons",0,&opflow->ineqconslogger);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OPFLOWEqConsJac",0,&opflow->eqconsjaclogger);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OPFLOWIneqConsJac",0,&opflow->ineqconsjaclogger);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OPFLOWHess",0,&opflow->hesslogger);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OPFLOWSolve",0,&opflow->solvelogger);CHKERRQ(ierr);

  /* Compute area participation factors */
  ierr = PSComputeParticipationFactors(ps);CHKERRQ(ierr);
  opflow->setupcalled = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/* OPFLOWSetInitialGuess - Sets the initial guess for OPFLOW

   Input Parameters:
.  opflow - the optimal power flow application object

   Output Parameters:
.  X - initial guess

*/
PetscErrorCode OPFLOWSetInitialGuess(OPFLOW opflow, Vec X)
{
  PetscErrorCode ierr;
  PetscBool      pflowconverged=PETSC_FALSE;

  PetscFunctionBegin;

  /* Set initial guess */
  switch (opflow->initializationtype) {
    case OPFLOWINIT_MIDPOINT:
    case OPFLOWINIT_FROMFILE:
    case OPFLOWINIT_FLATSTART:
      if(opflow->modelops.setinitialguess) {
	ierr = (*opflow->modelops.setinitialguess)(opflow,X);CHKERRQ(ierr);
      }
      break;
    case OPFLOWINIT_ACPF:
      ierr = OPFLOWComputePrePflow(opflow,&pflowconverged);CHKERRQ(ierr);
      if(!pflowconverged) {
	SETERRQ(PETSC_COMM_SELF,0,"AC power flow initialization did not converged\n");
      }
      ierr = (*opflow->modelops.setinitialguess)(opflow,X);CHKERRQ(ierr);

      break;
    default:
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Unknown OPFLOW initialization type\n");
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

  /* Solve */
  ierr = PetscLogEventBegin(opflow->solvelogger,0,0,0,0);CHKERRQ(ierr);
  ierr = (*opflow->solverops.solve)(opflow);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->solvelogger,0,0,0,0);CHKERRQ(ierr);

  //  ierr = VecView(opflow->X,0);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  OPFLOWGetObjective - Returns the objective function value

  Input Parameters:
- OPFLOW - the OPFLOW object

  Output Parameters:
- obj    - the objective function value

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode OPFLOWGetObjective(OPFLOW opflow,PetscReal *obj)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*opflow->solverops.getobjective)(opflow,obj);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetVariableBounds - Returns the variable bounds

  Input Parameters:
- OPFLOW - the OPFLOW object

  Output Parameters:
+ Xl     - lower bound on X
- Xu     - upper bound on X

  Notes: Should be called after the optimization finishes. This function merely
         returns the already computed lower and upper bound vectors.
*/
PetscErrorCode OPFLOWGetVariableBounds(OPFLOW opflow,Vec *Xl, Vec *Xu)
{
  PetscFunctionBegin;
  *Xl = opflow->Xl;
  *Xu = opflow->Xu;
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeObjective - Computes the objective function value

  Input Parameters:
+ OPFLOW - the OPFLOW object
- X      - the solution vector

  Output Parameters:
- obj    - the objective function value

*/
PetscErrorCode OPFLOWComputeObjective(OPFLOW opflow,Vec X,PetscReal *obj)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->objlogger,0,0,0,0);CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeobjective)(opflow,X,obj);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->objlogger,0,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeObjectiveArray - Computes the objective function value

  Input Parameters:
+ OPFLOW - the OPFLOW object
- x      - the solution array

  Output Parameters:
- obj    - the objective function value

*/
PetscErrorCode OPFLOWComputeObjectiveArray(OPFLOW opflow,const double* x,double *obj)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->objlogger,0,0,0,0);CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeobjectivearray)(opflow,x,obj);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->objlogger,0,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*
  OPFLOWComputeVariableBounds - Computes the bounds on X

  Input Parameters:
- OPFLOW - the OPFLOW object

  Output Parameters:
+ Xl     - lower bound on X
- Xu     - upper bound on X

  Notes: Should be called after the optimization is set up.
*/
PetscErrorCode OPFLOWComputeVariableBounds(OPFLOW opflow,Vec Xl, Vec Xu)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if(opflow->modelops.setvariablebounds) {
    ierr = (*opflow->modelops.setvariablebounds)(opflow,Xl,Xu);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeGradient - Computes the gradient of the objective function

  Input Parameters:
+ OPFLOW - the OPFLOW object
. X      - the solution vector
- grad    - the gradient of the objective function

  Notes: Should be called after the optimization is set up
*/
PetscErrorCode OPFLOWComputeGradient(OPFLOW opflow,Vec X,Vec grad)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->gradlogger,0,0,0,0);CHKERRQ(ierr);
  ierr = (*opflow->modelops.computegradient)(opflow,X,grad);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->gradlogger,0,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeGradientArray - Computes the gradient of the objective function (array version)

  Input Parameters:
+ OPFLOW - the OPFLOW object
. x      - the solution array
- grad    - the gradient of the objective function

  Notes: Should be called after the optimization is set up
*/
PetscErrorCode OPFLOWComputeGradientArray(OPFLOW opflow,const double *x,double *grad)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->gradlogger,0,0,0,0);CHKERRQ(ierr);
  ierr = (*opflow->modelops.computegradientarray)(opflow,x,grad);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->gradlogger,0,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*
  OPFLOWGetSolution - Returns the OPFLOW solution

  Input Parameters:
+ OPFLOW - the OPFLOW object
- X      - the opflow solution

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode OPFLOWGetSolution(OPFLOW opflow,Vec *X)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*opflow->solverops.getsolution)(opflow,X);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetConstraints - Returns the OPFLOW constraints

  Input Parameters:
+ OPFLOW - the OPFLOW object

  Output Parameters:
- G    - the opflow constraints

  Notes: Should be called after the optimization finishes.
         Equality constraints first followed by inequality constraints
*/
PetscErrorCode OPFLOWGetConstraints(OPFLOW opflow,Vec *G)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*opflow->solverops.getconstraints)(opflow,G);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeEqualityConstraints - Computes OPFLOW equality constraints
  Input Parameters:
+ OPFLOW - the OPFLOW object
. X      - the solution vector

  Output Parameters:
- Ge      - OPFLOW equallity constraints

  Notes: Should be called after the optimization is set up.

*/
PetscErrorCode OPFLOWComputeEqualityConstraints(OPFLOW opflow,Vec X,Vec Ge)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->eqconslogger,0,0,0,0);CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeequalityconstraints)(opflow,X,Ge);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->eqconslogger,0,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeInequalityConstraints - Computes OPFLOW inequality constraints
  Input Parameters:
+ OPFLOW - the OPFLOW object
. X      - the solution vector

  Output Parameters:
- Gi      - OPFLOW inequallity constraints

  Notes: Should be called after the optimization is set up.
*/
PetscErrorCode OPFLOWComputeInequalityConstraints(OPFLOW opflow,Vec X,Vec Gi)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->ineqconslogger,0,0,0,0);CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeinequalityconstraints)(opflow,X,Gi);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->ineqconslogger,0,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeEqualityConstraintsArray - Computes OPFLOW equality constraints (array version)

  Input Parameters:
+ OPFLOW - the OPFLOW object
. x      - the solution array

  Output Parameters:
- ge      - OPFLOW equality constraints array

  Notes: Should be called after the optimization is set up.
*/
PetscErrorCode OPFLOWComputeEqualityConstraintsArray(OPFLOW opflow,const double* x,double *ge)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->eqconslogger,0,0,0,0);CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeequalityconstraintsarray)(opflow,x,ge);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->eqconslogger,0,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeInequalityConstraintsArray - Computes OPFLOW inequality constraints (array version)

  Input Parameters:
+ OPFLOW - the OPFLOW object
. x      - the solution array

  Output Parameters:
- gi      - OPFLOW inequality constraints array

  Notes: Should be called after the optimization is set up.
*/
PetscErrorCode OPFLOWComputeInequalityConstraintsArray(OPFLOW opflow,const double* x,double *gi)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->ineqconslogger,0,0,0,0);CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeinequalityconstraintsarray)(opflow,x,gi);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->ineqconslogger,0,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeConstraints - Computes the OPFLOW constraints

  Input Parameters:
+ OPFLOW - the OPFLOW object
. X      - the solution vector

  Output Parameters:
- G    - the opflow constraints

  Notes: Should be called after the optimization is set up. G has
         equality constraints followed by inequality constraints
*/
PetscErrorCode OPFLOWComputeConstraints(OPFLOW opflow,Vec X,Vec G)
{
  PetscErrorCode ierr;
  PetscScalar    *g;

  PetscFunctionBegin;
  ierr = VecGetArray(G,&g);CHKERRQ(ierr);

  ierr = VecPlaceArray(opflow->Ge,g);CHKERRQ(ierr);
  ierr = OPFLOWComputeEqualityConstraints(opflow,X,opflow->Ge);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->Ge);CHKERRQ(ierr);

  if(opflow->Nconineq) {
    ierr = VecPlaceArray(opflow->Gi,g+opflow->nconeq);CHKERRQ(ierr);
    ierr = OPFLOWComputeInequalityConstraints(opflow,X,opflow->Gi);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Gi);CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(G,&g);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeHessian - Computes the OPFLOW Hessian

  Input Parameters:
+ OPFLOW - the OPFLOW object
. X      - the solution vector
. Lambda - Lagrange multiplier
- obj_factor - scaling factor for objective part of Hessian (used by IPOPT)

  Output Parameters:
- H    - the Hessian

  Notes: Should be called after the optimization is set up. Lambda has multipliers -
         equality constraints followed by inequality constraints
*/
PetscErrorCode OPFLOWComputeHessian(OPFLOW opflow,Vec X, Vec Lambda, PetscScalar obj_factor, Mat H)
{
  PetscErrorCode ierr;
  PetscScalar    *lambda;

  PetscFunctionBegin;

  opflow->obj_factor = obj_factor;

  /* IPOPT passes the scaled Lagrangian multipliers for Hessian calculation. The scaling
     factor is the Hessian objective factor it uses in the Hessian calculation
  */
  ierr = VecScale(Lambda,obj_factor);CHKERRQ(ierr);

  ierr = VecGetArray(Lambda,&lambda);CHKERRQ(ierr);

  ierr = VecPlaceArray(opflow->Lambdae,lambda);CHKERRQ(ierr);
  if(opflow->Nconineq) {
    ierr = VecPlaceArray(opflow->Lambdai,lambda+opflow->nconeq);CHKERRQ(ierr);
  }

  ierr = PetscLogEventBegin(opflow->hesslogger,0,0,0,0);CHKERRQ(ierr);
  /* Compute Hessian */
  ierr = (*opflow->modelops.computehessian)(opflow,X,opflow->Lambdae,opflow->Lambdai,H);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->hesslogger,0,0,0,0);CHKERRQ(ierr);

  ierr = VecResetArray(opflow->Lambdae);CHKERRQ(ierr);

  if(opflow->Nconineq) {
    ierr = VecResetArray(opflow->Lambdai);CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(Lambda,&lambda);CHKERRQ(ierr);

  /* Rescale Lambda back to its original value */
  ierr = VecScale(Lambda,1./obj_factor);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeConstraintJacobian - Computes the OPFLOW constraint Jacobians

  Input Parameters:
+ OPFLOW - the OPFLOW object
. X      - the solution vector

  Output Parameters:
+ Jeq    - equality constraint Jacobian
- Jineq  - Inequality constraint Jacobian

  Notes: Should be called after the optimization is set up. 
*/
PetscErrorCode OPFLOWComputeConstraintJacobian(OPFLOW opflow,Vec X,Mat Jeq, Mat Jineq)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MatZeroEntries(Jeq);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(opflow->eqconsjaclogger,0,0,0,0);CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeequalityconstraintjacobian)(opflow,X,Jeq);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->eqconsjaclogger,0,0,0,0);CHKERRQ(ierr);

  if(opflow->Nconineq) {
    ierr = MatZeroEntries(Jineq);CHKERRQ(ierr);
    ierr = PetscLogEventBegin(opflow->ineqconsjaclogger,0,0,0,0);CHKERRQ(ierr);
    ierr = (*opflow->modelops.computeinequalityconstraintjacobian)(opflow,X,Jineq);CHKERRQ(ierr);
    ierr = PetscLogEventEnd(opflow->ineqconsjaclogger,0,0,0,0);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}



/*
  OPFLOWGetConstraintBounds - Returns the OPFLOW constraint bounds

  Input Parameters:
+ OPFLOW - the OPFLOW object

  Output Parameters:
+ Gl    - lower bound on constraints
- Gu    - upper bound on constraints

  Notes: Should be called after the optimization finishes.
         Equality constraint bounds first followed by inequality constraints
*/
PetscErrorCode OPFLOWGetConstraintBounds(OPFLOW opflow,Vec *Gl,Vec *Gu)
{

  PetscFunctionBegin;
  *Gl = opflow->Gl;
  *Gu = opflow->Gu;
  PetscFunctionReturn(0);
}


/*
  OPFLOWGetConstraintJacobian - Returns the OPFLOW equality and inequality constraint jacobian

  Input Parameters:
+ OPFLOW - the OPFLOW object

  Output Parameters:
+ Jeq    - equality constraint Jacobian
- Jineq  - inequality constraint Jacobian

  Notes: Should be called after the optimization finishes.
         Jineq is set to NULL if no inequality constraints exists
*/
PetscErrorCode OPFLOWGetConstraintJacobian(OPFLOW opflow,Mat *Jeq,Mat *Jineq)
{

  PetscFunctionBegin;

  *Jeq = opflow->Jac_Ge;
  if(opflow->nconineq) *Jineq = opflow->Jac_Gi;
  else *Jineq = NULL;

  PetscFunctionReturn(0);
}

/*
  OPFLOWGetHessian - Returns the OPFLOW Hessian

  Input Parameters:
+ OPFLOW - the OPFLOW object

  Output Parameters:
+ Hess    - Hessian matrix
- obj_factor - factor used in scaling objective Hessian (used by IPOPT)
  Notes: Should be called after the optimization finishes.
*/
PetscErrorCode OPFLOWGetHessian(OPFLOW opflow,Mat *Hess,PetscScalar *obj_factor)
{

  PetscFunctionBegin;

  *Hess = opflow->Hes;
  *obj_factor = opflow->obj_factor;

  PetscFunctionReturn(0);
}


/*
  OPFLOWComputeConstraintBounds - Computes the bounds on constraints

  Input Parameters:
- OPFLOW - the OPFLOW object

  Output Parameters:
+ Gl     - lower bound on constraints
- Gu     - upper bound on constraints

  Notes: Should be called after the optimization is set up.
*/
PetscErrorCode OPFLOWComputeConstraintBounds(OPFLOW opflow,Vec Gl, Vec Gu)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if(opflow->modelops.setconstraintbounds) {
    ierr = (*opflow->modelops.setconstraintbounds)(opflow,Gl,Gu);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetConstraintMultipliers - Returns the OPFLOW constraint multipliers

  Input Parameters:
+ OPFLOW - the OPFLOW object

  Output Parameters:
- G    - the opflow constraint lagrange multipliers

  Notes: Should be called after the optimization finishes.
    Equality constraint multipliers first followed by inequality constraint multipliers
*/
PetscErrorCode OPFLOWGetConstraintMultipliers(OPFLOW opflow,Vec *Lambda)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*opflow->solverops.getconstraintmultipliers)(opflow,Lambda);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetConvergenceStatus - Did OPFLOW converge?

  Input Parameters:
+ OPFLOW - the OPFLOW object
- status - PETSC_TRUE if converged, PETSC_FALSE otherwise

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode OPFLOWGetConvergenceStatus(OPFLOW opflow,PetscBool *status)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*opflow->solverops.getconvergencestatus)(opflow,status);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWSolutionToPS - Updates the PS struct from OPFLOW solution

  Input Parameters:
. opflow - the OPFLOW object

  Notes: Updates the different fields in the PS struct from the OPFLOW solution
*/
PetscErrorCode OPFLOWSolutionToPS(OPFLOW opflow)
{
  PetscErrorCode     ierr;

  PetscFunctionBegin;

  ierr = (*opflow->modelops.solutiontops)(opflow);

  opflow->solutiontops = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/*
  OPFLOWSetObjectiveType - Sets the objective function for OPFLOW

  Input Parameters
+ opflow - the opflow object
- objtype - type of objective function

  Command-line option: -opflow_objective <obj_type>

  Notes: Must be called before OPFLOWSetUp
*/
PetscErrorCode OPFLOWSetObjectiveType(OPFLOW opflow,OPFLOWObjectiveType objtype)
{
  PetscFunctionBegin;
  opflow->objectivetype = objtype;
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetObjectiveType - Gets the objective function for OPFLOW

  Input Parameters
+ opflow - the opflow object
- objtype - type of objective function

  Notes: Must be called after OPFLOWSetUp
*/
PetscErrorCode OPFLOWGetObjectiveType(OPFLOW opflow,OPFLOWObjectiveType *objtype)
{
  PetscFunctionBegin;
  *objtype = opflow->objectivetype;
  PetscFunctionReturn(0);
}

