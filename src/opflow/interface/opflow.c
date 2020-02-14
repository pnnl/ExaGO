#include <private/pflowimpl.h>
#include <private/opflowimpl.h>
#include <petsc/private/dmnetworkimpl.h>

const char *const OPFLOWInitializationTypes[] = {"MIDPOINT","FROMFILE","ACPF","FLATSTART","OPFLOWInitializationType","OPFLOWINIT_",NULL};

void swap_dm(DM *dm1, DM *dm2)
{
  DM temp = *dm1;
  *dm1 = *dm2;
  *dm2 = temp;
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
  ierr = DMGetDefaultSection(plexdm,&opflow->defaultsection);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)opflow->defaultsection);CHKERRQ(ierr);

  ierr = DMGetDefaultGlobalSection(plexdm,&opflow->defaultglobalsection);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)opflow->defaultglobalsection);CHKERRQ(ierr);

  /* Set the new section created for initial power flow */
  ierr = DMSetDefaultSection(plexdm,opflow->initpflowpsection);CHKERRQ(ierr);
  ierr = DMGetDefaultGlobalSection(plexdm,&opflow->initpflowpglobsection);CHKERRQ(ierr);

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
  ierr = DMSetDefaultSection(plexdm,opflow->defaultsection);CHKERRQ(ierr);
  ierr = DMSetDefaultGlobalSection(plexdm,opflow->defaultglobalsection);CHKERRQ(ierr);

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
  ierr = DMSetDefaultSection(plexdm,opflow->initpflowpsection);CHKERRQ(ierr);
  ierr = DMSetDefaultGlobalSection(plexdm,opflow->initpflowpglobsection);CHKERRQ(ierr);

  ierr = DMCreateDefaultSF(plexdm,opflow->initpflowpsection,opflow->initpflowpglobsection);CHKERRQ(ierr);

  /* Reset the edge and bus starting locations of variables */
  ierr = PSSetEdgeandBusStartLoc(ps);CHKERRQ(ierr);

  /* Solve */
  ierr = PFLOWSolve(initpflow);CHKERRQ(ierr);
  ierr = PFLOWConverged(initpflow,converged);CHKERRQ(ierr);

  /* Update bus and gen structs in pflow->ps */
  ierr = PFLOWPostSolve(initpflow);CHKERRQ(ierr);

  /*
     Update the ref. counts for init pflow sections so that they do not
     get destroyed when DMSetDefaultSection is called
  */
  ierr = PetscObjectReference((PetscObject)opflow->initpflowpsection);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)opflow->initpflowpglobsection);CHKERRQ(ierr);

  /* Reset the sections */
  ierr = DMSetDefaultSection(plexdm,opflow->defaultsection);CHKERRQ(ierr);
  ierr = DMSetDefaultGlobalSection(plexdm,opflow->defaultglobalsection);CHKERRQ(ierr);

  /* Reset SF */
  ierr = DMCreateDefaultSF(plexdm,opflow->defaultsection,opflow->defaultglobalsection);CHKERRQ(ierr);

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
  opflow->compute_obj_hessian = PETSC_TRUE;

  opflow->solver   = NULL;
  opflow->formulation = NULL;

  opflow->initializationtype = OPFLOWINIT_MIDPOINT;

  opflow->nformulationsregistered = opflow->nsolversregistered = 0;
  opflow->OPFLOWFormulationRegisterAllCalled = opflow->OPFLOWSolverRegisterAllCalled = PETSC_FALSE;

  /* Register all formulations */
  ierr = OPFLOWFormulationRegisterAll(opflow);

  /* Register all solvers */
  ierr = OPFLOWSolverRegisterAll(opflow);

  /* Run-time options */
  opflow->ignore_lineflow_constraints = PETSC_FALSE;
  opflow->include_loadloss_variables = PETSC_FALSE;
  opflow->include_powerimbalance_variables = PETSC_FALSE;
  opflow->loadloss_penalty = 1e1;
  opflow->powerimbalance_penalty = 1e2;
  opflow->setupcalled = PETSC_FALSE;

  *opflowout = opflow;

  ierr = PetscPrintf(opflow->comm->type,"OPFLOW: Application created\n");
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

  if((*opflow)->solverops.destroy) {
    ierr = ((*opflow)->solverops.destroy)(*opflow);
  }

  if((*opflow)->formops.destroy) {
    ierr = ((*opflow)->formops.destroy)(*opflow);
  }

  ierr = PetscFree((*opflow)->eqconglobloc);CHKERRQ(ierr);

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
  ierr = PetscPrintf(opflow->comm->type,"OPFLOW: Finished reading network data file %s\n",netfile);CHKERRQ(ierr);
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
  ierr = (*opflow->formops.setnumconstraints)(opflow,branchnconeq,busnconeq,nconeq,nconineq);

  ierr = PetscSectionCreate(opflow->comm->type,&buseqconsection);CHKERRQ(ierr);

  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
  ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);

  ierr = PetscSectionSetChart(buseqconsection,eStart,vEnd);CHKERRQ(ierr);

  for(i=0; i < ps->nbranch; i++) {
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

  /* Update starting locations for variables at each line */
  for(i=0; i < ps->nbranch; i++) {
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
   OPFLOWSetFromOptions - Sets options for optimal power flow object

   Input Parameters:
.  opflow - the OPFLOW object
*/
PetscErrorCode OPFLOWSetFromOptions(OPFLOW opflow)
{
  char           formulationname[32],solvername[32];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr =  PetscOptionsBegin(opflow->comm->type,NULL,"OPFLOW options",NULL);CHKERRQ(ierr);
  PetscStrcpy(formulationname,"POWER_BALANCE_CARTESIAN");
  ierr = PetscOptionsString("-opflow_formulation","OPFLOW formulation type","",formulationname,formulationname,32,NULL);CHKERRQ(ierr);
  PetscStrcpy(solvername,"IPOPT");
  ierr = PetscOptionsString("-opflow_solver","OPFLOW solver type","",solvername,solvername,32,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnum("-opflow_initialization","Type of OPFLOW initialization","",OPFLOWInitializationTypes,(PetscEnum)opflow->initializationtype,(PetscEnum*)&opflow->initializationtype,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-opflow_ignore_lineflow_constraints","Ignore line flow constraints?","",opflow->ignore_lineflow_constraints,&opflow->ignore_lineflow_constraints,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-opflow_include_loadloss_variables","Include load loss?","",opflow->include_loadloss_variables,&opflow->include_loadloss_variables,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-opflow_loadloss_penalty","Penalty for load loss","",opflow->loadloss_penalty,&opflow->loadloss_penalty,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-opflow_include_powerimbalance_variables","Allow power imbalance?","",opflow->include_powerimbalance_variables,&opflow->include_powerimbalance_variables,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-opflow_powerimbalance_penalty","Power imbalance penalty","",opflow->powerimbalance_penalty,&opflow->powerimbalance_penalty,NULL);CHKERRQ(ierr);
  PetscOptionsEnd();

  /* Set formulation */
  ierr = OPFLOWSetFormulation(opflow,formulationname);CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"OPFLOW: Using %s formulation\n",formulationname);CHKERRQ(ierr);

  /* Set solver */
  ierr = OPFLOWSetSolver(opflow,solvername);CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"OPFLOW: Using %s solver\n",solvername);CHKERRQ(ierr);

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
  PetscInt       *busnvararray,*branchnvararray;
  PetscInt       *branchnconeq,*busnconeq;
  PetscInt       sendbuf[3],recvbuf[3];

  PetscFunctionBegin;

  ierr = OPFLOWSetFromOptions(opflow);CHKERRQ(ierr);

  /* Set up underlying PS object */
  ierr = PSSetUp(ps);CHKERRQ(ierr);

  ierr = PetscCalloc1(ps->nbus,&busnvararray);CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->nbranch,&branchnvararray);CHKERRQ(ierr);

  /* Set up number of variables for branches and buses */
  ierr = OPFLOWSetNumVariables(opflow,busnvararray,branchnvararray,&opflow->nx);CHKERRQ(ierr);

  ierr = PetscCalloc1(ps->nbranch,&branchnconeq);CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->nbus,&busnconeq);CHKERRQ(ierr);
  /* Set up number of equality and inequality constraints and 
     number of equality constraints at each bus */

  /* Set number of constraints */
  ierr = OPFLOWSetNumConstraints(opflow,branchnconeq,busnconeq,&opflow->nconeq,&opflow->nconineq);CHKERRQ(ierr);
  opflow->ncon = opflow->nconeq + opflow->nconineq;

  ierr = PetscFree(busnvararray);CHKERRQ(ierr); 
  ierr = PetscFree(branchnvararray);CHKERRQ(ierr);
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
  ierr = PetscPrintf(PETSC_COMM_SELF,"OPFLOW: Rank %d: nx = %d nconeq = %d, nconineq = %d, ncon = %d\n",opflow->comm->rank,opflow->nx,opflow->nconeq,opflow->nconineq,opflow->ncon);CHKERRQ(ierr);

  /* Set vertex local to global ordering */
  ierr = DMNetworkSetVertexLocalToGlobalOrdering(opflow->ps->networkdm);CHKERRQ(ierr);

  /* Create solution vector and upper/lower bounds */
  ierr = PSCreateGlobalVector(opflow->ps,&opflow->X);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->gradobj);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->Xl);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->Xu);CHKERRQ(ierr);
  /* Default variable bounds are -inf <= x <= +inf */
  ierr = VecSet(opflow->Xl,PETSC_NINFINITY);CHKERRQ(ierr);
  ierr = VecSet(opflow->Xu,PETSC_NINFINITY);CHKERRQ(ierr);

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

  ierr = (*opflow->solverops.setup)(opflow);CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"OPFLOW: Setup completed\n");CHKERRQ(ierr);

  if(opflow->initializationtype == OPFLOWINIT_ACPF) {
    ierr = OPFLOWSetUpInitPflow(opflow);CHKERRQ(ierr);
  }

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
      if(opflow->formops.setinitialguess) {
	ierr = (*opflow->formops.setinitialguess)(opflow,X);CHKERRQ(ierr);
      }
      break;
    case OPFLOWINIT_ACPF:
      ierr = OPFLOWComputePrePflow(opflow,&pflowconverged);CHKERRQ(ierr);
      if(!pflowconverged) {
	SETERRQ(PETSC_COMM_SELF,0,"AC power flow initialization did not converged\n");
      }
      ierr = (*opflow->formops.setinitialguess)(opflow,X);CHKERRQ(ierr);

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

  ierr = OPFLOWSetInitialGuess(opflow,opflow->X);CHKERRQ(ierr);

  ierr = VecSet(opflow->Lambda,1.0);CHKERRQ(ierr);
  ierr = VecSet(opflow->Lambdae,1.0);CHKERRQ(ierr);
  if(opflow->Nconineq) {
    ierr = VecSet(opflow->Lambdai,1.0);CHKERRQ(ierr);
  }

  /* Only for debugging */
  /*
  ierr = (*opflow->formops.computeequalityconstraints)(opflow,opflow->X,opflow->Ge);CHKERRQ(ierr);
  ierr = (*opflow->formops.computeinequalityconstraints)(opflow,opflow->X,opflow->Gi);CHKERRQ(ierr);
  ierr = VecView(opflow->Gi,0);
  exit(1);
  */

  /* Solve */
  ierr = (*opflow->solverops.solve)(opflow);CHKERRQ(ierr);

  //  ierr = VecView(opflow->X,0);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}




