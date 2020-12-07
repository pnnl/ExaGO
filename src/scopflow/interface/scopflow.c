#include <private/opflowimpl.h>
#include <private/tcopflowimpl.h>
#include <private/scopflowimpl.h>

/*
  SCOPFLOWCreate - Creates a security constrained optimal power flow application object

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

  scopflow->Nconeq   = scopflow->nconeq   = 0;
  scopflow->Nconineq = scopflow->nconineq = 0;
  scopflow->Ncon     = scopflow->ncon     = 0;
  scopflow->Nx       = scopflow->nx       = 0;
  scopflow->Nc       = scopflow->nc       = 0;
  scopflow->Gi       = NULL;
  scopflow->Lambdai  = NULL;
  
  scopflow->obj_factor = 1.0;
  scopflow->obj = 0.0;

  scopflow->ismultiperiod = PETSC_FALSE;

  scopflow->makeup_power_source = 1;

  scopflow->solver   = NULL;
  scopflow->model    = NULL;

  scopflow->nmodelsregistered = 0;
  scopflow->SCOPFLOWModelRegisterAllCalled = PETSC_FALSE;

  /* Register all models */
  ierr = SCOPFLOWModelRegisterAll(scopflow);

  scopflow->nsolversregistered = 0;
  scopflow->SCOPFLOWSolverRegisterAllCalled = PETSC_FALSE;

  /* Register all solvers */
  ierr = SCOPFLOWSolverRegisterAll(scopflow);

  /* Run-time options */
  scopflow->iscoupling = PETSC_FALSE;
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
  PetscInt       c;

  PetscFunctionBegin;
  ierr = COMMDestroy(&(*scopflow)->comm);CHKERRQ(ierr);

  /* Solution vector */
  ierr = VecDestroy(&(*scopflow)->X);CHKERRQ(ierr);
  ierr = VecDestroy(&(*scopflow)->gradobj);CHKERRQ(ierr);

  /* Lower and upper bounds on X */
  ierr = VecDestroy(&(*scopflow)->Xl);CHKERRQ(ierr);
  ierr = VecDestroy(&(*scopflow)->Xu);CHKERRQ(ierr);

  /* Constraints vector */
  ierr = VecDestroy(&(*scopflow)->G);CHKERRQ(ierr);

  ierr = VecDestroy(&(*scopflow)->Gl);CHKERRQ(ierr);
  ierr = VecDestroy(&(*scopflow)->Gu);CHKERRQ(ierr);
  ierr = VecDestroy(&(*scopflow)->Lambda);CHKERRQ(ierr);

  /* Jacobian of constraints */
  ierr = MatDestroy(&(*scopflow)->Jac);CHKERRQ(ierr);
  ierr = MatDestroy(&(*scopflow)->Hes);CHKERRQ(ierr);

  if((*scopflow)->solverops.destroy) {
    ierr = ((*scopflow)->solverops.destroy)(*scopflow);
  }

  if((*scopflow)->modelops.destroy) {
    ierr = ((*scopflow)->modelops.destroy)(*scopflow);
  }

  /* Destroy OPFLOW objects */
  for(c=0; c < (*scopflow)->nc; c++) {
    if(!(*scopflow)->ismultiperiod) {
      ierr = OPFLOWDestroy(&(*scopflow)->opflows[c]);CHKERRQ(ierr);
    } else {
      ierr = TCOPFLOWDestroy(&(*scopflow)->tcopflows[c]);CHKERRQ(ierr);
    }
  }

  ierr = PetscFree((*scopflow)->xstarti);CHKERRQ(ierr);
  ierr = PetscFree((*scopflow)->gstarti);CHKERRQ(ierr);
  ierr = PetscFree((*scopflow)->nxi);CHKERRQ(ierr);
  ierr = PetscFree((*scopflow)->ngi);CHKERRQ(ierr);
  ierr = PetscFree((*scopflow)->nconineqcoup);CHKERRQ(ierr);
  ierr = PetscFree((*scopflow)->ctgclist.cont);CHKERRQ(ierr);
  ierr = PetscFree((*scopflow)->opflows);CHKERRQ(ierr);
  ierr = PetscFree(*scopflow);CHKERRQ(ierr);
  //  *scopflow = 0;
  PetscFunctionReturn(0);
}

/*
  SCOPFLOWEnableMultiPeriod - Enable/Disable multi-period SCOPLOW

  Input Parameters:
+ scopflow - scopflow application object
- ismultperiod - PETSC_FALSE for single-period, PETSC_TRUE otherwise
*/
PetscErrorCode SCOPFLOWEnableMultiPeriod(SCOPFLOW scopflow,PetscBool ismultiperiod)
{
  PetscFunctionBegin;
  scopflow->ismultiperiod = ismultiperiod;
  PetscFunctionReturn(0);
}

/*
  SCOPFLOWSetModel - Sets the model for SCOPFLOW

  Input Parameters:
+ scopflow - scopflow application object
- modelname - name of the model
*/
PetscErrorCode SCOPFLOWSetModel(SCOPFLOW scopflow,const char* modelname)
{
  PetscErrorCode ierr,(*r)(SCOPFLOW)=NULL;
  PetscInt       i;
  PetscFunctionBegin;
  PetscBool match;
  for(i=0;i < scopflow->nmodelsregistered;i++) {
    ierr = PetscStrcmp(scopflow->SCOPFLOWModelList[i].name,modelname,&match);CHKERRQ(ierr);
    if(match) {
      r = scopflow->SCOPFLOWModelList[i].create;
      break;
    }
  }

  if(!r) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown type for SCOPFLOW Model %s",modelname);

  /* Null the function pointers */
  scopflow->modelops.destroy                        = 0;
  scopflow->modelops.setup                          = 0;
  scopflow->modelops.setnumvariablesandconstraints  = 0;
  scopflow->modelops.setvariablebounds              = 0;
  scopflow->modelops.setconstraintbounds            = 0;
  scopflow->modelops.setvariableandconstraintbounds = 0;
  scopflow->modelops.setinitialguess                = 0;
  scopflow->modelops.computeconstraints             = 0;
  scopflow->modelops.computejacobian                = 0;
  scopflow->modelops.computehessian                 = 0;
  scopflow->modelops.computeobjandgradient          = 0;
  scopflow->modelops.computeobjective               = 0;
  scopflow->modelops.computegradient                = 0;

  ierr = PetscStrcpy(scopflow->modelname,modelname);CHKERRQ(ierr);
  /* Call the underlying implementation constructor */
  ierr = (*r)(scopflow);CHKERRQ(ierr);

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

  ierr = PetscStrcpy(scopflow->solvername,solvername);CHKERRQ(ierr);
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
  PetscBool      scopflowsolverset;
  char           opflowmodelname[32]="POWER_BALANCE_POLAR";
  char           scopflowsolvername[32]="IPOPTNEW";
  PetscInt       c,i,j;
  char           scopflowmodelname[32]="GENRAMP";
  PS             ps;
  OPFLOW         opflow;
  char           ploadprofile[PETSC_MAX_PATH_LEN];
  char           qloadprofile[PETSC_MAX_PATH_LEN];
  char           windgenprofile[PETSC_MAX_PATH_LEN];
  PetscBool      flg=PETSC_FALSE,flg1=PETSC_FALSE,flg2=PETSC_FALSE,flg3=PETSC_FALSE;
  PetscReal      dT=5.0,duration=0.0; /* 5 minute time-step and single period */

  PetscFunctionBegin;

  ierr = PetscOptionsBegin(scopflow->comm->type,NULL,"SCOPFLOW options",NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-scopflow_model","SCOPFLOW model type","",scopflowmodelname,scopflowmodelname,32,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-scopflow_solver","SCOPFLOW solver type","",scopflowsolvername,scopflowsolvername,32,&scopflowsolverset);CHKERRQ(ierr);

  ierr = PetscOptionsBool("-scopflow_iscoupling","Include coupling between first stage and second stage","",scopflow->iscoupling,&scopflow->iscoupling,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-scopflow_Nc","Number of second stage scenarios","",scopflow->Nc,&scopflow->Nc,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-scopflow_replicate_basecase","Only for debugging: Replicate first stage for all second stage scenarios","",scopflow->replicate_basecase,&scopflow->replicate_basecase,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-scopflow_mode","Operation mode:Preventive (0) or Corrective (1)","",scopflow->mode,&scopflow->mode,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-scopflow_makeup_power_source","Make-up power source for contingencies","",scopflow->makeup_power_source,&scopflow->makeup_power_source,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-scopflow_enable_multiperiod","Multi-period SCOPFLOW?","",scopflow->ismultiperiod,&scopflow->ismultiperiod,NULL);CHKERRQ(ierr);
  if(scopflow->ismultiperiod) {
    /* Set loadp,loadq, and windgen files */
    ierr = PetscOptionsString("-scopflow_ploadprofile","Active power load profile","",ploadprofile,ploadprofile,100,&flg1);CHKERRQ(ierr);
    ierr = PetscOptionsString("-scopflow_qloadprofile","Reactive power load profile","",qloadprofile,qloadprofile,100,&flg2);CHKERRQ(ierr);
    ierr = PetscOptionsString("-scopflow_windgenprofile","Wind generation profile","",windgenprofile,windgenprofile,100,&flg3);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-scopflow_dT","Length of time-step (minutes)","",dT,&dT,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-scopflow_duration","Time horizon (hours)","",duration,&duration,NULL);CHKERRQ(ierr);

  }

  PetscOptionsEnd();

  if(scopflow->ctgcfileset && !scopflow->replicate_basecase) {
    if(scopflow->Nc < 0) scopflow->Nc = MAX_CONTINGENCIES;
    else scopflow->Nc += 1; 

    ierr = PetscCalloc1(scopflow->Nc,&scopflow->ctgclist.cont);CHKERRQ(ierr);
    for(c=0; c < scopflow->Nc; c++) scopflow->ctgclist.cont->noutages = 0;
    if(scopflow->Nc > 1) { 
      ierr = SCOPFLOWReadContingencyData(scopflow,scopflow->ctgcfileformat,scopflow->ctgcfile);CHKERRQ(ierr);
      scopflow->Nc = scopflow->ctgclist.Ncont+1;
    }
  } else {
    if(scopflow->Nc == -1) scopflow->Nc = 1;
  }

  int q = scopflow->Nc/scopflow->comm->size;
  int d = scopflow->Nc%scopflow->comm->size;
  if(d) {
    scopflow->nc = q + ((scopflow->comm->rank < d)?1:0); 
  } else {
    scopflow->nc = q;
  }
  ierr = MPI_Scan(&scopflow->nc,&scopflow->cend,1,MPIU_INT,MPI_SUM,scopflow->comm->type);CHKERRQ(ierr);
  scopflow->cstart = scopflow->cend - scopflow->nc;

  ierr = PetscPrintf(scopflow->comm->type,"SCOPFLOW running with %d contingencies (base case + %d contingencies)\n",scopflow->Nc,scopflow->Nc-1);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF,"Rank %d has %d contingencies, range [%d -- %d]\n",scopflow->comm->rank,scopflow->nc,scopflow->cstart,scopflow->cend);CHKERRQ(ierr);

  /* Set Model */
  if(!scopflow->ismultiperiod) {
    ierr = SCOPFLOWSetModel(scopflow,scopflowmodelname);CHKERRQ(ierr);
  } else {
    ierr = SCOPFLOWSetModel(scopflow,"GENRAMPT");CHKERRQ(ierr);
  }

  /* Set solver */
  if(scopflowsolverset) {
    if(scopflow->solver) ierr = (*scopflow->solverops.destroy)(scopflow);
    ierr = SCOPFLOWSetSolver(scopflow,scopflowsolvername);CHKERRQ(ierr);
    ierr = PetscPrintf(scopflow->comm->type,"SCOPFLOW: Using %s solver\n",scopflowsolvername);CHKERRQ(ierr);
  } else {
    if(!scopflow->solver) {
      ierr = SCOPFLOWSetSolver(scopflow,SCOPFLOWSOLVER_IPOPTNEW);CHKERRQ(ierr);
      ierr = PetscPrintf(scopflow->comm->type,"SCOPFLOW: Using %s solver\n",SCOPFLOWSOLVER_IPOPT);CHKERRQ(ierr); 
    }
  }
  
  if(!scopflow->ismultiperiod) {
    /* Create OPFLOW objects */
    ierr = PetscCalloc1(scopflow->nc,&scopflow->opflows);CHKERRQ(ierr);
    for(c=0; c < scopflow->nc; c++) {
      ierr = OPFLOWCreate(PETSC_COMM_SELF,&scopflow->opflows[c]);CHKERRQ(ierr);
      ierr = OPFLOWSetModel(scopflow->opflows[c],OPFLOWMODEL_PBPOL);CHKERRQ(ierr);
      //    ierr = OPFLOWSetSolver(scopflow->opflows[c],opflowsolvername);CHKERRQ(ierr);
      
      ierr = OPFLOWReadMatPowerData(scopflow->opflows[c],scopflow->netfile);CHKERRQ(ierr);
      /* Set up the PS object for opflow */
      ps = scopflow->opflows[c]->ps;
      ierr = PSSetUp(ps);CHKERRQ(ierr);
      
      /* Set contingencies */
      if(scopflow->ctgcfileset && !scopflow->replicate_basecase) {
	Contingency ctgc=scopflow->ctgclist.cont[scopflow->cstart+c];
	for(j=0; j < ctgc.noutages; j++) {
	  if(ctgc.outagelist[j].type == GEN_OUTAGE) {
	  PetscInt gbus=ctgc.outagelist[j].bus;
	  char     *gid = ctgc.outagelist[j].id;
	  PetscInt status = ctgc.outagelist[j].status;
	  ierr = PSSetGenStatus(ps,gbus,gid,status);CHKERRQ(ierr);
	  }
	  if(ctgc.outagelist[j].type == BR_OUTAGE) {
	    PetscInt fbus=ctgc.outagelist[j].fbus;
	    PetscInt tbus=ctgc.outagelist[j].tbus;
	    char     *brid = ctgc.outagelist[j].id;
	    PetscInt status = ctgc.outagelist[j].status;
	    ierr = PSSetLineStatus(ps,fbus,tbus,brid,status);CHKERRQ(ierr);
	  }
	}
      }
      
      /* Set up OPFLOW object */
      //    ierr = OPFLOWGenbusVoltageFixed(scopflow->opflows[c],PETSC_TRUE);CHKERRQ(ierr);
      //    if(scopflow->cstart+c > 0) scopflow->opflows[c]->obj_gencost = PETSC_FALSE; /* No gen. cost minimization for second stage */
      ierr = OPFLOWSetUp(scopflow->opflows[c]);CHKERRQ(ierr);
    }
  } else {
    TCOPFLOW          tcopflow;
    
    /* Create TCOPFLOW objects */
    ierr = PetscCalloc1(scopflow->nc,&scopflow->tcopflows);CHKERRQ(ierr);
    for(c=0; c < scopflow->nc; c++) {
      ierr = TCOPFLOWCreate(PETSC_COMM_SELF,&scopflow->tcopflows[c]);CHKERRQ(ierr);
      tcopflow = scopflow->tcopflows[c];

      ierr = TCOPFLOWSetNetworkData(tcopflow,scopflow->netfile);CHKERRQ(ierr);

      if(flg1 && flg2) {
	ierr = TCOPFLOWSetLoadProfiles(tcopflow,ploadprofile,qloadprofile);CHKERRQ(ierr);
      } else if(flg1) {
	ierr = TCOPFLOWSetLoadProfiles(tcopflow,ploadprofile,NULL);CHKERRQ(ierr);
      } else if(flg2) {
	ierr = TCOPFLOWSetLoadProfiles(tcopflow,NULL,qloadprofile);CHKERRQ(ierr);
      }
      if(flg3) {
	ierr = TCOPFLOWSetWindGenProfiles(tcopflow,windgenprofile);CHKERRQ(ierr);
      }

      ierr = TCOPFLOWSetTimeStepandDuration(tcopflow,dT,duration);CHKERRQ(ierr);

      /* Set contingencies */
      if(scopflow->ctgcfileset && !scopflow->replicate_basecase) {
	Contingency ctgc=scopflow->ctgclist.cont[scopflow->cstart+c];
	// Set this contingency with TCOPFLOW
	ierr = TCOPFLOWSetContingency(tcopflow,&ctgc);CHKERRQ(ierr);
      }
      
      ierr = TCOPFLOWSetUp(tcopflow);CHKERRQ(ierr);
    }
  }    
  
  ierr = PetscCalloc1(scopflow->nc,&scopflow->nconineqcoup);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->nc,&scopflow->nxi);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->nc,&scopflow->ngi);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->nc,&scopflow->xstarti);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->nc,&scopflow->gstarti);CHKERRQ(ierr);

  /* Set number of variables and constraints */
  ierr = (*scopflow->modelops.setnumvariablesandconstraints)(scopflow,scopflow->nxi,scopflow->ngi,scopflow->nconineqcoup);

  if(!scopflow->ismultiperiod) {
    scopflow->nx = scopflow->nxi[0];
    scopflow->ncon = scopflow->ngi[0];
    scopflow->nconcoup = scopflow->nconineqcoup[0];
    opflow = scopflow->opflows[0];
    scopflow->nconeq = opflow->nconeq;
    scopflow->nconineq = opflow->nconineq;
    
    for(i=1; i < scopflow->nc; i++) {
      scopflow->xstarti[i] = scopflow->xstarti[i-1] + scopflow->nxi[i-1];
      scopflow->gstarti[i] = scopflow->gstarti[i-1] + scopflow->ngi[i-1];
      scopflow->nx += scopflow->nxi[i];
      scopflow->ncon += scopflow->ngi[i];
      scopflow->nconcoup += scopflow->nconineqcoup[i];
      opflow = scopflow->opflows[i];
      scopflow->nconeq   += opflow->nconeq;
      scopflow->nconineq += opflow->nconineq;
    }
  } else {
    TCOPFLOW tcopflow;
    scopflow->nx = scopflow->nxi[0];
    scopflow->ncon = scopflow->ngi[0];
    scopflow->nconcoup = scopflow->nconineqcoup[0];
    tcopflow = scopflow->tcopflows[0];
    scopflow->nconeq = tcopflow->Nconeq;
    scopflow->nconineq = tcopflow->Nconineq + tcopflow->Nconcoup;
    
    for(i=1; i < scopflow->nc; i++) {
      scopflow->xstarti[i] = scopflow->xstarti[i-1] + scopflow->nxi[i-1];
      scopflow->gstarti[i] = scopflow->gstarti[i-1] + scopflow->ngi[i-1];
      scopflow->nx += scopflow->nxi[i];
      scopflow->ncon += scopflow->ngi[i];
      scopflow->nconcoup += scopflow->nconineqcoup[i];
      tcopflow = scopflow->tcopflows[i];
      scopflow->nconeq   += tcopflow->Nconeq;
      scopflow->nconineq += tcopflow->Nconineq + tcopflow->Nconcoup;
    }
  }

  /* Create vector X */
  ierr = VecCreate(scopflow->comm->type,&scopflow->X);CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->X,scopflow->nx,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(scopflow->X);CHKERRQ(ierr);
  ierr = VecGetSize(scopflow->X,&scopflow->Nx);CHKERRQ(ierr);

  ierr = VecDuplicate(scopflow->X,&scopflow->Xl);CHKERRQ(ierr);
  ierr = VecDuplicate(scopflow->X,&scopflow->Xu);CHKERRQ(ierr);
  ierr = VecDuplicate(scopflow->X,&scopflow->gradobj);CHKERRQ(ierr);

  /* vector for constraints */
  ierr = VecCreate(scopflow->comm->type,&scopflow->G);CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->G,scopflow->ncon,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(scopflow->G);CHKERRQ(ierr);
  ierr = VecGetSize(scopflow->G,&scopflow->Ncon);CHKERRQ(ierr);

  /* Constraint bounds vectors  */
  ierr = VecDuplicate(scopflow->G,&scopflow->Gl);CHKERRQ(ierr);
  ierr = VecDuplicate(scopflow->G,&scopflow->Gu);CHKERRQ(ierr);

  /* Constraint Jacobian */
  ierr = MatCreate(scopflow->comm->type,&scopflow->Jac);CHKERRQ(ierr);
  ierr = MatSetSizes(scopflow->Jac,scopflow->ncon,scopflow->nx,scopflow->Ncon,scopflow->Nx);CHKERRQ(ierr);
  ierr = MatSetUp(scopflow->Jac);CHKERRQ(ierr);
  ierr = MatSetFromOptions(scopflow->Jac);CHKERRQ(ierr);

  /* Hessian */
  ierr = MatCreate(scopflow->comm->type,&scopflow->Hes);CHKERRQ(ierr);
  ierr = MatSetSizes(scopflow->Hes,scopflow->nx,scopflow->nx,scopflow->Nx,scopflow->Nx);CHKERRQ(ierr);
  ierr = MatSetUp(scopflow->Hes);CHKERRQ(ierr);
  ierr = MatSetFromOptions(scopflow->Hes);CHKERRQ(ierr);

  /* Lagrangian multipliers */
  ierr = VecDuplicate(scopflow->G,&scopflow->Lambda);CHKERRQ(ierr);

  ierr = (*scopflow->solverops.setup)(scopflow);CHKERRQ(ierr);
  ierr = PetscPrintf(scopflow->comm->type,"SCOPFLOW: Setup completed\n");CHKERRQ(ierr);
  
  scopflow->setupcalled = PETSC_TRUE;
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
  PetscBool      issolver_empar;

  PetscFunctionBegin;

  if(!scopflow->setupcalled) {
    ierr = SCOPFLOWSetUp(scopflow);
  }

  ierr = PetscStrcmp(scopflow->solvername,"EMPAR",&issolver_empar);CHKERRQ(ierr);

  if(!issolver_empar) { /* Don't need to do all this for embarassingly parallel solver */
    /* Set bounds on variables */
    if(scopflow->modelops.setvariablebounds) {
      ierr = (*scopflow->modelops.setvariablebounds)(scopflow,scopflow->Xl,scopflow->Xu);CHKERRQ(ierr);
    }
    
    
    /* Set bounds on constraints */
    if(scopflow->modelops.setconstraintbounds) {
      ierr = (*scopflow->modelops.setconstraintbounds)(scopflow,scopflow->Gl,scopflow->Gu);CHKERRQ(ierr);
    }
    
    /* Set initial guess */
    if(scopflow->modelops.setinitialguess) {
      ierr = (*scopflow->modelops.setinitialguess)(scopflow,scopflow->X);CHKERRQ(ierr);
    }

    ierr = VecSet(scopflow->Lambda,1.0);CHKERRQ(ierr);
  }

  /* Solve */
  ierr = (*scopflow->solverops.solve)(scopflow);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  SCOPFLOWGetObjective - Returns the objective function value

  Input Parameters:
+ SCOPFLOW - the SCOPFLOW object
- obj    - the objective function value

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode SCOPFLOWGetObjective(SCOPFLOW scopflow,PetscReal *obj)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*scopflow->solverops.getobjective)(scopflow,obj);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SCOPFLOWGetBaseSolution - Returns the SCOPFLOW solution for a given contingency

  Input Parameters:
+ SCOPFLOW - the SCOPFLOW object
. contnum  - Contingency number (0 for base/no-contingency)
- X        - the scopflow solution for the given contingency

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode SCOPFLOWGetSolution(SCOPFLOW scopflow,PetscInt contnum,Vec *X)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*scopflow->solverops.getsolution)(scopflow,contnum,X);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SCOPFLOWGetConstraints - Returns the SCOPFLOW constraints for a given contingency

  Input Parameters:
+ SCOPFLOW - the SCOPFLOW object
. contnum  - Contingency number (0 for base/no-contingency)
- G    - the scopflow constraints

  Notes: Should be called after the optimization finishes.
         Equality constraints first followed by inequality constraints
*/
PetscErrorCode SCOPFLOWGetConstraints(SCOPFLOW scopflow,PetscInt contnum,Vec *G)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*scopflow->solverops.getconstraints)(scopflow,contnum,G);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SCOPFLOWGetConstraintMultipliers - Returns the SCOPFLOW constraint multipliers for a given contingency

  Input Parameters:
+ SCOPFLOW - the SCOPFLOW object
. contnum  - Contingency number (0 for base/no-contingency)
- G    - the scopflow constraint lagrange multipliers

  Notes: Should be called after the optimization finishes.
    Equality constraint multipliers first followed by inequality constraint multipliers
*/
PetscErrorCode SCOPFLOWGetConstraintMultipliers(SCOPFLOW scopflow,PetscInt contnum,Vec *Lambda)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*scopflow->solverops.getconstraintmultipliers)(scopflow,contnum,Lambda);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SCOPFLOWGetConvergenceStatus - Did SCOPFLOW converge?

  Input Parameters:
+ SCOPFLOW - the SCOPFLOW object
- status - PETSC_TRUE if converged, PETSC_FALSE otherwise

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode SCOPFLOWGetConvergenceStatus(SCOPFLOW scopflow,PetscBool *status)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*scopflow->solverops.getconvergencestatus)(scopflow,status);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
