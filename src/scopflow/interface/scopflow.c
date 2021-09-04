#include <private/opflowimpl.h>
#include <private/tcopflowimpl.h>
#include <private/scopflowimpl.h>
#include <utils.hpp>

/**
 * @brief Creates a security constrained optimal power flow application object
 * 
 * @param[in] mpicomm the MPI communicator
 * @param[in] scopflowout pointer to the security constrained optimal power flow application object
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

  scopflow->mode = 0;
  scopflow->tolerance = 1e-6;

  scopflow->scen = NULL;

  scopflow->ismultiperiod = PETSC_FALSE;

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

  scopflow->ctgclist    = NULL;
  scopflow->ctgcfileset = PETSC_FALSE;

  scopflow->setupcalled = PETSC_FALSE;
  *scopflowout = scopflow;

  ExaGOLog(EXAGO_LOG_INFO,"%s","SCOPFLOW: Application created\n");
  PetscFunctionReturn(0);
}

/**
 * @brief Destroys the security constrained optimal power flow application object
 * @param[in] scopflow the scopflow object to destroy
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

  /* Destroy TCOPFLOW or OPFLOW objects */
  if(!(*scopflow)->ismultiperiod) {
    ierr = OPFLOWDestroy(&(*scopflow)->opflow0);CHKERRQ(ierr);
    for(c=0; c < (*scopflow)->nc; c++) {
      ierr = OPFLOWDestroy(&(*scopflow)->opflows[c]);CHKERRQ(ierr);
    }
    ierr = PetscFree((*scopflow)->opflows);CHKERRQ(ierr);
  } else {
    for(c=0; c < (*scopflow)->nc; c++) {
      ierr = TCOPFLOWDestroy(&(*scopflow)->tcopflows[c]);CHKERRQ(ierr);
    }
    ierr = PetscFree((*scopflow)->tcopflows);CHKERRQ(ierr);
  }

  ierr = PetscFree((*scopflow)->xstarti);CHKERRQ(ierr);
  ierr = PetscFree((*scopflow)->gstarti);CHKERRQ(ierr);
  ierr = PetscFree((*scopflow)->nxi);CHKERRQ(ierr);
  ierr = PetscFree((*scopflow)->ngi);CHKERRQ(ierr);
  ierr = PetscFree((*scopflow)->nconeqcoup);CHKERRQ(ierr);
  ierr = PetscFree((*scopflow)->nconineqcoup);CHKERRQ(ierr);
  if((*scopflow)->Nc > 1) {
    ierr = ContingencyListDestroy(&(*scopflow)->ctgclist);CHKERRQ(ierr);
  }
  ierr = PetscFree(*scopflow);CHKERRQ(ierr);
  //  *scopflow = 0;
  PetscFunctionReturn(0);
}

/**
 * @brief Enable/Disable multi-period SCOPLOW
 * 
 * @param[in] scopflow the scopflow application object
 * @param[in] ismultperiod PETSC_FALSE for single-period, PETSC_TRUE otherwise
*/
PetscErrorCode SCOPFLOWEnableMultiPeriod(SCOPFLOW scopflow,PetscBool ismultiperiod)
{
  PetscFunctionBegin;
  scopflow->ismultiperiod = ismultiperiod;
  PetscFunctionReturn(0);
}

/**
 * @brief Sets the model for SCOPFLOW
 * 
 * @param[in] scopflow scopflow application object
 * @param[in] modelname name of the model
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

/**
 * @brief Sets the solver for SCOPFLOW
 * 
 * @param[in] scopflow the scopflow application object
 * @param[in] solvername name of the solver
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

/**
 * @brief Sets and reads the network data
 * 
 * @param[in] scopflow the scopflow object
 * @param[in] netfile the name of the network file
 * 
 * @note The input data must be in MATPOWER data format
*/
PetscErrorCode SCOPFLOWSetNetworkData(SCOPFLOW scopflow,const char netfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscMemcpy(scopflow->netfile,netfile,PETSC_MAX_PATH_LEN*sizeof(char));CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/**
 * @brief Sets the pload data
 * 
 * @param[in] scopflow the scopflow object
 * @param[in] profile the name of the pload file
 * 
 * @note The input data must be in MATPOWER data format
*/
PetscErrorCode SCOPFLOWSetPLoadData(SCOPFLOW scopflow,const char profile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscMemcpy(scopflow->pload,profile,PETSC_MAX_PATH_LEN*sizeof(char));CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/**
 * @brief Sets the qload data
 * 
 * @param[in] scopflow the scopflow object
 * @param[in] profile the name of the qload file
 * 
 * @note The input data must be in MATPOWER data format
*/
PetscErrorCode SCOPFLOWSetQLoadData(SCOPFLOW scopflow,const char profile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscMemcpy(scopflow->qload,profile,PETSC_MAX_PATH_LEN*sizeof(char));CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/**
 * @brief Sets the WindGen profile data
 * 
 * @param[in] scopflow the scopflow object
 * @param[in] profile the name of the profile 
 * 
 * @note The input data must be in MATPOWER data format
*/
PetscErrorCode SCOPFLOWSetWindGenProfile(SCOPFLOW scopflow,const char profile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscMemcpy(scopflow->windgen,profile,PETSC_MAX_PATH_LEN*sizeof(char));CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* This is kind of a hack to update the variable bounds for OPFLOW based on the mode SCOPFLOW uses
*/
PetscErrorCode SCOPFLOWUpdateOPFLOWVariableBounds(OPFLOW opflow, Vec Xl, Vec Xu,void* ctx)
{
  PetscErrorCode ierr;
  SCOPFLOW       scopflow=(SCOPFLOW)ctx;

  PetscFunctionBegin;
  if(opflow->has_gensetpoint) {
    /* Modify the bounds on ramping variables */
    PetscInt       j,k;
    PS             ps = opflow->ps;
    PSBUS          bus;
    PSGEN          gen;
    PetscScalar    *xl,*xu;
    
    ierr = VecGetArray(Xl,&xl);CHKERRQ(ierr);
    ierr = VecGetArray(Xu,&xu);CHKERRQ(ierr);
    for(j = 0; j < ps->nbus; j++) {
      bus = &ps->bus[j];
      for(k=0; k < bus->ngen; k++) {
	ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	if(!gen->status) continue;
	if(scopflow->mode == 0) {
	  /* Only ref. bus responsible for make-up power for contingencies */
	  if(bus->ide != REF_BUS) {
	    xl[opflow->idxn2sd_map[gen->startxpdevloc]]   = xu[opflow->idxn2sd_map[gen->startxpdevloc]] = 0.0;
	  }
	} else {
	  xl[opflow->idxn2sd_map[gen->startxpdevloc]] = gen->pb - gen->pgs; //-gen->ramp_rate_30min;
	  xu[opflow->idxn2sd_map[gen->startxpdevloc]] =  gen->pt - gen->pgs; //gen->ramp_rate_30min;
	}
      }
    }
    ierr = VecRestoreArray(Xl,&xl);CHKERRQ(ierr);
    ierr = VecRestoreArray(Xu,&xu);CHKERRQ(ierr);
  } 
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
  char           scopflowsolvername[32]="IPOPT";
  PetscInt       c,i,j;
  char           scopflowmodelname[32]="GENRAMP";
  PS             ps;
  OPFLOW         opflow;
  char           ploadprofile[PETSC_MAX_PATH_LEN];
  char           qloadprofile[PETSC_MAX_PATH_LEN];
  char           windgenprofile[PETSC_MAX_PATH_LEN];
  PetscBool      flg1=PETSC_FALSE,flg2=PETSC_FALSE,flg3=PETSC_FALSE;

  PetscFunctionBegin;

  ierr = PetscOptionsBegin(scopflow->comm->type,NULL,"SCOPFLOW options",NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-scopflow_model","SCOPFLOW model type","",scopflowmodelname,scopflowmodelname,32,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-scopflow_solver","SCOPFLOW solver type","",scopflowsolvername,scopflowsolvername,32,&scopflowsolverset);CHKERRQ(ierr);

  ierr = PetscOptionsBool("-scopflow_iscoupling","Include coupling between first stage and second stage","",scopflow->iscoupling,&scopflow->iscoupling,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-scopflow_Nc","Number of second stage scenarios","",scopflow->Nc,&scopflow->Nc,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-scopflow_mode","Operation mode:Preventive (0) or Corrective (1)","",scopflow->mode,&scopflow->mode,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-scopflow_enable_multiperiod","Multi-period SCOPFLOW?","",scopflow->ismultiperiod,&scopflow->ismultiperiod,NULL);CHKERRQ(ierr);
  if(scopflow->ismultiperiod) {
    /* Set loadp,loadq, and windgen files */
    ierr = PetscOptionsString("-scopflow_ploadprofile","Active power load profile","",ploadprofile,ploadprofile,200,&flg1);CHKERRQ(ierr);
    ierr = PetscOptionsString("-scopflow_qloadprofile","Reactive power load profile","",qloadprofile,qloadprofile,200,&flg2);CHKERRQ(ierr);
    ierr = PetscOptionsString("-scopflow_windgenprofile","Wind generation profile","",windgenprofile,windgenprofile,200,&flg3);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-scopflow_dT","Length of time-step (minutes)","",scopflow->dT,&scopflow->dT,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-scopflow_duration","Time horizon (hours)","",scopflow->duration,&scopflow->duration,NULL);CHKERRQ(ierr);

  }
  ierr = PetscOptionsReal("-scopflow_tolerance", "optimization tolerance", "", scopflow->tolerance,&scopflow->tolerance, NULL); CHKERRQ(ierr);
  PetscOptionsEnd();

  if(scopflow->ctgcfileset) {
    if(scopflow->Nc < 0) scopflow->Nc = MAX_CONTINGENCIES;
    else scopflow->Nc += 1; 


    if(scopflow->Nc > 1) { 
      /* Create contingency list object */
      ierr = ContingencyListCreate(scopflow->Nc,&scopflow->ctgclist);CHKERRQ(ierr);
      ierr = ContingencyListSetData(scopflow->ctgclist,scopflow->ctgcfileformat,scopflow->ctgcfile);CHKERRQ(ierr);
      ierr = ContingencyListReadData(scopflow->ctgclist,&scopflow->Nc);CHKERRQ(ierr);
      scopflow->Nc += 1;
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

  
  ExaGOLog(EXAGO_LOG_INFO,"SCOPFLOW running with %d contingencies (base case + %d contingencies)\n",scopflow->Nc,scopflow->Nc-1);
  //  ExaGOLogUseEveryRank(PETSC_TRUE);
  //  ExaGOLog(EXAGO_LOG_INFO,"Rank %d has %d contingencies, range [%d -- %d]\n",scopflow->comm->rank,scopflow->nc,scopflow->cstart,scopflow->cend);
  //  ExaGOLogUseEveryRank(PETSC_FALSE);

  /* Set model */
  if(!scopflow->ismultiperiod) {
    ierr = SCOPFLOWSetModel(scopflow,scopflowmodelname);CHKERRQ(ierr);
  } else {
    ierr = SCOPFLOWSetModel(scopflow,"GENRAMPT");CHKERRQ(ierr);
  }

  /* Set solver */
  if(scopflowsolverset) {
    if(scopflow->solver) ierr = (*scopflow->solverops.destroy)(scopflow);
    ierr = SCOPFLOWSetSolver(scopflow,scopflowsolvername);CHKERRQ(ierr);
    ExaGOLog(EXAGO_LOG_INFO,"SCOPFLOW: Using %s solver\n",scopflowsolvername);
  } else {
    if(!scopflow->solver) {
      ierr = SCOPFLOWSetSolver(scopflow,SCOPFLOWSOLVER_IPOPT);CHKERRQ(ierr);
      ExaGOLog(EXAGO_LOG_INFO,"SCOPFLOW: Using %s solver\n",SCOPFLOWSOLVER_IPOPT);
    }
  }
  
  if(!scopflow->ismultiperiod) {
    /* Create OPFLOW objects */
    ierr = PetscCalloc1(scopflow->nc,&scopflow->opflows);CHKERRQ(ierr);

    /* Create base-case OPFLOW */
    ierr = OPFLOWCreate(PETSC_COMM_SELF,&scopflow->opflow0);CHKERRQ(ierr);
    ierr = OPFLOWSetModel(scopflow->opflow0,OPFLOWMODEL_PBPOL);CHKERRQ(ierr);
    ierr = OPFLOWReadMatPowerData(scopflow->opflow0,scopflow->netfile);CHKERRQ(ierr);
    ierr = OPFLOWSetInitializationType(scopflow->opflow0, scopflow->type);CHKERRQ(ierr);
    ierr = OPFLOWSetGenBusVoltageType(scopflow->opflow0, scopflow->genbusvoltagetype);CHKERRQ(ierr);

    ierr = PSSetUp(scopflow->opflow0->ps);CHKERRQ(ierr);
    if(scopflow->scen) {
      // Scenario set by SOPFLOW, apply it */
      ierr = PSApplyScenario(scopflow->opflow0->ps,*scopflow->scen);CHKERRQ(ierr);
    }
    ierr = OPFLOWSetUp(scopflow->opflow0);CHKERRQ(ierr);

    ierr = OPFLOWSetObjectiveType(scopflow->opflow0,MIN_GEN_COST);CHKERRQ(ierr);


    for(c=0; c < scopflow->nc; c++) {
      ierr = OPFLOWCreate(PETSC_COMM_SELF,&scopflow->opflows[c]);CHKERRQ(ierr);
      //      ierr = OPFLOWSetModel(scopflow->opflows[c],OPFLOWMODEL_PBPOL);CHKERRQ(ierr);
      //    ierr = OPFLOWSetSolver(scopflow->opflows[c],opflowsolvername);CHKERRQ(ierr);
      ierr = OPFLOWSetInitializationType(scopflow->opflows[c], scopflow->type);CHKERRQ(ierr);
      ierr = OPFLOWSetGenBusVoltageType(scopflow->opflows[c], scopflow->genbusvoltagetype);CHKERRQ(ierr);
      
      ierr = OPFLOWReadMatPowerData(scopflow->opflows[c],scopflow->netfile);CHKERRQ(ierr);
      /* Set up the PS object for opflow */
      ps = scopflow->opflows[c]->ps;
      ierr = PSSetUp(ps);CHKERRQ(ierr);

      if(scopflow->scen) {
	// Scenario set by SOPFLOW, apply it */
	ierr = PSApplyScenario(ps,*scopflow->scen);CHKERRQ(ierr);
      }
      
      /* Set contingencies */
      if(scopflow->ctgcfileset &&scopflow->Nc > 1) {
	Contingency ctgc=scopflow->ctgclist->cont[scopflow->cstart+c];
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
      
      if(scopflow->cstart+c ==  0) { /* First stage */
	ierr = OPFLOWSetObjectiveType(scopflow->opflows[c],MIN_GEN_COST);CHKERRQ(ierr);
      } else { /* Second stages */
	ierr = OPFLOWHasGenSetPoint(scopflow->opflows[c],PETSC_TRUE);CHKERRQ(ierr); /* Activates ramping variables */
	ierr = OPFLOWSetModel(scopflow->opflows[c],OPFLOWMODEL_PBPOL);CHKERRQ(ierr);
	ierr = OPFLOWSetObjectiveType(scopflow->opflows[c],NO_OBJ);CHKERRQ(ierr);
	ierr = OPFLOWSetUpdateVariableBoundsFunction(scopflow->opflows[c],SCOPFLOWUpdateOPFLOWVariableBounds,(void*)scopflow);
      }

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
      // CLI option overrides any existing setting in scopflow
      if(flg1)
      {
        SCOPFLOWSetPLoadData(scopflow, ploadprofile);
      }
      if(flg2)
      {
        SCOPFLOWSetQLoadData(scopflow, qloadprofile);
      }
      if(flg3)
      {
        SCOPFLOWSetWindGenProfile(scopflow,windgenprofile);
      }

    	ierr = TCOPFLOWSetLoadProfiles(tcopflow,scopflow->pload,scopflow->qload);CHKERRQ(ierr);
	    ierr = TCOPFLOWSetWindGenProfiles(tcopflow,scopflow->windgen);CHKERRQ(ierr);
      ierr = TCOPFLOWSetTimeStepandDuration(tcopflow,scopflow->dT,scopflow->duration);CHKERRQ(ierr);

      if(scopflow->scen) {
	/* SOPFLOW has set scenario, pass it to TCOPFLOW */
	ierr = TCOPFLOWSetScenario(tcopflow,scopflow->scen);CHKERRQ(ierr);
      }

      /* Set contingencies */
      if(scopflow->ctgcfileset) {
	Contingency ctgc=scopflow->ctgclist->cont[scopflow->cstart+c];
	// Set this contingency with TCOPFLOW
	ierr = TCOPFLOWSetContingency(tcopflow,&ctgc);CHKERRQ(ierr);
      }
      
      ierr = TCOPFLOWSetUp(tcopflow);CHKERRQ(ierr);
    }
  }    

  ierr = PetscCalloc1(scopflow->nc,&scopflow->nconeqcoup);CHKERRQ(ierr);  
  ierr = PetscCalloc1(scopflow->nc,&scopflow->nconineqcoup);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->nc,&scopflow->nxi);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->nc,&scopflow->ngi);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->nc,&scopflow->xstarti);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->nc,&scopflow->gstarti);CHKERRQ(ierr);

  /* Set number of variables and constraints */
  ierr = (*scopflow->modelops.setnumvariablesandconstraints)(scopflow,scopflow->nxi,scopflow->ngi,scopflow->nconeqcoup,scopflow->nconineqcoup);

  if(!scopflow->ismultiperiod) {
    scopflow->nx = scopflow->nxi[0];
    scopflow->ncon = scopflow->ngi[0];
    scopflow->nconcoup = scopflow->nconeqcoup[0] + scopflow->nconineqcoup[0];
    opflow = scopflow->opflows[0];
    scopflow->nconeq = opflow->nconeq;
    scopflow->nconineq = opflow->nconineq;
    
    for(i=1; i < scopflow->nc; i++) {
      scopflow->xstarti[i] = scopflow->xstarti[i-1] + scopflow->nxi[i-1];
      scopflow->gstarti[i] = scopflow->gstarti[i-1] + scopflow->ngi[i-1];
      scopflow->nx += scopflow->nxi[i];
      scopflow->ncon += scopflow->ngi[i];
      scopflow->nconcoup += scopflow->nconeqcoup[i] + scopflow->nconineqcoup[i];
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

  /* Vector for constraints */
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
  ExaGOLog(EXAGO_LOG_INFO,"%s","SCOPFLOW: Setup completed\n");
  
  scopflow->setupcalled = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/**
 * @brief Solves the AC security constrained optimal power flow
 *
 * @param[in] scopflow the security constrained optimal power flow application object
*/
PetscErrorCode SCOPFLOWSolve(SCOPFLOW scopflow)
{
  PetscErrorCode ierr;
  PetscBool      issolver_ipopt;

  PetscFunctionBegin;

  if(!scopflow->setupcalled) {
    ierr = SCOPFLOWSetUp(scopflow);
  }

  ierr = PetscStrcmp(scopflow->solvername,"IPOPT",&issolver_ipopt);CHKERRQ(ierr);

  if(issolver_ipopt) { /* Don't need to do this if solver is not ipopt */
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

/**
 * @brief Returns the objective function value
 * 
 * @param[in] scopflow the scopflow object
 * @param[in] obj the objective function value
 * 
 * @note Should be called after the optimization finishes
*/
PetscErrorCode SCOPFLOWGetObjective(SCOPFLOW scopflow,PetscReal *obj)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*scopflow->solverops.getobjective)(scopflow,obj);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/**
 * @brief Returns the SCOPFLOW solution for a given contingency
 * 
 * @param[in] scopflow the scopflow object
 * @param[in] contnum Contingency number (0 for base/no-contingency)
 * @param[in] X the scopflow solution for the given contingency
 * 
 * @note Should be called after the optimization finishes
*/
PetscErrorCode SCOPFLOWGetSolution(SCOPFLOW scopflow,PetscInt contnum,Vec *X)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*scopflow->solverops.getsolution)(scopflow,contnum,X);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/**
 * @brief Returns the SCOPFLOW constraints for a given contingency
 * 
 * @param[in] scopflow the scopflow object
 * @param[in] contnum contingency number (0 for base/no-contingency)
 * @param[in] G the scopflow constraints
 * 
 * @note Should be called after the optimization finishes.
 * Equality constraints first followed by inequality constraints
*/
PetscErrorCode SCOPFLOWGetConstraints(SCOPFLOW scopflow,PetscInt contnum,Vec *G)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*scopflow->solverops.getconstraints)(scopflow,contnum,G);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/**
 * @brief Returns the SCOPFLOW constraint multipliers for a given contingency
 * 
 * @param[in] scopflow the scopflow object
 * @param[in] contnum  Contingency number (0 for base/no-contingency)
 * @param[in] Lambda pointer to vector object
 * 
 * @note Should be called after the optimization finishes. 
 * Equality constraint multipliers first followed by inequality constraint multipliers
*/
PetscErrorCode SCOPFLOWGetConstraintMultipliers(SCOPFLOW scopflow,PetscInt contnum,Vec *Lambda)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*scopflow->solverops.getconstraintmultipliers)(scopflow,contnum,Lambda);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/**
 * @brief Check SCOPFLOW convergence
 * 
 * @param[in] scopflow the scopflow object
 * @param[in] status PETSC_TRUE if converged, PETSC_FALSE otherwise
 * 
 * @note Should be called after the optimization finishes
*/
PetscErrorCode SCOPFLOWGetConvergenceStatus(SCOPFLOW scopflow,PetscBool *status)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*scopflow->solverops.getconvergencestatus)(scopflow,status);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/**
 * @brief Gets the number of variables and constraints
 *
 * @param[in] scopflow the scopflow object
 * @param[in] nx number of variables
 * @param[in] ncon number of constraints
 *
 * @note Must be called after SCOPFLOWSetUp
*/
PetscErrorCode SCOPFLOWGetNumVariablesandConstraints(SCOPFLOW scopflow,PetscInt *Nx,PetscInt *Ncon)
{
  PetscFunctionBegin;
  *Nx = scopflow->nx;
  *Ncon = scopflow->ncon;
  PetscFunctionReturn(0);
}

/**
 * @brief Get the number of iterations for a given solver
 * 
 * @param[in] scopflow application object 
 * @param [in] iter pointer to solver iteration variable
*/
PetscErrorCode SCOPFLOWGetNumIterations(SCOPFLOW scopflow,PetscInt *iter)
{
  PetscErrorCode ierr;
  *iter = scopflow->numiter;
  PetscFunctionReturn(0);
}

/**
 * @brief Set the time step for SCOPFLOW
 * 
 * @param[in] scopflow application object 
 * @param[in] dT solver timestep 
*/
PetscErrorCode SCOPFLOWSetTimeStep(SCOPFLOW scopflow,PetscReal dT)
{
  PetscFunctionBegin;
  scopflow->dT = dT;
  PetscFunctionReturn(0);
}

/**
 * @brief Set the problem duration for SCOPFLOW
 * 
 * @param[in] scopflow application object 
 * @param[in] duration application duration
*/
PetscErrorCode SCOPFLOWSetDuration(SCOPFLOW scopflow,PetscReal duration)
{
  PetscFunctionBegin;
  scopflow->duration = duration;
  PetscFunctionReturn(0);
}

/**
 * @brief Set the solver tolerance for SCOPFLOW
 * 
 * @param[in] scopflow application object 
 * @param[in] tol solver tolerance
*/
PetscErrorCode SCOPFLOWSetTolerance(SCOPFLOW scopflow,PetscReal tol)
{
  PetscFunctionBegin;
  scopflow->tolerance = tol;
  PetscFunctionReturn(0);
}

/**
 * @brief Get the solver tolerance for SCOPFLOW
 * 
 * @param[in] scopflow application object 
 * @param[in] tol pointer to solver tolerance variable
*/
PetscErrorCode SCOPFLOWGetTolerance(SCOPFLOW scopflow,PetscReal *tol)
{
  PetscFunctionBegin;
  *tol = scopflow->tolerance;
  PetscFunctionReturn(0);
}

/**
 * @brief Set SCOPFLOW initialization type
 *
 * @param[in] scopflow application object
 * @param[in] type initialization type for underlying OPFLOW structs 
 */
PetscErrorCode SCOPFLOWSetInitilizationType(SCOPFLOW scopflow, OPFLOWInitializationType type)
{
  PetscFunctionBegin;
  scopflow->type = type;
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSetScenario(SCOPFLOW scopflow,Scenario *scen)
{
  PetscFunctionBegin;
  scopflow->scen = scen;
  PetscFunctionReturn(0);
}

/**
 * @brief Set SCOPFLOW number of contingencies
 *
 * @param[in] scopflow application object
 * @param[in] num_cont number of contingencies
 */
PetscErrorCode SCOPFLOWSetNumContingencies(SCOPFLOW scopflow, PetscInt num_cont)
{
  PetscFunctionBegin;
  scopflow->Nc = num_cont;
  PetscFunctionReturn(0);
}

/**
 * @brief Set SCOPFLOW genbusvoltage type
 *
 * @param[in] scopflow application object
 * @param[in] type genbusvoltage type for underlying OPFLOW structs
 */
PetscErrorCode SCOPFLOWSetGenBusVoltageType(SCOPFLOW scopflow, OPFLOWGenBusVoltageType type)
{
  PetscFunctionBegin;
  scopflow->genbusvoltagetype = type;
  PetscFunctionReturn(0);
}

/**
 * @brief Set SCOPFLOW mode
 *
 * @param[in] scopflow application object
 * @param[in] mode for scopflow. 0 = preventive, 1 = corrective 
 */
PetscErrorCode SCOPFLOWSetMode(SCOPFLOW scopflow, PetscInt mode)
{
  PetscFunctionBegin;
  scopflow->mode = mode;
  PetscFunctionReturn(0);
}

/*
 * @brief Set the contingency data file for SCOPFLOW
 * 
 * @param[in] scopflow application object 
 * @param[in] contingency file format
 * @param[in] name of the contingency list file
*/
PetscErrorCode SCOPFLOWSetContingencyData(SCOPFLOW scopflow,ContingencyFileInputFormat ctgcfileformat,const char ctgcfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscMemcpy(scopflow->ctgcfile,ctgcfile,PETSC_MAX_PATH_LEN*sizeof(char));CHKERRQ(ierr);
  scopflow->ctgcfileformat = ctgcfileformat;
  scopflow->ctgcfileset = PETSC_TRUE;
  PetscFunctionReturn(0);
}
