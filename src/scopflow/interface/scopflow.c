#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>

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

  scopflow->Nconeq   = scopflow->nconeq  = 0;
  scopflow->Nconineq = scopflow->nconineq = 0;
  scopflow->Ncon     = scopflow->ncon     = 0;
  scopflow->Nx       = scopflow->nx       = 0;
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
  scopflow->first_stage_gen_cost_only = PETSC_TRUE;
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
  PetscInt       i;

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
  if((*scopflow)->Nconineq) {
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

  /* Destroy OPFLOW objects */
  for(i=0; i < (*scopflow)->Ns; i++) {
    ierr = OPFLOWDestroy(&(*scopflow)->opflows[i]);CHKERRQ(ierr);
  }

  ierr = PetscFree((*scopflow)->nconineqcoup);CHKERRQ(ierr);
  ierr = PetscFree((*scopflow)->ctgclist.cont);CHKERRQ(ierr);
  ierr = PetscFree((*scopflow)->opflows);CHKERRQ(ierr);
  ierr = PetscFree(*scopflow);CHKERRQ(ierr);
  //  *scopflow = 0;
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
  SCOPFLOWReadContingencyData - Reads the contingency list data file

  Input Parameters
+ scopflow - the scopflow object
- ctgcfile - the contingency file name

*/
PetscErrorCode SCOPFLOWReadContingencyData(SCOPFLOW scopflow,const char ctgcfile[])
{
  PetscErrorCode ierr;
  FILE           *fp;
  ContingencyList *ctgclist=&scopflow->ctgclist;
  Contingency    *cont;
  Outage         *outage;
  char           line[MAXLINE];
  char           *out;
  PetscInt       bus,fbus,tbus,type,num;
  char           equipid[3];
  PetscInt       status;
  PetscScalar    prob;

  PetscFunctionBegin;

  fp = fopen(ctgcfile,"r");
  if (fp == NULL) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open file %s",ctgcfile);CHKERRQ(ierr);
  }

  ctgclist->Ncont = -1;

  while((out = fgets(line,MAXLINE,fp)) != NULL) {
    if(strcmp(line,"\r\n") == 0 || strcmp(line,"\n") == 0) {
      continue; /* Skip blank lines */
    }
    sscanf(line,"%d,%d,%d,%d,%d,'%[^\t\']',%d,%lf",&num,&type,&bus,&fbus,&tbus,equipid,&status,&prob);

    if(num == scopflow->Ns) break;

    if(num == MAX_CONTINGENCIES) {
      SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,"Exceeding max. allowed contingencies = %d\n",num,MAX_CONTINGENCIES);
    }
    cont   = &ctgclist->cont[num];
    outage = &cont->outagelist[cont->noutages];
    outage->num  = num;
    outage->type = (OutageType)type;
    outage->bus  = bus;
    outage->fbus = fbus;
    outage->tbus = tbus;
    ierr = PetscMemcpy(outage->id,equipid,3*sizeof(char));CHKERRQ(ierr);
    outage->status = status;
    outage->prob   = prob;
    cont->noutages++;


    if(num > ctgclist->Ncont) ctgclist->Ncont = num;
  }
  fclose(fp);

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
  PetscBool      solverset;
  char           formulationname[32]="POWER_BALANCE_CARTESIAN";
  char           solvername[32]="IPOPT";
  PetscInt       i,j;
  PS             ps;
  PetscBool      flg;

  PetscFunctionBegin;

  ierr = PetscOptionsBegin(scopflow->comm->type,NULL,"SCOPFLOW options",NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-scopflow_formulation","SCOPFLOW formulation type","",formulationname,formulationname,32,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-scopflow_solver","SCOPFLOW solver type","",solvername,solvername,32,&solverset);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-scopflow_iscoupling","Include coupling between first stage and second stage","",scopflow->iscoupling,&scopflow->iscoupling,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-scopflow_first_stage_gen_cost_only","Include objective cost for first stage only","",scopflow->first_stage_gen_cost_only,&scopflow->first_stage_gen_cost_only,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-scopflow_Ns","Number of second stage scenarios","",scopflow->Ns,&scopflow->Ns,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-scopflow_replicate_basecase","Only for debugging: Replicate first stage for all second stage scenarios","",scopflow->replicate_basecase,&scopflow->replicate_basecase,NULL);CHKERRQ(ierr);
  PetscOptionsEnd();

  if(scopflow->ctgcfileset && !scopflow->replicate_basecase) {
    if(scopflow->Ns < 0) scopflow->Ns = MAX_CONTINGENCIES;
    else scopflow->Ns += 1; 

    ierr = PetscCalloc1(scopflow->Ns,&scopflow->ctgclist.cont);CHKERRQ(ierr);
    for(i=0; i < scopflow->Ns; i++) scopflow->ctgclist.cont->noutages = 0;
    ierr = SCOPFLOWReadContingencyData(scopflow,scopflow->ctgcfile);CHKERRQ(ierr);
    scopflow->Ns = scopflow->ctgclist.Ncont+1;
  } else {
    if(scopflow->Ns == -1) scopflow->Ns = 1;
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD,"SCOPFLOW running with %d scenarios (base case + %d scenarios)\n",scopflow->Ns,scopflow->Ns-1);CHKERRQ(ierr);

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

  /* Create OPFLOW objects */
  ierr = PetscCalloc1(scopflow->Ns,&scopflow->opflows);CHKERRQ(ierr);
  for(i=0; i < scopflow->Ns; i++) {
    ierr = OPFLOWCreate(PETSC_COMM_SELF,&scopflow->opflows[i]);CHKERRQ(ierr);
    ierr = OPFLOWSetFormulation(scopflow->opflows[i],formulationname);CHKERRQ(ierr);
    ierr = PetscStrcmp(solvername,SCOPFLOWSOLVER_PIPS,&flg);CHKERRQ(ierr);
    if(flg) {
      ierr = OPFLOWSetSolver(scopflow->opflows[i],OPFLOWSOLVER_IPOPT);CHKERRQ(ierr);
    } else {
      ierr = OPFLOWSetSolver(scopflow->opflows[i],solvername);CHKERRQ(ierr);
    }
    ierr = OPFLOWReadMatPowerData(scopflow->opflows[i],scopflow->netfile);CHKERRQ(ierr);
    /* Set up the PS object for opflow */
    ps = scopflow->opflows[i]->ps;
    ierr = PSSetUp(ps);CHKERRQ(ierr);
    /* Set contingencies */
    if(scopflow->ctgcfileset && !scopflow->replicate_basecase) {
      Contingency ctgc=scopflow->ctgclist.cont[i];
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
    ierr = OPFLOWSetUp(scopflow->opflows[i]);CHKERRQ(ierr);
    if(i > 0 && scopflow->first_stage_gen_cost_only) scopflow->opflows[i]->obj_gencost = PETSC_FALSE; /* No gen. cost minimization for second stage */
  }
  
  ierr = PetscCalloc1(scopflow->Ns,&scopflow->nconineqcoup);CHKERRQ(ierr);
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

  PetscFunctionBegin;

  if(!scopflow->setupcalled) {
    ierr = SCOPFLOWSetUp(scopflow);
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
  SCOPFLOWGetBaseCaseSolution - Returns the SCOPFLOW base case (x0) solution

  Input Parameters:
+ SCOPFLOW - the SCOPFLOW object
- X        - the scopflow base case solution

  Notes: Should be called after the optimization finishes
*/
PetscErrorCode SCOPFLOWGetBaseCaseSolution(SCOPFLOW scopflow,Vec *X)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*scopflow->solverops.getbasecasesolution)(scopflow,X);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SCOPFLOWGetConstraints - Returns the SCOPFLOW constraints

  Input Parameters:
+ SCOPFLOW - the SCOPFLOW object
- G    - the scopflow constraints

  Notes: Should be called after the optimization finishes.
         Equality constraints first followed by inequality constraints
*/
PetscErrorCode SCOPFLOWGetConstraints(SCOPFLOW scopflow,Vec *G)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*scopflow->solverops.getconstraints)(scopflow,G);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SCOPFLOWGetConstraintMultipliers - Returns the SCOPFLOW constraint multipliers

  Input Parameters:
+ SCOPFLOW - the SCOPFLOW object
- G    - the scopflow constraint lagrange multipliers

  Notes: Should be called after the optimization finishes.
    Equality constraint multipliers first followed by inequality constraint multipliers
*/
PetscErrorCode SCOPFLOWGetConstraintMultipliers(SCOPFLOW scopflow,Vec *Lambda)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*scopflow->solverops.getconstraintmultipliers)(scopflow,Lambda);CHKERRQ(ierr);
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


