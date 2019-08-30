#include <private/psimpl.h>
#include <private/opflowimpl.h>
#include <private/scopflowipoptimpl.h>
#include <IpStdCInterface.h>

/* Function Declarations */
Bool eval_scopflow_f(PetscInt n, PetscScalar* x, Bool new_x,
            PetscScalar* obj_value, UserDataPtr user_data);

Bool eval_scopflow_grad_f(PetscInt n, PetscScalar* x, Bool new_x,
                 PetscScalar* grad_f, UserDataPtr user_data);

Bool eval_scopflow_g(PetscInt n, PetscScalar* x, Bool new_x,
            PetscInt m, PetscScalar* g, UserDataPtr user_data);

Bool eval_scopflow_jac_g(PetscInt n, PetscScalar *x, Bool new_x,
                PetscInt m, PetscInt nele_jac,
                PetscInt *iRow, PetscInt *jCol, PetscScalar *values,
                UserDataPtr user_data);

Bool eval_scopflow_h(PetscInt n, PetscScalar *x, Bool new_x, PetscScalar obj_factor,
            PetscInt m, PetscScalar *lambda, Bool new_lambda,
            PetscInt nele_hess, PetscInt *iRow, PetscInt *jCol,
            PetscScalar *values, UserDataPtr user_data);


static int CCMatrixToMatrixMarketValuesOnly(struct CCMatrix *ccmatrix,PetscInt nz,PetscScalar *values)
{
  PetscErrorCode ierr;

  ierr = PetscMemcpy(values,ccmatrix->values,nz*sizeof(PetscScalar));

  return 0;
}

static int CCMatrixToMatrixMarketLocationsOnly(struct CCMatrix *ccmatrix,PetscInt ncol,PetscInt *iRow,PetscInt *jCol,PetscInt roffset,PetscInt coffset,PetscInt nval)
{
  PetscInt *rowidx;
  PetscInt *colptr;
  PetscInt j,k,ctr=0;
  
  rowidx = ccmatrix->rowidx;
  colptr = ccmatrix->colptr;

  /* Copy from compressed column to (row,col,val) format */
  for(j=0; j < ncol; j++) {
    for(k=colptr[j]; k < colptr[j+1]; k++) {
      iRow[ctr] = rowidx[k] + roffset;
      jCol[ctr] = j + coffset;
      ctr++;
    }
  }
  if(ctr != nval) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,"Incorrect number of entries ctr = %d given = %d\n",ctr,nval);

  return 0;
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
  ContingencyList ctgclist=scopflow->ctgclist;
  Contingency    *cont;
  Outage         *outage;
  char           line[MAXLINE];
  char           *out;
  PetscInt       bus,fbus,tbus,type,num;
  char           id[10];
  PetscInt       status;
  PetscScalar    prob;

  PetscFunctionBegin;

  fp = fopen(ctgcfile,"r");
  if (fp == NULL) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open file %s",ctgcfile);CHKERRQ(ierr);
  }

  while((out = fgets(line,MAXLINE,fp)) != NULL) {
    if(strcmp(line,"\r\n") == 0 || strcmp(line,"\n") == 0) {
      continue; /* Skip blank lines */
    }
    sscanf(line,"%d,%d,%d,%d,%d,'%[^\t\']',%d,%lf",&num,&type,&bus,&fbus,&tbus,id,&status,&prob);

    cont = ctgclist.cont+num;
    outage = &cont->outagelist[cont->noutages];
    outage->num  = num;
    outage->type = type;
    outage->bus  = bus;
    outage->fbus = fbus;
    outage->tbus = tbus;
    ierr = PetscMemcpy(outage->id,id,2*sizeof(char));CHKERRQ(ierr);
    outage->status = status;
    outage->prob   = prob;
  }
  fclose(fp);

  PetscFunctionReturn(0);
}

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
  SCOPFLOW       scopflow;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&scopflow);CHKERRQ(ierr);

  ierr = COMMCreate(mpicomm,&scopflow->comm);CHKERRQ(ierr);

  if(scopflow->comm->size > 1) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"SCOPF with IPOPT is not supported in parallel");

  scopflow->Ns              = -1;
  scopflow->Nx              = -1;
  scopflow->Nc              = -1;
  scopflow->Jcoup           = NULL;
  scopflow->iscoupling      = PETSC_FALSE;
  scopflow->first_stage_gen_cost_only = PETSC_TRUE;
  scopflow->ignore_line_flow_constraints = PETSC_TRUE;

  scopflow->nlp_ipopt       = NULL;
  scopflow->ctgcfileset     = PETSC_FALSE;
  scopflow->setupcalled = PETSC_FALSE;
  
  *scopflowout = scopflow;
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

  FreeIpoptProblem((*scopflow)->nlp_ipopt);CHKERRQ(ierr);
  ierr = PetscFree(*scopflow);CHKERRQ(ierr);
  *scopflow = 0;

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


extern int str_init_x0(double*,CallBackDataPtr);

/*
  SCOPFLOWSetInitialGuess - Sets the initial guess for the optimization

  Input Parameters:
. scopflow - the SCOPFLOW object

  Output Parameters:
+ X     - initial guess

  Notes:
   Sets X[i] = (Xl[i] + Xu[i])/2
*/
PetscErrorCode SCOPFLOWSetInitialGuess(SCOPFLOW scopflow, Vec X)
{
  PetscErrorCode ierr;
  struct CallBackData cbd;
  PetscInt i;
  PetscScalar *x,*xi;

  PetscFunctionBegin;
  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  
  for(i=0; i < scopflow->Ns; i++) {
    cbd.row_node_id = i;
    cbd.col_node_id = i;
    cbd.prob = (void*)scopflow;
    
    xi = x + scopflow->xstart[i];

    str_init_x0(xi,&cbd);
  }
  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
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
  PetscScalar    *xl,*xu,*gl,*gu;
  char           hessiantype[100];
  PetscBool      flg;

  PetscFunctionBegin;
  if(!scopflow->setupcalled) {
    ierr = SCOPFLOWSetUp(scopflow);CHKERRQ(ierr);
  }

  if(scopflow->nlp_ipopt) {
    FreeIpoptProblem(scopflow->nlp_ipopt);
  }

  ierr = VecGetArray(scopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Gu,&gu);CHKERRQ(ierr);

  /* Create IPOPT problem */
  scopflow->nlp_ipopt = CreateIpoptProblem(scopflow->Nx,xl,xu,scopflow->Nc,gl,gu,scopflow->nnz_jac_g,scopflow->nnz_hes,0,&eval_scopflow_f,
				   &eval_scopflow_g,&eval_scopflow_grad_f,
				   &eval_scopflow_jac_g,&eval_scopflow_h);

  ierr = VecRestoreArray(scopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Gu,&gu);CHKERRQ(ierr);

  /* Options for IPOPT. This need to go through PetscOptionsBegin later */
  
  AddIpoptNumOption(scopflow->nlp_ipopt, (char*)"tol", 1e-4);
  AddIpoptNumOption(scopflow->nlp_ipopt, (char*)"acceptable_tol", 1e-4);
  AddIpoptNumOption(scopflow->nlp_ipopt, (char*)"mu_init", 0.01);
  AddIpoptNumOption(scopflow->nlp_ipopt,(char*)"obj_scaling_factor",100);
  
  AddIpoptNumOption(scopflow->nlp_ipopt,(char*) "bound_frac",0.01);
  AddIpoptNumOption(scopflow->nlp_ipopt, (char*)"bound_push",0.01);
  AddIpoptNumOption(scopflow->nlp_ipopt, (char*)"dual_inf_tol", 1e-2);
  AddIpoptNumOption(scopflow->nlp_ipopt, (char*)"compl_inf_tol", 1e-2);
  AddIpoptNumOption(scopflow->nlp_ipopt, (char*)"constr_viol_tol", 5e-6);

  //  AddIpoptStrOption(scopflow->nlp_ipopt,(char*)"fixed_variable_treatment",(char*)"relax_bounds");
  // AddIpoptStrOption(scopflow->nlp_ipopt,(char*)"nlp_scaling_method",(char*)"none");
  AddIpoptIntOption(scopflow->nlp_ipopt,(char*)"max_iter",500);
  //AddIpoptStrOption(scopflow->nlp_ipopt, (char*)"mu_strategy", (char*)"adaptive");
  AddIpoptStrOption(scopflow->nlp_ipopt, (char*)"print_user_options", (char*)"yes");
  AddIpoptStrOption(scopflow->nlp_ipopt, (char*)"output_file", (char*)"ipopt.out");
  
  //  AddIpoptStrOption(scopflow->nlp_ipopt,"warm_start_init_point","yes");
  ierr = PetscOptionsGetString(NULL,NULL,"-scopflow_hessian_type",hessiantype,sizeof(hessiantype),&flg);CHKERRQ(ierr);
  if(flg) {
    ierr = PetscStrcmp(hessiantype,"L-BFGS",&flg);CHKERRQ(ierr);
    if(flg) {
      AddIpoptStrOption(scopflow->nlp_ipopt, (char*)"hessian_approximation", (char*)"limited-memory");
    } else {
      ierr = PetscStrcmp(hessiantype,"EXACT",&flg);CHKERRQ(ierr);
      if(!flg) {
	SETERRQ(PETSC_COMM_SELF,0,"Incorrect input to SCOPFLOW Hessian type option -scopflow_hessian_type\n\tAvailable types are L-BFGS and EXACT");
      }
    }
  }
  //AddIpoptStrOption(scopflow->nlp_ipopt, (char*)"derivative_test", (char*)"second-order");
  AddIpoptStrOption(scopflow->nlp_ipopt,(char*)"linear_solver",(char*)"mumps");
  // AddIpoptNumOption(scopflow->nlp_ipopt,(char*)"bound_relax_factor",1e-4);
  
  /* Set Initial Guess */
  ierr = SCOPFLOWSetInitialGuess(scopflow,scopflow->X);CHKERRQ(ierr);

  PetscScalar *x;
  ierr = VecGetArray(scopflow->X,&x);CHKERRQ(ierr);

  /* Solve */
  scopflow->solve_status = IpoptSolve(scopflow->nlp_ipopt,x,NULL,&scopflow->obj,NULL,NULL,NULL,scopflow);

  ierr = VecRestoreArray(scopflow->X,&x);CHKERRQ(ierr);

  //  ierr = VecView(scopflow->X,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

extern int str_prob_info(int* n, double* col_lb, double* col_ub, int* m,
			 double* row_lb, double* row_ub, CallBackDataPtr cbd);
extern int str_eval_jac_g(double* x0, double* x1, int* e_nz, double* e_elts,
		int* e_rowidx, int* e_colptr, int* i_nz, double* i_elts, int* i_rowidx,
			  int* i_colptr, CallBackDataPtr cbd);

extern int str_eval_h(double* x0, double* x1, double* lambda, int* nz, double* elts,
		      int* rowidx, int* colptr, CallBackDataPtr cbd);


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
  PetscInt       i;
  struct CallBackData cbd;
  PetscInt       n,m;
  PetscScalar     *xl,*xu,*gl,*gu;
  PetscScalar     *x,*x0,*x1,*lambda,*lambda1;

  PetscFunctionBegin;
  if(scopflow->Ns == -1) {
    SETERRQ(scopflow->comm->type,0,"Must call SCOPFLOWSetNumScenarios or SCOPFLOWSetScenariosFile before calling this function\n");
  }

  /* Options */
  ierr = PetscOptionsGetBool(NULL,NULL,"-scopflow_iscoupling",&scopflow->iscoupling,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-scopflow_first_stage_gen_cost_only",&scopflow->first_stage_gen_cost_only,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-scopflow_ignore_line_flow_constraints",&scopflow->ignore_line_flow_constraints,NULL);CHKERRQ(ierr);

  if(scopflow->ctgcfileset) {
    scopflow->ctgclist.Ncont = scopflow->Ns;
    ierr = PetscMalloc1(scopflow->Ns,&scopflow->ctgclist.cont);CHKERRQ(ierr);
    for(i=0; i < scopflow->Ns; i++) scopflow->ctgclist.cont->noutages = 0;
    ierr = SCOPFLOWReadContingencyData(scopflow,scopflow->ctgcfile);CHKERRQ(ierr);
  }

  ierr = PetscMalloc1(scopflow->Ns,&scopflow->opflows);CHKERRQ(ierr);
  /* Starting locations for x and g for each scenario */
  ierr = PetscCalloc1(scopflow->Ns,&scopflow->xstart);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->Ns,&scopflow->gstart);CHKERRQ(ierr);
  /* Number of variables and constraints for each scenario */
  ierr = PetscCalloc1(scopflow->Ns,&scopflow->Nxi);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->Ns,&scopflow->Nci);CHKERRQ(ierr);

  /* Equality and Inequality Jacobian blocks for each scenario */
  ierr = PetscCalloc1(scopflow->Ns,&scopflow->e_nz_jac_self);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->Ns,&scopflow->i_nz_jac_self);CHKERRQ(ierr);

  ierr = PetscCalloc1(scopflow->Ns,&scopflow->e_jac_self);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->Ns,&scopflow->i_jac_self);CHKERRQ(ierr);

  /* Equality and Inequality Jacobian blocks for coupling */
  ierr = PetscCalloc1(scopflow->Ns,&scopflow->e_nz_jac_coupled);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->Ns,&scopflow->i_nz_jac_coupled);CHKERRQ(ierr);

  ierr = PetscCalloc1(scopflow->Ns,&scopflow->e_jac_coupled);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->Ns,&scopflow->i_jac_coupled);CHKERRQ(ierr);

  /* Hessian block for each scenario */
  ierr = PetscCalloc1(scopflow->Ns,&scopflow->nz_hess_self);CHKERRQ(ierr);
  ierr = PetscCalloc1(scopflow->Ns,&scopflow->hess_self);CHKERRQ(ierr);


  scopflow->Nx = 0;
  scopflow->Nc = 0;
  scopflow->xstart[0] = 0;
  scopflow->gstart[0] = 0;

  /* Get the number of variables and constraints and
     set up vectors,matrices, etc. for scenarios
  */
  for(i=0; i < scopflow->Ns; i++) {
    cbd.row_node_id = i;
    cbd.col_node_id = i;
    cbd.prob = (void*)scopflow;

    str_prob_info(&n, NULL, NULL, &m,NULL,NULL,&cbd);

    scopflow->Nx += n;
    scopflow->Nc += m;
    scopflow->Nxi[i] = n;
    scopflow->Nci[i] = m;
    if(i < scopflow->Ns-1) {
      scopflow->xstart[i+1] = scopflow->Nx;
      scopflow->gstart[i+1] = scopflow->Nc;
    }
  }

  /* Create vector X */
  ierr = VecCreate(scopflow->comm->type,&scopflow->X);CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->X,scopflow->Nx,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(scopflow->X);CHKERRQ(ierr);

  ierr = VecDuplicate(scopflow->X,&scopflow->Xl);CHKERRQ(ierr);
  ierr = VecDuplicate(scopflow->X,&scopflow->Xu);CHKERRQ(ierr);

  /* vector for constraints */
  ierr = VecCreate(scopflow->comm->type,&scopflow->G);CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->G,scopflow->Nc,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(scopflow->G);CHKERRQ(ierr);

  /* Constraint bounds vectors  */
  ierr = VecDuplicate(scopflow->G,&scopflow->Gl);CHKERRQ(ierr);
  ierr = VecDuplicate(scopflow->G,&scopflow->Gu);CHKERRQ(ierr);
  
  /* Lagrangian multipliers */
  ierr = VecDuplicate(scopflow->G,&scopflow->Lambda);CHKERRQ(ierr);

  /* Set Bounds on variables and constraints */
  ierr = VecGetArray(scopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Gu,&gu);CHKERRQ(ierr);

  /* Set the variable and constraint bounds */
  for(i=0; i < scopflow->Ns; i++) {
    PetscScalar *xli,*xui,*gli,*gui;
    cbd.row_node_id = i;
    cbd.col_node_id = i;
    cbd.prob = (void*)scopflow;

    xli = xl + scopflow->xstart[i];
    xui = xu + scopflow->xstart[i];
    gli = gl + scopflow->gstart[i];
    gui = gu + scopflow->gstart[i];

    str_prob_info(&n,xli , xui, &m,gli,gui,&cbd);
  }

  ierr = VecRestoreArray(scopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Gu,&gu);CHKERRQ(ierr);

  /* Jacobian and Hessian non-zeros */
  scopflow->nnz_jac_g = 0;
  scopflow->nnz_hes = 0;

  ierr = VecSet(scopflow->X,1.0);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->X,&x);CHKERRQ(ierr);
  ierr = VecSet(scopflow->Lambda,1.0);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Lambda,&lambda);CHKERRQ(ierr);

  for(i=0; i < scopflow->Ns; i++) {
    cbd.row_node_id = i;
    cbd.col_node_id = i;
    cbd.prob = (void*)scopflow;

    x0 = x + scopflow->xstart[i];
    x1 = x + scopflow->xstart[i];

    /* Self part - Get Nonzeros */
    str_eval_jac_g(x0, x1, scopflow->e_nz_jac_self+i, NULL,NULL,NULL,
		   scopflow->i_nz_jac_self+i,NULL,NULL,NULL,&cbd);

    scopflow->nnz_jac_g += scopflow->e_nz_jac_self[i] + scopflow->i_nz_jac_self[i];

    ierr = PetscMalloc1(scopflow->opflows[i]->Nvar+1,&scopflow->e_jac_self[i].colptr);CHKERRQ(ierr);
    ierr = PetscMalloc1(scopflow->e_nz_jac_self[i],&scopflow->e_jac_self[i].rowidx);CHKERRQ(ierr);
    ierr = PetscMalloc1(scopflow->e_nz_jac_self[i],&scopflow->e_jac_self[i].values);CHKERRQ(ierr);

    if(scopflow->i_nz_jac_self[i]) {
      ierr = PetscMalloc1(scopflow->opflows[i]->Nvar+1,&scopflow->i_jac_self[i].colptr);CHKERRQ(ierr);
      ierr = PetscMalloc1(scopflow->i_nz_jac_self[i],&scopflow->i_jac_self[i].rowidx);CHKERRQ(ierr);
      ierr = PetscMalloc1(scopflow->i_nz_jac_self[i],&scopflow->i_jac_self[i].values);CHKERRQ(ierr);
    }

    lambda1 = lambda + scopflow->gstart[i];
    str_eval_h(x0, x1, lambda1, scopflow->nz_hess_self+i, NULL,
	       NULL, NULL, &cbd);

    ierr = PetscMalloc1(scopflow->opflows[i]->Nvar+1,&scopflow->hess_self[i].colptr);CHKERRQ(ierr);
    ierr = PetscMalloc1(scopflow->nz_hess_self[i],&scopflow->hess_self[i].rowidx);CHKERRQ(ierr);
    ierr = PetscMalloc1(scopflow->nz_hess_self[i],&scopflow->hess_self[i].values);CHKERRQ(ierr);

    scopflow->nnz_hes += scopflow->nz_hess_self[i];

    cbd.col_node_id = 0;
    /* Coupling part */
    if(i > 0) {
      x0 = x + scopflow->xstart[0];
      str_eval_jac_g(x0, x1, scopflow->e_nz_jac_coupled+i, NULL,NULL,NULL,
		     scopflow->i_nz_jac_coupled+i,NULL,NULL,NULL,&cbd);

      if(scopflow->i_nz_jac_coupled[i]) {
	ierr = PetscMalloc1(scopflow->opflows[i]->Nvar+1,&scopflow->i_jac_coupled[i].colptr);CHKERRQ(ierr);
	ierr = PetscMalloc1(scopflow->i_nz_jac_coupled[i],&scopflow->i_jac_coupled[i].rowidx);CHKERRQ(ierr);
	ierr = PetscMalloc1(scopflow->i_nz_jac_coupled[i],&scopflow->i_jac_coupled[i].values);CHKERRQ(ierr);
      }

      scopflow->nnz_jac_g += scopflow->e_nz_jac_coupled[i] + scopflow->i_nz_jac_coupled[i];

    }
  }
  ierr = VecRestoreArray(scopflow->X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Lambda,&lambda);CHKERRQ(ierr);
    
  scopflow->setupcalled = PETSC_TRUE;
  PetscFunctionReturn(0);
}

extern int str_eval_f(double*,double*,double*,CallBackDataPtr);

Bool eval_scopflow_f(PetscInt n, PetscScalar* x, Bool new_x,
            PetscScalar* obj_value, UserDataPtr user_data)
{
  SCOPFLOW       scopflow=(SCOPFLOW)user_data;
  PetscInt       i; 
  struct CallBackData cbd;
  PetscScalar    *x0,*x1;
  PetscScalar    obj_self=0.0;

  x0 = x + scopflow->xstart[0];
  
  *obj_value = 0.0;
  for(i=0; i < scopflow->Ns; i++) {
    cbd.row_node_id = i;
    cbd.col_node_id = i;
    cbd.prob = (void*)scopflow;
    
    x1 = x + scopflow->xstart[i];

    str_eval_f(x0,x1,&obj_self,&cbd);

    *obj_value += obj_self;
  }
  return TRUE;
}

extern int str_eval_grad_f(double*,double*,double*,CallBackDataPtr);

Bool eval_scopflow_grad_f(PetscInt n, PetscScalar* x, Bool new_x,
                 PetscScalar* grad_f, UserDataPtr user_data)
{
  SCOPFLOW       scopflow=(SCOPFLOW)user_data;
  PetscInt       i; 
  struct CallBackData cbd;
  PetscScalar    *x0,*x1,*df1;

  x0 = x + scopflow->xstart[0];
  
  for(i=0; i < scopflow->Ns; i++) {
    cbd.row_node_id = i;
    cbd.col_node_id = i;
    cbd.prob = (void*)scopflow;
    
    x1 = x + scopflow->xstart[i];
    df1 = grad_f + scopflow->xstart[i];

    str_eval_grad_f(x0,x1,df1,&cbd);
  }
  return TRUE;
}

extern int str_eval_g(double* x0, double* x1, double* eq_g, double* inq_g,
	       CallBackDataPtr cbd);

Bool eval_scopflow_g(PetscInt n, PetscScalar* x, Bool new_x,
             PetscInt m, PetscScalar* c, UserDataPtr user_data)
{
  SCOPFLOW       scopflow=(SCOPFLOW)user_data;
  PetscInt       i; 
  struct CallBackData cbd;
  PetscScalar    *x0,*x1;
  PetscScalar    *eq_g_self,*inq_g_self,*g1;

  x0 = x + scopflow->xstart[0];
  
  for(i=0; i < scopflow->Ns; i++) {
    cbd.row_node_id = i;
    cbd.col_node_id = i;
    cbd.prob = (void*)scopflow;
    
    x1 = x + scopflow->xstart[i];
    g1 = c + scopflow->gstart[i];
    eq_g_self = g1;
    inq_g_self = g1 + scopflow->opflows[i]->Nconeq;

    str_eval_g(x0,x1,eq_g_self,inq_g_self,&cbd);

  }
  return TRUE;
}

Bool eval_scopflow_jac_g(PetscInt n, PetscScalar *x, Bool new_x,
                PetscInt m, PetscInt nele_jac,
                PetscInt *iRow, PetscInt *jCol, PetscScalar *values,
                UserDataPtr user_data)
{
  PetscErrorCode ierr;
  SCOPFLOW       scopflow=(SCOPFLOW)user_data;
  PetscInt       i; 
  struct CallBackData cbd;
  PetscScalar    *x0,*x1;
  PetscInt       *iRowstart=iRow,*jColstart=jCol;
  PetscInt       roffset,coffset;
  PetscScalar    *xdup;
  Vec            Xdup;

  if(values == NULL) {
    /* x is null when values == NULL, so we create a duplicate vector to pass x0 and x1.
       This is used only for the getting the locations only 
    */
    ierr = VecDuplicate(scopflow->X,&Xdup);CHKERRQ(ierr);
    //    ierr = VecSet(Xdup,1.0);CHKERRQ(ierr);
    ierr = SCOPFLOWSetInitialGuess(scopflow,Xdup);CHKERRQ(ierr);
    ierr = VecGetArray(Xdup,&xdup);CHKERRQ(ierr);

    x0 = xdup + scopflow->xstart[0];

    for(i=0; i < scopflow->Ns; i++) {
      cbd.row_node_id = i;
      cbd.col_node_id = i;
      cbd.prob = (void*)scopflow;

      x1 = xdup + scopflow->xstart[i];

      /* Self part - Get Locations */
      str_eval_jac_g(x0, x1, scopflow->e_nz_jac_self+i,scopflow->e_jac_self[i].values,scopflow->e_jac_self[i].rowidx,scopflow->e_jac_self[i].colptr,scopflow->i_nz_jac_self+i,scopflow->i_jac_self[i].values,scopflow->i_jac_self[i].rowidx,scopflow->i_jac_self[i].colptr,&cbd);
      
      roffset = scopflow->gstart[i];
      coffset = scopflow->xstart[i];

      /* Insert equality constraints self part */
      CCMatrixToMatrixMarketLocationsOnly(&scopflow->e_jac_self[i],scopflow->opflows[i]->Nvar,iRowstart,jColstart,roffset,coffset,scopflow->e_nz_jac_self[i]);

      /* Increase iRow,jCol pointers */
      iRowstart += scopflow->e_nz_jac_self[i];
      jColstart += scopflow->e_nz_jac_self[i];

      PetscInt Nconineq = scopflow->Nci[i] - scopflow->opflows[i]->Nconeq;

      if(Nconineq) {
	roffset = scopflow->gstart[i] + scopflow->opflows[i]->Nconeq;

	/* Insert inequality constraints self part */
	CCMatrixToMatrixMarketLocationsOnly(&scopflow->i_jac_self[i],scopflow->opflows[i]->Nvar,iRowstart,jColstart,roffset,coffset,scopflow->i_nz_jac_self[i]);

	/* Increase iRow, jCol pointers */
	iRowstart += scopflow->i_nz_jac_self[i];
	jColstart += scopflow->i_nz_jac_self[i];
      }
      /* Coupling part */
      if(i > 0) {
	if(scopflow->i_nz_jac_coupled[i]) {
	  cbd.col_node_id = 0;

	  str_eval_jac_g(x0, x1, scopflow->e_nz_jac_coupled+i,scopflow->e_jac_coupled[i].values,scopflow->e_jac_coupled[i].rowidx,scopflow->e_jac_coupled[i].colptr,scopflow->i_nz_jac_coupled+i,scopflow->i_jac_coupled[i].values,scopflow->i_jac_coupled[i].rowidx,scopflow->i_jac_coupled[i].colptr,&cbd);  

	  roffset = scopflow->gstart[i] + scopflow->opflows[i]->Nconeq;
	  coffset = scopflow->xstart[0];

	  /* Insert inequality constraints coupled part */
	  CCMatrixToMatrixMarketLocationsOnly(&scopflow->i_jac_coupled[i],scopflow->opflows[0]->Nvar,iRowstart,jColstart,roffset,coffset,scopflow->i_nz_jac_coupled[i]);

	  iRowstart += scopflow->i_nz_jac_coupled[i];
	  jColstart += scopflow->i_nz_jac_coupled[i];
	}
      }
    }
    ierr = VecRestoreArray(Xdup,&xdup);CHKERRQ(ierr);
    ierr = VecDestroy(&Xdup);CHKERRQ(ierr);
  } else {
    PetscScalar *valuesstart = values;

    x0 = x + scopflow->xstart[0];

    for(i=0; i < scopflow->Ns; i++) {
      cbd.row_node_id = i;
      cbd.col_node_id = i;
      cbd.prob = (void*)scopflow;

      x1 = x + scopflow->xstart[i];

      /* Self part - Get Locations */
      str_eval_jac_g(x0, x1, scopflow->e_nz_jac_self+i,scopflow->e_jac_self[i].values,scopflow->e_jac_self[i].rowidx,scopflow->e_jac_self[i].colptr,scopflow->i_nz_jac_self+i,scopflow->i_jac_self[i].values,scopflow->i_jac_self[i].rowidx,scopflow->i_jac_self[i].colptr,&cbd);
      
      /* Insert equality constraints self part */
      CCMatrixToMatrixMarketValuesOnly(&scopflow->e_jac_self[i],scopflow->e_nz_jac_self[i],valuesstart);

      /* Increase valuesstart pointers */
      valuesstart += scopflow->e_nz_jac_self[i];

      PetscInt Nconineq = scopflow->Nci[i] - scopflow->opflows[i]->Nconeq;

      if(Nconineq) {
	roffset = scopflow->gstart[i] + scopflow->opflows[i]->Nconeq;

	/* Insert inequality constraints self part */
	CCMatrixToMatrixMarketValuesOnly(&scopflow->i_jac_self[i],scopflow->i_nz_jac_self[i],valuesstart);

	/* Increase valuesstart pointers */
	valuesstart += scopflow->i_nz_jac_self[i];

      }
      /* Coupling part */
      if(i > 0) {
	if(scopflow->i_nz_jac_coupled[i]) {
	  cbd.col_node_id = 0;

	  str_eval_jac_g(x0, x1, scopflow->e_nz_jac_coupled+i,scopflow->e_jac_coupled[i].values,scopflow->e_jac_coupled[i].rowidx,scopflow->e_jac_coupled[i].colptr,scopflow->i_nz_jac_coupled+i,scopflow->i_jac_coupled[i].values,scopflow->i_jac_coupled[i].rowidx,scopflow->i_jac_coupled[i].colptr,&cbd);  

	  /* Insert coupled constraints self part */
	  CCMatrixToMatrixMarketValuesOnly(&scopflow->i_jac_coupled[i],scopflow->i_nz_jac_coupled[i],valuesstart);

	  /* Increase valuesstart pointers */
	  valuesstart += scopflow->i_nz_jac_coupled[i];

	}
      }
    }    
  }
  return TRUE;
}

Bool eval_scopflow_h(PetscInt n, PetscScalar *x, Bool new_x, PetscScalar obj_factor,
                PetscInt m, PetscScalar *lambda, Bool new_lambda, 
                PetscInt nele_hess, PetscInt *iRow, PetscInt *jCol, 
                PetscScalar *values, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  SCOPFLOW       scopflow=(SCOPFLOW)user_data;
  PetscInt       i; 
  struct CallBackData cbd;
  PetscScalar    *x0,*x1;
  PetscInt       *rowidx,*colptr;
  PetscInt       *iRowstart=iRow,*jColstart=jCol;
  PetscInt       roffset,coffset;
  PetscScalar    *xdup,*lambdadup,*lambda1;
  Vec            Xdup;
  PetscInt       j,k,ctr=0;

  scopflow->obj_factor = obj_factor;
  if(values == NULL) {
    /* x is null when values == NULL, so we create a duplicate vector to pass x0 and x1.
       This is used only for the getting the locations only 
    */
    ierr = VecDuplicate(scopflow->X,&Xdup);CHKERRQ(ierr);
    //    ierr = VecSet(Xdup,1.0);CHKERRQ(ierr);
    ierr = SCOPFLOWSetInitialGuess(scopflow,Xdup);CHKERRQ(ierr);
    ierr = VecGetArray(Xdup,&xdup);CHKERRQ(ierr);
    ierr = VecSet(scopflow->Lambda,1.0);CHKERRQ(ierr);
    ierr = VecGetArray(scopflow->Lambda,&lambdadup);CHKERRQ(ierr);

    x0 = xdup + scopflow->xstart[0];

    for(i=0; i < scopflow->Ns; i++) {
      cbd.row_node_id = i;
      cbd.col_node_id = i;
      cbd.prob = (void*)scopflow;

      x1 = xdup + scopflow->xstart[i];
      lambda1 = lambdadup + scopflow->gstart[i];

      /* Self part - Get Locations */
      str_eval_h(x0, x1, lambda1, scopflow->nz_hess_self+i,scopflow->hess_self[i].values,scopflow->hess_self[i].rowidx, scopflow->hess_self[i].colptr, &cbd);
      
      roffset = scopflow->xstart[i];
      coffset = scopflow->xstart[i];
  
      rowidx = scopflow->hess_self[i].rowidx;
      colptr = scopflow->hess_self[i].colptr;

      ctr = 0;
      /* Copy from compressed column to (row,col,val) format */
      for(j=0; j < scopflow->opflows[i]->Nvar; j++) {
	for(k=colptr[j]; k < colptr[j+1]; k++) {
	  jColstart[ctr] = rowidx[k] + roffset;
	  iRowstart[ctr] = j + coffset;
	  ctr++;
	}
      }

      /* Hessian self part */
      //      CCMatrixToMatrixMarketLocationsOnly(&scopflow->hess_self[i],scopflow->opflows[i]->Nvar,iRowstart,jColstart,roffset,coffset,scopflow->nz_hess_self[i]);

      /* Increase iRow,jCol pointers */
      iRowstart += scopflow->nz_hess_self[i];
      jColstart += scopflow->nz_hess_self[i];

    }
    ierr = VecRestoreArray(Xdup,&xdup);CHKERRQ(ierr);
    ierr = VecDestroy(&Xdup);CHKERRQ(ierr);
    ierr = VecRestoreArray(scopflow->Lambda,&lambdadup);CHKERRQ(ierr);
  } else {
    PetscScalar *valuesstart = values;

    x0 = x + scopflow->xstart[0];

    for(i=0; i < scopflow->Ns; i++) {
      cbd.row_node_id = i;
      cbd.col_node_id = i;
      cbd.prob = (void*)scopflow;

      x1 = x + scopflow->xstart[i];
      lambda1 = lambda + scopflow->gstart[i];

      /* Self part - Values */
      str_eval_h(x0, x1, lambda1, scopflow->nz_hess_self+i,scopflow->hess_self[i].values,scopflow->hess_self[i].rowidx, scopflow->hess_self[i].colptr, &cbd);

      /* Hessian self part */
      CCMatrixToMatrixMarketValuesOnly(&scopflow->hess_self[i],scopflow->nz_hess_self[i],valuesstart);

      /* Increase valuesstart pointers */
      valuesstart += scopflow->nz_hess_self[i];
    }    
  }
  return TRUE;
}

/*
  SCOPFLOWSetNumScenarios - Sets the total number of scenarios in the SCOPF problem

  Input Parameters:
+ scopflow - the security constrained optimal power flow application object
- Ns       - the number of scenarios

  Notes: The total number of scenarios set by SCOPFLOW is actually Ns+1,
  i.e., Ns scenarios + 1 base-case
*/
PetscErrorCode SCOPFLOWSetNumScenarios(SCOPFLOW scopflow,PetscInt Ns)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  scopflow->Ns = Ns;
  ierr = PetscOptionsGetInt(NULL,NULL,"-scopflow_Ns",&scopflow->Ns,NULL);CHKERRQ(ierr);
  scopflow->Ns += 1; /* One addition for the base case */
  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"SCOPFLOW running with %d scenarios (base case + %d scenarios)\n",scopflow->Ns,scopflow->Ns-1);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


  
