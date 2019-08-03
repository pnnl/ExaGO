#include <private/psimpl.h>
#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>
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

/*
  SCOPFLOWGetConstraintJacobianNonzeros - Gets the number of nonzeros in the lagrangian hessian matrix

  Input Parameters:
. scopflow - the security constrained optimal power flow application object

  Output Parameters:
. nnz - number of nonzeros in the lagrangian hessian

  Notes:
  SCOPFLOWSetUp() must be called before calling this routine.
*/
PetscErrorCode SCOPFLOWGetLagrangianHessianNonzeros(SCOPFLOW scopflow,PetscInt *nnz)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  *nnz = 0;

  PetscFunctionReturn(0);
}

/*
  SCOPFLOWSetConstraintJacobianLocations - Sets the row and column nonzero locations for the
              lagrangian hessian

  Input Parameters:
+ scopflow - the SCOPFLOW object

  Output Parameters:
+ row - the row locations
- col - the col locations

  Notes:
   The row and column locations should be such that H(row[i],col[i]) = val
*/
PetscErrorCode SCOPFLOWSetLagrangianHessianLocations(SCOPFLOW scopflow, PetscInt *row, PetscInt *col)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  PetscFunctionReturn(0);
}



/*
  SCOPFLOWSetHessianValues - Sets the nonzero values for the
              Lagrangian Hessian

  Input Parameters:
+ scopflow - the SCOPFLOW object
- X      - the current iterate


  Output Parameters:
. values - the nonzero values in the Lagrangian Hessian

  Notes:
   The values should be in the same row and col order as set in SCOPFLOWSetLagrangianHessianLocations
*/
PetscErrorCode SCOPFLOWSetLagrangianHessianValues(SCOPFLOW scopflow, PetscScalar obj_factor, Vec X, Vec Lambda, PetscScalar *values)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

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

  //  if(scopflow->comm->size > 1) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"SCOPF with IPOPT is not supported in parallel");

  scopflow->ns              = -1;
  scopflow->Ns              = -1;
  scopflow->nc              =  0;
  scopflow->Nc              =  0;
  scopflow->nx              =  0;
  scopflow->Nx              =  0;

  scopflow->nlp_ipopt       = NULL;

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
  OPFLOW         *opflows;
  PetscInt       i;

  PetscFunctionBegin;
  if(scopflow->ns == -1) {
    SETERRQ(scopflow->comm->type,0,"Must call SCOPFLOWSetNumScenarios or SCOPFLOWSetScenariosFile before calling this function\n");
  }
  ierr = PetscMalloc1(scopflow->ns,&opflows);CHKERRQ(ierr);
  scopflow->opflow = opflows;
  for(i=0; i < scopflow->ns; i++) {
    /* Create OPFLOW object for each scenario */
    ierr = OPFLOWCreate(PETSC_COMM_SELF,&opflows[i]);CHKERRQ(ierr);
    /* Have the OPFLOW object read MatPower data file */
    ierr = OPFLOWReadMatPowerData(opflows[i],netfile);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  SCOPFLOWCreateGlobalVector - Returns a global vector of the appropriate size
  and distribution conforming to the distribution of the PS object.

  Input Paramereters:
. scopflow - the security constrained optimal power flow application object

  Output Parameters:
. vec - the global vector

  Notes:
  SCOPFLOWSetUp() must be called before calling this routine.
*/
PetscErrorCode SCOPFLOWCreateGlobalVector(SCOPFLOW scopflow,Vec *vec)
{
  PetscErrorCode ierr;
  PetscInt       nx;
  OPFLOW         *opflows=scopflow->opflow;

  PetscFunctionBegin;
  if(!scopflow->setupcalled) SETERRQ(scopflow->comm->type,0,"SCOPFLOWSetUp() must be called before calling SCOPFLOWCreateGlobalVector");
  
  nx = opflows[0]->nvar*scopflow->ns;
  ierr = VecCreate(scopflow->comm->type,&scopflow->X);CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->X,nx,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(scopflow->X);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  SCOPFLOWCreateConstraintJacobian - Returns a distributed matrix of appropriate size that can
                                   be used as the Jacobian


  Input Paramereters:
. scopflow - the security constrained optimal power flow application object

  Output Parameters:
. mat - the jacobian of the constraints

  Notes:
  SCOPFLOWSetUp() must be called before calling this routine.
*/
PetscErrorCode SCOPFLOWCreateConstraintJacobian(SCOPFLOW scopflow,Mat *mat)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

/*
  SCOPFLOWGetConstraintJacobianNonzeros - Gets the number of nonzeros in the constraint jacobian matrix

  Input Paramereters:
. scopflow - the security constrained optimal power flow application object

  Output Parameters:
. nnz - number of nonzeros in the constraint Jacobian

  Notes:
  SCOPFLOWSetUp() must be called before calling this routine.
*/
PetscErrorCode SCOPFLOWGetConstraintJacobianNonzeros(SCOPFLOW scopflow,PetscInt *nnz)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  *nnz = 0;

  PetscFunctionReturn(0);
}

/*
  SCOPFLOWSetVariableandConstraintBounds - Sets the bounds on variables and constraints

  Input Parameters:
. scopflow - the SCOPFLOW object

  Output Parameters:
+ Xl     - vector of lower bound on variables
. Xu     - vector of upper bound on variables
. Gl     - vector of lower bound on constraints
. Gu     - vector of lower bound on constraints
*/
PetscErrorCode SCOPFLOWSetVariableandConstraintBounds(SCOPFLOW scopflow, Vec Xl, Vec Xu, Vec Gl, Vec Gu)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

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

  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWConstraintFunction(SCOPFLOW,Vec,Vec);

/*
  SCOPFLOWSolve - Solves the AC security constrained optimal power flow

  Input Parameters:
. scopflow - the security constrained optimal power flow application object
*/
PetscErrorCode SCOPFLOWSolve(SCOPFLOW scopflow)
{
  PetscErrorCode ierr;
  PetscScalar    *xl,*xu,*cl,*cu;
  char           hessiantype[100];
  PetscBool      flg;

  PetscFunctionBegin;
  if(!scopflow->setupcalled) {
    ierr = SCOPFLOWSetUp(scopflow);CHKERRQ(ierr);
  }

  if(scopflow->nlp_ipopt) {
    FreeIpoptProblem(scopflow->nlp_ipopt);
  }

  /* Set bounds on variables and constraints */
  ierr = SCOPFLOWSetVariableandConstraintBounds(scopflow,scopflow->Xl,scopflow->Xu,scopflow->Cl,scopflow->Cu);CHKERRQ(ierr);

  ierr = VecGetArray(scopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Cl,&cl);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Cu,&cu);CHKERRQ(ierr);

  /* Create IPOPT problem */
  scopflow->nlp_ipopt = CreateIpoptProblem(scopflow->Nx,xl,xu,scopflow->Nc,cl,cu,scopflow->nnz_jac_g,scopflow->nnz_hes,0,&eval_scopflow_f,
				   &eval_scopflow_g,&eval_scopflow_grad_f,
				   &eval_scopflow_jac_g,&eval_scopflow_h);

  ierr = VecRestoreArray(scopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Cl,&cl);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Cu,&cu);CHKERRQ(ierr);

  /* Options for IPOPT. This need to go through PetscOptionsBegin later */
  
  AddIpoptNumOption(scopflow->nlp_ipopt, (char*)"tol", 1e-6);
  AddIpoptNumOption(scopflow->nlp_ipopt, (char*)"acceptable_tol", 1e-6);
  AddIpoptNumOption(scopflow->nlp_ipopt, (char*)"mu_init", 0.1);
  AddIpoptNumOption(scopflow->nlp_ipopt,(char*)"obj_scaling_factor",100);
  
  AddIpoptNumOption(scopflow->nlp_ipopt,(char*) "bound_frac",1);
  AddIpoptNumOption(scopflow->nlp_ipopt, (char*)"bound_push",1);
  AddIpoptNumOption(scopflow->nlp_ipopt, (char*)"dual_inf_tol", 1e-2);
  AddIpoptNumOption(scopflow->nlp_ipopt, (char*)"compl_inf_tol", 1e-2);
  AddIpoptNumOption(scopflow->nlp_ipopt, (char*)"constr_viol_tol", 5e-6);

  //  AddIpoptStrOption(scopflow->nlp_ipopt,(char*)"fixed_variable_treatment",(char*)"relax_bounds");
  // AddIpoptStrOption(scopflow->nlp_ipopt,(char*)"nlp_scaling_method",(char*)"none");
  AddIpoptIntOption(scopflow->nlp_ipopt,(char*)"max_iter",500);
  //  AddIpoptStrOption(scopflow->nlp_ipopt, (char*)"mu_strategy", (char*)"adaptive");
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
  //  AddIpoptStrOption(scopflow->nlp_ipopt, (char*)"derivative_test", (char*)"first-order");
  AddIpoptStrOption(scopflow->nlp_ipopt,(char*)"linear_solver",(char*)"mumps");
  // AddIpoptNumOption(scopflow->nlp_ipopt,(char*)"bound_relax_factor",1e-4);
  
  /* Set Initial Guess */
  ierr = SCOPFLOWSetInitialGuess(scopflow,scopflow->X);CHKERRQ(ierr);

  Vec C;
  ierr = VecDuplicate(scopflow->C,&C);CHKERRQ(ierr);
  ierr = SCOPFLOWConstraintFunction(scopflow,scopflow->X,C);CHKERRQ(ierr);

  PetscScalar *x;
  ierr = VecGetArray(scopflow->X,&x);CHKERRQ(ierr);

  /* Solve */
  scopflow->solve_status = IpoptSolve(scopflow->nlp_ipopt,x,NULL,&scopflow->obj,NULL,NULL,NULL,scopflow);

  ierr = VecRestoreArray(scopflow->X,&x);CHKERRQ(ierr);

  if(!scopflow->comm->rank) {
    ierr = PetscPrintf(PETSC_COMM_SELF,"Objective function value = %lf\n",scopflow->obj);CHKERRQ(ierr);
  }
  ierr = VecView(scopflow->X,0);CHKERRQ(ierr);
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
  OPFLOW         *opflows=scopflow->opflow;
  PetscInt       ngen; /* Number of generators */
  PetscInt       i,cloc_size;

  PetscFunctionBegin;
  scopflow->nx = 0;
  /* Set up OPFLOW objects for scenarios */
  for(i=0; i < scopflow->ns; i++) {
    ierr = OPFLOWSetUp(opflows[i]);CHKERRQ(ierr);
  }

  if(scopflow->ns) {
    scopflow->nx = opflows[0]->n;
    ierr = PSGetNumGenerators(opflows[0]->ps,&ngen,NULL);CHKERRQ(ierr);
    scopflow->nd = ngen;
  }

  /* Create vectors X,Xloc, and X0 */
  ierr = VecCreate(PETSC_COMM_SELF,&scopflow->Xloc);CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->Xloc,scopflow->ns*scopflow->nx,scopflow->ns*scopflow->nx);CHKERRQ(ierr);
  ierr = VecSetFromOptions(scopflow->Xloc);CHKERRQ(ierr);

  ierr = VecDuplicate(scopflow->opflow[0]->X,&scopflow->X0);CHKERRQ(ierr);

  ierr = VecCreate(scopflow->comm->type,&scopflow->X);CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->X,scopflow->ns*scopflow->nx,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(scopflow->X);CHKERRQ(ierr);

  ierr = VecDuplicate(scopflow->X,&scopflow->Xl);CHKERRQ(ierr);
  ierr = VecDuplicate(scopflow->X,&scopflow->Xu);CHKERRQ(ierr);

  /* Create vectors for constraints */
  /* Vector for coupling constraints */
  ierr = VecCreate(PETSC_COMM_SELF,&scopflow->Di);CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->Di,scopflow->nd,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(scopflow->Di);CHKERRQ(ierr);

  /* Vector for opflow + coupling constraints for one scenario */
  ierr = VecCreate(PETSC_COMM_SELF,&scopflow->Ci);CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->Ci,scopflow->nd+scopflow->opflow[0]->m,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(scopflow->Ci);CHKERRQ(ierr);

  /* Vector for opflow + coupling constraints for ns scenarios */
  if(!scopflow->comm->rank) cloc_size = scopflow->ns*scopflow->opflow[0]->m + (scopflow->ns-1)*scopflow->nd;
  else cloc_size =  scopflow->ns*(scopflow->opflow[0]->m + scopflow->nd);
  ierr = VecCreate(PETSC_COMM_SELF,&scopflow->Cloc);CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->Cloc,cloc_size,cloc_size);CHKERRQ(ierr);
  ierr = VecSetFromOptions(scopflow->Cloc);CHKERRQ(ierr);

  /* Global parallel vector for constraints for Ns scenarios */
  ierr = VecCreate(scopflow->comm->type,&scopflow->C);CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->C,cloc_size,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(scopflow->C);CHKERRQ(ierr);

  /* Constraint bounds vectors  */
  ierr = VecDuplicate(scopflow->C,&scopflow->Cl);CHKERRQ(ierr);
  ierr = VecDuplicate(scopflow->C,&scopflow->Cu);CHKERRQ(ierr);
  
  scopflow->setupcalled = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/*
  SCOPFLOWObjectiveFunction - The objective function for the security constrained optimal power flow

  Input Parameters:
+ scopflow - the SCOPFLOW object
. X      - the current iterate

  Output Parameters:
. obj - the objective function value (scalar)
*/
PetscErrorCode SCOPFLOWObjectiveFunction(SCOPFLOW scopflow,Vec X, PetscScalar* obj)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

Bool eval_scopflow_f(PetscInt n, PetscScalar* x, Bool new_x,
            PetscScalar* obj_value, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  SCOPFLOW         scopflow=(SCOPFLOW)user_data;
  
  ierr = VecPlaceArray(scopflow->X,x);CHKERRQ(ierr);
  ierr = SCOPFLOWObjectiveFunction(scopflow,scopflow->X,obj_value);CHKERRQ(ierr);
  ierr = VecResetArray(scopflow->X);CHKERRQ(ierr);
  return TRUE;
}

/*
  SCOPFLOWObjGradientFunction - The gradient of the objective function for the security constrained optimal power flow

  Input Parameters:
+ scopflow - the SCOPFLOW object
. X      - the current iterate

  Output Parameters:
. grad - the objective function gradient
*/
PetscErrorCode SCOPFLOWObjGradientFunction(SCOPFLOW scopflow,Vec X, Vec grad)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

Bool eval_scopflow_grad_f(PetscInt n, PetscScalar* x, Bool new_x,
                 PetscScalar* grad_f, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  SCOPFLOW         scopflow=(SCOPFLOW)user_data;
  
  ierr = VecPlaceArray(scopflow->X,x);CHKERRQ(ierr);
  ierr = VecPlaceArray(scopflow->gradobj,grad_f);CHKERRQ(ierr);
  ierr = SCOPFLOWObjGradientFunction(scopflow,scopflow->X,scopflow->gradobj);CHKERRQ(ierr);
  ierr = VecResetArray(scopflow->gradobj);CHKERRQ(ierr);
  ierr = VecResetArray(scopflow->X);CHKERRQ(ierr);
  return TRUE;
}

/*
  SCOPFLOWConstraintFunction - Evalulates the constraints for the security constrained optimal power flow

  Input Parameters:
+ scopflow - the SCOPFLOW object
. X      - the current iterate

  Output Parameters:
. G  - vector of constraints
*/
PetscErrorCode SCOPFLOWConstraintFunction(SCOPFLOW scopflow,Vec X,Vec C)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

Bool eval_scopflow_g(PetscInt n, PetscScalar* x, Bool new_x,
             PetscInt m, PetscScalar* c, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  SCOPFLOW       scopflow=(SCOPFLOW)user_data;

  ierr = VecPlaceArray(scopflow->X,x);CHKERRQ(ierr);
  ierr = VecPlaceArray(scopflow->C,c);CHKERRQ(ierr);
  ierr = SCOPFLOWConstraintFunction(scopflow,scopflow->X,scopflow->C);CHKERRQ(ierr);
  ierr = VecResetArray(scopflow->X);CHKERRQ(ierr);
  ierr = VecResetArray(scopflow->C);CHKERRQ(ierr);
  return TRUE;
}

/*
  SCOPFLOWSetConstraintJacobianLocations - Sets the row and column nonzero locations for the
              constraint Jacobian

  Input Parameters:
+ scopflow - the SCOPFLOW object

  Output Parameters:
+ row - the row locations
- col - the col locations

  Notes:
   The row and column locations should be such that J(row[i],col[i]) = val
*/
PetscErrorCode SCOPFLOWSetConstraintJacobianLocations(SCOPFLOW scopflow, PetscInt *row, PetscInt *col)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

/*
  SCOPFLOWSetConstraintJacobianValues - Sets the nonzero values for the
              constraint Jacobian

  Input Parameters:
+ scopflow - the SCOPFLOW object
- X      - the current iterate


  Output Parameters:
. values - the nonzero values in the constraint jacobian

  Notes:
   The values should be in the same row and col order as set in SCOPFLOWSetConstraintJacobianLocations
*/
PetscErrorCode SCOPFLOWSetConstraintJacobianValues(SCOPFLOW scopflow, Vec X,PetscScalar *values)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  PetscFunctionReturn(0);
}


Bool eval_scopflow_jac_g(PetscInt n, PetscScalar *x, Bool new_x,
                PetscInt m, PetscInt nele_jac,
                PetscInt *iRow, PetscInt *jCol, PetscScalar *values,
                UserDataPtr user_data)
{
  PetscErrorCode ierr;
  SCOPFLOW         scopflow=(SCOPFLOW)user_data;

  ierr = VecPlaceArray(scopflow->X,x);CHKERRQ(ierr);
  if(values == NULL) {
    ierr = SCOPFLOWSetConstraintJacobianLocations(scopflow,iRow,jCol);CHKERRQ(ierr);
  } else {
    ierr = SCOPFLOWSetConstraintJacobianValues(scopflow,scopflow->X,values);CHKERRQ(ierr);
  }
  ierr = VecResetArray(scopflow->X);CHKERRQ(ierr);
  return TRUE;
}

Bool eval_scopflow_h(PetscInt n, PetscScalar *x, Bool new_x, PetscScalar obj_factor,
                PetscInt m, PetscScalar *lambda, Bool new_lambda, 
                PetscInt nele_hess, PetscInt *iRow, PetscInt *jCol, 
                PetscScalar *values, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  SCOPFLOW         scopflow=(SCOPFLOW)user_data;

  ierr = VecPlaceArray(scopflow->X,x);CHKERRQ(ierr);
  ierr = VecPlaceArray(scopflow->lambda_g,lambda);CHKERRQ(ierr);
  
  if(values == NULL) {
    ierr = SCOPFLOWSetLagrangianHessianLocations(scopflow,iRow,jCol);CHKERRQ(ierr);
  } else {
    ierr = SCOPFLOWSetLagrangianHessianValues(scopflow,obj_factor, scopflow->X,scopflow->lambda_g,values);CHKERRQ(ierr);
  }
  ierr = VecResetArray(scopflow->X);CHKERRQ(ierr);
  ierr = VecResetArray(scopflow->lambda_g);CHKERRQ(ierr);
  
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
  COMM           comm=scopflow->comm;

  PetscFunctionBegin;
  scopflow->Ns = Ns;
  ierr = PetscOptionsGetInt(NULL,NULL,"-scopflow_Ns",&scopflow->Ns,NULL);CHKERRQ(ierr);
  scopflow->Ns += 1; /* One addition for the base case */
  scopflow->ns = scopflow->Ns / comm->size;
  if(comm->rank < (scopflow->Ns % comm->size)) scopflow->ns++;
  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"SCOPFLOW running with %d scenarios (base case + %d scenarios)\n",scopflow->Ns,scopflow->Ns-1);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF,"Rank[%d] has %d scenarios\n",comm->rank,scopflow->ns);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


  
