#include <private/psimpl.h>
#include <private/opflowimpl.h>
#include <private/scopflowpipsimpl.h>
#include <math.h>

extern int str_eval_g(double*,double*,double*,double*,CallBackDataPtr);
extern int str_eval_jac_g(double*,double*,int*,double*,int*,int*,int*,double*,int*,int*,CallBackDataPtr);
extern int str_eval_h(double*, double*, double*, int*, double*,
		      int*, int*, CallBackDataPtr);


/* Bounds on coupling constraints */
PetscErrorCode SCOPFLOWAddCouplingConstraintBounds(SCOPFLOW scopflow,int row,Vec Gl,Vec Gu)
{
  PetscErrorCode ierr;
  PetscScalar    *gl,*gu;
  PS             ps = scopflow->opflows[row]->ps; /* PS for the scenario */
  PetscInt       gloc; /* starting location for coupling constraints in G vector */
  PetscInt       i,k;
  PSBUS          bus;

  PetscFunctionBegin;

  gloc = scopflow->opflows[row]->Nconeq + scopflow->opflows[row]->Nconineq;

  ierr = VecGetArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(Gu,&gu);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    PetscInt k;

    bus = &ps->bus[i];

    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;

      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) {
	gl[gloc] = PETSC_NINFINITY;
	gu[gloc] = PETSC_INFINITY;
      } else {
	/* Generator can do a full ramp up to its max. capacity */
	gl[gloc] = -gen->pt;
	gu[gloc] =  gen->pt;
      }
      gloc += 1;
    }
  }

  ierr = VecRestoreArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gu,&gu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


/*
  OPFLOWSetVariableandConstraintBounds - Sets the bounds on variables and constraints

  Input Parameters:
. opflow - the OPFLOW object

  Output Parameters:
+ Xl     - vector of lower bound on variables
. Xu     - vector of upper bound on variables
. Gl     - vector of lower bound on constraints
. Gu     - vector of lower bound on constraints
*/
PetscErrorCode OPFLOWSetVariableandConstraintBounds(OPFLOW opflow, Vec Xl, Vec Xu, Vec Gl, Vec Gu)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscScalar    *xl,*xu,*gl,*gu;
  PetscInt       i;
  PSLINE         line;
  PSBUS          bus;
  PetscInt       loc=0,gloc=0;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(Xu,&xu);CHKERRQ(ierr);
  ierr = VecGetArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(Gu,&gu);CHKERRQ(ierr);
  
  for(i=0; i < ps->nbus; i++) {
    PetscInt k;

    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);

    /* Bounds on voltage angles and bounds on real power mismatch equality constraints */
    xl[loc] = PETSC_NINFINITY; xu[loc] = PETSC_INFINITY;
    gl[gloc] = 0.0;   gu[gloc] = 0.0;

    /* Bounds on voltage magnitudes and bounds on reactive power mismatch equality constraints */
    xl[loc+1] = bus->Vmin; xu[loc+1] = bus->Vmax;
    gl[gloc+1] = 0.0;       gu[gloc+1] = 0.0;

    if(bus->ide == REF_BUS || bus->ide == ISOLATED_BUS) xl[loc] = xu[loc] = bus->va*PETSC_PI/180.0;
    if(bus->ide == ISOLATED_BUS) xl[loc+1] = xu[loc+1] = bus->vm;
    
    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;

      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      loc = loc+2;
      if(!gen->status) xl[loc] = xu[loc] = xl[loc+1] = xu[loc+1] = 0.0;
      else {
	xl[loc] = gen->pb; /* PGmin */
	xu[loc] = gen->pt; /* PGmax */
	xl[loc+1] = gen->qb; /* QGmin */
	xu[loc+1] = gen->qt; /* QGmax */
	/* pb, pt, qb, qt are converted in p.u. in ps.c */
      }
    }
    loc  += 2;
    gloc += 2;
  }

  for(i=0; i < ps->Nbranch; i++) {

    line = &ps->line[i];

    /* Line flow inequality constraints */
    if(!line->status) gl[gloc] = gu[gloc] = gl[gloc+1] = gu[gloc+1] = 0.0;
    else {
      gl[gloc] = gl[gloc+1] = 0.0; 
      gu[gloc] = gu[gloc+1] = (line->rateA/ps->MVAbase)*(line->rateA/ps->MVAbase);
    }    
    gloc += 2;
  }


  ierr = VecRestoreArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Xu,&xu);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gu,&gu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* Sets the underlying OPFLOW object for each SCOPFLOW scenario
   and first-stage 
*/
PetscErrorCode SCOPFLOWSetUp_OPFLOW(OPFLOW opflow,PetscInt row)
{
  PetscErrorCode ierr;
  PetscInt       ngen;
  PetscInt       ncoup;

  PetscFunctionBegin;
  /* Set up PS object */
  ierr = PSSetUp(opflow->ps);CHKERRQ(ierr);
  ierr = PSGetNumGenerators(opflow->ps,&ngen,NULL);CHKERRQ(ierr);  

  ncoup = (row == 0) ? 0:ngen;

  opflow->nconeq   = 2*opflow->ps->Nbus;
  opflow->Nconeq   = 2*opflow->ps->Nbus;
  opflow->nconineq = 2*opflow->ps->nbranch;
  opflow->Nconineq = 2*opflow->ps->Nbranch; /* 0 <= Sf^2 <= Smax^2, 0 <= St^2 <= Smax^2 */

  /* Vector to hold solution */
  ierr = PSCreateGlobalVector(opflow->ps,&opflow->X);CHKERRQ(ierr);
  ierr = VecGetSize(opflow->X,&opflow->Nvar);CHKERRQ(ierr);

  /* Create the vector for upper and lower bounds on X */
  ierr = VecDuplicate(opflow->X,&opflow->Xl);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->Xu);CHKERRQ(ierr);

  /* Create the equality constraint vector */
  ierr = VecCreate(opflow->comm->type,&opflow->Ge);CHKERRQ(ierr);
  ierr = VecSetSizes(opflow->Ge,opflow->nconeq,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(opflow->Ge);CHKERRQ(ierr);

  /* Create the inequality constraint vector */
  ierr = VecCreate(opflow->comm->type,&opflow->Gi);CHKERRQ(ierr);
  ierr = VecSetSizes(opflow->Gi,opflow->nconineq,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(opflow->Gi);CHKERRQ(ierr);

  /* Create lower and upper bounds on constraint vector */
  /* Size Nconeq + Nconineq (+Ncoupling for row != 0 only) */
  ierr = VecCreate(opflow->comm->type,&opflow->Gl);CHKERRQ(ierr);
  ierr = VecSetSizes(opflow->Gl,opflow->nconineq+opflow->nconeq+ncoup,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(opflow->Gl);CHKERRQ(ierr);
  
  ierr = VecDuplicate(opflow->Gl,&opflow->Gu);CHKERRQ(ierr);

  /* Create Lagrangian multiplier vector */
  ierr = VecDuplicate(opflow->Gl,&opflow->Lambda);CHKERRQ(ierr);

  ierr = DMNetworkSetVertexLocalToGlobalOrdering(opflow->ps->networkdm);CHKERRQ(ierr);

  /* Create Equality Constraint Jacobian */
  ierr = MatCreate(opflow->comm->type,&opflow->Jac_Ge);CHKERRQ(ierr);
  ierr = MatSetSizes(opflow->Jac_Ge,opflow->nconeq,opflow->Nvar,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = MatSetType(opflow->Jac_Ge,MATSEQAIJ);CHKERRQ(ierr);
  ierr = MatSetUp(opflow->Jac_Ge);CHKERRQ(ierr);
  ierr = MatSetFromOptions(opflow->Jac_Ge);CHKERRQ(ierr);

  /* Create Inequality Constraint Jacobian */
  ierr = MatCreate(opflow->comm->type,&opflow->Jac_Gi);CHKERRQ(ierr);
  ierr = MatSetSizes(opflow->Jac_Gi,opflow->nconineq+ncoup,opflow->Nvar,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = MatSetType(opflow->Jac_Gi,MATSEQAIJ);CHKERRQ(ierr);
  ierr = MatSetUp(opflow->Jac_Gi);CHKERRQ(ierr);
  ierr = MatSetFromOptions(opflow->Jac_Gi);CHKERRQ(ierr);

  /* Create Hessian */
  ierr = MatCreate(opflow->comm->type,&opflow->Hes);CHKERRQ(ierr);
  ierr = MatSetSizes(opflow->Hes,opflow->Nvar,opflow->Nvar,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = MatSetType(opflow->Hes,MATSEQAIJ);CHKERRQ(ierr);
  ierr = MatSetUp(opflow->Hes);CHKERRQ(ierr);
  ierr = MatSetFromOptions(opflow->Hes);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

int str_init_x0(double* x0, CallBackDataPtr cbd) {
	int row = cbd->row_node_id;
	int col = cbd->col_node_id;
	if(row == 0)
	{
		x0[0] = 1.0;
		x0[1] = 1.0;
	}
	else if(row == 1)
	{
		x0[0] = 1.0;
	}
	else if(row == 2)
	{
		x0[0] = 1.0;
	}

	return 1;
}

int str_prob_info(int* n, double* col_lb, double* col_ub, int* m,
		  double* row_lb, double* row_ub, CallBackDataPtr cbd) {
  int row = cbd->row_node_id;
  int col = cbd->col_node_id;
  SCOPFLOW scopflow = (SCOPFLOW)cbd->prob;
  PetscInt rank=scopflow->comm->rank;
  int type = cbd->typeflag;
  PetscErrorCode ierr;
  PetscInt       ngen;
  Vec            Xl,Xu,Gl,Gu;

  if(type == 1) {
    if(row_lb == NULL){
      *m = 0;
    }
    return 1;
  }
  
  /* Set sizes of variables and constraints */
  if(col_lb == NULL){
    /* Create OPFLOW */
    ierr = OPFLOWCreate(PETSC_COMM_SELF,&scopflow->opflows[row]);CHKERRQ(ierr);
    ierr = OPFLOWReadMatPowerData(scopflow->opflows[row],scopflow->netfile);CHKERRQ(ierr);
    /* SCOPFLOWSetUp_OPFLOW handles creating the vectors and the sizes of the constraints */
    ierr = SCOPFLOWSetUp_OPFLOW(scopflow->opflows[row],row);CHKERRQ(ierr);
    if(row == 1) {
      /* Create the constraint Jacobian for coupling */
      ierr = MatDuplicate(scopflow->opflows[row]->Jac_Gi,MAT_DO_NOT_COPY_VALUES,&scopflow->Jcoup);CHKERRQ(ierr);
    
      PS    ps=scopflow->opflows[row]->ps;
      PSBUS bus;
      PetscInt locglob,genctr,k,i;
      PetscInt rowid[2],colid[2];
      PetscScalar val[2];
      PetscInt    gloc=scopflow->opflows[row]->nconineq;
      for(i=0; i < ps->nbus; i++) {
	bus = &ps->bus[i];
	
	ierr = PSBUSGetVariableGlobalLocation(bus,&locglob);CHKERRQ(ierr);
	genctr = 0;
	for(k=0; k < bus->ngen; k++) {
	  val[0] = -1;
	  rowid[0] = gloc;
	  colid[0] = locglob + 2 + genctr;
	  ierr = MatSetValues(scopflow->Jcoup,1,rowid,1,colid,val,ADD_VALUES);CHKERRQ(ierr);
	  genctr += 2;
	  gloc += 1;
	}
      }
      
      ierr = MatAssemblyBegin(scopflow->Jcoup,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(scopflow->Jcoup,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatTranspose(scopflow->Jcoup,MAT_INITIAL_MATRIX,&scopflow->JcoupT);CHKERRQ(ierr);
    }

    ierr = VecGetSize(scopflow->opflows[row]->X,n);CHKERRQ(ierr);
    ierr = VecGetSize(scopflow->opflows[row]->Gl,m);CHKERRQ(ierr);
    scopflow->ns++;

    ierr = PetscPrintf(PETSC_COMM_SELF,"Rank %d row %d col %d type %d Nvar=%d Ncon = %d\n",rank,row,col,type,*n,*m);CHKERRQ(ierr);
  } else {
    /* Bounds on variables and constraints */
    Xl = scopflow->opflows[row]->Xl;
    Xu = scopflow->opflows[row]->Xu;
    Gl = scopflow->opflows[row]->Gl;
    Gu = scopflow->opflows[row]->Gu;

    ierr = VecPlaceArray(Xl,col_lb);CHKERRQ(ierr);
    ierr = VecPlaceArray(Xu,col_ub);CHKERRQ(ierr);
    ierr = VecPlaceArray(Gl,row_lb);CHKERRQ(ierr);
    ierr = VecPlaceArray(Gu,row_ub);CHKERRQ(ierr);
    
    ierr = OPFLOWSetVariableandConstraintBounds(scopflow->opflows[row],Xl,Xu,Gl,Gu);CHKERRQ(ierr);

    if(row != 0) {
      /* Adding bounds on coupling constraints between first-stage and scenarios */
      ierr = SCOPFLOWAddCouplingConstraintBounds(scopflow,row,Gl,Gu);CHKERRQ(ierr);
    }
    
    ierr = VecResetArray(Xl);CHKERRQ(ierr);
    ierr = VecResetArray(Xu);CHKERRQ(ierr);
    ierr = VecResetArray(Gl);CHKERRQ(ierr);
    ierr = VecResetArray(Gu);CHKERRQ(ierr);
    
  }
  return 1;
}

int str_eval_f(double* x0, double* x1, double* obj, CallBackDataPtr cbd) {
	int row = cbd->row_node_id;
	int col = cbd->col_node_id;
	if(row == 0 )
	{   // (x0 + x1) ^ 2 + x0 + x1
		*obj =  ( x0[0] + x0[1] ) * ( x0[0] + x0[1] ) + x0[0] + x0[1];
	}
	else if(row == 1)
	{   // (x0 + x1)x3 + x3
		*obj = ( x0[0] + x0[1] ) * x1[0] + x1[0];
	}
	else if(row == 2)
	{  // (x0 + x1)x4 + x4
		*obj = ( x0[0] + x0[1] ) * x1[0] + x1[0];
	}
	return 1;
}

int str_eval_grad_f(double* x0, double* x1, double* grad, CallBackDataPtr cbd) {
	int row = cbd->row_node_id;
	int col = cbd->col_node_id;
	if(row == 0 && col == 0)
	{
		grad[0] = 2.0 * (x0[0] + x0[1]) + 1.0;
		grad[1] = 2.0 * (x0[0] + x0[1]) + 1.0;
	}
	else if(row == 1 && col == 1)
	{
		grad[0] = (x0[0] + x0[1]) + 1.0;
	}
	else if(row == 2 && col == 2)
	{
		grad[0] = (x0[0] + x0[1]) + 1.0;
	}
	else if(row == 1 && col == 0)
	{
		grad[0] = x1[0];
		grad[1] = x1[0];
	}
	else if(row == 2 && col == 0)
	{
		grad[0] = x1[0];
		grad[1] = x1[0];
	}

	return 1;
}

int str_write_solution(double* x, double* lam_eq, double* lam_ieq,CallBackDataPtr cbd)
{
	int row = cbd->row_node_id;
	int col = cbd->col_node_id;
	if(row == 0)
	{
	}
	else if(row == 1 || row == 2)
	{
	}
	return 1;
}

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

  scopflow->ns              =  0;
  scopflow->Ns              = -1;
  scopflow->nc              =  0;
  scopflow->Nc              =  0;
  scopflow->nx              =  0;
  scopflow->Nx              =  0;

  scopflow->nlp_pips       = NULL;

  //  gmyid = scopflow->comm->rank;

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

  //  FreeIpoptProblem((*scopflow)->nlp_ipopt);CHKERRQ(ierr);
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
  if(scopflow->Ns == -1) {
    SETERRQ(scopflow->comm->type,0,"Must call SCOPFLOWSetNumScenarios or SCOPFLOWSetScenariosFile before calling this function\n");
  }

  ierr = PetscMemcpy(scopflow->netfile,netfile,100*sizeof(char));CHKERRQ(ierr);

  if(scopflow->comm->size == 1) { /* Serial */
    ierr = PetscMalloc1(scopflow->Ns+1,&opflows);CHKERRQ(ierr);
  } else {
    ierr = PetscMalloc1(6,&opflows);CHKERRQ(ierr); /* Max. 6 (5 scenarios + 1 first-stage per rank) */
  }

  scopflow->opflows = opflows;
  for(i=0; i < scopflow->ns; i++) {
    ierr = OPFLOWCreate(PETSC_COMM_SELF,&opflows[i]);CHKERRQ(ierr);
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
  OPFLOW         *opflows=scopflow->opflows;

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

  PetscFunctionBegin;
  if(!scopflow->setupcalled) {
    ierr = SCOPFLOWSetUp(scopflow);CHKERRQ(ierr);
  }

  /* Solve with PIPS */
  PipsNlpSolveStruct(scopflow->nlp_pips);

  /*
  ierr = SCOPFLOWSetVariableandConstraintBounds(scopflow,scopflow->Xl,scopflow->Xu,scopflow->Cl,scopflow->Cu);CHKERRQ(ierr);

  ierr = VecGetArray(scopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Cl,&cl);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Cu,&cu);CHKERRQ(ierr);

  ierr = VecRestoreArray(scopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Cl,&cl);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Cu,&cu);CHKERRQ(ierr);

  ierr = SCOPFLOWSetInitialGuess(scopflow,scopflow->X);CHKERRQ(ierr);

  Vec C;
  ierr = VecDuplicate(scopflow->C,&C);CHKERRQ(ierr);
  ierr = SCOPFLOWConstraintFunction(scopflow,scopflow->X,C);CHKERRQ(ierr);

  PetscScalar *x;
  ierr = VecGetArray(scopflow->X,&x);CHKERRQ(ierr);

  //  scopflow->solve_status = IpoptSolve(scopflow->nlp_ipopt,x,NULL,&scopflow->obj,NULL,NULL,NULL,scopflow);

  ierr = VecRestoreArray(scopflow->X,&x);CHKERRQ(ierr);

  */
  if(!scopflow->comm->rank) {
    ierr = PetscPrintf(PETSC_COMM_SELF,"Objective function value = %lf\n",scopflow->obj);CHKERRQ(ierr);
  }

  //  ierr = VecView(scopflow->X,0);CHKERRQ(ierr);
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
  OPFLOW         *opflows=scopflow->opflows;
  PetscInt       ngen; /* Number of generators */
  PetscInt       i,cloc_size;

  PetscFunctionBegin;
  scopflow->nx = 0;
  /* Set up OPFLOW objects for scenarios */
  for(i=0; i < scopflow->ns; i++) {
    ierr = OPFLOWSetUp(opflows[i]);CHKERRQ(ierr);
  }

  str_init_x0_cb init_x0 = &str_init_x0;
  str_prob_info_cb prob_info = &str_prob_info;
  str_eval_f_cb eval_f = &str_eval_f;
  str_eval_g_cb eval_g = &str_eval_g;
  str_eval_grad_f_cb eval_grad_f = &str_eval_grad_f;
  str_eval_jac_g_cb eval_jac_g = &str_eval_jac_g;
  str_eval_h_cb eval_h = &str_eval_h;
  str_write_solution_cb write_solution = &str_write_solution;

  scopflow->nlp_pips = CreatePipsNlpProblemStruct(MPI_COMM_WORLD, scopflow->Ns,
							    init_x0, prob_info, eval_f, eval_g, eval_grad_f, eval_jac_g,
						  eval_h, str_write_solution, (UserDataPtr)scopflow);
  /*
  if(scopflow->ns) {
    scopflow->nx = opflows[0]->n;
    ierr = PSGetNumGenerators(opflows[0]->ps,&ngen,NULL);CHKERRQ(ierr);
    scopflow->nd = ngen;
  }
  ierr = VecCreate(PETSC_COMM_SELF,&scopflow->Xloc);CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->Xloc,scopflow->ns*scopflow->nx,scopflow->ns*scopflow->nx);CHKERRQ(ierr);
  ierr = VecSetFromOptions(scopflow->Xloc);CHKERRQ(ierr);

  ierr = VecDuplicate(scopflow->opflow[0]->X,&scopflow->X0);CHKERRQ(ierr);

  ierr = VecCreate(scopflow->comm->type,&scopflow->X);CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->X,scopflow->ns*scopflow->nx,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(scopflow->X);CHKERRQ(ierr);

  ierr = VecDuplicate(scopflow->X,&scopflow->Xl);CHKERRQ(ierr);
  ierr = VecDuplicate(scopflow->X,&scopflow->Xu);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_SELF,&scopflow->Di);CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->Di,scopflow->nd,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(scopflow->Di);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_SELF,&scopflow->Ci);CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->Ci,scopflow->nd+scopflow->opflow[0]->m,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(scopflow->Ci);CHKERRQ(ierr);

  if(!scopflow->comm->rank) cloc_size = scopflow->ns*scopflow->opflow[0]->m + (scopflow->ns-1)*scopflow->nd;
  else cloc_size =  scopflow->ns*(scopflow->opflow[0]->m + scopflow->nd);
  ierr = VecCreate(PETSC_COMM_SELF,&scopflow->Cloc);CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->Cloc,cloc_size,cloc_size);CHKERRQ(ierr);
  ierr = VecSetFromOptions(scopflow->Cloc);CHKERRQ(ierr);

  ierr = VecCreate(scopflow->comm->type,&scopflow->C);CHKERRQ(ierr);
  ierr = VecSetSizes(scopflow->C,cloc_size,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(scopflow->C);CHKERRQ(ierr);

  
  ierr = VecDuplicate(scopflow->C,&scopflow->Cl);CHKERRQ(ierr);
  ierr = VecDuplicate(scopflow->C,&scopflow->Cu);CHKERRQ(ierr);
  */

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
  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"SCOPFLOW running with %d scenarios\n",scopflow->Ns);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


  
