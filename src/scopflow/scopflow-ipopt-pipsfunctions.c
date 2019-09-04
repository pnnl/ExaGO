#include <private/scopflowipoptimpl.h>
#include <private/psimpl.h>
#include <private/opflowimpl.h>
#include <math.h>

/* This file contains functions from the scopflow pips
   implementation.
*/



/* Bounds on coupling constraints */
PetscErrorCode SCOPFLOWAddCouplingConstraintBounds(SCOPFLOW scopflow,int row,Vec Gl,Vec Gu)
{
  PetscErrorCode ierr;
  PetscScalar    *gl,*gu;
  PS             ps = scopflow->opflows[row]->ps; /* PS for the scenario */
  PetscInt       gloc; /* starting location for coupling constraints in G vector */
  PetscInt       i;
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
      /* Generator can do a full ramp up to its max. capacity */
      gl[gloc] = PETSC_NINFINITY;//-gen->pt;
      gu[gloc] = PETSC_INFINITY;//gen->pt;
      gloc += 1;
    }
  }

  ierr = VecRestoreArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gu,&gu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


/*
  SCOPFLOWSetVariableandConstraintBounds - Sets the bounds on variables and constraints

  Input Parameters:
+ scopflow - the SCOPFLOW object
- row      - scenario number

  Output Parameters:
+ Xl     - vector of lower bound on variables
. Xu     - vector of upper bound on variables
. Gl     - vector of lower bound on constraints
. Gu     - vector of lower bound on constraints
*/
PetscErrorCode SCOPFLOWSetVariableandConstraintBounds(SCOPFLOW scopflow, PetscInt row, Vec Xl, Vec Xu, Vec Gl, Vec Gu)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=scopflow->opflows[row];
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
    //    xl[loc] = PETSC_NINFINITY; xu[loc] = PETSC_INFINITY;
    xl[loc] = -1e6; xu[loc] = 1e6;
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
    gloc += 2;
  }

  if(!scopflow->ignore_line_flow_constraints) {
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
PetscErrorCode SCOPFLOWSetUp_OPFLOW(SCOPFLOW scopflow,PetscInt row)
{
  PetscErrorCode ierr;
  PetscInt       ngen;
  PetscInt       ncoup;
  OPFLOW         opflow=scopflow->opflows[row];
  PetscInt       i;

  PetscFunctionBegin;
  /* Set up PS object */
  ierr = PSSetUp(opflow->ps);CHKERRQ(ierr);

  /* Set contingencies */
  if(scopflow->ctgcfileset) {
    Contingency ctgc=scopflow->ctgclist.cont[row];
    for(i=0; i < ctgc.noutages; i++) {
      if(ctgc.outagelist[i].type == GEN_OUTAGE) {
	PetscInt gbus=ctgc.outagelist[i].bus;
	char     *gid = ctgc.outagelist[i].id;
	PetscInt status = ctgc.outagelist[i].status;
	ierr = PSSetGenStatus(scopflow->opflows[row]->ps,gbus,gid,status);CHKERRQ(ierr);
      }
      if(ctgc.outagelist[i].type == BR_OUTAGE) {
	PetscInt fbus=ctgc.outagelist[i].fbus;
	PetscInt tbus=ctgc.outagelist[i].tbus;
	char     *brid = ctgc.outagelist[i].id;
	PetscInt status = ctgc.outagelist[i].status;
	ierr = PSSetLineStatus(scopflow->opflows[row]->ps,fbus,tbus,brid,status);CHKERRQ(ierr);
      }
    }
  }

  ierr = PSGetNumGenerators(opflow->ps,&ngen,NULL);CHKERRQ(ierr);  

  if(scopflow->iscoupling) ncoup = (row == 0) ? 0:ngen;
  else ncoup = 0;

  opflow->nconeq   = 2*opflow->ps->Nbus;
  opflow->Nconeq   = 2*opflow->ps->Nbus;
  if(scopflow->ignore_line_flow_constraints) {
    opflow->nconineq = 0;
    opflow->Nconineq = 0;
  } else {
    opflow->nconineq = 2*opflow->ps->nbranch;
    opflow->Nconineq = 2*opflow->ps->Nbranch; /* 0 <= Sf^2 <= Smax^2, 0 <= St^2 <= Smax^2 */
  }

  /* Vector to hold solution */
  ierr = PSCreateGlobalVector(opflow->ps,&opflow->X);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->localX);CHKERRQ(ierr);
  ierr = VecGetSize(opflow->X,&opflow->Nvar);CHKERRQ(ierr);

  /* Vector to hold gradient */
  ierr = VecDuplicate(opflow->X,&opflow->gradobj);CHKERRQ(ierr);

  /* Create the vector for upper and lower bounds on X */
  ierr = VecDuplicate(opflow->X,&opflow->Xl);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->Xu);CHKERRQ(ierr);

  /* Create the equality constraint vector */
  ierr = VecCreate(opflow->comm->type,&opflow->Ge);CHKERRQ(ierr);
  ierr = VecSetSizes(opflow->Ge,opflow->nconeq,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(opflow->Ge);CHKERRQ(ierr);


  if(opflow->nconineq+ncoup) {
    /* Create the inequality constraint vector */
    ierr = VecCreate(opflow->comm->type,&opflow->Gi);CHKERRQ(ierr);
    ierr = VecSetSizes(opflow->Gi,opflow->nconineq+ncoup,PETSC_DETERMINE);CHKERRQ(ierr);
    ierr = VecSetFromOptions(opflow->Gi);CHKERRQ(ierr);
  }

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

  if(opflow->nconineq+ncoup) {
    /* Create Inequality Constraint Jacobian */
    ierr = MatCreate(opflow->comm->type,&opflow->Jac_Gi);CHKERRQ(ierr);
    ierr = MatSetSizes(opflow->Jac_Gi,opflow->nconineq+ncoup,opflow->Nvar,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
    ierr = MatSetType(opflow->Jac_Gi,MATSEQAIJ);CHKERRQ(ierr);
    ierr = MatSetUp(opflow->Jac_Gi);CHKERRQ(ierr);
    ierr = MatSetFromOptions(opflow->Jac_Gi);CHKERRQ(ierr);
  }

  /* Create Hessian */
  ierr = MatCreate(opflow->comm->type,&opflow->Hes);CHKERRQ(ierr);
  ierr = MatSetSizes(opflow->Hes,opflow->Nvar,opflow->Nvar,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = MatSetType(opflow->Hes,MATSEQAIJ);CHKERRQ(ierr);
  ierr = MatSetUp(opflow->Hes);CHKERRQ(ierr);
  ierr = MatSetFromOptions(opflow->Hes);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

extern PetscErrorCode OPFLOWSetInitialGuess(OPFLOW,Vec);

int str_init_x0(double* x0, CallBackDataPtr cbd) {
  int row = cbd->row_node_id;
  SCOPFLOW scopflow=cbd->prob;
  OPFLOW   opflow = scopflow->opflows[row];
  PetscErrorCode ierr;
  
  ierr = VecPlaceArray(opflow->X,x0);CHKERRQ(ierr);
  ierr = OPFLOWSetInitialGuess(opflow,opflow->X);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
	
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
  Vec            Xl,Xu,Gl,Gu;
  PetscInt       i;

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
    ierr = SCOPFLOWSetUp_OPFLOW(scopflow,row);CHKERRQ(ierr);

    if(scopflow->iscoupling && row > 0 && !scopflow->Jcoup) {
      /* Create the constraint Jacobian for coupling */
      ierr = MatDuplicate(scopflow->opflows[row]->Jac_Gi,MAT_DO_NOT_COPY_VALUES,&scopflow->Jcoup);CHKERRQ(ierr);
    
      PS    ps=scopflow->opflows[row]->ps;
      PSBUS bus;
      PetscInt locglob,genctr,k;
      PetscInt rowid[2],colid[2];
      PetscScalar val[2];
      PetscInt    gloc=scopflow->opflows[row]->nconineq;
      PSGEN       gen;
      for(i=0; i < ps->nbus; i++) {
	bus = &ps->bus[i];
	
	ierr = PSBUSGetVariableGlobalLocation(bus,&locglob);CHKERRQ(ierr);
	genctr = 0;
	for(k=0; k < bus->ngen; k++) {
	  ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	  //	  if(!gen->status) continue;
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
    
    ierr = SCOPFLOWSetVariableandConstraintBounds(scopflow,row,Xl,Xu,Gl,Gu);CHKERRQ(ierr);

    
    if(scopflow->iscoupling && row != 0) {
      /* Adding bounds on coupling constraints between first-stage and scenarios */
      ierr = SCOPFLOWAddCouplingConstraintBounds(scopflow,row,Gl,Gu);CHKERRQ(ierr);
    }
    
    ierr = VecResetArray(Xl);CHKERRQ(ierr);
    ierr = VecResetArray(Xu);CHKERRQ(ierr);
    ierr = VecResetArray(Gl);CHKERRQ(ierr);
    ierr = VecResetArray(Gu);CHKERRQ(ierr);
  
    /* Copy over col_lb and col_ub to vectors Xl and Xu so
       that initialization of X can work correctly as
       Xinital = (Xl + Xu)/2
    */
    PetscScalar *xl,*xu;
    ierr = VecGetArray(Xl,&xl);CHKERRQ(ierr);
    ierr = VecGetArray(Xu,&xu);CHKERRQ(ierr);
    ierr = PetscMemcpy(xl,col_lb,scopflow->opflows[row]->Nvar*sizeof(PetscScalar));CHKERRQ(ierr);
    ierr = PetscMemcpy(xu,col_ub,scopflow->opflows[row]->Nvar*sizeof(PetscScalar));CHKERRQ(ierr);
    ierr = VecRestoreArray(Xl,&xl);CHKERRQ(ierr);
    ierr = VecRestoreArray(Xu,&xu);CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_SELF,"Rank %d: row = %d, col = %d set up problem info\n",rank,row,col);CHKERRQ(ierr);
  }
  return 1;
}
