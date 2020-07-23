#include <scopflow_config.h>
#if defined(EXAGO_HAVE_HIOP)
#include <private/opflowimpl.h>
#include "opflow-hiop.hpp"
/* Converts an array xin in natural ordering to an array xout in sparse-dense
   ordering
*/
void OPFLOWHIOPInterface::naturaltospdense(const double *xin,double *xout)
{
  int i;

  for(i=0; i < opflow->nx; i++) {
    xout[idxn2sd_map[i]] = xin[i];
  }
}

/* Converts an array xin in sparse dense ordering to an array xout in natural
   ordering
*/
void OPFLOWHIOPInterface::spdensetonatural(const double *xin,double *xout)
{
  int i;

  for(i=0; i < opflow->nx; i++) {
    xout[i] = xin[idxn2sd_map[i]];
  }
}

OPFLOWHIOPInterface::OPFLOWHIOPInterface(OPFLOW opflowin) 
{
  PS ps;
  PSBUS bus;
  int   ngen;
  PSGEN gen;
  
  opflow   = opflowin;
  
  ps = opflow->ps;
  nxsparse = 2*ps->ngenON;
  nxdense  = 2*ps->nbus;
  
  PetscMalloc1(opflow->nx,&idxn2sd_map);
  
  int i,k;
  int spct=0,dnct=0;
  int loc;
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    PSBUSGetVariableLocation(bus,&loc);

    idxn2sd_map[loc] = nxsparse + dnct;
    idxn2sd_map[loc+1] = nxsparse + dnct+1;

    dnct += 2;
    loc += 2;
    PSBUSGetNGen(bus,&ngen);
    for(k=0; k < ngen; k++) {
      PSBUSGetGen(bus,k,&gen);
      if(!gen->status) continue;

      idxn2sd_map[loc] = spct;
      idxn2sd_map[loc+1] = spct + 1;

      spct += 2;
      loc += 2;
    }
  }
  
  /*  
  for(i=0; i < opflow->nx; i++) {
    PetscPrintf(PETSC_COMM_SELF,"Variable[%d]: natural =%d\tn2sd=%d\n",i,i,idxn2sd_map[i]);
  }
  */
  
}

bool OPFLOWHIOPInterface::get_prob_sizes(long long& n, long long& m)
{ 
  n = opflow->nx;
  m = opflow->ncon;
  return true; 
}

bool OPFLOWHIOPInterface::get_vars_info(const long long& n, double *xlow, double* xupp, NonlinearityType* type)
{
  PetscInt       i;
  PetscErrorCode ierr;
  const PetscScalar    *xl,*xu;

  ierr = (*opflow->modelops.setvariablebounds)(opflow,opflow->Xl,opflow->Xu);CHKERRQ(ierr);

  ierr = VecGetArrayRead(opflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArrayRead(opflow->Xu,&xu);CHKERRQ(ierr);


  /* Natural to sparse-dense conversion */
  naturaltospdense(xl,xlow);
  naturaltospdense(xu,xupp);

  for(i=0; i < n; i++) {
    type[i] = hiopNonlinear;
  }

  ierr = VecRestoreArrayRead(opflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(opflow->Xu,&xu);CHKERRQ(ierr);

  return true;
}

bool OPFLOWHIOPInterface::get_cons_info(const long long& m, double* clow, double* cupp, NonlinearityType* type)
{
  PetscInt i;
  PetscErrorCode ierr;

  ierr = VecPlaceArray(opflow->Gl,clow);CHKERRQ(ierr);
  ierr = VecPlaceArray(opflow->Gu,cupp);CHKERRQ(ierr);

  ierr = (*opflow->modelops.setconstraintbounds)(opflow,opflow->Gl,opflow->Gu);CHKERRQ(ierr);

  ierr = VecResetArray(opflow->Gl);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->Gu);CHKERRQ(ierr);

  for(i=0; i < m; i++) type[i] = hiopNonlinear;

  return true;
}

bool OPFLOWHIOPInterface::get_sparse_dense_blocks_info(int& nx_sparse, int& nx_dense,
				  int& nnz_sparse_Jace, int& nnz_sparse_Jaci,
				  int& nnz_sparse_Hess_Lagr_SS, int& nnz_sparse_Hess_Lagr_SD)
{
  nx_sparse = nxsparse;
  nx_dense  = nxdense;

  nnz_sparse_Jace = nnz_sparse_Hess_Lagr_SS = nxsparse;
  /* HIOP requires both the equality and inequality constraints to be dependent on both the sparse and dense variables. The inequality constraints not being functions of the sparse variables (Pg,Qg) cause a crash in HIOP.Hence, setting nnz_sparse_Jaci = 1 to add a fake value that will (hopefully) not make HIOP crash 
   */
  nnz_sparse_Jaci = opflow->nconineq?1:0; 
  nnz_sparse_Hess_Lagr_SD = 0;
  return true;
}

bool OPFLOWHIOPInterface::eval_f(const long long& n, const double* x, bool new_x, double& obj_value)
{
  PetscErrorCode ierr;
  PetscScalar    *xarr;

  obj_value = 0.0;
  ierr = VecGetArray(opflow->X,&xarr);CHKERRQ(ierr);
  /* Convert from sparse-dense to natural ordering */
  spdensetonatural(x,xarr);
  ierr = VecRestoreArray(opflow->X,&xarr);CHKERRQ(ierr);

  /* Compute objective */
  ierr = (*opflow->modelops.computeobjective)(opflow,opflow->X,&obj_value);CHKERRQ(ierr);

  return true;
}

bool OPFLOWHIOPInterface::eval_cons(const long long& n, const long long& m, 
	       const long long& num_cons, const long long* idx_cons,  
	       const double* x, bool new_x, double* cons)
{
  PetscErrorCode ierr;
  PetscScalar    *xarr;

  if(!num_cons) return true;

  ierr = VecGetArray(opflow->X,&xarr);CHKERRQ(ierr);
  /* Convert from sparse-dense to natural ordering */
  spdensetonatural(x,xarr);
  ierr = VecRestoreArray(opflow->X,&xarr);CHKERRQ(ierr);

  if(idx_cons[0] == 0) {
    /* Equality constaints */
    ierr = VecPlaceArray(opflow->Ge,cons+idx_cons[0]);CHKERRQ(ierr);
    ierr = (*opflow->modelops.computeequalityconstraints)(opflow,opflow->X,opflow->Ge);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Ge);CHKERRQ(ierr);
  }

  if(idx_cons[0] == opflow->nconeq && opflow->nconineq) {
    /* Inequality constraints */
    ierr = VecPlaceArray(opflow->Gi,cons);CHKERRQ(ierr);
    ierr = (*opflow->modelops.computeinequalityconstraints)(opflow,opflow->X,opflow->Gi);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Gi);CHKERRQ(ierr);
  }

  return true;
}

bool OPFLOWHIOPInterface::eval_grad_f(const long long& n, const double* x, bool new_x, double* gradf)
{
  PetscErrorCode ierr;
  PetscScalar    *xarr;
  const PetscScalar *gradarr;

  ierr = VecGetArray(opflow->X,&xarr);CHKERRQ(ierr);
  /* Convert from sparse-dense to natural ordering */
  spdensetonatural(x,xarr);
  ierr = VecRestoreArray(opflow->X,&xarr);CHKERRQ(ierr);

  ierr = VecSet(opflow->gradobj,0.0);CHKERRQ(ierr);
  ierr = (*opflow->modelops.computegradient)(opflow,opflow->X,opflow->gradobj);CHKERRQ(ierr);

  ierr = VecGetArrayRead(opflow->gradobj,&gradarr);CHKERRQ(ierr);
  /* Convert from natural to sparse-dense ordering */
  naturaltospdense(gradarr,gradf);
  ierr = VecRestoreArrayRead(opflow->gradobj,&gradarr);CHKERRQ(ierr);

  return true;
}

bool OPFLOWHIOPInterface::eval_Jac_cons(const long long& n, const long long& m, 
		   const long long& num_cons, const long long* idx_cons,
		   const double* x, bool new_x,
		   const long long& nsparse, const long long& ndense, 
		   const int& nnzJacS, int* iJacS, int* jJacS, double* MJacS, 
		   double** JacD)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k;
  PetscInt       ncols;
  const PetscInt *cols;
  const PetscScalar *vals;
  int               nnzs=0,gi=0;
  PS             ps=opflow->ps;
  PSBUS          bus;
  PetscScalar    *xarr;
  PetscInt       dcol;

  if(!num_cons) return true;

  if(idx_cons[0] == 0 && iJacS != NULL && jJacS!= NULL) {
    /* Sparse equality constraint Jacobian locations w.r.t Pg,Qg */
    for(i=0; i < ps->nbus; i++) {
      bus = &ps->bus[i];

      for(k=0; k < bus->ngenON; k++) {
	iJacS[nnzs + k] = 2*i;
	jJacS[nnzs + k] = 2*gi;

	iJacS[nnzs + bus->ngenON + k] = 2*i+1;
	jJacS[nnzs + bus->ngenON + k] = 2*gi+1;

	gi += 1;
      }
      nnzs += 2*bus->ngenON;
    }
    if(nnzs != nnzJacS) SETERRQ(PETSC_COMM_SELF,0,"Incorrect number of entries in sparse equality constraint Jacobian\n");
  }

  nnzs = 0;
  if(idx_cons[0] == opflow->nconeq && opflow->nconineq && iJacS != NULL && jJacS!= NULL) {
    /* Dummy entry to help avoid HIOP crash */
    iJacS[nnzs] = 0;
    jJacS[nnzs] = 0;
  }

  if(idx_cons[0] == opflow->nconeq && opflow->nconineq && MJacS != NULL) {
    MJacS[nnzs] = 0.0;
  }
  
  nnzs = 0;
  if(idx_cons[0] == 0 && MJacS != NULL) {
    /* Sparse equality constraint Jacobian values w.r.t Pg,Qg */
    for(i=0; i < ps->nbus; i++) {
      bus = &ps->bus[i];

      for(k=0; k < bus->ngenON; k++) {
	MJacS[nnzs] = -1;
	MJacS[nnzs+1] = -1;

	nnzs += 2;
      }
    }
  }

  if(JacD != NULL) {
    ierr = VecGetArray(opflow->X,&xarr);CHKERRQ(ierr);
    /* Convert from sparse-dense to natural ordering */
    spdensetonatural(x,xarr);
    ierr = VecRestoreArray(opflow->X,&xarr);CHKERRQ(ierr);

    if(idx_cons[0] == 0) {
      /* Equality constraints Jacobian */
      ierr = (*opflow->modelops.computeequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Ge);CHKERRQ(ierr);

      for(i=0; i < opflow->nconeq; i++) {
	for(j=0; j < nxdense; j++) JacD[i][j] = 0.0;
	ierr = MatGetRow(opflow->Jac_Ge,i,&ncols,&cols,&vals);CHKERRQ(ierr);
	for(k=0; k < ncols; k++) {
	  if(idxn2sd_map[cols[k]]-nxsparse >= 0) {
	    /* Dense variables */
	    dcol = idxn2sd_map[cols[k]] - nxsparse; /* Column number for dense variable in the dense block */
	    JacD[i][dcol] = vals[k];
	  }
	}
	ierr = MatRestoreRow(opflow->Jac_Ge,i,&ncols,&cols,&vals);CHKERRQ(ierr);
      }
    } else {
      
      /* Inequality constraints Jacobian */
      ierr = (*opflow->modelops.computeinequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Gi);CHKERRQ(ierr);

      for(i=0; i < opflow->nconineq; i++) {
	for(j=0; j < nxdense; j++) JacD[i][j] = 0.0;
	ierr = MatGetRow(opflow->Jac_Gi,i,&ncols,&cols,&vals);CHKERRQ(ierr);
	for(k=0; k < ncols; k++) {
	  if(idxn2sd_map[cols[k]]-nxsparse >= 0) {
	    /* Dense variables */
	    dcol = idxn2sd_map[cols[k]] - nxsparse; /* Column number for dense variable in the dense block */
	    JacD[i][dcol] = vals[k];
	  }
	}
	ierr = MatRestoreRow(opflow->Jac_Gi,i,&ncols,&cols,&vals);CHKERRQ(ierr);
      }
    }
  }
  return true;
}

bool OPFLOWHIOPInterface::eval_Hess_Lagr(const long long& n, const long long& m, 
		    const double* x, bool new_x, const double& obj_factor,
		    const double* lambda, bool new_lambda,
		    const long long& nsparse, const long long& ndense, 
		    const int& nnzHSS, int* iHSS, int* jHSS, double* MHSS, 
		    double** HDD,
		    int& nnzHSD, int* iHSD, int* jHSD, double* MHSD)
{
  PetscErrorCode ierr;
  PetscInt       i,k;
  PetscInt       ncols;
  const PetscInt *cols;
  const PetscScalar *vals;
  PetscInt       nnzs=0;
  PetscScalar    *xarr;
  PS             ps=opflow->ps;
  PetscInt       dcol,dnct=0;

  opflow->obj_factor = obj_factor;

  if(iHSS != NULL && jHSS!= NULL) {
    for(i=0; i < ps->ngenON; i++) {
      iHSS[nnzs] = 2*i;
      jHSS[nnzs] = 2*i;

      iHSS[nnzs+1] = 2*i+1;
      jHSS[nnzs+1] = 2*i+1;

      nnzs += 2;
    }
    if(nnzHSS != nnzs) SETERRQ2(PETSC_COMM_SELF,0,"Incorrect number of non-zeros for sparse Hessian %d != %d\n",nnzHSS,nnzs);
  }

  ierr = VecGetArray(opflow->X,&xarr);CHKERRQ(ierr);
  /* Convert from sparse-dense to natural ordering */
  spdensetonatural(x,xarr);
  ierr = VecRestoreArray(opflow->X,&xarr);CHKERRQ(ierr);

  ierr = VecPlaceArray(opflow->Lambdae,lambda);CHKERRQ(ierr);
  if(opflow->Nconineq) {
    ierr = VecPlaceArray(opflow->Lambdai,lambda+opflow->nconeq);CHKERRQ(ierr);
  }
    
  /* Compute Hessian */
  ierr = (*opflow->modelops.computehessian)(opflow,opflow->X,opflow->Lambdae,opflow->Lambdai,opflow->Hes);CHKERRQ(ierr);
    
  ierr = VecResetArray(opflow->Lambdae);CHKERRQ(ierr);
  if(opflow->Nconineq) {
    ierr = VecResetArray(opflow->Lambdai);CHKERRQ(ierr);
  }

  nnzs = 0;
  if(MHSS != NULL) {
    for(i=0; i < n; i++) {
      if(idxn2sd_map[i] < nxsparse) {
	/* Rows for sparse variables */
	ierr = MatGetRow(opflow->Hes,i,&ncols,&cols,&vals);CHKERRQ(ierr);
	for(k=0; k < ncols; k++) {
	  if(cols[k] == i) {
	    /* Diagonal element */
	    MHSS[nnzs] = vals[k];
	    nnzs += 1;
	  }
	}
      }
    }
  }

  if(HDD != NULL) {
    for(i=0; i < n; i++) {
      if(idxn2sd_map[i] >= nxsparse) {
	for(k=0; k < nxdense; k++) HDD[dnct][k] = 0.0;
	/* Rows for dense variables */
	ierr = MatGetRow(opflow->Hes,i,&ncols,&cols,&vals);CHKERRQ(ierr);
	for(k=0; k < ncols; k++) {
	  if(idxn2sd_map[cols[k]] >= nxsparse) {
	    dcol = idxn2sd_map[cols[k]] - nxsparse; /* Column number for dense variable in the dense block */
	    HDD[dnct][dcol] = vals[k];
	  }
	}
	ierr = MatRestoreRow(opflow->Hes,i,&ncols,&cols,&vals);CHKERRQ(ierr);

	dnct++;
      }
    }
  }

  return true;
}

bool OPFLOWHIOPInterface::get_starting_point(const long long& global_n, double* x0)
{
  PetscErrorCode ierr;
  const PetscScalar    *xarr;

  /* Set initial guess */
  ierr = (*opflow->modelops.setinitialguess)(opflow,opflow->X);CHKERRQ(ierr);

  ierr = VecGetArrayRead(opflow->X,&xarr);CHKERRQ(ierr);
  /* Convert from natural to sparse-dense ordering */
  naturaltospdense(xarr,x0);
  ierr = VecRestoreArrayRead(opflow->X,&xarr);CHKERRQ(ierr);

  return true;
}

void OPFLOWHIOPInterface::solution_callback(hiop::hiopSolveStatus status,
							 int n, const double* xsol,
							 const double* z_L,
							 const double* z_U,
							 int m, const double* gsol,
							 const double* lamsol,
							 double obj_value)
{
  PetscErrorCode    ierr;
  OPFLOWSolver_HIOP hiop=(OPFLOWSolver_HIOP)opflow->solver;
  PetscScalar       *x,*lam,*g;

  /* Copy over solution details */
  hiop->status = status;
  opflow->obj = obj_value;

  ierr = VecGetArray(opflow->X,&x);CHKERRV(ierr);
  spdensetonatural(xsol,x);
  ierr = VecRestoreArray(opflow->X,&x);CHKERRV(ierr);

  if(lamsol) { 
    /* HIOP returns a NULL for lamsol - probably lamsol needs to be added to HIOP. Need to
     remove this condition once it is fixed 
    */
    ierr = VecGetArray(opflow->Lambda,&lam);CHKERRV(ierr);
    ierr = PetscMemcpy((double*)lamsol,lam,opflow->nconeq*sizeof(PetscScalar));CHKERRV(ierr);
    if(opflow->Nconineq) {
      ierr = PetscMemcpy((double*)(lamsol+opflow->nconeq),lam+opflow->nconeq,opflow->nconineq*sizeof(PetscScalar));CHKERRV(ierr);
    }
    ierr = VecRestoreArray(opflow->Lambda,&lam);CHKERRV(ierr);
  } else {
    ierr = VecSet(opflow->Lambda,-9999.0);CHKERRV(ierr);
  }

  if(gsol) {
    /* Same situation as lamsol - gsol is NULL */
    ierr = VecGetArray(opflow->G,&g);CHKERRV(ierr);
    ierr = PetscMemcpy((double*)gsol,g,opflow->nconeq*sizeof(PetscScalar));CHKERRV(ierr);
    if(opflow->Nconineq) {
      ierr = PetscMemcpy((double*)(gsol+opflow->nconeq),g+opflow->nconeq,opflow->nconineq*sizeof(PetscScalar));CHKERRV(ierr);
    }
    ierr = VecRestoreArray(opflow->G,&g);CHKERRV(ierr);
  }
}

PetscErrorCode OPFLOWSolverSetUp_HIOP(OPFLOW opflow)
{
  PetscErrorCode ierr;
  OPFLOWSolver_HIOP hiop=(OPFLOWSolver_HIOP)opflow->solver;
  PetscBool         flg;
  PetscFunctionBegin;

  hiop->nlp = new OPFLOWHIOPInterface(opflow);
  hiop->mds = new hiop::hiopNlpMDS(*hiop->nlp);

  /* Options set in hiop.options file */

  hiop->solver = new hiop::hiopAlgFilterIPMNewton(hiop->mds);

  /* Error if model is not power balance */
  ierr = PetscStrcmp(opflow->modelname,OPFLOWMODEL_PBPOL,&flg);CHKERRQ(ierr);
  if(!flg) {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Only power balance model supported with HIOP solver\nUse -opflow_model POWER_BALANCE_POLAR\n");
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverSolve_HIOP(OPFLOW opflow)
{
  OPFLOWSolver_HIOP  hiop=(OPFLOWSolver_HIOP)opflow->solver;

  PetscFunctionBegin;

  hiop->status = hiop->solver->run();
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConvergenceStatus_HIOP(OPFLOW opflow,PetscBool *status)
{
  OPFLOWSolver_HIOP hiop = (OPFLOWSolver_HIOP)opflow->solver;

  PetscFunctionBegin;
  if(hiop->status < 3) *status = PETSC_TRUE; /* See hiopInterface.hpp. The first three denote convergence */
  else *status = PETSC_FALSE;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetObjective_HIOP(OPFLOW opflow,PetscReal *obj)
{
  PetscFunctionBegin;
  *obj = opflow->obj;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetSolution_HIOP(OPFLOW opflow,Vec *X)
{
  PetscFunctionBegin;
  *X = opflow->X;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConstraints_HIOP(OPFLOW opflow,Vec *G)
{
  PetscFunctionBegin;
  *G = opflow->G;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConstraintMultipliers_HIOP(OPFLOW opflow,Vec *Lambda)
{
  PetscFunctionBegin;
  *Lambda = opflow->Lambda;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverDestroy_HIOP(OPFLOW opflow)
{
  PetscErrorCode     ierr;
  OPFLOWSolver_HIOP  hiop=(OPFLOWSolver_HIOP)opflow->solver;

  PetscFunctionBegin;

  delete hiop->mds;
  delete hiop->nlp;

  ierr = PetscFree(hiop);CHKERRQ(ierr);

#ifdef HIOP_USE_MAGMA
  magma_finalize();
#endif

  PetscFunctionReturn(0);
}


extern "C" {

PetscErrorCode OPFLOWSolverCreate_HIOP(OPFLOW opflow)
{
  PetscErrorCode ierr;
  OPFLOWSolver_HIOP hiop;
  
  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&hiop);CHKERRQ(ierr);

  opflow->solver = hiop;
  
  opflow->solverops.setup = OPFLOWSolverSetUp_HIOP;
  opflow->solverops.solve = OPFLOWSolverSolve_HIOP;
  opflow->solverops.destroy = OPFLOWSolverDestroy_HIOP;
  opflow->solverops.getobjective = OPFLOWSolverGetObjective_HIOP;
  opflow->solverops.getconvergencestatus = OPFLOWSolverGetConvergenceStatus_HIOP;
  opflow->solverops.getsolution = OPFLOWSolverGetSolution_HIOP;
  opflow->solverops.getconstraints = OPFLOWSolverGetConstraints_HIOP;
  opflow->solverops.getconstraintmultipliers = OPFLOWSolverGetConstraintMultipliers_HIOP;

#ifdef HIOP_USE_MAGMA
  magma_init();
#endif
  PetscFunctionReturn(0);
}

} // End of extern "C"

#endif
