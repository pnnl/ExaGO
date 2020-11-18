/*
   This is the first HiOP interface for OPFLOW. It only works with power_balance_polar model. This
   implementation is deprecated. Use the other hiop implementation
*/
#include <exago_config.h>
#if defined(EXAGO_HAVE_HIOP)
#include <private/opflowimpl.h>
#include "opflow-hiopold.hpp"

typedef enum {
  AUTO = 0, CPU = 1,HYBRID = 2
}HIOPOLDComputeMode;
const char* HIOPOLDComputeModeChoices[] = {"auto","cpu","hybrid","HIOPOLDComputeModeChoices","",0};

/* Converts an array xin in natural ordering to an array xout in sparse-dense
   ordering
*/
void OPFLOWHIOPOLDInterface::naturaltospdense(const double *xin,double *xout)
{
  int i;

  for(i=0; i < opflow->nx; i++) {
    xout[idxn2sd_map[i]] = xin[i];
  }
}

/* Converts an array xin in sparse dense ordering to an array xout in natural
   ordering
*/
void OPFLOWHIOPOLDInterface::spdensetonatural(const double *xin,double *xout)
{
  int i;

  for(i=0; i < opflow->nx; i++) {
    xout[i] = xin[idxn2sd_map[i]];
  }
}

OPFLOWHIOPOLDInterface::OPFLOWHIOPOLDInterface(OPFLOW opflowin) 
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

bool OPFLOWHIOPOLDInterface::get_prob_sizes(long long& n, long long& m)
{ 
  n = opflow->nx;
  m = opflow->ncon;
  return true; 
}

bool OPFLOWHIOPOLDInterface::get_vars_info(const long long& n, double *xlow, double* xupp, NonlinearityType* type)
{
  PetscInt       i;
  PetscErrorCode ierr;
  const PetscScalar    *xl,*xu;

  ierr = (*opflow->modelops.setvariablebounds)(opflow,opflow->Xl,opflow->Xu);CHKERRQ(ierr);
    
  ierr = VecGetArrayRead(opflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArrayRead(opflow->Xu,&xu);CHKERRQ(ierr);

  naturaltospdense(xl,xlow);
  naturaltospdense(xu,xupp);

  for(i=0; i < n; i++) {
    type[i] = hiopNonlinear;
  }    

  ierr = VecRestoreArrayRead(opflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(opflow->Xu,&xu);CHKERRQ(ierr);  

  return true;
}

bool OPFLOWHIOPOLDInterface::get_cons_info(const long long& m, double* clow, double* cupp, NonlinearityType* type)
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

bool OPFLOWHIOPOLDInterface::get_sparse_dense_blocks_info(int& nx_sparse, int& nx_dense,
				  int& nnz_sparse_Jace, int& nnz_sparse_Jaci,
				  int& nnz_sparse_Hess_Lagr_SS, int& nnz_sparse_Hess_Lagr_SD)
{
  nx_sparse = nxsparse;
  nx_dense  = nxdense;

  nnz_sparse_Jace = nnz_sparse_Hess_Lagr_SS = nxsparse;
  nnz_sparse_Jaci = 0; 
  nnz_sparse_Hess_Lagr_SD = 0;
  return true;
}

bool OPFLOWHIOPOLDInterface::eval_f(const long long& n, const double* x, bool new_x, double& obj_value)
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

bool OPFLOWHIOPOLDInterface::eval_cons(const long long& n, const long long& m, 
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

bool OPFLOWHIOPOLDInterface::eval_grad_f(const long long& n, const double* x, bool new_x, double* gradf)
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

bool OPFLOWHIOPOLDInterface::eval_Jac_cons(const long long& n, const long long& m, 
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
	//	printf("%d:",i);
	for(k=0; k < ncols; k++) {
	  if(idxn2sd_map[cols[k]]-nxsparse >= 0) {
	    /* Dense variables */
	    dcol = idxn2sd_map[cols[k]] - nxsparse; /* Column number for dense variable in the dense block */
	    JacD[i][dcol] = vals[k];
	    //	    printf("(%d, %lf)",dcol,vals[k]);
	  }
	}
	//	printf("\n");
	ierr = MatRestoreRow(opflow->Jac_Ge,i,&ncols,&cols,&vals);CHKERRQ(ierr);
      }
    } else {
      
      /* Inequality constraints Jacobian */
      ierr = (*opflow->modelops.computeinequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Gi);CHKERRQ(ierr);

      for(i=0; i < opflow->nconineq; i++) {
	for(j=0; j < nxdense; j++) JacD[i][j] = 0.0;
	ierr = MatGetRow(opflow->Jac_Gi,i,&ncols,&cols,&vals);CHKERRQ(ierr);
	//	printf("%d:",i);
	for(k=0; k < ncols; k++) {
	  if(idxn2sd_map[cols[k]]-nxsparse >= 0) {
	    /* Dense variables */
	    dcol = idxn2sd_map[cols[k]] - nxsparse; /* Column number for dense variable in the dense block */
	    JacD[i][dcol] = vals[k];
	    //	    printf("(%d, %lf)",dcol,vals[k]);
	  }
	}
	//	printf("\n");
	ierr = MatRestoreRow(opflow->Jac_Gi,i,&ncols,&cols,&vals);CHKERRQ(ierr);
      }
    }
  }
  return true;
}

bool OPFLOWHIOPOLDInterface::eval_Hess_Lagr(const long long& n, const long long& m, 
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
	//	printf("%d:",dnct);
	for(k=0; k < ncols; k++) {
	  if(idxn2sd_map[cols[k]] >= nxsparse) {
	    dcol = idxn2sd_map[cols[k]] - nxsparse; /* Column number for dense variable in the dense block */
	    HDD[dnct][dcol] = vals[k];
	    //	    printf("(%d, %lf)",dcol, vals[k]);

	  }
	}
	//	printf("\n");
	ierr = MatRestoreRow(opflow->Hes,i,&ncols,&cols,&vals);CHKERRQ(ierr);

	dnct++;
      }
    }
  }

  return true;
}

bool OPFLOWHIOPOLDInterface::get_starting_point(const long long& global_n, double* x0)
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

void OPFLOWHIOPOLDInterface::solution_callback(hiop::hiopSolveStatus status,
							 int n, const double* xsol,
							 const double* z_L,
							 const double* z_U,
							 int m, const double* gsol,
							 const double* lamsol,
							 double obj_value)
{
  PetscErrorCode    ierr;
  OPFLOWSolver_HIOPOLD hiop=(OPFLOWSolver_HIOPOLD)opflow->solver;
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

PetscErrorCode OPFLOWSolverSetUp_HIOPOLD(OPFLOW opflow)
{
  PetscErrorCode    ierr;
  OPFLOWSolver_HIOPOLD hiop=(OPFLOWSolver_HIOPOLD)opflow->solver;
  PetscBool         flg1;
  PetscReal         tol=1e-6;
  HIOPOLDComputeMode   compute_mode=AUTO;
  int               verbose_level=0;

  PetscFunctionBegin;

  hiop->nlp = new OPFLOWHIOPOLDInterface(opflow);
  hiop->mds = new hiop::hiopNlpMDS(*hiop->nlp);

  ierr = PetscOptionsBegin(opflow->comm->type,NULL,"HIOP options",NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnum("-hiop_compute_mode","Type of compute mode","",HIOPOLDComputeModeChoices,(PetscEnum)compute_mode,(PetscEnum*)&compute_mode,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-hiop_tolerance","HIOP solver tolerance","",tol,&tol,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-hiop_verbosity_level","HIOP verbosity level (Integer 0 to 12)","",verbose_level,&verbose_level,NULL);CHKERRQ(ierr);
  PetscOptionsEnd();

  hiop->mds->options->SetStringValue("dualsUpdateType", "linear");
  hiop->mds->options->SetStringValue("dualsInitialization", "zero");
  hiop->mds->options->SetStringValue("fixed_var", "relax");

  hiop->mds->options->SetStringValue("Hessian", "analytical_exact");
  hiop->mds->options->SetStringValue("KKTLinsys", "xdycyd");
  hiop->mds->options->SetStringValue("compute_mode", HIOPOLDComputeModeChoices[compute_mode]);

  hiop->mds->options->SetIntegerValue("verbosity_level", verbose_level);
  hiop->mds->options->SetNumericValue("mu0", 1e-1);
  hiop->mds->options->SetNumericValue("tolerance", tol);

  hiop->solver = new hiop::hiopAlgFilterIPMNewton(hiop->mds);

  /* Error if model is not power balance hiop */
  ierr = PetscStrcmp(opflow->modelname,OPFLOWMODEL_PBPOL,&flg1);CHKERRQ(ierr);
  if(!flg1) {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Only power balance polar model allowed\n Run with -opflow_model POWER_BALANCE_POLAR\n");
    exit(1);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverSolve_HIOPOLD(OPFLOW opflow)
{
  OPFLOWSolver_HIOPOLD  hiop=(OPFLOWSolver_HIOPOLD)opflow->solver;

  PetscFunctionBegin;

  hiop->status = hiop->solver->run();
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConvergenceStatus_HIOPOLD(OPFLOW opflow,PetscBool *status)
{
  OPFLOWSolver_HIOPOLD hiop = (OPFLOWSolver_HIOPOLD)opflow->solver;

  PetscFunctionBegin;
  if(hiop->status < 3) *status = PETSC_TRUE; /* See hiopInterface.hpp. The first three denote convergence */
  else *status = PETSC_FALSE;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetObjective_HIOPOLD(OPFLOW opflow,PetscReal *obj)
{
  PetscFunctionBegin;
  *obj = opflow->obj;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetSolution_HIOPOLD(OPFLOW opflow,Vec *X)
{
  PetscFunctionBegin;
  *X = opflow->X;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConstraints_HIOPOLD(OPFLOW opflow,Vec *G)
{
  PetscFunctionBegin;
  *G = opflow->G;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConstraintMultipliers_HIOPOLD(OPFLOW opflow,Vec *Lambda)
{
  PetscFunctionBegin;
  *Lambda = opflow->Lambda;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverDestroy_HIOPOLD(OPFLOW opflow)
{
  PetscErrorCode     ierr;
  OPFLOWSolver_HIOPOLD  hiop=(OPFLOWSolver_HIOPOLD)opflow->solver;

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

PetscErrorCode OPFLOWSolverCreate_HIOPOLD(OPFLOW opflow)
{
  PetscErrorCode ierr;
  OPFLOWSolver_HIOPOLD hiop;
  
  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&hiop);CHKERRQ(ierr);

  opflow->solver = hiop;
  
  opflow->solverops.setup = OPFLOWSolverSetUp_HIOPOLD;
  opflow->solverops.solve = OPFLOWSolverSolve_HIOPOLD;
  opflow->solverops.destroy = OPFLOWSolverDestroy_HIOPOLD;
  opflow->solverops.getobjective = OPFLOWSolverGetObjective_HIOPOLD;
  opflow->solverops.getconvergencestatus = OPFLOWSolverGetConvergenceStatus_HIOPOLD;
  opflow->solverops.getsolution = OPFLOWSolverGetSolution_HIOPOLD;
  opflow->solverops.getconstraints = OPFLOWSolverGetConstraints_HIOPOLD;
  opflow->solverops.getconstraintmultipliers = OPFLOWSolverGetConstraintMultipliers_HIOPOLD;

#ifdef HIOP_USE_MAGMA
  magma_init();
#endif
  PetscFunctionReturn(0);
}

} // End of extern "C"

#endif
