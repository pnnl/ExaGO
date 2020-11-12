
#include <exago_config.h>
#if defined(EXAGO_HAVE_HIOP)
#include <private/opflowimpl.h>
#include "opflow-hiopnew.hpp"

typedef enum {
  AUTO = 0, CPU = 1,HYBRID = 2
}HIOPNEWComputeMode;
const char* HIOPNEWComputeModeChoices[] = {"auto","cpu","hybrid","HIOPNEWComputeModeChoices","",0};

/* Converts an array xin in natural ordering to an array xout in sparse-dense
   ordering
*/
void OPFLOWHIOPNEWInterface::naturaltospdense(const double *xin,double *xout)
{
  int i;

  for(i=0; i < opflow->nx; i++) {
    xout[idxn2sd_map[i]] = xin[i];
  }
}

/* Converts an array xin in sparse dense ordering to an array xout in natural
   ordering
*/
void OPFLOWHIOPNEWInterface::spdensetonatural(const double *xin,double *xout)
{
  int i;

  for(i=0; i < opflow->nx; i++) {
    xout[i] = xin[idxn2sd_map[i]];
  }
}

OPFLOWHIOPNEWInterface::OPFLOWHIOPNEWInterface(OPFLOW opflowin) 
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

bool OPFLOWHIOPNEWInterface::get_prob_sizes(long long& n, long long& m)
{ 
  n = opflow->nx;
  m = opflow->ncon;
  return true; 
}

bool OPFLOWHIOPNEWInterface::get_vars_info(const long long& n, double *xlow, double* xupp, NonlinearityType* type)
{
  PetscInt       i;
  PetscErrorCode ierr;
  //  PetscPrintf(MPI_COMM_SELF,"Enter get_vars_info\n");

  ierr = (*opflow->modelops.setvariableboundsarray)(opflow,xlow,xupp);CHKERRQ(ierr);
    
  for(i=0; i < n; i++) {
    type[i] = hiopNonlinear;
  }    
  //  PetscPrintf(MPI_COMM_SELF,"Exit get_vars_info\n");
  return true;
}

bool OPFLOWHIOPNEWInterface::get_cons_info(const long long& m, double* clow, double* cupp, NonlinearityType* type)
{
  PetscInt i;
  PetscErrorCode ierr;
  //  PetscPrintf(MPI_COMM_SELF,"Enter get_cons_info \n");

  ierr = (*opflow->modelops.setconstraintboundsarray)(opflow,clow,cupp);CHKERRQ(ierr);

  for(i=0; i < m; i++) type[i] = hiopNonlinear;

  //  PetscPrintf(MPI_COMM_SELF,"Exit get_cons_info \n");

  return true;
}

bool OPFLOWHIOPNEWInterface::get_sparse_dense_blocks_info(int& nx_sparse, int& nx_dense,
				  int& nnz_sparse_Jace, int& nnz_sparse_Jaci,
				  int& nnz_sparse_Hess_Lagr_SS, int& nnz_sparse_Hess_Lagr_SD)
{
  //  PetscPrintf(MPI_COMM_SELF,"Enter sparse_dense_blocks_info \n");

  nx_sparse = nxsparse;
  nx_dense  = nxdense;

  nnz_sparse_Jace = nnz_sparse_Hess_Lagr_SS = nxsparse;
  nnz_sparse_Jaci = 0; 
  nnz_sparse_Hess_Lagr_SD = 0;

  //  PetscPrintf(MPI_COMM_SELF,"Enter sparse_dense_blocks_info \n");

  return true;
}

bool OPFLOWHIOPNEWInterface::eval_f(const long long& n, const double* x, bool new_x, double& obj_value)
{
  PetscErrorCode ierr;

  //  PetscPrintf(MPI_COMM_SELF,"Enter eval_f \n");

  obj_value = 0.0;

  /* Compute objective */
  ierr = OPFLOWComputeObjectiveArray(opflow,x,&obj_value);CHKERRQ(ierr);

  //  PetscPrintf(MPI_COMM_SELF,"Exit eval_f \n");

  return true;
}

bool OPFLOWHIOPNEWInterface::eval_cons(const long long& n, const long long& m, 
	       const long long& num_cons, const long long* idx_cons,  
	       const double* x, bool new_x, double* cons)
{
  PetscErrorCode ierr;

  //  PetscPrintf(MPI_COMM_SELF,"Enter eval_cons \n");

  if(!num_cons) return true;

  if(idx_cons[0] == 0) {
    /* Equality constaints */
    ierr = OPFLOWComputeEqualityConstraintsArray(opflow,x,cons);CHKERRQ(ierr);
  }

  if(idx_cons[0] == opflow->nconeq && opflow->nconineq) {
    /* Inequality constraints */
    ierr = OPFLOWComputeInequalityConstraintsArray(opflow,x,cons);CHKERRQ(ierr);
  }

  //  PetscPrintf(MPI_COMM_SELF,"Exit eval_cons \n");
  return true;
}

bool OPFLOWHIOPNEWInterface::eval_grad_f(const long long& n, const double* x, bool new_x, double* gradf)
{
  PetscErrorCode ierr;
  //  PetscPrintf(MPI_COMM_SELF,"Enter eval_grad_f \n");

  ierr = OPFLOWComputeGradientArray(opflow,x,gradf);CHKERRQ(ierr);

  //  PetscPrintf(MPI_COMM_SELF,"Exit eval_grad_f \n");

  return true;
}

bool OPFLOWHIOPNEWInterface::eval_Jac_cons(const long long& n, const long long& m, 
		   const long long& num_cons, const long long* idx_cons,
		   const double* x, bool new_x,
		   const long long& nsparse, const long long& ndense, 
		   const int& nnzJacS, int* iJacS, int* jJacS, double* MJacS, 
		   double** JacD)
{
  PetscErrorCode ierr;
  //  PetscPrintf(MPI_COMM_SELF,"Enter eval_Jac_cons \n");

  if(!num_cons) return true;

  /* Equality constraints */
  if(idx_cons[0] == 0) {
    //    PetscPrintf(MPI_COMM_SELF,"Came here eq. \n");

    ierr = PetscLogEventBegin(opflow->eqconsjaclogger,0,0,0,0);CHKERRQ(ierr);
    /* Sparse Jacobian */
    ierr = (*opflow->modelops.computesparsejacobianhiop)(opflow,iJacS,jJacS,MJacS);CHKERRQ(ierr);

    /* Dense equality constraint Jacobian */
    ierr = (*opflow->modelops.computedenseequalityconstraintjacobianhiop)(opflow,x,JacD);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->eqconsjaclogger,0,0,0,0);CHKERRQ(ierr);

  } else {
    /* Dense inequality constraint Jacobian */
    ierr = PetscLogEventBegin(opflow->ineqconsjaclogger,0,0,0,0);CHKERRQ(ierr);
    if(opflow->nconineq) {
      ierr = (*opflow->modelops.computedenseinequalityconstraintjacobianhiop)(opflow,x,JacD);CHKERRQ(ierr);
    }
    ierr = PetscLogEventEnd(opflow->ineqconsjaclogger,0,0,0,0);CHKERRQ(ierr);
  }
  //  PetscPrintf(MPI_COMM_SELF,"Exit eval_Jac_cons \n");

  return true;
}

bool OPFLOWHIOPNEWInterface::eval_Hess_Lagr(const long long& n, const long long& m, 
		    const double* x, bool new_x, const double& obj_factor,
		    const double* lambda, bool new_lambda,
		    const long long& nsparse, const long long& ndense, 
		    const int& nnzHSS, int* iHSS, int* jHSS, double* MHSS, 
		    double** HDD,
		    int& nnzHSD, int* iHSD, int* jHSD, double* MHSD)
{
  PetscErrorCode ierr;
  //  PetscPrintf(MPI_COMM_SELF,"Enter eval_Hess_Lagr \n");

  ierr = PetscLogEventBegin(opflow->hesslogger,0,0,0,0);CHKERRQ(ierr);
  /* Compute sparse hessian */    
  ierr = (*opflow->modelops.computesparsehessianhiop)(opflow,x,iHSS,jHSS,MHSS);CHKERRQ(ierr);

  /* Compute dense hessian */
  ierr = (*opflow->modelops.computedensehessianhiop)(opflow,x,lambda,HDD);CHKERRQ(ierr);

  ierr = PetscLogEventEnd(opflow->hesslogger,0,0,0,0);CHKERRQ(ierr);

  //  PetscPrintf(MPI_COMM_SELF,"Exit eval_Hess_Lagr \n");

  return true;
}

bool OPFLOWHIOPNEWInterface::get_starting_point(const long long& global_n, double* x0)
{
  PetscErrorCode ierr;
  //  PetscPrintf(MPI_COMM_SELF,"Enter get_starting_point \n");

  ierr = (*opflow->modelops.setinitialguessarray)(opflow,x0);CHKERRQ(ierr);

  //  PetscPrintf(MPI_COMM_SELF,"Exit get_starting_point \n");

  return true;
}

void OPFLOWHIOPNEWInterface::solution_callback(hiop::hiopSolveStatus status,
							 int n, const double* xsol,
							 const double* z_L,
							 const double* z_U,
							 int m, const double* gsol,
							 const double* lamsol,
							 double obj_value)
{
  PetscErrorCode    ierr;
  OPFLOWSolver_HIOPNEW hiop=(OPFLOWSolver_HIOPNEW)opflow->solver;
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

PetscErrorCode OPFLOWSolverSetUp_HIOPNEW(OPFLOW opflow)
{
  PetscErrorCode ierr;
  OPFLOWSolver_HIOPNEW hiop=(OPFLOWSolver_HIOPNEW)opflow->solver;
  PetscBool         flg1,flg2;
  PetscReal         tol=1e-6;
  HIOPNEWComputeMode   compute_mode=AUTO;
  int               verbose_level=0;

  PetscFunctionBegin;

  hiop->nlp = new OPFLOWHIOPNEWInterface(opflow);
  hiop->mds = new hiop::hiopNlpMDS(*hiop->nlp);

  ierr = PetscOptionsBegin(opflow->comm->type,NULL,"HIOP options",NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnum("-hiop_compute_mode","Type of compute mode","",HIOPNEWComputeModeChoices,(PetscEnum)compute_mode,(PetscEnum*)&compute_mode,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-hiop_tolerance","HIOP solver tolerance","",tol,&tol,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-hiop_verbosity_level","HIOP verbosity level (Integer 0 to 12)","",verbose_level,&verbose_level,NULL);CHKERRQ(ierr);
  PetscOptionsEnd();

  hiop->mds->options->SetStringValue("dualsUpdateType", "linear");
  hiop->mds->options->SetStringValue("dualsInitialization", "zero");
  hiop->mds->options->SetStringValue("fixed_var", "relax");

  hiop->mds->options->SetStringValue("Hessian", "analytical_exact");
  hiop->mds->options->SetStringValue("KKTLinsys", "xdycyd");
  hiop->mds->options->SetStringValue("compute_mode", HIOPNEWComputeModeChoices[compute_mode]);

  hiop->mds->options->SetIntegerValue("verbosity_level", verbose_level);
  hiop->mds->options->SetNumericValue("mu0", 1e-1);
  hiop->mds->options->SetNumericValue("tolerance", tol);
#if defined(HIOP_USE_RAJA)
  hiop::LinearAlgebraFactory::set_mem_space("um");
  hiop->mds->options->SetStringValue("mem_space","um");
#endif
  //  ierr = PetscPrintf(MPI_COMM_SELF,"Came in OPFLOWSetUp\n");CHKERRQ(ierr);
  hiop->solver = new hiop::hiopAlgFilterIPMNewton(hiop->mds);

  /* Error if model is not power balance hiop or power balance raja hiop */
  ierr = PetscStrcmp(opflow->modelname,OPFLOWMODEL_PBPOLHIOP,&flg1);CHKERRQ(ierr);
#if defined(EXAGO_HAVE_RAJA)
  ierr = PetscStrcmp(opflow->modelname,OPFLOWMODEL_PBPOLRAJAHIOP,&flg2);CHKERRQ(ierr);
#endif
  if(!flg1 && !flg2) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"%s opflow model not supported with HIOP\n",opflow->modelname);
    PetscFunctionReturn(1);
    exit(0);
  }
  //  ierr = PetscPrintf(MPI_COMM_SELF,"Exit OPFLOWSetUp\n");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverSolve_HIOPNEW(OPFLOW opflow)
{
  OPFLOWSolver_HIOPNEW  hiop=(OPFLOWSolver_HIOPNEW)opflow->solver;

  PetscFunctionBegin;

  hiop->status = hiop->solver->run();
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConvergenceStatus_HIOPNEW(OPFLOW opflow,PetscBool *status)
{
  OPFLOWSolver_HIOPNEW hiop = (OPFLOWSolver_HIOPNEW)opflow->solver;

  PetscFunctionBegin;
  if(hiop->status < 3) *status = PETSC_TRUE; /* See hiopInterface.hpp. The first three denote convergence */
  else *status = PETSC_FALSE;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetObjective_HIOPNEW(OPFLOW opflow,PetscReal *obj)
{
  PetscFunctionBegin;
  *obj = opflow->obj;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetSolution_HIOPNEW(OPFLOW opflow,Vec *X)
{
  PetscFunctionBegin;
  *X = opflow->X;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConstraints_HIOPNEW(OPFLOW opflow,Vec *G)
{
  PetscFunctionBegin;
  *G = opflow->G;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverGetConstraintMultipliers_HIOPNEW(OPFLOW opflow,Vec *Lambda)
{
  PetscFunctionBegin;
  *Lambda = opflow->Lambda;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverDestroy_HIOPNEW(OPFLOW opflow)
{
  PetscErrorCode     ierr;
  OPFLOWSolver_HIOPNEW  hiop=(OPFLOWSolver_HIOPNEW)opflow->solver;

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

PetscErrorCode OPFLOWSolverCreate_HIOPNEW(OPFLOW opflow)
{
  PetscErrorCode ierr;
  OPFLOWSolver_HIOPNEW hiop;
  
  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&hiop);CHKERRQ(ierr);

  opflow->solver = hiop;
  
  opflow->solverops.setup = OPFLOWSolverSetUp_HIOPNEW;
  opflow->solverops.solve = OPFLOWSolverSolve_HIOPNEW;
  opflow->solverops.destroy = OPFLOWSolverDestroy_HIOPNEW;
  opflow->solverops.getobjective = OPFLOWSolverGetObjective_HIOPNEW;
  opflow->solverops.getconvergencestatus = OPFLOWSolverGetConvergenceStatus_HIOPNEW;
  opflow->solverops.getsolution = OPFLOWSolverGetSolution_HIOPNEW;
  opflow->solverops.getconstraints = OPFLOWSolverGetConstraints_HIOPNEW;
  opflow->solverops.getconstraintmultipliers = OPFLOWSolverGetConstraintMultipliers_HIOPNEW;

#ifdef HIOP_USE_MAGMA
  magma_init();
#endif
  PetscFunctionReturn(0);
}

} // End of extern "C"

#endif
