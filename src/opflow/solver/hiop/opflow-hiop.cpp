
#include <exago_config.h>
#if defined(EXAGO_HAVE_HIOP)
#include <private/opflowimpl.h>
#include "opflow-hiop.hpp"

typedef enum {
  AUTO = 0, CPU = 1,HYBRID = 2, GPU = 3
}HIOPComputeMode;
const char* HIOPComputeModeChoices[] = {"auto","cpu","hybrid","gpu","HIOPComputeModeChoices","",0};

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
  //  PetscPrintf(MPI_COMM_SELF,"Enter get_vars_info\n");

  ierr = (*opflow->modelops.setvariableboundsarray)(opflow,xlow,xupp);CHKERRQ(ierr);
    
  for(i=0; i < n; i++) {
    type[i] = hiopNonlinear;
  }    
  //  PetscPrintf(MPI_COMM_SELF,"Exit get_vars_info\n");
  return true;
}

bool OPFLOWHIOPInterface::get_cons_info(const long long& m, double* clow, double* cupp, NonlinearityType* type)
{
  PetscInt i;
  PetscErrorCode ierr;
  //  PetscPrintf(MPI_COMM_SELF,"Enter get_cons_info \n");

  ierr = (*opflow->modelops.setconstraintboundsarray)(opflow,clow,cupp);CHKERRQ(ierr);

  for(i=0; i < m; i++) type[i] = hiopNonlinear;

  //  PetscPrintf(MPI_COMM_SELF,"Exit get_cons_info \n");

  return true;
}

bool OPFLOWHIOPInterface::get_sparse_dense_blocks_info(int& nx_sparse, int& nx_dense,
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

bool OPFLOWHIOPInterface::eval_f(const long long& n, const double* x, bool new_x, double& obj_value)
{
  PetscErrorCode ierr;

  //  PetscPrintf(MPI_COMM_SELF,"Enter eval_f \n");

  obj_value = 0.0;

  /* Compute objective */
  ierr = (*opflow->modelops.computeobjectivearray)(opflow,x,&obj_value);CHKERRQ(ierr);

  //  PetscPrintf(MPI_COMM_SELF,"Exit eval_f \n");

  return true;
}

bool OPFLOWHIOPInterface::eval_cons(const long long& n, const long long& m, 
	       const long long& num_cons, const long long* idx_cons,  
	       const double* x, bool new_x, double* cons)
{
  PetscErrorCode ierr;

  //  PetscPrintf(MPI_COMM_SELF,"Enter eval_cons \n");

  if(!num_cons) return true;

  if(idx_cons[0] == 0) {
    /* Equality constaints */
    ierr = (*opflow->modelops.computeequalityconstraintsarray)(opflow,x,cons);CHKERRQ(ierr);
  }

  if(idx_cons[0] == opflow->nconeq && opflow->nconineq) {
    /* Inequality constraints */
    ierr = (*opflow->modelops.computeinequalityconstraintsarray)(opflow,x,cons);CHKERRQ(ierr);
  }

  //  PetscPrintf(MPI_COMM_SELF,"Exit eval_cons \n");
  return true;
}

bool OPFLOWHIOPInterface::eval_grad_f(const long long& n, const double* x, bool new_x, double* gradf)
{
  PetscErrorCode ierr;
  //  PetscPrintf(MPI_COMM_SELF,"Enter eval_grad_f \n");

  ierr = (*opflow->modelops.computegradientarray)(opflow,x,gradf);CHKERRQ(ierr);

  //  PetscPrintf(MPI_COMM_SELF,"Exit eval_grad_f \n");

  return true;
}

bool OPFLOWHIOPInterface::eval_Jac_cons(const long long& n, const long long& m, 
		   const long long& num_cons, const long long* idx_cons,
		   const double* x, bool new_x,
		   const long long& nsparse, const long long& ndense, 
		   const int& nnzJacS, int* iJacS, int* jJacS, double* MJacS, 
		   double* JacD)
{
  PetscErrorCode ierr;
  //  PetscPrintf(MPI_COMM_SELF,"Enter eval_Jac_cons \n");

  if(!num_cons) return true;

  /* Equality constraints */
  if(idx_cons[0] == 0) {
    //    PetscPrintf(MPI_COMM_SELF,"Came here eq. \n");

    /* Sparse Jacobian */
    ierr = (*opflow->modelops.computesparsejacobianhiop)(opflow,iJacS,jJacS,MJacS);CHKERRQ(ierr);

    /* Dense equality constraint Jacobian */
    ierr = (*opflow->modelops.computedenseequalityconstraintjacobianhiop)(opflow,x,JacD);CHKERRQ(ierr);
  } else {
    /* Dense inequality constraint Jacobian */
    //    PetscPrintf(MPI_COMM_SELF,"Came here ineq. \n");

    if(opflow->nconineq) {
      ierr = (*opflow->modelops.computedenseinequalityconstraintjacobianhiop)(opflow,x,JacD);CHKERRQ(ierr);
    }
  }
  //  PetscPrintf(MPI_COMM_SELF,"Exit eval_Jac_cons \n");

  return true;
}

bool OPFLOWHIOPInterface::eval_Hess_Lagr(const long long& n, const long long& m, 
		    const double* x, bool new_x, const double& obj_factor,
		    const double* lambda, bool new_lambda,
		    const long long& nsparse, const long long& ndense, 
		    const int& nnzHSS, int* iHSS, int* jHSS, double* MHSS, 
		    double* HDD,
		    int& nnzHSD, int* iHSD, int* jHSD, double* MHSD)
{
  PetscErrorCode ierr;
  //  PetscPrintf(MPI_COMM_SELF,"Enter eval_Hess_Lagr \n");

  /* Compute sparse hessian */    
  ierr = (*opflow->modelops.computesparsehessianhiop)(opflow,x,iHSS,jHSS,MHSS);CHKERRQ(ierr);

  /* Compute dense hessian */
  ierr = (*opflow->modelops.computedensehessianhiop)(opflow,x,lambda,HDD);CHKERRQ(ierr);

  //  PetscPrintf(MPI_COMM_SELF,"Exit eval_Hess_Lagr \n");

  return true;
}

bool OPFLOWHIOPInterface::get_starting_point(const long long& global_n, double* x0)
{
  PetscErrorCode ierr;
  //  PetscPrintf(MPI_COMM_SELF,"Enter get_starting_point \n");

  ierr = (*opflow->modelops.setinitialguessarray)(opflow,x0);CHKERRQ(ierr);

  //  PetscPrintf(MPI_COMM_SELF,"Exit get_starting_point \n");

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
  opflow->numits = n;

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
  PetscBool         ismodelpbpolhiop,ismodelpbpolrajahiop;
  HIOPComputeMode   compute_mode=AUTO;
  int               verbose_level=0;

  PetscFunctionBegin;

  hiop->nlp = new OPFLOWHIOPInterface(opflow);
  hiop->mds = new hiop::hiopNlpMDS(*hiop->nlp);

  ierr = PetscOptionsBegin(opflow->comm->type,NULL,"HIOP options",NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnum("-hiop_compute_mode","Type of compute mode","",HIOPComputeModeChoices,(PetscEnum)compute_mode,(PetscEnum*)&compute_mode,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-hiop_verbosity_level","HIOP verbosity level (Integer 0 to 12)","",verbose_level,&verbose_level,NULL);CHKERRQ(ierr);
  PetscOptionsEnd();

  hiop->mds->options->SetStringValue("dualsUpdateType", "linear");
  hiop->mds->options->SetStringValue("dualsInitialization", "zero");
  hiop->mds->options->SetStringValue("fixed_var", "relax");

  hiop->mds->options->SetStringValue("Hessian", "analytical_exact");
  hiop->mds->options->SetStringValue("KKTLinsys", "xdycyd");
  hiop->mds->options->SetStringValue("compute_mode", HIOPComputeModeChoices[compute_mode]);

  hiop->mds->options->SetIntegerValue("verbosity_level", verbose_level);
  hiop->mds->options->SetNumericValue("mu0", 1e-1);
  hiop->mds->options->SetNumericValue("tolerance", opflow->tolerance);

  /* Error if model is not power balance hiop or power balance raja hiop */
  ierr = PetscStrcmp(opflow->modelname,OPFLOWMODEL_PBPOLHIOP,&ismodelpbpolhiop);CHKERRQ(ierr);
#if defined(EXAGO_HAVE_RAJA)
  ierr = PetscStrcmp(opflow->modelname,OPFLOWMODEL_PBPOLRAJAHIOP,&ismodelpbpolrajahiop);CHKERRQ(ierr);
#endif
  if(!ismodelpbpolhiop && !ismodelpbpolrajahiop) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"%s opflow model not supported with HIOP\n",opflow->modelname);
    PetscFunctionReturn(1);
    exit(0);
  }

#if defined(HIOP_USE_RAJA)
  if(ismodelpbpolrajahiop) {
#ifndef EXAGO_HAVE_GPU
    hiop->mds->options->SetStringValue("mem_space","host");
#else
    // TODO: replace this with "device" when supported by HiOp
    hiop->mds->options->SetStringValue("mem_space","um");
#endif
  }
#endif

  //  ierr = PetscPrintf(MPI_COMM_SELF,"Came in OPFLOWSetUp\n");CHKERRQ(ierr);
  hiop->solver = new hiop::hiopAlgFilterIPMNewton(hiop->mds);

  //  ierr = PetscPrintf(MPI_COMM_SELF,"Exit OPFLOWSetUp\n");CHKERRQ(ierr);

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

#ifdef EXAGO_HAVE_GPU
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

#ifdef EXAGO_HAVE_GPU
  magma_init();
#endif
  PetscFunctionReturn(0);
}

} // End of extern "C"

#endif
