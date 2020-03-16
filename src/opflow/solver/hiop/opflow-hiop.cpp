#if defined(SCOPFLOW_HAVE_HIOP)
#include <private/opflowimpl.h>
#include "opflow-hiop.hpp"


bool OPFLOWSolverHIOP::get_prob_sizes(long long& n, long long& m)
{ 
  n = opflow->nx;
  m = opflow->nconeq;
  return true; 
}

bool OPFLOWSolverHIOP::get_vars_info(const long long& n, double *xlow, double* xupp, NonlinearityType* type)
{
  PetscInt       i;
  PetscErrorCode ierr;

  ierr = VecPlaceArray(opflow->Xl,xlow);CHKERRQ(ierr);
  ierr = VecPlaceArray(opflow->Xu,xupp);CHKERRQ(ierr);

  ierr = (*opflow->formops.setvariablebounds)(opflow,opflow->Xl,opflow->Xu);CHKERRQ(ierr);

  ierr = VecResetArray(opflow->Xl);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->Xu);CHKERRQ(ierr);

  for(i=0; i < n; i++) type[i] = hiopNonlinear;
  return true;
}

bool OPFLOWSolverHIOP::get_cons_info(const long long& m, double* clow, double* cupp, NonlinearityType* type)
{
  PetscInt i;
  PetscErrorCode ierr;

  ierr = VecPlaceArray(opflow->Gl,clow);CHKERRQ(ierr);
  ierr = VecPlaceArray(opflow->Gu,cupp);CHKERRQ(ierr);

  ierr = (*opflow->formops.setconstraintbounds)(opflow,opflow->Gl,opflow->Gu);CHKERRQ(ierr);

  ierr = VecResetArray(opflow->Gl);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->Gu);CHKERRQ(ierr);

  for(i=0; i < m; i++) type[i] = hiopNonlinear;

  return true;
}

bool OPFLOWSolverHIOP::get_sparse_dense_blocks_info(int& nx_sparse, int& nx_dense,
				  int& nnz_sparse_Jace, int& nnz_sparse_Jaci,
				  int& nnz_sparse_Hess_Lagr_SS, int& nnz_sparse_Hess_Lagr_SD)
{
  nx_sparse = 0;
  nx_dense  = opflow->nx;

  nnz_sparse_Jace = nnz_sparse_Jaci = nnz_sparse_Hess_Lagr_SS = nnz_sparse_Hess_Lagr_SD = 0;
  return true;
}

bool OPFLOWSolverHIOP::eval_f(const long long& n, const double* x, bool new_x, double& obj_value)
{
  PetscErrorCode ierr;

  obj_value = 0.0;
  ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
  ierr = (*opflow->formops.computeobjective)(opflow,opflow->X,&obj_value);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  return true;
}

bool OPFLOWSolverHIOP::eval_cons(const long long& n, const long long& m, 
	       const long long& num_cons, const long long* idx_cons,  
	       const double* x, bool new_x, double* cons)
{
  PetscErrorCode ierr;

  if(!num_cons) return true;
  ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
  ierr = VecPlaceArray(opflow->Ge,cons+idx_cons[0]);CHKERRQ(ierr);
  ierr = (*opflow->formops.computeequalityconstraints)(opflow,opflow->X,opflow->Ge);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->Ge);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->X);CHKERRQ(ierr);

  return true;
}

bool OPFLOWSolverHIOP::eval_grad_f(const long long& n, const double* x, bool new_x, double* gradf)
{
  PetscErrorCode ierr;

  ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
  ierr = VecPlaceArray(opflow->gradobj,gradf);CHKERRQ(ierr);
  ierr = VecSet(opflow->gradobj,0.0);CHKERRQ(ierr);
  ierr = (*opflow->formops.computegradient)(opflow,opflow->X,opflow->gradobj);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->gradobj);CHKERRQ(ierr);

  return true;
}

bool OPFLOWSolverHIOP::eval_Jac_cons(const long long& n, const long long& m, 
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

  if(!num_cons) return true;
  if(JacD != NULL) {
    ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
    
    /* Equality constraints Jacobian */
    ierr = (*opflow->formops.computeequalityconstraintjacobian)(opflow,opflow->X,opflow->Jac_Ge);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
    
    for(i=0; i < m; i++) {
      for(j=0; j < n; j++) JacD[i][j] = 0.0;
      ierr = MatGetRow(opflow->Jac_Ge,i,&ncols,&cols,&vals);CHKERRQ(ierr);
      for(k=0; k < ncols; k++) JacD[i][cols[k]] = vals[k];
      ierr = MatRestoreRow(opflow->Jac_Ge,i,&ncols,&cols,&vals);CHKERRQ(ierr);
    }
  }
  return true;
}

bool OPFLOWSolverHIOP::eval_Hess_Lagr(const long long& n, const long long& m, 
		    const double* x, bool new_x, const double& obj_factor,
		    const double* lambda, bool new_lambda,
		    const long long& nsparse, const long long& ndense, 
		    const int& nnzHSS, int* iHSS, int* jHSS, double* MHSS, 
		    double** HDD,
		    int& nnzHSD, int* iHSD, int* jHSD, double* MHSD)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k;
  PetscInt       ncols;
  const PetscInt *cols;
  const PetscScalar *vals;

  opflow->obj_factor = obj_factor;
  if(HDD != NULL) {
    ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Lambdae,lambda);CHKERRQ(ierr);
    if(opflow->Nconineq) {
      ierr = VecPlaceArray(opflow->Lambdai,lambda+opflow->nconeq);CHKERRQ(ierr);
    }
    
    /* Compute Hessian */
    ierr = (*opflow->formops.computehessian)(opflow,opflow->X,opflow->Lambdae,opflow->Lambdai,opflow->Hes);CHKERRQ(ierr);
    
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Lambdae);CHKERRQ(ierr);
    if(opflow->Nconineq) {
      ierr = VecResetArray(opflow->Lambdai);CHKERRQ(ierr);
    }

    for(i=0; i < n; i++) {
      for(j=0; j < n; j++) HDD[i][j] = 0.0;
      ierr = MatGetRow(opflow->Hes,i,&ncols,&cols,&vals);CHKERRQ(ierr);
      for(k=0; k < ncols; k++) HDD[i][cols[k]] = vals[k];
      ierr = MatRestoreRow(opflow->Hes,i,&ncols,&cols,&vals);CHKERRQ(ierr);
    }
  }

  return true;
}

bool OPFLOWSolverHIOP::get_starting_point(const long long& global_n, double* x0)
{
  PetscErrorCode ierr;

  ierr = VecPlaceArray(opflow->X,x0);CHKERRQ(ierr);
  /* Set initial guess */
  ierr = (*opflow->formops.setinitialguess)(opflow,opflow->X);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->X);CHKERRQ(ierr);

  return true;
}

PetscErrorCode OPFLOWSolverSetUp_HIOP(OPFLOW opflow)
{
  PetscErrorCode ierr;
  OPFLOWSolver_HIOP hiop=(OPFLOWSolver_HIOP)opflow->solver;
  PetscBool         flg;
  PetscFunctionBegin;

  hiop->nlp = new OPFLOWSolverHIOP();
  hiop->nlp->opflow = opflow;
  hiop->mds = new hiop::hiopNlpMDS(*hiop->nlp);

  /* Options set in hiop.options file */

  hiop->solver = new hiop::hiopAlgFilterIPMNewton(hiop->mds);

  /* Error if formulation is not power balance */
  ierr = PetscStrcmp(opflow->formulationname,OPFLOWFORMULATION_PBPOL,&flg);CHKERRQ(ierr);
  if(!flg) {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Only power balance form formulation supported with HIOP solver\nUse -opflow_formulation POWER_BALANCE_POLAR\n");
  }
  if(opflow->nconineq) {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Line flow inequality constraints are not supported with HIOP solver\n Run with option -opflow_ignore_lineflow_constraints\n");
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverSolve_HIOP(OPFLOW opflow)
{
  OPFLOWSolver_HIOP  hiop=(OPFLOWSolver_HIOP)opflow->solver;
  PetscScalar        obj_value;

  PetscFunctionBegin;

  hiop->status = hiop->solver->run();
  obj_value = hiop->solver->getObjective();
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverDestroy_HIOP(OPFLOW opflow)
{
  PetscErrorCode     ierr;
  OPFLOWSolver_HIOP  hiop=(OPFLOWSolver_HIOP)opflow->solver;

  PetscFunctionBegin;

  free(hiop->mds);
  free(hiop->nlp);

  ierr = PetscFree(hiop);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


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

  PetscFunctionReturn(0);
}

#endif
