#if defined(SCOPFLOW_HAVE_HIOP)
#include <private/opflowimpl.h>
#include "opflow-hiop.hpp"


bool OPFLOWSolverHIOP::get_prob_sizes(long long& n, long long& m)
{ 
  return true; 
}

bool OPFLOWSolverHIOP::get_vars_info(const long long& n, double *xlow, double* xupp, NonlinearityType* type)
{
  return true;
}

bool OPFLOWSolverHIOP::get_cons_info(const long long& m, double* clow, double* cupp, NonlinearityType* type)
{
  return true;
}

bool OPFLOWSolverHIOP::get_sparse_dense_blocks_info(int& nx_sparse, int& nx_dense,
				  int& nnz_sparse_Jace, int& nnz_sparse_Jaci,
				  int& nnz_sparse_Hess_Lagr_SS, int& nnz_sparse_Hess_Lagr_SD)
{
    return true;
}

bool OPFLOWSolverHIOP::eval_f(const long long& n, const double* x, bool new_x, double& obj_value)
{
  return true;
}

bool OPFLOWSolverHIOP::eval_cons(const long long& n, const long long& m, 
	       const long long& num_cons, const long long* idx_cons,  
	       const double* x, bool new_x, double* cons)
{
    return true;
}

bool OPFLOWSolverHIOP::eval_grad_f(const long long& n, const double* x, bool new_x, double* gradf)
{
  return true;
}

bool OPFLOWSolverHIOP::eval_Jac_cons(const long long& n, const long long& m, 
		   const long long& num_cons, const long long* idx_cons,
		   const double* x, bool new_x,
		   const long long& nsparse, const long long& ndense, 
		   const int& nnzJacS, int* iJacS, int* jJacS, double* MJacS, 
		   double** JacD)
{
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
  return true;
}

bool OPFLOWSolverHIOP::get_starting_point(const long long& global_n, double* x0)
{
  return true;
}

PetscErrorCode OPFLOWSolverSetUp_HIOP(OPFLOW opflow)
{
  PetscErrorCode ierr;
  OPFLOWSolver_HIOP hiop=(OPFLOWSolver_HIOP)opflow->solver;
  PetscFunctionBegin;

  hiop->mds = new hiop::hiopNlpMDS(hiop->nlp);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverCreate_HIOP(OPFLOW opflow)
{
  PetscErrorCode ierr;
  OPFLOWSolver_HIOP hiop;
  
  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&hiop);CHKERRQ(ierr);

  opflow->solver = hiop;
  hiop->nlp.opflow = opflow;

  //  opflow->solverops.setup = OPFLOWSolverSetUp_HIOP;
  //  opflow->solverops.solve = OPFLOWSolverSolve_HIOP;
  //  opflow->solverops.destroy = OPFLOWSolverDestroy_HIOP;

  PetscFunctionReturn(0);
}

#endif
