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

  hiop->nlp = new OPFLOWSolverHIOP();
  hiop->nlp->opflow = opflow;
  hiop->mds = new hiop::hiopNlpMDS(*hiop->nlp);

  /* Set options */
  hiop->mds->options->SetStringValue("dualsUpdateType", "linear");
  hiop->mds->options->SetStringValue("dualsInitialization", "zero");

  hiop->mds->options->SetStringValue("Hessian", "analytical_exact");
  hiop->mds->options->SetStringValue("KKTLinsys", "xdycyd");
  hiop->mds->options->SetStringValue("compute_mode", "cpu");

  hiop->mds->options->SetIntegerValue("verbosity_level", 3);
  hiop->mds->options->SetNumericValue("mu0", 1e-1);

  hiop->solver = new hiop::hiopAlgFilterIPMNewton(hiop->mds);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverSolve_HIOP(OPFLOW opflow)
{
  PetscErrorCode     ierr;
  OPFLOWSolver_HIOP  hiop=(OPFLOWSolver_HIOP)opflow->solver;

  PetscFunctionBegin;

  hiop->solver->run();
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
