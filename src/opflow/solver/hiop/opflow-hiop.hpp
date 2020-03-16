#if defined(SCOPFLOW_HAVE_HIOP)

#ifndef OPFLOWHIOP_H
#define OPFFLOWHIOP_H
#include <opflow.h>
#include <hiopInterface.hpp>
#include <hiopMatrix.hpp>
#include <hiopNlpFormulation.hpp>
#include <hiopAlgFilterIPM.hpp>
#include <cstdlib>
#include <string>
#include <cstring>
#include <cstdio>

class OPFLOWSolverHIOP : public hiop::hiopInterfaceMDS
{
public:
  OPFLOWSolverHIOP() {};
  OPFLOW opflow;

  bool get_prob_sizes(long long& n, long long& m);

  bool get_vars_info(const long long& n, double *xlow, double* xupp, NonlinearityType* type);


  bool get_cons_info(const long long& m, double* clow, double* cupp, NonlinearityType* type);


  bool get_sparse_dense_blocks_info(int& nx_sparse, int& nx_dense,
				    int& nnz_sparse_Jace, int& nnz_sparse_Jaci,
				    int& nnz_sparse_Hess_Lagr_SS, int& nnz_sparse_Hess_Lagr_SD);


  bool eval_f(const long long& n, const double* x, bool new_x, double& obj_value);


  bool eval_cons(const long long& n, const long long& m, 
		 const long long& num_cons, const long long* idx_cons,  
		 const double* x, bool new_x, double* cons);


  bool eval_grad_f(const long long& n, const double* x, bool new_x, double* gradf);

  bool eval_Jac_cons(const long long& n, const long long& m, 
		     const long long& num_cons, const long long* idx_cons,
		     const double* x, bool new_x,
		     const long long& nsparse, const long long& ndense, 
		     const int& nnzJacS, int* iJacS, int* jJacS, double* MJacS, 
		     double** JacD);

  bool eval_Hess_Lagr(const long long& n, const long long& m, 
		      const double* x, bool new_x, const double& obj_factor,
		      const double* lambda, bool new_lambda,
		      const long long& nsparse, const long long& ndense, 
		      const int& nnzHSS, int* iHSS, int* jHSS, double* MHSS, 
		      double** HDD,
		      int& nnzHSD, int* iHSD, int* jHSD, double* MHSD);

  bool get_starting_point(const long long& global_n, double* x0);

};

typedef struct _p_OPFLOWSolver_HIOP *OPFLOWSolver_HIOP;

struct _p_OPFLOWSolver_HIOP {
  
  OPFLOWSolverHIOP              *nlp;
  hiop::hiopSolveStatus         status;
  hiop::hiopNlpMDS              *mds;
  hiop::hiopAlgFilterIPMNewton  *solver;

};

#endif
#endif
