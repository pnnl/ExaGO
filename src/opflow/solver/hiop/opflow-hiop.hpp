#include <scopflow_config.h>

#if defined(SCOPFLOW_HAVE_HIOP)

#ifndef OPFLOWHIOP_HPP
#define OPFFLOWHIOP_HPP
#include <opflow.h>
#include <hiopInterface.hpp>
#include <hiopMatrix.hpp>
#include <hiopNlpFormulation.hpp>
#include <hiopAlgFilterIPM.hpp>
#include <cstdlib>
#include <string>
#include <cstring>
#include <cstdio>

class OPFLOWHIOPInterface : public hiop::hiopInterfaceMDS
{
public:
  OPFLOWHIOPInterface(OPFLOW);


  ~OPFLOWHIOPInterface() {
    PetscFree(idxn2sd_map);
  }

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

  void solution_callback(hiop::hiopSolveStatus status,
                                 int n, const double* x,
                                 const double* z_L,
                                 const double* z_U,
                                 int m, const double* g,
                                 const double* lambda,
			         double obj_value);

  void naturaltospdense(const double*,double*);

  void spdensetonatural(const double*,double*);
private:
  OPFLOW   opflow;
  PetscInt nxsparse,nxdense;
  /* HIOP uses an ordering of variables where the sparse variablesd
     are ordered first followed by the dense variables. OPFLOW uses
     a different (natural) ordering where variables are ordered by
     buses. idxn2sd_map is an index mapping to go back and forth between the sparse-dense ordering to natural ordering. 
     
     We assume that all generator powers (active and reactive power)
     are sparse variables and all voltage variables are dense variables
  */
  PetscInt *idxn2sd_map;
};

typedef struct _p_OPFLOWSolver_HIOP *OPFLOWSolver_HIOP;

struct _p_OPFLOWSolver_HIOP {
  
  OPFLOWHIOPInterface           *nlp;
  hiop::hiopSolveStatus         status;
  hiop::hiopNlpMDS              *mds;
  hiop::hiopAlgFilterIPMNewton  *solver;

};

#endif
#endif
